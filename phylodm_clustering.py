#!/usr/bin/env python3
"""
PhyloDM Clustering Pipeline

This script performs clustering of sequences based on phylogenetic distances.
It can also incorporate taxonomy information to ensure clusters are consistent at a specified taxonomic level.
"""

import sys
import os
import argparse
import pandas as pd
import numpy as np
import dendropy
from phylodm import PhyloDM
import multiprocessing
import glob
import re


def create_dir(directory):
    """Create a directory if it doesn't exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)


def process_matrix(df, cutoff, output_file):
    """Process the distance matrix and run MCL clustering."""
    temp_df = df.copy()
    temp_df = pd.melt(temp_df.reset_index(), id_vars='genome', value_vars=list(df.columns.values), var_name='genome2', value_name='sim')
    with pd.option_context('display.max_rows', None, 'display.max_columns', 3):
        temp_df['sim'] = temp_df['sim'].apply(lambda x: 1-x)
        temp_df = temp_df[(temp_df.sim > cutoff) & (temp_df.genome != temp_df.genome2)]
        temp_df.to_csv(output_file, sep="\t", index=None)

    # Use the absolute path for the MCL output file
    mcl_output = f"{output_file}.mcl"
    os.system(f"mcl {output_file} --abc -o {mcl_output} -I 1.5")


def run_phylodm_parallel(nwktre, out3c, cutoffs):
    """Run PhyloDM to generate distance matrices and perform clustering in parallel."""
    print(f"Loading tree from {nwktre}...")
    tree = dendropy.Tree.get_from_path(nwktre, schema='newick')
    pdm = PhyloDM.load_from_dendropy(tree)
    dm = pdm.dm(norm=False)
    labels = pdm.taxa()
    labels = [x.replace(" ", "_") for x in labels]
    df = pd.DataFrame(data=dm, columns=labels, index=labels)
    df.index.name = 'genome'

    # Create the output directory for PhyloDM results
    phylodm_out_dir = f"{out3c}_phylodm_out"
    create_dir(phylodm_out_dir)
    print(f"Running clustering with cutoffs: {', '.join([str(c) for c in cutoffs])}")

    # Use sequential processing instead of parallel to avoid file path issues
    for cutoff in cutoffs:
        # Use os.path.basename to get just the filename part
        base_name = os.path.basename(out3c)
        output_file = os.path.join(phylodm_out_dir, f"{base_name}_{cutoff:.2f}.3c")
        process_matrix(df, cutoff, output_file)


def read_taxonomy_file(taxonomy_file):
    """Read taxonomy file and return a dictionary mapping sequence IDs to taxonomy strings."""
    taxonomy_dict = {}
    if taxonomy_file:
        try:
            df = pd.read_csv(taxonomy_file, sep='\t', header=None)
            if len(df.columns) < 2:
                print(f"Warning: Taxonomy file {taxonomy_file} does not have at least 2 columns. Expected format: <sequence_id>\\t<taxonomy_string>")
                return {}

            for _, row in df.iterrows():
                seq_id = row[0]
                tax_string = row[1]
                taxonomy_dict[seq_id] = tax_string

            print(f"Loaded taxonomy information for {len(taxonomy_dict)} sequences")
        except Exception as e:
            print(f"Error reading taxonomy file {taxonomy_file}: {e}")
            return {}
    return taxonomy_dict


def get_taxonomy_level(tax_string, level):
    """Extract taxonomy at the specified level from a taxonomy string."""
    if not tax_string:
        return None

    parts = tax_string.split('|')
    # Level is 1-based, but we access 0-based array
    # Level 1 = species (rightmost), Level 7 = domain (leftmost)
    idx = len(parts) - level
    if idx >= 0 and idx < len(parts):
        return parts[idx]
    return None


def is_taxonomy_consistent(cluster, taxonomy_dict, tax_level):
    """Check if all sequences in a cluster with taxonomy information have consistent taxonomy at the specified level."""
    if not taxonomy_dict or tax_level <= 0:
        return True

    tax_values = set()
    for seq_id in cluster:
        if seq_id in taxonomy_dict:
            tax_value = get_taxonomy_level(taxonomy_dict[seq_id], tax_level)
            if tax_value:
                tax_values.add(tax_value)

    # If we have more than one unique taxonomy value at this level, the cluster is not consistent
    return len(tax_values) <= 1


def process_mcl_add_singletons(cluster_file, alignment_file, output_file, taxonomy_dict=None, tax_level=0):
    """Process MCL output, add singletons, and ensure taxonomy consistency if requested."""
    def read_cluster_file(file_path):
        clusters = []
        with open(file_path, 'r') as file:
            for line in file:
                clusters.append(line.strip().split('\t'))
        return clusters

    def read_alignment_file(file_path):
        alignment = []
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    genome_id = line[1:].split()[0]  # Remove ">" and get first field
                    alignment.append(genome_id)
        return alignment

    clusters = read_cluster_file(cluster_file)
    alignment = read_alignment_file(alignment_file)

    # If taxonomy information is provided, split clusters that are not consistent
    if taxonomy_dict and tax_level > 0:
        print(f"Checking taxonomy consistency at level {tax_level}...")
        consistent_clusters = []
        for cluster in clusters:
            if len(cluster) <= 1 or is_taxonomy_consistent(cluster, taxonomy_dict, tax_level):
                consistent_clusters.append(cluster)
            else:
                # Group sequences by taxonomy
                tax_groups = {}
                no_tax_seqs = []
                for seq_id in cluster:
                    if seq_id in taxonomy_dict:
                        tax_value = get_taxonomy_level(taxonomy_dict[seq_id], tax_level)
                        if tax_value:
                            if tax_value not in tax_groups:
                                tax_groups[tax_value] = []
                            tax_groups[tax_value].append(seq_id)
                        else:
                            no_tax_seqs.append(seq_id)
                    else:
                        no_tax_seqs.append(seq_id)

                # Add each taxonomy group as a separate cluster
                for tax_value, seqs in tax_groups.items():
                    consistent_clusters.append(seqs)

                # Add sequences without taxonomy as singletons
                for seq_id in no_tax_seqs:
                    consistent_clusters.append([seq_id])

        clusters = consistent_clusters

    # Add singletons
    flat_clusters = [item for sublist in clusters for item in sublist]
    singletons = [genome for genome in alignment if genome not in flat_clusters]

    for singleton in singletons:
        clusters.append([singleton])

    with open(output_file, 'w') as file:
        for cluster in clusters:
            file.write('\t'.join(cluster) + '\n')


def process_clusters_sortby_counts(count_file, id_file, output_file):
    """Sort clusters by sequence abundance based on count file."""
    if count_file and os.path.exists(count_file):
        try:
            df = pd.read_csv(count_file, sep='\t')
            # Check if the first column is unnamed
            if df.columns[0] == 'genome_id' or df.columns[0].lower() == 'genome_id':
                id_col = df.columns[0]
            else:
                # Rename the first column to genome_id
                df = df.rename(columns={df.columns[0]: 'genome_id'})
                id_col = 'genome_id'

            # Calculate the count (number of samples where the genome is present)
            df['count'] = df.drop(id_col, axis=1).apply(lambda row: sum(row > 0), axis=1)
            df = df[[id_col, 'count']]
            df.columns = ['genome_id', 'count']
        except Exception as e:
            # Silently create a default count file
            with open(id_file, 'r') as f:
                all_ids = set()
                for line in f:
                    all_ids.update(line.strip().split('\t'))

            df = pd.DataFrame({'genome_id': list(all_ids), 'count': [1] * len(all_ids)})
    else:
        # No count file provided, create a default one
        with open(id_file, 'r') as f:
            all_ids = set()
            for line in f:
                all_ids.update(line.strip().split('\t'))

        df = pd.DataFrame({'genome_id': list(all_ids), 'count': [1] * len(all_ids)})

    with open(id_file, 'r') as f:
        genome_ids = [line.strip().split('\t') for line in f.readlines()]

    reordered_ids = []
    for genome_id_list in genome_ids:
        sorted_id_list = sorted(genome_id_list, key=lambda x: df.loc[df['genome_id'] == x, 'count'].values[0] if x in df['genome_id'].tolist() and len(df.loc[df['genome_id'] == x, 'count'].values) > 0 else 0, reverse=True)
        reordered_ids.append(sorted_id_list)

    with open(output_file, 'w') as f:
        for reordered_id_list in reordered_ids:
            f.write('\t'.join(reordered_id_list) + '\n')


def generate_cluster_stats(cluster_dir, output_basename):
    """Generate statistics about the clustering results.

    Returns:
        pd.DataFrame: DataFrame containing the statistics, or None if no clusters were found.
    """
    stats = {'cutoff': [], 'num_clusters': [], 'num_singletons': [], 'avg_size': [], 'largest_cluster': []}

    found_clusters = False
    for cluster_file in os.listdir(cluster_dir):
        if cluster_file.endswith('.txt') and not cluster_file.startswith('temp'):
            # Extracting the cutoff value from the filename
            match = re.search(r'(\d+\.\d+)\.txt$', cluster_file)
            if match:
                found_clusters = True
                cutoff = float(match.group(1))
                with open(os.path.join(cluster_dir, cluster_file), 'r') as f:
                    clusters = [line.strip().split('\t') for line in f if line.strip()]
                    num_clusters = len(clusters)
                    num_singletons = sum(len(cluster) == 1 for cluster in clusters)
                    avg_size = np.mean([len(cluster) for cluster in clusters])
                    largest_cluster = max(len(cluster) for cluster in clusters)

                stats['cutoff'].append(cutoff)
                stats['num_clusters'].append(num_clusters)
                stats['num_singletons'].append(num_singletons)
                stats['avg_size'].append(avg_size)
                stats['largest_cluster'].append(largest_cluster)

    if not found_clusters:
        print(f"Warning: No cluster files found in {cluster_dir}")
        return None

    stats_df = pd.DataFrame(stats)
    stats_df.sort_values(by='cutoff', inplace=True)

    # Save individual stats file
    stats_df.to_csv(f'{output_basename}.clusterstats', sep='\t', index=False)

    return stats_df


def find_input_files(input_dir):
    """Find tree, alignment, and taxonomy files in the input directory."""
    tree_file = None
    alignment_file = None
    taxonomy_file = None

    # Find tree file
    tree_patterns = ['*.tree', '*.nwk', '*.tre', '*.treefile', '*.contree']
    for pattern in tree_patterns:
        tree_files = glob.glob(os.path.join(input_dir, pattern))
        if tree_files:
            tree_file = tree_files[0]
            break

    # Find alignment file
    aln_patterns = ['*.aln', '*.mafft', '*.mafft_t', '*.faa', '*.fa', '*.fasta']
    for pattern in aln_patterns:
        aln_files = glob.glob(os.path.join(input_dir, pattern))
        if aln_files:
            alignment_file = aln_files[0]
            break

    # Find taxonomy file
    tax_patterns = ['*.tsv', '*.tax', '*.taxonomy']
    for pattern in tax_patterns:
        tax_files = glob.glob(os.path.join(input_dir, pattern))
        if tax_files:
            # Check if it's a taxonomy file by looking at the content
            for tax_file in tax_files:
                try:
                    with open(tax_file, 'r') as f:
                        first_line = f.readline().strip()
                        if '|' in first_line and '\t' in first_line:
                            taxonomy_file = tax_file
                            break
                except:
                    continue
            if taxonomy_file:
                break

    return tree_file, alignment_file, taxonomy_file


def main():
    parser = argparse.ArgumentParser(description='PhyloDM Clustering Pipeline')
    parser.add_argument('--input_dir', type=str, required=True, help='Directory containing input files')
    parser.add_argument('--output_dir', type=str, help='Output directory (default: input_dir_results)')
    parser.add_argument('--tree', type=str, help='Tree file (if not specified, will look for *.tree, *.nwk, *.tre, *.treefile in input_dir)')
    parser.add_argument('--alignment', type=str, help='Alignment file (if not specified, will look for *.aln, *.mafft, *.mafft_t, *.faa in input_dir)')
    parser.add_argument('--taxonomy', type=str, help='Taxonomy file (if not specified, will look for *.tsv in input_dir)')
    parser.add_argument('--count', type=str, help='Count file (if not specified, will create a default count file)')
    parser.add_argument('--cutoffs', type=str, default='0.99,0.95,0.9,0.8,0.7', help='Comma-separated list of cutoff values (default: 0.99,0.95,0.9,0.8,0.7)')
    parser.add_argument('--tax_level', type=int, default=2, help='Taxonomy level for consistency check (1=species, 2=genus, 3=family, etc., 0=disable, default: 2)')

    args = parser.parse_args()

    # Set up input and output directories with absolute paths for internal use
    input_dir = os.path.abspath(args.input_dir)
    if args.output_dir:
        output_dir = os.path.abspath(args.output_dir)
    else:
        # Use relative path for default output directory name
        input_basename = os.path.basename(input_dir.rstrip('/'))
        output_dir = os.path.abspath(f"{input_basename}_results")

    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")

    create_dir(output_dir)

    # Find input files if not specified
    tree_file = args.tree
    alignment_file = args.alignment
    taxonomy_file = args.taxonomy
    count_file = args.count

    if not tree_file or not alignment_file:
        auto_tree, auto_aln, auto_tax = find_input_files(input_dir)
        if not tree_file and auto_tree:
            tree_file = auto_tree
        if not alignment_file and auto_aln:
            alignment_file = auto_aln
        if not taxonomy_file and auto_tax:
            taxonomy_file = auto_tax

    # Convert to absolute paths for internal use
    if tree_file and not os.path.isabs(tree_file):
        tree_file = os.path.abspath(tree_file)
    if alignment_file and not os.path.isabs(alignment_file):
        alignment_file = os.path.abspath(alignment_file)
    if taxonomy_file and not os.path.isabs(taxonomy_file):
        taxonomy_file = os.path.abspath(taxonomy_file)
    if count_file and not os.path.isabs(count_file):
        count_file = os.path.abspath(count_file)

    # Validate input files
    if not tree_file or not os.path.exists(tree_file):
        print(f"Error: Tree file {tree_file} not found")
        sys.exit(1)

    if not alignment_file or not os.path.exists(alignment_file):
        print(f"Error: Alignment file {alignment_file} not found")
        sys.exit(1)

    print(f"Tree file: {tree_file}")
    print(f"Alignment file: {alignment_file}")
    if taxonomy_file:
        print(f"Taxonomy file: {taxonomy_file}")
    if count_file:
        print(f"Count file: {count_file}")

    # Parse cutoffs
    try:
        cutoffs = [float(c) for c in args.cutoffs.split(',')]
    except:
        print(f"Error parsing cutoffs: {args.cutoffs}")
        print("Using default cutoffs: 0.99, 0.95, 0.9, 0.8, 0.7")
        cutoffs = [0.99, 0.95, 0.9, 0.8, 0.7]

    # Read taxonomy file if provided
    taxonomy_dict = {}
    if taxonomy_file and os.path.exists(taxonomy_file):
        taxonomy_dict = read_taxonomy_file(taxonomy_file)

    # Run the pipeline for each cutoff
    for cutoff in cutoffs:
        output_basename = os.path.join(output_dir, f"threshold_{cutoff:.2f}")

        print(f"Processing cutoff {cutoff}...")

        # Step 1: Run PhyloDM
        run_phylodm_parallel(tree_file, output_basename, [cutoff])

        # Step 2 and 3: Process the generated .mcl file
        clusters_dir = f"{output_basename}_clusters"
        create_dir(clusters_dir)
        create_dir(os.path.join(clusters_dir, "temp"))

        # Use os.path.join for proper path handling
        phylodm_out_dir = f"{output_basename}_phylodm_out"
        cluster_file = os.path.join(phylodm_out_dir, f"{os.path.basename(output_basename)}_{cutoff:.2f}.3c.mcl")
        temp_output_file = os.path.join(clusters_dir, "temp", f"temp_output_clusters_{cutoff:.2f}.txt")
        final_output_file = os.path.join(clusters_dir, f"{os.path.basename(output_basename)}_{cutoff:.2f}.txt")

        # Process MCL and add singletons
        process_mcl_add_singletons(cluster_file, alignment_file, temp_output_file, taxonomy_dict, args.tax_level)

        # Process clusters and sort by counts
        process_clusters_sortby_counts(count_file, temp_output_file, final_output_file)

    # Create directories for the final output
    clusters_dir = os.path.join(output_dir, "clusters")
    phylodm_out_dir = os.path.join(output_dir, "phylodm_out")
    create_dir(clusters_dir)
    create_dir(phylodm_out_dir)

    # Copy the final cluster files to the clusters directory
    for cutoff in cutoffs:
        # Copy cluster files
        src_cluster_dir = os.path.join(output_dir, f"threshold_{cutoff:.2f}_clusters")
        if os.path.exists(src_cluster_dir):
            # Copy the main cluster file
            src_file = os.path.join(src_cluster_dir, f"threshold_{cutoff:.2f}_{cutoff:.2f}.txt")
            dst_file = os.path.join(clusters_dir, f"threshold_{cutoff:.2f}.txt")
            if os.path.exists(src_file):
                with open(src_file, 'r') as src, open(dst_file, 'w') as dst:
                    dst.write(src.read())

            # Copy any other files in the cluster directory (except temp directory)
            for file in os.listdir(src_cluster_dir):
                if file != "temp" and os.path.isfile(os.path.join(src_cluster_dir, file)):
                    src_file = os.path.join(src_cluster_dir, file)
                    dst_file = os.path.join(clusters_dir, file)
                    with open(src_file, 'r') as src, open(dst_file, 'w') as dst:
                        dst.write(src.read())

        # Copy phylodm output files
        src_phylodm_dir = os.path.join(output_dir, f"threshold_{cutoff:.2f}_phylodm_out")
        if os.path.exists(src_phylodm_dir):
            for file in os.listdir(src_phylodm_dir):
                if os.path.isfile(os.path.join(src_phylodm_dir, file)):
                    src_file = os.path.join(src_phylodm_dir, file)
                    dst_file = os.path.join(phylodm_out_dir, file)
                    with open(src_file, 'r') as src, open(dst_file, 'w') as dst:
                        dst.write(src.read())

    # Generate cluster statistics for each cutoff
    all_stats = []
    for cutoff in cutoffs:
        cluster_dir = os.path.join(output_dir, f"threshold_{cutoff:.2f}_clusters")
        output_basename = os.path.join(output_dir, f"threshold_{cutoff:.2f}")
        stats = generate_cluster_stats(cluster_dir, output_basename)
        if stats is not None:
            all_stats.append(stats)

    # Create combined stats file
    if all_stats:
        combined_stats = pd.concat(all_stats)
        combined_stats.sort_values(by='cutoff', ascending=False, inplace=True)
        combined_stats.to_csv(os.path.join(output_dir, 'combined_stats.tsv'), sep='\t', index=False)

    # Clean up temporary directories
    import shutil
    for cutoff in cutoffs:
        # Remove threshold_*_clusters directories
        cluster_dir = os.path.join(output_dir, f"threshold_{cutoff:.2f}_clusters")
        if os.path.exists(cluster_dir):
            shutil.rmtree(cluster_dir)

        # Remove threshold_*_phylodm_out directories
        phylodm_dir = os.path.join(output_dir, f"threshold_{cutoff:.2f}_phylodm_out")
        if os.path.exists(phylodm_dir):
            shutil.rmtree(phylodm_dir)

        # Remove threshold_*.clusterstats files
        stats_file = os.path.join(output_dir, f"threshold_{cutoff:.2f}.clusterstats")
        if os.path.exists(stats_file):
            os.remove(stats_file)

    print(f"Clustering complete. Results are in {output_dir}/")
    print(f"See {output_dir}/README.md for more information.")

    # Create a README.md file
    create_readme(output_dir, tree_file, alignment_file, taxonomy_file, cutoffs, args.tax_level)


def create_readme(output_dir, tree_file, alignment_file, taxonomy_file, cutoffs, tax_level):
    """Create a README.md file in the output directory."""
    stats_file = os.path.join(output_dir, "combined_stats.tsv")
    if not os.path.exists(stats_file):
        print(f"Warning: Stats file {stats_file} not found")
        return

    stats_df = pd.read_csv(stats_file, sep='\t')

    # Get relative paths for display in README
    rel_tree_file = os.path.basename(tree_file)
    rel_alignment_file = os.path.basename(alignment_file)
    rel_taxonomy_file = os.path.basename(taxonomy_file) if taxonomy_file else None

    with open(os.path.join(output_dir, "README.md"), 'w') as f:
        f.write("# Phylogenetic Distance Matrix Clustering Results\n\n")
        f.write(f"This directory contains the results of clustering based on phylogenetic distances.\n\n")

        f.write("## Directory Structure\n\n")
        f.write("- `clusters/`: Contains the final clusters with sequences ordered by abundance for each threshold\n")
        f.write("- `phylodm_out/`: Contains the raw output from the PhyloDM clustering\n")
        f.write("- `combined_stats.tsv`: Contains statistics about the clustering for each threshold\n\n")

        f.write("## Input Files\n\n")
        f.write(f"- Tree file: `{rel_tree_file}`\n")
        f.write(f"- Alignment file: `{rel_alignment_file}`\n")
        if taxonomy_file:
            f.write(f"- Taxonomy file: `{rel_taxonomy_file}`\n")
            f.write(f"- Taxonomy level: {tax_level} ")
            if tax_level == 1:
                f.write("(species)\n")
            elif tax_level == 2:
                f.write("(genus)\n")
            elif tax_level == 3:
                f.write("(family)\n")
            elif tax_level == 4:
                f.write("(order)\n")
            elif tax_level == 5:
                f.write("(class)\n")
            elif tax_level == 6:
                f.write("(phylum)\n")
            elif tax_level == 7:
                f.write("(domain)\n")
            else:
                f.write("\n")
        f.write("\n")

        f.write("## Clustering Thresholds\n\n")
        f.write("The clustering was performed at the following thresholds:\n")
        for cutoff in cutoffs:
            if cutoff >= 0.99:
                f.write(f"- **{cutoff:.2f}**: Very stringent, only very closely related sequences are clustered\n")
            elif cutoff >= 0.95:
                f.write(f"- **{cutoff:.2f}**: Stringent clustering\n")
            elif cutoff >= 0.9:
                f.write(f"- **{cutoff:.2f}**: Moderate clustering\n")
            elif cutoff >= 0.8:
                f.write(f"- **{cutoff:.2f}**: Relaxed clustering\n")
            else:
                f.write(f"- **{cutoff:.2f}**: Very relaxed, distantly related sequences may be clustered\n")
        f.write("\n")

        f.write("## Clustering Statistics\n\n")
        f.write("| Threshold | Number of Clusters | Number of Singletons | Average Cluster Size | Largest Cluster Size |\n")
        f.write("|:---------:|:-----------------:|:-------------------:|:-------------------:|:-------------------:|\n")

        for _, row in stats_df.iterrows():
            f.write(f"| **{row['cutoff']:.2f}** | {row['num_clusters']} | {row['num_singletons']} | {row['avg_size']:.2f} | {row['largest_cluster']} |\n")

        f.write("\n")

        f.write("## Interpretation\n\n")
        f.write("As the threshold decreases from high to low values:\n")

        if len(stats_df) >= 2:
            highest_cutoff = stats_df.iloc[0]
            lowest_cutoff = stats_df.iloc[-1]

            f.write(f"- The number of clusters decreases (from {highest_cutoff['num_clusters']} to {lowest_cutoff['num_clusters']})\n")
            f.write(f"- The number of singletons decreases (from {highest_cutoff['num_singletons']} to {lowest_cutoff['num_singletons']})\n")
            f.write(f"- The average cluster size increases (from {highest_cutoff['avg_size']:.2f} to {lowest_cutoff['avg_size']:.2f})\n")
            f.write(f"- The largest cluster size increases (from {highest_cutoff['largest_cluster']} to {lowest_cutoff['largest_cluster']})\n")

        f.write("\nThis pattern is expected because:\n")
        f.write("- Lower thresholds (e.g., 0.7) result in fewer, larger clusters as more distantly related sequences are grouped together\n")
        f.write("- Higher thresholds (e.g., 0.99) result in more, smaller clusters as only very closely related sequences are grouped together\n\n")

        f.write("## Applications\n\n")
        f.write("The clustering results can be used to:\n")
        f.write("1. **Identify groups of closely related sequences**: Each cluster represents a group of sequences with similar evolutionary history\n")
        f.write("2. **Reduce redundancy in the dataset**: Select one representative sequence from each cluster\n")
        f.write("3. **Study the diversity and distribution of sequences**: Analyze the distribution of sequences across different clusters\n")
        f.write("4. **Select representative sequences for further analysis**: Choose sequences that represent the diversity of the dataset\n")

        if taxonomy_file and tax_level > 0:
            f.write("\n## Taxonomy-Based Clustering\n\n")
            f.write(f"The clustering was performed with taxonomy consistency checks at level {tax_level} ")
            if tax_level == 1:
                f.write("(species).\n")
            elif tax_level == 2:
                f.write("(genus).\n")
            elif tax_level == 3:
                f.write("(family).\n")
            elif tax_level == 4:
                f.write("(order).\n")
            elif tax_level == 5:
                f.write("(class).\n")
            elif tax_level == 6:
                f.write("(phylum).\n")
            elif tax_level == 7:
                f.write("(domain).\n")
            else:
                f.write(".\n")

            f.write("This means that clusters were split if they contained sequences with different taxonomy at this level.\n")
            f.write("Sequences without taxonomy information were allowed to cluster with any other sequences.\n")


if __name__ == "__main__":
    main()
