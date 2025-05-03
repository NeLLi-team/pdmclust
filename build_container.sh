#!/bin/bash

# This script builds the Apptainer container for the PhyloDM Clustering Pipeline

# Check if apptainer is installed
if command -v apptainer &> /dev/null; then
    CONTAINER_CMD="apptainer"
elif command -v singularity &> /dev/null; then
    CONTAINER_CMD="singularity"
else
    echo "Error: Neither apptainer nor singularity is installed"
    exit 1
fi

echo "Using $CONTAINER_CMD to build the container"

# Initialize build success flag
BUILD_SUCCESS=false

# Try different build approaches
echo "Building container with --fakeroot..."
$CONTAINER_CMD build --fakeroot phylodm-clustering.sif phylodm-clustering.def
if [ $? -eq 0 ]; then
    BUILD_SUCCESS=true
else
    echo "Build with --fakeroot failed, trying with --fakeroot and --ignore-fakeroot-command..."
    $CONTAINER_CMD build --fakeroot --ignore-fakeroot-command phylodm-clustering.sif phylodm-clustering.def
    if [ $? -eq 0 ]; then
        BUILD_SUCCESS=true
    else
        echo "Build with --fakeroot and --ignore-fakeroot-command failed, trying without fakeroot..."
        $CONTAINER_CMD build phylodm-clustering.sif phylodm-clustering.def
        if [ $? -eq 0 ]; then
            BUILD_SUCCESS=true
        else
            echo "Standard build failed, trying sandbox build method..."
            # Create a sandbox directory
            mkdir -p ./sandbox
            # Build into the sandbox
            $CONTAINER_CMD build --sandbox ./sandbox phylodm-clustering.def
            if [ $? -eq 0 ]; then
                # Convert sandbox to SIF
                $CONTAINER_CMD build phylodm-clustering.sif ./sandbox
                if [ $? -eq 0 ]; then
                    BUILD_SUCCESS=true
                    # Clean up sandbox
                    rm -rf ./sandbox
                else
                    echo "Failed to convert sandbox to SIF"
                fi
            else
                echo "Sandbox build failed"
            fi
        fi
    fi
fi

if [ "$BUILD_SUCCESS" = true ]; then
    echo "Container built successfully: phylodm-clustering.sif"
    echo "To run the container:"
    echo ""
    echo "Basic Usage:"
    echo "$CONTAINER_CMD run --bind /path/to/input:/input --bind /path/to/output:/output phylodm-clustering.sif --input_dir /input --output_dir /output"
    echo ""
    echo "With Taxonomy-Based Clustering:"
    echo "$CONTAINER_CMD run --bind /path/to/input:/input --bind /path/to/output:/output phylodm-clustering.sif --input_dir /input --output_dir /output --tax_level 2"
    echo ""
    echo "With Custom Cutoffs:"
    echo "$CONTAINER_CMD run --bind /path/to/input:/input --bind /path/to/output:/output phylodm-clustering.sif --input_dir /input --output_dir /output --cutoffs 0.95,0.9,0.85"
    echo ""
    echo "Extract Cluster Representatives:"
    echo "$CONTAINER_CMD run --bind /path/to/input:/input --bind /path/to/output:/output phylodm-clustering.sif --input_dir /input --output_dir /output --extract 0.7"
    echo ""
    echo "For more information, run:"
    echo "$CONTAINER_CMD run-help phylodm-clustering.sif"
else
    echo "Error building container"
    exit 1
fi
