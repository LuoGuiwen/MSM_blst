#!/bin/bash

# Parse the command-line arguments
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        benchmark=*)
        benchmark="${key#*=}"
        shift
        ;;
        config=*)
        configs="${key#*=}"
        shift
        ;;
        *)
        shift
        ;;
    esac
done

if [ -z "$benchmark" ]
then
    echo "Please enter a benchmark type (1 or 2)".
    exit 1
fi

# Set default name as "config" if not provided
if [ -z "$configs" ]
then
    if [ "$benchmark" -eq 1 ]; then
        configs=10
    elif [ "$benchmark" -eq 2 ]; then
        configs=10
    fi
fi

# Print the benchmark
echo "The benchmark is p$benchmark"

# Print each config
IFS=',' read -ra configs <<< "$configs"
for config in "${configs[@]}"
do
    echo "Execute with config_file_n_exp_$config.h"
    make benchmark=$benchmark config=$config
    ./main_test_p$benchmark
done