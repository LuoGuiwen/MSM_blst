#!/bin/bash

# Parse the command-line arguments
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        group=*)
        group="${key#*=}"
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

if [ -z "$group" ]
then
    group = 1
fi

# Set default name as "config" if not provided
if [ -z "$configs" ]
then
        configs=10
fi

# Print the group
echo "The group is p$group"

# Print each config
IFS=',' read -ra configs <<< "$configs"
for config in "${configs[@]}"
do
    echo "Execute with config_file_n_exp_$config.h"
    make group=$group config=$config
    ./main_test_p$group
done