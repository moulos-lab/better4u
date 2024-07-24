#!/bin/bash

# List of allowed executables
ALLOWED_COMMANDS=("bcftools" "vcftools" "PRSice.R" "plink" "regenie" "gcta64" "snptest")

# Function to print the list of allowed commands
print_allowed_commands() {
    echo "Allowed commands are:"
    for cmd in "${ALLOWED_COMMANDS[@]}"; do
        echo " - $cmd"
    done
}

# Check if the command is in the list of allowed commands
is_command_allowed() {
    local cmd=$1
    for allowed_cmd in "${ALLOWED_COMMANDS[@]}"; do
        if [ "$cmd" == "$allowed_cmd" ]; then
            return 0
        fi
    done
    return 1
}

# Check if at least one argument is provided
if [ $# -eq 0 ]; then
    echo "No command provided. Here is a list of allowed commands:"
    print_allowed_commands
    exit 1
fi

# Check if the provided command is in the allowed list
if is_command_allowed "$1"; then
    exec "$@"
else
    echo "Error: Command '$1' is not allowed."
    echo "Here is a list of allowed commands:"
    print_allowed_commands
    exit 1
fi