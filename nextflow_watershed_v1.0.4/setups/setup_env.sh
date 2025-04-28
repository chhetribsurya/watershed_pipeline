#!/bin/bash

# Function to check if a command exists
command_exists () {
    command -v "$1" >/dev/null 2>&1
}

# Set the installation directory
INSTALL_DIR=$(pwd)/tools
MINICONDA_DIR=$INSTALL_DIR/miniconda

# Check if conda is installed, if not, install Miniconda in the specified directory
if ! command_exists conda || [[ ! -d "$MINICONDA_DIR" ]]; then
    echo "Conda is not installed or not found in the specified directory. Installing Miniconda..."
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh -b -p $MINICONDA_DIR
        rm Miniconda3-latest-Linux-x86_64.sh
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
        bash Miniconda3-latest-MacOSX-x86_64.sh -b -p $MINICONDA_DIR
        rm Miniconda3-latest-MacOSX-x86_64.sh
    else
        echo "Unsupported OS. Please install Miniconda manually."
        exit 1
    fi
    export PATH="$MINICONDA_DIR/bin:$PATH"
else
    export PATH="$MINICONDA_DIR/bin:$PATH"
fi

# Check if mamba is installed, if not, install it
if ! command_exists mamba; then
    echo "Mamba is not installed. Installing mamba..."
    conda install -c conda-forge mamba -y
fi

# Remove existing environment if it exists
#if conda info --envs | grep -q 'watershed_env'; then
#    echo "Removing existing environment 'watershed_env'..."
#    mamba env remove --name watershed_env
#fi

# Check if the watershed_env environment exists
if conda info --envs | grep -q 'watershed_env'; then
    echo "Updating existing environment 'watershed_env'..."
    mamba env update --name watershed_env --file ./env_ymls/watershed_env.yml --prune
else
    # Create the new environment
    echo "Creating new environment 'watershed_env'..."
    mamba env create -f ./env_ymls/watershed_env.yml
fi

# Create the new environment
echo "Creating new environment 'watershed_env'..."
mamba env create -f ./env_ymls/watershed_env.yml

# Activate the environment
echo "Activating environment 'watershed_env'..."
eval "$(conda shell.bash hook)"
mamba activate watershed_env

echo "Environment 'watershed_env' setup is complete."
