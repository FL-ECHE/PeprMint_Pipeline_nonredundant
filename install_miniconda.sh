

# Detect the operating system
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    # Linux
    curl -o miniconda.sh -L https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
elif [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS
    curl -o miniconda.sh -L https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
else
    echo "Unsupported operating system."
    exit 1
fi

bash miniconda.sh -b
rm miniconda.sh
eval "$(~/miniconda3/bin/conda shell.bash hook)"
conda init
conda install -c conda-forge mamba -y
