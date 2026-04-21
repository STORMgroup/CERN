# Testing CERN

 This directory contains instructions to install the tools needed for processing datasets and generating results, as well as the following sub-directories:

 - [./data](./data/) contains the datasets used in the CERN paper and their full versions, as well as scripts to fully process datasets for read-mapping.
 - [./evaluation](./evaluation) contains the scripts to run Rawhash2 on each of the datasets used in the paper, and summarize the results.
 - [./scripts](./scripts) contains useful scripts for running tools and processing files.


# Installing necessary tools

Follow these steps to install the tools necessary to evaluate CERN.
To compile these tools, will require gcc version 11.2.0 or higher.

```bash
# Create the bin:
mkdir bin

#Re-adding the binary to your path is necessary after you start a new shell session.

# Step 0 Creating a conda environment:

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda create -n cern-env python=3.14
conda activate cern-env

# Step 1: Compiling tools

# Cloning and compiling RawHash2 v2.1
git clone --recursive https://github.com/STORMgroup/RawHash2.git rawhash2 \
  && cd rawhash2 \
  && git checkout a868b1e37b683d2c5b3f827d38a3931000647a7e \
  && git submodule update --init --recursive \
  && make cmake \
  && cp ./bin/rawhash2 ../bin/ && cd ..

# Downloading and compiling minimap2 v2.24
wget https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24.tar.bz2; tar -xf minimap2-2.24.tar.bz2; rm minimap2-2.24.tar.bz2; mv minimap2-2.24 minimap2; cd minimap2 && make && cp minimap2 ../bin/ && cd ..

# Downloading dorado v1.4.0
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.4.0-linux-x64.tar.gz; tar -xf dorado-1.4.0-linux-x64.tar.gz; rm dorado-1.4.0-linux-x64.tar.gz

# Downloading dorado v0.9.0 
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.0-linux-x64.tar.gz; tar -xf dorado-0.9.0-linux-x64.tar.gz; rm dorado-0.9.0-linux-x64.tar.gz

# Downloading Campolina
git clone https://github.com/lbcb-sci/Campolina; cd Campolina; mkdir weights; cd weights; wget https://zenodo.org/records/15626806/files/R10_model.pth; cd ../../

# Add binaries to PATH
export PATH=$PWD/bin:$PATH
```