# Installing necessary tools

Follow these steps to install the tools necessary to evaluate CERN.
To compile these tools, use gcc version 11.2.0 or higher.

```
#We will use cern-env directory to download and compile the tools from their repositories
mkdir -p cern-env/bin && cd cern-env

#Important: If you already completed the Step 0 and Step 1 as described, you can skip these steps and add the binaries to your PATH again
#Re-adding the binary to your path is necessary after you start a new shell session.

#Optional Step 0 Creating a conda environment (Note we highly recommend using conda for easy installation of dependencies).
#If not using conda, the packages with the specified versions below (e.g.,  python=3.6.15) must be installed manually in your environment
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda create -n cern-env python=3.6.15 pip=21.3.1 ont_vbz_hdf_plugin=1.0.1 seqkit=2.5.1
conda activate cern-env

#Installing SRA Toolkit
wget -qO- https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.1/sratoolkit.3.0.1-ubuntu64.tar.gz | tar xzv; cp -r ./sratoolkit.3.0.1-ubuntu64/bin/* bin/; rm -rf sratoolkit.3.0.1-ubuntu64

#Step 1 Compiling the tools
#Cloning and compiling RawHash2
#Recommended: using CMake (by default compiles with POD5 support only)
git clone --recursive https://github.com/STORMgroup/RawHash2.git rawhash2 && cd rawhash2 && make cmake && cp ./bin/rawhash2 ../bin/ && cd ..

#Alternative: using Make only (no CMake required)
# git clone --recursive https://github.com/STORMgroup/RawHash2.git rawhash2 && cd rawhash2 && make && cp ./bin/rawhash2 ../bin/ && cd ..


#Cloning and compiling UNCALLED v2.1 commit 58dbec69f625e0343739d821788d536b578ea41d
git clone --recursive https://github.com/skovaka/UNCALLED.git uncalled && cd uncalled && git checkout 58dbec69f625e0343739d821788d536b578ea41d && cd submods/bwa && git pull origin master && cd ../../ && pip3 install . && cd ..

#Downloading and compiling minimap2 v2.24
wget https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24.tar.bz2; tar -xf minimap2-2.24.tar.bz2; rm minimap2-2.24.tar.bz2; mv minimap2-2.24 minimap2; cd minimap2 && make && cp minimap2 ../bin/ && cd ..

#Downloading dorado v1.4.0
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.4.0-linux-x64.tar.gz; tar -xf dorado-1.4.0-linux-x64.tar.gz; rm dorado-1.4.0-linux-x64.tar.gz;

#Downloading dorado v0.9.0
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.0-linux-x64.tar.gz; tar -xf dorado-0.9.0-linux-x64.tar.gz; rm dorado-0.9.0-linux-x64.tar.gz;

#Downloading Campolina
git clone https://github.com/lbcb-sci/Campolina;

#Step 2 Adding binaries to PATH
#If you are skipping Step 0 and Step 1, uncomment the following line and execute:
# conda activate cern-env
export PATH=$PWD/bin:$PATH
```