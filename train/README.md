# Training HMMs

This directory contains the needed python scripts for training a nanopore event modeling HMM.

The training configuration can be edited at the top of each python file.

# Quickstart

First create a conda environment with the necessary python packages to train HMMs:

```bash
conda env create -n CERN_train -f environment.yml
```

Without editing any of the python files, the following command will create and train a 196-state HMM on Scrappie (R9.4)-segmented E. coli nanopore reads:

```bash
bash train_pipeline.sh
```

Note that the training parameters are in the python files, and editing those will change the output of the above command.

## Synthic Data Training

There are two steps CERN uses when training an HMM on synethic data:

1. Training on random sequences of DNA (in [./HMM_train.py](./HMM_train.py))
   1. This can be began from a provided starting HMM using [./HMM_continue_training.py](./HMM_train_continue_training.py)
2. Training on a debrujin sequence of DNA (in [./HMM_train_deBrujin.py](./HMM_train_deBrujin.py))

Each script will print the HMM to an output file specified in the training configuation.

## Experimental Data Training

In [./HMM_train_errors.py](./HMM_train_errors.py), the HMM is trained on real, segmented nanopore signals.

Note that in order to run this, the necessary event files must already be produced. (These can be created in [../test/data/cern_datasets](../test/data/cern_datasets))

## Pre-trained models

We provide already trained models in [./models](./models/) for the Scrappie (R9.4), Scrappie (R10.4.1), and Campolina segmentation algorithms.