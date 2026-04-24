# Training HMMs

This directory contains the code for training a nanopore event modeling HMM.

The training configuration can be edited at the top of each python file.


## Training

```bash
./train_hmm <n_states> <kmer_file> <k> <n_threads> > output.hmm
```


## Pre-trained models

We provide already trained models in [./models](./models/) for the Scrappie (R9.4), Scrappie (R10.4.1), and Campolina segmentation algorithms.

They are trained on a sequence which contains all 10-mers exactly once.