# Training HMMs

This directory contains the code for training a nanopore event modeling HMM.


## Training

This is an example of how to train an HMM, where n_states is the number of states, kmer_file is the pore model file with expected signal values, and the training sequence length is controlled by k, it will contain every k-mer exactly once.

```bash
make
# ./train_hmm <n_states> <kmer_file> <k> <n_threads> > output.tsv
./train_hmm 16 ../extern/kmer_models/uncalled_r1041_model_only_means 11 1 > output.tsv
```


## Pre-trained models

We provide already trained models in [./models](./models/) of various sizes.

They were trained on a sequence which contains all 11-mers exactly once.