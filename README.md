# CERN
Correcting Errors in Raw Nanopore Signals

## Features

### Segmenting and refining nanopore signals

./src contains the c++ code used to correct nanopore signals using HMMs, here is an example:

```
run_cern ./path/to/hmm_file ./path/to/event_file > corrected_events.txt
```

### Training HMMs

In ./train, you can train Hidden Markov Models to model nanopore dynamics. The models created here will be used to correct nanopore signals.

### Testing error correction performance

To evaluate CERN's error correction accuracy and runtime, clone this repository and visit ./test/evaluation for instructions to get started.

## Requirements

The requirements for each section of CERN are listed in their respective repositories. Get started by installing the necessary tools for evaluation in ./test/tools