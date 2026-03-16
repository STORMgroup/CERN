# CERN
Correcting Errors in Raw Nanopore Signals

## Features

### Training HMMs

In ./train, you can train Hidden Markov Models to model nanopore dynamics. The models created here will be used to correct nanopore signals.

### Segmenting and refining nanopore signals

./src contains the c++ code used to correct nanopore signals using HMMs, here is an example:

```
run_cern ./path/to/hmm_file ./path/to/event_file > corrected_events.txt
```

### Testing error correction performance

The general workflow to measure performace involves multiple steps:

1. Generate events, using the executables in ./src/segmentation/ for RawHash T-test based segmentation. These should be run the same way RawHash1/2 are run. You can also create events yourself using other tools, such as Campolina.
2. Correct events using the run_cern program in ./src/refinement/ Already trained models are available in ./src/train/models.
3. Run RawHash2 with and without the corrected events. RawHash2 has support for custom events built-in, using the `--events-file` flag.


## Requirements

In order to train the HMMs, you will need to set up your python environment using ./train/environment.yml

For evaluation and testing, you will need the following programs installed and in your $PATH:
 - RawHash2 (specifically, the version at https://github.com/STORMgroup/RawHash2) that has support for custom events.
 - minimap2
 - UNCALLED4
 - dorado

For instructions to install each of the above programs, check https://github.com/STORMgroup/RawHash2/tree/main/test

minimap2 is used to create a ground-truth mapping, UNCALLED4 is used to annotate and compare mappings, RawHash2 is used to evaluate the quality of corrected events, and dorado is used to basecall datasets if necessary.
