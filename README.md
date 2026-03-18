# CERN Overview
CERN is an HMM-based tool to correct raw nanopore signal events. To do this, it 1) Trains an HMM on synthetic and experimental nanopore event data and 2) uses the Viterbi algorithm and the learned HMM to correct novel signals.

In its current implementation, CERN takes in event sequences from a .tsv file, and outputs the corrected events to another .tsv.


# Installation

## Prequisites

 - C++ compiler
   - GCC 11+ (`g++` 11 or later)
 - GNU Make
 - Python 3.14.3+

## Quick Start

* Clone the code from its GitHub repository:

```bash
git clone https://github.com/STORMgroup/CERN.git
cd CERN
```

 * Make the executable

```bash
cd src
make
```

  * Then add it to your $PATH or bin:

```bash
cp run_cern /path/to/bin/
```

# Usage

Correct nanopore event sequences by providing the path to the trained HMM, as well as the .tsv containing event sequences:

```
run_cern ./path/to/hmm_file.txt ./path/to/event_file.tsv > corrected_events.tsv
```

# Training HMMs

In [train](./train/), you can train Hidden Markov Models to model nanopore dynamics. The models created here can be used to correct nanopore signals.

# Reproducing the results

Follow the steps in [test/evaluation/README.md](./test/evaluation/README.md).