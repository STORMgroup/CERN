# Parameter Search

This directory contains the code for doing an automatic parameter sweep over the CRANE parameters P(stay) and P(skip)


## Use

This is an example of how to search over parameters, where hmm_file is the .tsv file of a trained hmm, and segmenter must either be "scrappieR9", "scrappieR10", or "campolina".

Before running this, ensure that you have installed the necessary tools in [../test/](../test/), and downloaded and prepared the training data in [../test/data](../test/data)

```bash
# bash parameter_search.sh <n_threads> <hmm_file> <segmenter> 
bash parameter_search.sh 32 ../models/hmm_16.tsv scrappieR9
```