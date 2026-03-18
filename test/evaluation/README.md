# Error Correction Read Mapping Evaluation

In order to evaluate CERN, follow these steps:

1. Make the segmenter and CERN executable in [../../src/](../../src/)
2. Install the necessary tools by following the steps in [../tools/README.md](../tools/)
3. Download then prepare the datasets in [../data/cern_datasets/](../data/cern_datasets/)
4. Create the RawHash2 indices for each reference genome in this directory by running:
```bash
bash ./create_indexes.sh
```
5. Run RawHash2 with and without the corrected events, using the script map_all.sh in a dataset's evaluation directory
```bash
cd d1_ecoli_hpc_off # Navigate to the wanted directory
bash ./map_all.sh
```
6. Compare the results by running annotate_pafs.sh then summarize_annotations.sh in a dataset's evaluation directory
```bash
cd d1_ecoli_hpc_off # Navigate to the wanted directory
bash ./annotate_pafs.sh
bash ./summarize_annotations.sh
```