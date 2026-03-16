# Read Mapping Evaluation

In order to evaluate CERN, please follow these steps:

1. Download the dataset in /data/
2. Generate the ground truth, by basecalling in /data/, then running the minimap2 script in the dataset's evaluation directory
3. Run RawHash2 with and without the corrected events, using the script in the dataset's evaluation directory
4. Compare the results using the generate_results script in the dataset's evaluation directory

It is recommended to create a small subset of pod5 files for large datasets.