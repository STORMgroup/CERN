# CERN Evaluation Data

Get started by downloading the CERN evaluation datasets:

```bash
bash download_cern_data.sh
```


Then prepare the datasets using the preparation scripts below.

 - Before running these, we recommend inspecting them to ensure the needed tools have been installed and/or compiled.
 - Ensure that you have created the virtual environment in [../../../src/segmentation](../../../src/segmentation/).
   - This is necessary for the segmentation step to work correctly.
 - Running these will require access to GPUs due to the the dorado basecalling steps.

```bash
bash prepare_d1.sh
bash prepare_d2.sh
bash prepare_d3.sh
```

Currently, Campolina has issues running on some machines, so the Campolina segmentation and correction steps have been commented out in the above scripts. Feel free to uncomment and experiment.

If you want to train your own HMM models using the provided configuration, you will need to prepare the training dataset as well:

```bash
bash prepare_training.sh
```