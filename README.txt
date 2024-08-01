diffGEK for erythropoiesis

Folder run_models: To create the diffGEK input for erythropoiesis, we first create a matrix for spliced and unspliced counts per condition for the genes considered with the code prepare_input.ipynb. Since the adatas are difficult to upload, we provide the matrixes directly in the folder data_input.

Then we run the create_input.m code in MATLAB. This creates a data structure per gene. Then we run all the models per each gene. The command run_models.sh shows how to do it in the background. Since this is a long computation, we provide all the results in the results folder.

After the models have finished, we compute the statistics, running the codes of the stats folder in the given order (final tables already provided).

Finally, we can create the paper figures with the codes in the corresponding folder.
