diffGEK for simulated data

Folder generate_data: To create the diffGEK input for the simulated data,  we run the generate_data_model_n.m codes in MATLAB. They create a data structure per gene and per ground truth model. 

Then we run all the models in the folder run_models. The command run_models.sh shows how to do it in the background.

After the models have finished, we compute the statistics, running the codes of the stats folder in the given order.

Finally, we can create the paper figures with the codes in the corresponding folder.
