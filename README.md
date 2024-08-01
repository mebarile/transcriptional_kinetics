# diffGEK: a tool to compute transcriptional kinetics

This tool is meant to uncover differences in transcriptional kinetics among 2 conditions for any differentiating trajectory for which we have single cell resolution.

Since it is based on models such as scVelo (Bergen et al. 2020) and pseudodynamics (Fischer et al. 2019), it requires the same softwares as these tools: scanpy, scVelo, MATLAB, AMICI and PESTO. For our paper, we also used cellrank, tradeSeq and cellDanncer, but these latter are optional. 

We applied our tools to simulated data as well as erythropoiesis and myelopoiesis published datasets.

Before any analysis, you need to build the models with the corresponding codes in the model_building folder. Further instructions are provided in the respective folders.
