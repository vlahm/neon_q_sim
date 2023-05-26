
welcome to the dungeon. if you'd like to reproduce the LSTM results reported in the paper, you either need access to
an HPC cluster, or a desktop with a GPU and a couple of weeks for runtime

1. load a CUDA module (HPC) or install CUDA (local): https://developer.nvidia.com/cuda-downloads
     an example slurm batch script
1. create a python environment from environment.yml:
   conda env create -f environment.yml
1. run 



mention:
    need about 30GB ram to run lstms as-is
    need at least 62 to run full nhm
    also talk about gpu specs
        put any of this in paper?
    run_lstms.py was set up orig to do a lot more. sry 
    unzip lstm_runs.zip
    unzip input_data.zip
    need full camels dataset AND daymet isolate
    added custom eval param (pbias)
    install rgee
        setup config
    if you run everything from scratch, you'll have to update all the base_run_dirs for finetuning


things to hardcode:
    strategy courier output
    IS_PGDL
    ...
