welcome to the dungeon. here we've included two absolute monstrosities of R scripts that were used
to generate the configuration files required to specify NeuralHydrology models. these scripts were
"developed" over many moons, and in the most spaghettified fashion imaginable. they are not even spaghetti
scripts. they are capellini.

the good news is, you needn't pay them any mind. their resultant file structures and run_lstms.py are
all that's needed to reproduce our LSTM predictions. they are here only for completeness.

that's not to say run_lstms.py is any more elegant. it too received countless patches that were stitched
on as if by Captain Hook. Luckily it's only a few hundred lines, rather than 8000+. 

anyway, if you want to reproduce the LSTM results reported in the paper, here are the steps:

1. install CUDA on your machine. this process can be quite straightforward, or it can be an unforgettable nightmare. good luck.
   https://developer.nvidia.com/cuda-downloads
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
