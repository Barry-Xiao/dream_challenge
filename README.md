# dream_challenge
codes for dream challenge

## Install anaconda or miniconda

https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

https://docs.anaconda.com/anaconda/install/index.html


## Install jupyter notebook in base environment

Open your conda command prompt and install jupyter notebook in base environments 

```
conda install -c conda-forge notebook
conda install -c conda-forge nb_conda_kernels
```

## Create environments using the yml file in the repo

```
#repo-folder is the local path of this github repo

cd repo-folder

#create environment using dl.yml

conda env create -f dl.yml
```

## switch to dream_dl env

```
conda env activate dream_dl
```

## Launch jupyter notebook from this environment and specifying `Python [conda env:dream_dl]` kernel

```
jupyter notebook
```

And create new notebook specifying `Python [conda env:dream_dl]` kernel, or switch to `Python [conda env:dream_dl]` kernel in a existing notebook using `Change kernel` in `Kernel` tab.
