# MIAO (Multispectral Imaging Analysis Optimizer)

## Getting started with `LoadData`, the core function of MIAO
There are two main outputs of `LoadData`: `data` and `m`. `data` is self-explanatory: it contains your data! `m` is a structure that contains all the metadata for the run you are analyzing. Using LoadData is really straightforward once you understand how to set up your workspace.

### Step 1: Set path and initialize `m`
First, run this command to tell MATLAB where the load code is:

`addpath(genpath('/local_mount/space/juno/1/Software/MIAO'))`

To initialize m with the default load options, simply input `m = makem;` This outputs m with all of the neccesary default options to load a dataset. NOTE: _This tutorial assumes that you are working on a ssh session within the LFOI server cluster._

### Step 2: Find the folder where your data is located
LoadData only takes two inputs: path, and m. You can either directly input `path`, or find it using the handy command `findmousefolder` if you are working on the LFOI server cluster. Simply run the command `path = findmousefolder(mouse, run, stim)` and it should spit out the path to the dataset of interest.

### Step 3: Load your data
now you are ready to load your data! Simply run the command `[m,data] = LoadData(path,m);`. You should see some initialization messages, a loading bar indicating the load progress, and finally a `Done` message indicating completion.

Here's all this code with an example run to get you started.

```
addpath(genpath('/local_mount/space/juno/1/Software/MIAO'))
m = makem;
path = findmousefolder('cm62_2','runF','2');
[m,data] = LoadData(path,m);
```
