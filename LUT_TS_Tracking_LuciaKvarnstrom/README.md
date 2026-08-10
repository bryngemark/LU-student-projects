## Pattern-Based Electron Counting Algorithm for LDMX TS

This is the software developed for the thesis ["Pattern-Based Electron Counting Algorithm for LDMX"](https://lup.lub.lu.se/student-papers/search/publication/9232752). **NOTE: The files refered to as ClusterViewerAnalyzer and LUTAnalyzer in the paper have been renamed within *ldmx-sw* to ClusterTripletMaker and PatternLUTMaker, respectively.**

This software builds upon (and parts of it are included within) the GitHub repo [*ldmx-sw*](https://github.com/LDMX-Software/ldmx-sw). To get started with using *ldmx-sw*, please see the [tutorial](https://ldmx-software.github.io/using/analysis/ldmx-sw.html). You will need `denv` (see tutorial sec. 1) and if developing, see using `just` under sec. 11 of the tutorial. You will need to clone any version of *ldmx-sw* *after* v.4.8.1. I contributed my changes to *ldmx-sw* in [PR #2052](https://github.com/LDMX-Software/ldmx-sw/pull/2052/changes) and [PR #2070](https://github.com/LDMX-Software/ldmx-sw/pull/2070/changes).

Described below is the workflow to perform LUT-method TS tracking. The files included here under `/configs` are also in *ldmx-sw*.

## First, run simulations: 

You need a .ROOT file of events to work with. You can generate it using `runLUTeventgen.py` with the command:

```bash
$ denv fire TrigScint/exampleConfigs/runLUTeventgen.py
```

This will produce a file called `SimSamples.root`. I only used single-electron events and recommend `p.max_events = 10000`. `SimSamples.root` will then be used as input in the next step, which is to create a lookup table. 

## Create a LUT:

You will need to create a text file containing the clusters made in TS pads 1, 2, and 3. The command:

```bash
$ denv fire TrigScint/exampleConfigs/runLUTclusterstxt.py
```
will take the simulation samples file `SimSamples.root` as input, perform clustering and produce `Clusters.root`, and also create the output text file `clusters.txt`. In the file `runLUTclusterstxt.py`, you need to configure the boolean variable `truth_filtering` to decide if you would like to (before clustering is performed) use `TruthHitProducer` to isolate only those hits which occur from beam electrons. I performed these LUT steps twice in my thesis, once with `truth_filtering=True` and once with `truth_filtering=False` to compare. 

Once you have a `clusters.txt` file, you will use it as input in the next step, to create a lookup table text file from it. The command

```bash
$ denv fire TrigScint/exampleConfigs/runLUTana.py
```

will produce a file called `LUT.txt`. This is your lookup table. In the config `runLUTana.py`, you will need to configure the variable `lut_threshold`. I found that an optimal value is around 0.0008. It is also set to this value by default if you decide not to specify otherwise. 

## LUT-method tracking 

So, the files you have by now are `SimSamples.root`, `Clusters.root`, `clusters.txt`, and `LUT.txt`. Now, to perform LUT-method tracking on `Clusters.root`: 

```bash
$ denv fire TrigScint/exampleConfigs/runLUTtracking.py
```

The config file `runLUTtracking.py` will require `LUT.txt` as input of course, and it will transform `Clusters.root` (a .ROOT file with branches only up to clustering) into `Tracks.root` (a .ROOT file with two branches for tracking: one for lut-method and one for the max-delta method in `TrigScintTrackProducer`).

## To reproduce the thesis plots:

To reproduce the plots in my thesis you will need the jupyter notebook file `/plots/plots.ipynb`. You can launch jupyter notebook using 

```bash
$ jupyter notebook
```

(I have windows 10 and used [this tutorial](https://code.adamdimech.au/jupyter-notebook-in-windows-subsystem-for-linux-wsl/) to get jupyter notebook working in WSL 2). Assuming you have the following files: 

* clusters.txt
* LUT.txt
* Tracks.root

for both the truth and non-truth cases, and for multiple `lut_threshold`s, you should be able to reproduce all plots. **It is very important that these files are located in the same directory**.

## Other 

### Misalignment 

The TS geometry is specified in `runLUTeventgen.py`. To test a misaligned geometry, this would have to be changed. 

### Extra cluster in pad 2 

My paper mentions an extra cluster found in event number 7098 in pad 2. If you use `p.run_number=2` and `p.max_events=10000` (or atleast 7098) in `runLUTeventgen.py`, you can observe this for yourself. You can open root-files with the command 

```bash
$ denv rootbrowse examplefile.root
```

