# Disclaimer
ALL my code was run within the LUNARC HPC Desktop with access to both an ldmx container image and a Jupyter lab container image. I have no idea how the code reacts outside of HPC Desktop. Additionally, no simulation data files or accompanying sub-directories are provided, you will have to create your own work environment and tweak the names of files and directories to your liking. Especially the file `JobLauncher.sh` is highly personalized for the project and should be taken as a template rather than an actual file to run.

# Simulation and Analysis
The `singleNConfig.py`, `JobLauncher.sh`, `runAnalysisConfig.py` and `PNAnalyzer.cxx` can all be contained within the same directory with the container image initialized with
- `denv init ldmx/pro:v4.0.1`
I highly recommend creating a sub-directory to store the output files from the simulation runs. For me, this was the directory `OutputFiles` as can be seen in the remnants of the code. 

## Simulation file
The simulation file `singleNConfig.py` can be run directly from the terminal with the command
- `denv fire singleNConfig.py`
With up to 2 additional numeric arguments, such as
- `denv fire singleNConfig.py 1000 1`
The first of these arguments specifies the number of events to simulate. Within the file, its default value is 5000 events which takes around 2 hours to simulate and produces just shy of 500 MB of data. The second argument specifies the index of the output file to be created for the simulation data, as well as contributing to the rng seed within the file Within the file, its default value is 0. The output for the simulation is a single .root file named `singleNpn_simreco{index}.root` with whatever index was used. The data is stored in a TTree object, with subdirectories called "branches" containing lists of data in "leaves". One can view the contents of the .root file on the [online ROOT file reader by CERN](https://root.cern.ch/js/latest/). It can also be viewed from the terminal by loading the ROOT module and its parent directories with 
```
module load GCC/12.3.0
module load OpenMPI/4.1.5
module load ROOT/6.30.06
```
and opening a `TBrowser()` within the ROOT environment with
```
root
new TBrowser()
```
One exits the ROOT environment in the terminal by typing `.q`.

## Launching jobs
`JobLauncher.sh` is a simple file with commands for submitting a simulation job to the computing cluster. It can be run from the terminal with
- `denv fire JobLauncher.sh x`
and requires an additional numeric argument `x`. This argument is used as the index argument to `singleNConfig.py`, which is run from the file. It also jumps into the directory `OutputFiles` and runs the simulations from within, hence the need for such a directory or a substitute. It also contains commands to the cluster on where to send terminal output and error messages, as well as which user to notify upon conclusion of the job. These commands are as follows:
- `\#SBATCH -t hh:mm:ss` Specify time required for simulation. You should always overshoot the expected time so the cluster doesn't interrupt the job. It will end itself automatically once it's complete always.
- `\#SBATCH -A` Specify the working group you access the cluster from.
- `\#SBATCH -J` The name of the job. More useful if you're submitting multiple different jobs from the same account and want to keep track of which job finishes when, but also forms part of the output message file and error message file
- `\#SBATCH -o` The path to where to store terminal output messages from the job, such as echo or print messages. The template uses %x to copy the name of the job and %j to copy the index of the job, so the file name would be something like `PhotoNeutrons_123456.out`.
- `\#SBATCH -e` The path to where to store terminal error messages from the job. This uses the same template as the output files.
- `\#SBATCH --mail-user=someone@somewhere.something` The email adress of the user who gets notified when the job has finished. Please change this first, I do not want to be pranked with your computing cluster jobs.
- `\#SBATCH --mail-type=END` What type of message to convey to the aforementioned email adress. I have used the keyword `END` to get notified when the job has finished. Other options include `BEGIN`, `FAIL`, `REQUEUE` and `ALL`.
See [LUNARC Documentation site](https://lunarc-documentation.readthedocs.io/en/latest/manual/submitting_jobs/manual_specifying_requirements/) for more options when submitting jobs.


## Running Analysis Script
`runAnalysisConfig.py`is a simple file for creating an LDMX analyzer process from a separate analyzer file and giving it input files. It can be run from the terminal with 
- `denv fire runAnalysisConfig.py`
It fetches simulation files from a given directory, in my case `OutputFiles`, and places them in the `inputFiles` parameter of the analyzer process object. You can also directly place files in its input by specifying their name and path from the current directory.

## Analysis Script
`PNAnalyzer.cxx` is the analysis file used to separate various parameters from the simulation and place them into a new .root file. I do not know how to run it independently from the config script, see previous section. The input files is a list of .root files from the `singleNConfig.py` simulation file sent in with `runAnalysisConfig.py`. The analyzer processes each event in each file individually in the `analyze` method. The output of the analysis is a single .root file named `photoneutron\_8GeV.root` with a TTree containing whatever data I've ask it to store. I also manually add the date of creation to the name of the file to keep track of different versions of data (a remnant of which was left within the plotting notebook). If you want to create new variables to store, you need to:
- Create a variable name and initialize it among the `public` class variables.
- Add a branch to the TTree in the analyzer `onProcessStart()` method.
- Make sure the variable is reset and/or changes for each event within the analyzer `analyze()` method.
See CERN's [ROOT Manual](https://root.cern/manual/) for more functionalities within ROOT.

# Plotting
All the plots created for my thesis were made in the `AnalysisConfig.ipynb` Jupyter Lab notebook. It had no dependence on ldmx-sw and should be stored in a separate directory from the other files. It does, however have a python dependency in denv and needs the container image
- `denv init python:3.11`
in its directory. Jupyter Lab was installed _in denv_ with
- `denv pip install scikit-hep jupyterlab`.
Attempting to install this without denv will not give Jupyter lab access to the container image. Now, upon running
- `denv jupyter lab`,
(and after a boot-up sequence) a localhost link will appear which the user can copy into any browser. This is a Jupyter Lab browser which can run the .ipynb notebook. `AnalysisConfig` requires a .root data file to open from the analysis script, see previous section. I found it easiest to simply copy the .root file from the directory with the other files to the directory with the notebook. Be aware that the notebook version in this GitHub page requests a file with its date of creation appended. You will have to change this to your own file in order to run it properly.

