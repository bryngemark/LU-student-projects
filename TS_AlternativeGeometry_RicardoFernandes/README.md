# Alternative Trigger Scintillator

## Summary of my work

I worked on updating and optimizing the Track Producer (`TrigScint/src/TrigScint/TrigScintTrackProducer.cxx`), particularly for the alternative geometry (`ldmx-vertTS-v14-8gev`). Please see my commit from July 2, 2026:

https://github.com/LDMX-Software/ldmx-sw/commit/cd380014bc29c83ff0834eff87641f06c194ab4b

There is also an example configuration file for this geometry: `TrigScint/exampleConfigs/runTSvert.py`.

My main focus was comparing the Trigger Scintillator (TS) performance between `ldmx-vertTS-v14-8gev` and `ldmx-det-v14-8gev`.

## How to run

Follow the LDMX tutorial and generate TS data for both geometries. The `runTSvert.py` file provides an example of how to run the new geometry. The baseline geometry can be run by replacing the detector name with the original one and removing the parameters specific to the new geometry (the default parameters correspond to the baseline geometry).

I used version 4.7.1 of `ldmx-sw` for this work. However, any version after 4.7.6 should already include my changes to the Track Producer.

You can run the simulation using either `denv fire` or `just fire`, depending on whether you are using the container environment or a local Git installation.

I used Lunarc to generate my data. See the `Lunarc run` folder for more details.

## Analysis files

I performed the analysis using Jupyter Notebooks. You can follow the LDMX tutorial for instructions on how to set up the environment and then launch Jupyter with:

```bash
denv jupyter lab
```

For more details, see `analysis-code_example.ipynb` (1 electron/event) and `analysis-code_example2.ipynb` (2 electrons/event).

The notebooks are not organized in the most user-friendly way, as they are intended only as examples of how to run the analysis. Most of the code consists of Python plotting and data analysis routines.

