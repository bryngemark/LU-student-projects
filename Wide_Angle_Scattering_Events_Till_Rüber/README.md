I've been running my code using the latest development image (March 2025). To use LDMX code in general, one need to install "denv" and here is a descripton about how to do that: https://ldmx-software.github.io/using/getting-started.html . Using the "denv fire" command, scripts that can be found in the repository can be run. However, if one wants to run their personal developed code in the LDMX environment, one needs to use "just" instead of "denv". An installation guide to "just" can be found here: https://ldmx-software.github.io/developing/getting-started.html . There it is also described how to clone the current LDMX-SW repository so that one can make and test changes on the own computer.

The files used in this project work can be run through the following steps:

1. sim_GPS.py - "taskset -c 0,1,2 just fire sim_GPS.py 10000 60"

This produces 10000 wide angle scattering events and runs the reconstrcution
- "taskset -c 0,1,2" is used to restrict the job on 3 kernels as the code otherwise crashed my computer
- "just fire" says to run the code
- "10000" specifies the amount of events that we want to simulate
- "60" specifies the run number that we can set

2. do_pf.py - "taskset -c 0,1,2 just fire do_pf.py 10000 60"

This runs the particle flow algorithm and the parameters are the same as in the previous step

3. produceGraph.py - "taskset -c 0,1,2 just fire produceGraph.py 10000 60"

This uses the Analyzer to produce the graphs used in the thesis. It does so for single events like the reconstructed hits in a certain region and also plots for all events like the scoring plane hits. The produced file can be read via a root browser.

Notice that figures 10 and 11 show intermediate steps of the project work. To obtain these, the reconstructed hits must be modified or one uses the previous version of reconstructed hits without the fixes described in the project.
