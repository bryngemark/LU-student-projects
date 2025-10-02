# aerts_thesis25

The code consists of two analyzers, *may6Ana_1e.cxx* and *apr28Ana_2e.cxx*, which treat the one-electron and two-electron case respectively. The analyzers are run through their respective config file, *may6Config_1e.py* and *apr28Config_2e.py*, where input and output files must be specified. 

LDMX code is run using *denv*. https://ldmx-software.github.io/using/getting-started.html contains a usefull guide for downloading and setting up *denv*. In order to run *may6Ana_1e.cxx*, one writes  
```
denv fire may6Config_1e.py
```

input files:
- *may6Config_1e.py* requires an iput file containing single-electron events, which has been run through ParticleFlow. For my thesis this was generated using the *cfg_pfReco* generator from *ldmx-helpers* by *therwig*.  
- *apr28Config_2e.py* requires an input file with two-electron events, also run through ParticlFlow. This is done by overlaying two single-electron events into one file, before running it through ParticleFlow.

Both analyzers produce .root files with histograms for analyzing the performance of ParticleFlow depending on amount of electrons and signal energy.

Link to thesis: http://lup.lub.lu.se/student-papers/record/9210730
