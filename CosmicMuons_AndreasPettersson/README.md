# Thesis Paper Link

 http://lup.lub.lu.se/student-papers/record/9195143

# What the script does

This script analyzes cosmic muons in the back-HCal. An algorithm for the side-HCal (and ECal) is not implemented.

The whole analysis uses truth-level data, but should be implemented for reconstruction. We found that the positioning was off in reconstruction algorithm. Instead of giving a single straight muon path, it was split in two paths that pointed in different directions, in most cases. There is a simple reconstruction algorithm that takes a weighted energy-average for each position in a bar, altough this is the same method that I used in my script (the simple reco has some smearing as well I believe, which I didn't include). I didn't use simple reco as it had a bug that has since been fixed. An implementation in a reconstruction algorithm that uses TOA-positions will give a better representation of reality, and should be a next step (HcalRecProducer). 

A least squares method fits a line to the track of cosmic muons. A quadratic fit method was started but never finished (it finds a fit but it is not aligned to the data). I assume the problem lies in when I tried to plot the fit (curveFunc). Could be as simple as the equation I was using is wrong.

Three selections exist for MIM-like (minimally ionizing muons) events. These are angle_threshold_, R^2 and MinFCN. For explanation of these variables, see pages 10-12 in my thesis.

We can build a map that stores information for all scintillator bars in the back-HCal using the buildMap function.

With graphcheck enabled, the script plots the path of the muons in 3D plots.\
With estimatedcheck enabled, the script calculates the estimated energy in each bar.

We can analyse horisontal and vertical bars individually by changing chose_orientation_.

# How to run

Fork this repository and cd into it.

## Install denv and change image
install denv as according to https://ldmx-software.github.io/

Change image:\
denv init ldmx/pro:v4.2.12

I used ldmx/dev:4.2.2, but the script works for the pro version.\
See the ldmx-sw documentation for setting up the dev version.

## Generate an LHE-file

The LHE-file contains initial information about the muon.

run: (example command line with 10 events, choosing the backHcal)
denv python3 lheData/cosmic_muon_lhe_generator_v14.py --numEvents=10 --detector=backHcal

more options exist for the generator. See cosmic_muon_lhe_generator_v14.py.\
An LHE-file is now created in lheData/

## Simulate in Geant4

cosmics_config.py contains an option for max allowed events. Change this to the number of events you want to run. The run number needs to be changed if you want to have several files with unique events (altough I never tried this and can't guarantee that it works).

run:
denv fire cosmics_config.py lheData/name_of_file.lhe

This creates a ROOT-file in cosmicEvents/

## Analyse

config file uses CosmicsAnalysis.cxx & CosmicsAnalysis.h for analysis.

run:
denv fire Analysis/cosmics_ana_config.py cosmicEvents/name_of_file.root

Graphs and histograms are saved in Analysis/



