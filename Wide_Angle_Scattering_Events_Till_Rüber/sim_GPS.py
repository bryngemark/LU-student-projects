#sim_GPS.py

import os
import sys

#set up parameters for this run

if len(sys.argv) > 1 :
    nEvents=int(sys.argv[1])
else :
    nEvents=10

if len(sys.argv) > 2 :
    runNb=int(sys.argv[2])
else :
    runNb=1

# need the configuration python module
from LDMX.Framework import ldmxcfg

# Create the necessary process object (call this pass "sim")
p = ldmxcfg.Process( "sim" )

# Set the unique number that identifies this run
p.run = 9001 

# import a template simulator and change some of its parameters
from LDMX.SimCore import generators
from LDMX.SimCore import simulator

mySim = simulator.simulator( "mySim" )
mySim.setDetector( 'ldmx-det-v14-8gev', True)
mySim.generators = [ generators.single_e_wide_angle_downstream_target(minTheta = 30, maxTheta = 70, minPhi = 0, maxPhi = 360) ] # theta in deg
mySim.description = 'I am a wide angle scattering simulation!'

# import processor templates
from LDMX.Hcal import hcal
from LDMX.Hcal import digi
from LDMX.Ecal import digi as ecaldigi
from LDMX.DetDescr.HcalGeometry import HcalGeometry
from LDMX.Ecal import EcalGeometry
from LDMX.Hcal import HcalGeometry
from LDMX.Hcal import hcal_hardcoded_conditions
from LDMX.Ecal import ecal_hardcoded_conditions
from LDMX.Ecal import vetos

hcalDigis = digi.HcalDigiProducer()
ecalDigis = ecaldigi.EcalDigiProducer()
hcalDigis.hgcroc.noise = False
hcalrec = digi.HcalRecProducer()
ecalrec = ecaldigi.EcalRecProducer()

# add them to the sequence
p.sequence = [
    mySim,
    ecalDigis,  
    ecalrec, 
    hcalDigis, 
    hcalrec, 
    #ecal_digi.EcalDigiProducer(),
    #ecal_digi.EcalRecProducer(),
    # hcal_digi.HcalDigiProducer(),
    # hcal_digi.HcalRecProducer(),
    digi.HcalSimpleDigiAndRecProducer()
    ]

# During production (simulation), maxEvents is used as the number
# of events to simulate.
# Other times (like when analyzing a file), maxEvents is used as
# a cutoff to prevent fire from running over the entire file.
p.maxEvents = nEvents

# how frequently should the process print messages to the screen?
p.logFrequency = 1

# I want to see all of the information messages (the default is only warnings and errors)
p.termLogLevel = 1

# give name of output file
p.outputFiles = [ f"sim_wide_angle_events_{runNb}.root" ]

# print process object to make sure configuration is correct
# at beginning of run and wait for user to press enter to start
p.pause()
