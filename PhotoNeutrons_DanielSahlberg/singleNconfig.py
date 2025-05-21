from LDMX.Framework import ldmxcfg
# create my process object
p = ldmxcfg.Process( "testSingleNpn" )

# how many events to process?
import sys
p.maxEvents = 5000
outputIndex = 0
if len(sys.argv) > 1 :
    p.maxEvents = int(sys.argv[1])
    if len(sys.argv) >2 :
    	outputIndex = int(sys.argv[2])
 
p.totalEvents = p.maxEvents
# we want to see every event
p.logFrequency = int(p.totalEvents/5)
p.termLogLevel = 1

# Set a run number
p.run = 9003+outputIndex

# we also only have an output file
p.outputFiles = [ "singleNpn_simreco"+str(outputIndex)+".root" ]

from LDMX.SimCore import simulator as sim
import LDMX.Ecal.EcalGeometry
import LDMX.Ecal.ecal_hardcoded_conditions
import LDMX.Hcal.HcalGeometry
import LDMX.Hcal.hcal_hardcoded_conditions
import LDMX.Ecal.digi as ecal_digi
import LDMX.Ecal.vetos as ecal_vetos
import LDMX.Hcal.digi as hcal_digi

mySim = sim.simulator( "mySim" )
mySim.setDetector( 'ldmx-det-v14' , True )
mySim.description = 'Basic test simulation for single neutron PN events'
# Get a pre-written generator
from LDMX.SimCore import generators as gen
mySim.generators.append( gen.single_8gev_e_upstream_tagger() )
from LDMX.SimCore import bias_operators
mySim.biasing_operators = [ bias_operators.PhotoNuclear('ecal',550.,5000.,only_children_of_primary = True) ]


from LDMX.SimCore import photonuclear_models as pn 
from LDMX.Biasing import particle_filter 

myModel = pn.BertiniSingleNeutronModel()
# These are the default values of the parameters
myModel.hard_particle_threshold=200. # Count particles with >= 200 MeV as "hard"
#myModel.zmin = 0 # Apply the model for any nucleus
myModel.zmin = 5 
myModel.emin = 5000. # Apply the model for photonuclear reactions with > 5000 MeV photons
myModel.count_light_ions=True # Don't disregard deutrons, tritons, 3He, or 4He ions when counting hard particles 

# Change the default model to the single neutron model
mySim.photonuclear_model = myModel


singleNFilter = particle_filter.PhotoNuclearTopologyFilter.SingleNeutronFilter()
# Default values 
singleNFilter.hard_particle_threshold = 200. # MeV, use the same as for the model 
singleNFilter.count_light_ions = True # As above


# the filters are in a library that needs to be included
from LDMX.Biasing import filters
from LDMX.Biasing import util
from LDMX.Biasing import include as includeBiasing
includeBiasing.library()

 # Configure the sequence in which user actions should be called.
mySim.actions.extend([
    filters.TaggerVetoFilter(),
    # Only consider events where a hard brem occurs
    filters.TargetBremFilter(3000., 5000.), 
    # Only consider events where a PN reaction happnes in the ECal
    filters.EcalProcessFilter(),     
    # Tag all photo-nuclear tracks to persist them to the event.
    util.TrackProcessFilter.photo_nuclear(),
# Add the filter at the end of the current list of user actions.
    singleNFilter
]) 


# add your configured simulation to the sequence
p.sequence= [ mySim, 
	ecal_digi.EcalDigiProducer(),
        ecal_digi.EcalRecProducer(), 
        ecal_vetos.EcalVetoProcessor(),
        hcal_digi.HcalDigiProducer(),
        hcal_digi.HcalRecProducer()
        ]


#import json
#json.dumps(p.parameterDump())

#with open('parameterDump.json', 'w') as outfile:
#     json.dump(p.parameterDump(),  outfile, indent=4)