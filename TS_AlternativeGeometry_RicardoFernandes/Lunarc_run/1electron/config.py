from LDMX.Framework import ldmxcfg


p = ldmxcfg.Process('test')

from LDMX.SimCore import simulator as sim


my_sim = sim.Simulator( "my_sim" )
det = 'ldmx-vertTS-v14-8gev'
my_sim.set_detector(det, include_scoring_planes_minimal = True )
from LDMX.SimCore import generators as gen

my_sim.generators.append( gen.single_8gev_e_upstream_tagger() )
my_sim.description = 'Basic test Simulation'

p.sequence = [ my_sim ]

##################################################################
# Below should be the same for all sim scenarios

import os
import sys

# p.run = 2
# p.max_events = 1000000
# p.output_files = ['events.root']

p.run = int(sys.argv[1])
p.max_events = int(sys.argv[2])
filename = "1alternative"
#Filename will have the format (filename)(runnumber).root
#Feel free to modify filename with your own variables, or
#modify it to be more descriptive

p.output_files = ["data/{0}{1}.root".format(filename,int(sys.argv[1]))]

p.histogram_file = 'hist.root'
# Load the full tracking sequence
import LDMX.Ecal.digi as ecal_digi
import LDMX.Ecal.ecal_clusters as ecal_cluster

# Load the ECAL modules
import LDMX.Ecal.ecal_geometry
import LDMX.Ecal.ecal_hardcoded_conditions
import LDMX.Ecal.vetos as ecal_vetos
import LDMX.Hcal.digi as hcal_digi_and_reco

# Load the HCAL modules
import LDMX.Hcal.hcal_geometry
import LDMX.Hcal.hcal_hardcoded_conditions
#from LDMX.Tracking import full_tracking_sequence

# Load the TS modules
from LDMX.TrigScint.trig_scint import (
        TrigScintClusterProducer,
        TrigScintDigiProducer,
        trig_scint_track
)


ts_digis = [
        TrigScintDigiProducer.pad1(),
        TrigScintDigiProducer.pad2(),
        TrigScintDigiProducer.pad3(),
        ]

ts_clusters = [
        TrigScintClusterProducer.pad1(),
        TrigScintClusterProducer.pad2(),
        TrigScintClusterProducer.pad3(),
        ]
    
# trig_scint_track = TrigScintTrackProducer( "trig_scint_track" )
trig_scint_track.number_horizontal_bars = 16
trig_scint_track.horizontal_bar_gap = 2.1
trig_scint_track.number_vertical_bars = 8 
trig_scint_track.vertical_bar_gap = 0.1
#trig_scint_track.verbosity = 1000
trig_scint_track.seeding_collection="TriggerPad1Clusters"
trig_scint_track.further_input_collections=["TriggerPad2Clusters","TriggerPad3Clusters"]
trig_scint_track.vertical_bar_start_index = 52 
#trig_scint_track.delta_max=1.34
trig_scint_track.delta_vert_max=0.
trig_scint_track.horizontal_bar_length=30.


# Load the DQM modules
from LDMX.DQM import dqm

# Load electron counting and trigger
from LDMX.Recon.electron_counter import ElectronCounter


count = ElectronCounter(
    simulated_electron_number=1,
    instance_name="ElectronCounter",
    input_pass_name="",
)



p.logger.term_level = 1#0
# p.logger.custom(ecal_veto, level = -1)

# Add full tracking for both tagger and recoil trackers: digi, seeds, CFK, ambiguity resolution, GSF, DQM
#p.sequence.extend(full_tracking_sequence.sequence)
#p.sequence.extend(full_tracking_sequence.dqm_sequence)

p.sequence.extend([
        #ecal_digi.EcalDigiProducer(),
        #ecal_digi.EcalRecProducer(),
        #ecal_pres_skimmer,
        #ecal_cluster.EcalClusterProducer(),
        #ecal_veto,
        #ecal_mip,
        #ecal_veto_pnet,
        #hcal_digi,
        #hcal_reco,
        #hcal_veto,
        *ts_digis,
        *ts_clusters,
        trig_scint_track,
        count, 
        #TriggerProcessor('trigger', 8000.),
        #dqm.PhotoNuclearDQM(),
        #dqm.EcalClusterAnalyzer()
        ])

#p.sequence.extend(dqm.all_dqm)
