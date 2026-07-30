#!/usr/bin/python

import os
import sys
import glob

from LDMX.Framework import ldmxcfg


# first, we define the process, which must have a name which identifies this
# processing pass ("pass name").
sim_pass_name="ecal_pn"
pileup_file_pass_name="pileup"
this_pass_name="overlay"
p=ldmxcfg.Process(this_pass_name)

# p.run = 5
# p.max_events = 500
det = 'ldmx-vertTS-v14-8gev'


pileup_files = sorted(glob.glob('../pileup/data/*.root'))
main_files = sorted(glob.glob('../main/data/*.root'))
if not pileup_files:
    raise RuntimeError("No pileup ROOT files found in data/pileup/")
if not main_files:
    raise RuntimeError("No main ROOT files found in data/main/")

run_number = int(sys.argv[1])
file_index = (run_number - 1) % len(pileup_files)

run_file_number=run_number-4

p.run = int(sys.argv[1])
p.max_events = int(sys.argv[2])
filename = "events"
#Filename will have the format (filename)(runnumber).root
#Feel free to modify filename with your own variables, or
#modify it to be more descriptive

p.input_files = [main_files[(run_number - 1) % len(main_files)]]
# p.input_files = ['../main/data/*.root']
p.output_files = ["data/{0}{1}.root".format(filename,int(run_file_number))]


# Load the full tracking sequance
from LDMX.Recon.overlay import OverlayProducer
from LDMX.Tracking import full_tracking_sequence


overlay = OverlayProducer('OverlayProducer')
overlay.overlay_filename = pileup_files[file_index]
# overlay.overlay_filename = '../data/pileup/*.root'
overlay.sim_passname = sim_pass_name
overlay.overlay_passname = pileup_file_pass_name

overlay.calo_collections = [
    'TriggerPad1SimHits',
    'TriggerPad2SimHits',
    'TriggerPad3SimHits',
    'TargetSimHits',
    'EcalSimHits',
    'HcalSimHits'
]

overlay.tracker_collections = [
    "TaggerSimHits",
    "RecoilSimHits",
]

p.sequence = [overlay]

# ECal geometry nonsense
import LDMX.Ecal.ecal_clusters as ecal_cluster
import LDMX.Ecal.ecal_hardcoded_conditions
import LDMX.Hcal.hcal_hardcoded_conditions
from LDMX.Ecal import digi as ecal_digi_reco
from LDMX.Ecal import ecal_geometry
from LDMX.Ecal import vetos as ecal_vetos

# Hcal hardwired/geometry stuff
from LDMX.Hcal import hcal_geometry


# this is hardwired into the code to be appended to the sim hits collections
overlay_str="Overlay"

# Load the TS modules
from LDMX.TrigScint.trig_scint import (
    TrigScintClusterProducer,
    TrigScintDigiProducer,
    trig_scint_track,
)


ts_digis = [
        TrigScintDigiProducer.pad1(),
        TrigScintDigiProducer.pad2(),
        TrigScintDigiProducer.pad3(),
        ]
for digi in ts_digis :
    digi.input_collection += overlay_str

ts_clusters = [
        TrigScintClusterProducer.pad1(),
        TrigScintClusterProducer.pad2(),
        TrigScintClusterProducer.pad3(),
        ]
for clu in ts_clusters :
    clu.input_pass_name = this_pass_name

trig_scint_track.input_pass_name = this_pass_name

#trig_scint_track = TrigScintTrackProducer( "trig_scint_track" )
trig_scint_track.number_horizontal_bars = 16
trig_scint_track.horizontal_bar_gap = 2.1
trig_scint_track.number_vertical_bars = 8 
trig_scint_track.vertical_bar_gap = 0.1
# trig_scint_track.verbosity = 1000
trig_scint_track.seeding_collection="TriggerPad1Clusters"
trig_scint_track.further_input_collections=["TriggerPad2Clusters","TriggerPad3Clusters"]
trig_scint_track.delta_vert_max= 0.
trig_scint_track.horizontal_bar_length=30.

# Load ElectronCounter and Trigger
from LDMX.Recon.electron_counter import ElectronCounter
from LDMX.Recon.simple_trigger import TriggerProcessor


count = ElectronCounter(
    simulated_electron_number=2,
    instance_name="ElectronCounter",
    input_pass_name="",
)

count.input_pass_name = this_pass_name

# Load HCAL veto
import LDMX.Hcal.hcal as hcal

# Load and configure  particle flow sequence.
# Here we use PF "tracking" and CLUE Ecal clustering
from LDMX.Recon import pf_reco


# track_pf = pf_reco.PFTrackProducer()
# "EcalScoringPlaneHitsOverlay" #
# track_pf.input_track_coll_name=track_pf.input_track_coll_name+overlay_str
# track_pf.input_pass_name=this_pass_name
# track_pf.do_electron_tracking=True
# # reference info
# truth_pf = pf_reco.PFTruthProducer()

# # CLUE
# import LDMX.Ecal.ecal_clusters as cl


# cluster = cl.EcalClusterProducer()
# cluster.seed_threshold = 350.
# cluster.dc = 0.3
# cluster.nbr_of_layers = 1
# cluster.reclustering = True
# cluster.rec_hit_pass_name=this_pass_name #run on process+pileup

# # particle flow:
# pf_comb=pf_reco.PFProducer()
# pf_comb.input_ecal_coll_name = cluster.cluster_coll_name # use CLUE
# # pf_comb.input_ecal_pass_name = this_pass_name
# # trigger recasting existing CLUE to caloclusters
# pf_comb.use_existing_ecal_clusters = True

# # Load pileup finder
# from LDMX.Recon import pileup_finder


# pu_finder = pileup_finder.PileupFinder()
# pu_finder.rec_hit_pass_name=this_pass_name
# #needs recast caloclusters, not (CLUE) ecalclusters
# pu_finder.cluster_coll_name=pf_comb.input_ecal_coll_name+"Cast"
# pu_finder.pf_cand_coll_name=pf_comb.output_coll_name
# pu_finder.min_momentum=3000.

# # Load the DQM modules
from LDMX.DQM import dqm


trig_scint_sim_dqm = [
    dqm.TrigScintSimDQM(
        instance_name="TrigScintSimPad1",
        hit_collection="TriggerPad1SimHits",
        pad="pad1",
    ),
    dqm.TrigScintSimDQM(
        instance_name="TrigScintSimPad2",
        hit_collection="TriggerPad2SimHits",
        pad="pad2",
    ),
    dqm.TrigScintSimDQM(
        instance_name="TrigScintSimPad3",
        hit_collection="TriggerPad3SimHits",
        pad="pad3",
    ),
]

for ts_sim_dqm in trig_scint_sim_dqm:
    ts_sim_dqm.hit_collection += overlay_str

trig_scint_dqm = [
    dqm.TrigScintDigiDQM(
        instance_name="TrigScintDigiPad1",
        hit_collection="trigScintDigisPad1",
        pad="pad1",
    ),
    dqm.TrigScintDigiDQM(
        instance_name="TrigScintDigiPad2",
        hit_collection="trigScintDigisPad2",
        pad="pad2",
    ),
    dqm.TrigScintDigiDQM(
        instance_name="TrigScintDigiPad3",
        hit_collection="trigScintDigisPad3",
        pad="pad3",
    ),
    dqm.TrigScintClusterDQM(
        instance_name="TrigScintClusterPad1",
        cluster_collection="TriggerPad1Clusters",
        pad="pad1",
    ),
    dqm.TrigScintClusterDQM(
        instance_name="TrigScintClusterPad2",
        cluster_collection="TriggerPad2Clusters",
        pad="pad2",
    ),
    dqm.TrigScintClusterDQM(
        instance_name="TrigScintClusterPad3",
        cluster_collection="TriggerPad3Clusters",
        pad="pad3",
    ),
    dqm.TrigScintTrackDQM(
        instance_name="TrigScintTracks",
        track_collection="TriggerPadTracks",
    ),
]

for ts_dqm in trig_scint_dqm:
    ts_dqm.pass_name = this_pass_name


# Trigger DQM
trigger_dqm = dqm.Trigger()
trigger_dqm.trigger_pass = this_pass_name


dqm_with_overlay = (
    trig_scint_sim_dqm
    + trig_scint_dqm
    # + [
    #     trigger_dqm,
    #     ecal_digi_verify,
    #     ecal_shower_features,
    #     ecal_mip_tracking_features,
    #     ecal_veto_results,
    # ]
    # + hcal_dqm
)

p.logger.term_level = 0

# # Add full tracking for both tagger and recoil trackers: digi, seeds, CFK,
# # ambiguity resolution, GSF, DQM
# from LDMX.Tracking import full_tracking_sequence


# # append "Overlay" to sim collection names in tracking sequence
# full_tracking_sequence.set_overlay(this_pass_name)
# p.sequence.extend(full_tracking_sequence.sequence)
# p.sequence.extend(full_tracking_sequence.dqm_sequence)

p.sequence.extend([
        *ts_digis,
        *ts_clusters,
        trig_scint_track,
        count,
])

p.sequence.extend(dqm_with_overlay)

# # Add PFlow + pileup finding sequence
# p.sequence.extend(
#     [cluster, dqm.EcalClusterAnalyzer(), track_pf, truth_pf, pf_comb, pu_finder]
# )


# p.input_files = ['data2/main/*.root'] #ecal_pn
# p.output_files= ['data2/config/events.root']



p.histogram_file = 'hist.root'
