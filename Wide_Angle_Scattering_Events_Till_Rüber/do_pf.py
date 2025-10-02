from LDMX.Framework import ldmxcfg

import sys, math

p = ldmxcfg.Process('pf') # create new process named pf

if len(sys.argv)>1:
    p.maxEvents = int(sys.argv[1])
else :
    p.maxEvents = 100

if len(sys.argv) > 2 :
    p.inputFiles = [f'sim_wide_angle_events_{sys.argv[2]}.root']
    p.histogramFile = f'hist_pf_wide_angle_run_{sys.argv[2]}.root'
    p.outputFiles = [f'sim_wide_angle_events_pf_{sys.argv[2]}.root']
else :
    p.histogramFile = f'hist.root'
    # Set the output file name
    p.outputFiles = [f'sim_wide_angle_events_pf_4.root']

# Import the Hcal conditions
from LDMX.Hcal import hcal_hardcoded_conditions
from LDMX.Hcal import digi

# Geometry provider
import LDMX.Hcal.HcalGeometry

# HCal digi
hcalDigis = digi.HcalSimpleDigiAndRecProducer()
hcalDigis.position_resolution = 150 # 150 mm default

from LDMX.Recon import pfReco
     
ecalPF = pfReco.pfEcalClusterProducer()
hcalPF = pfReco.pfHcalClusterProducer()

trackEcalSPPF = pfReco.pfTrackProducer()
trackEcalSPPF.inputTrackCollName = 'EcalScoringPlaneHits'
trackHcalSPPF = pfReco.pfTrackProducer()
trackHcalSPPF.inputTrackCollName = 'HcalScoringPlaneHits'
truthPF = pfReco.pfTruthProducer()

# configure clustering options
ecalPF.doSingleCluster = False #True
ecalPF.logEnergyWeight = True

hcalPF.doSingleCluster = False
# hcalPF.clusterHitDist = 200. # mm
hcalPF.logEnergyWeight = True
hcalPF.hitCollName = 'HcalSimpleRecHits'

pfProducer = pfReco.pfProducer()
pfProducer.tkHadCaloMatchDist = 100.0
pfProducer.tkHadCaloMinEnergyRatio = 0.0
pfProducer.tkHadCaloMaxEnergyRatio = math.inf

# p.sequence.extend([
#      hcalDigis, # produced HcalSimpleRecHits
#     ])   
p.sequence = [
    ecalPF, hcalPF, trackEcalSPPF,
    # trackHcalSPPF, 
    pfProducer,
    truthPF,
    #ldmxcfg.Analyzer.from_file('Analyzer_creating_PFgraphs.cxx'),
    ]


p.keep = [ 
          # "drop .*SimHits.*",
          "drop EcalDigi.*",
          "drop HcalDigi.*",
          #"drop EcalRec.*",
          #"drop HcalRec.*",
          # slimmed into pfTruth collections instead
          "drop SimParticles.*",
          'drop TargetScoringPlaneHits.*',
          'drop TrackerScoringPlaneHits.*',
          'drop MagnetScoringPlaneHits.*',
          #'drop EcalScoringPlaneHits.*',
          #'drop HcalScoringPlaneHits.*',
          'drop TrigScintScoringPlaneHits.*',
     ]

allcols = [
     'stringParameters',
     'sampleOfInterest',
     'eventNumber',
     'TaggerSimHits',
     'EventHeader.floatParameters',
     'EventHeader.intParameters',
     'PFEcalClusters',
     'intParameters',
     'RecoilSimHits',
     'EventHeader.isRealData',
     'TriggerPad2SimHits',
     'TriggerPad1SimHits',
     'isRealData',
     'EventHeader.weight',
     'EventHeader.run',
     'TargetSimHits',
     'MagnetScoringPlaneHits',
     'EcalScoringPlaneHits',
     'HcalScoringPlaneHits',
     'TrigScintScoringPlaneHits',
     'version',
     'EventHeader',
     'HcalSimHits',
     'EcalDigis',
     'EventHeader.timestamp',
     'EventHeader.eventNumber',
     'timestamp',
     'numSamplesPerDigi',
     'HcalDigis',
     'EventHeader.stringParameters',
     'PFHcalClusters',
     'EcalSimHits',
     'PFTracks',
     'EcalRecHits',
     'weight',
     'samples',
     'floatParameters',
     'TriggerPad3SimHits',
     'run',
     'PFCandidates',
     'channelIDs',
     'SimParticles',
     'HcalRecHits',
     'HcalSimpleRecHits',
]

