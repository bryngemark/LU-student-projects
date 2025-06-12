from LDMX.Framework import ldmxcfg

import sys, math

p = ldmxcfg.Process('pf') # create new process named pf
p.inputFiles = ['sim_wide_angle_events_60_withBug.root']

if len(sys.argv)>1:
    p.maxEvents = int(sys.argv[1])
else :
    p.maxEvents = 50

if len(sys.argv)>2:
    # p.inputFiles = [f'sim_wide_angle_events_pf_{sys.argv[2]}.root']
    p.histogramFile = f'graphs_singleEvent/hist_single_event_run{sys.argv[2]}.root'
else :
    p.inputFiles = ['sim_wide_angle_events_pf_60.root']
    p.histogramFile = f'graphs_singleEvent/hist_single_event_run.root'
 
p.sequence = [
    ldmxcfg.Analyzer.from_file('Analyzer_single_Event.cxx'),
    ]

