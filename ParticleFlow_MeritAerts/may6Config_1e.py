from LDMX.Framework import ldmxcfg

p = ldmxcfg.Process('ana')
p.sequence = [ ldmxcfg.Analyzer.from_file('aerts_thesis25/may6Ana_1e.cxx') ]
p.inputFiles = [ '' ] #add input file here!
p.histogramFile = 'output/may6Hist_1e.root'