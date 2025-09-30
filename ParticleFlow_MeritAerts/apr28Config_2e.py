from LDMX.Framework import ldmxcfg

Amass='1000'

p = ldmxcfg.Process('ana')
p.sequence = [ ldmxcfg.Analyzer.from_file('aerts_thesis25/apr28Ana_2e.cxx') ]
p.inputFiles = [ '' ] #add input here!
p.histogramFile = 'output/apr28Hist_2e'+Amass+'MeV_test.root'