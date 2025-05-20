from LDMX.Framework import ldmxcfg
import sys, os

p = ldmxcfg.Process('ana')
p.sequence = [ldmxcfg.Analyzer.from_file('PNAnalyzer.cxx')]

# list content in your dir with sys or os
files=[f for f in os.listdir('./OutputFiles')]
p.inputFiles=[]
for filepart in files:
	p.inputFiles.append('./OutputFiles/'+filepart)

# How to set a fixed number of inputfiles (also debugging remnant)
#p.inputFiles=['singleNpn_simreco_test250.root']
