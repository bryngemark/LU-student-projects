from LDMX.Framework import ldmxcfg
import sys
import os

p = ldmxcfg.Process('pileup')


p.run = int(sys.argv[1])
p.max_events = int(sys.argv[2])
filename = "2pileup"
#Filename will have the format (filename)(runnumber).root
#Feel free to modify filename with your own variables, or
#modify it to be more descriptive

# p.run = 2
# # slightly less than the others to test wrapping
# p.max_events = 25
# p.output_files = ['data2/pileup/pileup.root']

p.logger.term_level = 4

from LDMX.SimCore import simulator as sim


my_sim = sim.Simulator( "my_sim" )
my_sim.set_detector( 'ldmx-vertTS-v14-8gev' )
from LDMX.SimCore import generators as gen


my_sim.generators.append( gen.single_8gev_e_upstream_tagger() )
my_sim.description = 'Basic test Simulation'

p.sequence = [ my_sim ]

import LDMX.Ecal.ecal_geometry
import LDMX.Hcal.hcal_geometry

run_number = int(sys.argv[1])-2

p.output_files = ["data/{0}{1}.root".format(filename,int(run_number))]