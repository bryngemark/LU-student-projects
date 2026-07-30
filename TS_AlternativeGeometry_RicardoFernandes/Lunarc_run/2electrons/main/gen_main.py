from LDMX.Framework import ldmxcfg
import os
import sys

p = ldmxcfg.Process('ecal_pn')



p.run = int(sys.argv[1])
p.max_events = int(sys.argv[2])
filename = "2main"
run_number=int((sys.argv[1]))-3
#Filename will have the format (filename)(runnumber).root
#Feel free to modify filename with your own variables, or
#modify it to be more descriptive

p.output_files = ["data/{0}{1}.root".format(filename,int(run_number))]

p.max_tries_per_event = 1000


# p.run = 3
# p.max_events = 20


p.logger.term_level = 4


from LDMX.Biasing import ecal
from LDMX.SimCore import generators as gen


my_sim = ecal.photo_nuclear('ldmx-vertTS-v14-8gev',gen.single_8gev_e_upstream_tagger())
my_sim.description = 'ECal PN Test Simulation'

p.sequence = [ my_sim ]

import LDMX.Ecal.ecal_geometry
import LDMX.Hcal.hcal_geometry


# p.output_files = ['data2/main/main.root']