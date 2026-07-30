import os
import sys
import time

#In this script, run numbers are set to increase sequentially
#You can do whatever you want to iterate run numbers, so long
#as they change to something unique in each iteration.

totaljobs = 20
eventsperjob = 1000
runnumber= 5
#The directory containing your justfile, with no / at the end
pathtojustfile='/home/ricardo/LDMX/ldmx-sw_custom'

#I made this script compatible with command line arguments
#If you run the script with no arguments like:
#
#>just fire submitjobs.py
#
#it will use the default values of totaljobs, eventsperjob,
#and runnumber. In this example script the default values
#are 100, 1000, and 983243, respectively
#
#Otherwise, if you provide 3 numerical arguments, like:
#
#> just fire submitjobs.py 50 100 203947
#
#The total jobs/starting run number/number of events per job
#will change accordingly
if(len(sys.argv)>3):
    print(sys.argv[0])
    totaljobs    = int(sys.argv[1])
    eventsperjob = int(sys.argv[2])
    runnumber    = int(sys.argv[3])


remainingjobs=totaljobs
    
while remainingjobs>0:
    os.system("sbatch --export=runnumber={0},eventsperjob={1},pathtojustfile={2} onejob.sh".format(runnumber,eventsperjob,pathtojustfile))
    remainingjobs=remainingjobs-1
    print("Job #{0}, run#{1} submitted.\n".format(totaljobs-remainingjobs,runnumber))
    runnumber=runnumber+1
    time.sleep(0.1)
