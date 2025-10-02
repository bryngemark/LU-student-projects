# Modified LDMX-SW code
The following changes have been made:
1. generators.py
- added General Particle source with incident angle of 30-70 degrees w.r.t. z-axis to simulate wide angle scattering events
2. digi.py
- added position resolution variable for side hcal to have side-hcal-specific smearing
3. HcalSimpleDigiAndRecProducer.h
- added position resolution variable for side hcal to have side-hcal-specific smearing
4. HcalSimpleDigiAndRecProducer.cxx
- fixed bug with energy weighted average
- added reconstruction for side hcal based on strip position with side-hcal-specific smearing along the bar
5. pfReco.py
- added configuration variables in pfProducer to make matching criteria between track and hcal cluster configurable
6. PFTrackProducer.cxx
- added extra for-loop through HcalScoringPlaneHits so that tracks are produced from both Ecal AND Hcal scoring plane
7. ParticleFlow.h
- added configuration parameters for matching criteria and distance variable to extract distance from matched cluster to track
8. ParticleFlow.cxx
- added matching algorithm between track and hcal cluster (similar to the one for ecal clusters)
- added arbitration
9. PFCandidate.h
- added methods to set and get distance parameters from matching algorithm
