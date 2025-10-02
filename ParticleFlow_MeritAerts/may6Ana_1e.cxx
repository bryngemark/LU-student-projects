#include "Framework/EventProcessor.h"
#include "Recon/ParticleFlow.h"
#include "Recon/Event/PFCandidate.h"
#include "Recon/Event/CaloCluster.h"
#include "Ecal/Event/EcalHit.h"
#include <iostream>
#include <fstream>
#include <string>
#include "TH3F.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLegend.h"

class may6Ana_1e: public framework::Analyzer {
    public:
    may6Ana_1e(const std::string& name, framework::Process& p)
       : framework::Analyzer(name, p) {}
     ~may6Ana_1e() override = default;
     void onProcessStart() override;
     void analyze(const framework::Event& event) override;

     TH2F* hist2D[10];
     TH1F* hist1D[10];
     int eventnr=0;
   };

   void may6Ana_1e::onProcessStart() {
    getHistoDirectory();

    hist1D[0]=new TH1F("trackE_pid3"," ",100,0,8500);
    hist1D[1]=new TH1F("trackE_!pid3"," ",100,0,8500);
    hist1D[2]=new TH1F("EcalE_pid3"," ",100,0,12000);
    hist1D[3]=new TH1F("EcalE_!pid3"," ",100,0,12000);
    hist1D[4]=new TH1F("dist_pid3"," ",100,0,20);
    hist1D[5]=new TH1F("dist_!pid3"," ",100,0,20);
    hist1D[6]=new TH1F("pid3_count"," ",3,0,3);
    hist1D[7]=new TH1F("low_energy_track"," ",100,0,700);

    hist2D[0]= new TH2F("clusVtrack_e"," ",100,0,12000,100,0,12000);

   }

   

   void may6Ana_1e::analyze(const framework::Event& event) {

    std::vector<ldmx::PFCandidate> cand{event.getCollection<ldmx::PFCandidate>("PFCandidates")};
    
    // pid count 
    int count_pid3=0;
    for (ldmx::PFCandidate hit : cand){
       if (hit.getPID()==3 || hit.getPID()==7){count_pid3++;};
      };
    hist1D[6]->Fill(count_pid3);



    // energy distributions
    for (ldmx::PFCandidate hit : cand){
      float trackE, clusE; 

      if (hit.getPID()==3){
        std::vector<float> p=hit.getTrackPxPyPz();
        float RawE=hit.getEcalRawEnergy(), ptot=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);

        hist1D[2]->Fill(RawE);
        hist1D[0]->Fill(ptot);
        hist2D[0]->Fill(RawE,ptot);
      };

      if (count_pid3==0){
        if(hit.getPID()==1){
          std::vector<float> p=hit.getTrackPxPyPz();
          float ptot=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);

          hist1D[1]->Fill(ptot);
        };

        if(hit.getPID()==2){
          hist1D[3]->Fill(hit.getEcalRawEnergy());
        };
      };
    };



    // distances
    if (count_pid3==1){
      for (ldmx::PFCandidate hit : cand){
        if (hit.getPID()==3){
          std::vector<float> t_xyz=hit.getEcalPositionXYZ(), e_xyz=hit.getEcalClusterXYZ();
          float dist=sqrt(pow(t_xyz[0]-e_xyz[0],2)+pow(t_xyz[1]-e_xyz[1],2));

          hist1D[4]->Fill(dist);
        };
      };
    };

    if (count_pid3==0){
      ldmx::PFCandidate x=cand[0];
      std::vector<float> t_xyz;
      t_xyz.push_back(0);

      for(ldmx::PFCandidate hit : cand){
        if (hit.getPID()==1){
          t_xyz=hit.getEcalPositionXYZ();
        }
        if (hit.getPID()==2){
          if (hit.getEcalRawEnergy()>x.getEcalRawEnergy()){x=hit;};
        };
      };
      std::vector<float> e_xyz=x.getEcalClusterXYZ();
      
      if (t_xyz[0]){
        float dist=sqrt(pow(t_xyz[0]-e_xyz[0],2)+pow(t_xyz[1]-e_xyz[1],2));
        hist1D[5]->Fill(dist);
      };
    };

    if(count_pid3==0){
      for (ldmx::PFCandidate hit : cand){
        if (hit.getPID()==1){
          std::vector<float> xyz=hit.getEcalPositionXYZ();
          float dist=sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
          hist1D[7]->Fill(dist);
        };
      };
    };

    eventnr++;

   };

DECLARE_ANALYZER(may6Ana_1e);