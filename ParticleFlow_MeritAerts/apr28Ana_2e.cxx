#include "Framework/EventProcessor.h"
#include "Recon/ParticleFlow.h"
#include "Recon/Event/PFCandidate.h"
#include "Recon/Event/CaloCluster.h"
#include "Ecal/Event/EcalHit.h"
#include "Event/TriggerResult.h"
#include <iostream>
#include <fstream>
#include <string>
#include "TH3F.h"
#include "TCanvas.h"
#include "TLegend.h"

class apr28Ana_2e: public framework::Analyzer {
    public:
    apr28Ana_2e(const std::string& name, framework::Process& p)
       : framework::Analyzer(name, p) {}
     ~apr28Ana_2e() override = default;
     void onProcessStart() override;
     void analyze(const framework::Event& event) override;

     int eventnr=0, plotnr=1; 

     TH1F* hist1D_;
     TH1F* hist1DE[2][3][3][3]; // [a] sorted by, [b] mom/(raw)energy, [c] 1st/2nd/3rd [d] trigger all/1/0
     TH3F* hist3D_s[5];
     TH3F* hist3D_o[5];
     TH3F* hist3D_c[5];
     TH2F* hist2D[10];
     TH1F* hist1Dp[10];
   };

   void apr28Ana_2e::onProcessStart() {
    getHistoDirectory();

    hist1D_=new TH1F("pid37count","am of pid3/7 hits per event",5,0,5);
    hist1Dp[0]=new TH1F("track_dist2","dist between tracks when 2 pid3",100,0,500);
    hist1Dp[1]=new TH1F("track_dist0","dist between tracks when 0 pid3",100,0,500);
    hist1Dp[2]=new TH1F("track_dist1","dist between tracks when 1 pid3",100,0,500);
    hist1Dp[3]=new TH1F("cluster_dist2","dist between clusters when 2 pid3",100,0,500);
    hist1Dp[4]=new TH1F("cluster_dist1","dist between clusters when 1 pid3",100,0,500);
    hist1Dp[5]=new TH1F("cluster_dist0","dist between clusters when 0 pid3",100,0,500);

    histograms_.create("clusterE"," ", 200,0,50); //checking for extreme low energy clustere (E~0 MeV)
    histograms_.create("pidcount1_match"," ",2,0,2);

    hist2D[2]=new TH2F("clusVStrack_dist_2"," ",100,0,500,100,0,500);
    hist2D[3]=new TH2F("clusVStrack_dist_1"," ",100,0,500,100,0,500);
    hist2D[4]=new TH2F("clusVStrack_dist_0", " ",100,0,500,100,0,500);
    hist2D[5]=new TH2F("(un)matched_dist"," ",100,0,60,100,0,60);
    hist2D[6]=new TH2F("clusVtrack_1signal", " ",100,0,500,100,0,500);
    hist2D[7]=new TH2F("clusVtrack_1pileup", " ",100,0,500,100,0,500);

    for (int i=0;i<plotnr;i++){
        std::string nr=std::to_string(i);
        hist3D_s[i]=new TH3F("3d_signal","signal",100,-200,200,100,-200,200,100,250,650);
        hist3D_o[i]=new TH3F("3d_overlay","signal",100,-200,200,100,-200,200,100,250,650);
        hist3D_c[i]=new TH3F("3d_clust","cluster centroid",100,-200,200,100,-200,200,100,250,650);
    };

    std::string name, desc=" ", tracknr, sort, trig;
    int high;
    for (int d=0;d<3;d++){
        for(int a=0;a<2;a++){
            for(int b=0;b<3;b++){
                for(int c=0;c<3;c++){
                    if(b==0){name="TrackE";high=8500;};
                    if(b==1){name="EcalE";high=8500;};
                    if(b==2){name="EcalRawE";high=12500;};

                    if(d==0){trig="_all";};
                    if(d==1){trig="_t0";};
                    if(d==2){trig="_t1";};

                    tracknr=std::to_string(c+1);
                    sort=(a?"E_":"T_");
                    name=sort+tracknr+name+trig;

                    hist1DE[a][b][c][d]=new TH1F(name.c_str(),desc.c_str(),200,0,high);
                };
            };
        };
    };

   };

   

   void apr28Ana_2e::analyze(const framework::Event& event) {

    std::vector<ldmx::PFCandidate> cand{event.getCollection<ldmx::PFCandidate>("PFCandidates")};
    std::vector<ldmx::SimTrackerHit> track{event.getCollection<ldmx::SimTrackerHit>("PFTracks")};
    ldmx::TriggerResult trig=event.getObject<ldmx::TriggerResult>("Trigger");
    int pass=trig.passed();

    int count_pid37=0;
    float m_e=0.51099895069; //MeV

    ldmx::PFCandidate hits_pid37[3];
    
    for (ldmx::PFCandidate hit : cand){
        histograms_.fill("clusterE",hit.getEcalRawEnergy());
    }

    //pid3 count
    for (ldmx::PFCandidate hit : cand){
        if(hit.getPID()==3 || hit.getPID()==7){
            hits_pid37[count_pid37]=hit;
            count_pid37++;
        };
    };
    hist1D_->Fill(count_pid37);



    //energy distributions
    if (count_pid37>1){

        std::sort(cand.begin(),cand.end(),
              [](ldmx::PFCandidate a, ldmx::PFCandidate b){
                  auto Pa=a.getTrackPxPyPz();
                  auto Pb=b.getTrackPxPyPz();
                  float Ea=sqrt(pow(Pa[0],2)+pow(Pa[1],2)+pow(Pa[2],2)), Eb=sqrt(pow(Pb[0],2)+pow(Pb[1],2)+pow(Pb[2],2));
                  return Ea>Eb;
              });

        for(int i=0;i<3;i++){

            std::vector<float> p=cand[i].getTrackPxPyPz();
            float track_E=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+m_e*m_e);

            hist1DE[0][0][i][0]->Fill(track_E);
            hist1DE[0][1][i][0]->Fill(cand[i].getEcalEnergy());
            hist1DE[0][2][i][0]->Fill(cand[i].getEcalRawEnergy());
            hist1DE[0][0][i][pass+1]->Fill(track_E);
            hist1DE[0][1][i][pass+1]->Fill(cand[i].getEcalEnergy());
            hist1DE[0][2][i][pass+1]->Fill(cand[i].getEcalRawEnergy());
        };

        std::sort(cand.begin(),cand.end(),
        [](ldmx::PFCandidate a, ldmx::PFCandidate b){return a.getEcalRawEnergy()>b.getEcalRawEnergy();});

        for(int i=0;i<3;i++){
            std::vector<float> p=cand[i].getTrackPxPyPz();
            float track_E=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+m_e*m_e);

            hist1DE[1][0][i][0]->Fill(track_E);
            hist1DE[1][1][i][0]->Fill(cand[i].getEcalEnergy());
            hist1DE[1][2][i][0]->Fill(cand[i].getEcalRawEnergy());
            hist1DE[1][0][i][pass+1]->Fill(track_E);
            hist1DE[1][1][i][pass+1]->Fill(cand[i].getEcalEnergy());
            hist1DE[1][2][i][pass+1]->Fill(cand[i].getEcalRawEnergy());
        };
        
    };


    
    //track distance
    if (track.size()==2){
        std::vector<float> pos1=track[0].getPosition(), pos2=track[1].getPosition();
        float track_dist=sqrt(pow(pos1[0]-pos2[0],2)+pow(pos1[1]-pos2[1],2)+pow(pos1[2]-pos2[2],2));

        if (count_pid37==2) {
            hist1Dp[0]->Fill(track_dist);
        } else if (count_pid37==1){
            hist1Dp[2]->Fill(track_dist);
        } else {
            hist1Dp[1]->Fill(track_dist);
        };
    };


    std::sort(cand.begin(),cand.end(),
    [](ldmx::PFCandidate a, ldmx::PFCandidate b){return a.getEcalRawEnergy()>b.getEcalRawEnergy();});
    
    //track vs cluster dist
    float clus_dist, track_dist;
    std::vector<float> xyz[2];
    if(track.size()==2 && cand[1].getEcalRawEnergy()>5){
        std::vector<float> pos1=track[0].getPosition(), pos2=track[1].getPosition();
        float track_dist=sqrt(pow(pos1[0]-pos2[0],2)+pow(pos1[1]-pos2[1],2));//+pow(pos1[2]-pos2[2],2));
        int ind=0;

        if (count_pid37==2){
            xyz[0]=hits_pid37[0].getEcalClusterXYZ();
            xyz[1]=hits_pid37[1].getEcalClusterXYZ();
        
            clus_dist=sqrt(pow(xyz[0][0]-xyz[1][0],2)+pow(xyz[0][1]-xyz[1][1],2));//+pow(xyz[0][2]-xyz[1][2],2));

            hist2D[2]->Fill(clus_dist,track_dist);
            hist1Dp[3]->Fill(clus_dist);
        } else {
            std::vector<float> clus1=cand[0].getEcalClusterXYZ(), clus2=cand[1].getEcalClusterXYZ();
            clus_dist=sqrt(pow(clus1[0]-clus2[0],2)+pow(clus1[1]-clus2[1],2));//+pow(clus1[2]-clus2[2],2));

            std::vector<double> p1=track[0].getMomentum(),p2=track[1].getMomentum();
            float E1=sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]+m_e*m_e), E2=sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]+m_e*m_e);

            if (E1>5 && E2>5 && count_pid37==0){
                hist2D[4]->Fill(clus_dist,track_dist);
                hist1Dp[5]->Fill(clus_dist);
            };
        }; 
    



    // ^signal v pileup
    if(track.size()==2 && cand[1].getEcalRawEnergy()>5 && count_pid37==1){
        float pid3TrackE, pid1TrackE;
        ldmx::PFCandidate hits[3]; 
        hits[2]=cand[0];
        for (ldmx::PFCandidate hit : cand){
            if (hit.getPID()==3 || hit.getPID()==7){
                std::vector<float> p=hit.getTrackPxPyPz();
                pid3TrackE=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
                hits[0]=hit;
            };
            if (hit.getPID()==1){
                std::vector<float> p=hit.getTrackPxPyPz();
                pid1TrackE=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
                hits[1]=hit;
            };
            if (hit.getPID()==2){
                if (hit.getEcalRawEnergy()>hits[2].getEcalRawEnergy()||hits[2].getPID()==3||hits[2].getPID()==7){hits[2]=hit;};
            }
        };
        float match=(pid3TrackE>pid1TrackE? 1:0);
        histograms_.fill("pidcount1_match",match);

        std::vector<float> Txyz1=hits[0].getEcalPositionXYZ(), Txyz2=hits[1].getEcalPositionXYZ();
        std::vector<float> Exyz1=hits[0].getEcalClusterXYZ(), Exyz2=hits[2].getEcalClusterXYZ();
        float dist_track=sqrt(pow(Txyz1[0]-Txyz2[0],2)+pow(Txyz1[1]-Txyz2[1],2));
        float dist_clus=sqrt(pow(Exyz1[0]-Exyz2[0],2)+pow(Exyz1[1]-Exyz2[1],2));

        if(pid3TrackE>pid1TrackE){
            hist2D[7]->Fill(dist_clus,dist_track);
        } else {
            hist2D[6]->Fill(dist_clus,dist_track);
        }
        hist2D[3]->Fill(dist_clus,dist_track);
        hist1Dp[4]->Fill(dist_clus);
    };



    // 1 pid3 distances
    if(count_pid37==1 && track.size()==2){
        std::vector<float> xyz=hits_pid37[0].getEcalPositionXYZ(), pxyz=hits_pid37[0].getTrackPxPyPz(), clusXYZ=hits_pid37[0].getEcalClusterXYZ();
        float tkXAtClus = xyz[0] + pxyz[0] / pxyz[2] * (clusXYZ[2] - xyz[2]), tkYAtClus = xyz[1] + pxyz[1] / pxyz[2] * (clusXYZ[2] - xyz[2]);
        float dist_pid3=sqrt(pow(tkXAtClus-clusXYZ[0],2)+pow(tkYAtClus-clusXYZ[1],2));

        bool foundClus=false, foundTrack=false;
        for (ldmx::PFCandidate hit:cand){ 
            if(!foundClus && hit.getPID()==2){
                clusXYZ=hit.getEcalPositionXYZ();
                foundClus=true;
            }
            if(!foundTrack && hit.getPID()==1){
                xyz=hit.getEcalPositionXYZ();
                pxyz=hit.getTrackPxPyPz();
                foundTrack=true;
            };
        }
        tkXAtClus = xyz[0] + pxyz[0] / pxyz[2] * (clusXYZ[2] - xyz[2]);
        tkYAtClus = xyz[1] + pxyz[1] / pxyz[2] * (clusXYZ[2] - xyz[2]);
        float dist_notpid3=sqrt(pow(tkXAtClus-clusXYZ[0],2)+pow(tkYAtClus-clusXYZ[1],2));
        
        hist2D[5]->Fill(dist_pid3,dist_notpid3);
    };



    //3d plots
    if(eventnr<plotnr){
        std::vector<ldmx::EcalHit> EcalO{event.getCollection<ldmx::EcalHit>("EcalRecHits","overlay")}, EcalS{event.getCollection<ldmx::EcalHit>("EcalRecHits","signal")};

        for(ldmx::EcalHit hit : EcalO){
            float x=hit.getXPos(), y=hit.getYPos(), z=hit.getZPos();
            hist3D_o[eventnr]-> Fill(x,y,z);
        };
        for(ldmx::EcalHit hit : EcalS){
            float x=hit.getXPos(), y=hit.getYPos(), z=hit.getZPos();
            hist3D_s[eventnr]-> Fill(x,y,z);
        };
        for(ldmx::PFCandidate hit : cand){
            std::vector<float> xyz=hit.getEcalClusterXYZ();
            hist3D_c[eventnr]-> Fill(xyz[0],xyz[1],xyz[2]);
        };
        hist3D_o[eventnr]->SetMarkerStyle(8);
        hist3D_s[eventnr]->SetMarkerStyle(8);
        hist3D_c[eventnr]->SetMarkerStyle(8);
    };

    eventnr++;

   };
};

DECLARE_ANALYZER(apr28Ana_2e);