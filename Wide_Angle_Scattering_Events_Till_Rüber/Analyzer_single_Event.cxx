// filename: Analyzer_single_Event.cxx

#include <iostream>
#include <iomanip>

#include "Framework/EventProcessor.h"

#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TObject.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TPolyLine3D.h>
#include <TView3D.h>
#include <TView.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TGaxis.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TStyle.h>

#include "Ecal/Event/EcalHit.h"
#include "Hcal/Event/HcalHit.h"
#include "SimCore/Event/SimTrackerHit.h"
#include "SimCore/Event/SimCalorimeterHit.h"

#include "Recon/Event/CaloCluster.h"
#include "Recon/Event/PFCandidate.h"

class Analyzer_single_Event : public framework::Analyzer {
  public:
  Analyzer_single_Event(const std::string& name, framework::Process& p)
      : framework::Analyzer(name, p) {}
    ~Analyzer_single_Event() override = default;
    
    void onProcessStart() override;
    void analyze(const framework::Event& event) override;
    void onProcessEnd() override;

  private:

  TGraph* graphEcalSPHits = new TGraph();
  TGraph* graphHcalSPHits = new TGraph();

  TGraph* graph = new TGraph();
  TGraph* graph1_2 = new TGraph();
  TGraph* graph2 = new TGraph();
  TGraph* graph3 = new TGraph();
  TGraph* graph4 = new TGraph();
  TGraph* graph5 = new TGraph();
  TGraph* graph5_2 = new TGraph();
  
  TH2D* graph6 = new TH2D("h2,","title,", 200, 0.0,5.0, 200, 0.0,5.0);
 

  THStack* hs = new THStack("hs","Stacked 1D x histograms");
  //histE = new TH1D("histE", "Something", 100 );
  TH1F* histXE_simhits = new TH1F("histXsimhits","Energy over X Position",100, -1100,1100);
  TH1F* histXE_simplerechits = new TH1F("histXsimplerechits","Energy over X Position",100, -1100,1100);

  TH1F* histXE_diff = new TH1F("histXE_diff","Energy over X Position",100, -1100,1100);
  
  TH1F* histE_clusters = new TH1F("histE_clusters","Energy over X Position",100, -1100,1100);

  TH1F* hist_nPID = new TH1F("hist_nPID", "PID entries", 7, 1, 8);
  TH1F* hist_PID_E = new TH1F("hist_PID_E", "Energy of PID entries", 7, 1, 8);

  TH1F* hist_EcalXPos_Track = new TH1F("hist_EcalXPos_Track","X Pos Track at Ecal front",100,-500,500); // ecal has dx 880

  TH1F* hist_EcalYPos_Track = new TH1F("hist_EcalYPos_Track","Y Pos Track at Ecal front",100,-500,500); // ecal has dx 600

  TGraph* graph7 = new TGraph();
  TGraph* graph8 = new TGraph(); // Ecal scoring plane - only electrons
  TGraph* graph9 = new TGraph(); // Hcal scoring plane - only electrons
  

  TH1F* hist_HcalScoringPlane_XYE = new TH1F("hist_HcalScoringPlane_XYE","HcalScoringPlane XY and Energy",100,-500,500);
  
  int pid5_idx = 0;

  int nPoints = 0;
  int nEcalSPe = 0;
  int nHcalSPe = 0;
  int nSPHit = 0;

  TDirectory* histo = getHistoDirectory();

};

void Analyzer_single_Event::onProcessStart() {
  std::cout << "onProcessStart" << std::endl;
}

void Analyzer_single_Event::analyze(const framework::Event& event) {
  std::cout << "Analyze Event:" << event.getEventNumber() << std::endl;

  TDirectory* dir = histo->mkdir(Form("Event_%d", event.getEventNumber()));
  dir->cd();

  TH2F* hist_simHits = new TH2F("histSimHits", ";Z [mm];X [mm]", 
    100, 200, 1000,  // Z bins and range
    100, -1500, 1500); // X bins and range

  TH2F* hist_simpleRecHits = new TH2F("histSimpleRecHits", ";Z [mm];X [mm]", 
        100, 200, 1000,  // Z bins and range
        100, -1500, 1500); // X bins and range

  TH2F* hist_recHits = new TH2F("histRecHits", "Rec Hits;Z [mm];X [mm]", 
        100, 200, 1000,  // Z bins and range
        100, -1500, 1500); // X bins and range
  
  TH2F* hist_simpleRecHits_stripCenter = new TH2F("histSimpleRecHits_stripCenter", "; Z [mm];X [mm]", 
        100, 200, 1000,  // Z bins and range
        100, -1500, 1500); // X bins and range

  TH2F* hist_simpleRecHits_withBug = new TH2F("histSimpleRecHits_withBug", ";Z [mm];X [mm]", 
        100, 200, 1000,  // Z bins and range
        100, -1500, 1500); // X bins and range

   TGraph* graph_cluster = new TGraph();

   TGraphErrors* graph_clusterE = new TGraphErrors();

  hist_simHits->SetDirectory(dir);
  hist_simpleRecHits_stripCenter->SetDirectory(dir);
  hist_simpleRecHits_withBug->SetDirectory(dir);
  hist_simpleRecHits->SetDirectory(dir);
  hist_recHits->SetDirectory(dir);

  const auto& hcal_sim_hits = event.getCollection<ldmx::SimCalorimeterHit>("HcalSimHits");
  for (const auto& hit : hcal_sim_hits) {
    hist_simHits->Fill(hit.getPosition()[2],hit.getPosition()[0],hit.getEdep());
  }

  const auto& simple_rec_hits_stripCenter = event.getCollection<ldmx::HcalHit>("HcalSimpleRecHitsStripCenter");
  for (const auto& hit :simple_rec_hits_stripCenter){
    hist_simpleRecHits_stripCenter->Fill(hit.getZPos(),hit.getXPos(),hit.getEnergy());
  }

  const auto& simple_rec_hits_withBug = event.getCollection<ldmx::HcalHit>("HcalSimpleRecHitsWithBug");
  for (const auto& hit :simple_rec_hits_withBug){
    hist_simpleRecHits_withBug->Fill(hit.getZPos(),hit.getXPos(),hit.getEnergy());
  }


  const auto& simple_rec_hits_1 = event.getCollection<ldmx::HcalHit>("HcalSimpleRecHits");
  for (const auto& hit :simple_rec_hits_1){
    hist_simpleRecHits->Fill(hit.getZPos(),hit.getXPos(),hit.getEnergy());
  }

  const auto& rec_hits = event.getCollection<ldmx::HcalHit>("HcalRecHits");
  for (const auto& hit :rec_hits){
    hist_recHits->Fill(hit.getZPos(),hit.getXPos(),hit.getEnergy());
  }

  int nHCalCluster = 0;
  const auto& pfHCalClusters = event.getCollection<ldmx::CaloCluster>("PFHcalClusters");
  for (const auto& cluster : pfHCalClusters){
    graph_cluster->SetPoint(nHCalCluster,cluster.getCentroidZ(),cluster.getCentroidX());
    nHCalCluster+=1;
  }

  int nHCalClusterE = 0;

  double x0 = 0.0;
  double x1 = 0.0;
  double z0 = 240.0;
  double z1 = 0.0;

  const auto& pfCandidates = event.getCollection<ldmx::PFCandidate>("PFCandidates");
  for (const auto& entry : pfCandidates){
    
    if(entry.getPID()==1){
      
    } else if (entry.getPID()==2){
      
    } else if (entry.getPID()==3){
      
    } else if (entry.getPID()==4){
      
    } else if (entry.getPID()==5){
        graph_clusterE->SetPoint(nHCalClusterE,entry.getHcalClusterXYZ()[2],entry.getHcalClusterXYZ()[0]);
        graph_clusterE->SetPointError(nHCalClusterE,entry.getHcalClusterEXYZ()[2], entry.getHcalClusterEXYZ()[0]);
        nHCalClusterE +=1;
        
        z1 = entry.getHcalClusterXYZ()[2];
        x0 = entry.getEcalPositionXYZ()[0];
        x1 = entry.getEcalPositionXYZ()[0] + (entry.getTrackPxPyPz()[0]) * (fabs(z1 - z0)/entry.getTrackPxPyPz()[2]);
    } else if (entry.getPID()==6){
      
    } else if (entry.getPID()==7){

    };
    
  }

  TLine* extrapTrack = new TLine(z0, x0, z1, x1);
  extrapTrack->SetLineStyle(9); // dashed
  extrapTrack->SetLineWidth(5);
  extrapTrack->SetLineColor(kBlack);

  hist_simpleRecHits->SetTitle("; Z [mm]; X [mm]");

  graph_cluster->SetMarkerColor(kRed);
  graph_cluster->SetMarkerStyle(8);

  graph_clusterE->SetMarkerColor(kRed);
  graph_clusterE->SetMarkerStyle(21);
  graph_clusterE->SetLineColor(kRed);
  graph_clusterE->SetLineWidth(3);

  TLegend* legend3 = new TLegend(0.65, 0.75, 0.88, 0.88);
  legend3->AddEntry(graph_cluster, "Cluster", "p");
  legend3->AddEntry(graph_clusterE, "Matched Cluster", "lep");
  legend3->AddEntry(extrapTrack,"Track","l");

  TCanvas c1_recHits_cluster("c1_recHits_cluster","",800,600);
  hist_simpleRecHits->Draw("COLZ");
  graph_cluster->Draw("P SAME");
  graph_clusterE->Draw("P SAME");
  extrapTrack->Draw("same");

  legend3->Draw();

  c1_recHits_cluster.Write();

  // Sim Hits With Cluster and Track

  TCanvas c2_simHits_cluster("c2_simHits_cluster","",800,600);

  hist_simHits->GetXaxis()->SetRangeUser(230, 800);
  hist_simHits->GetYaxis()->SetRangeUser(0, 800);

  hist_simHits->Draw("COLZ");
  graph_cluster->Draw("P SAME");
  graph_clusterE->Draw("P SAME");
  extrapTrack->Draw("same");

  legend3->Draw();

  c2_simHits_cluster.Write();

  // Sim Hits

  

  TCanvas c3_simHits("c3_simHits","",800,600);

  gStyle->SetOptTitle(0); // Hides the title box in all subsequent plots

  hist_simHits->GetXaxis()->SetRangeUser(230, 800);
  hist_simHits->GetYaxis()->SetRangeUser(0, 800);

  hist_simHits->Draw("COLZ");
  //extrapTrack->Draw("same");

  c3_simHits.Write();

  // Simple Rec Hits Strip Center

  TCanvas c4_simpleRecHits_stripCenter("c4_simpleRecHits_stripCenter","",800,600);

  hist_simpleRecHits_stripCenter->GetXaxis()->SetRangeUser(230, 800);
  hist_simpleRecHits_stripCenter->GetYaxis()->SetRangeUser(0, 800);

  hist_simpleRecHits_stripCenter->Draw("COLZ");
  //extrapTrack->Draw("same");

  c4_simpleRecHits_stripCenter.Write();

  // Simple Rec Hits with Bug

  TCanvas c5_simpleRecHits_withBug("c5_simpleRecHits_withBug","",800,600);

  hist_simpleRecHits_withBug->GetXaxis()->SetRangeUser(230, 800);
  hist_simpleRecHits_withBug->GetYaxis()->SetRangeUser(0, 800);

  hist_simpleRecHits_withBug->Draw("COLZ");
  //extrapTrack->Draw("same");

  c5_simpleRecHits_withBug.Write();

  // Simple Rec Hits

  TCanvas c6_simpleRecHits("c6_simpleRecHits","",800,600);

  hist_simpleRecHits->GetXaxis()->SetRangeUser(230, 800);
  hist_simpleRecHits->GetYaxis()->SetRangeUser(0, 800);

  hist_simpleRecHits->Draw("COLZ");

  c6_simpleRecHits.Write();

    
  int nPFCand = 0;
  int nPIDTrack = 0;
  int nPIDEcal = 0;
  int nPIDTrackEcal = 0;
  int nPIDHcal = 0;
  int nPIDTrackHcal = 0;
  int nPIDEcalHcal = 0;
  int nPIDTrackEcalHcal = 0;

  double simple_rec_hits_energy_sum = 0.0;

  


  //const auto& pfCandidates = event.getCollection<ldmx::PFCandidate>("PFCandidates");
  for (const auto& entry : pfCandidates){
    nPFCand += 1;
    if(entry.getPID()==1){
      hist_nPID->Fill(1,1);
      hist_PID_E->Fill(1,entry.getTruthEnergy());
      nPIDTrack += 1;
    } else if (entry.getPID()==2){
      hist_nPID->Fill(2,1);
      hist_PID_E->Fill(2,entry.getEcalRawEnergy());
      nPIDEcal += 1;
    } else if (entry.getPID()==3){
      nPIDTrackEcal += 1;
      hist_nPID->Fill(3,1);
      hist_PID_E->Fill(3,entry.getEcalRawEnergy());
    } else if (entry.getPID()==4){
      nPIDHcal += 1;
      hist_nPID->Fill(4,1);
      hist_PID_E->Fill(4,entry.getHcalRawEnergy());
    } else if (entry.getPID()==5){
      nPIDTrackHcal += 1;
      double p = sqrt( pow(entry.getTrackPxPyPz()[0],2) + pow(entry.getTrackPxPyPz()[1],2) +  pow(entry.getTrackPxPyPz()[2],2));
      graph->SetPoint(pid5_idx, p, entry.getHcalRawEnergy());
      graph1_2->SetPoint(pid5_idx, p, simple_rec_hits_energy_sum);
      
      graph2->SetPoint(pid5_idx,entry.getDistTkHcalMatch(),entry.getHcalRawEnergy());
      graph3->SetPoint(pid5_idx,entry.getDistTkHcalMatch(),entry.getHcalRawEnergy()/p);
      graph4->SetPoint(pid5_idx,entry.getHcalSecondEnergy(),entry.getHcalRawEnergy());
      graph5->SetPoint(pid5_idx,entry.getDistTkHcalMatch(),entry.getHcalSecondEnergy()/entry.getHcalRawEnergy());
      graph5_2->SetPoint(pid5_idx,entry.getDistTkSecondHcalCluster() - entry.getDistTkHcalMatch(),entry.getHcalSecondEnergy()/entry.getHcalRawEnergy());
      graph6->Fill(entry.getDistTkHcalMatch(),entry.getHcalRawEnergy()/p,entry.getHcalRawEnergy());
      pid5_idx += 1;
      hist_EcalXPos_Track->Fill(entry.getEcalPositionXYZ()[0]);
      hist_EcalYPos_Track->Fill(entry.getEcalPositionXYZ()[1]);
      
      std::cout << "Ecal Front Hit (pid=5) Z: " << entry.getEcalPositionXYZ()[2] << std::endl;
      hist_nPID->Fill(5,1);
      hist_PID_E->Fill(5,entry.getHcalRawEnergy());
    } else if (entry.getPID()==6){
      nPIDEcalHcal += 1;
      hist_nPID->Fill(6,1);
      hist_PID_E->Fill(6,entry.getEcalRawEnergy() + entry.getHcalRawEnergy());
    } else if (entry.getPID()==7){
      nPIDTrackEcalHcal += 1;
      hist_nPID->Fill(7,1);
      hist_PID_E->Fill(7,entry.getEcalRawEnergy() + entry.getHcalRawEnergy());
    };
    
  }

  const auto& ecal_scoring_hits = event.getCollection<ldmx::SimTrackerHit>("EcalScoringPlaneHits");
  for(const auto& hit :ecal_scoring_hits){
    if(hit.getTrackID() == 1 && fabs(hit.getPosition()[2]- 240) < 0.1){
      std::cout << "EcalScoringPlane Primary particle detected: " << std::endl;
      std::cout << "Electron at EcalScoringPlane X: " << hit.getPosition()[0] << " Y: " << hit.getPosition()[1] << " Z: " << hit.getPosition()[2] << " E: "<< hit.getEnergy() << std::endl;
      std::cout << "Electron at EcalScoringPlane px: " << hit.getMomentum()[0] << " py: " << hit.getMomentum()[1] << " pz: " << hit.getMomentum()[2] << " p: "<< sqrt( pow(hit.getMomentum()[0],2) + pow(hit.getMomentum()[1],2)+  pow(hit.getMomentum()[2],2)) << std::endl;
    graph8->SetPoint(nEcalSPe,hit.getPosition()[0],hit.getPosition()[1]);
    nEcalSPe+=1;
    nSPHit+=1;
    graphEcalSPHits->SetPoint(nSPHit,hit.getPosition()[0],hit.getPosition()[1]);
    }
  }

  const auto& hcal_scoring_hits = event.getCollection<ldmx::SimTrackerHit>("HcalScoringPlaneHits");
  for(const auto& hit :hcal_scoring_hits){
    if(hit.getTrackID() == 1 && fabs(hit.getPosition()[2]- 240) < 0.1){
    std::cout << "HcalScoringPlane Primary particle detected: " << std::endl;
    std::cout << "Electron at HcalScoringPlane X: " << hit.getPosition()[0] << " Y: " << hit.getPosition()[1] << " Z: " << hit.getPosition()[2] << " E: "<< hit.getEnergy() << std::endl;
    std::cout << "Electron at HcalScoringPlane px: " << hit.getMomentum()[0] << " py: " << hit.getMomentum()[1] << " pz: " << hit.getMomentum()[2] << " p: "<< sqrt( pow(hit.getMomentum()[0],2) + pow(hit.getMomentum()[1],2)+  pow(hit.getMomentum()[2],2)) << std::endl;  
    graph9->SetPoint(nHcalSPe,hit.getPosition()[0],hit.getPosition()[1]);
    nHcalSPe+=1;
    nSPHit+=1;
    graphHcalSPHits->SetPoint(nSPHit,hit.getPosition()[0],hit.getPosition()[1]);
    
    }
  }
  
  histXE_diff->Sumw2(kFALSE);
  histE_clusters->Sumw2(kFALSE);

  double sim_hits_energy_sum = 0.0;
  const auto& sim_hits = event.getCollection<ldmx::SimCalorimeterHit>("HcalSimHits");
  for (const auto& hit :sim_hits){
    histXE_simhits->Fill(hit.getPosition()[0],hit.getEdep());
    sim_hits_energy_sum += hit.getEdep();
    histXE_diff->Fill(hit.getPosition()[0],hit.getEdep());
  }

  
  const auto& simple_rec_hits = event.getCollection<ldmx::HcalHit>("HcalSimpleRecHits");
  for (const auto& hit :simple_rec_hits){
    histXE_simplerechits->Fill(hit.getXPos(),-hit.getEnergy());
    simple_rec_hits_energy_sum += hit.getEnergy();
    histXE_diff->Fill(hit.getXPos(),-hit.getEnergy());
  }

  double hcal_cluster_energy_sum = 0.0;
  const auto& hcal_clusters = event.getCollection<ldmx::CaloCluster>("PFHcalClusters");
  for (const auto& cluster :hcal_clusters){
    hcal_cluster_energy_sum += cluster.getEnergy();
    histE_clusters->Fill(cluster.getCentroidX(),cluster.getEnergy());
  }

  std::cout << "SimHitsEnergySum: " << sim_hits_energy_sum << std::endl;
  std::cout << "SimpleRecHitsEnergySum: " << simple_rec_hits_energy_sum << std::endl;
  std::cout << "HcalClusterEnergySum: " << hcal_cluster_energy_sum << std::endl;
  
}

void Analyzer_single_Event::onProcessEnd() {
  histo->cd();

  TCanvas c1_EcalSideHcalSP("c1_EcalSideHcalSP","Scoring Plane Ecal and Side Hcal",800,600);

  graphEcalSPHits->SetTitle("Primary Electrons at ECal and Side HCal Scoring Plane;X [mm];Y [mm]");
  graphEcalSPHits->SetMarkerStyle(7);
  graphEcalSPHits->SetMarkerColor(kBlue);
  

  graphHcalSPHits->SetTitle("Primary Electrons at ECal and Side HCal Scoring Plane;X [mm];Y [mm]");
  graphHcalSPHits->SetMarkerStyle(7);
  graphHcalSPHits->SetMarkerColor(kRed);

  graphHcalSPHits->GetXaxis()->SetLimits(-1500, 1500); // Optional: set X manually
  graphHcalSPHits->SetMinimum(-1500);                  // Y-axis lower limit
  graphHcalSPHits->SetMaximum(1500);                  // Y-axis upper limit

  graphHcalSPHits->Draw("AP"); // First graph with axes
  graphEcalSPHits->Draw("P SAME"); // Second graph on top

  auto legend = new TLegend(0.65, 0.75, 0.88, 0.88); // x1,y1,x2,y2 (normalized coordinates)
  legend->AddEntry(graphEcalSPHits, "ECal SP Hits", "p");
  legend->AddEntry(graphHcalSPHits, "HCal SP Hits", "p");
  legend->Draw();

  gPad->SetGrid();
  c1_EcalSideHcalSP.Write();

  /*

  TCanvas c1_p_clusterE_ratio("c1_p_clusterEnergy", "Matched Track Hcal Cluster Ratio", 800, 600);
  //TCanvas->SetDirectory(histo);
  graph->SetLineColor(0);
  graph->SetMarkerStyle(24);
  graph->SetTitle("Hcal Cluster energy vs Track Momentum; Momentum [MeV]; Cluster Energy [MeV]");
  graph->Draw("AP");
  c1_p_clusterE_ratio.Update();
  c1_p_clusterE_ratio.Write();
  //graph->Write();

  TCanvas c1_2_p_totalHcalEnergy("c1_2_p_totalHcalEnergy", "Matched Track over Total Hcal Energy", 800, 600);
  //TCanvas->SetDirectory(histo);
  graph1_2->SetLineColor(0);
  graph1_2->SetMarkerStyle(24);
  graph1_2->SetTitle("Total Hcal Cluster energy vs Track Momentum; Momentum [MeV]; Total Energy [MeV]");
  graph1_2->Draw("AP");
  c1_2_p_totalHcalEnergy.Update();
  c1_2_p_totalHcalEnergy.Write();
  //graph->Write();

  TCanvas c2_dist_clusterE("c2_dist_clusterEnergy", "Hcal Cluster Energy over distance to matched Track", 800, 600);
  //TCanvas->SetDirectory(histo);
  graph2->SetLineColor(0);
  graph2->SetMarkerStyle(24);
  graph2->SetTitle("Hcal Cluster Energy over distance to matched Track; Dist [RMS of cluster]; Cluster Energy [MeV]");
  graph2->Draw("AP");
  c2_dist_clusterE.Update();
  c2_dist_clusterE.Write();
  //graph->Write();

  TCanvas c3_dist_ERatio("c3_dist_ERatio", "Hcal Cluster Energy/Track Momentum Ratio over Distance to Matched Track", 800, 600);
  //TCanvas->SetDirectory(histo);
  graph3->SetLineColor(0);
  graph3->SetMarkerStyle(24);
  graph3->SetTitle("Hcal Cluster Energy/Track Momentum Ratio over Distance to Matched Track; Dist [RMS of cluster]; E_cluster/p_track");
  graph3->Draw("AP");
  c3_dist_ERatio.Update();
  c3_dist_ERatio.Write();
  //graph->Write();
  */
  TCanvas c4_Energy_first_second_cluster("c4_Energy_first_second_cluster", "First (matched) vs second Hcal Cluster Energy", 800, 600);
  graph4->SetLineColor(0);
  graph4->SetMarkerStyle(7);
  graph4->SetMarkerColor(kBlue);
  graph4->SetTitle("First (matched) vs second HCal Cluster Energy; Second Cluster Energy [MeV]; Matched Cluster Energy [MeV]");
  graph4->Draw("AP");
  c4_Energy_first_second_cluster.Update();
  c4_Energy_first_second_cluster.Write();
  
  //graph->Write();
  
  TCanvas c5_1st_2nd_ratio_dist("c5_1st_2nd_ratio_dist", "2nd to 1st (matched) Hcal Cluster Energy vs Track Distance", 800, 600);
  graph5->SetLineColor(0);
  graph5->SetMarkerStyle(7);
  graph5->SetMarkerColor(kBlue);
  graph5->SetTitle("2nd to 1st (matched) Hcal Cluster Energy vs Track Distance; Dist [RMS of cluster]; 2nd / 1st (matched) Cluster Energy");
  graph5->Draw("AP");
  c5_1st_2nd_ratio_dist.Update();
  c5_1st_2nd_ratio_dist.Write();

  
  //graph->Write();

  TCanvas c5_2_1st_2nd_ratio_dist("c5_2_1st_2nd_ratio_dist", "2nd to 1st (matched) Hcal Cluster Energy vs 2nd to 1st Distance from Track", 800, 600);
  //TCanvas->SetDirectory(histo);
  graph5_2->SetLineColor(0);
  graph5_2->SetMarkerStyle(7);
  graph5_2->SetMarkerColor(kBlue);
  graph5_2->SetTitle("2nd to 1st (matched) Hcal Cluster Energy vs 2nd to 1st Distance from Track; 2nd Dist - 1st Dist [RMS of cluster]; 2nd / 1st (matched) Cluster Energy");
  graph5_2->Draw("AP");
  gPad->SetGrid();
  c5_2_1st_2nd_ratio_dist.Update();
  c5_2_1st_2nd_ratio_dist.Write();
  
  /*
  TCanvas c6_dist_ERatio_Ecluster("c6_dist_ERatio_Ecluster", "Hcal Cluster Energy/Track Momentum Ratio over Distance to Matched Track with Cluster Energy (colour)", 800, 600);
  //TCanvas->SetDirectory(histo);
  //graph6->SetTitle("1st (matched) to 2nd Hcal Cluster Energy vs trak distance; Dist [RMS of cluster]; 2nd Cluster Energy / 1st (matched) Cluster Energy");
  graph6->Draw("COLZ");
  c6_dist_ERatio_Ecluster.Update();
  c6_dist_ERatio_Ecluster.Write();
  //graph->Write();
  */

  TCanvas c2_simHits_simpleRecHits_cluster_EnergyDiff("c2_energyDiffX","energy diff",800,600);
  histXE_simhits->SetLineColor(kBlue);
  histXE_simplerechits->SetLineColor(kGreen+3);

  histXE_simhits->SetTitle("Energy Distribution along x; X [mm]; E [MeV]");

  histXE_simhits->Draw("hist");
  histXE_simplerechits->Draw("hist same");
  histXE_diff->SetMarkerStyle(2);
  histXE_diff->SetMarkerColor(kRed);
  histXE_diff->SetLineColor(0);
  histXE_diff->Draw("P same");
  histE_clusters->SetMarkerStyle(20);
  histE_clusters->SetLineColor(0);
  histE_clusters->Draw("P same");

  TLegend* legend2 = new TLegend(0.65, 0.75, 0.88, 0.88);
  legend2->AddEntry(histXE_simhits, "SimHits", "l");
  legend2->AddEntry(histXE_simplerechits, "SimpleRecHits", "l");
  legend2->AddEntry(histXE_diff, "Difference", "p");
  legend2->AddEntry(histE_clusters, "Clusters", "p");

  legend2->Draw();

  c2_simHits_simpleRecHits_cluster_EnergyDiff.Write();

  // PID

  TCanvas c3_nPID("c3_nPID","PID entries",800,600);
  c3_nPID.SetLogy();

  hist_nPID->SetTitle("PID Entries; PID; Counts");

  hist_nPID->GetXaxis()->SetBinLabel(1, "1");
  hist_nPID->GetXaxis()->SetBinLabel(2, "2");
  hist_nPID->GetXaxis()->SetBinLabel(3, "3");
  hist_nPID->GetXaxis()->SetBinLabel(4, "4");
  hist_nPID->GetXaxis()->SetBinLabel(5, "5");
  hist_nPID->GetXaxis()->SetBinLabel(6, "6");
  hist_nPID->GetXaxis()->SetBinLabel(7, "7");

  hist_nPID->GetXaxis()->SetLabelSize(0.05);

  hist_nPID->SetLineColor(kBlue);
  hist_nPID->SetLineWidth(2);
  hist_nPID->Draw("hist");

  gPad->SetGrid();

  c3_nPID.Write();

  // PID Energy

  TCanvas c4_PID_E("c4_PID_E","Energy of PID entries",800,600);
  c4_PID_E.SetLogy();

  hist_PID_E->SetTitle("Energy of PID Entries; PID; Energy [MeV]");

  hist_PID_E->GetXaxis()->SetBinLabel(1, "1");
  hist_PID_E->GetXaxis()->SetBinLabel(2, "2");
  hist_PID_E->GetXaxis()->SetBinLabel(3, "3");
  hist_PID_E->GetXaxis()->SetBinLabel(4, "4");
  hist_PID_E->GetXaxis()->SetBinLabel(5, "5");
  hist_PID_E->GetXaxis()->SetBinLabel(6, "6");
  hist_PID_E->GetXaxis()->SetBinLabel(7, "7");

  hist_PID_E->GetXaxis()->SetLabelSize(0.05);

  hist_PID_E->SetLineColor(kBlue);
  hist_PID_E->SetLineWidth(2);
  hist_PID_E->Draw("hist");

  gPad->SetGrid();

  c4_PID_E.Write();
    
  /// Single Sim Hits Event


  std::cout << "Entries in EcalScoringPlane:"<< nEcalSPe << std::endl;
  std::cout << "Entries in HcalScoringPlane:"<< nHcalSPe << std::endl;
  std::cout << "Entries in Total Scoring Plane:"<< nSPHit << std::endl; 
  std::cout << "The analysis has ended." << std::endl;
}

DECLARE_ANALYZER(Analyzer_single_Event);
