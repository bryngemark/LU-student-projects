#include "Framework/EventProcessor.h"
#include "SimCore/Event/SimParticle.h"
#include "SimCore/Event/SimCalorimeterHit.h"
#include "SimCore/Event/SimTrackerHit.h"
#include "Ecal/Event/EcalHit.h"
#include "Hcal/Event/HcalHit.h"
#include "Ecal/Event/EcalVetoResult.h"

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TDirectory.h"
#include "TH2F.h"


class PNAnalyzer : public framework::Analyzer {
 public:
  PNAnalyzer(const std::string& name, framework::Process& p)
    : framework::Analyzer(name, p) {}
  ~PNAnalyzer() override = default;
  void onProcessStart() override;
  void analyze(const framework::Event& event) override;
  void onProcessEnd() override;
  
  TFile* file=new TFile("photoneutron_8GeV.root", "RECREATE");
  TTree* photoNeutrons=new TTree("photoneutronData","PN");
  float electronE=0, gammaE=0, neutronE=0;
  float EcalVeto_nReadoutHits=0, EcalVeto_deepestLayerHit=0, EcalVeto_summedDet=0, EcalVeto_summedTightIso=0, EcalVeto_ecalBackEnergy=0, EcalVeto_discValue=0, neutronEndZ=0, neutronEndX=0, neutronEndY=0, recoilE=0;

  float EcalSimTotE=0, HcalSimTotE=0, EcalRecTotE=0, HcalRecTotE=0, EcalSimWeightE=0;

  float electronRestEnergy=0.511, neutronRestEnergy=939.6; // MeV

  const float layerZ[35]={246, 253, 271, 279, 297, 305, 325, 335, 353, 363, 383, 393, 411, 421, 441, 451, 469, 481,  499, 509, 529, 539, 557, 567, 589, 603, 625, 639, 661, 675, 697, 711, 733, 747, 800};
  const int numLayers=sizeof(layerZ)/sizeof(layerZ[0]);
  const float weights[34]={2.312, 4.312, 6.522, 7.49, 8.595, 10.253, 10.915, 10.915, 10.915, 10.915,  10.915, 10.915, 10.915, 10.915, 10.915, 10.915, 10.915, 10.915, 10.915, 10.915, 10.915,  10.915, 10.915, 14.783, 18.539, 18.539, 18.539, 18.539, 18.539, 18.539, 18.539, 18.539, 18.539, 9.938};
  const float mip_si_energy=0.13, scndOEC=1.0150996066489024;

  int printEvery=5000;
};

void PNAnalyzer::onProcessStart() {
 std::cout << "Electron is 11.\nPhoton is 22.\nNeutron is 2112.\n";
 std::cout << "Photonnuclear reaction is 9\n";

 photoNeutrons->Branch("electronE", &electronE, "electronE/F");
 photoNeutrons->Branch("gammaE", &gammaE, "gammaE/F");
 photoNeutrons->Branch("neutronE", &neutronE, "neutronE/F");
 photoNeutrons->Branch("EcalVeto_nReadoutHits", &EcalVeto_nReadoutHits, "EcalVeto_nReadoutHits/F");
 photoNeutrons->Branch("EcalVeto_deepestLayerHit", &EcalVeto_deepestLayerHit, "EcalVeto_deepestLayerHit/F");
 photoNeutrons->Branch("EcalVeto_summedDet", &EcalVeto_summedDet, "EcalVeto_summedDet/F");
 photoNeutrons->Branch("EcalVeto_summedTightIso", &EcalVeto_summedTightIso, "EcalVeto_summedTightIso/F");
 photoNeutrons->Branch("EcalVeto_ecalBackEnergy", &EcalVeto_ecalBackEnergy, "EcalVeto_ecalBackEnergy/F");
 photoNeutrons->Branch("EcalVeto_discValue", &EcalVeto_discValue, "EcalVeto_discValue/F");
 photoNeutrons->Branch("neutronEndZ", &neutronEndZ, "neutronEndZ/F");
 photoNeutrons->Branch("neutronEndX", &neutronEndX, "neutronEndX/F");
 photoNeutrons->Branch("neutronEndY", &neutronEndY, "neutronEndY/F");
 photoNeutrons->Branch("recoilE", &recoilE, "recoilE/F");

 photoNeutrons->Branch("EcalSimTotE", &EcalSimTotE, "EcalSimTotE/F");
 photoNeutrons->Branch("EcalSimWeightE", &EcalSimWeightE, "EcalSimWeightE/F");
 photoNeutrons->Branch("HcalSimTotE", &HcalSimTotE, "HcalSimTotE/F");

 photoNeutrons->Branch("EcalRecTotE", &EcalRecTotE, "EcalRecTotE/F");
 photoNeutrons->Branch("HcalRecTotE", &HcalRecTotE, "HcalRecTotE/F");
}

//All events have SimParticle 1 as the source electron at exactly 8 GeV
void PNAnalyzer::analyze(const framework::Event& event) {
 const auto& particles{event.getMap<int, ldmx::SimParticle>("SimParticles", "")};
 const ldmx::EcalVetoResult& EcalVeto{event.getObject<ldmx::EcalVetoResult>("EcalVeto", "")};
 
 const ldmx::SimParticle *sourceElectron=nullptr, *highEnergyGamma=nullptr, *PNNeutron=nullptr;
 sourceElectron=&particles.at(1);
 
 if(event.getEventNumber()%printEvery==0) std::cout << "Source electron with energy " << sourceElectron->getEnergy() << " MeV creates\n";
 // Search its daughters for high energy photon
 ldmx::SimParticle daughter;
 for(auto daughterID : sourceElectron->getDaughters()) { // 1
   
   // Make sure the daughter is stored
   if(particles.find(daughterID) != std::end(particles)) { // 2
     daughter = particles.at(daughterID);
     if(event.getEventNumber()%printEvery==0) std::cout << "\tDaughter of type " << daughter.getPdgID() << " has energy " << daughter.getEnergy() << std::endl; // 2.1

     // Find high energy photon daughter
     if(daughter.getPdgID()==22 && daughter.getEnergy()>5000) { // 3
       highEnergyGamma = &particles.at(daughterID);
       ldmx::SimParticle grandDaughter;

       // Search photon for photoneutron daughter
       for(auto grandDaughterID : daughter.getDaughters()) { // 4

         //Make sure granddaughter is stored
         if(particles.find(grandDaughterID) != std::end(particles)) { // 5
           grandDaughter = particles.at(grandDaughterID);
           if(event.getEventNumber()%printEvery==0) std::cout << "\t\tGranddaughter of type " << grandDaughter.getPdgID() << " has energy " << grandDaughter.getEnergy() << std::endl; // 5.1

           // Find photoneutron granddaughter
           if(grandDaughter.getPdgID()==2112 && grandDaughter.getProcessType()==ldmx::SimParticle::ProcessType::photonNuclear) { // 6
	     
            // Make sure only highest energy neutron is stored (is always last listed daughter?)
            if(PNNeutron==nullptr)  PNNeutron = &particles.at(grandDaughterID);
            if(grandDaughter.getEnergy() > PNNeutron->getEnergy()) PNNeutron = &grandDaughter;
           } // 6
         } // 5
       } // 4
       //break;
     } // 3
   } // 2
 } // 1


 if(PNNeutron) {
   electronE=sourceElectron->getEnergy(); // Store this minus rest mass?
   gammaE=highEnergyGamma->getEnergy();
   neutronE=PNNeutron->getEnergy()-neutronRestEnergy; // Better to store kinetic energy (rename?)
   EcalVeto_nReadoutHits=EcalVeto.getNReadoutHits();
   EcalVeto_deepestLayerHit=EcalVeto.getDeepestLayerHit();
   EcalVeto_summedDet=EcalVeto.getSummedDet();
   EcalVeto_summedTightIso=EcalVeto.getSummedTightIso();
   EcalVeto_ecalBackEnergy=EcalVeto.getEcalBackEnergy();
   EcalVeto_discValue=EcalVeto.getDisc();
   neutronEndZ=(PNNeutron->getEndPoint())[2];
   neutronEndX=(PNNeutron->getEndPoint())[0];
   neutronEndY=(PNNeutron->getEndPoint())[1];

   const auto& ecal_sim_hits{event.getCollection<ldmx::SimCalorimeterHit>("EcalSimHits")};
   const auto& hcal_sim_hits{event.getCollection<ldmx::SimCalorimeterHit>("HcalSimHits")};
   EcalSimTotE=0;
   EcalSimWeightE=0;
   HcalSimTotE=0;

   const auto& ecal_rec_hits{event.getCollection<ldmx::EcalHit>("EcalRecHits")};
   const auto& hcal_rec_hits{event.getCollection<ldmx::HcalHit>("HcalRecHits")};
   EcalRecTotE=0; // Make sure to reset total energy
   HcalRecTotE=0;
   
   int ESimIndex=0, HSimIndex=0, ERecIndex=0, HRecIndex=0;

   for(const auto& hit : ecal_sim_hits) {
     EcalSimTotE+=hit.getEdep();
     int layer=-1;
     for(int i=numLayers-1; i>=0; i--) {
       if(layerZ[i]<hit.getPosition()[2]) {
        layer=i;
        break;
       }
     }

     //Assuming sim getEdep() is rec energy_deposited_in_si
     float energy_deposited_in_si=hit.getEdep();
     float num_mips_equivalent=energy_deposited_in_si/mip_si_energy;

     EcalSimWeightE+=(energy_deposited_in_si+num_mips_equivalent*weights[layer])*scndOEC;

     // Debugging remnant
     /*if(event.getEventNumber()==1) {
      std::cout << EcalSimTotE << " * " << weights[layer] << " = " << EcalSimWeightE << 
      " for layer " << layer << " at z = " << hit.getPosition()[2] << std::endl;
     }*/
   }
   for(const auto& hit : hcal_sim_hits) {
     HcalSimTotE+=hit.getEdep();
   }
   
   for(const auto& hit : ecal_rec_hits) {
     EcalRecTotE+=hit.getEnergy();
   }
   for(const auto& hit : hcal_rec_hits) {
     HcalRecTotE+=hit.getEnergy();
   }

   const auto& EcalSurfHit{event.getCollection<ldmx::SimTrackerHit>("EcalScoringPlaneHits")};
   recoilE=0;
   for(const auto& hit : EcalSurfHit) {
     if(hit.getTrackID()==1 && hit.getPosition()[2]<253.0 && hit.getPdgID()==11) {
      recoilE=hit.getEnergy();
      break;
     }
   }
   
   photoNeutrons->Fill();
 } else {
   std::cout << "Event " << event.getEventNumber() << " generated no photoneutron!\n";
 }
}

void PNAnalyzer::onProcessEnd() {
  file->cd();
  photoNeutrons->Write("photoNeutrons");
  file->Close();
  delete file;
}

DECLARE_ANALYZER(PNAnalyzer);
