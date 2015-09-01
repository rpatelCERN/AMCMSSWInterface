// -*- C++ -*-
//
// Package:    L1TkPartAnalyzer
// Class:      TkTriggerParticleAnalzer
// 
/**\class TkTriggerParticleAnalzer TkTriggerParticleAnalzer.cc SLHCUpgradeSimulations/TkTriggerParticleAnalzer/src/TkTriggerParticleAnalzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Emmanuelle Perez,40 1-A28,+41227671915,
//         Created:  Thu Nov 14 11:22:13 CET 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"


// Gen-level stuff:
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkPrimaryVertex.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkHTMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkHTMissParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTRoad.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace l1extra;



//
// class declaration
//

class TkTriggerParticleAnalzer : public edm::EDAnalyzer {
   public:

   typedef TTTrack< Ref_PixelDigi_ >  L1TkTrackType;
   typedef std::vector< L1TkTrackType >  L1TkTrackCollectionType;

      explicit TkTriggerParticleAnalzer(const edm::ParameterSet&);
      ~TkTriggerParticleAnalzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      //virtual void endRun(edm::Run const&, edm::EventSetup const&);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
	float tp_pt, tp_eta, tp_phi, tp_z;
  //      float tk_pt, tk_eta, tk_phi, tk_z;   
  //      float mu_pt, mu_eta, mu_phi, mu_z;
	std::vector<float>tk_pt; std::vector<float>tk_eta;  std::vector<float>tk_phi; std::vector<float>tk_z;
	std::vector<float>tk_chi2; std::vector<int>tk_truth;
        std::vector<float>mu_pt; std::vector<float>mu_eta;  std::vector<float>mu_phi; std::vector<float>mu_z;
	std::vector<float>muStandAlone_pt; std::vector<float>muIso; 
	std::vector<float>ele_pt; std::vector<float>ele_eta;  std::vector<float>ele_phi; std::vector<float>ele_z;
	std::vector<float>eleIso;
        std::vector<float>ele_tkpt; std::vector<float>ele_tketa;  std::vector<float>ele_tkphi;
        std::vector<float>ele_calopt; std::vector<float>ele_caloeta;  std::vector<float>ele_calophi;

        std::vector<float>pho_pt; std::vector<float>pho_eta;  std::vector<float>pho_phi;
        std::vector<float>phoIso;
        std::vector<float>pho_calopt; std::vector<float>pho_caloeta;  std::vector<float>pho_calophi;
	TTree*AmMuons;
     // for L1TkMuonParticle

        edm::InputTag L1TkMuonsInputTag;
        edm::InputTag GenPartInputTag;
        edm::InputTag TrackPartInputTag;
	edm::InputTag TTTracksInputTag;
        edm::InputTag TTTracksAssocInputTag;
	edm::InputTag L1TkElectronsInputTag;	
        edm::InputTag L1TkPhotonsInputTag;	
	int PDG;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TkTriggerParticleAnalzer::TkTriggerParticleAnalzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  L1TkMuonsInputTag = iConfig.getParameter<edm::InputTag>("L1TkMuonsInputTag");
  GenPartInputTag   = iConfig.getParameter<edm::InputTag>("GenPartInputTag");
  TrackPartInputTag =iConfig.getParameter<edm::InputTag>("TrackPartTag");
  TTTracksInputTag  =iConfig.getParameter<edm::InputTag>("TTTracksInputTag");
  TTTracksAssocInputTag=iConfig.getParameter<edm::InputTag>("inputTagMC");
  L1TkElectronsInputTag = iConfig.getParameter<edm::InputTag>("L1TkElectronsInputTag");
  L1TkPhotonsInputTag= iConfig.getParameter<edm::InputTag>("L1TkPhotonsInputTag");
  PDG=iConfig.getParameter< int >("ParticleType");
}


TkTriggerParticleAnalzer::~TkTriggerParticleAnalzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TkTriggerParticleAnalzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


        //  
        // ----------------------------------------------------------------------
        // retrieve the L1TkMuons
        //
 edm::Handle<L1TkMuonParticleCollection> L1TkMuonsHandle;
 iEvent.getByLabel(L1TkMuonsInputTag, L1TkMuonsHandle);
 std::vector<L1TkMuonParticle>::const_iterator muIter;

edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle;
iEvent.getByLabel(TrackPartInputTag, TrackingParticleHandle);

edm::Handle< std::vector< TTTrack< Ref_PixelDigi_ > > > TTTrackHandle;
iEvent.getByLabel(TTTracksInputTag,TTTrackHandle);

edm::Handle< TTTrackAssociationMap< Ref_PixelDigi_ > > MCTruthTTTrackHandle;
 iEvent.getByLabel(TTTracksAssocInputTag,MCTruthTTTrackHandle);

edm::Handle<L1TkElectronParticleCollection> L1TkElectronsHandle;
iEvent.getByLabel(L1TkElectronsInputTag, L1TkElectronsHandle);
std::vector<L1TkElectronParticle>::const_iterator eleIter ;


edm::Handle<L1TkEmParticleCollection> L1TkPhotonsHandle;
iEvent.getByLabel(L1TkPhotonsInputTag, L1TkPhotonsHandle);
std::vector<L1TkEmParticle>::const_iterator phoIter ;
/*
if(GenHandle.isValid()){
//std::cout<<"Gen Part Handle "<<std::endl;
    //pt, eta, phi, z
     
}
*/
tp_eta=999;
tp_phi=999;
tp_pt=-1;
tp_z=999;

tk_eta.resize(0);
tk_phi.resize(0);
tk_pt.resize(0);
tk_z.resize(0);
tk_chi2.resize(0);
tk_truth.resize(0);
mu_eta.resize(0);
mu_phi.resize(0);
mu_pt.resize(0);
muStandAlone_pt.resize(0);
mu_z.resize(0);
muIso.resize(0);

ele_pt.resize(0); 
ele_eta.resize(0);  
ele_phi.resize(0); 
ele_z.resize(0);
eleIso.resize(0);


pho_pt.resize(0);
pho_eta.resize(0);
pho_phi.resize(0);
phoIso.resize(0);
      pho_calopt.resize(0); pho_caloeta.resize(0);  pho_calophi.resize(0);
//std::vector<TTRoad> roads=matcher.makeRoads(iEvent);

if(TrackingParticleHandle.isValid()){
   std::vector< TrackingParticle >::const_iterator iterTP;
for(iterTP=TrackingParticleHandle->begin(); iterTP!=TrackingParticleHandle->end(); ++iterTP){
	if(abs(iterTP->pdgId())!=PDG)continue;
	//pt, eta, phi, z		
	tp_eta=iterTP->eta();
	tp_phi=iterTP->phi();
	tp_pt=iterTP->pt();
	tp_z=iterTP->vz();
	break;	
   }

}

if(TTTrackHandle.isValid()){
std::vector< TTTrack< Ref_PixelDigi_ > >::const_iterator iterL1Track;

//float dRMin=999;
int this_l1track = 0;
for ( iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++ ) {
edm::Ptr< TTTrack< Ref_PixelDigi_ > > l1track_ptr(TTTrackHandle, this_l1track);
++this_l1track;
tk_eta.push_back(iterL1Track->getMomentum().eta());
tk_phi.push_back(iterL1Track->getMomentum().phi());
tk_pt.push_back(iterL1Track->getMomentum().perp());
tk_z.push_back(iterL1Track->getPOCA().z());
tk_chi2.push_back(iterL1Track->getChi2());
int tmp_trk_genuine = 0;

if (MCTruthTTTrackHandle->isGenuine(l1track_ptr)) tmp_trk_genuine = 1;
else tmp_trk_genuine=0;
tk_truth.push_back(tmp_trk_genuine);
}
}
 if ( L1TkMuonsHandle.isValid() ) {
    for (muIter = L1TkMuonsHandle -> begin(); muIter != L1TkMuonsHandle->end(); ++muIter) {

        float pt = muIter -> pt();


        float eta = muIter -> eta();
        float phi = muIter -> phi();
        float zvtx = muIter -> getTrkzVtx();
	mu_pt.push_back(pt);
	muStandAlone_pt.push_back(muIter->getMuExtendedRef()->pt());
	mu_eta.push_back(eta);
	mu_phi.push_back(phi);
        mu_z.push_back(zvtx);
	//compute isolation from Tracks
	float isolation=0;
	for(unsigned t=0; t<tk_z.size(); ++t){
	   float dEta=fabs(tk_eta[t]-eta);
           float dPhi=fabs(tk_phi[t]-phi);
	   float dR=sqrt(dEta*dEta+dPhi*dPhi);
	   float dZ=fabs(zvtx-tk_z[t]);
	   if(dZ>1.0)continue;
	   if(dR>0.3)continue;
	   isolation+=tk_pt[t];
	}
	   muIso.push_back(isolation);	
	    
    }
 }
if ( L1TkElectronsHandle.isValid() ) {
for (eleIter = L1TkElectronsHandle -> begin(); eleIter != L1TkElectronsHandle->end(); ++eleIter) {
	 float et = eleIter -> pt();
	 float phi = eleIter -> phi();
	 float eta = eleIter -> eta();
	 float isolation=eleIter -> getTrkIsol() ;
        ele_pt.push_back(et);
	ele_eta.push_back(eta);
	ele_phi.push_back(phi);
	eleIso.push_back(isolation);

	 const edm::Ref< L1EmParticleCollection > EGref = eleIter -> getEGRef();

	 if ( EGref.isNonnull() ) {

	float et_L1Calo = EGref -> et();
	float eta_calo = EGref -> eta();
	float phi_calo = EGref -> phi();
	   ele_calopt.push_back(et_L1Calo);
           ele_calophi.push_back(phi_calo);
           ele_caloeta.push_back(eta_calo);
	}	
	const edm::Ptr< L1TkTrackType > TrkRef = eleIter -> getTrkPtr();
	if ( TrkRef.isNonnull() ) {
	   float pt_track = TrkRef -> getMomentum().perp();
           float phi_track = TrkRef -> getMomentum().phi();
           float eta_track = TrkRef -> getMomentum().eta();
           float ztrack = TrkRef -> getPOCA().z();
	   ele_tkpt.push_back(pt_track);
           ele_tkphi.push_back(phi_track);
           ele_tketa.push_back(eta_track);
	   ele_z.push_back(ztrack);
	}
     }
  }
if ( L1TkPhotonsHandle.isValid() ) {
for (phoIter = L1TkPhotonsHandle -> begin(); phoIter != L1TkPhotonsHandle->end(); ++phoIter) {
	float et = phoIter -> pt();
	float phi = phoIter -> phi();
	float eta = phoIter -> eta();
	float trkisol = phoIter -> getTrkIsol() ;
	const edm::Ref< L1EmParticleCollection > EGref = phoIter -> getEGRef();
	float et_L1Calo = EGref -> et();
	float eta_calo = EGref -> eta();
	float phi_calo = EGref -> phi();
	pho_pt.push_back(et);
        pho_eta.push_back(eta);
        pho_phi.push_back(phi);
        phoIso.push_back(trkisol);
        pho_calopt.push_back(et_L1Calo);
        pho_caloeta.push_back(eta_calo);
        pho_calophi.push_back(phi_calo);

   }
}

AmMuons->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
TkTriggerParticleAnalzer::beginJob()
{
   edm::Service<TFileService> fs;
    AmMuons = fs->make<TTree>("AmMuons", "");
    AmMuons->Branch("tp_eta", &tp_eta, "tp_eta/F");
    AmMuons->Branch("tp_phi", &tp_phi, "tp_phi/F");
    AmMuons->Branch("tp_pt", &tp_pt, "tp_pt/F");
    AmMuons->Branch("tp_z", &tp_z, "tp_z/F");

    AmMuons->Branch("tk_eta", &tk_eta);
    AmMuons->Branch("tk_phi", &tk_phi);
    AmMuons->Branch("tk_pt", &tk_pt);
    AmMuons->Branch("tk_z", &tk_z);
    AmMuons->Branch("tk_chi2", &tk_chi2);
    AmMuons->Branch("tk_truth", &tk_truth);
    AmMuons->Branch("mu_eta", &mu_eta);
    AmMuons->Branch("mu_phi", &mu_phi);
    AmMuons->Branch("mu_pt", &mu_pt);
    AmMuons->Branch("muStandAlone_pt", &muStandAlone_pt);
    AmMuons->Branch("mu_z", &mu_z);
    AmMuons->Branch("muIso", &muIso);

    AmMuons->Branch("ele_eta", &ele_eta);
    AmMuons->Branch("ele_phi", &ele_phi);
    AmMuons->Branch("ele_pt", &ele_pt);
    AmMuons->Branch("ele_z", &ele_z);
    AmMuons->Branch("eleIso", &eleIso);

    AmMuons->Branch("ele_tketa", &ele_tketa);
    AmMuons->Branch("ele_tkphi", &ele_tkphi);
    AmMuons->Branch("ele_tkpt", &ele_tkpt);

    AmMuons->Branch("ele_caloeta", &ele_caloeta);
    AmMuons->Branch("ele_calophi", &ele_calophi);
    AmMuons->Branch("ele_calopt", &ele_calopt);


    AmMuons->Branch("pho_eta", &pho_eta);
    AmMuons->Branch("pho_phi", &pho_phi);
    AmMuons->Branch("pho_pt", &pho_pt);
    AmMuons->Branch("phoIso", &phoIso);
    AmMuons->Branch("pho_caloeta", &pho_caloeta);
    AmMuons->Branch("pho_calophi", &pho_calophi);
    AmMuons->Branch("pho_calopt", &pho_calopt);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TkTriggerParticleAnalzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
TkTriggerParticleAnalzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
TkTriggerParticleAnalzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
TkTriggerParticleAnalzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
TkTriggerParticleAnalzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TkTriggerParticleAnalzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TkTriggerParticleAnalzer);
