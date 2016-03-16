#ifndef NTupleTools_AMTrackProducer_h_

#include <memory>
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/PatternMatcher.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitter.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "SLHCUpgradeSimulations/L1TrackTrigger/interface/StubPtConsistency.h"
#include "SLHCL1TrackTriggerSimulations/NTupleTools/interface/MapTTStubs.h" 
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TriggerTowerMap.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearizedTrackFitter.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitterAlgoLTF.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/CombinationFactory.h"
//#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
//#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include<TFile.h>
#include<TTree.h>
#include<TH2F.h>
#include<map>
class AMTrackProducer : public edm::EDProducer {
  public:
    explicit AMTrackProducer(const edm::ParameterSet&);

  private:
    virtual void beginJob();
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void initOptions(ProgramOption &option, int PBIndex);
double genTrackDistanceLongitudinal(const double &z0, const double &cotTheta, const double &pt, const double &d0,
const int charge, const double &B, const double &R, const double &z);
double genTrackDistanceTransverse(const double &pt, const double &phi0, const double &d0,
                                                 const int charge, const double &B, const double &phi, const double &R);
int PatternBankMap(float eta, float phi);

 virtual void endJob();

ProgramOption option_;
PatternMatcher* matcher_;

std::vector<float>tp_pt, tp_eta, tp_phi;
std::vector<float>tp_vz, tp_chgOverPt;
std::vector<float>matchtp_pt, matchtp_eta, matchtp_phi;
std::vector<float>matchtp_vz, matchtp_chgOverPt;

std::vector<float>trk_pt, trk_eta, trk_phi;
std::vector<float>trk_vz, trk_chgOverPt;
std::vector<float>trk_chi2, trk_cotheta;
std::vector<float>trk_rinv;
std::vector<int>trk_class;
std::vector<int>trk_ghost;
int NRoads;int NStubs;
int NStubs0, NStubs1, NStubs2, NStubs3, NStubs4, NStubs5;
int NgoodStubs0, NgoodStubs1, NgoodStubs2, NgoodStubs3, NgoodStubs4, NgoodStubs5;
std::vector<int> StubsinRoad;
std::vector<int> goodStubsinRoad;
TTree*Amtracks;
edm::InputTag RoadsTag_;
edm::InputTag StubsTag_;
TH2F*TrigTowerMap;
TriggerTowerMap * ttmap_;
TFile*fin;
std::shared_ptr<LinearizedTrackFitter> linearizedTrackFitter_;//(std::make_shared<LinearizedTrackFitter>("/fdata/hepx/store/user/rish/AMSIMULATION/Forked/CMSSW_6_2_0_SLHC25_patch3/src/LinearizedTrackFit/LinearizedTrackFit/python/ConstantsProduction/", true, true));
};
void AMTrackProducer::beginJob(){
   ttmap_ = new TriggerTowerMap();
   ttmap_->read("/fdata/hepx/store/user/rish/AMSIMULATION/CMSSW_6_2_0_SLHC25_patch3/src/SLHCL1TrackTriggerSimulations/AMSimulation/data/"); 
   fin=new TFile("TrigTowerMap.root", "READ");
   TrigTowerMap=(TH2F*)fin->Get("h2TrigMap");
   edm::Service<TFileService> fs;

    linearizedTrackFitter_=(std::make_shared<LinearizedTrackFitter>("/fdata/hepx/store/user/rish/AMSIMULATION/Forked/CMSSW_6_2_0_SLHC25_patch3/src/LinearizedTrackFit/LinearizedTrackFit/python/ConstantsProduction/", true, true));
    Amtracks = fs->make<TTree>("Amtracks", "");
/*
    Amtracks->Branch("tp_pt", &tp_pt);
    Amtracks->Branch("tp_eta", &tp_eta);
    Amtracks->Branch("tp_phi", &tp_phi);
    Amtracks->Branch("tp_vz", &tp_vz);
    Amtracks->Branch("tp_chgOverPt", &tp_chgOverPt);
    Amtracks->Branch("matchtp_pt", &matchtp_pt);
    Amtracks->Branch("matchtp_eta", &matchtp_eta);
    Amtracks->Branch("matchtp_phi", &matchtp_phi);
    Amtracks->Branch("matchtp_vz", &matchtp_vz);
    Amtracks->Branch("matchtp_chgOverPt", &matchtp_chgOverPt);
*/
    Amtracks->Branch("trk_pt", &trk_pt);
    Amtracks->Branch("trk_eta", &trk_eta);
    Amtracks->Branch("trk_phi", &trk_phi);
    Amtracks->Branch("trk_vz", &trk_vz);
    Amtracks->Branch("trk_chi2", &trk_chi2);
    Amtracks->Branch("trk_rinv", &trk_rinv);
    Amtracks->Branch("trk_cotheta", &trk_cotheta);
  /*
    Amtracks->Branch("trk_class", &trk_class);
    Amtracks->Branch("trk_ghost", &trk_ghost);
    Amtracks->Branch("NRoads", &NRoads, "NRoads/I");
    Amtracks->Branch("NStubs", &NStubs, "NStubs/I");
    Amtracks->Branch("NStubs0", &NStubs0, "NStubs0/I");
    Amtracks->Branch("NStubs1", &NStubs1, "NStubs1/I");
    Amtracks->Branch("NStubs2", &NStubs2, "NStubs2/I");
    Amtracks->Branch("NStubs3", &NStubs3, "NStubs3/I");
    Amtracks->Branch("NStubs4", &NStubs4, "NStubs4/I");
    Amtracks->Branch("NStubs5", &NStubs5, "NStubs5/I");
    Amtracks->Branch("NgoodStubs0", &NgoodStubs0, "NgoodStubs0/I");
    Amtracks->Branch("NgoodStubs1", &NgoodStubs1, "NgoodStubs1/I");
    Amtracks->Branch("NgoodStubs2", &NgoodStubs2, "NgoodStubs2/I");
    Amtracks->Branch("NgoodStubs3", &NgoodStubs3, "NgoodStubs3/I");
    Amtracks->Branch("NgoodStubs4", &NgoodStubs4, "NgoodStubs4/I");
    Amtracks->Branch("NgoodStubs5", &NgoodStubs5, "NgoodStubs5/I");
    Amtracks->Branch("StubsinRoad", &StubsinRoad);
    Amtracks->Branch("goodStubsinRoad", &goodStubsinRoad);
*/

}
void AMTrackProducer::endJob(){
delete ttmap_;
delete fin;
}
AMTrackProducer::AMTrackProducer(const edm::ParameterSet& iConfig) {
     StubsTag_ =(iConfig.getParameter<edm::InputTag>("inputTagStub"));
     RoadsTag_=(iConfig.getParameter<edm::InputTag>("RoadsInputTag"));
     produces< std::vector< TTTrack< Ref_PixelDigi_ > > >( "Level1TTTracks" ).setBranchAlias("Level1TTTracks");
}

void AMTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    //std::auto_ptr<std::vector<float> > dummies(new std::vector<float>());
//  typedef TTStubAssociationMap<Ref_PixelDigi_> assocmap_stub;

  std::auto_ptr< std::vector< TTTrack< Ref_PixelDigi_ > > > L1TkTracksForOutput( new std::vector< TTTrack< Ref_PixelDigi_ > > );
	 edm::Handle<edmNew::DetSetVector<TTStub<Ref_PixelDigi_> > > pixelDigiTTStubs;
        iEvent.getByLabel(StubsTag_, pixelDigiTTStubs);

    std::auto_ptr<std::vector< TTTrack<Ref_PixelDigi_> > > TTTrackVector(new std::vector<TTTrack<Ref_PixelDigi_> >());

 edm::Handle< std::vector< TTTrack< Ref_PixelDigi_ > > > TTRoadHandle;
   iEvent.getByLabel( RoadsTag_, TTRoadHandle );
  edm::ESHandle<StackedTrackerGeometry> stackedGeometryHandle;
    iSetup.get<StackedTrackerGeometryRecord>().get(stackedGeometryHandle);
 const StackedTrackerGeometry *theStackedGeometry = stackedGeometryHandle.product();
	trk_pt.resize(0);
        trk_phi.resize(0);
        trk_eta.resize(0);

if (TTRoadHandle->size() > 0 ){
//std::cout<<"ROADS FILLED IN MODULE "<<std::endl;
 unsigned int tkCnt = 0;
 std::vector< TTTrack< Ref_PixelDigi_ > >::const_iterator iterTTTrack;
for ( iterTTTrack = TTRoadHandle->begin();
	  iterTTTrack != TTRoadHandle->end();
	  ++iterTTTrack )
    {
      edm::Ptr< TTTrack< Ref_PixelDigi_ > > tempTrackPtr( TTRoadHandle, tkCnt++ );
      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > > trackStubs = tempTrackPtr->getStubRefs();
       std::vector<double> vars;
       for(unsigned int i=0;i<trackStubs.size();i++){
		edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > tempStubRef = trackStubs.at(i);
		GlobalPoint posStub = theStackedGeometry->findGlobalPosition( &(*tempStubRef) );
		vars.push_back(posStub.phi());
		vars.push_back(posStub.perp());
		vars.push_back(posStub.z());
	
	}
	//do fit for this road
	unsigned comb=trackStubs.size();
	if (comb>6)comb=6;
	double normChi2 = linearizedTrackFitter_->fit(vars, comb);	
	const std::vector<double>& pars = linearizedTrackFitter_->estimatedPars();
	std::cout<<pars[0]<<" "<<pars[1]<<" "<<pars[2]<<std::endl;
	trk_pt.push_back(pars[0]);
        trk_phi.push_back(pars[1]);
        trk_eta.push_back(pars[2]);
	//float cottheta=pars[2];
	//float theta=std::atan(1.0 /cottheta)
	//float eta= -std::log(tan(theta()/2.0));
	float pt=1.0/fabs(pars[0]);
	float px=pt*cos(fabs(pars[1]));
	float py=pt*sin(fabs(pars[1]));
	float pz=pt*pars[2];
	GlobalVector p3(px,py,pz);
	TTTrack<Ref_PixelDigi_> aTrack;
	aTrack.setMomentum(p3,4);
	aTrack.setRInv(0.003 * 3.8 * pars[0],4);
        aTrack.setChi2(normChi2,4);
	GlobalPoint POCA(0, 0,pars[3]);
	aTrack.setPOCA(POCA,4);
   	//if(ttracks[t].chi2Red()<5 && aTrack.getStubRefs().size()>=4 && !ttracks[t].isGhost())
	L1TkTracksForOutput->push_back( aTrack);
   }
   std::cout<<"Tracks "<<L1TkTracksForOutput->size()<<std::endl;
//iEvent.put( L1TkTracksForOutput, "Level1TTTracks");
	 //std::cout<<"chi2 "<<normChi2<<std::endl;
iEvent.put( L1TkTracksForOutput, "Level1TTTracks");
Amtracks->Fill();
    }
}



double AMTrackProducer::genTrackDistanceTransverse(const double &pt, const double &phi0, const double &d0,
                                                 const int charge, const double &B, const double &phi, const double &R) 
{
  //std::cout<<"Transverse input: "<<pt<<", "<<phi0<<", "<<d0<<", "<<charge<<", "<<B<<", "<<phi<<", "<<R<<std::endl;
    double rho = charge*pt/(B*0.003);
    double phiGen = phi0 - asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0)));
    double deltaPhi = (phi - phiGen);
    if (deltaPhi > M_PI) deltaPhi -= M_PI;
    else if (deltaPhi < -M_PI) deltaPhi += M_PI;
  //            //std::cout<<"Transverse output: "<<deltaPhi<<std::endl;
    return deltaPhi;
  }
double AMTrackProducer::genTrackDistanceLongitudinal(const double &z0, const double &cotTheta, const double &pt, const double &d0,
const int charge, const double &B, const double &R, const double &z) {
  //                                                                   //std::cout<<"Longitudinal input: "<<z0<<", "<<cotTheta<<", "<<pt<<", "<<d0<<", "<<charge<<", "<<B<<", "<<R<<", "<<z<<std::endl;
   double rho = charge*pt/(B*0.003);
 double zGen = z0 + charge*rho*cotTheta*acos((pow(rho,2) + pow(d0+rho,2) - pow(R,2))/(2*rho*(rho+d0)));
  //                                                                         //std::cout<<"Longitudinal output: "<<z-zGen<<std::endl;
  return (z - zGen);
  }

int  AMTrackProducer::PatternBankMap(float eta, float phi){
int pb=0;
int ib=TrigTowerMap->FindBin(eta,phi);
if(ib>0)pb=TrigTowerMap->GetBinContent(ib);
return pb;
}

void AMTrackProducer::initOptions(ProgramOption &option, int PBIndex){
//    option.bankfile=TString::Format("/fdata/hepx/store/user/rish/AMSIMULATION/CMSSW_6_2_0_SLHC25_patch3/src/SLHCL1TrackTriggerSimulations/AMSimulation/data/Patterns/patternBank_tt%d_sf1_nz1_pt3_100M.root", PBIndex);
    option.bankfile=TString::Format("/fdata/hepx/store/user/rish/AMSIMULATION/Forked/CMSSW_6_2_0_SLHC25_patch3/src/patternBank_tt%d_sf1_nz4_pt3_5M.root", PBIndex);
    option.verbose=2;
    option.tower=PBIndex;
    option.superstrip="sf1_nz4";
 //   option.maxPatterns=150327;
    option.nLayers=6;
    option.minFrequency=1;
     option.maxStubs=99999;
    option.maxMisses=0;
    option.maxRoads=99999;

    option.algo="LTF";
   option.maxChi2=99999;
    option.minNdof=1;
    option.maxCombs=99999;
    option.maxTracks=99999;
    option.hitBits=0;
    option.view="XYZ";
    //option.matrix="
}
DEFINE_FWK_MODULE(AMTrackProducer);

#endif
