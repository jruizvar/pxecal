// -*- C++ -*-
//
// Package:    Pxecal
// Class:      Pxecal
// 
/**\class Pxecal Pxecal.cc pixtrack/Pxecal/src/Pxecal.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jose Cupertino Ruiz Vargas,32 4-C14,+41227674949,
//         Created:  Wed Aug  7 11:57:33 CEST 2013
// $Id$
//
//


// system include files
#include <cmath>
// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Magnetic field
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

// Geometry
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

// Tracker data formats
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

// --- L1 Egamma dataFormats
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"

// SiPixelRecHit dataFormat
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

// Pileup
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// --- GenParticles
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// --- SimTracks & SimVertex
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

// --- Beam Spot
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// ROOT includes
#include "TTree.h"


class Pxecal : public edm::EDAnalyzer {
   public:
      explicit Pxecal(const edm::ParameterSet&);
      ~Pxecal();

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------Function taken from atani's code ---------------------------
      GlobalPoint getCalorimeterPosition(double phi, double eta, double e);
      
      int                         run;
      int                       event;
      float                        bz;
      int                     bunch_n;
      std::vector<int>         pileup;
      float                beamspot_x;
      float                beamspot_y;
      float           beamspot_widthX;
      float           beamspot_widthY;
      float           beamspot_sigmaZ;
      int                    simtrk_n;
      std::vector<float>    simtrk_pt;
      std::vector<float>   simtrk_eta;
      std::vector<float>   simtrk_phi;
      std::vector<int>      simtrk_id;
      std::vector<int>    simtrk_type;
      float                    gen_vx;
      float                    gen_vy;
      float                    gen_vz;
      int                   genpart_n;
      std::vector<float>    genpart_e;
      std::vector<float>   genpart_et;
      std::vector<float>   genpart_pt;
      std::vector<float>  genpart_eta;
      std::vector<float>  genpart_phi;
      std::vector<int> genpart_charge;
      std::vector<int>     genpart_id;
      int                     genel_n;
      std::vector<float>      genel_e;
      std::vector<float>     genel_et;
      std::vector<float>     genel_pt;
      std::vector<float>    genel_eta;
      std::vector<float>    genel_phi;
      std::vector<int>   genel_charge;
      std::vector<int>       genel_id;
      int                    egamma_n;
      std::vector<float>     egamma_e;
      std::vector<float>    egamma_et;
      std::vector<float>   egamma_eta;
      std::vector<float>   egamma_phi;
      std::vector<float>    egamma_gx;
      std::vector<float>    egamma_gy;
      std::vector<float>    egamma_gz;
      std::vector<int>  egamma_charge;
      int                    recHit_n;
      std::vector<float>    recHit_lx;
      std::vector<float>    recHit_ly;
      std::vector<float>    recHit_gx;
      std::vector<float>    recHit_gy;
      std::vector<float>    recHit_gz;
      std::vector<float>   recHit_rho;
      std::vector<int>   recHit_subid;
      std::vector<int>   recHit_layer;
      std::vector<int>    recHit_disk;
      std::vector<int>  recHit_spread;
      std::vector<int> recHit_spreadx;
      std::vector<int> recHit_spready;
      void InitializeVectors();
      TTree* t;
};

Pxecal::Pxecal(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  t = fs->make<TTree>("t","t");
  t->Branch("Run",            &run);
  t->Branch("Event",          &event);
  t->Branch("BZfield",        &bz);
  t->Branch("BunchN",         &bunch_n);
  t->Branch("Pileup",         &pileup);
  t->Branch("BeamSpotX",      &beamspot_x);
  t->Branch("BeamSpotY",      &beamspot_y);
  t->Branch("BeamSpotWidthX", &beamspot_widthX);
  t->Branch("BeamSpotWidthY", &beamspot_widthY);
  t->Branch("BeamSpotSigmaZ", &beamspot_sigmaZ);
  t->Branch("SimTrkN",        &simtrk_n);
  t->Branch("SimTrkPt",       &simtrk_pt);
  t->Branch("SimTrkEta",      &simtrk_eta);
  t->Branch("SimTrkPhi",      &simtrk_phi);
  t->Branch("SimTrkId",       &simtrk_id);
  t->Branch("SimTrkType",     &simtrk_type);
  t->Branch("GenVx",          &gen_vx);
  t->Branch("GenVy",          &gen_vy);
  t->Branch("GenVz",          &gen_vz);
  t->Branch("GenPartN",       &genpart_n);
  t->Branch("GenPartE",       &genpart_e);
  t->Branch("GenPartEt",      &genpart_et);
  t->Branch("GenPartPt",      &genpart_pt);
  t->Branch("GenPartEta",     &genpart_eta);
  t->Branch("GenPartPhi",     &genpart_phi);
  t->Branch("GenPartCharge",  &genpart_charge);
  t->Branch("GenPartId",      &genpart_id);
  t->Branch("GenElN",         &genel_n);
  t->Branch("GenElE",         &genel_e);
  t->Branch("GenElEt",        &genel_et);
  t->Branch("GenElPt",        &genel_pt);
  t->Branch("GenElEta",       &genel_eta);
  t->Branch("GenElPhi",       &genel_phi);
  t->Branch("GenElCharge",    &genel_charge);
  t->Branch("GenElId",        &genel_id);
  t->Branch("EgN",            &egamma_n);
  t->Branch("EgE",            &egamma_e);
  t->Branch("EgEt",           &egamma_et);
  t->Branch("EgEta",          &egamma_eta);
  t->Branch("EgPhi",          &egamma_phi);
  t->Branch("EgGx",           &egamma_gx);
  t->Branch("EgGy",           &egamma_gy);
  t->Branch("EgGz",           &egamma_gz);
  t->Branch("EgCharge",       &egamma_charge);
  t->Branch("ClN",            &recHit_n);
  t->Branch("ClLx",           &recHit_lx);
  t->Branch("ClLy",           &recHit_ly);
  t->Branch("ClGx",           &recHit_gx);
  t->Branch("ClGy",           &recHit_gy);
  t->Branch("ClGz",           &recHit_gz);
  t->Branch("ClRho",          &recHit_rho);
  t->Branch("ClSubid",        &recHit_subid);
  t->Branch("ClLayer",        &recHit_layer);
  t->Branch("ClDisk",         &recHit_disk);
  t->Branch("ClSize",         &recHit_spread);
  t->Branch("ClSizeX",        &recHit_spreadx);
  t->Branch("ClSizeY",        &recHit_spready);
}

Pxecal::~Pxecal()
{
}

// Taken from atanis' code
GlobalPoint Pxecal::getCalorimeterPosition(double phi, double eta, double e) {
  double x = 0;
  double y = 0;
  double z = 0;
  double depth = 0.89*(7.7+ log(e) );
  double theta = 2*atan(exp(-1*eta));
  double r = 0;
  if( fabs(eta) > 1.479 )
    {
      double ecalZ = 314*fabs(eta)/eta;

      r = ecalZ / cos( 2*atan( exp( -1*eta ) ) ) + depth;
      x = r * cos( phi ) * sin( theta );
      y = r * sin( phi ) * sin( theta );
      z = r * cos( theta );
    }
  else
    {
      double rperp = 129.0;
      double zface =  sqrt( cos( theta ) * cos( theta ) /
                            ( 1 - cos( theta ) * cos( theta ) ) *
                            rperp * rperp ) * abs( eta ) / eta;
      r = sqrt( rperp * rperp + zface * zface ) + depth;
      x = r * cos( phi ) * sin( theta );
      y = r * sin( phi ) * sin( theta );
      z = r * cos( theta );
    }
  GlobalPoint pos(x,y,z);
  return pos;
}

// ------------ method called for each event  ------------
void Pxecal::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  InitializeVectors();
  run = e.id().run();
  event = e.id().event();
  ///////////////////////////////////////////////////////////
  // Magnetic Field  
  //////////////////////////////////////////////////////////
  edm::ESHandle<MagneticField> theMagField;
  es.get<IdealMagneticFieldRecord>().get(theMagField);
  bz = theMagField.product()->inTesla(GlobalPoint(0,0,0)).z();
  //////////////////////////////////////////////////////////
  // Geometry 
  //////////////////////////////////////////////////////////
  edm::ESHandle<TrackerGeometry> geom;
  es.get<TrackerDigiGeometryRecord>().get(geom);
  const TrackerGeometry* theTracker = geom.product();
  //////////////////////////////////////////////////////////
  // Topology 
  //////////////////////////////////////////////////////////
  edm::ESHandle<TrackerTopology> topo;
  es.get<IdealGeometryRecord>().get(topo);
  //////////////////////////////////////////////////////////
  // RecHits
  //////////////////////////////////////////////////////////
  edm::Handle<SiPixelRecHitCollection> recHitColl;
  e.getByLabel( "siPixelRecHits", recHitColl);
  SiPixelRecHitCollection::const_iterator recHitIdIterator    = (recHitColl.product())->begin();
  SiPixelRecHitCollection::const_iterator recHitIdIteratorEnd = (recHitColl.product())->end();
  // Loop over Detector IDs
  recHit_n = 0;
  for ( ; recHitIdIterator != recHitIdIteratorEnd; recHitIdIterator++) {
    SiPixelRecHitCollection::DetSet detset = *recHitIdIterator;
    DetId detId = DetId(detset.detId()); // Get the Detid object
    const GeomDet* geomDet( theTracker->idToDet(detId) );
    int subid = detId.subdetId();
    SiPixelRecHitCollection::DetSet::const_iterator iterRecHit    = detset.begin();
    SiPixelRecHitCollection::DetSet::const_iterator iterRecHitEnd = detset.end();
    // Loop over rechits for this detid
    for ( ; iterRecHit != iterRecHitEnd; ++iterRecHit) {
      LocalPoint lp = iterRecHit->localPosition();
      GlobalPoint GP = geomDet->surface().toGlobal(lp);
      double rho = sqrt(GP.x()*GP.x() + GP.y()*GP.y());
      if(rho < 20){
        recHit_lx.push_back(lp.x());
        recHit_ly.push_back(lp.y());
        recHit_gx.push_back(GP.x());
        recHit_gy.push_back(GP.y());
        recHit_gz.push_back(GP.z());
        recHit_rho.push_back(rho);
        recHit_subid.push_back(subid);
        if(subid==PixelSubdetector::PixelBarrel){
          recHit_layer.push_back(topo->pxbLayer(detId.rawId()));
          recHit_disk.push_back(-99);
        }
        if(subid==PixelSubdetector::PixelEndcap){
          recHit_layer.push_back(-99);
          recHit_disk.push_back(topo->pxfDisk(detId()));
        }
        // Cluster size
        SiPixelRecHit::ClusterRef const& Cluster =  iterRecHit->cluster();
        recHit_spread.push_back(Cluster->size());
        recHit_spreadx.push_back(Cluster->sizeX());
        recHit_spready.push_back(Cluster->sizeY());
        recHit_n ++;
      }
    } // close loop over rechits for this detid
  } // close loop over detector IDs
  ///////////////////////////////////////////////////////////
  // L1 Ecal Trigger Primitives
  //////////////////////////////////////////////////////////
  edm::Handle< l1extra::L1EmParticleCollection > Egamma;
  e.getByLabel( "SLHCL1ExtraParticles","EGamma", Egamma );
  egamma_n = 0;
  for(l1extra::L1EmParticleCollection::const_iterator it = Egamma->begin(); it!=Egamma->end(); ++it){
    egamma_e.push_back(it->energy());
    egamma_et.push_back(it->et());
    egamma_eta.push_back(it->eta());
    egamma_phi.push_back(it->phi());
    egamma_charge.push_back(it->charge());
    GlobalPoint pos= getCalorimeterPosition(it->phi(), it->eta(), it->energy());
    egamma_gx.push_back(pos.x());
    egamma_gy.push_back(pos.y());
    egamma_gz.push_back(pos.z());
    egamma_n++;
  }
  ///////////////////////////////////////////////////////////
  // SimTracks & SimVertices
  //////////////////////////////////////////////////////////
  edm::Handle< edm::SimTrackContainer >   simTrackHandle;
  edm::Handle< edm::SimVertexContainer >  simVertex;
  e.getByLabel( "g4SimHits", simTrackHandle );
  e.getByLabel( "g4SimHits", simVertex );
  edm::SimTrackContainer::const_iterator iterSimTracks;
  simtrk_n = 0;
  for (iterSimTracks = simTrackHandle->begin(); iterSimTracks != simTrackHandle->end(); ++iterSimTracks) {
    simtrk_n++;
    simtrk_pt.push_back( iterSimTracks->momentum().pt() );   
    simtrk_eta.push_back( iterSimTracks->momentum().eta() );   
    simtrk_phi.push_back( iterSimTracks->momentum().phi() );   
    simtrk_id.push_back( iterSimTracks->trackId() );   
    simtrk_type.push_back( iterSimTracks->type() );   
    int index = iterSimTracks->vertIndex();
    // index==0 gets the primary vertex;
    if(index==0){
      gen_vx = simVertex->at(index).position().x();   
      gen_vy = simVertex->at(index).position().y();   
      gen_vz = simVertex->at(index).position().z();   
    }
  }
  ///////////////////////////////////////////////////////////
  // Gen Particles  
  //////////////////////////////////////////////////////////
  edm::Handle< reco::GenParticleCollection > genParticles;
  e.getByLabel( "genParticles", genParticles );
  genpart_n = 0;
  genel_n = 0;
  for(size_t i = 0; i < genParticles->size(); ++i){
    const reco::GenParticle & genParticle = genParticles->at(i);
    genpart_e.push_back( genParticle.energy() );
    genpart_et.push_back( genParticle.et() );
    genpart_pt.push_back( genParticle.pt() );
    genpart_eta.push_back( genParticle.eta() );
    genpart_phi.push_back( genParticle.phi() );
    genpart_charge.push_back( genParticle.charge() );
    genpart_id.push_back( genParticle.pdgId() );
    genpart_n++;
    if ( abs(genParticles->at(i).pdgId()) == 11  ) {
      const reco::GenParticle & genElectron = genParticles->at(i);
      genel_e.push_back( genElectron.energy() );
      genel_et.push_back( genElectron.et() );
      genel_pt.push_back( genElectron.pt() );
      genel_eta.push_back( genElectron.eta() );
      genel_phi.push_back( genElectron.phi() );
      genel_charge.push_back( genElectron.charge() );
      genel_id.push_back( genElectron.pdgId() );
      genel_n++;
    }
  }
  ///////////////////////////////////////////////////////////
  // Beam Spot  
  //////////////////////////////////////////////////////////
  edm::Handle<reco::BeamSpot> thebeamSpot;
  e.getByLabel("BeamSpotFromSim", "BeamSpot", thebeamSpot);
  beamspot_x = thebeamSpot->x0();
  beamspot_y = thebeamSpot->y0();
  beamspot_widthX = thebeamSpot->BeamWidthX();
  beamspot_widthY = thebeamSpot->BeamWidthY();
  beamspot_sigmaZ = thebeamSpot->sigmaZ();
  ///////////////////////////////////////////////////////////
  // Pileup  
  //////////////////////////////////////////////////////////
  edm::Handle< std::vector<PileupSummaryInfo> > puinfo;
  e.getByLabel("addPileupInfo", puinfo);
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  bunch_n=0;
  for(PVI = puinfo->begin(); PVI != puinfo->end(); ++PVI) {
     pileup.push_back(PVI->getPU_NumInteractions());
     bunch_n++;
  }
  // Fill tree
  t->Fill();

}// close Pxecal::analyze


// ------------ method called once each job just before starting event loop  ------------
void 
Pxecal::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Pxecal::endJob() 
{
}

void Pxecal::InitializeVectors()
{
          pileup.clear();
       simtrk_pt.clear();  
      simtrk_eta.clear();
      simtrk_phi.clear();
       simtrk_id.clear();
     simtrk_type.clear();
       genpart_e.clear();
      genpart_et.clear();
      genpart_pt.clear();
     genpart_eta.clear();
     genpart_phi.clear();
  genpart_charge.clear();
      genpart_id.clear();
         genel_e.clear();
        genel_et.clear();
        genel_pt.clear();
       genel_eta.clear();
       genel_phi.clear();
    genel_charge.clear();
        genel_id.clear();
        egamma_e.clear();
       egamma_et.clear();
      egamma_eta.clear();
      egamma_phi.clear();
       egamma_gx.clear();
       egamma_gy.clear();
       egamma_gz.clear();
   egamma_charge.clear();
       recHit_lx.clear();
       recHit_ly.clear();
       recHit_gx.clear();
       recHit_gy.clear();
       recHit_gz.clear();
      recHit_rho.clear();
    recHit_subid.clear();
    recHit_layer.clear();
     recHit_disk.clear();
   recHit_spread.clear();
  recHit_spreadx.clear();
  recHit_spready.clear();
}

//define this as a plug-in
DEFINE_FWK_MODULE(Pxecal);
