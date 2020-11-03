#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"

// #include "DataFormats/L1TCorrelator/interface/TkMuon.h"
// #include "DataFormats/L1TCorrelator/interface/TkMuonFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"

#include "DataFormats/L1Trigger/interface/Muon.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"
#include <iostream>

class GenMuAnalyzer : public edm::EDAnalyzer {
public:
  explicit GenMuAnalyzer(const edm::ParameterSet&);
  virtual ~GenMuAnalyzer() {};

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  edm::EDGetTokenT< reco::GenParticleCollection >     t_genParticle;
  edm::EDGetTokenT< reco::RecoChargedCandidateCollection >   t_L2Muon;
  edm::EDGetTokenT< reco::RecoChargedCandidateCollection >   t_L2MuonHJ;
  edm::EDGetTokenT< reco::RecoChargedCandidateCollection >   t_L2MuonHJTk;
  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > ttTrackToken_;
  edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > ttTrackMCTruthToken_;
  edm::EDGetTokenT< l1t::MuonBxCollection >                  t_L1Muon_;
  double _pt_min;
  double dRcone;

  TH1F *h_gen_pt;
  TH1F *h_gen_eta;
  TH1F *h_genQ_pt;
  TH1F *h_genQ_eta;
  TH1F *h_gen_matL1_pt;
  TH1F *h_gen_matL1_eta;  
  TH1F *h_gen_matL1Q_pt;
  TH1F *h_gen_matL1Q_eta; 
  TH1F *h_gen22_matL1_pt;
  TH1F *h_gen22_matL1_eta;    
  TH1F *h_gen_matL2_pt;
  TH1F *h_gen_matL2_eta;
  TH1F *h_gen22_matL2_pt;
  TH1F *h_gen22_matL2_eta;
  TH1F *h_gen_matL2HJ_pt;
  TH1F *h_gen_matL2HJ_eta;
  TH1F *h_gen_matL2HJTk_pt;
  TH1F *h_gen_matL2HJTk_eta;
  TH1F *h_gen22_matL2HJ_pt;
  TH1F *h_gen22_matL2HJ_eta;
  TH1F *h_gen22_matL2HJTk_pt;
  TH1F *h_gen22_matL2HJTk_eta;
  TH1F *h_gen_matTT_pt;
  TH1F *h_gen_matTT_eta;  
  TH1F *h_gen_matTTQ_pt;
  TH1F *h_gen_matTTQ_eta; 
  TH1F *h_L2_pt;
  TH1F *h_L2_eta;
  TH1F *h_L2Q_eta;
  TH1F *h_L2HJ_pt;
  TH1F *h_L2HJ_eta; 
  TH1F *h_L2HJQ_eta; 
  TH1F *h_L2HJTk_pt;
  TH1F *h_L2HJTk_eta;
  TH1F *h_L2_matgen_pt;
  TH1F *h_L2_matgen_eta;
  TH1F *h_L2Q_matgen_eta;
  TH1F *h_L2HJ_matgen_pt;
  TH1F *h_L2HJ_matgen_eta; 
  TH1F *h_L2HJQ_matgen_eta; 
  TH1F *h_L2HJTk_matgen_pt;
  TH1F *h_L2HJTk_matgen_eta;      
};


GenMuAnalyzer::GenMuAnalyzer(const edm::ParameterSet& iConfig)
  : t_genParticle( consumes< reco::GenParticleCollection >( iConfig.getParameter<edm::InputTag>("genParticle_src") ) ),
    t_L2Muon( consumes< reco::RecoChargedCandidateCollection >   (iConfig.getUntrackedParameter<edm::InputTag>("L2Muon")) ),
    t_L2MuonHJ( consumes< reco::RecoChargedCandidateCollection >   (iConfig.getUntrackedParameter<edm::InputTag>("L2MuonHJ")) ),
    t_L2MuonHJTk( consumes< reco::RecoChargedCandidateCollection >   (iConfig.getUntrackedParameter<edm::InputTag>("L2MuonHJTk")) ),
    ttTrackToken_        ( consumes<std::vector<TTTrack<Ref_Phase2TrackerDigi_> > >(iConfig.getParameter<edm::InputTag>("L1TrackInputTag"     )) ),
    ttTrackMCTruthToken_ ( consumes< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >(iConfig.getParameter<edm::InputTag>("MCTruthTrackInputTag"))),
    t_L1Muon_            ( consumes< l1t::MuonBxCollection  >                 (iConfig.getUntrackedParameter<edm::InputTag>("L1Muon"            )) ),

    _pt_min(       iConfig.getParameter< double >("pt_min")),
    dRcone(        iConfig.getParameter< double >("dRcone_"))
{
  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2(true);

  h_gen_pt              = fs->make<TH1F>("h_gen_pt",  "", 1000, 0, 1000);
  h_gen_eta             = fs->make<TH1F>("h_gen_eta", "", 60, -3, 3);
  h_genQ_pt              = fs->make<TH1F>("h_genQ_pt",  "", 1000, 0, 1000);
  h_genQ_eta             = fs->make<TH1F>("h_genQ_eta", "", 60, -3, 3);
  h_gen_matL1_pt  = fs->make<TH1F>("h_gen_matL1_pt",  "", 1000, 0, 1000);
  h_gen_matL1_eta = fs->make<TH1F>("h_gen_matL1_eta", "", 60, -3, 3);  
  h_gen_matL1Q_pt  = fs->make<TH1F>("h_gen_matL1Q_pt",  "", 1000, 0, 1000);
  h_gen_matL1Q_eta = fs->make<TH1F>("h_gen_matL1Q_eta", "", 60, -3, 3);  
  h_gen22_matL1_pt  = fs->make<TH1F>("h_gen22_matL1_pt",  "", 1000, 0, 1000);
  h_gen22_matL1_eta = fs->make<TH1F>("h_gen22_matL1_eta", "", 60, -3, 3);  
  h_gen_matL2_pt  = fs->make<TH1F>("h_gen_matL2_pt",  "", 1000, 0, 1000);
  h_gen_matL2_eta = fs->make<TH1F>("h_gen_matL2_eta", "", 60, -3, 3);
  h_gen22_matL2_pt  = fs->make<TH1F>("h_gen22_matL2_pt",  "", 1000, 0, 1000);
  h_gen22_matL2_eta = fs->make<TH1F>("h_gen22_matL2_eta", "", 60, -3, 3);
  h_gen_matL2HJ_pt  = fs->make<TH1F>("h_gen_matL2HJ_pt",  "", 1000, 0, 1000);
  h_gen_matL2HJ_eta = fs->make<TH1F>("h_gen_matL2HJ_eta", "", 60, -3, 3);
  h_gen_matL2HJTk_pt  = fs->make<TH1F>("h_gen_matL2HJTk_pt",  "", 1000, 0, 1000);
  h_gen_matL2HJTk_eta = fs->make<TH1F>("h_gen_matL2HJTk_eta", "", 60, -3, 3);
  h_gen22_matL2HJ_pt  = fs->make<TH1F>("h_gen22_matL2HJ_pt",  "", 1000, 0, 1000);
  h_gen22_matL2HJ_eta = fs->make<TH1F>("h_gen22_matL2HJ_eta", "", 60, -3, 3);
  h_gen22_matL2HJTk_pt  = fs->make<TH1F>("h_gen22_matL2HJTk_pt",  "", 1000, 0, 1000);
  h_gen22_matL2HJTk_eta = fs->make<TH1F>("h_gen22_matL2HJTk_eta", "", 60, -3, 3);
  h_gen_matTT_pt  = fs->make<TH1F>("h_gen_matTT_pt",  "", 1000, 0, 1000);
  h_gen_matTT_eta = fs->make<TH1F>("h_gen_matTT_eta", "", 60, -3, 3); 
  h_gen_matTTQ_pt  = fs->make<TH1F>("h_gen_matTTQ_pt",  "", 1000, 0, 1000);
  h_gen_matTTQ_eta = fs->make<TH1F>("h_gen_matTTQ_eta", "", 60, -3, 3); 

  h_L2_pt  = fs->make<TH1F>("h_L2_pt",  "", 1000, 0, 1000);
  h_L2_eta = fs->make<TH1F>("h_L2_eta",  "", 60, -3, 3);
  h_L2Q_eta = fs->make<TH1F>("h_L2Q_eta",  "", 60, -3, 3);
  h_L2HJ_pt  = fs->make<TH1F>("h_L2HJ_pt",  "", 1000, 0, 1000);
  h_L2HJ_eta = fs->make<TH1F>("h_L2HJ_eta",  "", 60, -3, 3);  
  h_L2HJQ_eta = fs->make<TH1F>("h_L2HJQ_eta",  "", 60, -3, 3);  
  h_L2HJTk_pt  = fs->make<TH1F>("h_L2HJTk_pt",  "", 1000, 0, 1000);
  h_L2HJTk_eta = fs->make<TH1F>("h_L2HJTk_eta",  "", 60, -3, 3); 
  h_L2_matgen_pt  = fs->make<TH1F>("h_L2_matgen_pt",  "", 1000, 0, 1000);
  h_L2_matgen_eta = fs->make<TH1F>("h_L2_matgen_eta",  "", 60, -3, 3);
  h_L2Q_matgen_eta = fs->make<TH1F>("h_L2Q_matgen_eta",  "", 60, -3, 3);
  h_L2HJ_matgen_pt  = fs->make<TH1F>("h_L2HJ_matgen_pt",  "", 1000, 0, 1000);
  h_L2HJ_matgen_eta = fs->make<TH1F>("h_L2HJ_matgen_eta",  "", 60, -3, 3);  
  h_L2HJQ_matgen_eta = fs->make<TH1F>("h_L2HJQ_matgen_eta",  "", 60, -3, 3);  
  h_L2HJTk_matgen_pt  = fs->make<TH1F>("h_L2HJTk_matgen_pt",  "", 1000, 0, 1000);
  h_L2HJTk_matgen_eta = fs->make<TH1F>("h_L2HJTk_matgen_eta",  "", 60, -3, 3); 
}

void GenMuAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup&) {

  // double dRcone = 0.3;

  edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon;
  bool hasL2 = iEvent.getByToken( t_L2Muon, h_L2Muon ); 
  edm::Handle<reco::RecoChargedCandidateCollection> h_L2MuonHJ;
  bool hasL2HJ = iEvent.getByToken( t_L2MuonHJ, h_L2MuonHJ );
  edm::Handle<reco::RecoChargedCandidateCollection> h_L2MuonHJTk;
  bool hasL2HJTk = iEvent.getByToken( t_L2MuonHJTk, h_L2MuonHJTk );
  edm::Handle<l1t::MuonBxCollection> h_L1Muon;
  bool hasL1 = iEvent.getByToken(t_L1Muon_, h_L1Muon);

  edm::Handle<reco::GenParticleCollection> h_genParticle;
  if(hasL2HJTk&& hasL2HJ&&hasL2 && iEvent.getByToken(t_genParticle,h_genParticle) ) {
  // if(hasL1&&hasL2HJTk&& hasL2HJ&&hasL2 && iEvent.getByToken(t_genParticle,h_genParticle) ) {
    reco::GenParticleCollection::const_iterator genp = h_genParticle->begin();
    for (; genp != h_genParticle->end(); genp++) {

      if( fabs(genp->pdgId()) != 13 )  continue;
      if( !genp->isPromptFinalState() )  continue;
      if( !genp->fromHardProcessFinalState() )  continue;
      if( fabs(genp->eta()) > 2.4 )  continue;

      h_gen_pt->Fill( genp->pt() );
      if( genp->pt() > _pt_min )  h_gen_eta->Fill( genp->eta() );

      if(genp->pt()>24){
        h_genQ_pt->Fill( genp->pt() );
        if( genp->pt() > _pt_min )  h_genQ_eta->Fill( genp->eta() );

      }

      bool L2matched = false;
      for(auto itL2 = h_L2Muon->begin(); itL2 != h_L2Muon->end(); itL2++) {
        // h_n_L2_pt->Fill( itL2->pt() );
        // h_n_L2_eta->Fill( itL2->eta() );

        if( reco::deltaR( itL2->momentum(), *genp ) < dRcone ) {
          L2matched = true;
          break;
        }
      }
      if( L2matched ) {
        h_gen_matL2_pt->Fill( genp->pt() );
        if( genp->pt() > _pt_min )  h_gen_matL2_eta->Fill( genp->eta() );
      }
      if( L2matched && genp->pt() >24 ) {
        h_gen22_matL2_pt->Fill( genp->pt() );
        if( genp->pt() > _pt_min )  h_gen22_matL2_eta->Fill( genp->eta() );
      }
      bool L2matchedHJ = false;
      for(auto itL2 = h_L2MuonHJ->begin(); itL2 != h_L2MuonHJ->end(); itL2++) {
        // h_n_L2HJ_pt->Fill( itL2->pt() );
        // h_n_L2HJ_eta->Fill( itL2->eta() );

        if( reco::deltaR( itL2->momentum(), *genp ) < dRcone ) {
          L2matchedHJ = true;
          break;
        }
      }
      if( L2matchedHJ ) {
        h_gen_matL2HJ_pt->Fill( genp->pt() );
        if( genp->pt() > _pt_min )  h_gen_matL2HJ_eta->Fill( genp->eta() );
      }
      if( L2matchedHJ &&genp->pt() >24) {
        h_gen22_matL2HJ_pt->Fill( genp->pt() );
        if( genp->pt() > _pt_min )  h_gen22_matL2HJ_eta->Fill( genp->eta() );
      }
      bool L2matchedHJTk = false;
      for(auto itL2 = h_L2MuonHJTk->begin(); itL2 != h_L2MuonHJTk->end(); itL2++) {
        // h_n_L2HJTk_pt->Fill( itL2->pt() );
        // h_n_L2HJTk_eta->Fill( itL2->eta() );

        if( reco::deltaR( itL2->momentum(), *genp ) < dRcone ) {
          L2matchedHJTk = true;
          break;
        }
      }
      if( L2matchedHJTk ) {
        h_gen_matL2HJTk_pt->Fill( genp->pt() );
        if( genp->pt() > _pt_min )  h_gen_matL2HJTk_eta->Fill( genp->eta() );
      }
      if( L2matchedHJTk && genp->pt() >24) {
        h_gen22_matL2HJTk_pt->Fill( genp->pt() );
        if( genp->pt() > _pt_min )  h_gen22_matL2HJTk_eta->Fill( genp->eta() );
      }
      bool L1genmatched = false;
      bool L1Qgenmatched = false;
      for(int ibx = h_L1Muon->getFirstBX(); ibx<=h_L1Muon->getLastBX(); ++ibx){
      if(ibx != 0) continue; // -- only take when ibx == 0 -- //
      for(auto it=h_L1Muon->begin(ibx); it!=h_L1Muon->end(ibx); it++){
        if( reco::deltaR( it->momentum(), *genp ) < dRcone ) {
          L1genmatched = true;
          if(L1genmatched&&it->hwQual()>7&&it->pt()>22) L1Qgenmatched = true;
          break;
          }

        }
      }
      if( L1genmatched ) {
        h_gen_matL1_pt->Fill( genp->pt() );
        if( genp->pt() > _pt_min )  h_gen_matL1_eta->Fill( genp->eta() );
      }  
      if( L1genmatched && genp->pt()>24) {
        h_gen22_matL1_pt->Fill( genp->pt() );
        if( genp->pt() > _pt_min )  h_gen22_matL1_eta->Fill( genp->eta() );
      }           
      if( L1Qgenmatched ) {
        h_gen_matL1Q_pt->Fill( genp->pt() );
        if( genp->pt() > _pt_min )  h_gen_matL1Q_eta->Fill( genp->eta() );
      }      

    }//genparticle loop

    for(auto itL2 = h_L2Muon->begin(); itL2 != h_L2Muon->end(); itL2++) {
      
      h_L2_pt->Fill( itL2->pt() );
      h_L2_eta->Fill( itL2->eta() );
      if(itL2->pt()>24) h_L2Q_eta->Fill( itL2->eta() );

      bool L2_genmatched = false;
      reco::GenParticleCollection::const_iterator genp = h_genParticle->begin();
      for (; genp != h_genParticle->end(); genp++) {

        if( fabs(genp->pdgId()) != 13 )  continue;
        if( !genp->isPromptFinalState() )  continue;
        if( !genp->fromHardProcessFinalState() )  continue;
        if( fabs(genp->eta()) > 2.4 )  continue;
      
        if( reco::deltaR( itL2->momentum(), *genp ) < dRcone ) {
          L2_genmatched = true;
          break;
          }
        }
      if( L2_genmatched ) {
        h_L2_matgen_pt->Fill( itL2->pt() );
        h_L2_matgen_eta->Fill( itL2->eta() );
        if(itL2->pt()>24) h_L2Q_matgen_eta->Fill( itL2->eta() );
      }

    }// L2 loop  

    for(auto itL2 = h_L2MuonHJ->begin(); itL2 != h_L2MuonHJ->end(); itL2++) {
      
      h_L2HJ_pt->Fill( itL2->pt() );
      h_L2HJ_eta->Fill( itL2->eta() );
      if(itL2->pt()>24) h_L2HJQ_eta->Fill( itL2->eta() );

      bool L2HJ_genmatched = false;
      reco::GenParticleCollection::const_iterator genp = h_genParticle->begin();
      for (; genp != h_genParticle->end(); genp++) {

        if( fabs(genp->pdgId()) != 13 )  continue;
        if( !genp->isPromptFinalState() )  continue;
        if( !genp->fromHardProcessFinalState() )  continue;
        if( fabs(genp->eta()) > 2.4 )  continue;
      
        if( reco::deltaR( itL2->momentum(), *genp ) < dRcone ) {
          L2HJ_genmatched = true;
          break;
          }
        }
      if( L2HJ_genmatched ) {
        h_L2HJ_matgen_pt->Fill( itL2->pt() );
        h_L2HJ_matgen_eta->Fill( itL2->eta() );
        if(itL2->pt()>24) h_L2HJQ_matgen_eta->Fill( itL2->eta() );
      }

    }// L2HJ loop  

    for(auto itL2 = h_L2MuonHJTk->begin(); itL2 != h_L2MuonHJTk->end(); itL2++) {
      
      h_L2HJTk_pt->Fill( itL2->pt() );
      h_L2HJTk_eta->Fill( itL2->eta() );

      bool L2HJTk_genmatched = false;
      reco::GenParticleCollection::const_iterator genp = h_genParticle->begin();
      for (; genp != h_genParticle->end(); genp++) {

        if( fabs(genp->pdgId()) != 13 )  continue;
        if( !genp->isPromptFinalState() )  continue;
        if( !genp->fromHardProcessFinalState() )  continue;
        if( fabs(genp->eta()) > 2.4 )  continue;
      
        if( reco::deltaR( itL2->momentum(), *genp ) < dRcone ) {
          L2HJTk_genmatched = true;
          break;
          }
        }
      if( L2HJTk_genmatched ) {
        h_L2HJTk_matgen_pt->Fill( itL2->pt() );
        h_L2HJTk_matgen_eta->Fill( itL2->eta() );
      }

    }// L2HJTk loop  
  } //if(hasL2HJTk&& hasL2HJ&&hasL2 && iEvent.getByToken(t_genParticle,h_genParticle) )

  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTTrackHandle;
  iEvent.getByToken(ttTrackToken_, TTTrackHandle);

  edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTTrackHandle;
  iEvent.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);

  int this_l1track = 0;
  std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator iterL1Track;
  for ( iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++ ) {    
    edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > l1track_ptr(TTTrackHandle, this_l1track);
    this_l1track++;

  edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(l1track_ptr);
  if (!my_tp.isNull()&&fabs(my_tp->pdgId())==13){
    // std::cout << "tp pdgid: "<< my_tp->pdgId() << std::endl;
    reco::GenParticleCollection::const_iterator genp = h_genParticle->begin();
    for (; genp != h_genParticle->end(); genp++) {  
      if( fabs(genp->pdgId()) != 13 )  continue;
      if( !genp->isPromptFinalState() )  continue;
      if( !genp->fromHardProcessFinalState() )  continue;
      if( fabs(genp->eta()) > 2.4 )  continue; 

      // std::cout << "gen eta: "<< genp->eta() << std::endl;
      // std::cout << "tp eta: "<< my_tp->eta() << std::endl;


      bool L1TTgenmatched = false;
      // std::cout << reco::deltaR( *my_tp, *genp ) << std::endl;
      if( reco::deltaR( *my_tp, *genp ) < dRcone ){
        L1TTgenmatched = true;
        if(L1TTgenmatched){
          h_gen_matTT_pt->Fill( genp->pt() );
          if( genp->pt() > _pt_min )  h_gen_matTT_eta->Fill( genp->eta() );
          if(genp->pt() > 24){
              h_gen_matTTQ_pt->Fill( genp->pt() );
              h_gen_matTTQ_eta->Fill( genp->eta() );
            }
          }
        break;
        }
       
      }  //gen loop        
    }
  }// l1 track loop
}

DEFINE_FWK_MODULE(GenMuAnalyzer);
