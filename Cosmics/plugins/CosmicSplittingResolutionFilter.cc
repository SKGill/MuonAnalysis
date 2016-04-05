#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/DataRecord/interface/SiStripCondDataRecords.h"
#include "CondFormats/SiStripObjects/interface/SiStripLatency.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "CosmicSplittingResolutionNtuple.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimGeneral/TrackingAnalysis/interface/SimHitTPAssociationProducer.h"

////////////////////////////////////////////////////////////////////////////////

track_type remapMuonTrackType(reco::Muon::MuonTrackType t) {
  switch (t) {
  case reco::Muon::InnerTrack: return tk_tkonly;
  case reco::Muon::CombinedTrack: return tk_global;
  case reco::Muon::TPFMS: return tk_tpfms;
  case reco::Muon::Picky: return tk_picky;
  case reco::Muon::DYT: return tk_dyt;
  default:
    throw cms::Exception("remapMuonTrackType") << "t value of " << t << " unexpected in cocktail function.";
  }
}


std::vector<reco::Muon::MuonTrackTypePair> tevOptimized_other(const reco::TrackRef& combinedTrack,
						  const reco::TrackRef& trackerTrack,
						  const reco::TrackRef& tpfmsTrack,
						  const reco::TrackRef& pickyTrack,
						  const reco::TrackRef& dytTrack,
						  const double ptThreshold,
						  const double tune1,
						  const double tune2,
						  double dptcut) {

  // Array for convenience below.
  const reco::Muon::MuonTrackTypePair refit[5] = { 
    make_pair(trackerTrack, reco::Muon::InnerTrack), 
    make_pair(combinedTrack,reco::Muon::CombinedTrack),
    make_pair(tpfmsTrack,   reco::Muon::TPFMS),
    make_pair(pickyTrack,   reco::Muon::Picky),
    make_pair(dytTrack,     reco::Muon::DYT)
  }; 
  
  // Calculate the log(tail probabilities). If there's a problem,
  // signify this with prob == 0. The current problems recognized are:
  // the track being not available, whether the (re)fit failed or it's
  // just not in the event, or if the (re)fit ended up with no valid
  // hits.
  double prob[5] = {0.,0.,0.,0.,0.};
  bool valid[5] = {0,0,0,0,0};

  double dptmin = 1.;

  if (dptcut>0) {  
    for (unsigned int i = 0; i < 4; ++i)
      if (refit[i].first.isNonnull())
        if (refit[i].first->ptError()/refit[i].first->pt()<dptmin) dptmin = refit[i].first->ptError()/refit[i].first->pt();
  
    if (dptmin>dptcut) dptcut = dptmin+0.15;
  }

  for (unsigned int i = 0; i < 4; ++i) 
    if (refit[i].first.isNonnull()){ 
      valid[i] = true;
      if (refit[i].first->numberOfValidHits() && (refit[i].first->ptError()/refit[i].first->pt()<dptcut || dptcut<0)) 
	prob[i] = muon::trackProbability(refit[i].first); 
    }

  
  // Start with picky.
  // Old tuneP
  int chosen = 3;
  // TuneR sigma 
  int chosen_tuneR_s = 3;
  // TuneR prob (no DYT)
  int chosen_tuneR_p = 3;
  // TuneR sigma2 (prob cut)
  int chosen_tuneR_s2 = 3;
  
  // If there's a problem with picky, make the default one of the
  // other tracks. Try TPFMS first, then global, then tracker-only.
  if (prob[3] == 0.) { 

    // split so that passing dptcut<0 recreates EXACTLY the old tuneP behavior
    if (dptcut>0) {
      if      (prob[0] > 0.) chosen = 0;
      else if (prob[2] > 0.) chosen = 2;
      else if (prob[1] > 0.) chosen = 1;
    } else {
      if      (prob[2] > 0.) chosen = 2;
      else if (prob[1] > 0.) chosen = 1;
      else if (prob[0] > 0.) chosen = 0;
    }
  } 
  
  // Now the old tuneP algorithm: switch from picky to tracker-only if the
  // difference, log(tail prob(picky)) - log(tail prob(tracker-only))
  // is greater than a tuned value. Then compare the
  // so-picked track to TPFMS in the same manner using another tuned
  // value.
  if (prob[0] > 0. && prob[3] > 0. && (prob[3] - prob[0]) > tune1)
    chosen = 0;
  if (prob[2] > 0. && (prob[chosen] - prob[2]) > tune2)
    chosen = 2;

  // Sanity checks 
  if (chosen == 3 && !valid[3] ) chosen = 2;
  if (chosen == 2 && !valid[2] ) chosen = 1;
  if (chosen == 1 && !valid[1] ) chosen = 0; 
  ////////////////////////////////////////////
  // TuneR sigma algorithm: TuneR + DYT
  ////////////////////////////////////////////

  dptmin = 1.;

  if (dptcut>0) {  
    for (unsigned int i = 0; i < 5; ++i)
      if (refit[i].first.isNonnull())
        if (refit[i].first->ptError()/refit[i].first->pt()<dptmin) dptmin = refit[i].first->ptError()/refit[i].first->pt();
  
    if (dptmin>dptcut) dptcut = dptmin+0.15;
  }

  for (unsigned int i = 0; i < 5; ++i) 
    if (refit[i].first.isNonnull()){ 
      valid[i] = true;
      if (refit[i].first->numberOfValidHits() && (refit[i].first->ptError()/refit[i].first->pt()<dptcut || dptcut<0)) 
	prob[i] = muon::trackProbability(refit[i].first); 
    }


  if (prob[4]>0. && prob[3]>0.) {
    if(refit[3].first->pt()>0 && refit[4].first->pt()>0 &&
       (refit[4].first->ptError()/refit[4].first->pt()-refit[3].first->ptError()/refit[3].first->pt())<=0.0)
      chosen_tuneR_s=4; // dyt
  }
  if (prob[2] > 0. && (prob[chosen_tuneR_s] - prob[2]) > 24) // Hardcoded tune1 for now
    chosen_tuneR_s = 2;
  if (prob[0] > 0. && (prob[chosen_tuneR_s] - prob[0]) > 35) // Hardcoded tune2 for now
    chosen_tuneR_s = 0;
  if ((refit[chosen_tuneR_s].first->pt() < 150 && prob[0] > 0.) || (trackerTrack->pt() < 150 && prob[0] > 0.) )
    if (prob[chosen_tuneR_s] - prob[0] > 1)
      chosen_tuneR_s = 0;
  ////////////////////////////////////////////
  // TuneR prob algorithm: TuneR NO DYT
  ////////////////////////////////////////////
  if (prob[2]>0. && prob[3]>0.) {
    //if(prob[3] - prob[2]>10)
    if(prob[3] - prob[2]>24)
      chosen_tuneR_p=2; // tpfms
  }
  if (prob[0] > 0. && (prob[chosen_tuneR_p] - prob[0]) > 35)
    chosen_tuneR_p = 0;
  if ((refit[chosen_tuneR_p].first->pt() < 150 && prob[0] > 0.) || (trackerTrack->pt() < 150 && prob[0] > 0.) )
    if (prob[chosen_tuneR_p] - prob[0] > 1)
      chosen_tuneR_p = 0;
  ////////////////////////////////////////////
  // TuneR sigma2 algorithm: TuneR + DYT with prob cut
  ////////////////////////////////////////////
  if (prob[4]>0. && prob[3]>0.) {
    if((prob[3]-prob[4])>23)
      chosen_tuneR_s2=4; // dyt
  }
  if (prob[2] > 0. && (prob[chosen_tuneR_s2] - prob[2]) > 10)
    chosen_tuneR_s2 = 2;
  if (prob[0] > 0. && (prob[chosen_tuneR_s2] - prob[0]) > 20)
    chosen_tuneR_s2 = 0;
        
  if ((refit[chosen_tuneR_s2].first->pt() < 150 && prob[0] > 0.) || (trackerTrack->pt() < 150 && prob[0] > 0.) )
    if (prob[chosen_tuneR_s2] - prob[0] > 1)
      chosen_tuneR_s2 = 0;
  
  std::vector<reco::Muon::MuonTrackTypePair> result = {refit[chosen],refit[chosen_tuneR_s],refit[chosen_tuneR_p],refit[chosen_tuneR_s2]};

  // Done. If pT of the chosen track (or pT of the tracker track) is below the threshold value, return the tracker track.
  if (valid[chosen] && refit[chosen].first->pt() < ptThreshold && prob[0] > 0.) result[0] = make_pair(trackerTrack,reco::Muon::InnerTrack);    
  if (trackerTrack->pt() < ptThreshold && prob[0] > 0.) result[0] = make_pair(trackerTrack,reco::Muon::InnerTrack);  
  
  // Return the chosen track (which can be the global track in
  // very rare cases).
  return result;
}

////////////////////////////////////////////////////////////////////////////////

// The cosmic MC values are not given at the PCA to the IP, but rather
// at their entry point to the CMS cavern. To do comparisons, need to
// propagate the MC values from there to the PCA (or wherever the
// track is defined).
// JMTBAD should move this into a library, make InputTags configurable...

template <typename T>
  static void dump_ref(const edm::Event& event, std::ostream& out, const edm::Ref<T>& ref) {
    out << "ref with product id " << ref.id().id();
    if (ref.id().id() == 0) {
      out << "\n";
      return;
    }
    if (&event) {
      edm::Provenance prov = event.getProvenance(ref.id());
      out << ", branch " << prov.branchName() << " (id " << prov.branchID().id() << "),";
    }
    out << " with index " << ref.index() << "\n";
  }


TrajectoryStateClosestToPoint* propagate_mc_to_point(const edm::Event& event, const edm::EventSetup& setup, const GlobalPoint& pos, std::pair<edm::EDGetToken, edm::EDGetToken> tokens, std::ostream* out=0) {

  if (out) *out << "trying to propagate MC truth track to position " << pos << "\n";

  edm::Handle<TrackingParticleCollection> gen_tracks;
  event.getByToken(tokens.first, gen_tracks);
  edm::Handle<SimHitTPAssociationProducer::SimHitTPAssociationList> simHitsTPAssoc;
  event.getByToken(tokens.second, simHitsTPAssoc);  

  if (out) *out << "found " << gen_tracks->size() << " truth tracks\n";

  // There should be exactly one TrackingParticle with |pdgId| =
  // 13. If there seems to be more than one, give up.
  const TrackingParticle* cosmic_tp = 0;
  int cosmic_tp_index = 0;
  for (TrackingParticleCollection::const_iterator it = gen_tracks->begin(), ite = gen_tracks->end(); it != ite; ++it) {    
    if (abs(it->pdgId()) == 13) {
      if (cosmic_tp != 0) {
	cosmic_tp = 0;
	break;
      }
      cosmic_tp = &*it;
      cosmic_tp_index = it - gen_tracks->begin();
    }
  }

  TrajectoryStateClosestToPoint* ret = 0;

  if (cosmic_tp == 0) {
    if (out) *out << "found more than one truth track with |pdgId| = 13! giving up\n";
  }
  else {
    if (out) *out << "found just one truth track with |pdgId| = 13. trying to find the simhit closest to " << pos << "\n";

    edm::ESHandle<MagneticField> b_field;
    edm::ESHandle<GlobalTrackingGeometry> geometry;
    setup.get<IdealMagneticFieldRecord>().get(b_field);
    setup.get<GlobalTrackingGeometryRecord>().get(geometry); 

    // Find the closest simHit to the point we want to propagate to, and
    // save its position and momentum.
    double min_mag2 = 1e99;
    GlobalPoint vtx;
    GlobalVector mom;

    edm::Ref<TrackingParticleCollection> ref(gen_tracks, cosmic_tp_index);
    std::pair<TrackingParticleRef, TrackPSimHitRef> clusterTPpairWithDummyTP(ref,TrackPSimHitRef());
    auto range = std::equal_range(simHitsTPAssoc->begin(), simHitsTPAssoc->end(), 
				  clusterTPpairWithDummyTP, SimHitTPAssociationProducer::simHitTPAssociationListGreater);

    dump_ref(event, *out, ref);

    for(auto ip = range.first; ip != range.second; ++ip) {
      TrackPSimHitRef TPhit = ip->second;
      dump_ref(event, *out, TPhit);
      if(abs(TPhit->particleType()) != 13)
	continue;
      const GeomDet* det = geometry->idToDet(DetId(TPhit->detUnitId()));
      const LocalVector& lv = TPhit->momentumAtEntry();
      const Local3DPoint& lp = TPhit->localPosition();
      const GlobalPoint pt = det->surface().toGlobal(lp);
      const GlobalVector delta = pt - pos;

      if (out) *out << "  its global position: " << pt << " distance^2 to point wanted: " << delta.mag2() << " current min distance^2: " << min_mag2 << "\n";
      
      if (delta.mag2() < min_mag2) {
	mom = det->surface().toGlobal(lv);
	vtx = pt;
	min_mag2 = delta.mag2();
      }

      if (out) *out << " vtx " << vtx << ", momentum " << mom << ", type " << TPhit->particleType() << "\n";

    }

    /* 
    //No trackPSimHit() member of TrackingParticle in CMSSW_7_3_2

    const std::vector<PSimHit>& simhits = cosmic_tp->trackPSimHit();
    for (std::vector<PSimHit>::const_iterator it = simhits.begin(), ite = simhits.end(); it != ite; ++it) {
      const PSimHit& hit = *it;
      if (out) *out << "  considering simhit in " << hit.detUnitId();

      const GeomDet* det = geometry->idToDet(DetId(hit.detUnitId()));
      const LocalVector& lv = hit.momentumAtEntry();
      const Local3DPoint& lp = hit.localPosition();
      const GlobalPoint pt = det->surface().toGlobal(lp);
      const GlobalVector delta = pt - pos;

      if (out) *out << "  its global position: " << pt << " distance^2 to point wanted: " << delta.mag2() << " current min distance^2: " << min_mag2 << "\n";

      if (delta.mag2() < min_mag2) {
	mom = det->surface().toGlobal(lv);
	vtx = pt;
	min_mag2 = delta.mag2();
      }
    }
    */
    if (out) *out << "making trajectory state with vtx " << vtx << ", momentum " << mom << ", charge " << cosmic_tp->charge() << "\n";
    
    // Propagate the momentum from the simHit's point as close to the
    // requested point as possible.
    FreeTrajectoryState fts(vtx, mom, TrackCharge(cosmic_tp->charge()), b_field.product());
    ret = new TrajectoryStateClosestToPoint(TSCPBuilderNoMaterial()(fts, pos));
  }

  return ret;
}

////////////////////////////////////////////////////////////////////////////////

class CosmicSplittingResolutionFilter : public edm::EDFilter {
public:
  explicit CosmicSplittingResolutionFilter(const edm::ParameterSet&);
  ~CosmicSplittingResolutionFilter() { delete nt; }

private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);

  // Whether running on data or simulation.
  bool is_mc;
  // Dataset identifier, different for data and the different binned
  // simulated samples (values defined in
  // CosmicSplittingResolutionNtuple.h).
  unsigned dataset_id;
  // Whether we're in "split tracks" mode (i.e. starting from 1-leg
  // reconstruction which is then split) or "unsplit tracks" mode
  // (i.e. starting from 2-leg reconstruction already).
  bool use_split_tracks;
  // Whether we're in pp reconstruction mode (i.e. not using the
  // specially-configured cosmics reco).
  bool use_pp_reco;
  // The GenParticles label.
  edm::EDGetToken gen_muon_label;
  // The global muon tracks label.
  edm::EDGetToken split_muon_label;
  // The reference track label (for e.g. binning by pT).
  edm::EDGetToken ref_track_label;
  // The label for the map which gets the ultimate global refit.
  edm::EDGetToken  global_map_label;
  // The label for the map which gets the ultimate TPFMS refit.
  edm::EDGetToken  tpfms_map_label;
  // The label for the map which gets the ultimate picky refit.
  edm::EDGetToken  picky_map_label;
  // The label for the map which gets the ultimate DYT refit.
  edm::EDGetToken  dyt_map_label;
  // The label for the map which gets the ultimate tracker-only refit.
  edm::EDGetToken  trackeronly_map_label;
  // Token for Digis
  edm::EDGetToken digisToken;
  // Token for SimTracks
  edm::EDGetToken simTrackToken;
  // Token for SimVertices
  edm::EDGetToken simVertexToken;
  // Token for TrackingParticles
  edm::EDGetToken trackingParticleToken;
  // Token for SimHit association to Tracking Particles
  edm::EDGetToken simHitTPAssocToken;
  // The "cut" value in the TMR algorithm.
  double tmr_cut;
  // The tune values for Tune P.
  double tunep_pt_threshold, tunep_tune1, tunep_tune2, tunep_dptcut;
  // Parameters controlling the sigma switch.
  double n_sigma_switch;
  double sigma_switch_pt_threshold;
  // Sanity check on the difference in angles between the upper and
  // lower tracks in each pair.
  double max_delta_phi;
  double max_delta_theta;
  // Whether to require at least one pixel hit before keeping the
  // event.
  bool require_pixels;
  // Whether we are not using refits.
  bool no_refits;

  // Some basic distributions, mostly for debugging.
  TH1F* h_track_multiplicity;
  TH1F* h_muon_multiplicity;

  // Keep track of how many events fail the cut sequence at which cut
  // ("error").
  TH1F* errors;
  enum error_codes { error_none, error_mc_trigger, error_collections, error_tracks_num, error_muons_num,
		     error_global_map, error_tpfms_map, error_picky_map, error_dyt_map, error_trackeronly_map,
		     error_track_shared, error_zig_zag, error_tracks_not_unsplit,
		     error_no_pixels, error_delta_phi, error_delta_theta, error_last };

  // The ntuple object, and our tree.
  CosmicSplittingResolutionNtuple* nt;
  TTree* tree;

  // We keep a separate simple ntuple to record the cuts by
  // run/lumi/event number.
  unsigned run, lum, evt;
  int error_code;
  TTree* error_tree;

  // Function to handle dealing with the above when we do cut an
  // event.
  void handle_error(int error) {
    errors->Fill(error);
    if (error != error_none) {
      error_code = error;
      error_tree->Fill();
    }
  }
};

////////////////////////////////////////////////////////////////////////////////

CosmicSplittingResolutionFilter::CosmicSplittingResolutionFilter(const edm::ParameterSet& cfg)
  : is_mc(cfg.getParameter<bool>("is_mc")),
    dataset_id(cfg.getParameter<unsigned>("dataset_id")),
    use_split_tracks(cfg.getParameter<bool>("use_split_tracks")),
    use_pp_reco(cfg.getParameter<bool>("use_pp_reco")),
    gen_muon_label(consumes<reco::GenParticleCollection >(cfg.getParameter<edm::InputTag>("gen_muon_label"))),
    split_muon_label(consumes<reco::MuonCollection >(cfg.getParameter<edm::InputTag>("split_muon_label"))),
    ref_track_label(consumes<reco::TrackCollection >(cfg.getParameter<edm::InputTag>("ref_track_label"))),
    global_map_label ( consumes<reco::TrackToTrackMap >(cfg.getParameter<edm::InputTag>("global_map_label"))),
    tpfms_map_label ( consumes<reco::TrackToTrackMap >(cfg.getParameter<edm::InputTag>("tpfms_map_label"))),
    picky_map_label ( consumes<reco::TrackToTrackMap >(cfg.getParameter<edm::InputTag>("picky_map_label"))),
    dyt_map_label ( consumes<reco::TrackToTrackMap >(cfg.getParameter<edm::InputTag>("dyt_map_label"))),
    trackeronly_map_label ( consumes<reco::TrackToTrackMap >(cfg.getParameter<edm::InputTag>("trackeronly_map_label"))),
    digisToken ( consumes<L1GlobalTriggerReadoutRecord >(edm::InputTag("gtDigis"))),
    simTrackToken ( consumes<std::vector<SimTrack >>(edm::InputTag("g4SimHits"))),
    simVertexToken ( consumes<std::vector<SimVertex >>(edm::InputTag("g4SimHits"))),
    trackingParticleToken ( consumes<std::vector<TrackingParticle >>(edm::InputTag("mix", "MergedTrackTruth"))),
    simHitTPAssocToken ( consumes<SimHitTPAssociationProducer::SimHitTPAssociationList >(edm::InputTag("simHitTPAssocProducer"))),
    tmr_cut(cfg.getParameter<double>("tmr_cut")),
    tunep_pt_threshold(cfg.getParameter<double>("tunep_pt_threshold")),
    tunep_tune1(cfg.getParameter<double>("tunep_tune1")),
    tunep_tune2(cfg.getParameter<double>("tunep_tune2")),
    tunep_dptcut(cfg.getParameter<double>("tunep_dptcut")),
    n_sigma_switch(cfg.getParameter<double>("n_sigma_switch")),
    sigma_switch_pt_threshold(cfg.getParameter<double>("sigma_switch_pt_threshold")),
    max_delta_phi(cfg.getParameter<double>("max_delta_phi")),
    max_delta_theta(cfg.getParameter<double>("max_delta_theta")),
    require_pixels(cfg.getParameter<bool>("require_pixels")),
    no_refits(cfg.getParameter<bool>("no_refits"))
{
  // Book 'em.

  edm::Service<TFileService> fs;
  
  h_track_multiplicity = fs->make<TH1F>("track_multiplicity", "", 10, 0, 10);
  h_muon_multiplicity  = fs->make<TH1F>("muon_multiplicity",  "", 10, 0, 10);

  errors = fs->make<TH1F>("errors", "", error_last, 0, error_last);
  errors->GetXaxis()->SetBinLabel(1 + error_none,                "none");
  errors->GetXaxis()->SetBinLabel(1 + error_mc_trigger,          "mc_trigger");
  errors->GetXaxis()->SetBinLabel(1 + error_collections,         "collections");
  errors->GetXaxis()->SetBinLabel(1 + error_tracks_num,          "tracks_num");
  errors->GetXaxis()->SetBinLabel(1 + error_muons_num,           "muons_num");
  errors->GetXaxis()->SetBinLabel(1 + error_global_map,          "global_map");
  errors->GetXaxis()->SetBinLabel(1 + error_tpfms_map,           "tpfms_map");
  errors->GetXaxis()->SetBinLabel(1 + error_dyt_map,             "dyt_map");
  errors->GetXaxis()->SetBinLabel(1 + error_picky_map,           "picky_map");
  errors->GetXaxis()->SetBinLabel(1 + error_trackeronly_map,     "trackeronly_map");
  errors->GetXaxis()->SetBinLabel(1 + error_track_shared,        "track_shared");
  errors->GetXaxis()->SetBinLabel(1 + error_zig_zag,             "zig_zag");
  errors->GetXaxis()->SetBinLabel(1 + error_tracks_not_unsplit,  "tracks_not_unsplit");
  errors->GetXaxis()->SetBinLabel(1 + error_no_pixels,           "no_pixels");
  errors->GetXaxis()->SetBinLabel(1 + error_delta_phi,           "delta_phi");
  errors->GetXaxis()->SetBinLabel(1 + error_delta_theta,         "delta_theta");

  nt = new CosmicSplittingResolutionNtuple;
  tree = fs->make<TTree>("t", "");
  write_to_tree(tree, nt);

  error_tree = fs->make<TTree>("error_t", "");
  error_tree->Branch("run", &run, "run/i");
  error_tree->Branch("lum", &lum, "lum/i");
  error_tree->Branch("event", &evt, "event/i");
  error_tree->Branch("error_code", &error_code, "error_code/I");
}

////////////////////////////////////////////////////////////////////////////////

bool CosmicSplittingResolutionFilter::filter(edm::Event& event, const edm::EventSetup& setup) {

  //std::cout<<"Event\n";

  // Clear out the ntuple, and store some basic per-event quantities.
  memset(nt, 0, sizeof(CosmicSplittingResolutionNtuple));
  nt->id = dataset_id;
  run = nt->run = event.id().run();
  lum = nt->lumi = event.luminosityBlock();
  evt = nt->event = event.id().event();
  error_code = error_none;

  // Cosmics taken during collisions runs can get on tape through a
  // special trigger, technical #25. Store whether this fired.
  edm::Handle<L1GlobalTriggerReadoutRecord> gtReadoutRecord;
  event.getByToken(digisToken, gtReadoutRecord);
  nt->tt25 = gtReadoutRecord->technicalTriggerWord().at(25);
  // JMTBAD need to add some other trigger selection, especially for MC.

  // In at least a few runs the magnetic field value in the database
  // is wrong. This should be dealt with by throwing out those runs,
  // but also keep track of it here for now.
  edm::ESHandle<MagneticField> magfield;
  setup.get<IdealMagneticFieldRecord>().get(magfield);
  static const GlobalPoint p(0,0,0);
  nt->bzat0 = magfield->inTesla(p).z();

  edm::ESHandle<SiStripLatency> apvlat;
  setup.get<SiStripLatencyRcd>().get(apvlat);
  nt->peak_mode = apvlat->singleReadOutMode();

  std::ostringstream out;

  ////////////////////////////////////////////////////////////////////////////////

  // Get the collections we need, depending on our configuration.

  edm::Handle<reco::GenParticleCollection> gen_muons;
  edm::Handle<std::vector<SimTrack>> sim_tracks;
  edm::Handle<std::vector<SimVertex>> sim_vertices;
  if (is_mc) {
    //event.getByToken(gen_muon_label, gen_muons);
    if (gen_muons.failedToGet()) {
      edm::LogInfo("CosmicSplittingResolution") << "Event failed: is_mc true and could not get the gen_muons collection ";// << gen_muon_label;
      handle_error(error_collections);
      //return false;
    }
    event.getByToken(simTrackToken, sim_tracks);
    event.getByToken(simVertexToken, sim_vertices);

    for (auto st:*sim_tracks) {
      if (abs(st.type()) == 13 && !st.noVertex()){
	edm::LogInfo("CosmicSplittingResolution") << "MOMX=" <<st.momentum().Px();
	edm::LogInfo("CosmicSplittingResolution") << "VX=" <<sim_vertices->at(st.vertIndex()).position().X();
      }
    }
  }
  
  edm::Handle<reco::MuonCollection> split_muons;
  event.getByToken(split_muon_label, split_muons);

  edm::Handle<reco::TrackCollection> ref_tracks;
  event.getByToken(ref_track_label, ref_tracks);

  edm::Handle<reco::TrackToTrackMap> global_map;
  event.getByToken(global_map_label,global_map);

  edm::Handle<reco::TrackToTrackMap> tpfms_map;
  event.getByToken(tpfms_map_label,tpfms_map);

  edm::Handle<reco::TrackToTrackMap> picky_map;
  event.getByToken(picky_map_label,picky_map);

  edm::Handle<reco::TrackToTrackMap> dyt_map;
  event.getByToken(dyt_map_label,dyt_map);

  edm::Handle<reco::TrackToTrackMap> trackeronly_map;
  const bool get_trackeronly_map = !no_refits;
  if (get_trackeronly_map)
    event.getByToken(trackeronly_map_label,trackeronly_map);

  // Check that we can get all the above collections.
  if (split_muons.failedToGet() || ref_tracks.failedToGet() || global_map.failedToGet() || tpfms_map.failedToGet() || picky_map.failedToGet() || dyt_map.failedToGet() || (get_trackeronly_map && trackeronly_map.failedToGet())) {
    edm::LogInfo("CosmicSplittingResolution") << "Event failed: could not get one of the collections: "
			   << " split_muons: " << split_muons.failedToGet()
			   << " ref_tracks: " << ref_tracks.failedToGet() 
			   << " global_map: " << global_map.failedToGet()
			   << " tpfms_map: " << tpfms_map.failedToGet()
			   << " picky_map: " << picky_map.failedToGet()
			   << " dyt_map: " << dyt_map.failedToGet()
			   << " trackeronly_map: " << (get_trackeronly_map && trackeronly_map.failedToGet());
    handle_error(error_collections);
    return false;
  }

  ////////////////////////////////////////////////////////////////////////////////

  // Get the one reference track we're going to use for e.g. binning
  // by pT, and store its values in the ntuple.

  unsigned ntracks = ref_tracks->size();
  h_track_multiplicity->Fill(ntracks);

  // Check that multiplicity of the source tracks is what we want: 1
  // unsplit track for split-tracks mode, 2 tracks otherwise. If not
  // in split-tracks mode, pick the lower-hemipshere track out of the
  // two as the reference track, otherwise the ref track is just the
  // unsplit track.

  if (ntracks != (use_split_tracks ? 1 : 2)) {
    edm::LogInfo("CosmicSplittingResolution") << "Event failed: tracks.size = " << ntracks;
    handle_error(error_tracks_num);
    return false;
  }

  unsigned which_ref = 0;
  if (!use_split_tracks)
    which_ref = ref_tracks->at(0).innerPosition().y() < ref_tracks->at(1).innerPosition().y() ? 0 : 1;

  const reco::Track& ref_track = ref_tracks->at(which_ref);

  nt->ref_charge = ref_track.charge();
  nt->ref_chi2   = ref_track.chi2();
  nt->ref_ndof   = (unsigned char)ref_track.ndof();
  nt->ref_pt     = ref_track.pt();
  nt->ref_theta  = ref_track.theta();
  nt->ref_phi    = ref_track.phi();
  nt->ref_dxy    = ref_track.dxy();
  nt->ref_dz     = ref_track.dz();

  nt->ref_error_pt    = ref_track.ptError();
  nt->ref_error_theta = ref_track.thetaError();
  nt->ref_error_phi   = ref_track.phiError();
  nt->ref_error_dxy   = ref_track.dxyError();
  nt->ref_error_dz    = ref_track.dzError();

  const reco::HitPattern& ref_hp = ref_track.hitPattern();
  nt->ref_pixel_layers = ref_hp.pixelLayersWithMeasurement();
  nt->ref_strip_layers = ref_hp.stripLayersWithMeasurement();
  nt->ref_pixel_hits   = ref_hp.numberOfValidPixelHits();
  nt->ref_strip_hits   = ref_hp.numberOfValidStripHits();
  nt->ref_muon_hits    = ref_hp.numberOfValidMuonHits();
  nt->ref_muon_hits_DT    = ref_hp.numberOfValidMuonDTHits();
  nt->ref_muon_hits_CSC    = ref_hp.numberOfValidMuonCSCHits();
  nt->ref_muon_hits_RPC    = ref_hp.numberOfValidMuonRPCHits();  

  ////////////////////////////////////////////////////////////////////////////////

  // Ensure we have exactly two global muons reconstructed, and get
  // them. Then use them to get all the tracks we consider, from the
  // composed refit maps.

  int num_global_muons = 0;
  const reco::Muon* global_muons[2];
  for (reco::MuonCollection::const_iterator mu = split_muons->begin(); mu != split_muons->end(); ++mu) {
    if (mu->isGlobalMuon()) {
      if (num_global_muons < 2)
	global_muons[num_global_muons] = &*mu;
      ++num_global_muons;
    }
  }

  h_muon_multiplicity->Fill(num_global_muons);

  if (num_global_muons != 2) {
    edm::LogInfo("CosmicSplittingResolution") << "Event failed: global muons found: " << num_global_muons;
    handle_error(error_muons_num);
    return false;
  }

  // Figure out which track is on the upper and which is on the lower
  // half of the detector.
  int which_upper = global_muons[0]->globalTrack()->innerPosition().y() > global_muons[1]->globalTrack()->innerPosition().y() ? 0 : 1;
  edm::LogInfo("CosmicSplittingResolution") << "found 2 tracks, track #" << which_upper << " is on upper\n";

  // Provide uniform access to all the tracks for the two muon
  // halves. Store in the order {track type: global, standalone,
  // tracker-only, ... } x {upper, lower}. These will be replaced
  // below with the values from the composed refit maps.
  reco::TrackRef tracks[n_tracks][2]  = { 
    { global_muons[which_upper]->globalTrack(), global_muons[!which_upper]->globalTrack() },
    { global_muons[which_upper]->outerTrack(),  global_muons[!which_upper]->outerTrack() },
    { global_muons[which_upper]->innerTrack(),  global_muons[!which_upper]->innerTrack() },
    { reco::TrackRef(), reco::TrackRef() }, // Null refs, fill next in a separate step since we have to use the map to get TPFMS tracks.
    { reco::TrackRef(), reco::TrackRef() }, // Ditto for Picky tracks.
    { reco::TrackRef(), reco::TrackRef() }, // Ditto for DYT tracks.
    { reco::TrackRef(), reco::TrackRef() }, // Ditto for TMR tracks.
    { reco::TrackRef(), reco::TrackRef() }, // Ditto for SS tracks.
    { reco::TrackRef(), reco::TrackRef() }, // Ditto for TuneP tracks.
    { reco::TrackRef(), reco::TrackRef() }, // Ditto for old TuneP tracks.
    { reco::TrackRef(), reco::TrackRef() }, // Ditto for TuneR sigma tracks.
    { reco::TrackRef(), reco::TrackRef() }, // Ditto for TuneR prob tracks.
    { reco::TrackRef(), reco::TrackRef() }, // Ditto for TuneR sigma2 tracks.
    { global_muons[which_upper]->globalTrack(), global_muons[!which_upper]->globalTrack() }, // Duplicate the global track here for the one where we propagate to the stand-alone track's position.
    { global_muons[which_upper]->innerTrack(), global_muons[!which_upper]->innerTrack() }    // Ditto for the tracker track.
  };

  // Save a copy of the global TrackRefs so they can be used if we're
  // not using refits.
  reco::TrackRef orig_glb[2] = { tracks[tk_global][0], tracks[tk_global][1] };

  // Tracks are now swapped so that the upper track is at index
  // upper=0. Set which_upper to this so we don't use it later and get
  // confused.
  which_upper = upper;

  // We used the TevMuonProducer to refit the global track, so use the
  // refit tracks for the global level instead.
  for (int j = 0; j < 2; ++j) {
    reco::TrackToTrackMap::const_iterator it = global_map->find(tracks[tk_global][j]);
    if (it == global_map->end()) {
      edm::LogInfo("CosmicSplittingResolution") << "Event failed: muon " << j << "'s global refit not found in map!";
      handle_error(error_global_map);
      return false;
    }
    tracks[tk_global][j] = it->val;
  }

  // If not doing refits, use the "keys" stored above.
  reco::TrackRef tev_keys[2];
  if (no_refits) { 
    tev_keys[0] = orig_glb[0];
    tev_keys[1] = orig_glb[1];
  }
  else {
    tev_keys[0] = tracks[tk_global][0];
    tev_keys[1] = tracks[tk_global][1];
  }

  // TPFMS comes from the last refit Global -- get it now.
  for (int j = 0; j < 2; ++j) {
    reco::TrackToTrackMap::const_iterator it = tpfms_map->find(tev_keys[j]);
    if (it == tpfms_map->end()) {
      edm::LogInfo("CosmicSplittingResolution") << "Event failed: muon " << j << "'s TPFMS fit not found in map!";
      handle_error(error_tpfms_map);
      return false;
    }
    tracks[tk_tpfms][j] = it->val;
  }

  // Picky comes from the last refit Global -- get it now.
  for (int j = 0; j < 2; ++j) {
    reco::TrackToTrackMap::const_iterator it = picky_map->find(tev_keys[j]);
    if (it == picky_map->end()) {
      edm::LogInfo("CosmicSplittingResolution") << "Event failed: muon " << j << "'s Picky fit not found in map!";
      handle_error(error_picky_map);
      return false;
    }
    tracks[tk_picky][j] = it->val;
  }

  // DYT comes from the last refit Global -- get it now.
  for (int j = 0; j < 2; ++j) {
    reco::TrackToTrackMap::const_iterator it = dyt_map->find(tev_keys[j]);
    if (it == dyt_map->end()) {
      edm::LogInfo("CosmicSplittingResolution") << "Event failed: muon " << j << "'s DYT fit not found in map!";
      handle_error(error_dyt_map);
      return false;
    }
    tracks[tk_dyt][j] = it->val;
  }

  // Use refit tracker-only tracks, too, if specified to do so.
  if (get_trackeronly_map) {
    for (int j = 0; j < 2; ++j) {
      reco::TrackToTrackMap::const_iterator it = trackeronly_map->find(tracks[tk_tkonly][j]);
      if (it == trackeronly_map->end()) {
	edm::LogInfo("CosmicSplittingResolution") << "Event failed: muon " << j << "'s tracker-only refit not found in map!";
	handle_error(error_trackeronly_map);
	return false;
      }
      tracks[tk_tkonly][j] = it->val;
    }
  }

  // Now that all the tracks are found, do cocktails using those
  // tracks.
  for (int j = 0; j < 2; ++j) {
    const reco::TrackRef& glb = tracks[tk_global][j];
    const reco::TrackRef& tko = tracks[tk_tkonly][j];
    const reco::TrackRef& fms = tracks[tk_tpfms] [j];
    const reco::TrackRef& pmr = tracks[tk_picky] [j];
    const reco::TrackRef& dyt = tracks[tk_dyt] [j];

    // TMR: use tracker-only if the chi^2,dof tail prob difference
    // between TPFMS and tracker-only is bigger than the cut value.
    reco::Muon::MuonTrackTypePair tmr = muon::TMR(tko, fms, tmr_cut);
    if (tmr.second == reco::Muon::None) { // JMTBAD this should be in the function itself...
      tmr.first = glb;
      tmr.second = reco::Muon::CombinedTrack;
    }
    tracks[tk_tmr][j] = tmr.first;
    nt->choice_tmr[j] = remapMuonTrackType(tmr.second);

    // N-sigma switch: if both global and tracker-only pT are above
    // threshold, use global if |1/pT global - 1/pT tracker-only| < N
    // * sigma_{1/pT} tracker-only, else use tracker-only.
    reco::Muon::MuonTrackTypePair sigsw = muon::sigmaSwitch(glb, tko, n_sigma_switch, sigma_switch_pt_threshold);
    tracks[tk_sigsw][j] = sigsw.first;
    nt->choice_sigsw[j] = remapMuonTrackType(sigsw.second);

    // Tune P (see MuonCocktails.h/cc).
    // Added dyt track ref to the function
    reco::Muon::MuonTrackTypePair tunep = muon::tevOptimized(glb, tko, fms, pmr, dyt, tunep_pt_threshold, tunep_tune1, tunep_tune2, tunep_dptcut);
    std::vector<reco::Muon::MuonTrackTypePair> tunep_other = tevOptimized_other(glb, tko, fms, pmr, dyt, tunep_pt_threshold, tunep_tune1, tunep_tune2, tunep_dptcut);
    // 0 - old tuneP
    // 1 - tuneR sigma (tuneR + DYT)
    // 2 - tuneR prob (tuneR no DYT)
    // 3 - tuneR sigma2 (tuneR + DYT prob cut)

    tracks[tk_tunep][j] = tunep.first;
    nt->choice_tunep[j] = remapMuonTrackType(tunep.second);
    tracks[tk_tunep_old][j] = tunep_other[0].first;
    nt->choice_tunep_old[j] = remapMuonTrackType(tunep_other[0].second);
    tracks[tk_tuner_sigma][j] = tunep_other[1].first;
    nt->choice_tuner_sigma[j] = remapMuonTrackType(tunep_other[1].second);
    tracks[tk_tuner_prob][j] = tunep_other[2].first;
    nt->choice_tuner_prob[j] = remapMuonTrackType(tunep_other[2].second);
    tracks[tk_tuner_sigma2][j] = tunep_other[3].first;
    nt->choice_tuner_sigma2[j] = remapMuonTrackType(tunep_other[3].second);
  }

  ////////////////////////////////////////////////////////////////////////////////

  // Perform sanity checks.

  // Ensure the two halves do not share tracks for any of the
  // reconstructors.
  for (int i = 0; i < n_tracks; ++i) {
    if (tracks[i][0] == tracks[i][1]) {
      edm::LogInfo("CosmicSplittingResolution") << "Event failed: muon 0 and muon 1 share the same " << track_nicks[i] << " track!";
      handle_error(error_track_shared);
      return false;
    }
  }

  // Do all the sets of upper, lower agree in order? Otherwise,
  // something fishy going on. JMTBAD revisit. In the loop, skip
  // global, which we are using to define the order, and also skip
  // stand-alone.
  for (int i = 2; i < n_tracks; ++i) {
    if (tracks[i][0]->innerPosition().y() < tracks[i][1]->innerPosition().y()) {
      edm::LogInfo("CosmicSplittingResolution") << "Event failed: " << track_nicks[i] << " zig zags, #0 is at lower y than #1!";
      handle_error(error_zig_zag);
      return false;
    }
  }

  // If we're using split tracks, check if split tracker tracks came
  // from the one unsplit track. Do so by comparing hit raw detids.
  if (use_split_tracks) {
    std::map<unsigned, int> unsplit_hit_ids;
    int num_unsplit_hits = int(ref_track.recHitsSize());
    int num_split_hits = 0;

    // Store the detids of valid hits of the unsplit track in a map so
    // we can do a two-way check (i.e. check if the hits from the
    // unsplit track all wound up on the split tracks, and if all the
    // hits from the split tracks are from the unsplit one).
    for (int k = 0; k < num_unsplit_hits; ++k) {
      if (ref_track.recHit(k)->isValid()) {
	unsigned id = ref_track.recHit(k)->geographicalId().rawId();
	if (unsplit_hit_ids.find(id) == unsplit_hit_ids.end())
	  unsplit_hit_ids[id] = 1;
	else
	  unsplit_hit_ids[id]++;
      }
    }

    out << "Unsplit map:\n";
    for (std::map<unsigned, int>::const_iterator it = unsplit_hit_ids.begin(); it != unsplit_hit_ids.end(); ++it)
      out << it->first << " -> " << it->second << "\n";

    // Check the split hits.
    int num_unsplit_not_on_split = 0;
    int num_split_not_on_unsplit = 0;
    int num_unsplit_duplicated = 0;
    for (int j = 0; j < 2; ++j) {
      const reco::TrackRef& tk = tracks[tk_tkonly][j];
      int num_hits = int(tk->recHitsSize());
      num_split_hits += num_hits;

      for (int k = 0; k < num_hits; ++k) {
	if (tk->recHit(k)->isValid()) {
	  unsigned id = tk->recHit(k)->geographicalId().rawId();
	  if (unsplit_hit_ids.find(id) == unsplit_hit_ids.end()) {
	    out << "Split track " << tk.index() << " has hit with id " << id << " not on unsplit track!\n";
	    num_split_not_on_unsplit++;
	  }
	  else
	    unsplit_hit_ids[id]--;
	}
      }
    }

    // Check the unsplit hits.
    for (std::map<unsigned, int>::const_iterator it = unsplit_hit_ids.begin(); it != unsplit_hit_ids.end(); ++it) {
      if (it->second > 0) {
	out << "Unsplit track has hit with id " << it->first << " not on either split track!\n";
	num_unsplit_not_on_split++;
      }
      else if (it->second < 0) {
	int duplicated_times = -it->second;
	out << "Unsplit track has hit with id " << it->first << " on split tracks " << duplicated_times << " extra times!\n";
	num_unsplit_duplicated += duplicated_times;
      }
    }

    if (num_unsplit_not_on_split/double(num_unsplit_hits) > 0.2 ||
	num_split_not_on_unsplit/double(num_split_hits)   > 0.2 ||
	num_unsplit_duplicated > 0) {
      edm::LogInfo("CosmicSplittingResolution") << "Event failed: unsplit hits did not satisfactorily find their way to the split tracks:\n"
					    << out.str()
					    << "  num_unsplit_hits: " << num_unsplit_hits << " num_split_hits: " << num_split_hits
					    << " num_unsplit_not_on_split: " << num_unsplit_not_on_split
					    << " num_split_not_on_unsplit: " << num_split_not_on_unsplit
					    << " num_unsplit_duplicated: " << num_unsplit_duplicated;
      handle_error(error_tracks_not_unsplit);
      return false;
    }
  }

  // Don't cut on numbers of hits or other things at this stage,
  // except throw out any events that have absolutely no pixel hits.
  if (require_pixels) {
    for (int i = 0; i < n_tracks; ++i) {
      // Don't check stand-alone, which by definition has no pixel hits.
      if (i == tk_stalone)
	continue;
      for (int j = 0; j < 2; ++j) {
	if (tracks[i][j]->hitPattern().pixelLayersWithMeasurement() < 1) {
	  edm::LogInfo("CosmicSplittingResolution") << "Event failed selection on tracks " << track_nicks[i] << " #" << j << ": no pixel hits!";
	  handle_error(error_no_pixels);
	  return false;
	}
      }
    }
  }

  // Check the delta phi and theta between the upper and lower tracks
  // to guard against trying to compute resolution of two different
  // muons (say one of the legs of each were not reconstructed so that
  // we compare the upper leg of one to the lower leg of another).
  // Don't do this for pp_reco, since here the phi will come out back
  // to back. (JMTBAD use epsilon close to pi in this case.)
  if (!use_pp_reco) {
    for (int i = 0; i < n_tracks; ++i) {
      // Don't check stand-alone (may check delta phi(tracker-only,
      // stand-alone) later).
      if (i == tk_stalone)
	continue;

      double delta;
      // JMTBAD check that this shouldn't be using reco::deltaPhi()
      if ((delta = fabs(tracks[i][0]->phi() - tracks[i][1]->phi())) > max_delta_phi) {
	edm::LogInfo("CosmicSplittingResolution") << "Event failed selection on tracks " << track_nicks[i] << ": delta phi: " << delta;
	handle_error(error_delta_phi);
	return false;
      }
      if ((delta = fabs(tracks[i][0]->theta() - tracks[i][1]->theta())) > max_delta_theta) {
	edm::LogInfo("CosmicSplittingResolution") << "Event failed selection on tracks " << track_nicks[i] << ": delta theta: " << delta;
	handle_error(error_delta_theta);
	return false;
      }
    }
  }

  // Check for the pairs of tracks sharing hits. (In unsplit tracks
  // mode, this should always be 0.) Don't cut, but just save the
  // count. JMTBAD can't check when doing no_refits for standard
  // EventContent since the rechits for tevMuons branches aren't
  // there!
  if (!no_refits) {
    for (int i = 0; i < n_tracks; ++i) {
      nt->shared_hits[i] = 0;
      std::vector<unsigned> hits_seen;
      for (int j = 0; j < 2; ++j) {
	for (int k = 0; k < int(tracks[i][j]->recHitsSize()); ++k) {
	  unsigned id = tracks[i][j]->recHit(k)->geographicalId().rawId();
	  if (j == 0)
	    hits_seen.push_back(id);
	  else if (std::find(hits_seen.begin(), hits_seen.end(), id) != hits_seen.end())
	    nt->shared_hits[i] += 1;
	}
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////

  // Passed selection -- fill the rest of the ntuple.
  out.str("");
  for (int i = 0; i < n_tracks; ++i) {
    for (int j = 0; j < 2; ++j) {
      nt->inner_pos[i][j][0] = tracks[i][j]->innerPosition().x();
      nt->inner_pos[i][j][1] = tracks[i][j]->innerPosition().y();
      nt->inner_pos[i][j][2] = tracks[i][j]->innerPosition().z();

      nt->outer_pos[i][j][0] = tracks[i][j]->outerPosition().x();
      nt->outer_pos[i][j][1] = tracks[i][j]->outerPosition().y();
      nt->outer_pos[i][j][2] = tracks[i][j]->outerPosition().z();

      nt->chi2[i][j] = tracks[i][j]->chi2();
      assert(tracks[i][j]->ndof() < 256);
      nt->ndof[i][j] = (unsigned char)tracks[i][j]->ndof();
      nt->prob[i][j] = muon::trackProbability(tracks[i][j]);

      nt->qoverp      [i][j] = tracks[i][j]->qoverp();
      nt->error_qoverp[i][j] = tracks[i][j]->qoverpError();

      nt->unprop_charge[i][j] = tracks[i][j]->charge();
      nt->unprop_pt    [i][j] = tracks[i][j]->pt();
      nt->unprop_theta [i][j] = tracks[i][j]->theta();
      nt->unprop_phi   [i][j] = tracks[i][j]->phi();
      nt->unprop_dxy   [i][j] = tracks[i][j]->dxy();
      nt->unprop_dz    [i][j] = tracks[i][j]->dz();

      nt->unprop_error_pt   [i][j] = tracks[i][j]->ptError();
      nt->unprop_error_theta[i][j] = tracks[i][j]->thetaError();
      nt->unprop_error_phi  [i][j] = tracks[i][j]->phiError();
      nt->unprop_error_dxy  [i][j] = tracks[i][j]->dxyError();
      nt->unprop_error_dz   [i][j] = tracks[i][j]->dzError();

      const reco::HitPattern& hp = tracks[i][j]->hitPattern();
      nt->pixel_layers[i][j] = hp.pixelLayersWithMeasurement();
      nt->strip_layers[i][j] = hp.stripLayersWithMeasurement();
      nt->pixel_hits  [i][j] = hp.numberOfValidPixelHits();
      nt->strip_hits  [i][j] = hp.numberOfValidStripHits();
      nt->muon_hits   [i][j] = hp.numberOfValidMuonHits();
      nt->muon_hits_DT[i][j] = hp.numberOfValidMuonDTHits();
      nt->muon_hits_CSC[i][j]= hp.numberOfValidMuonCSCHits();
      nt->muon_hits_RPC[i][j]= hp.numberOfValidMuonRPCHits(); 

      out << track_nicks[i] << " #" << j << ": q/pt = " << nt->unprop_charge[i][j] << "/" << nt->unprop_pt[i][j] << ";  chi2/ndf = " << nt->chi2[i][j] << "/" << nt->ndof[i][j] << ";  prob = " << nt->prob[i][j] << "\n";
    }
  }

  edm::LogInfo("CosmicSplittingResolution") << "Cocktail info:\n" << out.str() << "choice_tunep[0] = " << track_nicks[nt->choice_tunep[0]] << "  [1] = " << track_nicks[nt->choice_tunep[1]];
  edm::LogInfo("CosmicSplittingResolution") << "Cocktail info:\n" << out.str() << "choice_tunep_old[0] = " << track_nicks[nt->choice_tunep_old[0]] << "  [1] = " << track_nicks[nt->choice_tunep_old[1]];

  if (is_mc) {
    // TRACKING PARTICLE
    event.getByToken(simTrackToken, sim_tracks);
    std::cout<<"SIMTRACKS SIZE="<<sim_tracks->size()<<std::endl;
    event.getByToken(simVertexToken, sim_vertices);

    edm::Handle<TrackingParticleCollection> gen_tracks;
    event.getByToken(trackingParticleToken, gen_tracks);
    std::cout<<"ASSOCIATED SIMTRACKS SIZE="<<gen_tracks->size()<<std::endl;

    for (auto st:*sim_tracks) {
      bool not_seen = true;
      std::cout<<"type="<<st.type()<<" MOMX=" <<st.momentum().Px()<<" MOMY=" <<st.momentum().Py()<<" MOMZ=" <<st.momentum().Pz()<<" E=" <<st.momentum().E()<<std::endl;
      std::cout<< "VX=" <<sim_vertices->at(st.vertIndex()).position().X()<< "VY=" <<sim_vertices->at(st.vertIndex()).position().Y()<< "VZ=" <<sim_vertices->at(st.vertIndex()).position().Z()<<std::endl;;

      if (abs(st.type()) == 13 && !st.noVertex()){
	if (not_seen) {
	  nt->mc_vertex[0] = sim_vertices->at(st.vertIndex()).position().X();
	  nt->mc_vertex[1] = sim_vertices->at(st.vertIndex()).position().Y();
	  nt->mc_vertex[2] = sim_vertices->at(st.vertIndex()).position().Z();

	  nt->unprop_mc_charge = st.charge();
	  nt->unprop_mc_pt     = st.momentum().Pt();
	  nt->unprop_mc_theta  = st.momentum().Theta();
	  nt->unprop_mc_phi    = st.momentum().Phi();  
	  not_seen = false;
	}
	else
	  break;
      }
    }

  //   // Get the status-1 generated muon (the one propagated already to
  //   // the CMS cavern). JMTBAD true?
  //   size_t i = 0;
  //   for (; i < gen_muons->size(); ++i)
  //     if (gen_muons->at(i).status() == 1)
  // 	break;
  //   const reco::GenParticle& gen_mu(gen_muons->at(i));

    edm::LogInfo("CosmicSplittingResolution") << "unprop. MC vertex: " << nt->mc_vertex[0] << ", " << nt->mc_vertex[1] << ", " << nt->mc_vertex[2] << "  pT: " << nt->unprop_mc_pt << " eta: " << -log(tan(nt->unprop_mc_theta/2)) << " phi: " << nt->unprop_mc_phi;
   }

  // Store whether the track had any hits in CSCs or in DTs, and what
  // the TPFMS first station was.
  nt->hit_csc = tracks[tk_global][0]->hitPattern().cscStationsWithAnyHits() > 0 || tracks[tk_global][1]->hitPattern().cscStationsWithAnyHits() > 0;
  nt->hit_dt  = tracks[tk_global][0]->hitPattern().dtStationsWithAnyHits()  > 0 || tracks[tk_global][1]->hitPattern().dtStationsWithAnyHits()  > 0;
  nt->tpfms_first_station[0] = tracks[tk_global][0]->hitPattern().innermostMuonStationWithHits(0);
  nt->tpfms_first_station[1] = tracks[tk_global][1]->hitPattern().innermostMuonStationWithHits(0);
  edm::LogInfo("CosmicSplittingResolution") << "hit_csc: " << nt->hit_csc << " hit_dt: " << nt->hit_dt << " tpfms_first_station: " << nt->tpfms_first_station[0] << ", " << nt->tpfms_first_station[1];

  // Need transient tracks so we can propagate to the same point.
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
  setup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);

  for (int i = 0; i < n_tracks; ++i) {
    // Propagate the lower track to the upper track's bottom position
    // first. First, need to find that bottom position, since the tracks
    // don't all point the same direction (globalCosmicMuons seem always
    // to point down, and the tracker tracks seem always to point
    // radially outward, but who knows if that will change).
    //
    // JMTBAD am I using the right TransientTrackBuilder?  I.e. is it
    // using the right propagator/etc.?
   
    const math::XYZPoint& pos = (tracks[i][upper]->outerPosition().y() < tracks[i][upper]->innerPosition().y()) ? tracks[i][upper]->outerPosition() : tracks[i][upper]->innerPosition();
    GlobalPoint upper_bottom_pos(pos.x(), pos.y(), pos.z());

    for (int j = 0; j < 2; ++j) {
      nt->prop_valid[i][j] = true; // until proven otherwise...

      reco::TransientTrack ttk;
      TrajectoryStateClosestToPoint tscp;

      try {
	ttk = transientTrackBuilder->build(tracks[i][j]);

	// If this is the track set tk_glb_for_sta, propagate each
	// global track to the innerPosition of the appropriate
	// stand-alone track; otherwise do like we describe
	// above. Propagate to innerPosition so that we can just use
	// the unpropagated stand-alone values in the analysis (may
	// change this later).
	if (i == tk_staglb || i == tk_statko) {
	  const math::XYZPoint& pos = tracks[tk_stalone][j]->innerPosition();
	  GlobalPoint prop_to(pos.x(), pos.y(), pos.z());
	  tscp = ttk.trajectoryStateClosestToPoint(prop_to);
	}
	else
	  tscp = ttk.trajectoryStateClosestToPoint(upper_bottom_pos);

	nt->prop_valid[i][j] = tscp.isValid();
      }
      catch (const cms::Exception& e) {
	nt->prop_valid[i][j] = false;
      }
      
      if (!nt->prop_valid[i][j]) {
	edm::LogInfo("CosmicSplittingResolution") << "propagation for track " << i << ", " << j << " failed, skipping filling it";
	continue;
      }

      nt->prop_pos[i][j][0] = tscp.position().x();
      nt->prop_pos[i][j][1] = tscp.position().y();
      nt->prop_pos[i][j][2] = tscp.position().z();

      nt->charge[i][j] = tscp.charge();
      nt->pt    [i][j] = tscp.momentum().perp();
      nt->theta [i][j] = tscp.momentum().theta();
      nt->phi   [i][j] = tscp.momentum().phi();
      nt->dxy   [i][j] = tscp.perigeeParameters().transverseImpactParameter();
      nt->dz    [i][j] = tscp.perigeeParameters().longitudinalImpactParameter();

      double kappa           = tscp.perigeeParameters().transverseCurvature();
      double kappa_error     = tscp.perigeeError().transverseCurvatureError();
      double Bz              = tscp.theState().parameters().magneticField().inTesla(tscp.position()).z();
      double q_over_pt       = kappa/2.99792458e-3/Bz;
      double q_over_pt_error = q_over_pt/kappa * kappa_error;
      double inv_pt_error    = q_over_pt_error/tscp.charge();

      nt->kappa      [i][j] = kappa;
      nt->error_kappa[i][j] = kappa_error;
      nt->bz         [i][j] = Bz;

      nt->error_pt   [i][j] = 1./q_over_pt/q_over_pt * inv_pt_error;
      nt->error_theta[i][j] = tscp.perigeeError().thetaError();
      nt->error_phi  [i][j] = tscp.perigeeError().phiError();
      nt->error_dxy  [i][j] = tscp.perigeeError().transverseImpactParameterError();
      nt->error_dz   [i][j] = tscp.perigeeError().longitudinalImpactParameterError();

      edm::LogInfo("CosmicSplittingResolution") 
	<< "propagated values for " << (j == upper ? "upper " : "lower ") << " " << track_nicks[i] << " track:"
	<< "\n  charge:    " << int(nt->charge[i][j])
	<< "\n  pt:        " << nt->pt[i][j]    << " +/- " << nt->error_pt[i][j]
	<< "\n  theta:     " << nt->theta[i][j] << " +/- " << nt->error_theta[i][j]
	<< "\n  eta:       " << -log(tan(nt->theta[i][j]/2))
	<< "\n  phi:       " << nt->phi[i][j]   << " +/- " << nt->error_phi[i][j]
	<< "\n  dxy:       " << nt->dxy[i][j]   << " +/- " << nt->error_dxy[i][j]
	<< "\n  dz:        " << nt->dz[i][j]    << " +/- " << nt->error_dz[i][j]
	<< "\n  kappa:     " << kappa           << " +/- " << kappa_error
	<< "\n  q_over_pt: " << q_over_pt       << " +/- " << q_over_pt_error
	<< "\n  Bz:        " << Bz;
    }

    if (is_mc) {
      TrajectoryStateClosestToPoint* mc_tscp = 0;
      try {
	std::pair <edm::EDGetToken, edm::EDGetToken> tokens(trackingParticleToken,simHitTPAssocToken);
	mc_tscp = propagate_mc_to_point(event, setup, upper_bottom_pos, tokens, &std::cout);
      }
      catch (const cms::Exception& e) {
	mc_tscp = 0;
	nt->mc_prop_valid[i] = false;
      }

      if (mc_tscp != 0) {
	nt->mc_prop_valid[i] = true;

	nt->mc_prop_pos[i][0] = mc_tscp->position().x();
	nt->mc_prop_pos[i][1] = mc_tscp->position().y();
	nt->mc_prop_pos[i][2] = mc_tscp->position().z();

	nt->mc_charge[i] = mc_tscp->charge();
	nt->mc_pt    [i] = mc_tscp->momentum().perp();
	nt->mc_theta [i] = mc_tscp->momentum().theta();
	nt->mc_phi   [i] = mc_tscp->momentum().phi();
	nt->mc_dxy   [i] = mc_tscp->perigeeParameters().transverseImpactParameter();
	nt->mc_dz    [i] = mc_tscp->perigeeParameters().longitudinalImpactParameter();

	delete mc_tscp; // Done with it, and we own the pointer, so delete it.

	edm::LogInfo("CosmicSplittingResolution") 
	  << "propagated MC truth for " << track_nicks[i] << " track: "
	  << "\n  vertex: " << nt->mc_prop_pos[i][0] << ", " << nt->mc_prop_pos[i][1] << ", " << nt->mc_prop_pos[i][2] << ", "
	  << "\n  q*pt:   " << nt->mc_charge[i]*nt->mc_pt[i]
	  << "\n  theta:  " << nt->mc_theta[i]
	  << "\n  eta:    " << -log(tan(nt->mc_theta[i]/2))
	  << "\n  phi:    " << nt->mc_phi[i]
	  << "\n  dxy:    " << nt->mc_dxy[i]
	  << "\n  dz:     " << nt->mc_dz[i];
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////

  // Done. Store this event into the ntuple, and since we got here,
  // return true for this event.

  handle_error(error_none);
  tree->Fill();
  return true;
}

DEFINE_FWK_MODULE(CosmicSplittingResolutionFilter);
