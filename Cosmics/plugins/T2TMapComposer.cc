#include <map>
#include <vector>
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// JMTBAD move these three to a library.

template <typename T>
void dump_ref(std::ostream& out, const edm::Ref<T>& ref, const edm::Event* event=0) {
  out << "ref with product id " << ref.id().id();
  if (ref.id().id() == 0) {
    out << "\n";
    return;
  }
  if (event) {
    edm::Provenance prov = event->getProvenance(ref.id());
    out << ", branch " << prov.branchName() << " (id " << prov.branchID().id() << "),";
  }
  out << " with index " << ref.index() << "\n";
}

template <typename T>
void dump_handle(std::ostream& out, const edm::Handle<T>& h) {
  out << "handle with product id " << h.id().id();
  if (h.id().id() == 0) {
    out << "\n";
    return;
  }

  const edm::Provenance* prov = h.provenance();
  out << ", branch " << prov->branchName() << " (id " << prov->branchID().id() << ")\n";
}

void dump_t2tmap(std::ostream& out, const reco::TrackToTrackMap& map, const edm::Event* event) {
  out << "map size: " << map.size() << "\n";
  int ipair = 0; // subtracting iterators below doesn't work for some reason
  for (reco::TrackToTrackMap::const_iterator it = map.begin(); it != map.end(); ++it, ++ipair) {
    out << "pair #" << ipair << ":\n  key ref: ";
    dump_ref(out, it->key, event);
    out << "  val ref: ";
    dump_ref(out, it->val, event);
  }
}

class T2TMapComposer : public edm::EDProducer {
public:
  explicit T2TMapComposer(const edm::ParameterSet&);

private:
  typedef std::vector<edm::Handle<reco::TrackToTrackMap> > map_vector;
  reco::TrackToTrackMap* compose(edm::Event&, edm::Handle<reco::TrackCollection>& first_tracks, edm::Handle<reco::TrackCollection>& last_tracks, const map_vector&);
  virtual void produce(edm::Event&, const edm::EventSetup&);

  const bool debug;
  std::vector<std::string> new_map_names;
  std::vector<edm::InputTag> first_track_tags;
  std::vector<edm::InputTag> last_track_tags;
  std::vector<std::vector<edm::InputTag> > map_tags;
};

T2TMapComposer::T2TMapComposer(const edm::ParameterSet& cfg)
  : debug(true),
    first_track_tags(cfg.getParameter<std::vector<edm::InputTag> >("first_track_tags")),
    last_track_tags(cfg.getParameter<std::vector<edm::InputTag> >("last_track_tags"))
{
  if (cfg.existsAs<std::vector<std::string> >("new_map_names")) {
    new_map_names = cfg.getParameter<std::vector<std::string> >("new_map_names");
    for (std::vector<std::string>::const_iterator i = new_map_names.begin(), e = new_map_names.end(); i != e; ++i) {
      map_tags.push_back(cfg.getParameter<std::vector<edm::InputTag> >(*i));
      produces<reco::TrackToTrackMap>(*i);
    }
  }
  else {
    new_map_names.push_back("");
    map_tags.push_back(cfg.getParameter<std::vector<edm::InputTag> >("map_tags"));
    produces<reco::TrackToTrackMap>();
  }
}

reco::TrackToTrackMap* T2TMapComposer::compose(edm::Event& event, edm::Handle<reco::TrackCollection>& first_tracks, edm::Handle<reco::TrackCollection>& last_tracks, const map_vector& maps) {
  std::ostringstream out;
  reco::TrackToTrackMap* t2tmap = new reco::TrackToTrackMap(first_tracks, last_tracks);

  if (debug) {
    out << "first tracks handle: ";
    dump_handle(out, first_tracks);
    out << "last  tracks handle: ";
    dump_handle(out, last_tracks);
  }

  // For each key->value in the first map, follow the mappings through
  // to the end of the maps vector and store key->ultimate_value in
  // the output map.
  for (reco::TrackToTrackMap::const_iterator orig = maps.at(0)->begin(); orig != maps.at(0)->end(); ++orig) {
    if (debug) out << "orig key pt " << orig->key->pt() << " eta " << orig->key->eta() << " phi " << orig->key->phi() << "\n";
    reco::TrackRef curr = orig->val;
    if (debug) out << "first curr pt " << curr->pt() << " eta " << curr->eta() << " phi " << curr->phi() << "\n";

    for (map_vector::const_iterator map = maps.begin() + 1; map != maps.end(); ++map) {
      reco::TrackToTrackMap::const_iterator it = (*map)->find(curr);
      if (it == (*map)->end()) {
	curr = reco::TrackRef();
	break;
      }
      else {
	curr = it->val;
        if (debug) out << "curr now pt " << curr->pt() << " eta " << curr->eta() << " phi " << curr->phi() << "\n";
      }
    }

    // If the link is broken somewhere, don't store anything for this
    // key.
    if (!curr.isNull()) {
      if (debug) {
        out << "before insert:\n"
                  << "orig key pt " << orig->key->pt() << " eta " << orig->key->eta() << " phi " << orig->key->phi() << "\n";
        dump_ref(out, orig->key, &event);
        out << "curr final pt " << curr->pt() << " eta " << curr->eta() << " phi " << curr->phi() << "\n";
        dump_ref(out, curr, &event);
      }

      t2tmap->insert(orig->key, curr);
    }
  }

  if (debug) edm::LogInfo("T2TMapComposer") << out.str();

  return t2tmap;
}

void T2TMapComposer::produce(edm::Event& event, const edm::EventSetup&) {
  std::ostringstream out;

  for (size_t i = 0; i < map_tags.size(); ++i) {
    // Need refs to first and last track collections with changes to
    // AssociationMap in 750.
    edm::Handle<reco::TrackCollection> first_tracks, last_tracks;
    event.getByLabel(first_track_tags[i], first_tracks);
    event.getByLabel(last_track_tags [i], last_tracks);

    // Get all the maps we're to compose and put them in the vector in order.
    map_vector maps(map_tags[i].size());
    int itag = 0;
    for (std::vector<edm::InputTag>::const_iterator tag = map_tags[i].begin(); tag != map_tags[i].end(); ++tag, ++itag) {
      event.getByLabel(*tag, maps.at(itag));
      
      if (debug) {
	out << "src #" << itag << ": " << tag->encode() << " t2t map\n";
	dump_t2tmap(out, *maps.at(itag), &event);
      }
    }

    std::auto_ptr<reco::TrackToTrackMap> t2tmap(compose(event, first_tracks, last_tracks, maps));

    if (debug) {
      out << "composed map with name '" << new_map_names[i] << "':\n";
      dump_t2tmap(out, *t2tmap, &event);
    }

    event.put(t2tmap, new_map_names[i]);
  }

  edm::LogInfo("T2TMapComposer") << out.str();
}

DEFINE_FWK_MODULE(T2TMapComposer);
