#ifndef CosmicSplittingResolutionNtuple_h
#define CosmicSplittingResolutionNtuple_h

#include "TTree.h"

// The possible values of dataset_id in
// CosmicSplittingResolutionFilter. Used to distinguish data and the three
// different bins of simulation.
enum data_id_type { data, mc_10, mc_100, mc_500 };

// Below in the ntuple definition we store e.g. pT for each of the
// track types separately in arrays. These two enums define the layout
// of the arrays.
//
// tk_staglb contains global tracks propagated to the stand-alone
// tracks' positions, and tk_statko is the same for tracker-only
// tracks. This is for the comparison global vs. standalone or
// tracker-only vs. standalone. (Therefore there are double the
// statistics for this comparison: upper and lower pairs.)
enum track_type { tk_global, tk_stalone, tk_tkonly, tk_tpfms, tk_picky, tk_dyt, tk_tmr, tk_sigsw, tk_tunep, tk_tunep_old, tk_tuner_sigma, tk_tuner_prob, tk_tuner_sigma2, tk_staglb, tk_statko, n_tracks };
enum track_pos { upper, lower, n_track_pos };

extern const char* track_nicks[n_tracks];

struct CosmicSplittingResolutionNtuple {
  typedef unsigned short ushort;

  unsigned id;
  //unsigned new_id;

  unsigned run;
  unsigned lumi;
  unsigned event;

  ushort shared_hits[n_tracks];

  float prop_pos[n_tracks][n_track_pos][3];
  float inner_pos[n_tracks][n_track_pos][3];
  float outer_pos[n_tracks][n_track_pos][3];

  short  charge[n_tracks][n_track_pos];
  float  chi2[n_tracks][n_track_pos];
  ushort ndof[n_tracks][n_track_pos];
  float  prob[n_tracks][n_track_pos];
  ushort pixel_layers[n_tracks][n_track_pos];
  ushort strip_layers[n_tracks][n_track_pos];
  ushort pixel_hits[n_tracks][n_track_pos];
  ushort strip_hits[n_tracks][n_track_pos];
  ushort muon_hits[n_tracks][n_track_pos];
  ushort muon_hits_DT[n_tracks][n_track_pos];  
  ushort muon_hits_CSC[n_tracks][n_track_pos];
  ushort muon_hits_RPC[n_tracks][n_track_pos];

  float qoverp[n_tracks][n_track_pos];
  float error_qoverp[n_tracks][n_track_pos];

  bool prop_valid[n_tracks][n_track_pos];

  float pt[n_tracks][n_track_pos];
  float theta[n_tracks][n_track_pos];
  float phi[n_tracks][n_track_pos];
  float dxy[n_tracks][n_track_pos];
  float dz[n_tracks][n_track_pos];

  float bz[n_tracks][n_track_pos];
  float kappa[n_tracks][n_track_pos];
  float error_kappa[n_tracks][n_track_pos];

  float error_pt[n_tracks][n_track_pos];
  float error_theta[n_tracks][n_track_pos];
  float error_phi[n_tracks][n_track_pos];
  float error_dxy[n_tracks][n_track_pos];
  float error_dz[n_tracks][n_track_pos];

  short unprop_charge[n_tracks][n_track_pos];
  float unprop_pt[n_tracks][n_track_pos];
  float unprop_theta[n_tracks][n_track_pos];
  float unprop_phi[n_tracks][n_track_pos];
  float unprop_dxy[n_tracks][n_track_pos];
  float unprop_dz[n_tracks][n_track_pos];

  float unprop_error_pt[n_tracks][n_track_pos];
  float unprop_error_theta[n_tracks][n_track_pos];
  float unprop_error_phi[n_tracks][n_track_pos];
  float unprop_error_dxy[n_tracks][n_track_pos];
  float unprop_error_dz[n_tracks][n_track_pos];

  float mc_vertex[3];
  short unprop_mc_charge;
  float unprop_mc_pt;
  float unprop_mc_theta;
  float unprop_mc_phi;

  float mc_prop_pos[n_tracks][3];

  bool mc_prop_valid[n_tracks];

  short mc_charge[n_tracks];
  float mc_pt[n_tracks];
  float mc_theta[n_tracks];
  float mc_phi[n_tracks];
  float mc_dxy[n_tracks];
  float mc_dz[n_tracks];

  // L1Particle quantities
  ushort nl1;
  float l1_pt[4];
  float l1_eta[4];
  float l1_phi[4];
  short l1_charge[4];
  unsigned int l1_quality[4];
  int l1_bx[4];
  bool l1_isol[4];
  bool l1_isFwd[4];
  bool l1_isRPC[4];

  // Muon Stations
  ushort muon_stations[n_track_pos];

  short  ref_charge;
  float  ref_chi2;
  ushort ref_ndof;
  ushort ref_pixel_layers;
  ushort ref_strip_layers;
  ushort ref_pixel_hits;
  ushort ref_strip_hits;
  ushort ref_muon_hits;
  ushort ref_muon_hits_DT;
  ushort ref_muon_hits_CSC;
  ushort ref_muon_hits_RPC;

  float ref_pt;
  float ref_theta;
  float ref_phi;
  float ref_dxy;
  float ref_dz;

  float ref_error_pt;
  float ref_error_theta;
  float ref_error_phi;
  float ref_error_dxy;
  float ref_error_dz;

  bool   hit_dt;
  bool   hit_csc;
  ushort tpfms_first_station[n_track_pos];

  ushort choice_tmr[n_track_pos];
  ushort choice_sigsw[n_track_pos];
  ushort choice_tunep[n_track_pos];
  ushort choice_tunep_old[n_track_pos];
  ushort choice_tuner_sigma[n_track_pos];
  ushort choice_tuner_prob[n_track_pos];
  ushort choice_tuner_sigma2[n_track_pos];

  bool tt25;

  float bzat0;

  bool peak_mode; // false == deconvolution
};

void write_to_tree (TTree*, CosmicSplittingResolutionNtuple*);
void read_from_tree(TTree*, CosmicSplittingResolutionNtuple*);

#endif
