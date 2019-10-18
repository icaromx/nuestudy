//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 16 11:58:01 2019 by ROOT version 6.18/00
// from TTree NeutrinoSelectionFilter/Neutrino Selection TTree
// found on file: /uboone/data/users/davidc/searchingfornues/v08_00_00_25/cc0pinp/1013/prodgenie_bnb_intrinsic_nue_uboone_overlay_mcc9.1_run3_G_reco2.root
//////////////////////////////////////////////////////////

#ifndef NeutrinoSelectionFilter_h
#define NeutrinoSelectionFilter_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
//#include <da_cos6-linux-gnu/include/c++/7.3.0/vector>
//#include <da_cos6-linux-gnu/include/c++/7.3.0/vector>
#include "vector"
#include "vector"
#include "string"
#include "vector"

class NeutrinoSelectionFilter {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           selected;
   UInt_t          trk_id;
   UInt_t          shr_id;
   Float_t         shr_energy_tot;
   Float_t         shr_energy;
   Float_t         shr_energy_tot_cali;
   Float_t         shr_energy_cali;
   Float_t         shr_theta;
   Float_t         shr_phi;
   Float_t         shr_pca_0;
   Float_t         shr_pca_1;
   Float_t         shr_pca_2;
   Float_t         shr_px;
   Float_t         shr_py;
   Float_t         shr_pz;
   Float_t         shr_openangle;
   Float_t         shr_tkfit_start_x;
   Float_t         shr_tkfit_start_y;
   Float_t         shr_tkfit_start_z;
   Float_t         shr_tkfit_theta;
   Float_t         shr_tkfit_phi;
   Float_t         shr_start_x;
   Float_t         shr_start_y;
   Float_t         shr_start_z;
   Float_t         shr_dedx_Y;
   Float_t         shr_dedx_V;
   Float_t         shr_dedx_U;
   Float_t         shr_dedx_Y_cali;
   Float_t         shr_dedx_V_cali;
   Float_t         shr_dedx_U_cali;
   Float_t         shr_tkfit_dedx_Y;
   Float_t         shr_tkfit_dedx_V;
   Float_t         shr_tkfit_dedx_U;
   UInt_t          shr_tkfit_nhits_Y;
   UInt_t          shr_tkfit_nhits_V;
   UInt_t          shr_tkfit_nhits_U;
   Float_t         shr_chipr;
   Float_t         shr_chimu;
   Float_t         shr_bragg_p;
   Float_t         shr_bragg_mu;
   Float_t         shr_bragg_mip;
   Float_t         shr_bragg_kaon;
   Float_t         shr_bragg_pion;
   Float_t         tksh_distance;
   Float_t         tksh_angle;
   Float_t         shr_distance;
   Float_t         shr_score;
   Int_t           shr_bkt_pdg;
   Float_t         shr_bkt_purity;
   Float_t         shr_bkt_completeness;
   Float_t         shr_bkt_E;
   Float_t         trk_len;
   Float_t         trk_theta;
   Float_t         trk_phi;
   Float_t         trk_energy;
   Float_t         trk_energy_muon;
   Float_t         trk_energy_tot;
   Float_t         trk_energy_muon_tot;
   Float_t         trk_distance;
   Float_t         trk_score;
   Int_t           trk_bkt_pdg;
   Float_t         trk_bkt_purity;
   Float_t         trk_bkt_completeness;
   Float_t         trk_bkt_E;
   Float_t         trk_chipr_best;
   Float_t         trk_chipr_worst;
   Float_t         trk_chimu_best;
   Float_t         trk_chimu_worst;
   Float_t         trk_chipr;
   Float_t         trk_chimu;
   Float_t         trk_pida;
   Float_t         trk_bragg_p;
   Float_t         trk_bragg_mu;
   Float_t         trk_bragg_mip;
   Float_t         trk_bragg_kaon;
   Float_t         trk_bragg_pion;
   UInt_t          trk_hits_max;
   UInt_t          shr_hits_max;
   UInt_t          total_hits_y;
   Float_t         extra_energy_y;
   Float_t         trk_energy_hits_tot;
   UInt_t          shr_hits_tot;
   UInt_t          shr_hits_y_tot;
   UInt_t          shr_hits_u_tot;
   UInt_t          shr_hits_v_tot;
   UInt_t          trk_hits_tot;
   UInt_t          trk_hits_y_tot;
   UInt_t          trk_hits_u_tot;
   UInt_t          trk_hits_v_tot;
   UInt_t          n_tracks_contained;
   UInt_t          n_showers_contained;
   Float_t         matched_E;
   Float_t         hits_ratio;
   Float_t         contained_fraction;
   Float_t         pt;
   Float_t         p;
   Float_t         pt_assume_muon;
   Float_t         p_assume_muon;
   Float_t         dvtx;
   Float_t         dtrk;
   vector<double>  *dtrk_x_boundary;
   vector<double>  *dtrk_y_boundary;
   vector<double>  *dtrk_z_boundary;
   vector<double>  *dshr_x_boundary;
   vector<double>  *dshr_y_boundary;
   vector<double>  *dshr_z_boundary;
   vector<double>  *dvtx_x_boundary;
   vector<double>  *dvtx_y_boundary;
   vector<double>  *dvtx_z_boundary;
   vector<vector<double> > *dtrk_boundary;
   vector<vector<double> > *dvtx_boundary;
   vector<vector<double> > *dshr_boundary;
   vector<vector<double> > *dmc_boundary;
   Float_t         CosmicIP;
   Int_t           _run;
   Int_t           _sub;
   Int_t           _evt;
   Int_t           _nAssnCosmics;
   Int_t           _kMaxCosm;
   Double_t        _t0Cosm_v[200];   //[_kMaxCosm]
   Double_t        _t0timesVCosm_v[200];   //[_kMaxCosm]
   Double_t        _xStartCosm_v[200];   //[_kMaxCosm]
   Double_t        _yStartCosm_v[200];   //[_kMaxCosm]
   Double_t        _zStartCosm_v[200];   //[_kMaxCosm]
   Double_t        _xEndCosm_v[200];   //[_kMaxCosm]
   Double_t        _yEndCosm_v[200];   //[_kMaxCosm]
   Double_t        _zEndCosm_v[200];   //[_kMaxCosm]
   Double_t        _recoNu_vtx_x;
   Double_t        _recoNu_vtx_y;
   Double_t        _recoNu_vtx_z;
   Double_t        _t0_nu_cosmic;
   Double_t        _nu_cosmic_x;
   Double_t        _nu_cosmic_y;
   Double_t        _nu_cosmic_z;
   Double_t        _nu_cosmic_Length;
   Double_t        _nu_cosmic_Start_x;
   Double_t        _nu_cosmic_Start_y;
   Double_t        _nu_cosmic_Start_z;
   Double_t        _nu_cosmic_End_x;
   Double_t        _nu_cosmic_End_y;
   Double_t        _nu_cosmic_End_z;
   Double_t        _nu_cosmic_TrackID;
   Double_t        _closestNuCosmicDist;
   Double_t        _rand_vtx_x;
   Double_t        _rand_vtx_y;
   Double_t        _rand_vtx_z;
   Double_t        _t0_rand_cosmic;
   Double_t        _rand_cosmic_x;
   Double_t        _rand_cosmic_y;
   Double_t        _rand_cosmic_z;
   Double_t        _rand_cosmic_Length;
   Double_t        _rand_cosmic_Start_x;
   Double_t        _rand_cosmic_Start_y;
   Double_t        _rand_cosmic_Start_z;
   Double_t        _rand_cosmic_End_x;
   Double_t        _rand_cosmic_End_y;
   Double_t        _rand_cosmic_End_z;
   Double_t        _rand_cosmic_TrackID;
   Double_t        _closestRandCosmicDist;
   Float_t         leeweight;
   Float_t         nu_vtx_x;
   Float_t         nu_vtx_y;
   Float_t         nu_vtx_z;
   Float_t         nu_sce_x;
   Float_t         nu_sce_y;
   Float_t         nu_sce_z;
   Float_t         true_pt;
   Float_t         true_pt_visible;
   Float_t         true_p;
   Float_t         true_p_visible;
   Float_t         true_e_visible;
   Int_t           nu_pdg;
   Int_t           ccnc;
   Int_t           interaction;
   Float_t         nu_e;
   Float_t         nu_pt;
   Float_t         theta;
   Float_t         vtx_x;
   Float_t         vtx_y;
   Float_t         vtx_z;
   Bool_t          isVtxInActive;
   Bool_t          isVtxInFiducial;
   Float_t         true_nu_vtx_t;
   Float_t         true_nu_vtx_x;
   Float_t         true_nu_vtx_y;
   Float_t         true_nu_vtx_z;
   Float_t         true_nu_vtx_sce_x;
   Float_t         true_nu_vtx_sce_y;
   Float_t         true_nu_vtx_sce_z;
   Float_t         reco_nu_vtx_x;
   Float_t         reco_nu_vtx_y;
   Float_t         reco_nu_vtx_z;
   Float_t         reco_nu_vtx_sce_x;
   Float_t         reco_nu_vtx_sce_y;
   Float_t         reco_nu_vtx_sce_z;
   Int_t           nmuon;
   Float_t         muon_e;
   Float_t         muon_c;
   Float_t         muon_p;
   Int_t           nelec;
   Float_t         elec_e;
   Float_t         elec_c;
   Float_t         elec_p;
   Float_t         elec_vx;
   Float_t         elec_vy;
   Float_t         elec_vz;
   Int_t           npi0;
   Float_t         pi0_e;
   Float_t         pi0_c;
   Float_t         pi0_p;
   Int_t           nneutron;
   Int_t           nproton;
   Float_t         proton_e;
   Float_t         proton_c;
   Float_t         proton_p;
   Int_t           npion;
   Float_t         pion_e;
   Float_t         pion_c;
   Float_t         pion_p;
   Int_t           nslice;
   Int_t           crtveto;
   Float_t         crthitpe;
   vector<int>     *pfp_slice_idx;
   Int_t           category;
   vector<int>     *backtracked_pdg;
   vector<float>   *backtracked_e;
   vector<float>   *backtracked_purity;
   vector<float>   *backtracked_completeness;
   vector<float>   *backtracked_overlay_purity;
   vector<float>   *backtracked_px;
   vector<float>   *backtracked_py;
   vector<float>   *backtracked_pz;
   vector<float>   *backtracked_start_x;
   vector<float>   *backtracked_start_y;
   vector<float>   *backtracked_start_z;
   vector<float>   *backtracked_start_t;
   vector<float>   *backtracked_start_U;
   vector<float>   *backtracked_start_V;
   vector<float>   *backtracked_start_Y;
   vector<float>   *backtracked_sce_start_x;
   vector<float>   *backtracked_sce_start_y;
   vector<float>   *backtracked_sce_start_z;
   vector<float>   *backtracked_sce_start_U;
   vector<float>   *backtracked_sce_start_V;
   vector<float>   *backtracked_sce_start_Y;
   Float_t         lep_e;
   Int_t           pass;
   Int_t           run;
   Int_t           sub;
   Int_t           evt;
   Int_t           swtrig;
   Float_t         xtimeoffset;
   Float_t         xsceoffset;
   Float_t         ysceoffset;
   Float_t         zsceoffset;
   Int_t           evnhits;
   Int_t           slpdg;
   Int_t           slnhits;
   Int_t           n_pfps;
   Int_t           n_tracks;
   Int_t           n_showers;
   vector<float>   *trk_score_v;
   vector<int>     *pfpdg;
   vector<int>     *pfnhits;
   vector<int>     *pfnplanehits_U;
   vector<int>     *pfnplanehits_V;
   vector<int>     *pfnplanehits_Y;
   UInt_t          hits_u;
   UInt_t          hits_v;
   UInt_t          hits_y;
   Float_t         topological_score;
   Float_t         slclustfrac;
   vector<int>     *mc_pdg;
   vector<float>   *mc_E;
   vector<float>   *mc_vx;
   vector<float>   *mc_vy;
   vector<float>   *mc_vz;
   vector<float>   *mc_endx;
   vector<float>   *mc_endy;
   vector<float>   *mc_endz;
   vector<float>   *mc_px;
   vector<float>   *mc_py;
   vector<float>   *mc_pz;
   vector<float>   *mc_completeness;
   vector<float>   *mc_purity;
   string          *endmuonprocess;
   Float_t         endmuonmichel;
   map<string,vector<double> > *weights;
   vector<double>  *weightsFlux;
   vector<double>  *weightsGenie;
   vector<double>  *weightsReint;
   Float_t         weightSpline;
   Float_t         nu_flashmatch_score;
   Float_t         best_cosmic_flashmatch_score;
   Float_t         best_obviouscosmic_flashmatch_score;
   vector<float>   *cosmic_flashmatch_score_v;
   Float_t         NeutrinoEnergy0;
   Float_t         NeutrinoEnergy1;
   Float_t         NeutrinoEnergy2;
   Float_t         SliceCaloEnergy0;
   Float_t         SliceCaloEnergy1;
   Float_t         SliceCaloEnergy2;
   Float_t         gamma_edep;
   Float_t         gamma_etot;
   Float_t         gamma_dist;
   Int_t           gamma_parent;
   Float_t         elec_edep;
   Float_t         elec_etot;
   Float_t         elec_dist;
   Int_t           elec_parent;
   Float_t         gamma1_edep;
   Float_t         gamma1_etot;
   Float_t         gamma2_edep;
   Float_t         gamma2_etot;
   Int_t           nflag_pl1;
   Int_t           nnoise_pl1;
   Int_t           nslhits_pl1;
   Int_t           nslnoise_pl1;
   Int_t           nhits_pl1;
   Float_t         frac_slnoise_pl1;
   vector<float>   *shr_dedx_u_v;
   vector<float>   *shr_dedx_v_v;
   vector<float>   *shr_dedx_y_v;
   vector<float>   *shr_energy_u_v;
   vector<float>   *shr_energy_v_v;
   vector<float>   *shr_energy_y_v;
   vector<unsigned long> *shr_pfp_id_v;
   vector<float>   *shr_start_x_v;
   vector<float>   *shr_start_y_v;
   vector<float>   *shr_start_z_v;
   vector<float>   *shr_start_U_v;
   vector<float>   *shr_start_V_v;
   vector<float>   *shr_dist_v;
   vector<int>     *shr_nclus_v;
   vector<float>   *shr_clushitfrac_v;
   vector<float>   *shr_px_v;
   vector<float>   *shr_py_v;
   vector<float>   *shr_pz_v;
   vector<float>   *shr_openangle_v;
   vector<float>   *shr_theta_v;
   vector<float>   *shr_phi_v;
   vector<int>     *shr_tkfit_nhits_v;
   vector<float>   *shr_tkfit_start_x_v;
   vector<float>   *shr_tkfit_start_y_v;
   vector<float>   *shr_tkfit_start_z_v;
   vector<float>   *shr_tkfit_start_U_v;
   vector<float>   *shr_tkfit_start_V_v;
   vector<float>   *shr_tkfit_theta_v;
   vector<float>   *shr_tkfit_phi_v;
   vector<float>   *shr_tkfit_dedx_u_v;
   vector<float>   *shr_tkfit_dedx_v_v;
   vector<float>   *shr_tkfit_dedx_y_v;
   vector<int>     *shr_tkfit_dedx_nhits_u_v;
   vector<int>     *shr_tkfit_dedx_nhits_v_v;
   vector<int>     *shr_tkfit_dedx_nhits_y_v;
   vector<float>   *shr_spacepoint_start_x_v;
   vector<float>   *shr_spacepoint_start_y_v;
   vector<float>   *shr_spacepoint_start_z_v;
   vector<float>   *shr_spacepoint_start_U_v;
   vector<float>   *shr_spacepoint_start_V_v;
   vector<float>   *shr_hits_start_U_wire_v;
   vector<float>   *shr_hits_start_U_x_v;
   vector<float>   *shr_hits_start_V_wire_v;
   vector<float>   *shr_hits_start_V_x_v;
   vector<float>   *shr_hits_start_Y_wire_v;
   vector<float>   *shr_hits_start_Y_x_v;
   Int_t           evnunhits;
   Int_t           evlepnhits;
   Int_t           evpronhits;
   Int_t           evpi1nhits;
   Int_t           evpi0nhits;
   Int_t           evneunhits;
   Int_t           evgamnhits;
   Int_t           evothnhits;
   Int_t           slnunhits;
   Int_t           sllepnhits;
   Int_t           slpronhits;
   Int_t           slpi1nhits;
   Int_t           slpi0nhits;
   Int_t           slneunhits;
   Int_t           slgamnhits;
   Int_t           slothnhits;
   vector<int>     *pfnunhits;
   vector<int>     *pflepnhits;
   vector<int>     *pfpronhits;
   vector<int>     *pfpi1nhits;
   vector<int>     *pfpi0nhits;
   vector<int>     *pfneunhits;
   vector<int>     *pfgamnhits;
   vector<int>     *pfothnhits;
   vector<float>   *trk_bragg_p_v;
   vector<float>   *trk_bragg_mu_v;
   vector<float>   *trk_bragg_mip_v;
   vector<float>   *trk_pida_v;
   vector<float>   *trk_pid_chipr_v;
   vector<float>   *trk_pid_chipi_v;
   vector<float>   *trk_pid_chika_v;
   vector<float>   *trk_pid_chimu_v;
   vector<float>   *trk_bragg_p_u_v;
   vector<float>   *trk_bragg_mu_u_v;
   vector<float>   *trk_bragg_mip_u_v;
   vector<float>   *trk_pida_u_v;
   vector<float>   *trk_pid_chipr_u_v;
   vector<float>   *trk_pid_chipi_u_v;
   vector<float>   *trk_pid_chika_u_v;
   vector<float>   *trk_pid_chimu_u_v;
   vector<float>   *trk_bragg_p_v_v;
   vector<float>   *trk_bragg_mu_v_v;
   vector<float>   *trk_bragg_mip_v_v;
   vector<float>   *trk_pida_v_v;
   vector<float>   *trk_pid_chipr_v_v;
   vector<float>   *trk_pid_chipi_v_v;
   vector<float>   *trk_pid_chika_v_v;
   vector<float>   *trk_pid_chimu_v_v;
   vector<unsigned long> *trk_pfp_id_v;
   vector<float>   *trk_dir_x_v;
   vector<float>   *trk_dir_y_v;
   vector<float>   *trk_dir_z_v;
   vector<float>   *trk_start_x_v;
   vector<float>   *trk_start_y_v;
   vector<float>   *trk_start_z_v;
   vector<float>   *trk_end_x_v;
   vector<float>   *trk_end_y_v;
   vector<float>   *trk_end_z_v;
   vector<float>   *trk_distance_v;
   vector<float>   *trk_theta_v;
   vector<float>   *trk_phi_v;
   vector<float>   *trk_len_v;
   vector<float>   *trk_mcs_muon_mom_v;
   vector<float>   *trk_energy_proton_v;
   vector<float>   *trk_energy_muon_v;
   Float_t         bdt_nuNCpi0;
   Float_t         bdt_numuCCpi0;
   Float_t         bdt_numuCC;
   Float_t         bdt_ext;
   Float_t         bdt_cosmic;
   Float_t         bdt_global;

   // List of branches
   TBranch        *b_selected;   //!
   TBranch        *b_trk_pfp_id;   //!
   TBranch        *b_shr_pfp_id;   //!
   TBranch        *b_shr_energy_tot;   //!
   TBranch        *b_shr_energy;   //!
   TBranch        *b_shr_energy_tot_cali;   //!
   TBranch        *b_shr_energy_cali;   //!
   TBranch        *b_shr_theta;   //!
   TBranch        *b_shr_phi;   //!
   TBranch        *b_shr_pca_0;   //!
   TBranch        *b_shr_pca_1;   //!
   TBranch        *b_shr_pca_2;   //!
   TBranch        *b_shr_px;   //!
   TBranch        *b_shr_py;   //!
   TBranch        *b_shr_pz;   //!
   TBranch        *b_shr_openangle;   //!
   TBranch        *b_shr_tkfit_start_x;   //!
   TBranch        *b_shr_tkfit_start_y;   //!
   TBranch        *b_shr_tkfit_start_z;   //!
   TBranch        *b_shr_tkfit_theta;   //!
   TBranch        *b_shr_tkfit_phi;   //!
   TBranch        *b_shr_start_x;   //!
   TBranch        *b_shr_start_y;   //!
   TBranch        *b_shr_start_z;   //!
   TBranch        *b_shr_dedx_Y;   //!
   TBranch        *b_shr_dedx_V;   //!
   TBranch        *b_shr_dedx_U;   //!
   TBranch        *b_shr_dedx_Y_cali;   //!
   TBranch        *b_shr_dedx_V_cali;   //!
   TBranch        *b_shr_dedx_U_cali;   //!
   TBranch        *b_shr_tkfit_dedx_Y;   //!
   TBranch        *b_shr_tkfit_dedx_V;   //!
   TBranch        *b_shr_tkfit_dedx_U;   //!
   TBranch        *b_shr_tkfit_nhits_Y;   //!
   TBranch        *b_shr_tkfit_nhits_V;   //!
   TBranch        *b_shr_tkfit_nhits_U;   //!
   TBranch        *b_shr_chipr;   //!
   TBranch        *b_shr_chimu;   //!
   TBranch        *b_shr_bragg_p;   //!
   TBranch        *b_shr_bragg_mu;   //!
   TBranch        *b_shr_bragg_mip;   //!
   TBranch        *b_shr_bragg_kaon;   //!
   TBranch        *b_shr_bragg_pion;   //!
   TBranch        *b_tksh_distance;   //!
   TBranch        *b_tksh_angle;   //!
   TBranch        *b_shr_distance;   //!
   TBranch        *b_shr_score;   //!
   TBranch        *b_shr_bkt_pdg;   //!
   TBranch        *b_shr_bkt_purity;   //!
   TBranch        *b_shr_bkt_completeness;   //!
   TBranch        *b_shr_bkt_E;   //!
   TBranch        *b_trk_len;   //!
   TBranch        *b_trk_theta;   //!
   TBranch        *b_trk_phi;   //!
   TBranch        *b_trk_energy;   //!
   TBranch        *b_trk_energy_muon;   //!
   TBranch        *b_trk_energy_tot;   //!
   TBranch        *b_trk_energy_muon_tot;   //!
   TBranch        *b_trk_distance;   //!
   TBranch        *b_trk_score;   //!
   TBranch        *b_trk_bkt_pdg;   //!
   TBranch        *b_trk_bkt_purity;   //!
   TBranch        *b_trk_bkt_completeness;   //!
   TBranch        *b_trk_bkt_E;   //!
   TBranch        *b_trk_chipr_best;   //!
   TBranch        *b_trk_chipr_worst;   //!
   TBranch        *b_trk_chimu_best;   //!
   TBranch        *b_trk_chimu_worst;   //!
   TBranch        *b_trk_chipr;   //!
   TBranch        *b_trk_chimu;   //!
   TBranch        *b_trk_pida;   //!
   TBranch        *b_trk_bragg_p;   //!
   TBranch        *b_trk_bragg_mu;   //!
   TBranch        *b_trk_bragg_mip;   //!
   TBranch        *b_trk_bragg_kaon;   //!
   TBranch        *b_trk_bragg_pion;   //!
   TBranch        *b_trk_hits_max;   //!
   TBranch        *b_shr_hits_max;   //!
   TBranch        *b_total_hits_y;   //!
   TBranch        *b_extra_energy_y;   //!
   TBranch        *b_trk_energy_hits_tot;   //!
   TBranch        *b_shr_hits_tot;   //!
   TBranch        *b_shr_hits_y_tot;   //!
   TBranch        *b_shr_hits_u_tot;   //!
   TBranch        *b_shr_hits_v_tot;   //!
   TBranch        *b_trk_hits_tot;   //!
   TBranch        *b_trk_hits_y_tot;   //!
   TBranch        *b_trk_hits_u_tot;   //!
   TBranch        *b_trk_hits_v_tot;   //!
   TBranch        *b_n_tracks_contained;   //!
   TBranch        *b_n_showers_contained;   //!
   TBranch        *b_matched_E;   //!
   TBranch        *b_hits_ratio;   //!
   TBranch        *b_contained_fraction;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_p;   //!
   TBranch        *b_pt_assume_muon;   //!
   TBranch        *b_p_assume_muon;   //!
   TBranch        *b_dvtx;   //!
   TBranch        *b_dtrk;   //!
   TBranch        *b_dtrk_x_boundary;   //!
   TBranch        *b_dtrk_y_boundary;   //!
   TBranch        *b_dtrk_z_boundary;   //!
   TBranch        *b_dshr_x_boundary;   //!
   TBranch        *b_dshr_y_boundary;   //!
   TBranch        *b_dshr_z_boundary;   //!
   TBranch        *b_dvtx_x_boundary;   //!
   TBranch        *b_dvtx_y_boundary;   //!
   TBranch        *b_dvtx_z_boundary;   //!
   TBranch        *b_dtrk_boundary;   //!
   TBranch        *b_dvtx_boundary;   //!
   TBranch        *b_dshr_boundary;   //!
   TBranch        *b_dmc_boundary;   //!
   TBranch        *b_CosmicIP;   //!
   TBranch        *b_run;   //!
   TBranch        *b_sub;   //!
   TBranch        *b_evt;   //!
   TBranch        *b__nAssnCosmics;   //!
   TBranch        *b__kMaxCosm;   //!
   TBranch        *b__t0Cosm_v;   //!
   TBranch        *b__t0timesVCosm_v;   //!
   TBranch        *b__xStartCosm_v;   //!
   TBranch        *b__yStartCosm_v;   //!
   TBranch        *b__zStartCosm_v;   //!
   TBranch        *b__xEndCosm_v;   //!
   TBranch        *b__yEndCosm_v;   //!
   TBranch        *b__zEndCosm_v;   //!
   TBranch        *b_reco_vtx_x;   //!
   TBranch        *b_reco_vtx_y;   //!
   TBranch        *b_reco_vtx_z;   //!
   TBranch        *b__t0_nu_cosmic;   //!
   TBranch        *b__nu_cosmic_x;   //!
   TBranch        *b__nu_cosmic_y;   //!
   TBranch        *b__nu_cosmic_z;   //!
   TBranch        *b__nu_cosmic_Length;   //!
   TBranch        *b__nu_cosmic_Start_x;   //!
   TBranch        *b__nu_cosmic_Start_y;   //!
   TBranch        *b__nu_cosmic_Start_z;   //!
   TBranch        *b__nu_cosmic_End_x;   //!
   TBranch        *b__nu_cosmic_End_y;   //!
   TBranch        *b__nu_cosmic_End_z;   //!
   TBranch        *b__nu_cosmic_TrackID;   //!
   TBranch        *b__closestNuCosmicDist;   //!
   TBranch        *b_rand_vtx_x;   //!
   TBranch        *b_rand_vtx_y;   //!
   TBranch        *b_rand_vtx_z;   //!
   TBranch        *b__t0_rand_cosmic;   //!
   TBranch        *b__rand_cosmic_x;   //!
   TBranch        *b__rand_cosmic_y;   //!
   TBranch        *b__rand_cosmic_z;   //!
   TBranch        *b__rand_cosmic_Length;   //!
   TBranch        *b__rand_cosmic_Start_x;   //!
   TBranch        *b__rand_cosmic_Start_y;   //!
   TBranch        *b__rand_cosmic_Start_z;   //!
   TBranch        *b__rand_cosmic_End_x;   //!
   TBranch        *b__rand_cosmic_End_y;   //!
   TBranch        *b__rand_cosmic_End_z;   //!
   TBranch        *b__rand_cosmic_TrackID;   //!
   TBranch        *b__closestRandCosmicDist;   //!
   TBranch        *b_leeweight;   //!
   TBranch        *b_nu_vtx_x;   //!
   TBranch        *b_nu_vtx_y;   //!
   TBranch        *b_nu_vtx_z;   //!
   TBranch        *b_nu_sce_x;   //!
   TBranch        *b_nu_sce_y;   //!
   TBranch        *b_nu_sce_z;   //!
   TBranch        *b_true_pt;   //!
   TBranch        *b_true_pt_visible;   //!
   TBranch        *b_true_p;   //!
   TBranch        *b_true_p_visible;   //!
   TBranch        *b_true_e_visible;   //!
   TBranch        *b_nu_pdg;   //!
   TBranch        *b_ccnc;   //!
   TBranch        *b_interaction;   //!
   TBranch        *b_nu_e;   //!
   TBranch        *b_nu_pt;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_vtx_x;   //!
   TBranch        *b_vtx_y;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_isVtxInActive;   //!
   TBranch        *b_isVtxInFiducial;   //!
   TBranch        *b_true_nu_vtx_t;   //!
   TBranch        *b_true_nu_vtx_x;   //!
   TBranch        *b_true_nu_vtx_y;   //!
   TBranch        *b_true_nu_vtx_z;   //!
   TBranch        *b_true_nu_vtx_sce_x;   //!
   TBranch        *b_true_nu_vtx_sce_y;   //!
   TBranch        *b_true_nu_vtx_sce_z;   //!
   TBranch        *b_reco_nu_vtx_x;   //!
   TBranch        *b_reco_nu_vtx_y;   //!
   TBranch        *b_reco_nu_vtx_z;   //!
   TBranch        *b_reco_nu_vtx_sce_x;   //!
   TBranch        *b_reco_nu_vtx_sce_y;   //!
   TBranch        *b_reco_nu_vtx_sce_z;   //!
   TBranch        *b_nmuon;   //!
   TBranch        *b_muon_e;   //!
   TBranch        *b_muon_c;   //!
   TBranch        *b_muon_p;   //!
   TBranch        *b_nelec;   //!
   TBranch        *b_elec_e;   //!
   TBranch        *b_elec_c;   //!
   TBranch        *b_elec_p;   //!
   TBranch        *b_elec_vx;   //!
   TBranch        *b_elec_vy;   //!
   TBranch        *b_elec_vz;   //!
   TBranch        *b_npi0;   //!
   TBranch        *b_pi0_e;   //!
   TBranch        *b_pi0_c;   //!
   TBranch        *b_pi0_p;   //!
   TBranch        *b_nneutron;   //!
   TBranch        *b_nproton;   //!
   TBranch        *b_proton_e;   //!
   TBranch        *b_proton_c;   //!
   TBranch        *b_proton_p;   //!
   TBranch        *b_npion;   //!
   TBranch        *b_pion_e;   //!
   TBranch        *b_pion_c;   //!
   TBranch        *b_pion_p;   //!
   TBranch        *b_nslice;   //!
   TBranch        *b_crtveto;   //!
   TBranch        *b_crthitpe;   //!
   TBranch        *b_pfp_slice_idx;   //!
   TBranch        *b_category;   //!
   TBranch        *b_backtracked_pdg;   //!
   TBranch        *b_backtracked_e;   //!
   TBranch        *b_backtracked_purity;   //!
   TBranch        *b_backtracked_completeness;   //!
   TBranch        *b_backtracked_overlay_purity;   //!
   TBranch        *b_backtracked_px;   //!
   TBranch        *b_backtracked_py;   //!
   TBranch        *b_backtracked_pz;   //!
   TBranch        *b_backtracked_start_x;   //!
   TBranch        *b_backtracked_start_y;   //!
   TBranch        *b_backtracked_start_z;   //!
   TBranch        *b_backtracked_start_t;   //!
   TBranch        *b_backtracked_start_U;   //!
   TBranch        *b_backtracked_start_V;   //!
   TBranch        *b_backtracked_start_Y;   //!
   TBranch        *b_backtracked_sce_start_x;   //!
   TBranch        *b_backtracked_sce_start_y;   //!
   TBranch        *b_backtracked_sce_start_z;   //!
   TBranch        *b_backtracked_sce_start_U;   //!
   TBranch        *b_backtracked_sce_start_V;   //!
   TBranch        *b_backtracked_sce_start_Y;   //!
   TBranch        *b_lep_e;   //!
   TBranch        *b_pass;   //!
   TBranch        *b_swtrig;   //!
   TBranch        *b_xtimeoffset;   //!
   TBranch        *b_xsceoffset;   //!
   TBranch        *b_ysceoffset;   //!
   TBranch        *b_zsceoffset;   //!
   TBranch        *b_evnhits;   //!
   TBranch        *b_slpdg;   //!
   TBranch        *b_slnhits;   //!
   TBranch        *b_n_pfps;   //!
   TBranch        *b_n_tracks;   //!
   TBranch        *b_n_showers;   //!
   TBranch        *b_trk_score_v;   //!
   TBranch        *b_pfpdg;   //!
   TBranch        *b_pfnhits;   //!
   TBranch        *b_pfnplanehits_U;   //!
   TBranch        *b_pfnplanehits_V;   //!
   TBranch        *b_pfnplanehits_Y;   //!
   TBranch        *b_hits_u;   //!
   TBranch        *b_hits_v;   //!
   TBranch        *b_hits_y;   //!
   TBranch        *b_topological_score;   //!
   TBranch        *b_slclustfrac;   //!
   TBranch        *b_mc_pdg;   //!
   TBranch        *b_mc_E;   //!
   TBranch        *b_mc_vx;   //!
   TBranch        *b_mc_vy;   //!
   TBranch        *b_mc_vz;   //!
   TBranch        *b_mc_endx;   //!
   TBranch        *b_mc_endy;   //!
   TBranch        *b_mc_endz;   //!
   TBranch        *b_mc_px;   //!
   TBranch        *b_mc_py;   //!
   TBranch        *b_mc_pz;   //!
   TBranch        *b_mc_completeness;   //!
   TBranch        *b_mc_purity;   //!
   TBranch        *b_endmuonprocess;   //!
   TBranch        *b_endmuonmichel;   //!
   TBranch        *b_weights;   //!
   TBranch        *b_weightsFlux;   //!
   TBranch        *b_weightsGenie;   //!
   TBranch        *b_weightsReint;   //!
   TBranch        *b_weightSpline;   //!
   TBranch        *b_nu_flashmatch_score;   //!
   TBranch        *b_best_cosmic_flashmatch_score;   //!
   TBranch        *b_best_obviouscosmic_flashmatch_score;   //!
   TBranch        *b_cosmic_flashmatch_score_v;   //!
   TBranch        *b_NeutrinoEnergy0;   //!
   TBranch        *b_NeutrinoEnergy1;   //!
   TBranch        *b_NeutrinoEnergy2;   //!
   TBranch        *b_SliceCaloEnergy0;   //!
   TBranch        *b_SliceCaloEnergy1;   //!
   TBranch        *b_SliceCaloEnergy2;   //!
   TBranch        *b_gamma_edep;   //!
   TBranch        *b_gamma_etot;   //!
   TBranch        *b_gamma_dist;   //!
   TBranch        *b_gamma_parent;   //!
   TBranch        *b_elec_edep;   //!
   TBranch        *b_elec_etot;   //!
   TBranch        *b_elec_dist;   //!
   TBranch        *b_elec_parent;   //!
   TBranch        *b_gamma1_edep;   //!
   TBranch        *b_gamma1_etot;   //!
   TBranch        *b_gamma2_edep;   //!
   TBranch        *b_gamma2_etot;   //!
   TBranch        *b_nflag_pl1;   //!
   TBranch        *b_nnoise_pl1;   //!
   TBranch        *b_nslhits_pl1;   //!
   TBranch        *b_nslnoise_pl1;   //!
   TBranch        *b_nhits_pl1;   //!
   TBranch        *b_frac_slnoise_pl1;   //!
   TBranch        *b_shr_dedx_u_v;   //!
   TBranch        *b_shr_dedx_v_v;   //!
   TBranch        *b_shr_dedx_y_v;   //!
   TBranch        *b_shr_energy_u_v;   //!
   TBranch        *b_shr_energy_v_v;   //!
   TBranch        *b_shr_energy_y_v;   //!
   TBranch        *b_shr_pfp_id_v;   //!
   TBranch        *b_shr_start_x_v;   //!
   TBranch        *b_shr_start_y_v;   //!
   TBranch        *b_shr_start_z_v;   //!
   TBranch        *b_shr_start_U_v;   //!
   TBranch        *b_shr_start_V_v;   //!
   TBranch        *b_shr_dist_v;   //!
   TBranch        *b_shr_nclus_v;   //!
   TBranch        *b_shr_clushitfrac_v;   //!
   TBranch        *b_shr_px_v;   //!
   TBranch        *b_shr_py_v;   //!
   TBranch        *b_shr_pz_v;   //!
   TBranch        *b_shr_openangle_v;   //!
   TBranch        *b_shr_theta_v;   //!
   TBranch        *b_shr_phi_v;   //!
   TBranch        *b_shr_tkfit_nhits_v;   //!
   TBranch        *b_shr_tkfit_start_x_v;   //!
   TBranch        *b_shr_tkfit_start_y_v;   //!
   TBranch        *b_shr_tkfit_start_z_v;   //!
   TBranch        *b_shr_tkfit_start_U_v;   //!
   TBranch        *b_shr_tkfit_start_V_v;   //!
   TBranch        *b_shr_tkfit_theta_v;   //!
   TBranch        *b_shr_tkfit_phi_v;   //!
   TBranch        *b_shr_tkfit_dedx_u_v;   //!
   TBranch        *b_shr_tkfit_dedx_v_v;   //!
   TBranch        *b_shr_tkfit_dedx_y_v;   //!
   TBranch        *b_shr_tkfit_dedx_nhits_u_v;   //!
   TBranch        *b_shr_tkfit_dedx_nhits_v_v;   //!
   TBranch        *b_shr_tkfit_dedx_nhits_y_v;   //!
   TBranch        *b_shr_spacepoint_start_x_v;   //!
   TBranch        *b_shr_spacepoint_start_y_v;   //!
   TBranch        *b_shr_spacepoint_start_z_v;   //!
   TBranch        *b_shr_spacepoint_start_U_v;   //!
   TBranch        *b_shr_spacepoint_start_V_v;   //!
   TBranch        *b_shr_hits_start_U_wire_v;   //!
   TBranch        *b_shr_hits_start_U_x_v;   //!
   TBranch        *b_shr_hits_start_V_wire_v;   //!
   TBranch        *b_shr_hits_start_V_x_v;   //!
   TBranch        *b_shr_hits_start_Y_wire_v;   //!
   TBranch        *b_shr_hits_start_Y_x_v;   //!
   TBranch        *b_evnunhits;   //!
   TBranch        *b_evlepnhits;   //!
   TBranch        *b_evpronhits;   //!
   TBranch        *b_evpi1nhits;   //!
   TBranch        *b_evpi0nhits;   //!
   TBranch        *b_evneunhits;   //!
   TBranch        *b_evgamnhits;   //!
   TBranch        *b_evothnhits;   //!
   TBranch        *b_slnunhits;   //!
   TBranch        *b_sllepnhits;   //!
   TBranch        *b_slpronhits;   //!
   TBranch        *b_slpi1nhits;   //!
   TBranch        *b_slpi0nhits;   //!
   TBranch        *b_slneunhits;   //!
   TBranch        *b_slgamnhits;   //!
   TBranch        *b_slothnhits;   //!
   TBranch        *b_pfnunhits;   //!
   TBranch        *b_pflepnhits;   //!
   TBranch        *b_pfpronhits;   //!
   TBranch        *b_pfpi1nhits;   //!
   TBranch        *b_pfpi0nhits;   //!
   TBranch        *b_pfneunhits;   //!
   TBranch        *b_pfgamnhits;   //!
   TBranch        *b_pfothnhits;   //!
   TBranch        *b_trk_bragg_p_v;   //!
   TBranch        *b_trk_bragg_mu_v;   //!
   TBranch        *b_trk_bragg_mip_v;   //!
   TBranch        *b_trk_pida_v;   //!
   TBranch        *b_trk_pid_chipr_v;   //!
   TBranch        *b_trk_pid_chipi_v;   //!
   TBranch        *b_trk_pid_chika_v;   //!
   TBranch        *b_trk_pid_chimu_v;   //!
   TBranch        *b_trk_bragg_p_u_v;   //!
   TBranch        *b_trk_bragg_mu_u_v;   //!
   TBranch        *b_trk_bragg_mip_u_v;   //!
   TBranch        *b_trk_pida_u_v;   //!
   TBranch        *b_trk_pid_chipr_u_v;   //!
   TBranch        *b_trk_pid_chipi_u_v;   //!
   TBranch        *b_trk_pid_chika_u_v;   //!
   TBranch        *b_trk_pid_chimu_u_v;   //!
   TBranch        *b_trk_bragg_p_v_v;   //!
   TBranch        *b_trk_bragg_mu_v_v;   //!
   TBranch        *b_trk_bragg_mip_v_v;   //!
   TBranch        *b_trk_pida_v_v;   //!
   TBranch        *b_trk_pid_chipr_v_v;   //!
   TBranch        *b_trk_pid_chipi_v_v;   //!
   TBranch        *b_trk_pid_chika_v_v;   //!
   TBranch        *b_trk_pid_chimu_v_v;   //!
   TBranch        *b_trk_pfp_id_v;   //!
   TBranch        *b_trk_dir_x_v;   //!
   TBranch        *b_trk_dir_y_v;   //!
   TBranch        *b_trk_dir_z_v;   //!
   TBranch        *b_trk_start_x_v;   //!
   TBranch        *b_trk_start_y_v;   //!
   TBranch        *b_trk_start_z_v;   //!
   TBranch        *b_trk_end_x_v;   //!
   TBranch        *b_trk_end_y_v;   //!
   TBranch        *b_trk_end_z_v;   //!
   TBranch        *b_trk_distance_v;   //!
   TBranch        *b_trk_theta_v;   //!
   TBranch        *b_trk_phi_v;   //!
   TBranch        *b_trk_len_v;   //!
   TBranch        *b_trk_mcs_muon_mom_v;   //!
   TBranch        *b_trk_energy_proton_v;   //!
   TBranch        *b_trk_energy_muon_v;   //!
   TBranch        *b_bdt_nuNCpi0;   //!
   TBranch        *b_bdt_numuCCpi0;   //!
   TBranch        *b_bdt_numuCC;   //!
   TBranch        *b_bdt_ext;   //!
   TBranch        *b_bdt_cosmic;   //!
   TBranch        *b_bdt_global;   //!

   NeutrinoSelectionFilter(TString filename, string sample, string namedir);
   virtual ~NeutrinoSelectionFilter();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(string sample, float potweight, string namedir);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef NeutrinoSelectionFilter_cxx
NeutrinoSelectionFilter::NeutrinoSelectionFilter(TString filename, string sample, string namedir) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
	float potweight;
	cout << Form("%s:/nuselection",filename.Data()) << endl;
	TTree *tree = 0;
	TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(Form("%s",filename.Data()));
	if (!f || !f->IsOpen()) {
		f = new TFile(Form("%s",filename.Data()));
	}
	TDirectory * dir = (TDirectory*)f->Get(Form("%s:/nuselection",filename.Data()));
	if (sample == "nue")       potweight=0.000055;
	if (sample == "numu")      potweight=0.00733;
	if (sample == "databnb")   potweight=1;
	if (sample == "dirt")      potweight=0.0234;
	if (sample == "ext")       potweight=0.05355;
	
	cout << "dir = " << dir << endl;
	cout << Form("%s.root:nuselection",filename.Data()) << Form(" potweight = %f",potweight) << endl;
	dir->GetObject("NeutrinoSelectionFilter",tree);
	Init(tree);
	Loop(sample, potweight, namedir);

}

NeutrinoSelectionFilter::~NeutrinoSelectionFilter()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NeutrinoSelectionFilter::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NeutrinoSelectionFilter::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void NeutrinoSelectionFilter::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   dtrk_x_boundary = 0;
   dtrk_y_boundary = 0;
   dtrk_z_boundary = 0;
   dshr_x_boundary = 0;
   dshr_y_boundary = 0;
   dshr_z_boundary = 0;
   dvtx_x_boundary = 0;
   dvtx_y_boundary = 0;
   dvtx_z_boundary = 0;
   dtrk_boundary = 0;
   dvtx_boundary = 0;
   dshr_boundary = 0;
   dmc_boundary = 0;
   pfp_slice_idx = 0;
   backtracked_pdg = 0;
   backtracked_e = 0;
   backtracked_purity = 0;
   backtracked_completeness = 0;
   backtracked_overlay_purity = 0;
   backtracked_px = 0;
   backtracked_py = 0;
   backtracked_pz = 0;
   backtracked_start_x = 0;
   backtracked_start_y = 0;
   backtracked_start_z = 0;
   backtracked_start_t = 0;
   backtracked_start_U = 0;
   backtracked_start_V = 0;
   backtracked_start_Y = 0;
   backtracked_sce_start_x = 0;
   backtracked_sce_start_y = 0;
   backtracked_sce_start_z = 0;
   backtracked_sce_start_U = 0;
   backtracked_sce_start_V = 0;
   backtracked_sce_start_Y = 0;
   trk_score_v = 0;
   pfpdg = 0;
   pfnhits = 0;
   pfnplanehits_U = 0;
   pfnplanehits_V = 0;
   pfnplanehits_Y = 0;
   mc_pdg = 0;
   mc_E = 0;
   mc_vx = 0;
   mc_vy = 0;
   mc_vz = 0;
   mc_endx = 0;
   mc_endy = 0;
   mc_endz = 0;
   mc_px = 0;
   mc_py = 0;
   mc_pz = 0;
   mc_completeness = 0;
   mc_purity = 0;
   endmuonprocess = 0;
   weights = 0;
   weightsFlux = 0;
   weightsGenie = 0;
   weightsReint = 0;
   cosmic_flashmatch_score_v = 0;
   shr_dedx_u_v = 0;
   shr_dedx_v_v = 0;
   shr_dedx_y_v = 0;
   shr_energy_u_v = 0;
   shr_energy_v_v = 0;
   shr_energy_y_v = 0;
   shr_pfp_id_v = 0;
   shr_start_x_v = 0;
   shr_start_y_v = 0;
   shr_start_z_v = 0;
   shr_start_U_v = 0;
   shr_start_V_v = 0;
   shr_dist_v = 0;
   shr_nclus_v = 0;
   shr_clushitfrac_v = 0;
   shr_px_v = 0;
   shr_py_v = 0;
   shr_pz_v = 0;
   shr_openangle_v = 0;
   shr_theta_v = 0;
   shr_phi_v = 0;
   shr_tkfit_nhits_v = 0;
   shr_tkfit_start_x_v = 0;
   shr_tkfit_start_y_v = 0;
   shr_tkfit_start_z_v = 0;
   shr_tkfit_start_U_v = 0;
   shr_tkfit_start_V_v = 0;
   shr_tkfit_theta_v = 0;
   shr_tkfit_phi_v = 0;
   shr_tkfit_dedx_u_v = 0;
   shr_tkfit_dedx_v_v = 0;
   shr_tkfit_dedx_y_v = 0;
   shr_tkfit_dedx_nhits_u_v = 0;
   shr_tkfit_dedx_nhits_v_v = 0;
   shr_tkfit_dedx_nhits_y_v = 0;
   shr_spacepoint_start_x_v = 0;
   shr_spacepoint_start_y_v = 0;
   shr_spacepoint_start_z_v = 0;
   shr_spacepoint_start_U_v = 0;
   shr_spacepoint_start_V_v = 0;
   shr_hits_start_U_wire_v = 0;
   shr_hits_start_U_x_v = 0;
   shr_hits_start_V_wire_v = 0;
   shr_hits_start_V_x_v = 0;
   shr_hits_start_Y_wire_v = 0;
   shr_hits_start_Y_x_v = 0;
   pfnunhits = 0;
   pflepnhits = 0;
   pfpronhits = 0;
   pfpi1nhits = 0;
   pfpi0nhits = 0;
   pfneunhits = 0;
   pfgamnhits = 0;
   pfothnhits = 0;
   trk_bragg_p_v = 0;
   trk_bragg_mu_v = 0;
   trk_bragg_mip_v = 0;
   trk_pida_v = 0;
   trk_pid_chipr_v = 0;
   trk_pid_chipi_v = 0;
   trk_pid_chika_v = 0;
   trk_pid_chimu_v = 0;
   trk_bragg_p_u_v = 0;
   trk_bragg_mu_u_v = 0;
   trk_bragg_mip_u_v = 0;
   trk_pida_u_v = 0;
   trk_pid_chipr_u_v = 0;
   trk_pid_chipi_u_v = 0;
   trk_pid_chika_u_v = 0;
   trk_pid_chimu_u_v = 0;
   trk_bragg_p_v_v = 0;
   trk_bragg_mu_v_v = 0;
   trk_bragg_mip_v_v = 0;
   trk_pida_v_v = 0;
   trk_pid_chipr_v_v = 0;
   trk_pid_chipi_v_v = 0;
   trk_pid_chika_v_v = 0;
   trk_pid_chimu_v_v = 0;
   trk_pfp_id_v = 0;
   trk_dir_x_v = 0;
   trk_dir_y_v = 0;
   trk_dir_z_v = 0;
   trk_start_x_v = 0;
   trk_start_y_v = 0;
   trk_start_z_v = 0;
   trk_end_x_v = 0;
   trk_end_y_v = 0;
   trk_end_z_v = 0;
   trk_distance_v = 0;
   trk_theta_v = 0;
   trk_phi_v = 0;
   trk_len_v = 0;
   trk_mcs_muon_mom_v = 0;
   trk_energy_proton_v = 0;
   trk_energy_muon_v = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("selected", &selected, &b_selected);
   fChain->SetBranchAddress("trk_id", &trk_id, &b_trk_pfp_id);
   fChain->SetBranchAddress("shr_id", &shr_id, &b_shr_pfp_id);
   fChain->SetBranchAddress("shr_energy_tot", &shr_energy_tot, &b_shr_energy_tot);
   fChain->SetBranchAddress("shr_energy", &shr_energy, &b_shr_energy);
   fChain->SetBranchAddress("shr_energy_tot_cali", &shr_energy_tot_cali, &b_shr_energy_tot_cali);
   fChain->SetBranchAddress("shr_energy_cali", &shr_energy_cali, &b_shr_energy_cali);
   fChain->SetBranchAddress("shr_theta", &shr_theta, &b_shr_theta);
   fChain->SetBranchAddress("shr_phi", &shr_phi, &b_shr_phi);
   fChain->SetBranchAddress("shr_pca_0", &shr_pca_0, &b_shr_pca_0);
   fChain->SetBranchAddress("shr_pca_1", &shr_pca_1, &b_shr_pca_1);
   fChain->SetBranchAddress("shr_pca_2", &shr_pca_2, &b_shr_pca_2);
   fChain->SetBranchAddress("shr_px", &shr_px, &b_shr_px);
   fChain->SetBranchAddress("shr_py", &shr_py, &b_shr_py);
   fChain->SetBranchAddress("shr_pz", &shr_pz, &b_shr_pz);
   fChain->SetBranchAddress("shr_openangle", &shr_openangle, &b_shr_openangle);
   fChain->SetBranchAddress("shr_tkfit_start_x", &shr_tkfit_start_x, &b_shr_tkfit_start_x);
   fChain->SetBranchAddress("shr_tkfit_start_y", &shr_tkfit_start_y, &b_shr_tkfit_start_y);
   fChain->SetBranchAddress("shr_tkfit_start_z", &shr_tkfit_start_z, &b_shr_tkfit_start_z);
   fChain->SetBranchAddress("shr_tkfit_theta", &shr_tkfit_theta, &b_shr_tkfit_theta);
   fChain->SetBranchAddress("shr_tkfit_phi", &shr_tkfit_phi, &b_shr_tkfit_phi);
   fChain->SetBranchAddress("shr_start_x", &shr_start_x, &b_shr_start_x);
   fChain->SetBranchAddress("shr_start_y", &shr_start_y, &b_shr_start_y);
   fChain->SetBranchAddress("shr_start_z", &shr_start_z, &b_shr_start_z);
   fChain->SetBranchAddress("shr_dedx_Y", &shr_dedx_Y, &b_shr_dedx_Y);
   fChain->SetBranchAddress("shr_dedx_V", &shr_dedx_V, &b_shr_dedx_V);
   fChain->SetBranchAddress("shr_dedx_U", &shr_dedx_U, &b_shr_dedx_U);
   fChain->SetBranchAddress("shr_dedx_Y_cali", &shr_dedx_Y_cali, &b_shr_dedx_Y_cali);
   fChain->SetBranchAddress("shr_dedx_V_cali", &shr_dedx_V_cali, &b_shr_dedx_V_cali);
   fChain->SetBranchAddress("shr_dedx_U_cali", &shr_dedx_U_cali, &b_shr_dedx_U_cali);
   fChain->SetBranchAddress("shr_tkfit_dedx_Y", &shr_tkfit_dedx_Y, &b_shr_tkfit_dedx_Y);
   fChain->SetBranchAddress("shr_tkfit_dedx_V", &shr_tkfit_dedx_V, &b_shr_tkfit_dedx_V);
   fChain->SetBranchAddress("shr_tkfit_dedx_U", &shr_tkfit_dedx_U, &b_shr_tkfit_dedx_U);
   fChain->SetBranchAddress("shr_tkfit_nhits_Y", &shr_tkfit_nhits_Y, &b_shr_tkfit_nhits_Y);
   fChain->SetBranchAddress("shr_tkfit_nhits_V", &shr_tkfit_nhits_V, &b_shr_tkfit_nhits_V);
   fChain->SetBranchAddress("shr_tkfit_nhits_U", &shr_tkfit_nhits_U, &b_shr_tkfit_nhits_U);
   fChain->SetBranchAddress("shr_chipr", &shr_chipr, &b_shr_chipr);
   fChain->SetBranchAddress("shr_chimu", &shr_chimu, &b_shr_chimu);
   fChain->SetBranchAddress("shr_bragg_p", &shr_bragg_p, &b_shr_bragg_p);
   fChain->SetBranchAddress("shr_bragg_mu", &shr_bragg_mu, &b_shr_bragg_mu);
   fChain->SetBranchAddress("shr_bragg_mip", &shr_bragg_mip, &b_shr_bragg_mip);
   fChain->SetBranchAddress("shr_bragg_kaon", &shr_bragg_kaon, &b_shr_bragg_kaon);
   fChain->SetBranchAddress("shr_bragg_pion", &shr_bragg_pion, &b_shr_bragg_pion);
   fChain->SetBranchAddress("tksh_distance", &tksh_distance, &b_tksh_distance);
   fChain->SetBranchAddress("tksh_angle", &tksh_angle, &b_tksh_angle);
   fChain->SetBranchAddress("shr_distance", &shr_distance, &b_shr_distance);
   fChain->SetBranchAddress("shr_score", &shr_score, &b_shr_score);
   fChain->SetBranchAddress("shr_bkt_pdg", &shr_bkt_pdg, &b_shr_bkt_pdg);
   fChain->SetBranchAddress("shr_bkt_purity", &shr_bkt_purity, &b_shr_bkt_purity);
   fChain->SetBranchAddress("shr_bkt_completeness", &shr_bkt_completeness, &b_shr_bkt_completeness);
   fChain->SetBranchAddress("shr_bkt_E", &shr_bkt_E, &b_shr_bkt_E);
   fChain->SetBranchAddress("trk_len", &trk_len, &b_trk_len);
   fChain->SetBranchAddress("trk_theta", &trk_theta, &b_trk_theta);
   fChain->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
   fChain->SetBranchAddress("trk_energy", &trk_energy, &b_trk_energy);
   fChain->SetBranchAddress("trk_energy_muon", &trk_energy_muon, &b_trk_energy_muon);
   fChain->SetBranchAddress("trk_energy_tot", &trk_energy_tot, &b_trk_energy_tot);
   fChain->SetBranchAddress("trk_energy_muon_tot", &trk_energy_muon_tot, &b_trk_energy_muon_tot);
   fChain->SetBranchAddress("trk_distance", &trk_distance, &b_trk_distance);
   fChain->SetBranchAddress("trk_score", &trk_score, &b_trk_score);
   fChain->SetBranchAddress("trk_bkt_pdg", &trk_bkt_pdg, &b_trk_bkt_pdg);
   fChain->SetBranchAddress("trk_bkt_purity", &trk_bkt_purity, &b_trk_bkt_purity);
   fChain->SetBranchAddress("trk_bkt_completeness", &trk_bkt_completeness, &b_trk_bkt_completeness);
   fChain->SetBranchAddress("trk_bkt_E", &trk_bkt_E, &b_trk_bkt_E);
   fChain->SetBranchAddress("trk_chipr_best", &trk_chipr_best, &b_trk_chipr_best);
   fChain->SetBranchAddress("trk_chipr_worst", &trk_chipr_worst, &b_trk_chipr_worst);
   fChain->SetBranchAddress("trk_chimu_best", &trk_chimu_best, &b_trk_chimu_best);
   fChain->SetBranchAddress("trk_chimu_worst", &trk_chimu_worst, &b_trk_chimu_worst);
   fChain->SetBranchAddress("trk_chipr", &trk_chipr, &b_trk_chipr);
   fChain->SetBranchAddress("trk_chimu", &trk_chimu, &b_trk_chimu);
   fChain->SetBranchAddress("trk_pida", &trk_pida, &b_trk_pida);
   fChain->SetBranchAddress("trk_bragg_p", &trk_bragg_p, &b_trk_bragg_p);
   fChain->SetBranchAddress("trk_bragg_mu", &trk_bragg_mu, &b_trk_bragg_mu);
   fChain->SetBranchAddress("trk_bragg_mip", &trk_bragg_mip, &b_trk_bragg_mip);
   fChain->SetBranchAddress("trk_bragg_kaon", &trk_bragg_kaon, &b_trk_bragg_kaon);
   fChain->SetBranchAddress("trk_bragg_pion", &trk_bragg_pion, &b_trk_bragg_pion);
   fChain->SetBranchAddress("trk_hits_max", &trk_hits_max, &b_trk_hits_max);
   fChain->SetBranchAddress("shr_hits_max", &shr_hits_max, &b_shr_hits_max);
   fChain->SetBranchAddress("total_hits_y", &total_hits_y, &b_total_hits_y);
   fChain->SetBranchAddress("extra_energy_y", &extra_energy_y, &b_extra_energy_y);
   fChain->SetBranchAddress("trk_energy_hits_tot", &trk_energy_hits_tot, &b_trk_energy_hits_tot);
   fChain->SetBranchAddress("shr_hits_tot", &shr_hits_tot, &b_shr_hits_tot);
   fChain->SetBranchAddress("shr_hits_y_tot", &shr_hits_y_tot, &b_shr_hits_y_tot);
   fChain->SetBranchAddress("shr_hits_u_tot", &shr_hits_u_tot, &b_shr_hits_u_tot);
   fChain->SetBranchAddress("shr_hits_v_tot", &shr_hits_v_tot, &b_shr_hits_v_tot);
   fChain->SetBranchAddress("trk_hits_tot", &trk_hits_tot, &b_trk_hits_tot);
   fChain->SetBranchAddress("trk_hits_y_tot", &trk_hits_y_tot, &b_trk_hits_y_tot);
   fChain->SetBranchAddress("trk_hits_u_tot", &trk_hits_u_tot, &b_trk_hits_u_tot);
   fChain->SetBranchAddress("trk_hits_v_tot", &trk_hits_v_tot, &b_trk_hits_v_tot);
   fChain->SetBranchAddress("n_tracks_contained", &n_tracks_contained, &b_n_tracks_contained);
   fChain->SetBranchAddress("n_showers_contained", &n_showers_contained, &b_n_showers_contained);
   fChain->SetBranchAddress("matched_E", &matched_E, &b_matched_E);
   fChain->SetBranchAddress("hits_ratio", &hits_ratio, &b_hits_ratio);
   fChain->SetBranchAddress("contained_fraction", &contained_fraction, &b_contained_fraction);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("p", &p, &b_p);
   fChain->SetBranchAddress("pt_assume_muon", &pt_assume_muon, &b_pt_assume_muon);
   fChain->SetBranchAddress("p_assume_muon", &p_assume_muon, &b_p_assume_muon);
   fChain->SetBranchAddress("dvtx", &dvtx, &b_dvtx);
   fChain->SetBranchAddress("dtrk", &dtrk, &b_dtrk);
   fChain->SetBranchAddress("dtrk_x_boundary", &dtrk_x_boundary, &b_dtrk_x_boundary);
   fChain->SetBranchAddress("dtrk_y_boundary", &dtrk_y_boundary, &b_dtrk_y_boundary);
   fChain->SetBranchAddress("dtrk_z_boundary", &dtrk_z_boundary, &b_dtrk_z_boundary);
   fChain->SetBranchAddress("dshr_x_boundary", &dshr_x_boundary, &b_dshr_x_boundary);
   fChain->SetBranchAddress("dshr_y_boundary", &dshr_y_boundary, &b_dshr_y_boundary);
   fChain->SetBranchAddress("dshr_z_boundary", &dshr_z_boundary, &b_dshr_z_boundary);
   fChain->SetBranchAddress("dvtx_x_boundary", &dvtx_x_boundary, &b_dvtx_x_boundary);
   fChain->SetBranchAddress("dvtx_y_boundary", &dvtx_y_boundary, &b_dvtx_y_boundary);
   fChain->SetBranchAddress("dvtx_z_boundary", &dvtx_z_boundary, &b_dvtx_z_boundary);
   fChain->SetBranchAddress("dtrk_boundary", &dtrk_boundary, &b_dtrk_boundary);
   fChain->SetBranchAddress("dvtx_boundary", &dvtx_boundary, &b_dvtx_boundary);
   fChain->SetBranchAddress("dshr_boundary", &dshr_boundary, &b_dshr_boundary);
   fChain->SetBranchAddress("dmc_boundary", &dmc_boundary, &b_dmc_boundary);
   fChain->SetBranchAddress("CosmicIP", &CosmicIP, &b_CosmicIP);
   fChain->SetBranchAddress("_run", &_run, &b_run);
   fChain->SetBranchAddress("_sub", &_sub, &b_sub);
   fChain->SetBranchAddress("_evt", &_evt, &b_evt);
   fChain->SetBranchAddress("_nAssnCosmics", &_nAssnCosmics, &b__nAssnCosmics);
   fChain->SetBranchAddress("_kMaxCosm", &_kMaxCosm, &b__kMaxCosm);
   fChain->SetBranchAddress("_t0Cosm_v", _t0Cosm_v, &b__t0Cosm_v);
   fChain->SetBranchAddress("_t0timesVCosm_v", _t0timesVCosm_v, &b__t0timesVCosm_v);
   fChain->SetBranchAddress("_xStartCosm_v", _xStartCosm_v, &b__xStartCosm_v);
   fChain->SetBranchAddress("_yStartCosm_v", _yStartCosm_v, &b__yStartCosm_v);
   fChain->SetBranchAddress("_zStartCosm_v", _zStartCosm_v, &b__zStartCosm_v);
   fChain->SetBranchAddress("_xEndCosm_v", _xEndCosm_v, &b__xEndCosm_v);
   fChain->SetBranchAddress("_yEndCosm_v", _yEndCosm_v, &b__yEndCosm_v);
   fChain->SetBranchAddress("_zEndCosm_v", _zEndCosm_v, &b__zEndCosm_v);
   fChain->SetBranchAddress("_recoNu_vtx_x", &_recoNu_vtx_x, &b_reco_vtx_x);
   fChain->SetBranchAddress("_recoNu_vtx_y", &_recoNu_vtx_y, &b_reco_vtx_y);
   fChain->SetBranchAddress("_recoNu_vtx_z", &_recoNu_vtx_z, &b_reco_vtx_z);
   fChain->SetBranchAddress("_t0_nu_cosmic", &_t0_nu_cosmic, &b__t0_nu_cosmic);
   fChain->SetBranchAddress("_nu_cosmic_x", &_nu_cosmic_x, &b__nu_cosmic_x);
   fChain->SetBranchAddress("_nu_cosmic_y", &_nu_cosmic_y, &b__nu_cosmic_y);
   fChain->SetBranchAddress("_nu_cosmic_z", &_nu_cosmic_z, &b__nu_cosmic_z);
   fChain->SetBranchAddress("_nu_cosmic_Length", &_nu_cosmic_Length, &b__nu_cosmic_Length);
   fChain->SetBranchAddress("_nu_cosmic_Start_x", &_nu_cosmic_Start_x, &b__nu_cosmic_Start_x);
   fChain->SetBranchAddress("_nu_cosmic_Start_y", &_nu_cosmic_Start_y, &b__nu_cosmic_Start_y);
   fChain->SetBranchAddress("_nu_cosmic_Start_z", &_nu_cosmic_Start_z, &b__nu_cosmic_Start_z);
   fChain->SetBranchAddress("_nu_cosmic_End_x", &_nu_cosmic_End_x, &b__nu_cosmic_End_x);
   fChain->SetBranchAddress("_nu_cosmic_End_y", &_nu_cosmic_End_y, &b__nu_cosmic_End_y);
   fChain->SetBranchAddress("_nu_cosmic_End_z", &_nu_cosmic_End_z, &b__nu_cosmic_End_z);
   fChain->SetBranchAddress("_nu_cosmic_TrackID", &_nu_cosmic_TrackID, &b__nu_cosmic_TrackID);
   fChain->SetBranchAddress("_closestNuCosmicDist", &_closestNuCosmicDist, &b__closestNuCosmicDist);
   fChain->SetBranchAddress("_rand_vtx_x", &_rand_vtx_x, &b_rand_vtx_x);
   fChain->SetBranchAddress("_rand_vtx_y", &_rand_vtx_y, &b_rand_vtx_y);
   fChain->SetBranchAddress("_rand_vtx_z", &_rand_vtx_z, &b_rand_vtx_z);
   fChain->SetBranchAddress("_t0_rand_cosmic", &_t0_rand_cosmic, &b__t0_rand_cosmic);
   fChain->SetBranchAddress("_rand_cosmic_x", &_rand_cosmic_x, &b__rand_cosmic_x);
   fChain->SetBranchAddress("_rand_cosmic_y", &_rand_cosmic_y, &b__rand_cosmic_y);
   fChain->SetBranchAddress("_rand_cosmic_z", &_rand_cosmic_z, &b__rand_cosmic_z);
   fChain->SetBranchAddress("_rand_cosmic_Length", &_rand_cosmic_Length, &b__rand_cosmic_Length);
   fChain->SetBranchAddress("_rand_cosmic_Start_x", &_rand_cosmic_Start_x, &b__rand_cosmic_Start_x);
   fChain->SetBranchAddress("_rand_cosmic_Start_y", &_rand_cosmic_Start_y, &b__rand_cosmic_Start_y);
   fChain->SetBranchAddress("_rand_cosmic_Start_z", &_rand_cosmic_Start_z, &b__rand_cosmic_Start_z);
   fChain->SetBranchAddress("_rand_cosmic_End_x", &_rand_cosmic_End_x, &b__rand_cosmic_End_x);
   fChain->SetBranchAddress("_rand_cosmic_End_y", &_rand_cosmic_End_y, &b__rand_cosmic_End_y);
   fChain->SetBranchAddress("_rand_cosmic_End_z", &_rand_cosmic_End_z, &b__rand_cosmic_End_z);
   fChain->SetBranchAddress("_rand_cosmic_TrackID", &_rand_cosmic_TrackID, &b__rand_cosmic_TrackID);
   fChain->SetBranchAddress("_closestRandCosmicDist", &_closestRandCosmicDist, &b__closestRandCosmicDist);
   fChain->SetBranchAddress("leeweight", &leeweight, &b_leeweight);
   fChain->SetBranchAddress("nu_vtx_x", &nu_vtx_x, &b_nu_vtx_x);
   fChain->SetBranchAddress("nu_vtx_y", &nu_vtx_y, &b_nu_vtx_y);
   fChain->SetBranchAddress("nu_vtx_z", &nu_vtx_z, &b_nu_vtx_z);
   fChain->SetBranchAddress("nu_sce_x", &nu_sce_x, &b_nu_sce_x);
   fChain->SetBranchAddress("nu_sce_y", &nu_sce_y, &b_nu_sce_y);
   fChain->SetBranchAddress("nu_sce_z", &nu_sce_z, &b_nu_sce_z);
   fChain->SetBranchAddress("true_pt", &true_pt, &b_true_pt);
   fChain->SetBranchAddress("true_pt_visible", &true_pt_visible, &b_true_pt_visible);
   fChain->SetBranchAddress("true_p", &true_p, &b_true_p);
   fChain->SetBranchAddress("true_p_visible", &true_p_visible, &b_true_p_visible);
   fChain->SetBranchAddress("true_e_visible", &true_e_visible, &b_true_e_visible);
   fChain->SetBranchAddress("nu_pdg", &nu_pdg, &b_nu_pdg);
   fChain->SetBranchAddress("ccnc", &ccnc, &b_ccnc);
   fChain->SetBranchAddress("interaction", &interaction, &b_interaction);
   fChain->SetBranchAddress("nu_e", &nu_e, &b_nu_e);
   fChain->SetBranchAddress("nu_pt", &nu_pt, &b_nu_pt);
   fChain->SetBranchAddress("theta", &theta, &b_theta);
   fChain->SetBranchAddress("vtx_x", &vtx_x, &b_vtx_x);
   fChain->SetBranchAddress("vtx_y", &vtx_y, &b_vtx_y);
   fChain->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z);
   fChain->SetBranchAddress("isVtxInActive", &isVtxInActive, &b_isVtxInActive);
   fChain->SetBranchAddress("isVtxInFiducial", &isVtxInFiducial, &b_isVtxInFiducial);
   fChain->SetBranchAddress("true_nu_vtx_t", &true_nu_vtx_t, &b_true_nu_vtx_t);
   fChain->SetBranchAddress("true_nu_vtx_x", &true_nu_vtx_x, &b_true_nu_vtx_x);
   fChain->SetBranchAddress("true_nu_vtx_y", &true_nu_vtx_y, &b_true_nu_vtx_y);
   fChain->SetBranchAddress("true_nu_vtx_z", &true_nu_vtx_z, &b_true_nu_vtx_z);
   fChain->SetBranchAddress("true_nu_vtx_sce_x", &true_nu_vtx_sce_x, &b_true_nu_vtx_sce_x);
   fChain->SetBranchAddress("true_nu_vtx_sce_y", &true_nu_vtx_sce_y, &b_true_nu_vtx_sce_y);
   fChain->SetBranchAddress("true_nu_vtx_sce_z", &true_nu_vtx_sce_z, &b_true_nu_vtx_sce_z);
   fChain->SetBranchAddress("reco_nu_vtx_x", &reco_nu_vtx_x, &b_reco_nu_vtx_x);
   fChain->SetBranchAddress("reco_nu_vtx_y", &reco_nu_vtx_y, &b_reco_nu_vtx_y);
   fChain->SetBranchAddress("reco_nu_vtx_z", &reco_nu_vtx_z, &b_reco_nu_vtx_z);
   fChain->SetBranchAddress("reco_nu_vtx_sce_x", &reco_nu_vtx_sce_x, &b_reco_nu_vtx_sce_x);
   fChain->SetBranchAddress("reco_nu_vtx_sce_y", &reco_nu_vtx_sce_y, &b_reco_nu_vtx_sce_y);
   fChain->SetBranchAddress("reco_nu_vtx_sce_z", &reco_nu_vtx_sce_z, &b_reco_nu_vtx_sce_z);
   fChain->SetBranchAddress("nmuon", &nmuon, &b_nmuon);
   fChain->SetBranchAddress("muon_e", &muon_e, &b_muon_e);
   fChain->SetBranchAddress("muon_c", &muon_c, &b_muon_c);
   fChain->SetBranchAddress("muon_p", &muon_p, &b_muon_p);
   fChain->SetBranchAddress("nelec", &nelec, &b_nelec);
   fChain->SetBranchAddress("elec_e", &elec_e, &b_elec_e);
   fChain->SetBranchAddress("elec_c", &elec_c, &b_elec_c);
   fChain->SetBranchAddress("elec_p", &elec_p, &b_elec_p);
   fChain->SetBranchAddress("elec_vx", &elec_vx, &b_elec_vx);
   fChain->SetBranchAddress("elec_vy", &elec_vy, &b_elec_vy);
   fChain->SetBranchAddress("elec_vz", &elec_vz, &b_elec_vz);
   fChain->SetBranchAddress("npi0", &npi0, &b_npi0);
   fChain->SetBranchAddress("pi0_e", &pi0_e, &b_pi0_e);
   fChain->SetBranchAddress("pi0_c", &pi0_c, &b_pi0_c);
   fChain->SetBranchAddress("pi0_p", &pi0_p, &b_pi0_p);
   fChain->SetBranchAddress("nneutron", &nneutron, &b_nneutron);
   fChain->SetBranchAddress("nproton", &nproton, &b_nproton);
   fChain->SetBranchAddress("proton_e", &proton_e, &b_proton_e);
   fChain->SetBranchAddress("proton_c", &proton_c, &b_proton_c);
   fChain->SetBranchAddress("proton_p", &proton_p, &b_proton_p);
   fChain->SetBranchAddress("npion", &npion, &b_npion);
   fChain->SetBranchAddress("pion_e", &pion_e, &b_pion_e);
   fChain->SetBranchAddress("pion_c", &pion_c, &b_pion_c);
   fChain->SetBranchAddress("pion_p", &pion_p, &b_pion_p);
   fChain->SetBranchAddress("nslice", &nslice, &b_nslice);
   fChain->SetBranchAddress("crtveto", &crtveto, &b_crtveto);
   fChain->SetBranchAddress("crthitpe", &crthitpe, &b_crthitpe);
   fChain->SetBranchAddress("pfp_slice_idx", &pfp_slice_idx, &b_pfp_slice_idx);
   fChain->SetBranchAddress("category", &category, &b_category);
   fChain->SetBranchAddress("backtracked_pdg", &backtracked_pdg, &b_backtracked_pdg);
   fChain->SetBranchAddress("backtracked_e", &backtracked_e, &b_backtracked_e);
   fChain->SetBranchAddress("backtracked_purity", &backtracked_purity, &b_backtracked_purity);
   fChain->SetBranchAddress("backtracked_completeness", &backtracked_completeness, &b_backtracked_completeness);
   fChain->SetBranchAddress("backtracked_overlay_purity", &backtracked_overlay_purity, &b_backtracked_overlay_purity);
   fChain->SetBranchAddress("backtracked_px", &backtracked_px, &b_backtracked_px);
   fChain->SetBranchAddress("backtracked_py", &backtracked_py, &b_backtracked_py);
   fChain->SetBranchAddress("backtracked_pz", &backtracked_pz, &b_backtracked_pz);
   fChain->SetBranchAddress("backtracked_start_x", &backtracked_start_x, &b_backtracked_start_x);
   fChain->SetBranchAddress("backtracked_start_y", &backtracked_start_y, &b_backtracked_start_y);
   fChain->SetBranchAddress("backtracked_start_z", &backtracked_start_z, &b_backtracked_start_z);
   fChain->SetBranchAddress("backtracked_start_t", &backtracked_start_t, &b_backtracked_start_t);
   fChain->SetBranchAddress("backtracked_start_U", &backtracked_start_U, &b_backtracked_start_U);
   fChain->SetBranchAddress("backtracked_start_V", &backtracked_start_V, &b_backtracked_start_V);
   fChain->SetBranchAddress("backtracked_start_Y", &backtracked_start_Y, &b_backtracked_start_Y);
   fChain->SetBranchAddress("backtracked_sce_start_x", &backtracked_sce_start_x, &b_backtracked_sce_start_x);
   fChain->SetBranchAddress("backtracked_sce_start_y", &backtracked_sce_start_y, &b_backtracked_sce_start_y);
   fChain->SetBranchAddress("backtracked_sce_start_z", &backtracked_sce_start_z, &b_backtracked_sce_start_z);
   fChain->SetBranchAddress("backtracked_sce_start_U", &backtracked_sce_start_U, &b_backtracked_sce_start_U);
   fChain->SetBranchAddress("backtracked_sce_start_V", &backtracked_sce_start_V, &b_backtracked_sce_start_V);
   fChain->SetBranchAddress("backtracked_sce_start_Y", &backtracked_sce_start_Y, &b_backtracked_sce_start_Y);
   fChain->SetBranchAddress("lep_e", &lep_e, &b_lep_e);
   fChain->SetBranchAddress("pass", &pass, &b_pass);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("sub", &sub, &b_sub);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("swtrig", &swtrig, &b_swtrig);
   fChain->SetBranchAddress("xtimeoffset", &xtimeoffset, &b_xtimeoffset);
   fChain->SetBranchAddress("xsceoffset", &xsceoffset, &b_xsceoffset);
   fChain->SetBranchAddress("ysceoffset", &ysceoffset, &b_ysceoffset);
   fChain->SetBranchAddress("zsceoffset", &zsceoffset, &b_zsceoffset);
   fChain->SetBranchAddress("evnhits", &evnhits, &b_evnhits);
   fChain->SetBranchAddress("slpdg", &slpdg, &b_slpdg);
   fChain->SetBranchAddress("slnhits", &slnhits, &b_slnhits);
   fChain->SetBranchAddress("n_pfps", &n_pfps, &b_n_pfps);
   fChain->SetBranchAddress("n_tracks", &n_tracks, &b_n_tracks);
   fChain->SetBranchAddress("n_showers", &n_showers, &b_n_showers);
   fChain->SetBranchAddress("trk_score_v", &trk_score_v, &b_trk_score_v);
   fChain->SetBranchAddress("pfpdg", &pfpdg, &b_pfpdg);
   fChain->SetBranchAddress("pfnhits", &pfnhits, &b_pfnhits);
   fChain->SetBranchAddress("pfnplanehits_U", &pfnplanehits_U, &b_pfnplanehits_U);
   fChain->SetBranchAddress("pfnplanehits_V", &pfnplanehits_V, &b_pfnplanehits_V);
   fChain->SetBranchAddress("pfnplanehits_Y", &pfnplanehits_Y, &b_pfnplanehits_Y);
   fChain->SetBranchAddress("hits_u", &hits_u, &b_hits_u);
   fChain->SetBranchAddress("hits_v", &hits_v, &b_hits_v);
   fChain->SetBranchAddress("hits_y", &hits_y, &b_hits_y);
   fChain->SetBranchAddress("topological_score", &topological_score, &b_topological_score);
   fChain->SetBranchAddress("slclustfrac", &slclustfrac, &b_slclustfrac);
   fChain->SetBranchAddress("mc_pdg", &mc_pdg, &b_mc_pdg);
   fChain->SetBranchAddress("mc_E", &mc_E, &b_mc_E);
   fChain->SetBranchAddress("mc_vx", &mc_vx, &b_mc_vx);
   fChain->SetBranchAddress("mc_vy", &mc_vy, &b_mc_vy);
   fChain->SetBranchAddress("mc_vz", &mc_vz, &b_mc_vz);
   fChain->SetBranchAddress("mc_endx", &mc_endx, &b_mc_endx);
   fChain->SetBranchAddress("mc_endy", &mc_endy, &b_mc_endy);
   fChain->SetBranchAddress("mc_endz", &mc_endz, &b_mc_endz);
   fChain->SetBranchAddress("mc_px", &mc_px, &b_mc_px);
   fChain->SetBranchAddress("mc_py", &mc_py, &b_mc_py);
   fChain->SetBranchAddress("mc_pz", &mc_pz, &b_mc_pz);
   fChain->SetBranchAddress("mc_completeness", &mc_completeness, &b_mc_completeness);
   fChain->SetBranchAddress("mc_purity", &mc_purity, &b_mc_purity);
   fChain->SetBranchAddress("endmuonprocess", &endmuonprocess, &b_endmuonprocess);
   fChain->SetBranchAddress("endmuonmichel", &endmuonmichel, &b_endmuonmichel);
   fChain->SetBranchAddress("weights", &weights, &b_weights);
   fChain->SetBranchAddress("weightsFlux", &weightsFlux, &b_weightsFlux);
   fChain->SetBranchAddress("weightsGenie", &weightsGenie, &b_weightsGenie);
   fChain->SetBranchAddress("weightsReint", &weightsReint, &b_weightsReint);
   fChain->SetBranchAddress("weightSpline", &weightSpline, &b_weightSpline);
   fChain->SetBranchAddress("nu_flashmatch_score", &nu_flashmatch_score, &b_nu_flashmatch_score);
   fChain->SetBranchAddress("best_cosmic_flashmatch_score", &best_cosmic_flashmatch_score, &b_best_cosmic_flashmatch_score);
   fChain->SetBranchAddress("best_obviouscosmic_flashmatch_score", &best_obviouscosmic_flashmatch_score, &b_best_obviouscosmic_flashmatch_score);
   fChain->SetBranchAddress("cosmic_flashmatch_score_v", &cosmic_flashmatch_score_v, &b_cosmic_flashmatch_score_v);
   fChain->SetBranchAddress("NeutrinoEnergy0", &NeutrinoEnergy0, &b_NeutrinoEnergy0);
   fChain->SetBranchAddress("NeutrinoEnergy1", &NeutrinoEnergy1, &b_NeutrinoEnergy1);
   fChain->SetBranchAddress("NeutrinoEnergy2", &NeutrinoEnergy2, &b_NeutrinoEnergy2);
   fChain->SetBranchAddress("SliceCaloEnergy0", &SliceCaloEnergy0, &b_SliceCaloEnergy0);
   fChain->SetBranchAddress("SliceCaloEnergy1", &SliceCaloEnergy1, &b_SliceCaloEnergy1);
   fChain->SetBranchAddress("SliceCaloEnergy2", &SliceCaloEnergy2, &b_SliceCaloEnergy2);
   fChain->SetBranchAddress("gamma_edep", &gamma_edep, &b_gamma_edep);
   fChain->SetBranchAddress("gamma_etot", &gamma_etot, &b_gamma_etot);
   fChain->SetBranchAddress("gamma_dist", &gamma_dist, &b_gamma_dist);
   fChain->SetBranchAddress("gamma_parent", &gamma_parent, &b_gamma_parent);
   fChain->SetBranchAddress("elec_edep", &elec_edep, &b_elec_edep);
   fChain->SetBranchAddress("elec_etot", &elec_etot, &b_elec_etot);
   fChain->SetBranchAddress("elec_dist", &elec_dist, &b_elec_dist);
   fChain->SetBranchAddress("elec_parent", &elec_parent, &b_elec_parent);
   fChain->SetBranchAddress("gamma1_edep", &gamma1_edep, &b_gamma1_edep);
   fChain->SetBranchAddress("gamma1_etot", &gamma1_etot, &b_gamma1_etot);
   fChain->SetBranchAddress("gamma2_edep", &gamma2_edep, &b_gamma2_edep);
   fChain->SetBranchAddress("gamma2_etot", &gamma2_etot, &b_gamma2_etot);
   fChain->SetBranchAddress("nflag_pl1", &nflag_pl1, &b_nflag_pl1);
   fChain->SetBranchAddress("nnoise_pl1", &nnoise_pl1, &b_nnoise_pl1);
   fChain->SetBranchAddress("nslhits_pl1", &nslhits_pl1, &b_nslhits_pl1);
   fChain->SetBranchAddress("nslnoise_pl1", &nslnoise_pl1, &b_nslnoise_pl1);
   fChain->SetBranchAddress("nhits_pl1", &nhits_pl1, &b_nhits_pl1);
   fChain->SetBranchAddress("frac_slnoise_pl1", &frac_slnoise_pl1, &b_frac_slnoise_pl1);
   fChain->SetBranchAddress("shr_dedx_u_v", &shr_dedx_u_v, &b_shr_dedx_u_v);
   fChain->SetBranchAddress("shr_dedx_v_v", &shr_dedx_v_v, &b_shr_dedx_v_v);
   fChain->SetBranchAddress("shr_dedx_y_v", &shr_dedx_y_v, &b_shr_dedx_y_v);
   fChain->SetBranchAddress("shr_energy_u_v", &shr_energy_u_v, &b_shr_energy_u_v);
   fChain->SetBranchAddress("shr_energy_v_v", &shr_energy_v_v, &b_shr_energy_v_v);
   fChain->SetBranchAddress("shr_energy_y_v", &shr_energy_y_v, &b_shr_energy_y_v);
   fChain->SetBranchAddress("shr_pfp_id_v", &shr_pfp_id_v, &b_shr_pfp_id_v);
   fChain->SetBranchAddress("shr_start_x_v", &shr_start_x_v, &b_shr_start_x_v);
   fChain->SetBranchAddress("shr_start_y_v", &shr_start_y_v, &b_shr_start_y_v);
   fChain->SetBranchAddress("shr_start_z_v", &shr_start_z_v, &b_shr_start_z_v);
   fChain->SetBranchAddress("shr_start_U_v", &shr_start_U_v, &b_shr_start_U_v);
   fChain->SetBranchAddress("shr_start_V_v", &shr_start_V_v, &b_shr_start_V_v);
   fChain->SetBranchAddress("shr_dist_v", &shr_dist_v, &b_shr_dist_v);
   fChain->SetBranchAddress("shr_nclus_v", &shr_nclus_v, &b_shr_nclus_v);
   fChain->SetBranchAddress("shr_clushitfrac_v", &shr_clushitfrac_v, &b_shr_clushitfrac_v);
   fChain->SetBranchAddress("shr_px_v", &shr_px_v, &b_shr_px_v);
   fChain->SetBranchAddress("shr_py_v", &shr_py_v, &b_shr_py_v);
   fChain->SetBranchAddress("shr_pz_v", &shr_pz_v, &b_shr_pz_v);
   fChain->SetBranchAddress("shr_openangle_v", &shr_openangle_v, &b_shr_openangle_v);
   fChain->SetBranchAddress("shr_theta_v", &shr_theta_v, &b_shr_theta_v);
   fChain->SetBranchAddress("shr_phi_v", &shr_phi_v, &b_shr_phi_v);
   fChain->SetBranchAddress("shr_tkfit_nhits_v", &shr_tkfit_nhits_v, &b_shr_tkfit_nhits_v);
   fChain->SetBranchAddress("shr_tkfit_start_x_v", &shr_tkfit_start_x_v, &b_shr_tkfit_start_x_v);
   fChain->SetBranchAddress("shr_tkfit_start_y_v", &shr_tkfit_start_y_v, &b_shr_tkfit_start_y_v);
   fChain->SetBranchAddress("shr_tkfit_start_z_v", &shr_tkfit_start_z_v, &b_shr_tkfit_start_z_v);
   fChain->SetBranchAddress("shr_tkfit_start_U_v", &shr_tkfit_start_U_v, &b_shr_tkfit_start_U_v);
   fChain->SetBranchAddress("shr_tkfit_start_V_v", &shr_tkfit_start_V_v, &b_shr_tkfit_start_V_v);
   fChain->SetBranchAddress("shr_tkfit_theta_v", &shr_tkfit_theta_v, &b_shr_tkfit_theta_v);
   fChain->SetBranchAddress("shr_tkfit_phi_v", &shr_tkfit_phi_v, &b_shr_tkfit_phi_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_u_v", &shr_tkfit_dedx_u_v, &b_shr_tkfit_dedx_u_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_v_v", &shr_tkfit_dedx_v_v, &b_shr_tkfit_dedx_v_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_y_v", &shr_tkfit_dedx_y_v, &b_shr_tkfit_dedx_y_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_nhits_u_v", &shr_tkfit_dedx_nhits_u_v, &b_shr_tkfit_dedx_nhits_u_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_nhits_v_v", &shr_tkfit_dedx_nhits_v_v, &b_shr_tkfit_dedx_nhits_v_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_nhits_y_v", &shr_tkfit_dedx_nhits_y_v, &b_shr_tkfit_dedx_nhits_y_v);
   fChain->SetBranchAddress("shr_spacepoint_start_x_v", &shr_spacepoint_start_x_v, &b_shr_spacepoint_start_x_v);
   fChain->SetBranchAddress("shr_spacepoint_start_y_v", &shr_spacepoint_start_y_v, &b_shr_spacepoint_start_y_v);
   fChain->SetBranchAddress("shr_spacepoint_start_z_v", &shr_spacepoint_start_z_v, &b_shr_spacepoint_start_z_v);
   fChain->SetBranchAddress("shr_spacepoint_start_U_v", &shr_spacepoint_start_U_v, &b_shr_spacepoint_start_U_v);
   fChain->SetBranchAddress("shr_spacepoint_start_V_v", &shr_spacepoint_start_V_v, &b_shr_spacepoint_start_V_v);
   fChain->SetBranchAddress("shr_hits_start_U_wire_v", &shr_hits_start_U_wire_v, &b_shr_hits_start_U_wire_v);
   fChain->SetBranchAddress("shr_hits_start_U_x_v", &shr_hits_start_U_x_v, &b_shr_hits_start_U_x_v);
   fChain->SetBranchAddress("shr_hits_start_V_wire_v", &shr_hits_start_V_wire_v, &b_shr_hits_start_V_wire_v);
   fChain->SetBranchAddress("shr_hits_start_V_x_v", &shr_hits_start_V_x_v, &b_shr_hits_start_V_x_v);
   fChain->SetBranchAddress("shr_hits_start_Y_wire_v", &shr_hits_start_Y_wire_v, &b_shr_hits_start_Y_wire_v);
   fChain->SetBranchAddress("shr_hits_start_Y_x_v", &shr_hits_start_Y_x_v, &b_shr_hits_start_Y_x_v);
   fChain->SetBranchAddress("evnunhits", &evnunhits, &b_evnunhits);
   fChain->SetBranchAddress("evlepnhits", &evlepnhits, &b_evlepnhits);
   fChain->SetBranchAddress("evpronhits", &evpronhits, &b_evpronhits);
   fChain->SetBranchAddress("evpi1nhits", &evpi1nhits, &b_evpi1nhits);
   fChain->SetBranchAddress("evpi0nhits", &evpi0nhits, &b_evpi0nhits);
   fChain->SetBranchAddress("evneunhits", &evneunhits, &b_evneunhits);
   fChain->SetBranchAddress("evgamnhits", &evgamnhits, &b_evgamnhits);
   fChain->SetBranchAddress("evothnhits", &evothnhits, &b_evothnhits);
   fChain->SetBranchAddress("slnunhits", &slnunhits, &b_slnunhits);
   fChain->SetBranchAddress("sllepnhits", &sllepnhits, &b_sllepnhits);
   fChain->SetBranchAddress("slpronhits", &slpronhits, &b_slpronhits);
   fChain->SetBranchAddress("slpi1nhits", &slpi1nhits, &b_slpi1nhits);
   fChain->SetBranchAddress("slpi0nhits", &slpi0nhits, &b_slpi0nhits);
   fChain->SetBranchAddress("slneunhits", &slneunhits, &b_slneunhits);
   fChain->SetBranchAddress("slgamnhits", &slgamnhits, &b_slgamnhits);
   fChain->SetBranchAddress("slothnhits", &slothnhits, &b_slothnhits);
   fChain->SetBranchAddress("pfnunhits", &pfnunhits, &b_pfnunhits);
   fChain->SetBranchAddress("pflepnhits", &pflepnhits, &b_pflepnhits);
   fChain->SetBranchAddress("pfpronhits", &pfpronhits, &b_pfpronhits);
   fChain->SetBranchAddress("pfpi1nhits", &pfpi1nhits, &b_pfpi1nhits);
   fChain->SetBranchAddress("pfpi0nhits", &pfpi0nhits, &b_pfpi0nhits);
   fChain->SetBranchAddress("pfneunhits", &pfneunhits, &b_pfneunhits);
   fChain->SetBranchAddress("pfgamnhits", &pfgamnhits, &b_pfgamnhits);
   fChain->SetBranchAddress("pfothnhits", &pfothnhits, &b_pfothnhits);
   fChain->SetBranchAddress("trk_bragg_p_v", &trk_bragg_p_v, &b_trk_bragg_p_v);
   fChain->SetBranchAddress("trk_bragg_mu_v", &trk_bragg_mu_v, &b_trk_bragg_mu_v);
   fChain->SetBranchAddress("trk_bragg_mip_v", &trk_bragg_mip_v, &b_trk_bragg_mip_v);
   fChain->SetBranchAddress("trk_pida_v", &trk_pida_v, &b_trk_pida_v);
   fChain->SetBranchAddress("trk_pid_chipr_v", &trk_pid_chipr_v, &b_trk_pid_chipr_v);
   fChain->SetBranchAddress("trk_pid_chipi_v", &trk_pid_chipi_v, &b_trk_pid_chipi_v);
   fChain->SetBranchAddress("trk_pid_chika_v", &trk_pid_chika_v, &b_trk_pid_chika_v);
   fChain->SetBranchAddress("trk_pid_chimu_v", &trk_pid_chimu_v, &b_trk_pid_chimu_v);
   fChain->SetBranchAddress("trk_bragg_p_u_v", &trk_bragg_p_u_v, &b_trk_bragg_p_u_v);
   fChain->SetBranchAddress("trk_bragg_mu_u_v", &trk_bragg_mu_u_v, &b_trk_bragg_mu_u_v);
   fChain->SetBranchAddress("trk_bragg_mip_u_v", &trk_bragg_mip_u_v, &b_trk_bragg_mip_u_v);
   fChain->SetBranchAddress("trk_pida_u_v", &trk_pida_u_v, &b_trk_pida_u_v);
   fChain->SetBranchAddress("trk_pid_chipr_u_v", &trk_pid_chipr_u_v, &b_trk_pid_chipr_u_v);
   fChain->SetBranchAddress("trk_pid_chipi_u_v", &trk_pid_chipi_u_v, &b_trk_pid_chipi_u_v);
   fChain->SetBranchAddress("trk_pid_chika_u_v", &trk_pid_chika_u_v, &b_trk_pid_chika_u_v);
   fChain->SetBranchAddress("trk_pid_chimu_u_v", &trk_pid_chimu_u_v, &b_trk_pid_chimu_u_v);
   fChain->SetBranchAddress("trk_bragg_p_v_v", &trk_bragg_p_v_v, &b_trk_bragg_p_v_v);
   fChain->SetBranchAddress("trk_bragg_mu_v_v", &trk_bragg_mu_v_v, &b_trk_bragg_mu_v_v);
   fChain->SetBranchAddress("trk_bragg_mip_v_v", &trk_bragg_mip_v_v, &b_trk_bragg_mip_v_v);
   fChain->SetBranchAddress("trk_pida_v_v", &trk_pida_v_v, &b_trk_pida_v_v);
   fChain->SetBranchAddress("trk_pid_chipr_v_v", &trk_pid_chipr_v_v, &b_trk_pid_chipr_v_v);
   fChain->SetBranchAddress("trk_pid_chipi_v_v", &trk_pid_chipi_v_v, &b_trk_pid_chipi_v_v);
   fChain->SetBranchAddress("trk_pid_chika_v_v", &trk_pid_chika_v_v, &b_trk_pid_chika_v_v);
   fChain->SetBranchAddress("trk_pid_chimu_v_v", &trk_pid_chimu_v_v, &b_trk_pid_chimu_v_v);
   fChain->SetBranchAddress("trk_pfp_id_v", &trk_pfp_id_v, &b_trk_pfp_id_v);
   fChain->SetBranchAddress("trk_dir_x_v", &trk_dir_x_v, &b_trk_dir_x_v);
   fChain->SetBranchAddress("trk_dir_y_v", &trk_dir_y_v, &b_trk_dir_y_v);
   fChain->SetBranchAddress("trk_dir_z_v", &trk_dir_z_v, &b_trk_dir_z_v);
   fChain->SetBranchAddress("trk_start_x_v", &trk_start_x_v, &b_trk_start_x_v);
   fChain->SetBranchAddress("trk_start_y_v", &trk_start_y_v, &b_trk_start_y_v);
   fChain->SetBranchAddress("trk_start_z_v", &trk_start_z_v, &b_trk_start_z_v);
   fChain->SetBranchAddress("trk_end_x_v", &trk_end_x_v, &b_trk_end_x_v);
   fChain->SetBranchAddress("trk_end_y_v", &trk_end_y_v, &b_trk_end_y_v);
   fChain->SetBranchAddress("trk_end_z_v", &trk_end_z_v, &b_trk_end_z_v);
   fChain->SetBranchAddress("trk_distance_v", &trk_distance_v, &b_trk_distance_v);
   fChain->SetBranchAddress("trk_theta_v", &trk_theta_v, &b_trk_theta_v);
   fChain->SetBranchAddress("trk_phi_v", &trk_phi_v, &b_trk_phi_v);
   fChain->SetBranchAddress("trk_len_v", &trk_len_v, &b_trk_len_v);
   fChain->SetBranchAddress("trk_mcs_muon_mom_v", &trk_mcs_muon_mom_v, &b_trk_mcs_muon_mom_v);
   fChain->SetBranchAddress("trk_energy_proton_v", &trk_energy_proton_v, &b_trk_energy_proton_v);
   fChain->SetBranchAddress("trk_energy_muon_v", &trk_energy_muon_v, &b_trk_energy_muon_v);
   fChain->SetBranchAddress("bdt_nuNCpi0", &bdt_nuNCpi0, &b_bdt_nuNCpi0);
   fChain->SetBranchAddress("bdt_numuCCpi0", &bdt_numuCCpi0, &b_bdt_numuCCpi0);
   fChain->SetBranchAddress("bdt_numuCC", &bdt_numuCC, &b_bdt_numuCC);
   fChain->SetBranchAddress("bdt_ext", &bdt_ext, &b_bdt_ext);
   fChain->SetBranchAddress("bdt_cosmic", &bdt_cosmic, &b_bdt_cosmic);
   fChain->SetBranchAddress("bdt_global", &bdt_global, &b_bdt_global);
   Notify();
}

Bool_t NeutrinoSelectionFilter::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NeutrinoSelectionFilter::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NeutrinoSelectionFilter::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NeutrinoSelectionFilter_cxx
