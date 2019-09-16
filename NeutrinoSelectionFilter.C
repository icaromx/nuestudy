#define NeutrinoSelectionFilter_cxx
#include "NeutrinoSelectionFilter.h"
#include <TTree.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <vector>       // std::vector
#include <algorithm>    // std::find
#include <iostream>     // std::cout

bool Found(string label, int nCuts, std::vector<string> v_cuts);
void MakeHistos(TH1F *h_hist[11], string label, int nbins, double minbin, double maxbin);
void vMakeHistos(std::vector<TH1F *> h_hist, string label, int nbins, double minbin, double maxbin);
void FillHistos(TH1F *h_hist[11], double fillvar, int category, float potweight);
void SaveHistos(TFile *outROOT, TH1F *h_hist[11], TDirectory *dir, string label);
//void SaveHistos(TFile *outROOT, std::vector<TH1F> *h_hist, TDirectory *dir, string label);


std::vector<TH1F *> vHistos(string label, string htitle, int nbins, double minbin, double maxbin);

void NeutrinoSelectionFilter::Loop(string sample, float potweight){
//   In a ROOT session, you can do:
//      root> .L NeutrinoSelectionFilter.C
//      root> NeutrinoSelectionFilter t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

	string cutvar;
	ifstream myfile ("HistosToMake.list");
	std::vector<string> v_cuts;
	while ( getline (myfile,cutvar) ){
		v_cuts.push_back(cutvar);
	}
	int nCuts = v_cuts.size();
	if (fChain == 0) return;

	//////////////////////////////////////
	// DEFINE HISTOGRAMS
	//////////////////////////////////////


	TDirectory *dir_C_nslice_P_shr_energy;
	TH1F *h_C_nslice_P_shr_energy[11];
	string L_C_nslice_P_shr_energy = "C_nslice_P_shr_energy";
	MakeHistos(h_C_nslice_P_shr_energy,L_C_nslice_P_shr_energy, 100, 0, 6);
	//std::vector<TH1F *> vh_C_nslice_P_shr_energy = vHistos(L_C_nslice_P_shr_energy, "nslice", 50, 0, 3);

	TDirectory *dir_C_nslice_P_n_showers_contained;
	TH1F *h_C_nslice_P_n_showers_contained[11];
	string L_C_nslice_P_n_showers_contained = "C_nslice_P_n_showers_contained";
	MakeHistos(h_C_nslice_P_n_showers_contained,L_C_nslice_P_n_showers_contained, 20, 0, 20);

	TDirectory *dir_C_nslice_P_n_tracks_contained;
	TH1F *h_C_nslice_P_n_tracks_contained[11];
	string L_C_nslice_P_n_tracks_contained = "C_nslice_P_n_tracks_contained";
	MakeHistos(h_C_nslice_P_n_tracks_contained,L_C_nslice_P_n_tracks_contained, 20, 0, 20);

	TDirectory *dir_C_nslice_P_slnhits;
	TH1F *h_C_nslice_P_slnhits[11];
	string L_C_nslice_P_slnhits = "C_nslice_P_slnhits";
	MakeHistos(h_C_nslice_P_slnhits,L_C_nslice_P_slnhits, 10000, 0, 10000);

	TDirectory *dir_C_nslice_P_shr_distance;
	TH1F *h_C_nslice_P_shr_distance[11];
	string L_C_nslice_P_shr_distance = "C_nslice_P_shr_distance";
	MakeHistos(h_C_nslice_P_shr_distance,L_C_nslice_P_shr_distance, 500, 0, 50);

  TDirectory *dir_C_nslice_P_trk_score;
	TH1F *h_C_nslice_P_trk_score[11];
	string L_C_nslice_P_trk_score = "C_nslice_P_trk_score";
	MakeHistos(h_C_nslice_P_trk_score,L_C_nslice_P_trk_score, 100, 0, 1);

  TDirectory *dir_C_nslice_P_shr_score;
  TH1F *h_C_nslice_P_shr_score[11];
  string L_C_nslice_P_shr_score = "C_nslice_P_shr_score";
  MakeHistos(h_C_nslice_P_shr_score,L_C_nslice_P_shr_score, 100, 0, 1);

  TDirectory *dir_C_nslice_P_score;
  TH1F *h_C_nslice_P_score[11];
  string L_C_nslice_P_score = "C_nslice_P_score";
  MakeHistos(h_C_nslice_P_score,L_C_nslice_P_score, 100, 0, 1);

	///////////////////////////////////////////////////////////////////////////////////////
  TDirectory *dir_C_slpdg_P_shr_energy;
	TH1F *h_C_slpdg_P_shr_energy[11];
	string L_C_slpdg_P_shr_energy = "C_slpdg_P_shr_energy";
	MakeHistos(h_C_slpdg_P_shr_energy, L_C_slpdg_P_shr_energy, 100, 0, 6);

	TDirectory *dir_C_n_showers_contained_P_shr_energy;
	TH1F *h_C_n_showers_contained_P_shr_energy[11];
	string L_C_n_showers_contained_P_shr_energy = "C_n_showers_contained_P_shr_energy";
	MakeHistos(h_C_n_showers_contained_P_shr_energy, L_C_n_showers_contained_P_shr_energy, 100, 0, 6);

	TDirectory *dir_C_n_tracks_contained_P_shr_energy;
	TH1F *h_C_n_tracks_contained_P_shr_energy[11];
	string L_C_n_tracks_contained_P_shr_energy = "C_n_tracks_contained_P_shr_energy";
	MakeHistos(h_C_n_tracks_contained_P_shr_energy, L_C_n_tracks_contained_P_shr_energy, 100, 0, 6);

	TDirectory *dir_C_cal_shr_energy_P_shr_energy;
	TH1F *h_C_cal_shr_energy_P_shr_energy[11];
	string L_C_cal_shr_energy_P_shr_energy = "C_cal_shr_energy_P_shr_energy";
	MakeHistos(h_C_cal_shr_energy_P_shr_energy, L_C_cal_shr_energy_P_shr_energy, 100, 0, 6);

	TDirectory *dir_C_frac_shr_energy_P_shr_energy;
	TH1F *h_C_frac_shr_energy_P_shr_energy[11];
	string L_C_frac_shr_energy_P_shr_energy = "C_frac_shr_energy_P_shr_energy";
	MakeHistos(h_C_frac_shr_energy_P_shr_energy, L_C_frac_shr_energy_P_shr_energy, 100, 0, 6);

  //////////////////////////////////////////////////////////////////////
  TDirectory *dir_C_nslice_P_nu_e;
  TH1F *h_C_nslice_P_nu_e[11];
  string L_C_nslice_P_nu_e = "C_nslice_P_nu_e";
  MakeHistos(h_C_nslice_P_nu_e,L_C_nslice_P_nu_e, 100, 0, 6);

  TDirectory *dir_C_slpdg_P_nu_e;
	TH1F *h_C_slpdg_P_nu_e[11];
	string L_C_slpdg_P_nu_e = "C_slpdg_P_nu_e";
	MakeHistos(h_C_slpdg_P_nu_e, L_C_slpdg_P_nu_e, 100, 0, 6);

	TDirectory *dir_C_n_showers_contained_P_nu_e;
	TH1F *h_C_n_showers_contained_P_nu_e[11];
	string L_C_n_showers_contained_P_nu_e = "C_n_showers_contained_P_nu_e";
	MakeHistos(h_C_n_showers_contained_P_nu_e, L_C_n_showers_contained_P_nu_e, 100, 0, 6);

	TDirectory *dir_C_n_tracks_contained_P_nu_e;
	TH1F *h_C_n_tracks_contained_P_nu_e[11];
	string L_C_n_tracks_contained_P_nu_e = "C_n_tracks_contained_P_nu_e";
	MakeHistos(h_C_n_tracks_contained_P_nu_e, L_C_n_tracks_contained_P_nu_e, 100, 0, 6);

	TDirectory *dir_C_cal_nu_e_P_nu_e;
	TH1F *h_C_cal_nu_e_P_nu_e[11];
	string L_C_cal_nu_e_P_nu_e = "C_cal_nu_e_P_nu_e";
	MakeHistos(h_C_cal_nu_e_P_nu_e, L_C_cal_nu_e_P_nu_e, 100, 0, 6);

	TDirectory *dir_C_frac_nu_e_P_nu_e;
	TH1F *h_C_frac_nu_e_P_nu_e[11];
	string L_C_frac_nu_e_P_nu_e = "C_frac_nu_e_P_nu_e";
	MakeHistos(h_C_frac_nu_e_P_nu_e, L_C_frac_nu_e_P_nu_e, 100, 0, 6);



	//////////////////////////////////////////////////////////////////////

	std::cout << "Sample is: " << sample << std::endl;
	Long64_t nentries = fChain->GetEntriesFast();
	std::cout << "nentries = " << nentries << std::endl;
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry); nbytes += nb;
		if(sample == "numu" && nu_pdg != 14) continue;
		//if (jentry == 500) break;

		///////////////////////////////////////////////////
		//CUTS
		///////////////////////////////////////////////////
		if(!(nslice == 1 && Found("nslice", nCuts, v_cuts))) continue;
		FillHistos(h_C_nslice_P_shr_energy, shr_energy, category, potweight);
    FillHistos(h_C_nslice_P_nu_e, nu_e, category, potweight);
		FillHistos(h_C_nslice_P_n_showers_contained, n_showers_contained, category, potweight);
		FillHistos(h_C_nslice_P_n_tracks_contained, n_tracks_contained, category, potweight);
		FillHistos(h_C_nslice_P_slnhits, slnhits, category, potweight);
		FillHistos(h_C_nslice_P_shr_distance, shr_distance, category, potweight);
    FillHistos(h_C_nslice_P_trk_score, trk_score, category, potweight);
    FillHistos(h_C_nslice_P_shr_score, shr_score, category, potweight);
    FillHistos(h_C_nslice_P_score, shr_score, category, potweight); FillHistos(h_C_nslice_P_score, trk_score, category, potweight);



    if(!(slpdg == 12 && Found("slpdg", nCuts, v_cuts))) continue;
    FillHistos(h_C_slpdg_P_shr_energy, shr_energy, category, potweight);


		if(!(n_showers_contained > 0 && Found("n_showers_contained", nCuts, v_cuts))) continue;
		//if(n_showers_contained > 0 && Found("n_showers_contained", nCuts, v_cuts))
		FillHistos(h_C_n_showers_contained_P_shr_energy, shr_energy, category, potweight);

		if(!(n_tracks_contained == 0 && Found("n_tracks_contained", nCuts, v_cuts))) continue;
		//if(n_tracks_contained == 0 && Found("n_tracks_contained", nCuts, v_cuts))
		FillHistos(h_C_n_tracks_contained_P_shr_energy, shr_energy, category, potweight);

		if(!((shr_energy_tot+0.02)/0.8 > 0.06 && Found("shr_energy_tot", nCuts, v_cuts))) continue;
		FillHistos(h_C_cal_shr_energy_P_shr_energy, shr_energy, category, potweight);

		if(!(shr_energy/shr_energy_tot > 0.8 && Found("shr_energy_tot", nCuts, v_cuts))) continue;
		FillHistos(h_C_frac_shr_energy_P_shr_energy, shr_energy, category, potweight);

	}


	TFile *outROOT = new TFile(Form("%s.root",sample.c_str()),"RECREATE");
	cout<<"writing out"<<endl;
	TTree *tInfo = new TTree("info","");
	tInfo->Branch("Sample", "sample", &sample);
	tInfo->Fill();



	////////////////////////////////////////
	//Save Histograms Here
	////////////////////////////////////////
	SaveHistos(outROOT, h_C_nslice_P_shr_energy, 	dir_C_nslice_P_shr_energy,	L_C_nslice_P_shr_energy);
  SaveHistos(outROOT, h_C_slpdg_P_shr_energy, 	dir_C_slpdg_P_shr_energy,	L_C_slpdg_P_shr_energy);
	SaveHistos(outROOT, h_C_n_showers_contained_P_shr_energy, 			dir_C_n_showers_contained_P_shr_energy, 		L_C_n_showers_contained_P_shr_energy);
	SaveHistos(outROOT, h_C_n_tracks_contained_P_shr_energy, 			dir_C_n_tracks_contained_P_shr_energy, 		L_C_n_tracks_contained_P_shr_energy);
	SaveHistos(outROOT, h_C_cal_shr_energy_P_shr_energy, 			dir_C_cal_shr_energy_P_shr_energy, 		L_C_cal_shr_energy_P_shr_energy);
	SaveHistos(outROOT, h_C_frac_shr_energy_P_shr_energy, 			dir_C_frac_shr_energy_P_shr_energy, 		L_C_frac_shr_energy_P_shr_energy);

	SaveHistos(outROOT, h_C_nslice_P_n_showers_contained, 	dir_C_nslice_P_n_showers_contained,	L_C_nslice_P_n_showers_contained);
	SaveHistos(outROOT, h_C_nslice_P_n_tracks_contained, 	dir_C_nslice_P_n_tracks_contained,	L_C_nslice_P_n_tracks_contained);
	SaveHistos(outROOT, h_C_nslice_P_slnhits, 	dir_C_nslice_P_slnhits,	L_C_nslice_P_slnhits);
	SaveHistos(outROOT, h_C_nslice_P_shr_distance, 	dir_C_nslice_P_shr_distance,	L_C_nslice_P_shr_distance);
  SaveHistos(outROOT, h_C_nslice_P_trk_score, 	dir_C_nslice_P_trk_score,	L_C_nslice_P_trk_score);
  SaveHistos(outROOT, h_C_nslice_P_shr_score, 	dir_C_nslice_P_shr_score,	L_C_nslice_P_shr_score);
  SaveHistos(outROOT, h_C_nslice_P_score, 	dir_C_nslice_P_score,	L_C_nslice_P_score);
  SaveHistos(outROOT, h_C_nslice_P_nu_e, 	dir_C_nslice_P_nu_e,	L_C_nslice_P_nu_e);




	/////////////////////////////////////////////////////////////////////

	outROOT->Write();
  	cout<<"finished writing"<<endl;
}

bool Found(string label, int nCuts, std::vector<string> v_cuts){
	std::vector<string>::iterator it;
	it = find (v_cuts.begin(), v_cuts.end(), label);
	int here = it - v_cuts.begin();
	bool found = here < nCuts;
	return found;
}
void FillHistos(TH1F *h_hist[11], double fillvar, int category, float potweight){
	int categories[11] = {1,10,11,2,21,3,31,4,5,6,0};
	for (int j = 0; j < 11; ++j){
		if(category==categories[j]) h_hist[j]->Fill(fillvar,potweight);
	}
}
void MakeHistos(TH1F *h_hist[11], string label, int nbins, double minbin, double maxbin){
	string channel_labels[11] = {"1eother","1e0p0pi","1enp0pi","1muother","1mu1pi0","ncother","ncpi0","cosmic","outoffv","other","data"};
	for(int i = 0; i < 11; i++){
		TString hTitle = Form("h_%s_%s",channel_labels[i].c_str(), label.c_str());
		h_hist[i] = new TH1F(hTitle,hTitle,nbins,minbin,maxbin);
	}
}
/*
void vMakeHistos(std::vector<TH1F *> h_hist, string label, int nbins, double minbin, double maxbin){
	string channel_labels[11] = {"1eother","1e0p0pi","1enp0pi","1muother","1mu1pi0","ncother","ncpi0","cosmic","outoffv","other","data"};
	for(int i = 0; i < 11; i++){
		TString hTitle = Form("vh_%s_%s",channel_labels[i].c_str(), label.c_str());
		TH1F *h = new TH1F(hTitle,hTitle,nbins,minbin,maxbin);
		//h_hist[i] = new TH1F(hTitle,hTitle,nbins,minbin,maxbin);
		h_hist.push_back(h);
		delete h;
	}
}
*/

void SaveHistos(TFile *outROOT, TH1F *h_hist[11], TDirectory *dir, string label){
	dir = outROOT->mkdir(label.c_str());
	dir -> cd();
	for (int i = 0; i < 11; ++i){
		dir -> Add(h_hist[i]);
	}
}
/*
std::vector<TH1F *> vHistos(string label, string htitle, int nbins, double minbin, double maxbin){
	std::vector<TH1F *> h_histos;
	string channel_labels[11] = {"1eother","1e0p0pi","1enp0pi","1muother","1mu1pi0","ncother","ncpi0","cosmic","outoffv","other","data"};
	for(int i = 0; i < 11; i++){
		TString hTitle = Form("vh_%s_%s",channel_labels[i].c_str(), label.c_str());
		TH1F *h_hist = new TH1F(hTitle,hTitle,nbins,minbin,maxbin);
		h_histos.push_back(h_hist);
		//h_histos[i] = new TH1F(hTitle,hTitle,nbins,minbin,maxbin);
	}
	return h_histos;
}
*/
