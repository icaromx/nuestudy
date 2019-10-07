#define NeutrinoSelectionFilter_cxx
#include "NeutrinoSelectionFilter.h"
#include "Histos.h"
#include <TTree.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <vector>       // std::vector
#include <algorithm>    // std::find
#include <iostream>     // std::cout

using namespace std;

void vmuFillHistos(std::vector<TH1F *> vhists, double fillvar, int nproton, int npion, float potweight);
void vSaveCutHistos(TFile *outROOT, std::vector<std::vector<TH1F *> > vhists, string label);

void vFillHistos(std::vector<TH1F *> vhists, double fillvar, int category, float potweight);
void vSaveHistos(TFile *outROOT, std::vector<TH1F *> vhists, string label);

TString getDir( const std::string& subdir );
void NeutrinoSelectionFilter::Loop(string sample, float potweight_nosw, string namedir){
//	TString workdir = getDir( "/Users/ivan/Work/eLEE/" );
//	TString outdir = getDir( Form("%s/results/%s",workdir.Data(),namedir.c_str()) );
	TString workdir = getDir( "$APP/work/eLEE/" );
	TString outdir = getDir( Form("%s/results/%s",workdir.Data(),namedir.c_str()) );


	std::vector<string> muchanLabels = {"1mu0p0pi","1mu1p0pi","1muNp0pi","1mu0p1pi","1mu1p1pi","1mu1pNpi","1mu0pNpi","1muNp1pi","1muNpNpi","1muother"};
	std::vector<string> chanLabels   = {"1e0p0pi","1enp0pi","1eother","1muother","1mu1pi0","ncother","ncpi0","cosmic","outoffv","other","data"};

	string cutvar;
	ifstream myfile ("HistosToMake.list");
	std::vector<string> v_cuts;
	while ( getline (myfile,cutvar) ){
		v_cuts.push_back(cutvar);
	}
	int nCuts = v_cuts.size();
	if (fChain == 0) return;
	double Edepo, shrStart, potweight;
	int mucat;
	double Low_nuE = 0.2; double High_nuE = 0.6;
	//////////////////////////////////////
	// DEFINE HISTOGRAMS
	//////////////////////////////////////
	CutHistosGen C_Edepo("P_Edepo","Deposited Energy [GeV]",8,100,0,6, chanLabels);
	std::vector< std::vector< TH1F *> > hC_Edepo = C_Edepo.GetHists();

	CutHistosGen C_LowEdepo("P_LowEdepo","Deposited Energy [GeV]",8,100,0,6, chanLabels);
	std::vector< std::vector< TH1F *> > hC_LowEdepo = C_LowEdepo.GetHists();

	CutHistosGen muC_Edepo("muC_Edepo","Deposited Energy [GeV]",8,100,0,6, muchanLabels);
	std::vector< std::vector< TH1F *> > muhC_Edepo = muC_Edepo.GetHists();

	CutHistosGen C_shr_dedx_Y("P_shr_dedx_Y","Shower dE/dx (Y Plane) [GeV/cm]",8,100,0,20, chanLabels);
	std::vector< std::vector< TH1F *> > hC_shr_dedx_Y = C_shr_dedx_Y.GetHists();

	CutHistosGen C_shr_dedx_V("P_shr_dedx_V","Shower dE/dx (V Plane) [GeV/cm]",8,100,0,20, chanLabels);
	std::vector< std::vector< TH1F *> > hC_shr_dedx_V = C_shr_dedx_V.GetHists();

	CutHistosGen C_shr_dedx_U("P_shr_dedx_U","Shower dE/dx (U Plane) [GeV/cm]",8,100,0,20, chanLabels);
	std::vector< std::vector< TH1F *> > hC_shr_dedx_U = C_shr_dedx_U.GetHists();

	CutHistosGen C_shr_start_z("P_shr_start_z","Shower Start Z [cm]",8,500,0,1050, chanLabels);
	std::vector< std::vector< TH1F *> > hC_shr_start_z = C_shr_start_z.GetHists();
	
	CutHistosGen C_shr_start_y("P_shr_start_y","Shower Start Y [cm]",8,500,-120,120, chanLabels);
	std::vector< std::vector< TH1F *> > hC_shr_start_y = C_shr_start_y.GetHists();
	
	CutHistosGen C_shr_start_x("P_shr_start_x","Shower Start X [cm]",8,500,0,250, chanLabels);
	std::vector< std::vector< TH1F *> > hC_shr_start_x = C_shr_start_x.GetHists();
	
	CutHistosGen C_vtx_z("P_vtx_z","Vertex Z [cm]",8,500,0,1050, chanLabels);
	std::vector< std::vector< TH1F *> > hC_vtx_z = C_vtx_z.GetHists();
	
	CutHistosGen C_vtx_y("P_vtx_y","Vertex Y [cm]",8,500,-120,120, chanLabels);
	std::vector< std::vector< TH1F *> > hC_vtx_y = C_vtx_y.GetHists();
	
	CutHistosGen C_vtx_x("P_vtx_x","Vertex X [cm]",8,500,0,250, chanLabels);
	std::vector< std::vector< TH1F *> > hC_vtx_x = C_vtx_x.GetHists();

	CutHistosGen C_nu_sce_x("P_nu_sce_x","nu_sce_x [cm]",8,500,0,250, chanLabels);
	std::vector< std::vector< TH1F *> > hC_nu_sce_x = C_nu_sce_x.GetHists();
	//////////////////////////////////////////////////////

	std::cout << "Sample is: " << sample << std::endl;
	Long64_t nentries = fChain->GetEntriesFast();
	std::cout << "nentries = " << nentries << std::endl;
	Long64_t nbytes = 0, nb = 0;
	cout << "VECTOR SIZE: " << hC_Edepo.size() << endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++){

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry); nbytes += nb;
		///////VARIABLES///////
		Edepo = shr_energy_tot + 0.02 + trk_energy_tot;
		if(sample == "databnb" || sample == "ext") {
			potweight = potweight_nosw;
		}else{
			//potweight = potweight_nosw*weightSpline;
			potweight = potweight_nosw;
		}
		//////////////////////
		if(sample == "numu" && nu_pdg != 14) continue;
		///////////////////////////////////////////////////
		//CHAIN CUTS
		///////////////////////////////////////////////////
		if (category == 2) vmuFillHistos(muhC_Edepo[0], Edepo, nproton, npion, potweight);
		vFillHistos(hC_Edepo[0], Edepo, category, potweight);
		vFillHistos(hC_shr_dedx_Y[0], shr_dedx_Y, category, potweight);
		vFillHistos(hC_shr_dedx_V[0], shr_dedx_V, category, potweight);
		vFillHistos(hC_shr_dedx_U[0], shr_dedx_U, category, potweight);
		if (nu_e >= Low_nuE && nu_e <= High_nuE) vFillHistos(hC_LowEdepo[0], Edepo, category, potweight);
		vFillHistos(hC_shr_start_z[0], shr_start_z, category, potweight);
		vFillHistos(hC_shr_start_y[0], shr_start_y, category, potweight);
		vFillHistos(hC_shr_start_x[0], shr_start_x, category, potweight);
		vFillHistos(hC_vtx_z[0], vtx_z, category, potweight);
		vFillHistos(hC_vtx_y[0], vtx_y, category, potweight);
		vFillHistos(hC_vtx_x[0], vtx_x, category, potweight);
		vFillHistos(hC_nu_sce_x[0],nu_sce_x,category,potweight);
		///////////////////////////////////////////////////
		if(!(nslice == 1)) continue;

		vFillHistos(hC_Edepo[1], Edepo, category, potweight);
		vFillHistos(hC_shr_dedx_Y[1], shr_dedx_Y, category, potweight);
		vFillHistos(hC_shr_dedx_V[1], shr_dedx_V, category, potweight);
		vFillHistos(hC_shr_dedx_U[1], shr_dedx_U, category, potweight);
		if (category == 2) vmuFillHistos(muhC_Edepo[1], Edepo, nproton, npion, potweight);
		if (nu_e >= Low_nuE && nu_e <= High_nuE) vFillHistos(hC_LowEdepo[1], Edepo, category, potweight);
		vFillHistos(hC_shr_start_z[1], shr_start_z, category, potweight);
		vFillHistos(hC_shr_start_y[1], shr_start_y, category, potweight);
		vFillHistos(hC_shr_start_x[1], shr_start_x, category, potweight);
		vFillHistos(hC_vtx_z[1], vtx_z, category, potweight);
		vFillHistos(hC_vtx_y[1], vtx_y, category, potweight);
		vFillHistos(hC_vtx_x[1], vtx_x, category, potweight);
		vFillHistos(hC_nu_sce_x[1],nu_sce_x,category,potweight);

		///////////////////////////////////////////////////
		if(!(slpdg == 12)) continue;
		vFillHistos(hC_shr_dedx_Y[2], shr_dedx_Y, category, potweight);
		vFillHistos(hC_shr_dedx_V[2], shr_dedx_V, category, potweight);
		vFillHistos(hC_shr_dedx_U[2], shr_dedx_U, category, potweight);
		vFillHistos(hC_Edepo[2], Edepo, category, potweight);
		if (category == 2) vmuFillHistos(muhC_Edepo[2], Edepo, nproton, npion, potweight);
		if (nu_e >= Low_nuE && nu_e <= High_nuE) vFillHistos(hC_LowEdepo[2], Edepo, category, potweight);
		vFillHistos(hC_shr_start_z[2], shr_start_z, category, potweight);
		vFillHistos(hC_shr_start_y[2], shr_start_y, category, potweight);
		vFillHistos(hC_shr_start_x[2], shr_start_x, category, potweight);
		vFillHistos(hC_vtx_z[2], vtx_z, category, potweight);
		vFillHistos(hC_vtx_y[2], vtx_y, category, potweight);
		vFillHistos(hC_vtx_x[2], vtx_x, category, potweight);
		vFillHistos(hC_nu_sce_x[2],nu_sce_x,category,potweight);

		///////////////////////////////////////////////////
		if(!(contained_fraction > 0.9)) continue;
		vFillHistos(hC_shr_dedx_Y[3], shr_dedx_Y, category, potweight);
		vFillHistos(hC_shr_dedx_V[3], shr_dedx_V, category, potweight);
		vFillHistos(hC_shr_dedx_U[3], shr_dedx_U, category, potweight);
		if (category == 2) vmuFillHistos(muhC_Edepo[3], Edepo, nproton, npion, potweight);
		vFillHistos(hC_Edepo[3], Edepo, category, potweight);
		if (nu_e >= Low_nuE && nu_e <= High_nuE) vFillHistos(hC_LowEdepo[3], Edepo, category, potweight);
		vFillHistos(hC_shr_start_z[3], shr_start_z, category, potweight);
		vFillHistos(hC_shr_start_y[3], shr_start_y, category, potweight);
		vFillHistos(hC_shr_start_x[3], shr_start_x, category, potweight);
		vFillHistos(hC_vtx_z[3], vtx_z, category, potweight);
		vFillHistos(hC_vtx_y[3], vtx_y, category, potweight);
		vFillHistos(hC_vtx_x[3], vtx_x, category, potweight);
		vFillHistos(hC_nu_sce_x[3],nu_sce_x,category,potweight);

		///////////////////////////////////////////////////
		if(!(n_showers_contained > 0)) continue;
		vFillHistos(hC_shr_dedx_Y[4], shr_dedx_Y, category, potweight);
		vFillHistos(hC_shr_dedx_V[4], shr_dedx_V, category, potweight);
		vFillHistos(hC_shr_dedx_U[4], shr_dedx_U, category, potweight);
		if (category == 2) vmuFillHistos(muhC_Edepo[4], Edepo, nproton, npion, potweight);
		vFillHistos(hC_Edepo[4], Edepo, category, potweight);
		if (nu_e >= Low_nuE && nu_e <= High_nuE) vFillHistos(hC_LowEdepo[4], Edepo, category, potweight);
		vFillHistos(hC_shr_start_z[4], shr_start_z, category, potweight);
		vFillHistos(hC_shr_start_y[4], shr_start_y, category, potweight);
		vFillHistos(hC_shr_start_x[4], shr_start_x, category, potweight);
		vFillHistos(hC_vtx_z[4], vtx_z, category, potweight);
		vFillHistos(hC_vtx_y[4], vtx_y, category, potweight);
		vFillHistos(hC_vtx_x[4], vtx_x, category, potweight);
		vFillHistos(hC_nu_sce_x[4],nu_sce_x,category,potweight);

		///////////////////////////////////////////////////
		if(!(n_tracks_contained == 0)) continue;
		vFillHistos(hC_shr_dedx_Y[5], shr_dedx_Y, category, potweight);
		vFillHistos(hC_shr_dedx_V[5], shr_dedx_V, category, potweight);
		vFillHistos(hC_shr_dedx_U[5], shr_dedx_U, category, potweight);
		if (category == 2) vmuFillHistos(muhC_Edepo[5], Edepo, nproton, npion, potweight);
		vFillHistos(hC_Edepo[5], Edepo, category, potweight);
		if (nu_e >= Low_nuE && nu_e <= High_nuE) vFillHistos(hC_LowEdepo[5], Edepo, category, potweight);
		vFillHistos(hC_shr_start_z[5], shr_start_z, category, potweight);
		vFillHistos(hC_shr_start_y[5], shr_start_y, category, potweight);
		vFillHistos(hC_shr_start_x[5], shr_start_x, category, potweight);
		vFillHistos(hC_vtx_z[5], vtx_z, category, potweight);
		vFillHistos(hC_vtx_y[5], vtx_y, category, potweight);
		vFillHistos(hC_vtx_x[5], vtx_x, category, potweight);
		vFillHistos(hC_nu_sce_x[5],nu_sce_x,category,potweight);

		///////////////////////////////////////////////////
		if(!((shr_energy_tot+0.02)/0.8 > 0.06)) continue;
		vFillHistos(hC_shr_dedx_Y[6], shr_dedx_Y, category, potweight);
		vFillHistos(hC_shr_dedx_V[6], shr_dedx_V, category, potweight);
		vFillHistos(hC_shr_dedx_U[6], shr_dedx_U, category, potweight);
		if (category == 2) vmuFillHistos(muhC_Edepo[6], Edepo, nproton, npion, potweight);
		vFillHistos(hC_Edepo[6], Edepo, category, potweight);
		if (nu_e >= Low_nuE && nu_e <= High_nuE) vFillHistos(hC_LowEdepo[6], Edepo, category, potweight);
		vFillHistos(hC_shr_start_z[6], shr_start_z, category, potweight);
		vFillHistos(hC_shr_start_y[6], shr_start_y, category, potweight);
		vFillHistos(hC_shr_start_x[6], shr_start_x, category, potweight);
		vFillHistos(hC_vtx_z[6], vtx_z, category, potweight);
		vFillHistos(hC_vtx_y[6], vtx_y, category, potweight);
		vFillHistos(hC_vtx_x[6], vtx_x, category, potweight);
		vFillHistos(hC_nu_sce_x[6], nu_sce_x, category,potweight);

		///////////////////////////////////////////////////
		if(!(shr_energy/shr_energy_tot > 0.8)) continue;
		vFillHistos(hC_shr_dedx_Y[7], shr_dedx_Y, category, potweight);
		vFillHistos(hC_shr_dedx_V[7], shr_dedx_V, category, potweight);
		vFillHistos(hC_shr_dedx_U[7], shr_dedx_U, category, potweight);
		if (category == 2) vmuFillHistos(muhC_Edepo[7], Edepo, nproton, npion, potweight);
		vFillHistos(hC_Edepo[7], Edepo, category, potweight);
		if (nu_e >= Low_nuE && nu_e <= High_nuE) vFillHistos(hC_LowEdepo[7], Edepo, category, potweight);
		vFillHistos(hC_shr_start_z[7], shr_start_z, category, potweight);
		vFillHistos(hC_shr_start_y[7], shr_start_y, category, potweight);
		vFillHistos(hC_shr_start_x[7], shr_start_x, category, potweight);
		vFillHistos(hC_vtx_z[7], vtx_z, category, potweight);
		vFillHistos(hC_vtx_y[7], vtx_y, category, potweight);
		vFillHistos(hC_vtx_x[7], vtx_x, category, potweight);
		vFillHistos(hC_nu_sce_x[7],nu_sce_x,category,potweight);

	}

	TFile *outROOT = new TFile(Form("%s/%s.root",outdir.Data(),sample.c_str()),"RECREATE");
	//TFile *LeTable_outROOT = new TFile(Form("%s/LeTable_%s.root",outdir.Data(),sample.c_str()),"RECREATE");
	//TFile *muLeTable_outROOT = new TFile(Form("%s/LeTable_%s_1muother.root",outdir.Data(),sample.c_str()),"RECREATE");
	//TFile *outROOT = new TFile(Form("%s/%s.root",outdir.Data(),sample.c_str()),"RECREATE");
	TFile *muoutROOT = new TFile(Form("%s/%s_1muother.root",outdir.Data(),sample.c_str()),"RECREATE");
	cout<<"writing out"<<endl;
//	t1muother->Write();
	TTree *tInfo = new TTree("info","");
	tInfo->Branch("Sample", "sample", &sample);
	tInfo->Fill();

	////////////////////////////////////////
	//Save Histograms Here
	////////////////////////////////////////
	//cout << "my name is: " << HC_P_shr_energy.GetLabel() << endl;
	vSaveCutHistos(outROOT, hC_Edepo, C_Edepo.GetLabel());
	vSaveCutHistos(outROOT, hC_LowEdepo, C_LowEdepo.GetLabel());
	vSaveCutHistos(outROOT, hC_shr_dedx_Y, C_shr_dedx_Y.GetLabel());
	vSaveCutHistos(outROOT, hC_shr_dedx_V, C_shr_dedx_V.GetLabel());
	vSaveCutHistos(outROOT, hC_shr_dedx_U, C_shr_dedx_U.GetLabel());
	vSaveCutHistos(outROOT, hC_shr_start_z, C_shr_start_z.GetLabel());
	vSaveCutHistos(outROOT, hC_shr_start_y, C_shr_start_y.GetLabel());
	vSaveCutHistos(outROOT, hC_shr_start_x, C_shr_start_x.GetLabel());
	vSaveCutHistos(outROOT, hC_vtx_z, C_vtx_z.GetLabel());
	vSaveCutHistos(outROOT, hC_vtx_y, C_vtx_y.GetLabel());
	vSaveCutHistos(outROOT, hC_vtx_x, C_vtx_x.GetLabel());
	vSaveCutHistos(outROOT, hC_nu_sce_x, C_nu_sce_x.GetLabel());

	vSaveCutHistos(muoutROOT, muhC_Edepo, muC_Edepo.GetLabel());
	//vSaveCutHistos(muoutROOT, muhC_Edepo, muC_Edepo.GetLabel());
	/////////////////////////////////////////////////////////////////////
//	LeTable_outROOT->Write(); LeTable_outROOT->Close();
//	muLeTable_outROOT->Write(); muLeTable_outROOT->Close();
	outROOT->Write(); outROOT->Close();
	muoutROOT->Write(); muoutROOT->Close();
	cout<<"finished writing"<<endl;
}

void vFillHistos(std::vector<TH1F *> vhists, double fillvar, int category, float potweight){
	int categories[11] = {10,11,1,2,21,3,31,4,5,6,0};
	for (unsigned j = 0; j < vhists.size(); ++j){
		if(category==categories[j]) {
			vhists[j]->Fill(fillvar,potweight);
			break;
		}
	}
}

void vmuFillHistos(std::vector<TH1F *> vhists, double fillvar, int nproton, int npion, float potweight){
	int mucat; //00,10,N0,01,11,1N,0N,N1,NN,other
	if(nproton == 0 && npion == 0) 		{mucat = 0;}
	else if(nproton == 1 && npion == 0) {mucat = 1;} //10
	else if(nproton >  1 && npion == 0)	{mucat = 2;} //N0
	else if(nproton == 0 && npion == 1)	{mucat = 3;} //01
	else if(nproton == 1 && npion == 1) {mucat = 4;} //11
	else if(nproton == 1 && npion >  1)	{mucat = 5;} //1N
	else if(nproton == 0 && npion >  1)	{mucat = 6;} //0N
	else if(nproton >  1 && npion == 1)	{mucat = 7;} //N1
	else if(nproton >  1 && npion >  1)	{mucat = 8;} //NN
	else {mucat = 9;};															 //other
	vhists[mucat]->Fill(fillvar,potweight);
}

void vSaveHistos(TFile *outROOT, std::vector<TH1F *> vhists, string label){
	TDirectory *dir;
	dir = outROOT->mkdir(label.c_str());
	dir -> cd();
	for (unsigned i = 0; i < vhists.size(); ++i){
		dir -> Add(vhists[i]);
	}
}

void vSaveCutHistos(TFile *outROOT, std::vector<std::vector<TH1F *> > vhists, string label){
	for (unsigned k = 0; k < vhists.size(); k++) {
		TDirectory *dir;
		dir = outROOT->mkdir(Form("C%d_%s",k,label.c_str()));
		dir -> cd();
		for (unsigned i = 0; i < vhists[k].size(); ++i){
			dir -> Add(vhists[k][i]);
		}
	}
}

TString getDir( const std::string& subdir ){
	TString getdir( subdir );
 	if( 0 != system( Form( "test -d %s", getdir.Data()))){
		std::cout << "histos directory does not exist, making one now.... " << std::endl;
    int madedir = system( Form( "mkdir -m 755 -p %s", getdir.Data() ) );
    if( 0 != madedir ) std::cout << "HistoDir, Could not make plot directory, " << getdir << std::endl;
   }
 return getdir;
}
