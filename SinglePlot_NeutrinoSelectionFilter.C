#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>

#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TMath.h>
#include <THStack.h>
#include <TAttFill.h>
#include <TRatioPlot.h>

using namespace std;

std::vector<string> get_labels(TFile *f);
TString getDir( const std::string& subdir );


void SinglePlot_NeutrinoSelectionFilter(string dirname, string varname, int cutnum, double xmin, double xmax, double ymin, double ymax){
  //TString workdir = getDir( "/Users/ivan/Work/eLEE");
  TString workdir = getDir( "$APP/work/eLEE");
  TString histdir = getDir( Form("%s/results/%s",workdir.Data(),dirname.c_str()));
  TString plotdir = getDir(Form("%s/PLOTS/%s",workdir.Data(),dirname.c_str()));

  gROOT -> SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetHistLineWidth(4);

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.05);

  //title format
  gStyle->SetTitleFont(22, "XYZ");
  gStyle->SetTitleSize(0.06, "XYZ");

  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1);

  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(22, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");


  gROOT->ForceStyle();
  string filepath = "results";

  TFile *f_nue = new TFile(Form("%s/nue.root",histdir.Data()));
  TFile *f_numu = new TFile(Form("%s/numu.root",histdir.Data()));
  TFile *f_databnb = new TFile(Form("%s/databnb.root",histdir.Data()));
  TFile *f_dirt = new TFile(Form("%s/dirt.root",histdir.Data()));
  TFile *f_ext = new TFile(Form("%s/ext.root",histdir.Data()));

  std::vector<string> histoLabels = get_labels(f_nue);

  vector<string> channel_labels = {"1e0p0pi","1enp0pi","1eother","1muother","1mu1pi0","ncother","ncpi0","cosmic","outoffv","other","data"};
  vector<int> hist_Fill = {kGreen,kGreen+3, kGreen-3, kAzure+2, kAzure+3, kAzure+6, kAzure+7, kRed-4, kRed+3, kGray+1, kBlack};
  vector<int> hist_FillStyle = {kGreen,kGreen+3, kGreen-3, kAzure+2, kAzure+3, kAzure+6, kAzure+7, kRed-4, kRed+3, kGray+1, 3544};

  //std::vector<double> xmax, ymax, logymax, xmin;
  
 // cout << "xmax.size() = " << xmax.size() << ", ymax.size() = " << ymax.size() <<", histoLabels.size() = " << histoLabels.size() << endl;	
  
    //cout << histoLabels[i] << endl;
    std::vector<TH1F *> histos_nue;
    std::vector<TH1F *> histos_numu;
    std::vector<TH1F *> histos_databnb;
    std::vector<TH1F *> histos_dirt;
    std::vector<TH1F *> histos_ext;
    //C0_P_closestNuCosmicDist/h_1e0p0pi_C0_P_closestNuCosmicDist
    for (int j = 0; j < channel_labels.size(); ++j){
      TString hname = Form("C%d_P_%s/h_%s_C%d_P_%s", cutnum, varname.c_str(), channel_labels[j].c_str(), cutnum, varname.c_str());
      TH1F *h = (TH1F*)f_nue->Get(hname);
      histos_nue.push_back(h);
      TH1F *h1 = (TH1F*)f_numu->Get(hname);
      histos_numu.push_back(h1);
      TH1F *h2 = (TH1F*)f_databnb->Get(hname);
      histos_databnb.push_back(h2);
      TH1F *h3 = (TH1F*)f_dirt->Get(hname);
      histos_dirt.push_back(h3);
      TH1F *h4 = (TH1F*)f_ext->Get(hname);
      histos_ext.push_back(h4);
    }
    for (int j = 0; j < histos_nue.size()-1; ++j){
      histos_nue[j]->Add(histos_numu[j]);
      histos_nue[j]->Add(histos_dirt[j]);
    }
    THStack *hs = new THStack("hs","");
    auto legend = new TLegend(0.898,0.6,0.6,0.926);
    for (int j = 0; j < histos_nue.size()-1; ++j){
      histos_nue[j]->SetFillColor(hist_Fill[j]);
      histos_nue[j]->SetLineColor(hist_FillStyle[j]);
      //histos_nue[j]->Rebin(10);
      hs->Add(histos_nue[j]);
      legend->AddEntry(histos_nue[j],Form("%s: %f",channel_labels[j].c_str(),histos_nue[j]->Integral() ),"f");
      //cout << histos_nue[j]->Integral() << ", name = " << histos_nue[j]->GetTitle() << endl;
    }
    TString picname = Form("C%d_%s.png", cutnum, varname.c_str());
    histos_ext[10]->SetFillColor(kBlack);
    histos_ext[10]->SetFillStyle(3544);
    //histos_ext[10]->Rebin(10);
    hs->Add(histos_ext[10]);
    legend->AddEntry(histos_ext[10],"EXT","f");
    TCanvas *g = new TCanvas("g","g",1500,1100);
    //histos_databnb[10]->Rebin(10);
    auto rpg = new TRatioPlot(hs, histos_databnb[10]);
    //auto rp = new TRatioPlot(histos_databnb[10],hs);
    rpg->SetSeparationMargin(0);
    g->SetLogy();
    //rpg->GetXaxis()->SetRangeUser(0.0,xmax[i]);
    //histos_databnb[10]->GetYaxis()->SetRangeUser(ymin[i],10);
    //rp->GetUpperRefYaxis()->SetRangeUser(ymin[i],10);
    rpg->Draw();
    //hs->GetYaxis()->SetRangeUser(ymin[i],logymax[i]);
    hs->SetMinimum(ymin);//,logymax[i]);
    //hs->SetMaximum(ymax[i]);
    hs->GetXaxis()->SetTitle(Form("%s",histos_databnb[10]->GetXaxis()->GetTitle()));
    rpg->GetLowYaxis()->SetNdivisions(505);
    //rp->GetUpperPad()->Range(0.0, ymin[i], 0.001,logymax[i]);
    //rp->GetUpperRefYaxis()->SetRangeUser(0,logymax[i]);
    //rp->GetUpperRefYaxis()->SetRangeUser(ymin[i],10);

    rpg->GetLowerRefYaxis()->SetTitle("(MC + EXT)/BNB");
    rpg->GetLowerRefYaxis()->SetTitleSize(0.035);
    rpg->GetLowerRefYaxis()->SetTitleOffset(1.2);
    rpg->GetLowerRefXaxis()->SetTitleSize(0.035);
    rpg->GetLowerRefXaxis()->SetTitleOffset(1.2);
    rpg->GetLowerRefGraph()->SetMaximum(2);

    //legend->Draw("same");
    g->Update();
    g->SaveAs(Form("%s/RP_LOG_%s",plotdir.Data(), picname.Data()));
    delete g;

    TCanvas *e = new TCanvas("e","e",1500,1100);
    auto rp = new TRatioPlot(hs, histos_databnb[10]);
    //auto rp = new TRatioPlot(histos_databnb[10],hs);

    rp->SetSeparationMargin(0);
    //e->SetLogy();
    //rp->GetXaxis()->SetRangeUser(0.0,xmax[i]);
    //histos_databnb[10]->GetYaxis()->SetRangeUser(ymin[i],10);
    //rp->GetUpperRefYaxis()->SetRangeUser(ymin[i],10);
    rp->Draw();
    //hs->GetYaxis()->SetRangeUser(ymin[i],logymax[i]);

    //hs->SetMinimum(ymin[i]);//,logymax[i]);//
   // hs->SetMaximum(ymax[i]);////

    //hs->GetXaxis()->SetTitle(Form("%s",xAxisLabels[i].c_str()));
    hs->GetXaxis()->SetTitle(Form("%s",histos_databnb[10]->GetXaxis()->GetTitle()));
    rp->GetLowYaxis()->SetNdivisions(505);
    //rp->GetUpperPad()->Range(0.0, ymin[i], 0.001,logymax[i]);
    //rp->GetUpperRefYaxis()->SetRangeUser(0,logymax[i]);
    //rp->GetUpperRefYaxis()->SetRangeUser(ymin[i],10);

    rp->GetLowerRefYaxis()->SetTitle("(MC + EXT)/BNB");
    rp->GetLowerRefYaxis()->SetTitleSize(0.035);
    rp->GetLowerRefYaxis()->SetTitleOffset(1.2);
    rp->GetLowerRefXaxis()->SetTitleSize(0.035);
    rp->GetLowerRefXaxis()->SetTitleOffset(1.2);
    rp->GetLowerRefGraph()->SetMaximum(2);

    //legend->Draw("same");
    e->Update();
    e->SaveAs(Form("%s/RP_LINE_%s",plotdir.Data(), picname.Data()));
    delete e;
}

std::vector<string> get_labels(TFile *f) {
  std::vector<string> labels;
  TKey *key;
  Int_t total = 0;
  TIter next((TList *)f->GetListOfKeys());
  while ((key = (TKey *)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TDirectory")) {
      TDirectory *dir = (TDirectory *)key->ReadObj();
      labels.push_back(dir->GetName());
      total++;
    }
  }
  return labels;
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
