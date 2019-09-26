#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

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


void muPlotter_NeutrinoSelectionFilter(string dirname){
  TString workdir = getDir( "/Users/ivan/Work/eLEE");
  TString histdir = getDir( Form("%s/results/%s",workdir.Data(),dirname.c_str()));
  TString plotdir = getDir( Form("%s/PLOTS/%s/1muother",workdir.Data(),dirname.c_str()));

  gROOT -> SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(2);
  //gStyle->SetLineWidth(2);
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

//  TFile *f_nue = new TFile(Form("%s/nue.root",histdir.Data()));
  TFile *f_numu = new TFile(Form("%s/numu_1muother.root",histdir.Data()));
//  TFile *f_databnb = new TFile(Form("%s/databnb.root",histdir.Data()));
//  TFile *f_dirt = new TFile(Form("%s/dirt.root",histdir.Data()));
//  TFile *f_ext = new TFile(Form("%s/ext.root",histdir.Data()));

  std::vector<string> histoLabels = get_labels(f_numu);
  cout << histoLabels.size() << endl;
  std::vector<string> channel_labels = {"1mu0p0pi","1mu1p0pi","1muNp0pi","1mu0p1pi","1mu1p1pi","1mu1pNpi","1mu0pNpi","1muNp1pi","1muNpNpi","1muother"};
  vector<int> hist_Fill = {kGreen,kGreen+3, kGreen-3, kAzure+2, kAzure+3, kAzure+6, kAzure+7, kRed-4, kRed+3, kGray+1, kBlack};
  vector<int> hist_FillStyle = {kGreen,kGreen+3, kGreen-3, kAzure+2, kAzure+3, kAzure+6, kAzure+7, kRed-4, kRed+3, kGray+1, 3544};

  std::vector<double> xmax = {2, 2, 2, 2, 2, 2, 2, 2, 10, 12, 3500, 10, 1, 1, 1, 3};
  std::vector<double> ymax = {1000,1000, 400, 400, 400, 400, 400, 400, 17, 500, 140, 400, 400, 400, 400, 2000};
  std::vector<double> logymax = {40000,40000,40000, 40000, 40000, 40000, 40000, 40000, 500, 1000, 140, 400, 400, 400, 400, 2000};
  std::vector<double> ymin = {0.01,0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.0001, 0.01, 0.01, 0.01, 0.01, 0.01};



  for (int i = 0; i < histoLabels.size(); ++i){
    cout << histoLabels[i] << endl;
    std::vector<TH1F *> histos_numu;
    for (int j = 0; j < channel_labels.size(); ++j){
      TH1F *h1 = (TH1F*)f_numu->Get(Form("%s/h_%s_%s",histoLabels[i].c_str(),channel_labels[j].c_str(),histoLabels[i].c_str()));
      cout << Form("%s/h_%s_%s",histoLabels[i].c_str(),channel_labels[j].c_str(),histoLabels[i].c_str()) << endl;
      histos_numu.push_back(h1);
    }
    THStack *hs = new THStack("hs","");
    auto legend = new TLegend(0.898,0.6,0.6,0.926);
    double sumEntries = 0;
    cout << "histos_numu.size() = " << histos_numu.size() << endl;
    for (int j = 0; j < histos_numu.size(); ++j){
      histos_numu[j]->SetFillColor(hist_Fill[j]);
      histos_numu[j]->SetLineColor(hist_FillStyle[j]);
      hs->Add(histos_numu[j]);
      legend->AddEntry(histos_numu[j],Form("%s: %f",channel_labels[j].c_str(),histos_numu[j]->Integral()),"f");
      sumEntries += histos_numu[j]->Integral();
      cout << histos_numu[j]->Integral() << ", name = " << histos_numu[j]->GetTitle() << endl;
    }
    cout << sumEntries << endl;
    TCanvas *c = new TCanvas("c","c",1500,1000);
    c->cd();
    c->SetLogy();
    hs->SetTitle("");
    hs->Draw("hist");
    hs->GetXaxis()->SetTitle(histos_numu[0]->GetXaxis()->GetTitle());
    hs->GetXaxis()->SetLimits(0.0,2.0);
    hs->GetYaxis()->SetLimits(0.0,6000);
    legend->Draw("same");
    c->SaveAs(Form("%s/LOG_%s.png",plotdir.Data(), histoLabels[i].c_str()));
    delete c;

    TCanvas *d = new TCanvas("d","d",1500,1000);
    d->cd();
    hs->SetTitle("");
    hs->Draw("hist");
    hs->GetXaxis()->SetTitle(histos_numu[0]->GetXaxis()->GetTitle());
    hs->GetXaxis()->SetLimits(0.0,2.0);
    hs->GetYaxis()->SetLimits(0.0,6000);
    legend->Draw("same");
    d->SaveAs(Form("%s/muLINE_%s.png",plotdir.Data(), histoLabels[i].c_str()));
    delete d;
  }

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
