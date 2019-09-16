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


void Plotter_NeutrinoSelectionFilter(){
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

  TFile *f_nue = new TFile(Form("%s/nue.root",filepath.c_str()));
  TFile *f_numu = new TFile(Form("%s/numu.root",filepath.c_str()));
  TFile *f_databnb = new TFile(Form("%s/databnb.root",filepath.c_str()));
  TFile *f_dirt = new TFile(Form("%s/dirt.root",filepath.c_str()));
  TFile *f_ext = new TFile(Form("%s/ext.root",filepath.c_str()));

  std::vector<string> histoLabels = get_labels(f_nue);

  vector<string> channel_labels = {"1eother","1e0p0pi","1enp0pi","1muother","1mu1pi0","ncother","ncpi0","cosmic","outoffv","other","data"};
  //vector<string> legend_labels = {"1eother","1e0p0pi","1enp0pi","1muother","1mu1pi0","ncother","ncpi0","cosmic","outoffv","dataEXT"};
  std::vector<unsigned> vOrder = {1,2,0,3,4,5,6,7,8,9,10};
  vector<int> hist_Fill = {kGreen,kGreen+3, kGreen-3, kAzure+2, kAzure+3, kAzure+6, kAzure+7, kRed-4, kRed+3, kGray+1, kBlack};
  vector<int> hist_FillStyle = {kGreen,kGreen+3, kGreen-3, kAzure+2, kAzure+3, kAzure+6, kAzure+7, kRed-4, kRed+3, kGray+1, 3544};
  vector<string> sample_labels = {"nue","numu","databnb","dirt","ext"};

  std::vector<string> xAxisLabels = { "Shower Energy [GeV]",
                                      "Shower Energy [GeV]",
                                      "Shower Energy [GeV]",
                                      "Shower Energy [GeV]",
                                      "Shower Energy [GeV]",
                                      "Shower Energy [GeV]",
                                      "n_showers_contained",
                                      "n_tracks_contained",
                                      "slnhits",
                                      "shr_distance",
                                      "trk_score",
                                      "shr_score",
                                      "score",
                                      "True #{nu_e} Energy [GeV]"};
  std::vector<double> xmax = {3, 3, 3, 3, 3, 3, 10, 12, 3500, 10, 1, 1, 1, 3};
  std::vector<double> ymax = {200, 200, 200, 200, 200, 200, 17, 500, 140, 400, 400, 400, 400, 2000};
  std::vector<double> logymax = {40000, 40000, 40000, 40000, 40000, 40000, 500, 1000, 140, 400, 400, 400, 400, 2000};
  std::vector<double> ymin = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.0001, 0.01, 0.01, 0.01, 0.01, 0.01};



  for (int i = 0; i < histoLabels.size(); ++i){
    cout << histoLabels[i] << endl;
    std::vector<TH1F *> histos_nue;
    std::vector<TH1F *> histos_numu;
    std::vector<TH1F *> histos_databnb;
    std::vector<TH1F *> histos_dirt;
    std::vector<TH1F *> histos_ext;
    for (int j = 0; j < channel_labels.size(); ++j){
      TH1F *h = (TH1F*)f_nue->Get(Form("%s/h_%s_%s", histoLabels[i].c_str(), channel_labels[j].c_str(), histoLabels[i].c_str()));
      histos_nue.push_back(h);
      TH1F *h1 = (TH1F*)f_numu->Get(Form("%s/h_%s_%s",histoLabels[i].c_str(),channel_labels[j].c_str(),histoLabels[i].c_str()));
      histos_numu.push_back(h1);
      TH1F *h2 = (TH1F*)f_databnb->Get(Form("%s/h_%s_%s",histoLabels[i].c_str(),channel_labels[j].c_str(),histoLabels[i].c_str()));
      histos_databnb.push_back(h2);
      TH1F *h3 = (TH1F*)f_dirt->Get(Form("%s/h_%s_%s",histoLabels[i].c_str(),channel_labels[j].c_str(),histoLabels[i].c_str()));
      histos_dirt.push_back(h3);
      TH1F *h4 = (TH1F*)f_ext->Get(Form("%s/h_%s_%s",histoLabels[i].c_str(),channel_labels[j].c_str(),histoLabels[i].c_str()));
      histos_ext.push_back(h4);
    }

    for (int j = 0; j < histos_nue.size()-1; ++j){
      histos_nue[j]->Add(histos_numu[j]);
      histos_nue[j]->Add(histos_dirt[j]);
    }
    THStack *hs = new THStack("hs","");
    auto legend = new TLegend(0.898,0.6,0.77,0.926);
    for (int j = 0; j < histos_nue.size()-1; ++j){
      histos_nue[vOrder[j]]->SetFillColor(hist_Fill[vOrder[j]]);
      histos_nue[vOrder[j]]->SetLineColor(hist_FillStyle[vOrder[j]]);

      hs->Add(histos_nue[vOrder[j]]);
      legend->AddEntry(histos_nue[vOrder[j]],channel_labels[vOrder[j]].c_str(),"f");
      //if(channel_labels[vOrder[j]] == "1eother") cout << "1eother = " << histos_nue[vOrder[j]]->Integral() << ", name = " << histos_nue[vOrder[j]]->GetTitle() << endl;
      //if(channel_labels[vOrder[j]] == "1enp0pi") cout << "1eNp0pi = " << histos_nue[vOrder[j]]->Integral() << ", name = " << histos_nue[vOrder[j]]->GetTitle() << endl;
      //if(channel_labels[vOrder[j]] == "1e0p0pi") cout << "1e0p0pi = " << histos_nue[vOrder[j]]->Integral() << ", name = " << histos_nue[vOrder[j]]->GetTitle() << endl;
    }

    histos_ext[10]->SetFillColor(kBlack);
    histos_ext[10]->SetFillStyle(3544);
    hs->Add(histos_ext[10]);


    legend->AddEntry(histos_ext[10],"data_ext","f");

    //auto rp = new TRatioPlot(histos_databnb[10],hs);
    /*
    TCanvas *c = new TCanvas("c","c",1500,1000);
    c->cd();
    c->SetLogy();
    histos_databnb[10]->SetTitle("");
    histos_databnb[10]->GetXaxis()->SetRangeUser(0.0,xmax[i]);
    histos_databnb[10]->GetYaxis()->SetRangeUser(ymin[i],logymax[i]);
    //histos_databnb[10]->GetYaxis()->SetTitle("Events/0.06 GeV");
    histos_databnb[10]->GetXaxis()->SetTitle(Form("%s",xAxisLabels[i].c_str()));
    histos_databnb[10]->Draw("e1");
    hs->Draw("hist same");
    legend->Draw("same");
    c->SaveAs(Form("LOG_%s.png",histoLabels[i].c_str()));
    delete c;

    TCanvas *d = new TCanvas("d","d",1500,1000);
    d->cd();
    histos_databnb[10]->SetTitle("");
    //histos_databnb[10]->GetYaxis()->SetTitle("Events/0.06 GeV");
    histos_databnb[10]->GetYaxis()->SetRangeUser(ymin[i],ymax[i]);
    histos_databnb[10]->GetXaxis()->SetTitle(Form("%s",xAxisLabels[i].c_str()));
    histos_databnb[10]->Draw("e1");
    hs->Draw("hist same");
    legend->Draw("same");
    d->SaveAs(Form("LINE_%s.png",histoLabels[i].c_str()));
    delete d;


    TCanvas *f = new TCanvas("f","f",1500,1000);
    f->cd();
    histos_databnb[10]->SetTitle("");
    //histos_databnb[10]->GetYaxis()->SetTitle("Events/0.06 GeV");
    histos_databnb[10]->GetYaxis()->SetRangeUser(ymin[i],2);
    histos_databnb[10]->GetXaxis()->SetRangeUser(0,xmax[i]);
    histos_databnb[10]->GetXaxis()->SetTitle(Form("%s",xAxisLabels[i].c_str()));
    histos_databnb[10]->Draw("e1");
    hs->Draw("hist same");
    legend->Draw("same");
    f->SaveAs(Form("ZOOMED_LINE_%s.png",histoLabels[i].c_str()));
    delete f;
    */

    //hs->GetYaxis()->SetRangeUser(ymin[i],logymax[i]);
    TCanvas *g = new TCanvas("e","e",1500,1100);
    auto rpg = new TRatioPlot(hs, histos_databnb[10]);
    //auto rp = new TRatioPlot(histos_databnb[10],hs);

    rpg->SetSeparationMargin(0);
    g->SetLogy();
    rpg->GetXaxis()->SetRangeUser(0.0,xmax[i]);
    //histos_databnb[10]->GetYaxis()->SetRangeUser(ymin[i],10);
    //rp->GetUpperRefYaxis()->SetRangeUser(ymin[i],10);
    rpg->Draw();
    //hs->GetYaxis()->SetRangeUser(ymin[i],logymax[i]);
    hs->SetMinimum(ymin[i]);//,logymax[i]);
    hs->SetMaximum(logymax[i]);
    hs->GetXaxis()->SetTitle(Form("%s",xAxisLabels[i].c_str()));
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

    legend->Draw("same");
    g->Update();
    g->SaveAs(Form("RP_LOG_%s.png",histoLabels[i].c_str()));
    delete g;


    TCanvas *e = new TCanvas("e","e",1500,1100);
    auto rp = new TRatioPlot(hs, histos_databnb[10]);
    //auto rp = new TRatioPlot(histos_databnb[10],hs);

    rp->SetSeparationMargin(0);
    //e->SetLogy();
    rp->GetXaxis()->SetRangeUser(0.0,xmax[i]);
    //histos_databnb[10]->GetYaxis()->SetRangeUser(ymin[i],10);
    //rp->GetUpperRefYaxis()->SetRangeUser(ymin[i],10);
    rp->Draw();
    //hs->GetYaxis()->SetRangeUser(ymin[i],logymax[i]);
    hs->SetMinimum(ymin[i]);//,logymax[i]);
    hs->SetMaximum(ymax[i]);
    hs->GetXaxis()->SetTitle(Form("%s",xAxisLabels[i].c_str()));
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

    legend->Draw("same");
    e->Update();
    e->SaveAs(Form("RP_LINE_%s.png",histoLabels[i].c_str()));
    delete e;


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
