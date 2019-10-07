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


void Plotter_NeutrinoSelectionFilter(string dirname){
  //TString workdir = getDir( "/Users/ivan/Work/eLEE");
  TString workdir = getDir( "$APP/work/eLEE");
  TString histdir = getDir( Form("%s/results/%s",workdir.Data(),dirname.c_str()));
  TString plotdir = getDir(Form("%s/PLOTS/%s",workdir.Data(),dirname.c_str()));

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

  TFile *f_nue = new TFile(Form("%s/nue.root",histdir.Data()));
  TFile *f_numu = new TFile(Form("%s/numu.root",histdir.Data()));
  TFile *f_databnb = new TFile(Form("%s/databnb.root",histdir.Data()));
  TFile *f_dirt = new TFile(Form("%s/dirt.root",histdir.Data()));
  TFile *f_ext = new TFile(Form("%s/ext.root",histdir.Data()));

  std::vector<string> histoLabels = get_labels(f_nue);

  vector<string> channel_labels = {"1e0p0pi","1enp0pi","1eother","1muother","1mu1pi0","ncother","ncpi0","cosmic","outoffv","other","data"};
  //vector<string> legend_labels = {"1eother","1e0p0pi","1enp0pi","1muother","1mu1pi0","ncother","ncpi0","cosmic","outoffv","dataEXT"};
  //std::vector<unsigned> vOrder = {1,2,0,3,4,5,6,7,8,9,10};
  std::vector<unsigned> vOrderdavid = {3, 4, 5, 7, 8, 6, 0, 1, 2, 9, 10};
  vector<int> hist_Fill = {kGreen,kGreen+3, kGreen-3, kAzure+2, kAzure+3, kAzure+6, kAzure+7, kRed-4, kRed+3, kGray+1, kBlack};
  vector<int> hist_FillStyle = {kGreen,kGreen+3, kGreen-3, kAzure+2, kAzure+3, kAzure+6, kAzure+7, kRed-4, kRed+3, kGray+1, 3544};

  std::vector<double> xmax, ymax, logymax, xmin;
  double ymin = 0.01; 
  for(unsigned i = 0 ; i < 64; i++){
	if(i < 8){//EDepo 8
		xmax.push_back(2.0);
		xmin.push_back(0.0);
		ymax.push_back(400000);
	}else if(i>=8 && i < 16){//LowEDepo 8
		xmax.push_back(0.8);
		xmin.push_back(0.0);
		ymax.push_back(400000);
	}else if(i >= 16 && i < 40){//shr_dedx 3x8=24
		xmax.push_back(20);
		xmin.push_back(0.0);
		ymax.push_back(600);
	}else if(i >= 40 && i < 48){//shr_start 3x8 = 24
		xmax.push_back(1050.0);
		xmin.push_back(0.0);
		ymax.push_back(100);
	}else if(i >= 48 && i < 56){//shr_start 3x8 = 24
		xmax.push_back(120.0);
		xmin.push_back(-120.0);
		ymax.push_back(100);
	}else if(i >= 56 && i < 64){//shr_start 3x8 = 24
		xmax.push_back(250.0);
		xmin.push_back(0.0);
		ymax.push_back(100);
	}else if(i >= 64 && i < 72){//shr_start 3x8 = 24
		xmax.push_back(1050.0);
		xmin.push_back(0.0);
		ymax.push_back(100);
	}else if(i >= 72 && i < 80){//shr_start 3x8 = 24
		xmax.push_back(120.0);
		xmin.push_back(-120.0);
		ymax.push_back(100);
	}else if(i >= 80 && i < 88){//shr_start 3x8 = 24
		xmax.push_back(250.0);
		xmin.push_back(0.0);
		ymax.push_back(100);
	}else{//nu_sce_x 8
		xmax.push_back(250);
		xmin.push_back(0.0);
		ymax.push_back(20);
	}
  }
  cout << "xmax.size() = " << xmax.size() << ", ymax.size() = " << ymax.size() <<", histoLabels.size() = " << histoLabels.size() << endl;	

//  std::vector<double> xmax = {2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 3500, 10, 1, 1, 1, 3};
//  std::vector<double> ymax = {6000,1000, 400, 400, 400, 400, 400, 400, 17, 500, 140, 400, 400, 400, 400, 2000};
 // std::vector<double> logymax = {40000,40000, 40000, 40000, 40000, 40000, 40000, 40000, 500, 1000, 140, 400, 400, 400, 400, 2000};
//  std::vector<double> ymin = {0.01,0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.0001, 0.01, 0.01, 0.01, 0.01, 0.01};



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
    THStack *hsdavid = new THStack("hsdavid","");
    auto legend = new TLegend(0.898,0.6,0.6,0.926);
    auto legenddavid = new TLegend(0.898,0.6,0.6,0.926);
    for (int j = 0; j < histos_nue.size()-1; ++j){
      histos_nue[j]->SetFillColor(hist_Fill[j]);
      histos_nue[j]->SetLineColor(hist_FillStyle[j]);
      hs->Add(histos_nue[j]);
      hsdavid->Add(histos_nue[vOrderdavid[j]]);
      legend->AddEntry(histos_nue[j],Form("%s: %f",channel_labels[j].c_str(),histos_nue[j]->Integral() ),"f");
      legenddavid->AddEntry(histos_nue[vOrderdavid[j]],Form("%s: %f",channel_labels[vOrderdavid[j]].c_str(),histos_nue[vOrderdavid[j]]->Integral(0,34) ),"f");
      cout << histos_nue[j]->Integral() << ", name = " << histos_nue[j]->GetTitle() << endl;
    }

    cout << histos_ext[10]->Integral() << ", name = EXT" << histos_ext[10]->GetTitle() << endl;
    cout << histos_databnb[10]->Integral() << ", name = " << histos_databnb[10]->GetTitle() << endl;
    histos_ext[10]->SetFillColor(kBlack);
    histos_ext[10]->SetFillStyle(3544);
    hs->Add(histos_ext[10]);
    hsdavid->Add(histos_ext[10]);
    legend->AddEntry(histos_ext[10],"EXT","f");
    legenddavid->AddEntry(histos_ext[10],Form("EXT: %f",histos_ext[10]->Integral(0,34)),"f");
    legenddavid->AddEntry(histos_databnb[10],Form("BNB: %f",histos_databnb[10]->Integral(0,34)));
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
    TCanvas *g = new TCanvas("g","g",1500,1100);
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
    hs->SetMaximum(ymax[i]);
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

    legend->Draw("same");
    g->Update();
    g->SaveAs(Form("%s/RP_LOG_%s.png",plotdir.Data(), histoLabels[i].c_str()));
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

    legend->Draw("same");
    e->Update();
    e->SaveAs(Form("%s/RP_LINE_%s.png",plotdir.Data(), histoLabels[i].c_str()));
    delete e;

    continue;
    TCanvas *dav = new TCanvas("dav","dav",1500,1100);
    auto rpdav = new TRatioPlot(hsdavid, histos_databnb[10]);
    //auto rp = new TRatioPlot(histos_databnb[10],hs);

    rpdav->SetSeparationMargin(0);
    //e->SetLogy();
    rpdav->GetXaxis()->SetRangeUser(0.0,2.0);
    //histos_databnb[10]->GetYaxis()->SetRangeUser(ymin[i],10);
    //rp->GetUpperRefYaxis()->SetRangeUser(ymin[i],10);
    rpdav->Draw();
    //hs->GetYaxis()->SetRangeUser(ymin[i],logymax[i]);
    hsdavid->SetMinimum(ymin);//,logymax[i]);
    hsdavid->SetMaximum(500);
    //hsdavid->GetXaxis()->SetTitle(Form("%s",xAxisLabels[i].c_str()));
    hsdavid->GetXaxis()->SetTitle(Form("%s",histos_databnb[10]->GetXaxis()->GetTitle()));
    rpdav->GetLowYaxis()->SetNdivisions(505);
    //rp->GetUpperPad()->Range(0.0, ymin[i], 0.001,logymax[i]);
    //rp->GetUpperRefYaxis()->SetRangeUser(0,logymax[i]);
    //rp->GetUpperRefYaxis()->SetRangeUser(ymin[i],10);

    rpdav->GetLowerRefYaxis()->SetTitle("(MC + EXT)/BNB");
    rpdav->GetLowerRefYaxis()->SetTitleSize(0.035);
    rpdav->GetLowerRefYaxis()->SetTitleOffset(1.2);
    rpdav->GetLowerRefXaxis()->SetTitleSize(0.035);
    rpdav->GetLowerRefXaxis()->SetTitleOffset(1.2);
    rpdav->GetLowerRefGraph()->SetMaximum(2);

    legenddavid->Draw("same");
    dav->Update();
    dav->SaveAs(Form("%s/DAV_RP_LINE_%s.png",plotdir.Data(), histoLabels[i].c_str()));
    delete dav;
    
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
