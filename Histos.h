//begin Hist class
#include <TTree.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <vector>       // std::vector
#include <algorithm>    // std::find
#include <iostream>     // std::cout

using namespace std;

class Histos{
  private:
    string flabel, fXlabel;
    int fnbins;
    double fminbin;
    double fmaxbin;
  public:
    Histos (string label, string Xlabel, int nbins, double minbin, double maxbin){
        flabel = label; fXlabel = Xlabel; fnbins = nbins; fminbin = minbin; fmaxbin = maxbin;
    }
    ~Histos (){};

		string GetLabel()  {return flabel;}
    string GetXLabel() {return fXlabel;}

		std::vector< TH1F *> GetHists(){
      std::vector< TH1F *> h_histos;
    	string channel_labels[11] = {"1eother","1e0p0pi","1enp0pi","1muother","1mu1pi0","ncother","ncpi0","cosmic","outoffv","other","data"};
    	for(int i = 0; i < 11; i++){
    		TString hTitle = Form("h_%s_%s",channel_labels[i].c_str(), flabel.c_str());
    		TH1F *h_hist = new TH1F(hTitle,hTitle,fnbins,fminbin,fmaxbin);
        h_hist->GetXaxis()->SetTitle(fXlabel.c_str());
    		h_histos.push_back(h_hist);
    	}
    	return h_histos;
    };

};

class muHistos{
  private:
    //std::vector< TH1F * > hists;
    //TDirectory *dir;
    string flabel, fXlabel;
    int fnbins;
    double fminbin;
    double fmaxbin;
  public:
    muHistos (string label, string Xlabel, int nbins, double minbin, double maxbin){
        flabel = label; fXlabel = Xlabel; fnbins = nbins; fminbin = minbin; fmaxbin = maxbin;
    }
    ~muHistos (){};

		string GetLabel()  {return flabel;}
    string GetXLabel() {return fXlabel;}

		std::vector< TH1F * > GetHists(){
      std::vector<TH1F *> h_histos;
    	                                                    //00,10,N0,01,11,1N,0N,N1,NN,other
      std::vector<string> channel_labels = {"1mu0p0pi","1mu1p0pi","1muNp0pi","1mu0p1pi","1mu1p1pi","1mu1pNpi","1mu0pNpi","1muNp1pi","1muNpNpi","1muother"};
      //00,10,01,11,NN,other
      //std::vector<string> channel_labels = {"1mu0p0pi","1mu1p0pi","1mu0p1pi","1mu1p1pi","1muNpNpi","1muother"};
    	for(unsigned i = 0; i < channel_labels.size(); i++){
    		TString hTitle = Form("h_%s_%s",channel_labels[i].c_str(), flabel.c_str());
    		TH1F *h_hist = new TH1F(hTitle,hTitle,fnbins,fminbin,fmaxbin);
        h_hist->GetXaxis()->SetTitle(fXlabel.c_str());
    		h_histos.push_back(h_hist);
    	}
    	return h_histos;
    };

};

class CutHistos{
  private:
    string flabel, fXlabel;
    int fnbins, fnumCuts;
    double fminbin;
    double fmaxbin;
  public:
    CutHistos (string label, string Xlabel, int numCuts, int nbins, double minbin, double maxbin){
        flabel = label; fXlabel = Xlabel; fnbins = nbins; fminbin = minbin; fmaxbin = maxbin, fnumCuts = numCuts;
    }
    ~CutHistos (){};

		string GetLabel()  {return flabel;}
    string GetXLabel() {return fXlabel;}

		std::vector< std::vector< TH1F *> > GetHists(){
      std::vector< std::vector< TH1F *> > h_histos;
    	string channel_labels[11] = {"1eother","1e0p0pi","1enp0pi","1muother","1mu1pi0","ncother","ncpi0","cosmic","outoffv","other","data"};
      for(int k = 0; k < fnumCuts; k++){
        std::vector< TH1F *> vhist;
        for(int i = 0; i < 11; i++){
      		TString hTitle = Form("h_%s_C%d_%s",channel_labels[i].c_str(),k, flabel.c_str());
      		TH1F *h_hist = new TH1F(hTitle,hTitle,fnbins,fminbin,fmaxbin);
          h_hist->GetXaxis()->SetTitle(fXlabel.c_str());
      		vhist.push_back(h_hist);
      	}
        h_histos.push_back(vhist);
      }
    	return h_histos;
    };
};

class CutHistosGen{
  private:
    string flabel, fXlabel;
    int fnbins, fnumCuts;
    double fminbin;
    double fmaxbin;
    std::vector<string> fchanLabels;
  public:
    CutHistosGen (string label, string Xlabel, int numCuts, int nbins, double minbin, double maxbin, std::vector<string> chanLabels){
        flabel = label; fXlabel = Xlabel; fnbins = nbins; fminbin = minbin; fmaxbin = maxbin, fnumCuts = numCuts; fchanLabels=chanLabels;
    }
    ~CutHistosGen (){};

		string GetLabel()  {return flabel;}
    string GetXLabel() {return fXlabel;}

		std::vector< std::vector< TH1F *> > GetHists(){
      std::vector< std::vector< TH1F *> > h_histos;
    	//string channel_labels[11] = {"1eother","1e0p0pi","1enp0pi","1muother","1mu1pi0","ncother","ncpi0","cosmic","outoffv","other","data"};
      for(int k = 0; k < fnumCuts; k++){
        std::vector< TH1F *> vhist;
        for(unsigned i = 0; i < fchanLabels.size(); i++){
      		TString hTitle = Form("h_%s_C%d_%s",fchanLabels[i].c_str(),k, flabel.c_str());
      		TH1F *h_hist = new TH1F(hTitle,hTitle,fnbins,fminbin,fmaxbin);
		h_hist->GetXaxis()->SetTitle(fXlabel.c_str());
      		vhist.push_back(h_hist);
      	}
        h_histos.push_back(vhist);
      }
    	return h_histos;
    };
};
