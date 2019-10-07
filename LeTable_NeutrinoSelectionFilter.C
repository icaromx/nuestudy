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

using namespace std;

std::vector<string> get_labels(TFile *f);
TString getDir( const std::string& subdir );


void LeTable_NeutrinoSelectionFilter(string dirname, int kcuts, string histonames){

  TString workdir = getDir( "$APP/work/eLEE" );

  //TString workdir = getDir( "/Users/ivan/Work/eLEE");
  TString histdir = getDir( Form("%s/results/%s",workdir.Data(),dirname.c_str()));
  TString plotdir = getDir(Form("%s/PLOTS/%s",workdir.Data(),dirname.c_str()));

  TFile *f_nue =      new TFile(Form("%s/nue.root",histdir.Data()));
  TFile *f_numu =     new TFile(Form("%s/numu.root",histdir.Data()));
  TFile *f_databnb =  new TFile(Form("%s/databnb.root",histdir.Data()));
  TFile *f_dirt =     new TFile(Form("%s/dirt.root",histdir.Data()));
  TFile *f_ext =      new TFile(Form("%s/ext.root",histdir.Data()));

  std::vector<string> histoLabels = get_labels(f_nue);

//  vector<string> channel_labels = {"1eother","1e0p0pi","1enp0pi","1muother","1mu1pi0","ncother","ncpi0","cosmic","outoffv","other","data"};
  vector<string> channel_labels = {"1e0p0pi","1enp0pi","1eother","1muother","1mu1pi0","ncother","ncpi0","cosmic","outoffv","other","data"};
  //std::vector<unsigned> vOrder = {1,2,0,3,4,5,6,7,8,9,10};
  vector<string> sample_labels = {"nue","numu","databnb","dirt","ext"};
  std::vector<string> xAxisLabels = {"Shower Energy [GeV]","Shower Energy [GeV]", "Shower Energy [GeV]", "Shower Energy [GeV]","Shower Energy [GeV]","n_showers_contained","n_tracks_contained"};
  std::vector<double> xmax = {3, 3, 3, 3, 10, 12};
  std::vector<string> cutLabel = {"No Cut","nslice$==$1","slpdg$==$12","containedFraction$>$0.9","nShrContained$>$0","nTrkContained$==$0","(ShrEnergyTot+0.02)/0.8$>$0.06","ShrEnergy/ShrEnergyTot$>$0.8"};


  std::vector<std::vector<double> > numEntries;
  double pur, eff;
  cout << histoLabels.size() << endl;
  for (int i = 0; i < kcuts; ++i){
    std::vector<double> vEntries;
    //cout << histoLabels[i] << endl;
    std::vector<TH1F *> histos_nue;
    std::vector<TH1F *> histos_numu;
    std::vector<TH1F *> histos_databnb;
    std::vector<TH1F *> histos_dirt;
    std::vector<TH1F *> histos_ext;
    for (int j = 0; j < channel_labels.size(); ++j){
    //  cout << Form("%s/h_%s_%s", histoLabels[i].c_str(), channel_labels[j].c_str(), histoLabels[i].c_str()) << endl;
    /*
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
    */
	TH1F *h = (TH1F*)f_nue->Get(Form("C%d_P_%s/h_%s_C%d_P_%s", i, histonames.c_str(), channel_labels[j].c_str(), i, histonames.c_str()));
	cout << Form("C%d_P_%s/h_%s_C%d_P_%s", i, histonames.c_str(), channel_labels[j].c_str(), i, histonames.c_str()) << endl;
	cout << h->Integral() << endl;
	histos_nue.push_back(h);
	TH1F *h1 = (TH1F*)f_numu->Get(Form("C%d_P_%s/h_%s_C%d_P_%s", i, histonames.c_str(), channel_labels[j].c_str(), i, histonames.c_str()));
	histos_numu.push_back(h1);	
	TH1F *h2 = (TH1F*)f_databnb->Get(Form("C%d_P_%s/h_%s_C%d_P_%s", i, histonames.c_str(), channel_labels[j].c_str(), i, histonames.c_str()));
	histos_databnb.push_back(h2);
	TH1F *h3 = (TH1F*)f_dirt->Get(Form("C%d_P_%s/h_%s_C%d_P_%s", i, histonames.c_str(), channel_labels[j].c_str(), i, histonames.c_str()));
	histos_dirt.push_back(h3);
	TH1F *h4 = (TH1F*)f_ext->Get(Form("C%d_P_%s/h_%s_C%d_P_%s", i, histonames.c_str(), channel_labels[j].c_str(), i, histonames.c_str()));
	histos_ext.push_back(h4);
    }
    for (int j = 0; j < histos_nue.size()-1; ++j){
      histos_nue[j]->Add(histos_numu[j]);
      histos_nue[j]->Add(histos_dirt[j]);
    }
    for (int j = 0; j < histos_nue.size()-1; ++j){
      cout << histos_nue[j]->GetTitle() << " = " << histos_nue[j]->Integral() << endl;
      //vEntries.push_back(histos_nue[vOrder[j]]->Integral());
      vEntries.push_back(histos_nue[j]->Integral());
    }
    double total = 0;


    vEntries.push_back(histos_ext[10]->Integral());
    for (int i = 0; i < vEntries.size(); ++i) total += vEntries[i];
    vEntries.push_back(total);

    vEntries.push_back(histos_databnb[10]->Integral());

    if(i == 0) {
      eff = -1;
      pur = -1;
    }else{
      pur = vEntries[0]/total*100;
      eff = vEntries[0]/numEntries[0][0]*100;
    }
    vEntries.push_back(pur);
    vEntries.push_back(eff);
    numEntries.push_back(vEntries);
  }

  vector<string> tableL = {" ","1e0p0pi","1enp0pi","1eother","1muother","1mu1pi0","ncother","ncpi0","cosmic","outoffv","other","EXT","MC+EXT","BNB","Purity [\\%]","Efficiency[\\%]"};
  ofstream myfile;
  cout << Form("%s",plotdir.Data()) << endl; 
  //myfile.open (Form("%s/LeTable.txt",plotdir.Data()));
  myfile.open ("LeTable.txt");
  myfile << "\\begin{center} \n";
  myfile << "\\begin{tabular}{|l|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}\n";
  myfile << "\\hline\n";

  for (int i = 0; i < tableL.size(); ++i){
    if(i == tableL.size() -1) {
      myfile << tableL[i].c_str();
      myfile <<"\\\\ \n";
      myfile <<"\\hline \n";
    }else{
      myfile << tableL[i].c_str() << " & ";
    }
  }
  for (int i = 0; i < numEntries.size(); ++i){
    myfile << cutLabel[i].c_str() << " & ";
    for (int j = 0; j < numEntries[i].size(); ++j){

      if(j == numEntries[i].size()-1) {
        myfile << numEntries[i][j];
        myfile <<"\\\\ \n";
        myfile <<"\\hline \n";
      }else{
        myfile << numEntries[i][j] << " & ";
      }
    }
  }
  myfile << "\\end{tabular} \n";
  myfile << "\\end{center} \n";
  myfile.close();

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
