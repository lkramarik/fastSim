#include<iostream>
#include"TFile.h"
#include"TProfile.h"
#include"TChain.h"
#include"TH1.h"
#include"TH2.h"
#include"TMath.h"
#include"TStyle.h"
#include"TPad.h"
#include"TCanvas.h"
#include"TFractionFitter.h"
#include"TObjArray.h"
#include"TLegend.h"
#include"TPaveText.h"
#include"TList.h"
#include"TAxis.h"
#include"TNtuple.h"
#include"TGraph.h"
#include <fstream>

using namespace std;
using namespace TMath;

void getEfficiency() {
    TString input = "D0.toyMc_c3_job0.root";

    Float_t eff[10], bincenter[10];

    const int nPtBins = 3;
    Double_t ptmin[] = {1, 2, 3, 4};
    Double_t ptmax[] = {2, 3, 5, 5};
    float decayL[] = {0.01, 0.01, 0.01, 0.01};
    float dcaDaughters[] = {0.013, 0.0046, 0.013, 0.013};
    float dcaD0[] = {0.0095, 0.0065, 0.0076, 0.0076};
    float cos[] = {0.5, 0.7, 0.7, 0.5};
    float kdca[] = {0.007, 0.0074, 0.007, 0.007};
    float pidca[] = {0.007, 0.006, 0.0079, 0.0079};

    for (int i = 0; i < nPtBins; ++i) {
        makeEff(input, ptmin[i], ptmax[i], decayL[i], dcaDaughters[i], dcaD0[i], cos[i], kdca[i], pidca[i]);
    }

}

void makeEff(TString input, float ptmin, float ptmax, float decayL, float dcaDaughters, float dcaD0, float cos, float kdca, float pidca){
    TFile *fileSim = new TFile(input, "READ");
    TNtuple *ntp = (TNtuple *) fileSim->Get("nt");

    TH1D* hS = new TH1D("hS", "hS", 1, 0, 4);
    TH1D* hAll = new TH1D("hAll", "hAll", 1, 0, 4);

    TString preCut = "abs(rEta<1)&&abs(kREta<1)&&abs(pREta<1)";

    TString ptCut = Form("rPt>%f && rPt<%f", ptmin, ptmax);
    TString decayLCut = Form("decayLength>%f",decayL);
    TString dcaDaughtersCut = Form("dca12<%f", dcaDaughters);
    TString dcaD0Cut = Form("dcaD0ToPv<%f", dcaD0);
    TString cosCut = Form("cosTheta>%f", cos);
    TString kdcaCut = Form("kRDca>%f",kdca);
    TString pidcaCut = Form("pRDca>%f", pidca);
    TString TOFmatching = "kHft>0 && pHft>0";
    TString setCuts = ptCut+" && "+preCut+" && "+decayLCut+" && "+dcaDaughtersCut+" && "+dcaD0Cut+" && "+cosCut+" && "+kdcaCut+" && "+pidcaCut+" && "+TOFmatching;
    cout<<"Cuts same as in data are:"<<endl;
    cout<<setCuts<<endl;

    Long64_t all=ntp->GetEntries(ptCut);
    Long64_t signal=ntp->GetEntries(setCuts);

    cout<<"Before reco: "<<all<<endl;
    cout<<"After reco: "<<signal<<endl;
    cout<<"Eff: "<<(double)signal/(double)all<<endl;

    ofstream yieldsErr;
    yieldsErr.open("raw_yields_eff.txt", std::ios::app);
    yieldsErr<<ptmin<<" "<<ptmax<<" "<<(double)signal/(double)all<<endl;
    yieldsErr.close();

//    bincenter[i] = (ptmax[i] + ptmin[i])/2;
//
//
//        TGraph* gr = new TGraph(nbins,bincenter,eff);
//        gr -> Draw("AC*");
//
//        cut = preCut+"&&"+ptcut;
//        ntp -> Project("hAll", "rM", cut);
//        hAll -> Draw();
//        hS -> Draw("same");
//
//        eff[i] = hS->GetEntries()/hAll->GetEntries();
//        cout<<"Efficiency for "<<ptmin[i]<<" to "<<ptmax[i]<<" is: "<<eff[i]<<endl;
    }

//    TGraph* gr = new TGraph(nbins,bincenter,eff);
//    gr -> Draw("AC*");
    //    ntpS -> Project("hS", "D_mass", "(D_pt>=1)&&(D_pt<2)&&(k_dca>0.0087)&&(pi1_dca>0.0099)&&(dcaD0ToPv<0.0075)&&(dcaDaughters<0.0093)&&(D_decayL>0.0232)");


//    TFile *outFile = new TFile("./outputs/Run16QA_histos_adc_av.root", "RECREATE");
//    TList *histoList = (TList*) fileData->Get("picoD0AnaMaker");

