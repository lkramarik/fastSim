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

void getEfficiency(TString input = "D0.toyMc_c3_job0.root"){
    TFile *fileSim = new TFile(input, "READ");
    TNtuple* ntp = (TNtuple*)fileSim -> Get("nt");

    Float_t ptmin[10] = {0,1,2,3,4,5,6,7,8,9};
    Float_t ptmax[10] = {1,2,3,4,5,6,7,8,9,10};
    Float_t kdca[10] = {0,0,0,0,0,0,0,0,0,0};
    Float_t pdca[10] = {0,0,0,0,0,0,0,0,0,0};
    Float_t dcaD0[10] = {0,0,0,0,0,0,0,0,0,0};
    Float_t dcaDaughters[10] = {0,0,0,0,0,0,0,0,0,0};
    Float_t decayLength[10] = {0,0,0,0,0,0,0,0,0,0};
    Float_t cosTheta[10] = {0,0,0,0,0,0,0,0,0,0};
    Float_t eff[10], bincenter[10];
    TString preCut = "(rEta<1)";

    TH1D* hS = new TH1D("hS", "hS", 2000, 0.4, 2.4);
    TH1D* hAll = new TH1D("hAll", "hAll", 2000, 0.4, 2.4);

    TString cut, ptcut, cutSim;
    int nbins = 5;
    for (int i = 0; i < nbins; ++i) {
        bincenter[i] = (ptmax[i] + ptmin[i])/2;
        cout<<bincenter[i]<<endl;
        cutSim = Form("(kRDca>%f)", kdca[i]);
        ptcut = Form("(rPt>%f)&&(rPt<%f)", ptmin[i], ptmax[i]);
        cut = preCut+"&&"+cutSim+"&&"+ptcut;
        cout<<"cut is: "<<cut<<endl;

        ntp -> Project("hS", "rM", cut);
        cut = preCut+"&&"+ptcut;
        ntp -> Project("hAll", "rM", cut);
        hAll -> Draw();
        hS -> Draw("same");

        eff[i] = hS->GetEntries()/hAll->GetEntries();
        cout<<"Efficiency for "<<ptmin[i]<<" to "<<ptmax[i]<<" is: "<<eff[i]<<endl;
    }

    TGraph* gr = new TGraph(nbins,bincenter,eff);
    gr -> Draw("AC*");
    //    ntpS -> Project("hS", "D_mass", "(D_pt>=1)&&(D_pt<2)&&(k_dca>0.0087)&&(pi1_dca>0.0099)&&(dcaD0ToPv<0.0075)&&(dcaDaughters<0.0093)&&(D_decayL>0.0232)");


//    TFile *outFile = new TFile("./outputs/Run16QA_histos_adc_av.root", "RECREATE");
//    TList *histoList = (TList*) fileData->Get("picoD0AnaMaker");

}