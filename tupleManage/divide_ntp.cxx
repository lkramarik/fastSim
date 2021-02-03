// Divide NTP to signalLike and backgroundLike

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TString.h"
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TNtuple.h"

using namespace std;

void divide_ntp(TString input="D0.toyMc.1605.root") {
    //void divide_ntp(TString input="D0.toyMc.root") {	
//    TString folder = "/media/lukas/376AD6A434B7392F/work/sim/";
//    input="D0.toyMC.0910.root";

    TString folder = "";
//    input="D0.toyMc.0303.95perc.root";

    TFile* data = new TFile(folder+input ,"r");
    TFile *fileOut = new TFile("ntp_full_"+input, "RECREATE");  // output root file

    TNtuple *ntpOut= new TNtuple("ntp_signal","D Meson Tree","D_mass:D_decayL:D_cosThetaStar:cosTheta:D_pt:D_ptSIM:pi1_pt:k_pt:pi1_dca:k_dca:pReta:kReta:dcaDaughters:dcaD0ToPv:hft:pid:etas:mcEtas:tpc:weight:rapidities:MCrapidities:refMult");
    //    TNtuple *ntpOut= new TNtuple("ntp_sideband","D Meson Tree","D_mass:D_decayL:D_theta:D_cosThetaStar:cosTheta:D_pt:pi1_pt:k_pt:pi1_dca:k_dca:dcaDaughters:dcaD0ToPv");

    TNtuple *ntp = (TNtuple*) data->Get("nt");
    Float_t w, flag, D_theta, D_mass, D_pt, D_decayL, k_pt, pi1_pt, pi1_dca, k_dca, k_nSigma, pi1_nSigma, pi1_TOFinvbeta, k_TOFinvbeta, dcaDaughters, pi1_eventId, k_eventId, dca_d0, dcaD0ToPv, cosTheta, D_cosThetaStar, D_ptSIM, y, rY, kY, kRY, pY, pRY;
    Float_t pPID, kPID, kHft, pHft, pTpc, kTpc, kREta, pREta, dReta, eta, kEta, pEta, rapidities, MCrapidities,refMult;

    ntp->SetBranchAddress("rM", &D_mass);
    ntp->SetBranchAddress("decayLength", &D_decayL);
    ntp->SetBranchAddress("angle12", &D_theta);
    ntp->SetBranchAddress("cosThetaStar", &D_cosThetaStar);
    ntp->SetBranchAddress("cosTheta", &cosTheta);
    ntp->SetBranchAddress("rPt", &D_pt);
    ntp->SetBranchAddress("pt", &D_ptSIM);
    ntp->SetBranchAddress("pRPt", &pi1_pt);

    ntp->SetBranchAddress("y", &y);
    ntp->SetBranchAddress("rY", &rY);
    ntp->SetBranchAddress("kY", &kY);
    ntp->SetBranchAddress("kRY", &kRY);
    ntp->SetBranchAddress("pY", &pY);
    ntp->SetBranchAddress("pRY", &pRY);

    ntp->SetBranchAddress("kREta", &kREta);
    ntp->SetBranchAddress("pREta", &pREta);
    ntp->SetBranchAddress("rEta", &dReta);

    ntp->SetBranchAddress("kEta", &kEta);
    ntp->SetBranchAddress("pEta", &pEta);
    ntp->SetBranchAddress("eta", &eta);

    ntp->SetBranchAddress("kRPt", &k_pt);
    ntp->SetBranchAddress("pRDca", &pi1_dca);
    ntp->SetBranchAddress("kRDca", &k_dca);
    ntp->SetBranchAddress("dca12", &dcaDaughters);
    ntp->SetBranchAddress("dcaD0ToPv", &dcaD0ToPv);

    ntp->SetBranchAddress("kPID", &kPID);
    ntp->SetBranchAddress("pPID", &pPID);
    ntp->SetBranchAddress("kHft", &kHft);
    ntp->SetBranchAddress("pHft", &pHft);
    ntp->SetBranchAddress("pTpc", &pTpc);
    ntp->SetBranchAddress("kTpc", &kTpc);
    ntp->SetBranchAddress("refMult", &refMult);

    ntp->SetBranchAddress("w", &w);

    //     ntp->SetBranchAddress("k_nSigma", &k_nSigma);
    //     ntp->SetBranchAddress("pi1_nSigma", &pi1_nSigma);
    //     ntp->SetBranchAddress("pi1_TOFinvbeta", &pi1_TOFinvbeta);
    //     ntp->SetBranchAddress("k_TOFinvbeta", &k_TOFinvbeta);

    TH1D* hpt = new TH1D("hMcPt", "hMcPt", 600, 0., 6);
    ntp -> Project("hMcPt", "pt","");

    const int nNtVars = ntpOut->GetNvar();
    float ntVar[30];

    float pid, hft, etas, tpc, mcEtas;
    cout<<"lets do this"<<endl;
    cout<<ntp->GetEntries()<<endl;
    for (long int i = 0; i < ntp->GetEntries(); i++) {
        ntp->GetEntry(i);

        tpc=0;
        if ((pTpc>0) && (kTpc>0)) tpc=1;

        pid=0;
        if ((kPID>0) && (pPID>0)) pid=1;

        hft=0;
        if ((kHft>0) && (pHft>0)) hft =1;

        etas=0;
        if ((abs(kREta)<1) && (abs(pREta)<1)) etas =1;

        mcEtas=0;
        if ((abs(kEta)<1) && (abs(pEta)<1)) mcEtas =1;

        rapidities=0;
        if ((abs(pRY)<1) && (abs(kRY)<1)) rapidities =1;

        MCrapidities=0;
        if ((abs(kY)<1) && (abs(pY)<1)) MCrapidities =1;


        int ii = 0;
        ntVar[ii++]=D_mass;
        ntVar[ii++]=D_decayL;
        ntVar[ii++]=D_cosThetaStar;
        ntVar[ii++]=cosTheta;
        ntVar[ii++]=D_pt;
        ntVar[ii++]=D_ptSIM;
        ntVar[ii++]=pi1_pt;
        ntVar[ii++]=k_pt;

        ntVar[ii++]=pi1_dca;
        ntVar[ii++]=k_dca;

        ntVar[ii++]=pREta;
        ntVar[ii++]=kREta;

        ntVar[ii++]=dcaDaughters;
        ntVar[ii++]=dcaD0ToPv;
        ntVar[ii++]=hft;
        ntVar[ii++]=pid;
        ntVar[ii++]=etas;
        ntVar[ii++]=mcEtas;
        ntVar[ii++]=tpc;
        ntVar[ii++]=w;
        ntVar[ii++]=rapidities;
        ntVar[ii++]=MCrapidities;
        ntVar[ii++]=refMult;

        ntpOut->Fill(ntVar);

    }
    fileOut->cd();
    ntpOut->Write(ntpOut->GetName(), TObject::kOverwrite);
    hpt->Write("hMcPt");
    fileOut->Close();

    data->Close();
    cout<<"Done."<<endl;
}


void project() {

    TFile *data = new TFile("ntpTMVA_D0.toyMc.Large.root", "OPEN");  // output root file
    TNtuple* ntpS = (TNtuple*)data -> Get("ntp_signal");
    TH1D* hS = new TH1D("signal", "signal", 2000, 0.4, 2.4);
    ntpS -> Project("signal", "D_mass");
    hS -> Draw();



}


