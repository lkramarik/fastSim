
#include "TH3.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TLegend.h"


void cTauTest(TString inputName="/home/lukas/work/tmva_d0/sim/D0.toyMc.test.root"){
    TFile* fIn = new TFile(inputName, "READ");
    TNtuple* ntp = (TNtuple*)fIn -> Get("nt");

    Float_t rcMass, rcpt, rcy, rcDecayL, mcDecayL, mcMass, mcPt, mcEta;
    ntp -> SetBranchAddress("rM",&rcMass);
    ntp -> SetBranchAddress("rPt",&rcpt);
    ntp -> SetBranchAddress("rY",&rcy);
    ntp -> SetBranchAddress("decayLength",&rcDecayL);
//    ntp -> SetBranchAddress("",&mcDecayL);
    ntp -> SetBranchAddress("m",&mcMass);
    ntp -> SetBranchAddress("pt",&mcPt);
    ntp -> SetBranchAddress("eta",&mcEta);


//    float mcCtau = mcDecayL * 1.e4 * tD0->mcMass / (mcPt * cosh(mcEta));
    TH1F* hRcCtau = new TH1F("hRcCtau","",300,0,3000);
    TH1F* hRcCtau1 = new TH1F("hRcCtau1","",300,0,3000);

    for (int i = 0; i < ntp->GetEntries()/5; ++i) {
        ntp -> GetEntry(i);
        float rcMt = sqrt(rcMass*rcMass + rcpt*rcpt);
        float rcPz = rcMt*sinh(rcy);
        float rcP = sqrt(rcpt*rcpt+rcPz*rcPz);
        float rcCtau = rcDecayL * mcMass / rcP * 1e4;
        float rcCtau1 = rcDecayL * rcMass / rcP * 1e4;

        hRcCtau->Fill(rcCtau);
        hRcCtau1->Fill(rcCtau1);

    }


    TF1* fitExp=new TF1("expCtau", "[0]*exp(-x/[1])", 50, 500);
    fitExp->SetParameter(1,124);
    hRcCtau->Fit("expCtau","R");
    cout<<fitExp->GetParameter(1)<<" "<<fitExp->GetParError(1)<<endl;
    hRcCtau->Draw();
    hRcCtau1->Draw("same");


//    fIn->Close();



}