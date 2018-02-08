/* *********************************************************************
 *
 *  Analysis code to read Lc *.toyMc.root files.
 *
 *  Authors:
 *            **Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  ** Code Maintainer
 *
 * *********************************************************************
*/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PlotFile;
#endif

#ifndef __CINT__
#include "iostream"
#include <string>
#include <cmath>
#include <vector>
#include <string>

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"

#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#endif

#include "dataDrivenFastSimulator.h"
#include "LcNt.h"
#include "anaCuts.h"

using namespace std;

TGraphErrors* gHftRatioCorrection = NULL;
TF1*          gf1AuAu010Weight = NULL;

TF1*          fpionNsig_eff = NULL;
TF1*          fpionNsigTof_eff = NULL; //divided 0.013
TF1*          fkaonNsig_eff = NULL;
TF1*          fkaonNsigTof_eff    = NULL;//divided 0.013
TGraphErrors*          gprotonNsig_eff = NULL;
TGraphErrors*          gprotonNsigTof_eff    = NULL;//divided 0.013

TF1* ftofpi[9];
TF1* ftofk[9];  
TF1* ftofp[9];

// TF1 *myGaus = NULL;

int getLcPtIndex(float const pt)
{
   int bin = -1;
   for (int i = 0; i < anaCuts::nPtBins; i++)
   {
      if ((pt >= anaCuts::PtEdge[i]) && (pt < anaCuts::PtEdge[i + 1]))
         bin = i;
   }
   return bin;
}

int getLcCentIndex(float const cent)
{
   int bin = -1;

   for (int i = 0; i < anaCuts::physNCentralities; ++i)
   {
      // if (cent <= anaCuts::physCentralityEdges[i + 1])//this is for combined MB
      if (cent < anaCuts::physCentralityEdges[i + 1])//this is for fine centrality bins
      {
         bin = i;
         break;
      }
   }

   return bin;
}

bool isGoodTrack(float const ptCut, float const pt, float const eta)
{
   return pt > ptCut && fabs(eta) < anaCuts::eta;
}

bool isGoodPair(float const pt, float const y, float const cosTheta, float const kDca, float const piDca, float const pDca,
                float const dcaDaughters, float const decayLength, float const dcaV0ToPv)
{
   int tmpIndex = getLcPtIndex(pt);
   if (tmpIndex < 0) return false;

   return fabs(y) < anaCuts::rapidity &&
          cosTheta > anaCuts::cosTheta[tmpIndex] &&
          kDca > anaCuts::kDca[tmpIndex] &&
          piDca > anaCuts::piDca[tmpIndex] &&
          pDca > anaCuts::pDca[tmpIndex] &&
          dcaDaughters < anaCuts::dcaDaughters[tmpIndex] &&
          decayLength > anaCuts::decayLength[tmpIndex] &&
          dcaV0ToPv < anaCuts::dcaV0ToPv[tmpIndex];
}

int trkHalf(float const phi)
{
   if (phi > anaCuts::rightHalfLowEdge && phi < anaCuts::rightHalfHighEdge) return +1; //right side
   else return -1;//lest side
}

bool gausRandom(float const mean, float const sigma)
{
   return gRandom->Gaus(mean, sigma);
}


struct Hists
{
   int centrality;
   float minPtCut;
   TH1D* hNoCuts;
   TH1D* hNoCutsPhysBinning;
   TH1D* hTopoCuts;
   TH1D* hHftMatchingOnly;
   TH1D* hPIDOnly;
   TH1D* hTpcOnly;
   TH1D* hTpcHftTopo;
   TH2D* h2MassPt;

   Hists(float ptCut, int cent)
   {

      centrality = cent;
      minPtCut = ptCut;

      int nBins = 300;
      float minPt = 0.;
      float maxPt = 15.;
      hNoCuts = new TH1D(Form("hNoCuts_minPt%i_%i", (int)(minPtCut * 1e3), cent), Form("No Cuts %s minPt>%1.2f", anaCuts::physCentralityName[cent].Data(), minPtCut), nBins, minPt, maxPt);
      hNoCutsPhysBinning = new TH1D(Form("hNoCutsPhysBinning_minPt%i_%i", (int)(minPtCut * 1e3), cent), Form("No Cuts Physics Binning %s minPt>%1.2f", anaCuts::physCentralityName[cent].Data(), minPtCut), anaCuts::physNPtBins, anaCuts::physPtEdge);
      hTopoCuts = new TH1D(Form("hTopoCuts_minPt%i_%i", (int)(minPtCut * 1e3), cent), Form("Topo Cuts %s minPt>%1.2f", anaCuts::physCentralityName[cent].Data(), minPtCut), nBins, minPt, maxPt);
      hHftMatchingOnly = new TH1D(Form("hHftMatchingOnly_minPt%i_%i", (int)(minPtCut * 1e3), cent), Form("HFT Matching Only %s minPt>%1.2f", anaCuts::physCentralityName[cent].Data(), minPtCut), nBins, minPt, maxPt);
      hPIDOnly= new TH1D(Form("hPIDOnly_minPt%i_%i", (int)(minPtCut * 1e3), cent), Form("PID Only %s minPt>%1.2f", anaCuts::physCentralityName[cent].Data(), minPtCut), nBins, minPt, maxPt);
      hTpcOnly = new TH1D(Form("hTpcOnly_minPt%i_%i", (int)(minPtCut * 1e3), cent), Form("TPC Only %s minPt>%1.2f", anaCuts::physCentralityName[cent].Data(), minPtCut), nBins, minPt, maxPt);
      hTpcHftTopo = new TH1D(Form("hTpcHftTopo_minPt%i_%i", (int)(minPtCut * 1e3), cent), Form("TPC + HFT + Topo %s minPt>%1.2f", anaCuts::physCentralityName[cent].Data(), minPtCut), nBins, minPt, maxPt);
      h2MassPt = new TH2D(Form("h2MassPt_minPt%i_%i", (int)(minPtCut * 1e3), cent), Form("Invariant Mass vs. Pt %s minPt>%1.2f", anaCuts::physCentralityName[cent].Data(), minPtCut), nBins, minPt, maxPt, 210, 0, 2.1);

      hNoCuts->Sumw2();
      hNoCutsPhysBinning->Sumw2();
      hTopoCuts->Sumw2();
      hHftMatchingOnly->Sumw2();
      hPIDOnly->Sumw2();
      hTpcOnly->Sumw2();
      hTpcHftTopo->Sumw2();
      h2MassPt->Sumw2();
   }

   // void fill(LcNt const* const t, TF1* myGaus, TF1* ftofpi, TF1* ftofk)
   void fill(LcNt const* const t)
   {
      if (t->rPt > 10.) return;
      bool passTopologicalCuts = isGoodPair(t->rPt, t->rY, t->cosTheta,  t->kRDca, t->piRDca, t->pRDca, t->dcaDaughters, t->decayLength, t->dcaLcToPv);
      // bool passHft = t->kHft > 0 && t->piHft > 0 && t->pHft > 0;
      // bool passTpc = t->kTpc > 0 && t->piTpc > 0 && t->pTpc > 0;
      // float weight = t->pt * t->w;
      float weight = gf1AuAu010Weight->Eval(t->pt);
      // float weight = 1.;

      hNoCuts->Fill(t->pt, weight);
      hNoCutsPhysBinning->Fill(t->pt, weight);

      // if (!isGoodTrack(minPtCut, t->kRPt, t->kREta) || !isGoodTrack(minPtCut, t->piRPt, t->piREta) || !isGoodTrack(minPtCut, t->pRPt, t->pREta)) return;

      bool passTpc = t->kTpc > 0 && t->piTpc > 0 && t->pTpc > 0 && isGoodTrack(minPtCut, t->kRPt, t->kREta) && isGoodTrack(minPtCut, t->piRPt, t->piREta) && isGoodTrack(minPtCut, t->pRPt, t->pREta);

      float hftRatioWeight = matchHft(0, t->vz, t->cent, t->piRPt, t->piRPhi, t->piREta);
      hftRatioWeight      *= matchHft(1, t->vz, t->cent, t->kRPt, t->kRPhi, t->kREta);
      hftRatioWeight      *= matchHft(2, t->vz, t->cent, t->pRPt, t->pRPhi, t->pREta);

      float hftRatioCorrection = 1.0;

      float piontpcpid = fpionNsig_eff->Eval(t->piRPt);// typo Oct25
      // cout << " pi tpc PID = " << piontpcpid << endl;
      float kaontpcpid = fkaonNsig_eff->Eval(t->kRPt);
      // cout << " k tpc PID = " << kaontpcpid << endl;
      float protontpcpid = gprotonNsig_eff->Eval(t->pRPt);//typo Oct25
      // cout << " p tpc PID = " << protontpcpid << endl;


      float piontofpid = fpionNsigTof_eff->Eval(t->piRPt);
      // cout << " pi tof PID = " << piontofpid << endl;
      float kaontofpid = fkaonNsigTof_eff->Eval(t->kRPt);
      // cout << " k tof PID = " << kaontofpid << endl;
      float protontofpid = gprotonNsigTof_eff->Eval(t->pRPt);
      // cout << " p tof PID = " << protontofpid << endl;

      float piontofmatch = ftofpi[int(t->cent)]->Eval(t->piRPt);
      // float piontofmatch = ftofpi0->Eval(t->piRPt);
      // cout << " pi tof match = " << piontofmatch << endl;
      float kaontofmatch = ftofk[int(t->cent)]->Eval(t->kRPt);
      // float kaontofmatch = ftofk0->Eval(t->kRPt);
      // cout << " k tof match = " << kaontofmatch << endl;
      float protontofmatch = ftofp[int(t->cent)]->Eval(t->pRPt);
      // float protontofmatch = ftofp0->Eval(t->pRPt);
      // cout << " p tof match = " << protontofmatch << endl;

      float PID = (piontofmatch * piontofpid * piontpcpid + (1 - piontofmatch) * piontpcpid) * (kaontofmatch * kaontofpid * kaontpcpid) * (protontofmatch * protontofpid * protontpcpid);


      // cout << " PID = " << PID << endl;
      // if (gRandom->Rndm() < PID) return;
      // weight *= PID;



      /*if(t->pRPt>0.6 && t->pRPt<1.99)
      {
        hftRatioCorrection *= gHftRatioCorrection->Eval(t->pRPt);
      }

      if(t->kRPt>0.6 && t->kRPt<1.99)
      {
        hftRatioCorrection *= gHftRatioCorrection->Eval(t->kRPt);
      }
      */

      if (passTopologicalCuts) hTopoCuts->Fill(t->rPt, weight);
      if (passTpc) hTpcOnly->Fill(t->rPt, weight);

      hHftMatchingOnly->Fill(t->rPt, weight * hftRatioCorrection * hftRatioWeight);

      hPIDOnly->Fill(t->rPt, weight * PID);

      // weight *= PID;

      if (passTpc && passTopologicalCuts)
      {
         hTpcHftTopo->Fill(t->rPt, weight * PID * hftRatioCorrection * hftRatioWeight);
         h2MassPt->Fill(t->rPt, t->rM, weight * PID);
      }
   }

   void makeEffciency(TDirectory* fOut, TH1D* hPass, TH1D* hTotal, TH1D* hTotalPhysBinning = NULL, bool graphAsym = false)
   {
      fOut->cd();
      TString name = hPass->GetName();
      name.Replace(0, 1, "");
      TH1D* hEff = (TH1D*)hPass->Clone(Form("hEff%s", name.Data()));
      hEff->SetTitle(Form("%s Eff.", hPass->GetTitle()));
      hEff->Divide(hTotal);
      hEff->Write();

      if (hTotalPhysBinning)
      {
         TH1* hEffPhysBinning = hPass->Rebin(anaCuts::physNPtBins, Form("hEffPhysBinning%s", name.Data()), anaCuts::physPtEdge);
         hEffPhysBinning->SetTitle(Form("%s Eff. - Physics Binning", hPass->GetTitle()));
         hEffPhysBinning->Divide(hTotalPhysBinning);
         hEffPhysBinning->Write();
      }

      if (graphAsym)
      {
         TGraphAsymmErrors* grEff = new TGraphAsymmErrors(hPass, hTotal, "n");
         grEff->SetName(Form("grEff%s", name.Data()));
         grEff->SetTitle(Form("%s Eff.", hPass->GetTitle()));
         grEff->Write();

         if (hTotalPhysBinning)
         {
            TH1* hPassPhysBinning = hPass->Rebin(anaCuts::physNPtBins, Form("hPassPhysBinning%s", name.Data()), anaCuts::physPtEdge);
            TGraphAsymmErrors* grEffPhysBinning = new TGraphAsymmErrors(hPassPhysBinning, hTotalPhysBinning, "n");
            grEffPhysBinning->SetName(Form("grEffPhysBinning%s", name.Data()));
            grEffPhysBinning->SetTitle(Form("%s Eff. Physics Binning", hPass->GetTitle()));
            grEffPhysBinning->Write();
            delete hPassPhysBinning;
         }
      }
   }

   void write(TFile* fOut)
   {
      TDirectory* dir = NULL;
      if (!(dir = (TDirectory*)fOut->Get(Form("minPt%iMeV", (int)(minPtCut * 1e3)))))
         dir = (TDirectory*)fOut->mkdir(Form("minPt%iMeV", (int)(minPtCut * 1e3)));

      dir->cd();
      hNoCuts->Write();
      hTopoCuts->Write();
      hHftMatchingOnly->Write();
      hPIDOnly->Write();
      hTpcOnly->Write();
      hTpcHftTopo->Write();
      h2MassPt->Write();
      makeEffciency(dir, hTopoCuts, hNoCuts, hNoCutsPhysBinning, false);
      makeEffciency(dir, hHftMatchingOnly, hNoCuts, hNoCutsPhysBinning, false);
      makeEffciency(dir, hPIDOnly, hNoCuts, hNoCutsPhysBinning, false);
      makeEffciency(dir, hTpcOnly, hNoCuts, hNoCutsPhysBinning, false);
      makeEffciency(dir, hTpcHftTopo, hNoCuts, hNoCutsPhysBinning, false);
   }
};

struct TopoHists
{
   float minPtCut;
   TH3F* mcPointingAngle;
   TH3F* mcDecayL;
   TH3F* mcDcaDaughters;
   TH3F* mcProtonDca2Vtx;
   TH3F* mcPionDca2Vtx;
   TH3F* mcKaonDca2Vtx;
   TH3F* mcLcDca2Vtx;

   TopoHists(float const minPt)
   {
      float const maxDca = 0.3;

      int const nDcaBins = 300;
      float const maxDca12 = 0.1;
      int const nDca12Bins = 200;

      minPtCut = minPt;
      mcPointingAngle = new TH3F(Form("%s_se_us_pointingangle_minPt%i", "mc", (int)(minPtCut * 1.e3)), "Same Event US pointing angle; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 1000, 0.9, 1.0);
      mcDecayL = new TH3F(Form("%s_se_us_decayL_minPt%i", "mc", (int)(minPtCut * 1.e3)), "Same Event US Decay Length; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins, 0, maxDca);
      mcDcaDaughters = new TH3F(Form("%s_se_us_dcaDaughters_minPt%i", "mc", (int)(minPtCut * 1.e3)), "Same Event US dca daughters; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDca12Bins, 0, maxDca12);
      mcProtonDca2Vtx = new TH3F(Form("%s_se_us_protonDca_minPt%i", "mc", (int)(minPtCut * 1.e3)), "Same Event p dca 2 vertex; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins, 0, maxDca);
      mcPionDca2Vtx = new TH3F(Form("%s_se_us_pionDca_minPt%i", "mc", (int)(minPtCut * 1.e3)), "Same Event #pi dca 2 vertex; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins, 0, maxDca);
      mcKaonDca2Vtx = new TH3F(Form("%s_se_us_kaonDca_minPt%i", "mc", (int)(minPtCut * 1.e3)), "Same Event US K dca 2 vertex; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins, 0, maxDca);
      mcLcDca2Vtx = new TH3F(Form("%s_se_us_LcDca2Vtx_minPt%i", "mc", (int)(minPtCut * 1.e3)), "SameEvent US Lc dca 2 vertex; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 100, 0, 0.05);

      mcPointingAngle->Sumw2();
      mcDecayL->Sumw2();
      mcDcaDaughters->Sumw2();
      mcProtonDca2Vtx->Sumw2();
      mcPionDca2Vtx->Sumw2();
      mcKaonDca2Vtx->Sumw2();
      mcLcDca2Vtx->Sumw2();
   }

   void fill(LcNt* t)
   {
      if (!isGoodTrack(minPtCut, t->kRPt, t->kREta) || !isGoodTrack(minPtCut, t->piRPt, t->piREta) || !isGoodTrack(minPtCut, t->pRPt, t->pREta)) return;
      if (t->rM <  anaCuts::massMin || t->rM > anaCuts::massMax) return;
      bool passHft = t->kHft > 0 && t->piHft > 0 && t->pHft > 0;
      bool passTpc = t->kTpc > 0  && t->piTpc > 0 && t->pTpc > 0;
      if (!passHft || !passTpc) return;

      float weight = t->pt * t->w;
      int ptIndex = getLcPtIndex(t->rPt);
      if (ptIndex < 0) return;

      //Cos theta
      if (t->pRDca > anaCuts::pDca[ptIndex] &&
            t->piRDca > anaCuts::piDca[ptIndex] &&
            t->kRDca > anaCuts::kDca[ptIndex] &&
            t->dcaDaughters < anaCuts::dcaDaughters[ptIndex] &&
            t->decayLength > anaCuts::decayLength[ptIndex] &&
            // t->cosTheta > anaCuts::cosTheta[ptIndex] &&
            t->dcaLcToPv < anaCuts::dcaV0ToPv[ptIndex])
         mcPointingAngle->Fill(t->rPt, t->cent,  t->cosTheta, weight);

      //DecayL
      if (t->pRDca > anaCuts::pDca[ptIndex] &&
            t->piRDca > anaCuts::piDca[ptIndex] &&
            t->kRDca > anaCuts::kDca[ptIndex] &&
            t->dcaDaughters < anaCuts::dcaDaughters[ptIndex] &&
            // t->decayLength > anaCuts::decayLength[ptIndex] &&
            t->cosTheta > anaCuts::cosTheta[ptIndex] &&
            t->dcaLcToPv < anaCuts::dcaV0ToPv[ptIndex])
         mcDecayL->Fill(t->rPt, t->cent, t->decayLength / 1.e4, weight);

      //DcaDaughter
      if (t->pRDca > anaCuts::pDca[ptIndex] &&
            t->piRDca > anaCuts::piDca[ptIndex] &&
            t->kRDca > anaCuts::kDca[ptIndex] &&
            // t->dcaDaughters < anaCuts::dcaDaughters[ptIndex] &&
            t->decayLength > anaCuts::decayLength[ptIndex] &&
            t->cosTheta > anaCuts::cosTheta[ptIndex] &&
            t->dcaLcToPv < anaCuts::dcaV0ToPv[ptIndex])
         mcDcaDaughters->Fill(t->rPt, t->cent, t->dcaDaughters / 1.e4, weight);

      //ProtonDca
      if (//t->pRDca > anaCuts::pDca[ptIndex] &&
         t->piRDca > anaCuts::piDca[ptIndex] &&
         t->kRDca > anaCuts::kDca[ptIndex] &&
         t->dcaDaughters < anaCuts::dcaDaughters[ptIndex] &&
         t->decayLength > anaCuts::decayLength[ptIndex] &&
         t->cosTheta > anaCuts::cosTheta[ptIndex] &&
         t->dcaLcToPv < anaCuts::dcaV0ToPv[ptIndex])
         mcProtonDca2Vtx->Fill(t->rPt, t->cent, t->pRDca / 1.e4, weight);

      //PionDca
      if (t->pRDca > anaCuts::pDca[ptIndex] &&
            //t->piRDca > anaCuts::piDca[ptIndex] &&
            t->kRDca > anaCuts::kDca[ptIndex] &&
            t->dcaDaughters < anaCuts::dcaDaughters[ptIndex] &&
            t->decayLength > anaCuts::decayLength[ptIndex] &&
            t->cosTheta > anaCuts::cosTheta[ptIndex] &&
            t->dcaLcToPv < anaCuts::dcaV0ToPv[ptIndex])
         mcPionDca2Vtx->Fill(t->rPt, t->cent, t->piRDca / 1.e4, weight);

      //Kaon Dca
      if (t->pRDca > anaCuts::pDca[ptIndex] &&
            t->piRDca > anaCuts::piDca[ptIndex] &&
            // t->kRDca > anaCuts::kDca[ptIndex] &&
            t->dcaDaughters < anaCuts::dcaDaughters[ptIndex] &&
            t->decayLength > anaCuts::decayLength[ptIndex] &&
            t->cosTheta > anaCuts::cosTheta[ptIndex] &&
            t->dcaLcToPv < anaCuts::dcaV0ToPv[ptIndex])
         mcKaonDca2Vtx->Fill(t->rPt, t->cent,  t->kRDca / 1.e4, weight);

      //Lc dca
      if (t->pRDca > anaCuts::pDca[ptIndex] &&
            t->piRDca > anaCuts::piDca[ptIndex] &&
            t->kRDca > anaCuts::kDca[ptIndex] &&
            t->dcaDaughters < anaCuts::dcaDaughters[ptIndex] &&
            t->decayLength > anaCuts::decayLength[ptIndex] &&
            t->cosTheta > anaCuts::cosTheta[ptIndex] &&
            // t->dcaLcToPv < anaCuts::dcaV0ToPv[ptIndex])
            1)
         mcLcDca2Vtx->Fill(t->rPt, t->cent, t->dcaLcToPv / 1.e4, weight);

   }

   void write(TFile* fOut)
   {
      fOut->cd();
      mcPointingAngle->Write();
      mcDecayL->Write();
      mcDcaDaughters->Write();
      mcProtonDca2Vtx->Write();
      mcPionDca2Vtx->Write();
      mcKaonDca2Vtx->Write();
      mcLcDca2Vtx->Write();
   }
};

int main(int argc, char **argv)
{
   TStopwatch*   stopWatch = new TStopwatch();
   stopWatch->Start();

   loadHftRatio();
   gRandom->SetSeed();
   TFile* fHftRatioCorrection = new TFile("hftRatioCorrection_v1.root");
   gHftRatioCorrection = (TGraphErrors*)fHftRatioCorrection->Get("Graph");
   // gHftRatioCorrection->SetDirectory(0);
   // fHftRatioCorrection->Close();

   TFile* fAuAu010Weight = new TFile("AuAu010_weight.root");
   gf1AuAu010Weight = (TF1*)fAuAu010Weight->Get("f1Levy010");
   // gf1AuAu010Weight->SetDirectory(0);
   // fAuAu010Weight->Close();

   TFile* fpionPid = new TFile("pion_PidEff_FromXiaolong.root");
   fpionNsig_eff = (TF1*)fpionPid->Get("fNsig_eff")->Clone("fpionNsig_eff");
   // fpionNsig_eff->SetDirectory(0);
   fpionNsigTof_eff = (TF1*)fpionPid->Get("fNsigTof_eff")->Clone("fpionNsigTof_eff");
   // fpionNsigTof_eff->SetDirectory(0);
   // fpionPid->Close();

   TFile* fkaonPid = new TFile("kaon_PidEff_FromXiaolong.root");
   fkaonNsig_eff = (TF1*)fkaonPid->Get("fNsig_eff")->Clone("fkaonNsig_eff");
   // fkaonNsig_eff->SetDirectory(0);
   fkaonNsigTof_eff = (TF1*)fkaonPid->Get("fNsigTof_eff")->Clone("fkaonNsigTof_eff");
   // fkaonNsigTof_eff->SetDirectory(0);
   // fkaonPid->Close();

   TFile* fprotonPid = new TFile("proton_PidEff_Lambda0.root");
   gprotonNsig_eff = (TGraphErrors*)fprotonPid->Get("gNsig_eff")->Clone("gprotonNsig_eff");
   // gprotonNsig_eff->SetDirectory(0);
   gprotonNsigTof_eff = (TGraphErrors*)fprotonPid->Get("gNsigTof_eff")->Clone("gprotonNsigTof_eff");
   // gprotonNsigTof_eff->SetDirectory(0);
   // fprotonPid->Close();

   // TF1 *myGaus = new TF1("myGaus", "gaus(0)", -10, 10); //for tpc
   // myGaus->SetName("myGaus");
   // myGaus->SetNpx(100);

   TFile* fTofMatchFile = new TFile("tofMatch_fit_Run14_17Jan19.root");
   for (int i = 0; i < 9; i++)
   {
      ftofpi[i] = (TF1*)fTofMatchFile->Get(Form("funpip_%d", i))->Clone(Form("funpi_%d", i));
      // ftofpi[i]->SetDirectory(0);
      ftofk[i] = (TF1*)fTofMatchFile->Get(Form("funkm_%d", i))->Clone(Form("funk_%d", i));
      // ftofk[i]->SetDirectory(0);
      ftofp[i] = (TF1*)fTofMatchFile->Get(Form("funpp_%d", i))->Clone(Form("funp_%d", i));
      // ftofp[i]->SetDirectory(0);
   }
      // ftofpi0 = (TF1*)fTofMatchFile->Get(Form("funpi_0"))->Clone(Form("funpi_0"));
      // ftofk0 = (TF1*)fTofMatchFile->Get(Form("funk_0"))->Clone(Form("funk_0"));
      // ftofp0 = (TF1*)fTofMatchFile->Get(Form("funp_0"))->Clone(Form("funp_0"));
   // fTofMatchFile->Close();

   //
   std::string file = argv[1];
   LcNt* t = new LcNt(file);

   TFile* fOut = new TFile("eff.root", "recreate");

   Long64_t nEntries = t->GetEntries();
   cout << "nEntries = " << nEntries << endl;

   // std::vector<Hists> hists300MeV;
   std::vector<Hists> hists500MeV;
   // std::vector<Hists> hists1000MeV;

   TopoHists          topoHists(0.5);

   for (int iCent = 0; iCent < anaCuts::physNCentralities; ++iCent)
   {
      // hists300MeV.push_back(Hists(0.3,iCent));
      hists500MeV.push_back(Hists(0.5, iCent));
      // hists1000MeV.push_back(Hists(1.0,iCent));
   }

   // for (Long64_t i = 0; i < 1000000; ++i)
   for (Long64_t i = 0; i < t->GetEntries(); ++i)
   {
      t->GetEntry(i);

      if (i && i % 1000000 == 0) cout << static_cast<float>(i) / nEntries << endl;

      if (fabs(t->y) > anaCuts::rapidity) continue;

      // if (t->pid > 0) continue;//D+ 411
      // if (t->pid < 0) continue;//D- -411

      int dpmCentBin = getLcCentIndex(t->cent);
      if (dpmCentBin < 0) continue;
      // hists300MeV[dpmCentBin].fill(t);
      hists500MeV[dpmCentBin].fill(t);
      // hists1000MeV[dpmCentBin].fill(t);

      topoHists.fill(t);
   }

   for (int iCent = 0; iCent < anaCuts::physNCentralities; ++iCent)
   {
      // hists300MeV[iCent].write(fOut);
      hists500MeV[iCent].write(fOut);
      // hists1000MeV[iCent].write(fOut);
   }

   topoHists.write(fOut);
   fOut->Close();
   stopWatch->Stop();
   stopWatch->Print();
}
