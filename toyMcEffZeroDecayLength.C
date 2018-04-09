/* *********************************************************************
 *  ROOT macro - Toy Monte Carlo Simulation for D0 decay
 *  Includes Momentum Resolution, DCA, hft ration, TPC efficiency ...
 *  Example for D0 --> Kpi
 *
 *  Authors:
 *            Guannan Xie (guannanxie@lbl.gov)
 *            **Mustafa Mustafa (mmustafa@lbl.gov)
 *            Hao Qiu (hqiu@lbl.gov)
 *
 *  ** Code Maintainer
 *
 * *********************************************************************
*/

#include <iostream>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TF1.h"
#include "TClonesArray.h"
#include "TPythia6.h"
#include "TPythia6Decayer.h"
#include "TRandom3.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TMath.h"
#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "TStopwatch.h"
// #include "TSystem.h"
// #include "TMemStat.h"

using namespace std;

void setDecayChannels(int const mdme);
void decayAndFill(int const kf, TLorentzVector* b, double const weight, TClonesArray& daughters);
void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& piMom, TVector3 v00);
void getKinematics(TLorentzVector& b, double const mass);
TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution);
TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos);
TVector3 smearPosData(int iParticleIndex, double vz, int zdcb, TLorentzVector const& rMom, TVector3 const& pos);
float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0);
TVector3 getVertex(int centrality);
int getZdcBin(int const centrality);
bool matchHft(int iParticleIndex, double vz, int zdcb, TLorentzVector const& mom);
bool tpcReconstructed(int iParticleIndex, float charge, int cent, TLorentzVector const& mom);
bool reconstructD0(int const centrality, TLorentzVector const& mom);
bool matchTOF(int const iParticleIndex, TLorentzVector const& mom);
void bookObjects();
void write();
int getPtIndexDca(double);
int getEtaIndexDca(double);
int getVzIndexDca(double);
int getPhiIndexDca(double);

int getPtIndexHftRatio(double);
int getEtaIndexHftRatio(double);
int getVzIndexHftRatio(double);
int getPhiIndexHftRatio(double);

TPythia6Decayer* pydecay;
TNtuple* nt;
TFile* result;

TF1* fKaonMomResolution = NULL;
TF1* fPionMomResolution = NULL;
TF1* fWeightFunction = NULL;
TF1* fWeightFunctionAuAu = NULL;
TGraph* grEff[3];
const Int_t nParticles = 2;
const Int_t nCentHftRatio = 9;

const int nmultEdge = 7;
float const multEdge[nmultEdge+1] = {0, 4, 8, 12, 16, 20, 24, 200};


const Int_t nZdcX = 5;
const Double_t zdcxBins[] = {0,40,50,60,70,200};


const int m_nZdc = 10;
float const m_zdcEdge[m_nZdc+1] = {0,50,70,90,110,130,150,170,190,210,250};

// HFT ratio binning
const Int_t nEtasHftRatio = 10;
const Int_t nVzsHftRatio = 6;
const Int_t nPtBinsHftRatio = 15;
const Int_t nPhisHftRatio = 11;
const Double_t EtaEdgeHftRatio[nEtasHftRatio + 1] = //ok
{
   -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0
};
const Double_t VzEdgeHftRatio[nVzsHftRatio + 1] = //ok
{
   -6.0e4, -4.0e4, -2.0e4, 0.0, 2.0e4, 4.0e4, 6.0e4  
};
const Double_t ptEdgeHftRatio[nPtBinsHftRatio + 1] = { //ok
   0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,2.0,2.5,3.0,4.0,6.0,12.0
};
const Double_t PhiEdgeHftRatio[nPhisHftRatio + 1] = //ok
{
   -3.14159 , -2.80359 , -2.17527 , -1.54696 , -0.918637 , -0.290319 , 0.338 , 0.966319 , 1.59464 , 2.22296 , 2.85127 , 3.14159
};

int const nVzsDca = 4;
const Int_t nPhisDca = 11;
int const nEtasDca = 3;
const Int_t nPtBinsDca = 12;
float const VzEdgeDca[nVzsDca + 1] = {   -6.e4, -3.e4, 0, 3.e4, 6.e4};
float const EtaEdgeDca[nEtasDca + 1] = {0, 0.4, 0.8, 1};
const Double_t PhiEdgeDca[nPhisDca + 1] =
{
   -3.14159 , -2.80359 , -2.17527 , -1.54696 , -0.918637 , -0.290319 , 0.338 , 0.966319 , 1.59464 , 2.22296 , 2.85127 , 3.14159 //Sector by Sector  // sector number 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3
};
const Double_t ptEdgeDca[nPtBinsDca + 1] =
{
        0.2,0.3,0.4,0.5,0.6,0.8,1.,1.25,1.5,2.,3.0,5.,12.
};

int centralitySelect=0;

TH1D* h1Vz[nmultEdge+1];
TH1D* h1ZdcX[nmultEdge+1];

TH1D* hHftRatio1[nParticles][nEtasHftRatio][nVzsHftRatio][nPhisHftRatio][m_nZdc];
int const nCentDca = 9;
TH2D* h2Dca[nParticles][nEtasDca][nVzsDca][nZdcX][nPtBinsDca];

TH1D* hTpcPiPlus[nmultEdge]; //embedding
TH1D* hTpcPiMinus[nmultEdge]; //embedding
TH1D* hTpcKPlus[nmultEdge]; //embedding
TH1D* hTpcKMinus[nmultEdge]; //embedding
TH1D* h_pi_tof_eff;//embedding
TH1D* h_k_tof_eff;//embedding

string outFileName = "D0.toyMc";
std::pair<int, int> const decayChannels(747, 807);
std::pair<float, float> const momentumRange(0, 12);

float const gVzCut = 6.0e4;
float const acceptanceRapidity = 1.0;
float const sigmaPos0 = 15.2;
float const pxlLayer1Thickness = 0.00486;
float const sigmaVertexCent[nCentHftRatio] = {31., 18.1, 12.8, 9.3, 7.2, 5.9, 5., 4.6, 4.}; //not using

//============== main  program ==================
int jobindx;
void toyMcEffZeroDecayLength(int npart = 1e8, int jobId=0)
//void toyMcEffZeroDecayLength(int npart = 500)
{
   // TMemStat mem;
   // mem.Enable();
//   jobindx = jobId/9;
//   centralitySelect = jobId%9;

    jobindx = jobId/7;
    centralitySelect = jobId%7;

   TStopwatch*   stopWatch = new TStopwatch();
   stopWatch->Start();
   gRandom->SetSeed();
   bookObjects();

   pydecay = TPythia6Decayer::Instance();
   pydecay->Init();
    cout<<"let decay"<<endl;
   setDecayChannels(763); // D0 --> Kpi
   TLorentzVector* b_d = new TLorentzVector;
   TClonesArray ptl("TParticle", 10);
   cout<<"going to the for cycle of desiref npart"<<endl;
   for (int ipart = 0; ipart < npart; ipart++)
   {
      if (!(ipart % 1000))
         cout << "____________ ipart = " << ipart / static_cast<float>(npart) << " ________________" << endl;

      getKinematics(*b_d, M_D_0); //random pt, y, phi, return b vector with correct prop.
    cout<<"decay and fill"<<endl;
//      decayAndFill(421, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);
//      decayAndFill(-421, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);
      decayAndFill(421, b_d, fWeightFunctionAuAu->Eval(b_d->Perp()), ptl);  //421 = D0, ptl = daughters array
      decayAndFill(-421, b_d, fWeightFunctionAuAu->Eval(b_d->Perp()), ptl);
      cout<<"autosave"<<endl;
       if (ipart % 1000 == 1) nt->AutoSave("SaveSelf");
   }

   write();
   // mem.Show();
   stopWatch->Stop();
   stopWatch->Print();
}

void setDecayChannels(int const mdme)
{
   for (int idc = decayChannels.first; idc < decayChannels.second + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0); //747-807 switching off
   TPythia6::Instance()->SetMDME(mdme, 1, 1);
}

void decayAndFill(int const kf, TLorentzVector* b, double const weight, TClonesArray& daughters)
{
   pydecay->Decay(kf, b);
   pydecay->ImportParticles(&daughters);

   TLorentzVector kMom;
   TLorentzVector pMom;
   TVector3 v00;

   int nTrk = daughters.GetEntriesFast();
   cout<<"going into daughters"<<endl;
   for (int iTrk = 0; iTrk < nTrk; ++iTrk)
   {
      TParticle* ptl0 = (TParticle*)daughters.At(iTrk);

      switch (abs(ptl0->GetPdgCode()))
      {
         case 321: //kaonplus
            ptl0->Momentum(kMom); //seting momentum to kMom
            // v00.SetXYZ(0,0,0);
            v00.SetXYZ(ptl0->Vx() * 1000., ptl0->Vy() * 1000., ptl0->Vz() * 1000.); // converted to μm, production vertex
            break;
         case 211:   // pionplus
            ptl0->Momentum(pMom);
            break;
         default:
            break;
      }
   }
   daughters.Clear();
    cout<<"going to fill"<<endl;
   fill(kf, b, weight, kMom, pMom, v00);
}

void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& pMom, TVector3 v00)
{
   int const centrality = floor(nmultEdge * gRandom->Rndm());
    cout<<centrality<<endl;
   TVector3 const vertex = getVertex(centrality); //from hVz
    cout<<"after get vertex"<<endl;
   int zdcb = getZdcBin(centrality); //from data
    cout<<"vz and zdc taken"<<endl;
   // smear primary vertex
   // float const sigmaVertex = sigmaVertexCent[cent];
   // TVector3 const vertex(gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex));

   v00 += vertex; //vertex smear

   // smear momentum
   TLorentzVector const kRMom = smearMom(kMom, fKaonMomResolution); //fKaonMomResolution is TF1
   TLorentzVector const pRMom = smearMom(pMom, fPionMomResolution);

   // smear position
   TVector3 const kRPos = smearPosData(1, vertex.z(), zdcb, kRMom, v00);
   TVector3 const pRPos = smearPosData(0, vertex.z(), zdcb, pRMom, v00);
   // TVector3 const kRPos = smearPos(kMom, kRMom, v00);
   // TVector3 const pRPos = smearPos(pMom, pRMom, v00);

   // reconstruct
   TLorentzVector const rMom = kRMom + pRMom;
   float const kDca = dca(kMom.Vect(), v00, vertex);
   float const pDca = dca(pMom.Vect(), v00, vertex);
   float kRDca = dca(kRMom.Vect(), kRPos, vertex);
   float const kRSDca = dcaSigned(kRMom.Vect(), kRPos, vertex);
   float const kRDcaXY = dcaXY(kRMom.Vect(), kRPos, vertex);
   float const kRDcaZ = dcaZ(kRMom.Vect(), kRPos, vertex);
   float pRDca = dca(pRMom.Vect(), pRPos, vertex);
   float const pRSDca = dcaSigned(pRMom.Vect(), pRPos, vertex);
   float const pRDcaXY = dcaXY(pRMom.Vect(), pRPos, vertex);
   float const pRDcaZ = dcaZ(pRMom.Vect(), pRPos, vertex);

   TVector3 v0;
   float dca12 = dca1To2(kRMom.Vect(), kRPos, pRMom.Vect(), pRPos, v0);
   float decayLength = (v0 - vertex).Mag();
   float dcaD0ToPv = dca(rMom.Vect(), v0, vertex);
   float const cosTheta = (v0 - vertex).Unit().Dot(rMom.Vect().Unit());
   float const angle12 = kRMom.Vect().Angle(pRMom.Vect());

   TLorentzVector kRMomRest = kRMom;
   TVector3 beta;
   beta.SetMagThetaPhi(rMom.Beta(), rMom.Theta(), rMom.Phi());
   kRMomRest.Boost(-beta);
   float const cosThetaStar = rMom.Vect().Unit().Dot(kRMomRest.Vect().Unit());

   int const charge = kf > 0 ? 1 : -1;
   
   // save
   
   dca12 /= 1e4;
   decayLength /= 1e4;
   dcaD0ToPv /= 1e4;
   pRDca /= 1e4;
   kRDca /= 1e4;

   float arr[112];
   int iArr = 0;
   arr[iArr++] = centrality;
   arr[iArr++] = zdcb;
   arr[iArr++] = vertex.X();
   arr[iArr++] = vertex.Y();
   arr[iArr++] = vertex.Z();
   arr[iArr++] = getVzIndexDca(vertex.Z());

   arr[iArr++] = kf;
   arr[iArr++] = weight;
   arr[iArr++] = b->M();
   arr[iArr++] = b->Perp();
   arr[iArr++] = b->PseudoRapidity();
   arr[iArr++] = b->Rapidity();
   arr[iArr++] = b->Phi();
   arr[iArr++] = v00.X();
   arr[iArr++] = v00.Y();
   arr[iArr++] = v00.Z();

   arr[iArr++] = rMom.M();
   arr[iArr++] = rMom.Perp();
   arr[iArr++] = rMom.PseudoRapidity();
   arr[iArr++] = rMom.Rapidity();
   arr[iArr++] = rMom.Phi();
   arr[iArr++] = v0.X();
   arr[iArr++] = v0.Y();
   arr[iArr++] = v0.Z();

   arr[iArr++] = dca12;
   arr[iArr++] = decayLength;
   arr[iArr++] = dcaD0ToPv;
   arr[iArr++] = cosTheta;
   arr[iArr++] = angle12;
   arr[iArr++] = cosThetaStar;

   arr[iArr++] = kMom.M();
   arr[iArr++] = kMom.Perp();
   arr[iArr++] = kMom.PseudoRapidity();
   arr[iArr++] = kMom.Rapidity();
   arr[iArr++] = kMom.Phi();
   arr[iArr++] = kDca;

   arr[iArr++] = kRMom.M();
   arr[iArr++] = kRMom.Perp();
   arr[iArr++] = kRMom.PseudoRapidity();
   arr[iArr++] = kRMom.Rapidity();
   arr[iArr++] = kRMom.Phi();
   arr[iArr++] = kRPos.X();
   arr[iArr++] = kRPos.Y();
   arr[iArr++] = kRPos.Z();
   arr[iArr++] = kRDca;
   arr[iArr++] = kRSDca;
   arr[iArr++] = kRDcaXY;
   arr[iArr++] = kRDcaZ;
   arr[iArr++] = tpcReconstructed(1, -1 * charge, centrality, kRMom);

   arr[iArr++] = pMom.M();
   arr[iArr++] = pMom.Perp();
   arr[iArr++] = pMom.PseudoRapidity();
   arr[iArr++] = pMom.Rapidity();
   arr[iArr++] = pMom.Phi();
   arr[iArr++] = pDca;

   arr[iArr++] = pRMom.M();
   arr[iArr++] = pRMom.Perp();
   arr[iArr++] = pRMom.PseudoRapidity();
   arr[iArr++] = pRMom.Rapidity();
   arr[iArr++] = pRMom.Phi();
   arr[iArr++] = pRPos.X();
   arr[iArr++] = pRPos.Y();
   arr[iArr++] = pRPos.Z();
   arr[iArr++] = pRDca;
   arr[iArr++] = pRSDca;
   arr[iArr++] = pRDcaXY;
   arr[iArr++] = pRDcaZ;
   arr[iArr++] = tpcReconstructed(0, charge, centrality, pRMom);
   arr[iArr++] = matchTOF(1, kRMom);
   arr[iArr++] = matchTOF(0, pRMom);
   arr[iArr++] = matchHft(1, vertex.z(), zdcb,  kRMom); //kaon = 1, pion = 0
   arr[iArr++] = matchHft(0, vertex.z(), zdcb, pRMom);

   nt->Fill(arr);
}

void getKinematics(TLorentzVector& b, double const mass)
{
   float const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
   float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
   float const phi = TMath::TwoPi() * gRandom->Rndm();

   float const mT = sqrt(mass * mass + pt * pt);
   float const pz = mT * sinh(y);
   float const E = mT * cosh(y);

   b.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
}

float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
   TVector3 posDiff = pos - vertex;
   return fabs(p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff));
}

float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
   TVector3 posDiff = pos - vertex;
   float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;

   return sign * p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff);
}

float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
   TVector3 newPos(pos);
   newPos.SetZ(0);

   TVector3 newP(p);
   newP.SetZ(0);

   TVector3 newVertex(vertex);
   newVertex.SetZ(0);

   TVector3 posDiff = newPos - newVertex;
   float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;
   return sign * newP.Cross(posDiff.Cross(newP)).Unit().Dot(posDiff);
}

float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
   TVector3 posDiff = pos - vertex;
   if (sin(p.Theta()) == 0) return 0;
   else return (-posDiff.x() * cos(p.Phi()) - posDiff.y() * sin(p.Phi())) * cos(p.Theta()) / sin(p.Theta()) + posDiff.z();
}

float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0)
{
   TVector3 posDiff = pos2 - pos1;
   TVector3 pu1 = p1.Unit();
   TVector3 pu2 = p2.Unit();
   double pu1Pu2 = pu1.Dot(pu2);
   double g = posDiff.Dot(pu1);
   double k = posDiff.Dot(pu2);
   double s2 = (k - pu1Pu2 * g) / (pu1Pu2 * pu1Pu2 - 1.);
   double s1 = g + s2 * pu1Pu2;
   TVector3 posDca1 = pos1 + pu1 * s1;
   TVector3 posDca2 = pos2 + pu2 * s2;
   v0 = 0.5 * (posDca1 + posDca2);
   return (posDca1 - posDca2).Mag();
}

TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution)
{
   float const pt = b.Perp();
   float const sPt = gRandom->Gaus(pt, pt * fMomResolution->Eval(pt));

   TLorentzVector sMom;
   sMom.SetXYZM(sPt * cos(b.Phi()), sPt * sin(b.Phi()), sPt * sinh(b.PseudoRapidity()), b.M());
   return sMom;
}

TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos)
{
   float thetaMCS = 13.6 / mom.Beta() / rMom.P() / 1000 * sqrt(pxlLayer1Thickness / fabs(sin(mom.Theta())));
   float sigmaMCS = thetaMCS * 28000 / fabs(sin(mom.Theta()));
   float sigmaPos = sqrt(pow(sigmaMCS, 2) + pow(sigmaPos0, 2));

   return TVector3(gRandom->Gaus(pos.X(), sigmaPos), gRandom->Gaus(pos.Y(), sigmaPos), gRandom->Gaus(pos.Z(), sigmaPos));
}

int getPtIndexDca(double pT)
{
   for (int i = 0; i < nPtBinsDca; i++)
   {
      if ((pT >= ptEdgeDca[i]) && (pT < ptEdgeDca[i + 1]))
         return i;
   }
   if (pT < ptEdgeDca[0]) return 0;
   return nPtBinsDca - 1 ;
}

int getEtaIndexDca(double Eta)
{
   for (int i = 0; i < nEtasDca; i++)
   {
      if ((fabs(Eta) >= EtaEdgeDca[i]) && (fabs(Eta) < EtaEdgeDca[i + 1]))
         return i;
   }
   return nEtasDca - 1 ;
}

int getVzIndexDca(double Vz)
{
   for (int i = 0; i < nVzsDca; i++)
   {
      if ((Vz >= VzEdgeDca[i]) && (Vz < VzEdgeDca[i + 1]))
         return i;
   }
   return nVzsDca - 1 ;
}

int getPhiIndexDca(double Phi)
{
   for (int i = 0; i < nPhisDca; i++)
   {
      if ((Phi >= PhiEdgeDca[i]) && (Phi < PhiEdgeDca[i + 1]))
         return i;
   }
   return nPhisDca - 1 ;
}

int getPtIndexHftRatio(double pT)
{
   for (int i = 0; i < nPtBinsHftRatio; i++)
   {
      if ((pT >= ptEdgeHftRatio[i]) && (pT < ptEdgeHftRatio[i + 1]))
         return i;
   }
   if (pT<ptEdgeHftRatio[0]) return 0;
   return nPtBinsHftRatio - 1 ;
}

int getEtaIndexHftRatio(double Eta)
{
   for (int i = 0; i < nEtasHftRatio; i++)
   {
      if ((Eta >= EtaEdgeHftRatio[i]) && (Eta < EtaEdgeHftRatio[i + 1]))
         return i;
   }
   return nEtasHftRatio - 1 ;
}

int getVzIndexHftRatio(double Vz)
{
   for (int i = 0; i < nVzsHftRatio; i++)
   {
      if ((Vz >= VzEdgeHftRatio[i]) && (Vz < VzEdgeHftRatio[i + 1]))
         return i;
   }
   return nVzsHftRatio - 1 ;
}

int getPhiIndexHftRatio(double Phi)
{
   for (int i = 0; i < nPhisHftRatio; i++)
   {
      if ((Phi >= PhiEdgeHftRatio[i]) && (Phi < PhiEdgeHftRatio[i + 1]))
         return i;
   }
   return nPhisHftRatio - 1 ;
}

TVector3 smearPosData(int const iParticleIndex, double const vz, int zdcb, TLorentzVector const& rMom, TVector3 const& pos)
{
   int const iEtaIndex = getEtaIndexDca(rMom.PseudoRapidity());
   int const iVzIndex = getVzIndexDca(vz);
   // int const iPhiIndex = getPhiIndexDca(rMom.Phi());
   int const iPtIndex = getPtIndexDca(rMom.Perp());

   double sigmaPosZ = 0;
   double sigmaPosXY = 0;

   h2Dca[iParticleIndex][iEtaIndex][iVzIndex][zdcb][iPtIndex]->GetRandom2(sigmaPosXY,sigmaPosZ);
   sigmaPosZ *= 1.e4;
   sigmaPosXY *= 1.e4;

   TVector3 newPos(pos);
   newPos.SetZ(0);
   TVector3 momPerp(-rMom.Vect().Y(), rMom.Vect().X(), 0.0);
   newPos -= momPerp.Unit() * sigmaPosXY;

   return TVector3(newPos.X(), newPos.Y(), pos.Z() + sigmaPosZ);
}

TVector3 getVertex(int const centrality)
{
   double rdmVz;

   if (h1Vz[centrality]->GetEntries() == 0) rdmVz = 0.;
   else
   {
      do {
          rdmVz = h1Vz[centrality]->GetRandom() * 1e4;
          cout<<rdmVz<<endl;
      }
      while (fabs(rdmVz) > gVzCut);
   }

   return TVector3(0., 0., rdmVz);
}

int getZdcBin(int const centrality) 
{
    float zdc;
   int zdcbin=-1;
   while (zdcbin<0 || zdcbin>=nZdcX) {
       cout<<zdcbin<<endl;
      zdc = h1ZdcX[centrality]->GetRandom();
      zdcbin = h1ZdcX[centrality]->FindBin(zdc)-1;
   }
   return zdcbin;
}

bool reconstructD0(int const centrality, TLorentzVector const& mom)
{
   /*TGraph* gr = NULL;

   if (centrality < 4) gr = grEff[0];
   else if (centrality < 7) gr = grEff[1];
   else gr = grEff[2];

   return gRandom->Rndm() < gr->Eval(mom.Perp());*/
   return 1;
}

bool tpcReconstructed(int iParticleIndex, float charge, int cent, TLorentzVector const& mom)
{
   TH1D* h = NULL;

   if (iParticleIndex == 0)
   {
      if (charge > 0) h = hTpcPiPlus[cent];
      else h = hTpcPiMinus[cent];
   }
   else
   {
      if (charge > 0) h = hTpcKPlus[cent];
      else h = hTpcKMinus[cent];
   }

   int const bin = h->FindBin(mom.Perp());

   return gRandom->Rndm() < h->GetBinContent(bin);
}

bool matchHft(int const iParticleIndex, double const vz, int const zdcb, TLorentzVector const& mom)
{
   int const iEtaIndex = getEtaIndexHftRatio(mom.PseudoRapidity());
   int const iVzIndex = getVzIndexHftRatio(vz);
   int const iPhiIndex = getPhiIndexHftRatio(mom.Phi());

   int const bin = hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][zdcb]->FindBin(mom.Perp());
   return gRandom->Rndm() < hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][zdcb]->GetBinContent(bin);
}

bool matchTOF(int const iParticleIndex, TLorentzVector const& mom)
{
   if (iParticleIndex == 0) { // pion
      return gRandom->Rndm() < h_pi_tof_eff->GetBinContent(h_pi_tof_eff->FindBin(mom.Perp())); //from histogram
   }
   else if (iParticleIndex == 1) { // kaon
      return gRandom->Rndm() < h_k_tof_eff->GetBinContent(h_k_tof_eff->FindBin(mom.Perp())); //from histogram
   } else {
      return false;
   }
}

//___________
void bookObjects()
{
   cout << "Loading input momentum resolution ..." << endl;
   TFile fPionMom("pion_momentum_resolution.root");
   TFile fKaonMom("kaon_momentum_resolution.root");
   fPionMomResolution = (TF1*)fPionMom.Get("pion_MomResFit")->Clone("pion_MomResFit");
   fKaonMomResolution = (TF1*)fKaonMom.Get("kaon_MomResFit")->Clone("kaon_MomResFit");
   fPionMom.Close();
   fKaonMom.Close();

   cout << "Loading input spectra ..." << endl;
   TFile fPP("pp200_spectra.root");
   fWeightFunction = (TF1*)fPP.Get("run12/f1Levy")->Clone("f1Levy");
   fPP.Close();
   
   TFile fAuAu("Run14_D0_MyRaa_pT1.0.root");
   fWeightFunctionAuAu = (TF1*)fAuAu.Get("myLevyFcn_9")->Clone("f1LevyAuAu");
   fAuAu.Close();

   TFile fHftRatio1Pion("hftratio_vs_pt_dAu_pion.root");
   TFile fHftRatio1Kaon("hftratio_vs_pt_dAu_kaon.root");
   TFile fDca1("Dca2D_AuAu2016_lumiprod.root");

   TFile fEvent("inputs.event.root");

   TH3F* mh3VzZdcMult = (TH3F*)fEvent.Get("mh3VzZdcMult");

   TFile f_tof("eff_tof.root");
   h_pi_tof_eff = (TH1D*)f_tof.Get("eff_TOF_p0_nsigma1");
   h_pi_tof_eff->SetDirectory(0);
   h_k_tof_eff = (TH1D*)f_tof.Get("eff_TOF_p1_nsigma1");
   h_k_tof_eff->SetDirectory(0);
   f_tof.Close();

   char name[500];
   //getting VZ histogram for each multiplicity
   for (int ii = 0; ii < nmultEdge; ++ii)   {
       int binVzmin = 0;
       int binVzup = mh3VzZdcMult->GetXaxis()->GetNbins();
       int binZDCmin = 0;
       int binZDCmax = mh3VzZdcMult->GetYaxis()->GetNbins();
       int binMultmin = mh3VzZdcMult->GetZaxis()->FindBin(multEdge[ii]);
       int binMultmax = mh3VzZdcMult->GetZaxis()->FindBin(multEdge[ii+1]);
       h1Vz[ii] = mh3VzZdcMult -> ProjectionX("_px",binZDCmin, binZDCmax, binMultmin, binMultmax, ""); //vz zdc
       h1Vz[ii]->SetDirectory(0);
       h1ZdcX[ii] = mh3VzZdcMult -> ProjectionY("_py",binVzmin, binVzup, binMultmin, binMultmax, ""); //vz zdc
       h1ZdcX[ii]->SetDirectory(0);
   }
   fEvent.Close();

   cout << "Loading input HFT ratios and DCA ..." << endl;

   for (int iParticle = 0; iParticle < nParticles; ++iParticle)
   {
      for (int iZdc = 0; iZdc < m_nZdc; ++iZdc)
      {
         for (int iEta = 0; iEta < nEtasHftRatio; ++iEta)
         {
            for (int iVz = 0; iVz < nVzsHftRatio; ++iVz)
            {
               for (int iPhi = 0; iPhi < nPhisHftRatio; ++iPhi)
               {
                  if (iParticle==0)  hHftRatio1[iParticle][iEta][iVz][iPhi][iZdc] = (TH1D*)(fHftRatio1Pion.Get(Form("h_hftratio_p%d_eta%d_vz%d_phi%d_z%d", iParticle, iEta, iVz, iPhi, iZdc)));
                  if (iParticle==1)  hHftRatio1[iParticle][iEta][iVz][iPhi][iZdc] = (TH1D*)(fHftRatio1Kaon.Get(Form("h_hftratio_p%d_eta%d_vz%d_phi%d_z%d", iParticle, iEta, iVz, iPhi, iZdc)));

                  hHftRatio1[iParticle][iEta][iVz][iPhi][iZdc]->SetDirectory(0);
               }
            }
         }
      }
      cout << "Finished loading HFT Ratio: " <<  endl;

      for(int iZdc = 0; iZdc < nZdcX; ++iZdc)
      {
      // DCA
      for (int iEta = 0; iEta < nEtasDca; ++iEta)
      {
         for (int iVz = 0; iVz < nVzsDca; ++iVz)
         {
            // for (int iPhi = 0; iPhi < nPhisDca; ++iPhi)
               for (int iPt = 0; iPt < nPtBinsDca; ++iPt)
               {
                  h2Dca[iParticle][iEta][iVz][iZdc][iPt] = (TH2D*)((fDca1.Get(Form("h_dcaXYvsZ_c%i_p%i_z%i_vz%i_eta%i_pt%i", centralitySelect, iParticle, iZdc, iVz, iEta, iPt))));
                  h2Dca[iParticle][iEta][iVz][iZdc][iPt]->SetDirectory(0);
               }
            }
         }
      }
      // cout << "Finished loading centrality: " << iCent << endl;
      }
      cout << "Finished loading Dca: " <<  endl;

   fHftRatio1Pion.Close();
   fHftRatio1Kaon.Close();
   fDca1.Close();

   cout << " Loading TPC tracking efficiencies " << endl;

   TFile fTpcPiPlus("Eff_PionPlus_embedding.root");
   TFile fTpcPiMinus("Eff_PionMinus_embedding.root");
   TFile fTpcKPlus("Eff_KaonPlus_embedding.root");
   TFile fTpcKMinus("Eff_KaonMinus_embedding.root");

   for (int iCent = 0; iCent < nmultEdge; ++iCent)
   {
      hTpcPiPlus[iCent] = (TH1D*)fTpcPiPlus.Get(Form("h1Ratiocent_%i", iCent));
      hTpcPiPlus[iCent]->SetDirectory(0);
      hTpcPiMinus[iCent] = (TH1D*)fTpcPiMinus.Get(Form("h1Ratiocent_%i", iCent));
      hTpcPiMinus[iCent] ->SetDirectory(0);
      hTpcKPlus[iCent] = (TH1D*)fTpcKPlus.Get(Form("h1Ratiocent_%i", iCent));
      hTpcKPlus[iCent]->SetDirectory(0);
      hTpcKMinus[iCent] = (TH1D*)fTpcKMinus.Get(Form("h1Ratiocent_%i", iCent));
      hTpcKMinus[iCent]->SetDirectory(0);
   }

   fTpcPiPlus.Close();
   fTpcPiMinus.Close();
   fTpcKPlus.Close();
   fTpcKMinus.Close();

   cout << "Done with loading all files ..." << endl;

   /*grEff[0] = new TGraph("eff_4080.csv", "%lg %lg", ",");
   grEff[1] = new TGraph("eff_1040.csv", "%lg %lg", ",");
   grEff[2] = new TGraph("eff_010.csv", "%lg %lg", ",");*/

   std::stringstream ss_cent, ss_indx; ss_cent << centralitySelect; ss_indx << jobindx;
   result = new TFile((outFileName + "_c" + ss_cent.str() + "_job" + ss_indx.str() + ".root").c_str(), "recreate");
   result->SetCompressionLevel(1);
   result->cd();

   int BufSize = (int)pow(2., 16.);
   // int Split = 1;
   nt = new TNtuple("nt", "", "cent:zdcxbin:vx:vy:vz:vzIdx:"
                    "pid:w:m:pt:eta:y:phi:v0x:v0y:v0z:" // MC D0
                    "rM:rPt:rEta:rY:rPhi:rV0x:rV0y:rV0z:" // Rc D0
                    "dca12:decayLength:dcaD0ToPv:cosTheta:angle12:cosThetaStar:" // Rc pair
                    "kM:kPt:kEta:kY:kPhi:kDca:" // MC Kaon
                    "kRM:kRPt:kREta:kRY:kRPhi:kRVx:kRVy:kRVz:kRDca:kRSDca:kRDcaXY:kRDcaZ:kTpc:" // Rc Kaon
                    "pM:pPt:pEta:pY:pPhi:pDca:" // MC Pion1
                    "pRM:pRPt:pREta:pRY:pRPhi:pRVx:pRVy:pRVz:pRDca:pRSDca:pRDcaXY:pRDcaZ:pTpc:" // Rc Pion1
                    "kHft:pHft", BufSize);
   // nt->SetAutoSave(-500000); // autosave every 1 Mbytes
}
//___________
void write()
{
   result->cd();
   nt->Write();
   result->Close();
}
