#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH3.h"
#include "TH2.h"
#include "TH1.h"
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
#include "StAnaCutsHijing.h"
// #include "TSystem.h"
// #include "TMemStat.h"

using namespace std;

void setDecayChannels(int const mdme);
void decayAndFill(int const kf, TLorentzVector* b, double const weight, TClonesArray& daughters);
void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& piMom, TVector3 v00);
void getKinematics(TLorentzVector& b, double const mass);
TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution);
TVector3 smearPosData(int iParticleIndex, double vz, int zdcb, TLorentzVector const& rMom, TVector3 const& pos, int const centrality);
float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0);
TVector3 getVertex(int centrality);
TVector3 getVertexWithError(int centrality);
int getZdcBin(int const centrality);
bool matchHft(int iParticleIndex, double vz, int zdcb, TLorentzVector const& mom);
bool tpcReconstructed(int iParticleIndex, float charge, int cent, TLorentzVector const& mom);
bool matchTOF(int const iParticleIndex, TLorentzVector const& mom);
bool goodPID(int const iParticleIndex, TLorentzVector const& mom);
void bookObjects();
void write();

int getPtIndexDca(double);
int getEtaIndexDca(double);
int getVzIndexDca(double);
int getPhiIndexDca(double);

int getVzIndexHftRatio(double);
int getPhiIndexHftRatio(double);
int getZdcBinRatio(float);
int getMultiplicityBin(double);
int getMultiplicityBinTPC(double);
int getMultiplicityBinDCA(double);

TVector3 const smearVertex(TVector3);

TPythia6Decayer* pydecay;
TNtuple* nt;
TFile* result;

TF1* fKaonMomResolution = NULL;
TF1* fPionMomResolution = NULL;
TF1* fWeightFunction = NULL;
TF1* f1PidPi = NULL;
TF1* f1PidK = NULL;
TGraph* grEff[3];
const Int_t nParticles = 2;
const Int_t nCentHftRatio = 9;

const int nmultEdge = 7;
float const multEdge[nmultEdge+1] = {0, 4, 8, 12, 16, 20, 24, 200};

//tpc Eff.
const int nmultEdgeTPC = 1;
float const multEdgeTPC[nmultEdgeTPC+1] = {0, 200};




// HFT ratio binning
const Int_t m_nEtasRatio = 6;
const Int_t nVzsHftRatio = 3;
const Int_t nPtBinsHftRatio = 15;
const Int_t nPhisHftRatio = 11;

const Double_t VzEdgeHftRatio[nVzsHftRatio + 1] = //ok
        {
                -6.0, -2.0, 2.0, 6.0
        };

const Double_t PhiEdgeHftRatio[nPhisHftRatio + 1] = //ok
        {
                -3.14159 , -2.80359 , -2.17527 , -1.54696 , -0.918637 , -0.290319 , 0.338 , 0.966319 , 1.59464 , 2.22296 , 2.85127 , 3.14159
        };

int const nVzsDca = 4;
const Int_t nPhisDca = 11;
int const nEtasDca = 4;
const Int_t nPtBinsDca = 7;
float const VzEdgeDca[nVzsDca + 1] = {   -6., -3., 0, 3., 6.};
float const EtaEdgeDca[nEtasDca + 1] = {0, 0.2, 0.4, 0.8, 1};

const Double_t PhiEdgeDca[nPhisDca + 1] =
        {
                -3.14159 , -2.80359 , -2.17527 , -1.54696 , -0.918637 , -0.290319 , 0.338 , 0.966319 , 1.59464 , 2.22296 , 2.85127 , 3.14159 //Sector by Sector
                // sector number 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3
        };
const Double_t ptEdgeDca[nPtBinsDca + 1] =
        {
                0.15, 0.4, 0.8, 1., 1.5, 2., 4., 12.
        };

const int nmultEdgeDCA = 1;
float const multEdgeDCA[nmultEdgeDCA+1] = {0, 200};

const int nZdcDCA = 1;
float const zdcxBinsDCA[nZdcDCA+1] = {0,210};

const int m_nZdc = 1;
float const m_zdcEdge[m_nZdc+1] = {0,210};

//

TH1D* h1Vz[nmultEdge+1];
TH1D* h1ZdcX[nmultEdge+1];
TH1D* h1VxError[nmultEdge+1];
TH1D* h1VzError[nmultEdge+1];

TH1F* hD0ptHIJING;
TH1F* hD0yHIJING;

TH1D* hHftRatio1[vars::m_nParticles][vars::m_nEtasRatio][vars::m_nVzsRatio][vars::m_nPhisRatio][vars::m_nZdc];
int const nCentDca = 9;
//TH2D* h2Dca[nParticles][nEtasDca][nVzsDca][nmultEdgeDCA][nPtBinsDca];
//TH2F* h2Dca[vars::m_nParticles][vars::m_nEtasDca][vars::m_nVzsDca][vars::m_nmultEdgeDCA][vars::m_nPtsDca];
TH2F* h2Dca[2][3][4][1][11];

TH1D* hTpcPiPlus[nmultEdgeTPC]; //embedding
TH1D* hTpcPiMinus[nmultEdgeTPC]; //embedding
TH1D* hTpcKPlus[nmultEdgeTPC]; //embedding
TH1D* hTpcKMinus[nmultEdgeTPC]; //embedding
TH1D* hRefMult;
TH1D* h_k_tof_eff;//embedding
TH1D* h_pi_tof_eff;//embedding

string outFileName = "D0.toyMc";
std::pair<int, int> const decayChannels(747, 807);
std::pair<float, float> const momentumRange(0, 10);

float const gVzCut = 6.0;
float const acceptanceRapidity = 1.0;
float const sigmaVertexCent[nCentHftRatio] = {31., 18.1, 12.8, 9.3, 7.2, 5.9, 5., 4.6, 4.}; //not using

//============== main  program ==================
int jobindx;

//_______________________________________________________________________________________________________________
void toyMcEffZeroDecayLength_hijingInputs(int npart = 1e5, int jobId=0)
//void toyMcEffZeroDecayLength(int npart = 500)
{
    jobindx = jobId;

    TStopwatch*   stopWatch = new TStopwatch();
    stopWatch->Start();

    gRandom->SetSeed(jobId);
    bookObjects();

    pydecay = TPythia6Decayer::Instance();
    pydecay->Init();
    setDecayChannels(763); // D0 --> Kpi
    TLorentzVector* b_d = new TLorentzVector;
    TClonesArray ptl("TParticle", 10);
    for (int ipart = 0; ipart < npart; ipart++)
    {
//      if (i((float)part/(float)npart*10) % 10 == 0)
//         cout << "____________ ipart = " << ipart / static_cast<float>(npart) << " ________________" << endl;

        getKinematics(*b_d, M_D_0); //random pt, y, phi, return b vector with correct prop.
        decayAndFill(421, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);  //421 = D0, ptl = daughters array
        decayAndFill(-421, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);
        if (ipart % 1000 == 1) nt->AutoSave("SaveSelf");
    }
    cout<<"lets write"<<endl;
    write();
    // mem.Show();
    stopWatch->Stop();
    stopWatch->Print();
}

//_______________________________________________________________________________________________________________
void setDecayChannels(int const mdme)
{
    for (int idc = decayChannels.first; idc < decayChannels.second + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0); //747-807 switching off
    TPythia6::Instance()->SetMDME(mdme, 1, 1);
}

//_______________________________________________________________________________________________________________
void decayAndFill(int const kf, TLorentzVector* b, double const weight, TClonesArray& daughters)
{
//    cout<<"decayAndFill start"<<endl;
    pydecay->Decay(kf, b);
    pydecay->ImportParticles(&daughters);

    TLorentzVector kMom;
    TLorentzVector pMom;
    TVector3 v00;

    int nTrk = daughters.GetEntriesFast();
    for (int iTrk = 0; iTrk < nTrk; ++iTrk)
    {
        TParticle* ptl0 = (TParticle*)daughters.At(iTrk);

        switch (abs(ptl0->GetPdgCode()))
        {
            case 321: //kaonplus kaon minus
                ptl0->Momentum(kMom); //seting momentum to kMom
                // v00.SetXYZ(0,0,0);
                v00.SetXYZ(ptl0->Vx() * 0.1, ptl0->Vy() * 0.1, ptl0->Vz() * 0.1); // converted to μm, production vertex
                break;
            case 211:   // pionplus
                ptl0->Momentum(pMom);
                break;
            default:
                break;
        }
    }
    daughters.Clear();
    fill(kf, b, weight, kMom, pMom, v00);
}

//_______________________________________________________________________________________________________________
void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& pMom, TVector3 v00)
{
//    cout<<"Fill() start"<<endl;
    Double_t refMult = hRefMult->GetRandom();
    int centrality = getMultiplicityBin(refMult);
    int centralityTPC = getMultiplicityBinTPC(refMult);
    int centralityDCA = getMultiplicityBinDCA(refMult);

//    int const centrality = floor(nmultEdge * gRandom->Rndm());

    TVector3 const vertex = getVertex(centrality); //from hVz,  converted to um
//    TVector3 const vertex = getVertexWithError(centrality); //from hVz,  converted to um
//    TVector3 const vertexR = smearVertex(vertex); not implemented yet
    int zdcb = getZdcBin(centrality); //from data

    // smear primary vertex
    // float const sigmaVertex = sigmaVertexCent[cent];
    // TVector3 const vertex(gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex));

    v00 += vertex; //SV + z of vertex from data, in um => Z position of the SV, nothing changed in xy

    // smear momentum
    TLorentzVector const kRMom = smearMom(kMom, fKaonMomResolution); //fKaonMomResolution is TF1
    TLorentzVector const pRMom = smearMom(pMom, fPionMomResolution);

    // smear position
    TVector3 const kRPos = smearPosData(1, vertex.z(), zdcb, kRMom, v00, centralityDCA); //particle dca smearing , transverse to its vector (why not to just change xy and z according to the dcaxy and dcaz from data without transverse position)
    TVector3 const pRPos = smearPosData(0, vertex.z(), zdcb, pRMom, v00, centralityDCA);

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
    float dca12 = dca1To2(kRMom.Vect(), kRPos, pRMom.Vect(), pRPos, v0); // v0 is reconstructed SV position
    float decayLength = (v0 - vertex).Mag();
    float dcaD0ToPv = dca(rMom.Vect(), v0, vertex);
    float const cosTheta = (v0 - vertex).Unit().Dot(rMom.Vect().Unit());
    float const angle12 = kRMom.Vect().Angle(pRMom.Vect());

    float cosThetaStar = kRMom.Vect().Unit().Dot(rMom.Vect().Unit());
    if (cosThetaStar!=cosThetaStar) cosThetaStar=-999;
    int const charge = kf > 0 ? 1 : -1;

    // save; change units from um to cm
//    dca12 /= 1e4;
//    decayLength /= 1e4;
//    dcaD0ToPv /= 1e4;
//    pRDca /= 1e4;
//    kRDca /= 1e4;

    const int nNtVars = nt->GetNvar();
    float arr[nNtVars];
    int iArr = 0;
    arr[iArr++] = centrality;
    arr[iArr++] = refMult;
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
    arr[iArr++] = tpcReconstructed(1, -1 * charge, centralityTPC, kRMom);

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
    arr[iArr++] = tpcReconstructed(0, charge, centralityTPC, pRMom);

    arr[iArr++] = goodPID(1, kRMom);
    arr[iArr++] = goodPID(0, pRMom);
    arr[iArr++] = matchHft(1, vertex.z(), zdcb, kRMom); //kaon = 1, pion = 0
    arr[iArr++] = matchHft(0, vertex.z(), zdcb, pRMom);

    nt->Fill(arr);
}

//_______________________________________________________________________________________________________________
void getKinematics(TLorentzVector& b, double const mass)
{
//    float const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
//    float const pt = fWeightFunction->GetRandom();
//    float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);

    float pt = 999;
    while (pt > 8)  pt = hD0ptHIJING->GetRandom();

    float y = 999;
    while (fabs(y)>1) y = hD0yHIJING->GetRandom();

    float const phi = TMath::TwoPi() * gRandom->Rndm(); //rndm [0,1]

    float const mT = sqrt(mass * mass + pt * pt);
    float const pz = mT * sinh(y);
    float const E = mT * cosh(y);

    b.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
}

//_______________________________________________________________________________________________________________
float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
    TVector3 posDiff = pos - vertex;
    return fabs(p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff));
}

//_______________________________________________________________________________________________________________
float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
    TVector3 posDiff = pos - vertex;
    float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;

    return sign * p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff);
}

//_______________________________________________________________________________________________________________
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

//_______________________________________________________________________________________________________________
float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
    TVector3 posDiff = pos - vertex;
    if (sin(p.Theta()) == 0) return 0;
    else return (-posDiff.x() * cos(p.Phi()) - posDiff.y() * sin(p.Phi())) * cos(p.Theta()) / sin(p.Theta()) + posDiff.z();
}

//_______________________________________________________________________________________________________________
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

//_______________________________________________________________________________________________________________
TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution)
{
    float const pt = b.Perp();
    float const sPt = gRandom->Gaus(pt, pt * fMomResolution->Eval(pt));

    TLorentzVector sMom;
    sMom.SetXYZM(sPt * cos(b.Phi()), sPt * sin(b.Phi()), sPt * sinh(b.PseudoRapidity()), b.M());
    return sMom;
}

//_______________________________________________________________________________________________________________
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

//_______________________________________________________________________________________________________________
int getEtaIndexDca(double Eta)
{
    for (int i = 0; i < nEtasDca; i++)
    {
        if ((fabs(Eta) >= EtaEdgeDca[i]) && (fabs(Eta) < EtaEdgeDca[i + 1]))
            return i;
    }
    return nEtasDca - 1 ;
}

//_______________________________________________________________________________________________________________
int getVzIndexDca(double Vz)
{
    for (int i = 0; i < nVzsDca; i++)
    {
        if ((Vz >= VzEdgeDca[i]) && (Vz < VzEdgeDca[i + 1]))
            return i;
    }
    return nVzsDca - 1 ;
}

//_______________________________________________________________________________________________________________
int getPhiIndexDca(double Phi)
{
    for (int i = 0; i < nPhisDca; i++)
    {
        if ((Phi >= PhiEdgeDca[i]) && (Phi < PhiEdgeDca[i + 1]))
            return i;
    }
    return nPhisDca - 1 ;
}

//_______________________________________________________________________________________________________________
int getVzIndexHftRatio(double Vz)
{
    for (int i = 0; i < nVzsHftRatio; i++)
    {
        if ((Vz >= VzEdgeHftRatio[i]) && (Vz < VzEdgeHftRatio[i + 1]))
            return i;
    }
    return -1 ;
}

//_______________________________________________________________________________________________________________
int getPhiIndexHftRatio(double Phi)
{
    for (int i = 0; i < nPhisHftRatio; i++)
    {
        if ((Phi >= PhiEdgeHftRatio[i]) && (Phi < PhiEdgeHftRatio[i + 1]))
            return i;
    }
    return -1 ;
}

//_______________________________________________________________________________________________________________
int getMultiplicityBin(double mult)
{
    for (int i = 0; i < nmultEdge; i++) {
        if ((mult >= multEdge[i]) && (mult < multEdge[i+1]))
            return i;
    }
    return -1 ;
}

//_______________________________________________________________________________________________________________
int getMultiplicityBinTPC(double mult)
{
    for (int i = 0; i < nmultEdgeTPC; i++) {
        if ((mult >= multEdgeTPC[i]) && (mult < multEdgeTPC[i+1]))
            return i;
    }
    return -1 ;
}

//_______________________________________________________________________________________________________________
int getMultiplicityBinDCA(double mult)
{
    for (int i = 0; i < nmultEdgeDCA; i++) {
        if ((mult >= multEdgeDCA[i]) && (mult < multEdgeDCA[i+1]))
            return i;
    }
    return -1 ;
}

//_______________________________________________________________________________________________________________
TVector3 smearPosData(int const iParticleIndex, double const vz, int zdcb, TLorentzVector const& rMom, TVector3 const& pos, int const centrality) //pos is SV
{
//    int const iEtaIndex = getEtaIndexDca(rMom.PseudoRapidity());
//    int const iVzIndex = getVzIndexDca(vz);
//    int const iPtIndex = getPtIndexDca(rMom.Perp());
//
    int const iEtaIndex = getIndex(rMom.PseudoRapidity(), vars::m_EtaEdgeDca, vars::m_nEtasDca);
    int const iVzIndex = getIndex(vz, vars::m_VzEdgeDca, vars::m_nVzsDca);
    int const iPtIndex = getIndex(rMom.Perp(), vars::m_PtEdgeDca, vars::m_nPtsDca);

    double sigmaPosZ = 0;
    double sigmaPosXY = 0;

    if (h2Dca[iParticleIndex][iEtaIndex][iVzIndex][centrality][iPtIndex]->GetEntries()>0) h2Dca[iParticleIndex][iEtaIndex][iVzIndex][centrality][iPtIndex]->GetRandom2(sigmaPosXY,sigmaPosZ);

    TVector3 newPos(pos);
    newPos.SetZ(0);
    TVector3 momPerp(-rMom.Vect().Y(), rMom.Vect().X(), 0.0); //momentum perp. to the original one in xy (transverse) plane
    newPos -= momPerp.Unit() * sigmaPosXY;

    return TVector3(newPos.X(), newPos.Y(), pos.Z() + sigmaPosZ);
}

//_______________________________________________________________________________________________________________
TVector3 getVertex(int const centrality)
{
    double rdmVz;

    if (h1Vz[centrality]->GetEntries() == 0) rdmVz = 0.;
    else
    {
        do {
            rdmVz = h1Vz[centrality]->GetRandom(); //um
        }
        while (fabs(rdmVz) > gVzCut);
    }

    return TVector3(0., 0., rdmVz);
}

//_______________________________________________________________________________________________________________
TVector3 getVertexWithError(int const centrality)
{
    double rdmVz;

    if (h1Vz[centrality]->GetEntries() == 0) rdmVz = 0.;
    else
    {
        do {
            rdmVz = h1Vz[centrality]->GetRandom();
        }
        while (fabs(rdmVz) > gVzCut);
    }

    double xError = h1VxError[centrality]->GetRandom();
    double yError = h1VxError[centrality]->GetRandom();
    double zError = h1VzError[centrality]->GetRandom();

    float rand;
    do {
        rand = gRandom->Uniform(-1, 1);
    } while ( rand==0 );
    if (rand<0) rand=-1;
    if (rand>0) rand=1;
    return TVector3(0.+rand*xError, 0.+rand*yError, rdmVz+rand*zError);
}

//_______________________________________________________________________________________________________________
int getZdcBin(int const centrality)
{
    float zdc;
    int zdcbin=-1;
    zdc = h1ZdcX[centrality]->GetRandom();
//    zdcbin = getZdcBinRatio(zdc);
    zdcbin = getIndex(zdc, vars::m_zdcEdge, vars::m_nZdc);

    return zdcbin;
}

//_______________________________________________________________________________________________________________
int getZdcBinRatio(float zdc){
    for (int i = 0; i < m_nZdc; i++){
        if ((zdc >= m_zdcEdge[i]) && (zdc < m_zdcEdge[i + 1]))
            return i;
    }
    return -1;
}

//_______________________________________________________________________________________________________________
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

//_______________________________________________________________________________________________________________
bool matchHft(int const iParticleIndex, double const vz, int const zdcb, TLorentzVector const& mom)
{
//    int const iEtaIndex = getEtaIndexHftRatio(mom.PseudoRapidity());
//    int const iVzIndex = getVzIndexHftRatio(vz);
//    int const iPhiIndex = getPhiIndexHftRatio(mom.Phi());

    int const iEtaIndex = getIndex(mom.PseudoRapidity(), vars::m_EtaEdgeRatio, m_nEtasRatio);
    int const iVzIndex = getIndex(vz, vars::m_VzEdgeRatio, vars::m_nVzsRatio);
    int const iPhiIndex = getIndex(mom.Phi(), vars::m_PhiEdgeRatio, vars::m_nPhisRatio);

    if (iEtaIndex<0 || iVzIndex<0 || iPhiIndex<0) return false;
    if (mom.Perp()>12) return false;
    if (mom.Perp()<0.15) return false;
    int bin = -1;
    if (hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][zdcb]->GetEntries()>0) bin = hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][zdcb]->FindBin(mom.Perp());
    else return false;
    if (bin<1) return false;
    return gRandom->Rndm() < hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][zdcb]->GetBinContent(bin);
}

//_______________________________________________________________________________________________________________
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

//_______________________________________________________________________________________________________________
bool goodPID(int const iParticleIndex, TLorentzVector const& mom)
{
//    cout<<"pid"<<endl;

    if (iParticleIndex == 0) { // pion
        return gRandom->Rndm() < f1PidPi->Eval(mom.Perp()); //from histogram
    }
    else if (iParticleIndex == 1) { // kaon
        return gRandom->Rndm() < f1PidK->Eval(mom.Perp()); //from histogram
    } else {
        return false;
    }
}

//_______________________________________________________________________________________________________________
TVector3 const smearVertex(TVector3 vertex){
    return vertex;
}


//_______________________________________________________________________________________________________________
void bookObjects()
{
    cout << "Loading input momentum resolution ..." << endl;
    TFile fPionMom("pion_momentum_resolution.root");
    TFile fKaonMom("kaon_momentum_resolution.root");
    fPionMomResolution = (TF1*)fPionMom.Get("fct_gaus_HFT_1.00")->Clone();
    fKaonMomResolution = (TF1*)fKaonMom.Get("fct_gaus_HFT_1.00")->Clone();
    fPionMom.Close();
    fKaonMom.Close();

    cout << "Loading input spectra ..." << endl;

    TFile fPP("published_run10_D0_AuAu_data.root");
    fWeightFunction = (TF1*)fPP.Get("Levy_pp")->Clone("f1Levy");
    fPP.Close();

    TFile filePidK("totalEff_K.root");
    f1PidK = (TF1*)filePidK.Get("fTotalGraphEffPid_K")->Clone("f1PidK");
    filePidK.Close();

    TFile filePidPi("totalEff_pi.root");
    f1PidPi = (TF1*)filePidPi.Get("fTotalGraphEffPid_pi")->Clone("f1PidPi");
    filePidPi.Close();

    TFile fileHijing("HIJING_D0_pt_y.root");
    hD0ptHIJING = (TH1F*)fileHijing.Get("hMcD0Pt")->Clone("hD0ptHIJING");
    hD0ptHIJING->SetDirectory(0);
    hD0yHIJING = (TH1F*)fileHijing.Get("hMcD0Y")->Clone("hD0yHIJING");
    hD0yHIJING->SetDirectory(0);
    fileHijing.Close();

    cout<<"Loading TOF eff."<<endl;
    TFile f_tof("eff_tof.root");
    h_pi_tof_eff = (TH1D*)f_tof.Get("eff_TOF_p0_nsigma1");
    h_pi_tof_eff->SetDirectory(0);
    h_k_tof_eff = (TH1D*)f_tof.Get("eff_TOF_p1_nsigma1");
    h_k_tof_eff->SetDirectory(0);
    f_tof.Close();

    char name[500];
    //getting VZ histogram for each multiplicity
    cout<<"Loading Vz and ZDCs..."<<endl;
    TFile fEvent("inputs.event_hijing.root");
    TH3D* mh3VzZdcMult = new TH3D();
    mh3VzZdcMult = (TH3D*)fEvent.Get("mh3VzZdcMult");
    mh3VzZdcMult->SetDirectory(0);

    int binVzmin = 1;
    int binVzup = mh3VzZdcMult->GetXaxis()->GetNbins();
    int binZDCmin = 1;
    int binZDCmax = mh3VzZdcMult->GetYaxis()->GetNbins();

    for (int ii = 0; ii < nmultEdge; ++ii)   {
        int binMultmin = mh3VzZdcMult->GetZaxis()->FindBin(multEdge[ii]);
        int binMultmax = mh3VzZdcMult->GetZaxis()->FindBin(multEdge[ii+1]);
        h1Vz[ii] = mh3VzZdcMult -> ProjectionX("_px",binZDCmin, binZDCmax, binMultmin, binMultmax, ""); //vz zdc
        h1Vz[ii]->SetDirectory(0);
        h1ZdcX[ii] = mh3VzZdcMult -> ProjectionY("_py",binVzmin, binVzup, binMultmin, binMultmax, ""); //vz zdc
        h1ZdcX[ii]->SetDirectory(0);
    }

    hRefMult = (TH1D*)fEvent.Get("hrefMult");
    hRefMult->SetDirectory(0);
    fEvent.Close();

    cout << "Loading input HFT ratios and DCA ...HIJING..." << endl;
    TFile fDca1("dcaxy_vs_dcaz_hijing.root");
    TFile fHftRatio1Pion("hftratio_vs_pt_dAu_pion_hijing.root");
    TFile fHftRatio1Kaon("hftratio_vs_pt_dAu_kaon_hijing.root");
    for (int iParticle = 0; iParticle < vars::m_nParticles; ++iParticle) {
        for (int iZdc = 0; iZdc < vars::m_nZdc; ++iZdc) {
            for (int iEta = 0; iEta < vars::m_nEtasRatio; ++iEta) {
                for (int iVz = 0; iVz < vars::m_nVzsRatio; ++iVz) {
                    for (int iPhi = 0; iPhi < vars::m_nPhisRatio; ++iPhi) {
                        if (iParticle==0)  hHftRatio1[iParticle][iEta][iVz][iPhi][iZdc] = (TH1D*)(fHftRatio1Pion.Get(Form("h_hftratio_p%d_eta%d_vz%d_phi%d_z%d", iParticle, iEta, iVz, iPhi, iZdc)));
                        if (iParticle==1)  hHftRatio1[iParticle][iEta][iVz][iPhi][iZdc] = (TH1D*)(fHftRatio1Kaon.Get(Form("h_hftratio_p%d_eta%d_vz%d_phi%d_z%d", iParticle, iEta, iVz, iPhi, iZdc)));
                        hHftRatio1[iParticle][iEta][iVz][iPhi][iZdc]->SetDirectory(0);
                    }
                }
            }
        }

        //DCA
        int iZdc=0;
        for (int iEta = 0; iEta < vars::m_nEtasDca; ++iEta) {
            for (int iVz = 0; iVz < vars::m_nVzsDca; ++iVz) {
                for (int iCent = 0; iCent < vars::m_nmultEdgeDCA; ++iCent) {
                    for (int iPt = 0; iPt < vars::m_nPtsDca; ++iPt) {
                        TString h2dName=Form("mh3DcaXyZPt_p%d_eta%d_vz%d_m%d_pt%d", iParticle, iEta, iVz, iCent, iPt);
                        h2Dca[iParticle][iEta][iVz][iCent][iPt] = (TH2F* )((fDca1.Get(h2dName)));
                        h2Dca[iParticle][iEta][iVz][iCent][iPt]->SetDirectory(0);
                    }
                }
            }
        }
    }
    cout << "Finished loading Dca: " <<  endl;

    fHftRatio1Pion.Close();
    fHftRatio1Kaon.Close();
    fDca1.Close();

    cout << " Loading TPC tracking efficiencies " << endl;
    TFile fTpcPiPlus("piplus_tpc_eff_embedding.root");
    TFile fTpcPiMinus("piminus_tpc_eff_embedding.root");
    TFile fTpcKPlus("kplus_tpc_eff_embedding.root");
    TFile fTpcKMinus("kminus_tpc_eff_embedding.root");

    for (int iCent = 0; iCent < nmultEdgeTPC; ++iCent) {
        hTpcPiPlus[iCent] = (TH1D*)fTpcPiPlus.Get(Form("TrackEffMult%i", iCent));
        hTpcPiPlus[iCent]->SetDirectory(0);
        hTpcPiMinus[iCent] = (TH1D*)fTpcPiMinus.Get(Form("TrackEffMult%i", iCent));
        hTpcPiMinus[iCent] ->SetDirectory(0);
        hTpcKPlus[iCent] = (TH1D*)fTpcKPlus.Get(Form("TrackEffMult%i", iCent));
        hTpcKPlus[iCent]->SetDirectory(0);
        hTpcKMinus[iCent] = (TH1D*)fTpcKMinus.Get(Form("TrackEffMult%i", iCent));
        hTpcKMinus[iCent]->SetDirectory(0);
    }

    fTpcPiPlus.Close();
    fTpcPiMinus.Close();
    fTpcKPlus.Close();
    fTpcKMinus.Close();

//   cout << "Done with loading all files ..." << endl;

    std::stringstream ss_indx; ss_indx << jobindx;
    result = new TFile((outFileName + "_c" + "_job" + ss_indx.str() + ".root").c_str(), "recreate");
    result->SetCompressionLevel(1);
    result->cd();

    int BufSize = (int)pow(2., 16.);
// int Split = 1;
    nt = new TNtuple("nt", "", "cent:refMult:zdcxbin:vx:vy:vz:vzIdx:"
                               "pid:w:m:pt:eta:y:phi:v0x:v0y:v0z:" // MC D0
                               "rM:rPt:rEta:rY:rPhi:rV0x:rV0y:rV0z:" // Rc D0
                               "dca12:decayLength:dcaD0ToPv:cosTheta:angle12:cosThetaStar:" // Rc pair
                               "kM:kPt:kEta:kY:kPhi:kDca:" // MC Kaon
                               "kRM:kRPt:kREta:kRY:kRPhi:kRVx:kRVy:kRVz:kRDca:kRSDca:kRDcaXY:kRDcaZ:kTpc:" // Rc Kaon
                               "pM:pPt:pEta:pY:pPhi:pDca:" // MC Pion1
                               "pRM:pRPt:pREta:pRY:pRPhi:pRVx:pRVy:pRVz:pRDca:pRSDca:pRDcaXY:pRDcaZ:pTpc:" // Rc Pion1
                               "kPID:pPID:kHft:pHft", BufSize);
// nt->SetAutoSave(-500000); // autosave every 1 Mbytes

    cout<<"files loaded"<<endl;
}

//___________
void write() {
    result->cd();
    nt->Write(nt->GetName(), TObject::kOverwrite);
    result->Close();
}
