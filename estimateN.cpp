//
// Created by lukas on 20.6.2019.
//

void estimateN() {
//    TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,10,700,900);

    TFile *inFile = new TFile("pp200_spectra.root", "READ");

    TF1* f =  (TF1*)inFile->Get("run12/f1Levy")->Clone("f1Levy");
//    f->SetName("f1Levy");

    TH1F* h1f = new TH1F("h1f","Test",6,0.5,6);
    TH1F* h2 = new TH1F("h2","Test2",10,0.5,10);

    h1f->FillRandom("f1Levy",2000000);
//    h1f->FillRandom("sqroot",10000);

    TCanvas *c1 = new TCanvas("c1","c1",900,900);
//    f->Draw();
    h1f->Draw();

    cout<<"pt12:"<<endl;
    cout<<0.002*h1f->GetBinContent(2)<<" "<<h1f->GetBinContent(2)<<endl;
    float pt12Reco = 0.002*h1f->GetBinContent(2);
    float pt12Orig = h1f->GetBinContent(2);

    cout<<"pt45:"<<endl;
    cout<<0.025*h1f->GetBinContent(h1f->FindBin(4.5))<<" "<<h1f->GetBinContent(h1f->FindBin(4.5))<<endl;
    float pt45Reco = 0.025*h1f->GetBinContent(h1f->FindBin(4.5));
    float pt45Orig = h1f->GetBinContent(h1f->FindBin(4.5));

    float pt12 = sqrt( pow(sqrt(pt12Reco)/(pt12Orig),2) + pow(sqrt(pt12Orig)*pt12Reco/((pow(pt12Orig,2))),2) )/0.002;
    float pt45 = sqrt( pow(sqrt(pt45Reco)/(pt45Orig),2) + pow(sqrt(pt45Orig)*pt45Reco/((pow(pt45Orig,2))),2) )/0.025;

    cout<<"Stat. uncer. pT 1-2: "<<pt12*100<<endl;
    cout<<"Stat. uncer. pT 4-5: "<<pt45*100<<endl;

    TF1* f2 = new TF1("f2", "exp(-1.45-1.73*x)", 0, 10);

    TCanvas *c2 = new TCanvas("c2","c2",900,900);
    f->Draw();
    f2->Draw("same");




}