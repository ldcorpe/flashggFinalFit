{
TFile *_file0 = TFile::Open("/afs/cern.ch/user/m/musella/public/workspace/exo/full_analysis_spring15_7415v2_sync_v5_data_ecorr_cic2_final_ws.root");
t_EBEE=_file0->Get("tree_data_cic2_EBEE");
t_EBEB=_file0->Get("tree_data_cic2_EBEB");
//TCanvas *t = new TCanvas ("functionLab","functionLab",800,400);
//t->Divide(2,1);
gStyle->SetOptStat(11111);
gStyle->SetOptFit(11111);

//t->cd(1);
t_EBEE->Draw("mgg>>h_EBEE(100,320,1600)");
TH1F *h_EBEE = (TH1F*)gPad->GetPrimitive("h_EBEE");
h_EBEE->SetMarkerColor(kBlack);
h_EBEE->SetMarkerSize(4);
h_EBEE->Draw("EP");

//t->cd(2);
t_EBEB->Draw("mgg>>h_EBEB(100,230,1600)");
TH1F *h_EBEB = (TH1F*)gPad->GetPrimitive("h_EBEB");
h_EBEB->SetMarkerColor(kBlack);
h_EBEB->Draw("EP");

std::cout << " --------------------------------------------------------------------------- " << std::endl;
std::cout << " Usage : The two categories considered in this analysis are EBEE and EBEB.   " << std::endl;
std::cout << " The mgg dists are loaded with correct ranges as TH1Fs : h_EBEB and h_EBEE   " << std::endl;
std::cout << " You can try to fit them using CINT the usual format for ROOT fitting:       " << std::endl;
std::cout << " h_EBEE->Fit(\"expo\")                                                       " << std::endl;
std::cout << " or, define arbitrary functions, using any number of params [0],[1],[2]...:  " << std::endl;
std::cout << " TF1 *f1 = new TF1(\"f1\",\"[0]*exp(-x*[1])\", 0, 1600);                     " << std::endl;
std::cout << " h_EBEE->Fit(\"f1\")                                                         " << std::endl;
std::cout << " --------------------------------------------------------------------------- " << std::endl;


}


