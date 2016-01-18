#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"

#include "TFile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"
#include "RooMsgService.h"
#include "RooDataHist.h"
#include "RooExtendPdf.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TMacro.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TArrow.h"
#include "TKey.h"

#include "RooCategory.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"

#include "../interface/PdfModelBuilder.h"
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFunc.h>
#include <iomanip>
#include "boost/program_options.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"

using namespace std;
using namespace RooFit;
using namespace boost;

namespace po = program_options;

bool BLIND = true;
bool runFtestCheckWithToys=false;
int nBinsForMass = 100;

TRandom3 *RandomGen = new TRandom3();

RooAbsPdf* getPdf(PdfModelBuilder &pdfsModel, string type, int order, const char* ext=""){
  if (type=="DijetSimple") return pdfsModel.getDijetSimple(Form("%s_dijetsimple%d",ext,order),order); 
  else if (type=="Dijet") return pdfsModel.getDijet(Form("%s_dijet%d",ext,order),order); 
  else if (type=="Atlas") return pdfsModel.getAtlas(Form("%s_atlas%d",ext,order),order); 
  else if (type=="Expow") return pdfsModel.getExpow(Form("%s_expow%d",ext,order),order); 
  else if (type=="Chebychev") return pdfsModel.getChebychev(Form("%s_cheb%d",ext,order),order); 
  else if (type=="Exponential") return pdfsModel.getExponentialSingle(Form("%s_exp%d",ext,order),order); 
  else if (type=="Expow") return pdfsModel.getExpow(Form("%s_expow%d",ext,order),order); 
  else if (type=="PowerLaw") return pdfsModel.getPowerLawSingle(Form("%s_pow%d",ext,order),order); 
  else if (type=="Laurent") return pdfsModel.getLaurentSeries(Form("%s_lau%d",ext,order),order); 
  else {
    cerr << "[ERROR] -- getPdf() -- type " << type << " not recognised." << endl;
    return NULL;
  }
}

void runFit(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries){

	int ntries=0;
  	RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
	//std::cout << " BEFORE ITERATIONS-------------------------------" << std::endl;
	//params_test->Print("v");
	int stat=1;
	double minnll=10e8;
	while (stat!=0){
	  if (ntries>=MaxTries) break;
	  //std::cout << "----------------------------- BEFORE FIT-------------------------------" << std::endl;
	  //params_test->Print("v");
	  //std::cout << "-----------------------------------------------------------------------" << std::endl;
	  RooFitResult *fitTest = pdf->fitTo(*data,RooFit::Save(1)
		,RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME
          stat = fitTest->status();
	  minnll = fitTest->minNll();
	  if (stat!=0) params_test->assignValueOnly(fitTest->randomizePars());
	  ntries++; 
	}
	*stat_t = stat;
	*NLL = minnll;
}
double getProbabilityFtest(double chi2, int ndof,RooAbsPdf *pdfNull, RooAbsPdf *pdfTest, RooRealVar *mass, RooDataSet *data, std::string name){
 
  double prob_asym = TMath::Prob(chi2,ndof);
  if (!runFtestCheckWithToys) return prob_asym;

  int ndata = data->sumEntries();
  
  // fit the pdfs to the data and keep this fit Result (for randomizing)
  RooFitResult *fitNullData = pdfNull->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
		,RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE),RooFit::PrintLevel(-1)); //FIXME
  RooFitResult *fitTestData = pdfTest->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
		,RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE),RooFit::PrintLevel(-1)); //FIXME

  // Ok we want to check the distribution in toys then 
  // Step 1, cache the parameters of each pdf so as not to upset anything 
  RooArgSet *params_null = pdfNull->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_null;
  params_null->snapshot(preParams_null);
  RooArgSet *params_test = pdfTest->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_test;
  params_test->snapshot(preParams_test);
 
  int ntoys =5000;
  TCanvas *can = new TCanvas();
  TH1F toyhist(Form("toys_fTest_%s.pdf",pdfNull->GetName()),";Chi2;",60,-2,10);
  TH1I toyhistStatN(Form("Status_%s.pdf",pdfNull->GetName()),";FitStatus;",8,-4,4);
  TH1I toyhistStatT(Form("Status_%s.pdf",pdfTest->GetName()),";FitStatus;",8,-4,4);

  TGraph *gChi2 = new TGraph();
  gChi2->SetLineColor(kGreen+2);
  double w = toyhist.GetBinWidth(1);

  int ipoint=0;

  for (int b=0;b<toyhist.GetNbinsX();b++){
	double x = toyhist.GetBinCenter(b+1);
	if (x>0){
	  gChi2->SetPoint(ipoint,x,(ROOT::Math::chisquared_pdf(x,ndof)));
	  ipoint++;
	}
  }
  int npass =0; int nsuccesst =0;
  mass->setBins(nBinsForMass);
  for (int itoy = 0 ; itoy < ntoys ; itoy++){

        params_null->assignValueOnly(preParams_null);
        params_test->assignValueOnly(preParams_test);
  	RooDataHist *binnedtoy = pdfNull->generateBinned(RooArgSet(*mass),ndata,0,1);

	int stat_n=1;
        int stat_t=1;
	int ntries = 0;
	double nllNull,nllTest;
	// Iterate on the fit 
	int MaxTries = 2;
	while (stat_n!=0){
	  if (ntries>=MaxTries) break;
	  RooFitResult *fitNull = pdfNull->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1),RooFit::SumW2Error(kTRUE) //FIXME
		,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1));
		//,RooFit::Optimize(0));

	  nllNull = fitNull->minNll();
          stat_n = fitNull->status();
	  if (stat_n!=0) params_null->assignValueOnly(fitNullData->randomizePars());
	  ntries++; 
	}
	
	ntries = 0;
	while (stat_t!=0){
	  if (ntries>=MaxTries) break;
	  RooFitResult *fitTest = pdfTest->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1),RooFit::SumW2Error(kTRUE) //FIXME
		,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1));
	  nllTest = fitTest->minNll();
          stat_t = fitTest->status();
	  if (stat_t!=0) params_test->assignValueOnly(fitTestData->randomizePars()); 
	  ntries++; 
	}
       
	toyhistStatN.Fill(stat_n);
	toyhistStatT.Fill(stat_t);

        if (stat_t !=0 || stat_n !=0) continue;
	nsuccesst++;
	double chi2_t = 2*(nllNull-nllTest);
	if (chi2_t >= chi2) npass++;
        toyhist.Fill(chi2_t);
  }

  double prob=0;
  if (nsuccesst!=0)  prob = (double)npass / nsuccesst;
  toyhist.Scale(1./(w*toyhist.Integral()));
  toyhist.Draw();
  TArrow lData(chi2,toyhist.GetMaximum(),chi2,0);
  lData.SetLineWidth(2);
  lData.Draw();
  gChi2->Draw("L");
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->DrawLatex(0.1,0.91,Form("Prob (asymptotic) = %.4f (%.4f)",prob,prob_asym));
  can->SaveAs(name.c_str());

  TCanvas *stas =new TCanvas();
  toyhistStatT.SetLineColor(1); 
  TLegend *leg = new TLegend(0.2,0.6,0.4,0.89); leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(&toyhistStatN,"Null Hyp","L");
  leg->AddEntry(&toyhistStatT,"Test Hyp","L");
  toyhistStatN.Draw();
  toyhistStatT.Draw("same");
  leg->Draw();
  stas->SaveAs(Form("%s_fitstatus.pdf",name.c_str()));
  //reassign params
  params_null->assignValueOnly(preParams_null);
  params_test->assignValueOnly(preParams_test);

  delete can; delete stas;
  delete gChi2;
  delete leg;
  delete lat;

  // Still return the asymptotic prob (usually its close to the toys one)
  return prob_asym;

}

double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooDataSet *data, std::string name){
	
  double prob;
  int ntoys = 500;

  // Routine to calculate the goodness of fit. 
  name+="_gofTest.pdf";
  RooRealVar norm("norm","norm",data->sumEntries(),0,10E6);
  //norm.removeRange();

  RooExtendPdf *pdf = new RooExtendPdf("ext","ext",*mpdf,norm);

  // get The Chi2 value from the data
  RooPlot *plot_chi2 = mass->frame();
  data->plotOn(plot_chi2,Binning(nBinsForMass),Name("data"));

  pdf->plotOn(plot_chi2,Name("pdf"));
  int np = pdf->getParameters(*data)->getSize();

  double chi2 = plot_chi2->chiSquare("pdf","data",np);
  std::cout << "[INFO] Calculating GOF for pdf " << pdf->GetName() << ", using " <<np << " fitted parameters" <<std::endl;

  // The first thing is to check if the number of entries in any bin is < 5 
  // if so, we don't rely on asymptotic approximations
 
  if ((double)data->sumEntries()/nBinsForMass < 5 ){

    std::cout << "[INFO] Running toys for GOF test " << std::endl;
    // store pre-fit params 
    RooArgSet *params = pdf->getParameters(*data);
    RooArgSet preParams;
    params->snapshot(preParams);
    int ndata = data->sumEntries();
 
    int npass =0;
    std::vector<double> toy_chi2;
    for (int itoy = 0 ; itoy < ntoys ; itoy++){
    //  std::cout << "[INFO] " <<Form("\t.. %.1f %% complete\r",100*float(itoy)/ntoys) << std::flush;
      params->assignValueOnly(preParams);
      int nToyEvents = RandomGen->Poisson(ndata);
      RooDataHist *binnedtoy = pdf->generateBinned(RooArgSet(*mass),nToyEvents,0,1);
      pdf->fitTo(*binnedtoy,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Strategy(0),RooFit::SumW2Error(kTRUE)); //FIXME

      RooPlot *plot_t = mass->frame();
      binnedtoy->plotOn(plot_t);
      pdf->plotOn(plot_t);//,RooFit::NormRange("fitdata_1,fitdata_2"));

      double chi2_t = plot_t->chiSquare(np);
      if( chi2_t>=chi2) npass++;
      toy_chi2.push_back(chi2_t*(nBinsForMass-np));
      delete plot_t;
    }
    std::cout << "[INFO] complete" << std::endl;
    prob = (double)npass / ntoys;

    TCanvas *can = new TCanvas();
    can->SetLogy();
    can->SetLogx();
	double medianChi2 = toy_chi2[(int)(((float)ntoys)/2)];
    double rms = TMath::Sqrt(medianChi2);

    TH1F toyhist(Form("gofTest_%s.pdf",pdf->GetName()),";Chi2;",50,medianChi2-5*rms,medianChi2+5*rms);
    for (std::vector<double>::iterator itx = toy_chi2.begin();itx!=toy_chi2.end();itx++){
      toyhist.Fill((*itx));
    }
    toyhist.Draw();

    TArrow lData(chi2*(nBinsForMass-np),toyhist.GetMaximum(),chi2*(nBinsForMass-np),0);
    lData.SetLineWidth(2);
    lData.Draw();
    can->SaveAs(name.c_str());

    // back to best fit 	
    params->assignValueOnly(preParams);
  } else {
    prob = TMath::Prob(chi2*(nBinsForMass-np),nBinsForMass-np);
  }
  std::cout << "[INFO] Chi2 in Observed =  " << chi2*(nBinsForMass-np) << std::endl;
  std::cout << "[INFO] p-value  =  " << prob << std::endl;
  delete pdf;
  return prob;

}

void plot(RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, string name,vector<string> diphotonCats_, int status, double *prob){
  
  // Chi2 taken from full range fit
  RooPlot *plot_chi2 = mass->frame();
  data->plotOn(plot_chi2,Binning(nBinsForMass));
  pdf->plotOn(plot_chi2);

  int np = pdf->getParameters(*data)->getSize()+1; //Because this pdf has no extend
  double chi2 = plot_chi2->chiSquare(np);
 
  *prob = getGoodnessOfFit(mass,pdf,data,name);
  RooPlot *plot = mass->frame();
  mass->setRange("unblindReg_1",0,500);
  //mass->setRange("unblindReg_2",150,180);
  if (BLIND) {
    data->plotOn(plot,Binning(nBinsForMass),CutRange("unblindReg_1"));
    //data->plotOn(plot,Binning(nBinsForMass),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(nBinsForMass),Invisible());
  }
  else data->plotOn(plot,Binning(nBinsForMass));

 // data->plotOn(plot,Binning(80));
  TCanvas *canv = new TCanvas();
  canv->SetLogy();
  canv->SetLogx();
  pdf->plotOn(plot);//,RooFit::NormRange("fitdata_1,fitdata_2"));
  pdf->paramOn(plot,RooFit::Layout(0.34,0.96,0.89),RooFit::Format("NEA",AutoPrecision(1)));
  if (BLIND) plot->SetMinimum(0.0001);
  plot->SetTitle("");
  plot->Draw();

  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->DrawLatex(0.1,0.92,Form("#chi^{2} = %.3f, Prob = %.2f, Fit Status = %d ",chi2*(nBinsForMass-np),*prob,status));
  canv->SaveAs(name.c_str());
 	
	//plot_chi2->Draw();
  //canv->SaveAs((name+"debug").c_str());

  delete canv;
  delete lat;
}
void plot(RooRealVar *mass, RooMultiPdf *pdfs, RooCategory *catIndex, RooDataSet *data, string name, vector<string> diphotonCats_, int cat, int bestFitPdf=-1){
  
  int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
  TCanvas *canv = new TCanvas();
  canv->SetLogy();
  canv->SetLogx();
  TLegend *leg = new TLegend(0.5,0.55,0.92,0.92);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  RooPlot *plot = mass->frame();
		
		//data->plotOn(plot);
		//plot->Draw();
		///canv->Print(Form("test_LC.pdf"));

  mass->setRange("unblindReg_1",0,500);
  //mass->setRange("unblindReg_2",150,180);
  if (BLIND) {
    data->plotOn(plot,Binning(nBinsForMass),CutRange("unblindReg_1"));
    //data->plotOn(plot,Binning(80),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(nBinsForMass),Invisible());
  }
  else data->plotOn(plot,Binning(nBinsForMass)); 

  int currentIndex = catIndex->getIndex();
  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  leg->AddEntry(datLeg,Form("Data - %s",diphotonCats_[cat].c_str()),"LEP");
  int style=1;
  for (int icat=0;icat<catIndex->numTypes();icat++){
    int col;
    if (icat<=6) col=color[icat];
    else {col=kBlack; style++;}
    catIndex->setIndex(icat);
    pdfs->getCurrentPdf()->fitTo(*data,RooFit::Minos(0),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE));	 //FIXME
    pdfs->getCurrentPdf()->plotOn(plot,LineColor(col),LineStyle(style));//,RooFit::NormRange("fitdata_1,fitdata_2"));
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    if (bestFitPdf==icat) ext=" (Best Fit Pdf) ";
    leg->AddEntry(pdfLeg,Form("%s%s",pdfs->getCurrentPdf()->GetName(),ext.c_str()),"L");
  }
  plot->SetTitle(Form("Category %s",diphotonCats_[cat].c_str()));
  if (BLIND) plot->SetMinimum(0.0001);
  plot->Draw();
  leg->Draw("same");
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  canv->SaveAs(Form("%s.root",name.c_str()));
  catIndex->setIndex(currentIndex);
  delete canv;
}

void plot(RooRealVar *mass, map<string,RooAbsPdf*> pdfs, RooDataSet *data, string name, vector<string> diphotonCats_, int cat, int bestFitPdf=-1){
  
  int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
  TCanvas *canv = new TCanvas();
  canv->SetLogy();
  canv->SetLogx();
  TLegend *leg = new TLegend(0.6,0.65,0.89,0.89);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  RooPlot *plot = mass->frame();

  mass->setRange("unblindReg_1",0,500);
  //mass->setRange("unblindReg_2",150,180);
  if (BLIND) {
    data->plotOn(plot,Binning(nBinsForMass),CutRange("unblindReg_1"));
    //data->plotOn(plot,Binning(nBinsForMass),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(nBinsForMass),Invisible());
  }
  else data->plotOn(plot,Binning(nBinsForMass));

  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
	if(diphotonCats_.size() >0){
  leg->AddEntry(datLeg,Form("Data - %s",diphotonCats_[cat].c_str()),"LEP");
	} else {
  leg->AddEntry(datLeg,Form("Data - %d",cat),"LEP");
	}
  int i=0;
  int style=1;
  for (map<string,RooAbsPdf*>::iterator it=pdfs.begin(); it!=pdfs.end(); it++){
    int col;
    if (i<=6) col=color[i];
    else {col=kBlack; style++;}
    it->second->plotOn(plot,LineColor(col),LineStyle(style));//,RooFit::NormRange("fitdata_1,fitdata_2"));
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    if (bestFitPdf==i) ext=" (Best Fit Pdf) ";
    leg->AddEntry(pdfLeg,Form("%s%s",it->first.c_str(),ext.c_str()),"L");
    i++;
  }
  plot->SetTitle(Form(" %s",diphotonCats_[cat].c_str()));
  if (BLIND) plot->SetMinimum(0.0001);
  plot->Draw();
  leg->Draw("same");
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  canv->SaveAs(Form("%s.root",name.c_str()));
  delete canv;
}

void transferMacros(TFile *inFile, TFile *outFile){
  
  TIter next(inFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())){
    if (string(key->ReadObj()->ClassName())=="TMacro") {
      //cout << key->ReadObj()->ClassName() << " : " << key->GetName() << endl;
      TMacro *macro = (TMacro*)inFile->Get(key->GetName());
      outFile->cd();
      macro->Write();
    }
  }
}
int getBestFitFunction(RooMultiPdf *bkg, RooDataSet *data, RooCategory *cat, bool silent=false){


	double global_minNll = 1E10;
	int best_index = 0;
	int number_of_indeces = cat->numTypes();
		
	RooArgSet snap,clean;
	RooArgSet *params = bkg->getParameters((const RooArgSet*)0);
	params->remove(*cat);
	params->snapshot(snap);
	params->snapshot(clean);
	if (!silent) {
		std::cout << "[INFO] CLEAN SET OF PARAMETERS" << std::endl;
		//params->Print("V");
		std::cout << "-----------------------" << std::endl;
	}
	
	//bkg->setDirtyInhibit(1);
	//RooAbsReal *nllm = bkg->createNLL(*data);
	//RooMinimizer minim(*nllm);
	//minim.setStrategy(1);
	
	for (int id=0;id<number_of_indeces;id++){		
		params->assignValueOnly(clean);
		cat->setIndex(id);

		//RooAbsReal *nllm = bkg->getCurrentPdf()->createNLL(*data);

		if (!silent) {
			/*
			std::cout << "BEFORE  MAKING FIT" << std::endl;
			params->Print("V");
			std::cout << "-----------------------" << std::endl;		
			*/
		}
		
		//minim.minimize("Minuit2","minimize");
		double minNll=0; //(nllm->getVal())+bkg->getCorrection();
		int fitStatus=1;		
		runFit(bkg->getCurrentPdf(),data,&minNll,&fitStatus,/*max iterations*/3);
		// Add the penalty

		minNll=minNll+bkg->getCorrection();

		if (!silent) {
			/*
			std::cout << "After Minimization ------------------  " <<std::endl;
			std::cout << bkg->getCurrentPdf()->GetName() << " " << minNll <<std::endl;
			bkg->Print("v");
			bkg->getCurrentPdf()->getParameters(*data)->Print("V");
			std::cout << " ------------------------------------  " << std::endl;
	
			params->Print("V");
			*/
			std::cout << "[INFO] AFTER FITTING" << std::endl;
			std::cout << "[INFO] Function was " << bkg->getCurrentPdf()->GetName() <<std::endl;
			std::cout << "[INFO] Correction Applied is " << bkg->getCorrection() <<std::endl;
			std::cout << "[INFO] NLL + c = " <<  minNll << std::endl;
			std::cout << "-----------------------" << std::endl;
		}
			
		if (minNll < global_minNll){
        		global_minNll = minNll;
			snap.assignValueOnly(*params);
        		best_index=id;
		}
	}
    	cat->setIndex(best_index);
	params->assignValueOnly(snap);
	
	if (!silent) {
		std::cout << "[INFO] Best fit Function -- " << bkg->getCurrentPdf()->GetName() << " " << cat->getIndex() <<std::endl;
		//bkg->getCurrentPdf()->getParameters(*data)->Print("v");
	}
	return best_index;
}

int main(int argc, char* argv[]){
 
  string fileName;
  int ncats;
  string datfile;
  string outDir;
  string outfilename;
  bool verbose=false;
  bool saveMultiPdf=false;
string diphotonCatsStr_;
vector<string> diphotonCats_;
 bool isData_ =0;
 string sqrts_ = "13TeV";
 string analysisType_ = "";

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("infilename,i", po::value<string>(&fileName),                                              "In file name")
    ("ncats,c", po::value<int>(&ncats)->default_value(5),                                       "Number of categories")
    ("datfile,d", po::value<string>(&datfile)->default_value("dat/fTest.dat"),                  "Right results to datfile for BiasStudy")
    ("outDir,D", po::value<string>(&outDir)->default_value("plots/fTest"),                      "Out directory for plots")
    ("analysisType,a", po::value<string>(&analysisType_)->default_value("cic2"),                 "What kind of analysis you are rinning (cic, cic2...)... <Eventually this might also be able to directly populate the diphotonCats ?>")
    ("saveMultiPdf", po::value<string>(&outfilename),         					"Save a MultiPdf model with the appropriate pdfs")
    ("runFtestCheckWithToys", 									"When running the F-test, use toys to calculate pvals (and make plots) ")
    ("unblind",  									        "Dont blind plots")
    ("isData",    								    	        "Use Data not MC-based pseudodata ")
		("diphotonCats,f", po::value<string>(&diphotonCatsStr_)->default_value("EBEB,EBEE"),       "Flashgg category names to consider")
    ("verbose,v",                                                                               "Run with more output")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")) { cout << desc << endl; exit(1); }
	if (vm.count("unblind")) BLIND=false;
	if (vm.count("isData")) isData_=true;
  saveMultiPdf = vm.count("saveMultiPdf");

  if (vm.count("verbose")) verbose=true;
  if (vm.count("runFtestCheckWithToys")) runFtestCheckWithToys=true;

  if (!verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
    gErrorIgnoreLevel=kWarning;
  }
	split(diphotonCats_,diphotonCatsStr_,boost::is_any_of(","));
  
	int startingCategory=0;
	ncats= diphotonCats_.size();

  if(verbose) std::cout << "[INFO] SaveMultiPdf? " << saveMultiPdf << std::endl;
  TFile *outputfile;
  RooWorkspace *outputws;

  if (saveMultiPdf){
	outputfile = new TFile(outfilename.c_str(),"RECREATE");
	outputws = new RooWorkspace(); outputws->SetName("wtemplates");
  }

  system(Form("mkdir -p %s",outDir.c_str()));
  TFile *inFile = TFile::Open(fileName.c_str());
  RooWorkspace *inWS;
	if (saveMultiPdf){
		transferMacros(inFile,outputfile);

		RooRealVar *intL; 
		RooRealVar *sqrts;
		
		/*std::cout << "DEBUG LC inNT entries  " << inNT->GetEntries() << std::endl;
		//inNT->Print();
		float mgg;
		inNT->SetBranchAddress("mgg",&mgg);
		for ( int i =0 ; i < inNT->GetEntries() ; i++){
		inNT->GetEntry(i);
		
		} */
	}
	vector<string> functionClasses;
//	functionClasses.push_back("DijetSimple");
    functionClasses.push_back("Dijet");
	functionClasses.push_back("Exponential");
//	functionClasses.push_back("Expow");
// 	functionClasses.push_back("PowerLaw");
	functionClasses.push_back("Laurent");
	functionClasses.push_back("Atlas");
	map<string,string> namingMap;
//	namingMap.insert(pair<string,string>("DijetSimple","dijetsimp"));
	namingMap.insert(pair<string,string>("Dijet","dijet"));
	namingMap.insert(pair<string,string>("Atlas","atlas"));
	namingMap.insert(pair<string,string>("Exponential","exp"));
//	namingMap.insert(pair<string,string>("Expow","expow"));
//   	namingMap.insert(pair<string,string>("PowerLaw","pow"));
	namingMap.insert(pair<string,string>("Laurent","lau"));
	// store results here

	FILE *resFile ;
	resFile = fopen(Form("%s/fTestResults.txt",outDir.c_str()),"w");
	vector<map<string,int> > choices_vec;
	vector<map<string,std::vector<int> > > choices_envelope_vec;
	vector<map<string,RooAbsPdf*> > pdfs_vec;

	PdfModelBuilder pdfsModel;
	
	double upperEnvThreshold = 0.1; // upper threshold on delta(chi2) to include function in envelope (looser than truth function)

	fprintf(resFile,"Truth Model & d.o.f & $\\Delta NLL_{N+1}$ & $p(\\chi^{2}>\\chi^{2}_{(N\\rightarrow N+1)})$ \\\\\n");
	fprintf(resFile,"\\hline\n");

	for (int cat=startingCategory; cat<ncats; cat++){
	string catname;
	catname = Form("%s",diphotonCats_[cat].c_str());
    RooRealVar *mass = new RooRealVar (Form("mgg%s",catname.c_str()),Form("mgg%s",catname.c_str()), 200, 1600);
//	RooRealVar *mass = (RooRealVar*)inWS->var(Form("mgg%s",catname.c_str()));
	pdfsModel.setObsVar(mass);

  TNtuple *inNT;
		if (isData_){
			inNT = (TNtuple*)inFile->Get(Form("tree_data_%s_%s",analysisType_.c_str(),diphotonCats_[cat].c_str()));
		} else {
			inNT = (TNtuple*)inFile->Get(Form("tree_data_%s_%s",analysisType_.c_str(),diphotonCats_[cat].c_str()));
		}
	if(diphotonCats_[cat]=="EBEE") mass->setRange(320,1600); //FIXME Need a more configurable method to set range
	if(diphotonCats_[cat]=="EBEB") mass->setRange(230,1600); //FIXME Need a more configurable method to set range
	mass->setBins(nBinsForMass);
	cout << "-------------------------------------------------------------------------------" << endl;
	cout << "right mass loaded? " << endl;
	mass->Print("v");

	if (verbose) std::cout << "[INFO]  considering nTuple " << inNT->GetName() << std::endl;
		map<string,int> choices;
		map<string,std::vector<int> > choices_envelope;
		map<string,RooAbsPdf*> pdfs;
		map<string,RooAbsPdf*> allPdfs;
		RooDataSet *dataFull = new RooDataSet( Form("data_%s",catname.c_str()),"data_unbinned" ,inNT, RooArgSet(*mass) );

 		if (dataFull){
		std::cout << "[INFO] opened data for  "  << dataFull->GetName()  <<" - " << dataFull <<std::endl;
		if (verbose) dataFull->Print("V");
		}
	 	
		//TCanvas *t = new TCanvas();
		//RooPlot *frame = mass->frame();
		//dataFull->plotOn(frame);


		RooDataSet *data;
		string thisdataBinned_name;
		thisdataBinned_name =Form("data_binned_%s",diphotonCats_[cat].c_str());
		thisdataBinned_name =Form("data_binned_%s",diphotonCats_[cat].c_str());
		RooDataHist thisdataBinned(thisdataBinned_name.c_str(),"data",*mass,*dataFull);
		data = (RooDataSet*)&thisdataBinned;

		RooArgList storedPdfs("store");
		fprintf(resFile,"\\multicolumn{4}{|c|}{\\textbf{Category %d}} \\\\\n",cat);
		fprintf(resFile,"\\hline\n");

		double MinimimNLLSoFar=1e10;
		int simplebestFitPdfIndex = 0;

		// Standard F-Test to find the truth functions
		for (vector<string>::iterator funcType=functionClasses.begin(); 
				funcType!=functionClasses.end(); funcType++){

			double thisNll=0.; double prevNll=0.; double chi2=0.; double prob=0.; 
			int order=0; int prev_order=0; int cache_order=0;

			RooAbsPdf *prev_pdf=NULL;
			RooAbsPdf *cache_pdf=NULL;
			std::vector<int> pdforders;

			int counter =0;
			//	while (prob<0.05){
			while (prob<0.05 && order < 5){ //FIXME should be around order 3
				cout << "ORDER      " << order   << endl;
				RooAbsPdf *bkgPdf = getPdf(pdfsModel,*funcType,order,Form("ftest_pdf_%d_%s",cat,sqrts_.c_str()));
				if (!bkgPdf){
					// assume this order is not allowed
					order++;
				}
				else {

					//RooFitResult *fitRes = bkgPdf->fitTo(*data,Save(true),RooFit::Minimizer("Minuit2","minimize"));
					int fitStatus = 0;
					//thisNll = fitRes->minNll();
					runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/3);//bkgPdf->fitTo(*data,Save(true),RooFit::Minimizer("Minuit2","minimize"));
					if (fitStatus!=0) std::cout << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;
					chi2 = 2.*(prevNll-thisNll);
					if (chi2<0. && order>1) chi2=0.;
					if (prev_pdf!=NULL){
						prob = getProbabilityFtest(chi2,order-prev_order,prev_pdf,bkgPdf,mass,data
								,Form("%s/Ftest_from_%s%d_cat%d.pdf",outDir.c_str(),funcType->c_str(),order,cat));
						std::cout << "[INFO]  F-test Prob(chi2>chi2(data)) == " << prob << std::endl;
					} else {
						prob = 0;
					}
					double gofProb=0;
					// otherwise we get it later ...
					if (!saveMultiPdf) plot(mass,bkgPdf,data,Form("%s/%s%d_cat%d.png",outDir.c_str(),funcType->c_str(),order,cat),diphotonCats_,fitStatus,&gofProb);
					cout << "[INFO]\t " << *funcType << " " << order << " " << prevNll << " " << thisNll << " " << chi2 << " " << prob << endl;
					//fprintf(resFile,"%15s && %d && %10.2f && %10.2f && %10.2f \\\\\n",funcType->c_str(),order,thisNll,chi2,prob);
					prevNll=thisNll;
					cache_order=prev_order;
					cache_pdf=prev_pdf;
					prev_order=order;
					prev_pdf=bkgPdf;
					order++;
				}
				counter++;
			}

			fprintf(resFile,"%15s & %d & %5.2f & %5.2f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
			choices.insert(pair<string,int>(*funcType,cache_order));
			pdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),cache_order),cache_pdf));

			int truthOrder = cache_order;

			// Now run loop to determine functions inside envelope
			if (saveMultiPdf){
				chi2=0.;
				thisNll=0.;
				prevNll=0.;
				prob=0.;
				order=1;
				prev_order=0;
				cache_order=0;
				std::cout << "[INFO] Determining Envelope Functions for Family " << *funcType << ", cat " << cat << std::endl;
				std::cout << "[INFO] Upper end Threshold for highest order function " << upperEnvThreshold <<std::endl;

				while (prob<upperEnvThreshold){
					RooAbsPdf *bkgPdf = getPdf(pdfsModel,*funcType,order,Form("env_pdf_%d_%s",cat,sqrts_.c_str()));
					if (!bkgPdf ){
						// assume this order is not allowed
						if (order >6) { std::cout << " [WARNING] could not add ] " << std::endl; break ;}
						order++;
					}

					else {
						//RooFitResult *fitRes;
						int fitStatus=0;
						runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/3);//bkgPdf->fitTo(*data,Save(true),RooFit::Minimizer("Minuit2","minimize"));
						//thisNll = fitRes->minNll();
						if (fitStatus!=0) std::cout << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;
						double myNll = 2.*thisNll;
						chi2 = 2.*(prevNll-thisNll);
						if (chi2<0. && order>1) chi2=0.;
						prob = TMath::Prob(chi2,order-prev_order);

						cout << "[INFO] \t " << *funcType << " " << order << " " << prevNll << " " << thisNll << " " << chi2 << " " << prob << endl;
						prevNll=thisNll;
						cache_order=prev_order;
						cache_pdf=prev_pdf;

						// Calculate goodness of fit for the thing to be included (will use toys for lowstats)!
						double gofProb =0; 
						plot(mass,bkgPdf,data,Form("%s/%s%d_cat%d.pdf",outDir.c_str(),funcType->c_str(),order,cat),diphotonCats_,fitStatus,&gofProb);

						if ((prob < upperEnvThreshold) ) { // Looser requirements for the envelope

							if (gofProb > 0.01 || order == truthOrder ) {  // Good looking fit or one of our regular truth functions

								std::cout << "[INFO] Adding to Envelope " << bkgPdf->GetName() << " "<< gofProb 
									<< " 2xNLL + c is " << myNll + bkgPdf->getVariables()->getSize() <<  std::endl;
								allPdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),order),bkgPdf));
								storedPdfs.add(*bkgPdf);
								pdforders.push_back(order);

								// Keep track but we shall redo this later
								if ((myNll + bkgPdf->getVariables()->getSize()) < MinimimNLLSoFar) {
									simplebestFitPdfIndex = storedPdfs.getSize()-1;
									MinimimNLLSoFar = myNll + bkgPdf->getVariables()->getSize();
								}
							}
						}

						prev_order=order;
						prev_pdf=bkgPdf;
						order++;
					}
				}

				fprintf(resFile,"%15s & %d & %5.2f & %5.2f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
				choices_envelope.insert(pair<string,std::vector<int> >(*funcType,pdforders));
			}
		}

		fprintf(resFile,"\\hline\n");
		choices_vec.push_back(choices);
		choices_envelope_vec.push_back(choices_envelope);
		pdfs_vec.push_back(pdfs);

		plot(mass,pdfs,data,Form("%s/truths_cat%d",outDir.c_str(),cat),diphotonCats_,cat);

		if (saveMultiPdf){

			// Put selectedModels into a MultiPdf
			string catindexname;
			string catname;
				catindexname = Form("pdfindex_%s",diphotonCats_[cat].c_str());
				catname = Form("%s",diphotonCats_[cat].c_str());
			RooCategory catIndex(catindexname.c_str(),"c");
			RooMultiPdf *pdf = new RooMultiPdf(Form("model_bkg_%s",catname.c_str()),"all pdfs",catIndex,storedPdfs);
			RooRealVar nBackground(Form("model_bkg_%s_norm",catname.c_str()),"nbkg",data->sumEntries(),0,10E8);
			//nBackground.removeRange(); // bug in roofit will break combine until dev branch brought in
			//double check the best pdf!
			int bestFitPdfIndex = getBestFitFunction(pdf,data,&catIndex,!verbose);
			catIndex.setIndex(bestFitPdfIndex);
			std::cout << "// ------------------------------------------------------------------------- //" <<std::endl; 
			std::cout << "[INFO] Created MultiPdf " << pdf->GetName() << ", in Category " << cat << " with a total of " << catIndex.numTypes() << " pdfs"<< std::endl;
			storedPdfs.Print();
			std::cout << "[INFO] Best Fit Pdf = " << bestFitPdfIndex << ", " << storedPdfs.at(bestFitPdfIndex)->GetName() << std::endl;
			std::cout << "// ------------------------------------------------------------------------- //" <<std::endl;
			std::cout << "[INFO] Simple check of index "<< simplebestFitPdfIndex <<std::endl;

			mass->setBins(nBinsForMass);
			RooDataHist dataBinned(Form("data_binned_%s",catname.c_str()),"data",*mass,*dataFull);

			// Save it (also a binned version of the dataset
			outputws->import(*pdf);
			outputws->import(nBackground);
			outputws->import(catIndex);
			outputws->import(dataBinned);
			outputws->import(*data);
			outputws->import(*dataFull);
			plot(mass,pdf,&catIndex,data,Form("%s/multipdf_%s",outDir.c_str(),catname.c_str()),diphotonCats_,cat,bestFitPdfIndex);

		}

		}
		if (saveMultiPdf){
			outputfile->cd();
			outputws->Write();
			outputfile->Close();	
		}

		FILE *dfile = fopen(datfile.c_str(),"w");
		cout << "[RESULT] Recommended options" << endl;

		for (int cat=startingCategory; cat<ncats; cat++){
			cout << "Cat: " << diphotonCats_[cat] << endl;
			fprintf(dfile,"cat=%s\n",diphotonCats_[cat].c_str()); 
			for (map<string,int>::iterator it=choices_vec[cat-startingCategory].begin(); it!=choices_vec[cat-startingCategory].end(); it++){
				cout << "\t" << it->first << " - " << it->second << endl;
				fprintf(dfile,"truth=%s:%d:%s%d\n",it->first.c_str(),it->second,namingMap[it->first].c_str(),it->second);
			}
			for (map<string,std::vector<int> >::iterator it=choices_envelope_vec[cat-startingCategory].begin(); it!=choices_envelope_vec[cat-startingCategory].end(); it++){
				std::vector<int> ords = it->second;
				for (std::vector<int>::iterator ordit=ords.begin(); ordit!=ords.end(); ordit++){
					fprintf(dfile,"paul=%s:%d:%s%d\n",it->first.c_str(),*ordit,namingMap[it->first].c_str(),*ordit);
				}
			}
			fprintf(dfile,"\n");
		}
		inFile->Close();

		return 0;
	}
