  #include <TFile.h>
#include <TLine.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TLegend.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TFrame.h>
#include <TLatex.h>
//#include "../macro/mkPlotsLivia/CMS_lumi.C"
#include <iostream>
#include <vector>

double makeInterpolation(TGraph* g, bool isFirst);
void makePlots2HDM(TString, TString);
void makePlotsBaryo(TString, TString);
void makePlotsScalar(TString, TString);
void getLimits(TFile*, Double_t &, Double_t); 
void getXsec(TFile*, int, int, Double_t &); 
void getXsecBaryo(TFile*, int, int, Double_t &); 
void getXsecScalar(TFile*, int, int, Double_t &); 
void getEff(TFile*, int, int, Double_t &, Double_t &); 
void make2Dlimitplots_76X();

void make2Dlimitplots_76X(){
 
  cout << "Making 2D limit plots" << endl;

  //TString inDir	= "~soffi/public/4Margaret/2Dinputs/";
  TString inDir		= "ntuples4fit_pho_met0_met130_cic_default_shapes_lumi_35.9/";
  TString outDir	= "~/www/plotsMonoH/FitLimits/2HDM-2017";
  // SPECIFY LUMI in mkPlotsLivia/CMS_lumi.C

  makePlotsScalar(inDir, outDir);
  makePlotsBaryo(inDir, outDir);
  makePlots2HDM(inDir, outDir);

}



double makeInterpolation(TGraph* g, bool isFirst){

  double*  xvals=g->GetX();
  double*  yvals=g->GetY();

  double    deltaylow=99999999;
  double   deltayup=99999999;
  double yobs=1.;
  
  double yup=999.;
  double xup=999.;
  double ylow=999.;
  double xlow=999.;
  if(isFirst){
  for(int i=0;i<2; i++){
    //std::cout<<yvals[i]<<" "<<xvals[i]<<std::endl;
    double deltay=yvals[i]-yobs;
    if(deltay>=0 && deltay<deltayup){
      deltayup=deltay;
      yup=yvals[i];
      xup=xvals[i];
    }else if(deltay<0 && abs(deltay)<abs(deltaylow)){
      deltaylow=deltay;
      xlow=xvals[i];
      ylow=yvals[i];
      
    }
  }
  }else{
    for(int i=2;i<g->GetN(); i++){
    //std::cout<<yvals[i]<<" "<<xvals[i]<<std::endl;
    double deltay=yvals[i]-yobs;
    if(deltay>=0 && deltay<deltayup){
      deltayup=deltay;
      yup=yvals[i];
      xup=xvals[i];
    }else if(deltay<0 && abs(deltay)<abs(deltaylow)){
      deltaylow=deltay;
      xlow=xvals[i];
      ylow=yvals[i];
      
    }
  }

  }
  //  std::cout<<deltayup<<" "<<deltaylow<<" "<<xlow<<" "<<xup<<" "<<ylow<<" "<<yup<<std::endl;
  double m=(yup-ylow)/(xup-xlow);
  //  std::cout<<" m: "<< m << std::endl;
  double    xobs=(1-ylow+m*xlow)/m;
  
  return xobs;

}


void makePlotsScalar(TString inDir, TString outDir){

  // lumi
  Float_t lumi = 2.3;

  // mZp masses 
  Double_t mass[8] = {10,50,100,200,300,500,1000,10000};
  Int_t nMasses = 8;

  // pick up the higgsCombine files for each A0 mass
  std::vector<TFile* > higgsCombineFiles_Mchi1;
  std::vector<TFile* > higgsCombineFiles_Mchi10;
  std::vector<TFile* > higgsCombineFiles_Mchi50;
  std::vector<TFile* > higgsCombineFiles_Mchi150;
  std::vector<TFile* > higgsCombineFiles_Mchi500;
  std::vector<TFile* > higgsCombineFiles_Mchi1000;

  higgsCombineFiles_Mchi1.resize(nMasses);
  higgsCombineFiles_Mchi10.resize(nMasses);
  higgsCombineFiles_Mchi50.resize(nMasses);
  higgsCombineFiles_Mchi150.resize(nMasses);
  higgsCombineFiles_Mchi500.resize(nMasses);
  higgsCombineFiles_Mchi1000.resize(nMasses);

  for (int n=0; n<nMasses; n++){                 
    higgsCombineFiles_Mchi1[n] = new TFile(Form("%s/higgsCombineMonoHgg_sig_ScalarZp_mZP%d_mChi1.Asymptotic.mH1.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
    higgsCombineFiles_Mchi10[n] = new TFile(Form("%s/higgsCombineMonoHgg_sig_ScalarZp_mZP%d_mChi10.Asymptotic.mH10.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
    higgsCombineFiles_Mchi50[n] = new TFile(Form("%s/higgsCombineMonoHgg_sig_ScalarZp_mZP%d_mChi50.Asymptotic.mH50.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
    higgsCombineFiles_Mchi150[n] = new TFile(Form("%s/higgsCombineMonoHgg_sig_ScalarZp_mZP%d_mChi150.Asymptotic.mH150.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
    higgsCombineFiles_Mchi500[n] = new TFile(Form("%s/higgsCombineMonoHgg_sig_ScalarZp_mZP%d_mChi500.Asymptotic.mH500.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
    higgsCombineFiles_Mchi1000[n] = new TFile(Form("%s/higgsCombineMonoHgg_sig_ScalarZp_mZP%d_mChi1000.Asymptotic.mH1000.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
  }

 // pick up theory xsec
 TFile* theory = new TFile("theory_Scalar.root");

 // pick up efficiencies 
 // std::cout<<"debug1"<<std::endl;

 // make TLatex label
 TString latexCMSname = "CMS";
 TLatex *l1 = new TLatex(0.12,0.99,latexCMSname);
 l1->SetTextSize(0.036);
 l1->SetTextAlign(12);
 l1->SetNDC(kTRUE);
 l1->SetTextFont(62);
 TLatex *l1b = new TLatex(0.11,0.99,latexCMSname);
 l1b->SetTextSize(0.040);
 l1b->SetTextAlign(12);
 l1b->SetNDC(kTRUE);
 l1b->SetTextFont(62);

 TString latexlumi = Form("%1.1f fb^{-1}",lumi);
 TString latexenergy = " (13 TeV)";
 TString latexname = latexlumi+latexenergy;
 TLatex *l2 = new TLatex(0.65,0.99,latexname);
 l2->SetTextSize(0.034);
 l2->SetTextAlign(12);
 l2->SetNDC(kTRUE);
 l2->SetTextFont(42);
 //TLatex *l2b = new TLatex(0.74,0.90,latexname);
 //TLatex *l2b = new TLatex(0.725,0.95,latexname);
 TLatex *l2b = new TLatex(0.68,0.95,latexname);
 l2b->SetTextSize(0.034);
 l2b->SetTextAlign(12);
 l2b->SetNDC(kTRUE);
 l2b->SetTextFont(42);

 TString thestring = "Z'#rightarrow DM+h(#gamma#gamma)";
 TString add2hdm   = "(Baryonic)";
 //latex.SetTextSize(0.036);
 //latex.SetTextAlign(12); // centered
 //const char *thestring = "Z'#rightarrow DM+H(#gamma#gamma)";
 TLatex *l3 = new TLatex(0.11,0.82,thestring);
 l3->SetTextSize(0.036);
 l3->SetTextAlign(12);
 l3->SetNDC(kTRUE);
 l3->SetTextFont(42);
 TLatex *l3b = new TLatex(0.11,0.86,thestring);
 l3b->SetTextSize(0.031);
 l3b->SetTextAlign(12);
 l3b->SetNDC(kTRUE);
 l3b->SetTextFont(42);
 TLatex *l4 = new TLatex(0.11,0.85,add2hdm);
 l4->SetTextSize(0.036);
 l4->SetTextAlign(12);
 l4->SetNDC(kTRUE);
 l4->SetTextFont(42);
 //TLatex *l4b = new TLatex(0.115,0.73,add2hdm);
 TLatex *l4b = new TLatex(0.11,0.83,add2hdm);
 l4b->SetTextSize(0.031);
 l4b->SetTextAlign(12);
 l4b->SetNDC(kTRUE);
 l4b->SetTextFont(42);


 // setup 1D plots - expected
 TGraph* limit1;
 TGraph* limit10;
 TGraph* limit50;
 TGraph* limit150;
 TGraph* limit500;
 TGraph* limit1000;

 // setup 1D plots - observed
 TGraph* limit1_obs;
 TGraph* limit10_obs;
 TGraph* limit50_obs;
 TGraph* limit150_obs;
 TGraph* limit500_obs;
 TGraph* limit1000_obs;

 // setup output plot
 
 TH2D * limits = new TH2D("limits","limits",8,0,8,6,0,6);
 limits->GetXaxis()->SetTitle("m_{Z'} [GeV]");
 limits->GetYaxis()->SetTitle("m_{#chi} [GeV]");
 limits->GetZaxis()->SetTitle("#sigma_{95\% CL} / #sigma_{th}");
 limits->SetTitle("");
 limits->GetZaxis()->SetRangeUser(0.01,1e04);
 limits->SetMarkerSize(1.6);
 //limits->GetZaxis()->SetLimits(0.,10000);
 gStyle->SetHistMinimumZero(kFALSE);

 // limits->GetYaxis()->SetTitleOffset(1.3);
 //limits->GetZaxis()->SetTitleOffset(0.75);
 // size of axis labels
 limits->GetXaxis()->SetTitleSize(0.04);
 limits->GetYaxis()->SetTitleSize(0.04);
 limits->GetZaxis()->SetTitleSize(0.035);
 limits->GetXaxis()->SetLabelSize(0.05);
 limits->GetYaxis()->SetLabelSize(0.05); 
 limits->GetZaxis()->SetLabelSize(0.025);

 // set the lables for the Xaxis (mZp)
 limits->GetXaxis()->SetBinLabel(1,"10");
 limits->GetXaxis()->SetBinLabel(2,"50");
 limits->GetXaxis()->SetBinLabel(3,"100");
 limits->GetXaxis()->SetBinLabel(4,"200");
 limits->GetXaxis()->SetBinLabel(5,"300");
 limits->GetXaxis()->SetBinLabel(6,"500");
 limits->GetXaxis()->SetBinLabel(7,"1000");
 limits->GetXaxis()->SetBinLabel(8,"10000");

 // set the lables for the Yaxis (mA0)
 limits->GetYaxis()->SetBinLabel(1,"1");
 limits->GetYaxis()->SetBinLabel(2,"10");
 limits->GetYaxis()->SetBinLabel(3,"50");
 limits->GetYaxis()->SetBinLabel(4,"150");
 limits->GetYaxis()->SetBinLabel(5,"500");
 limits->GetYaxis()->SetBinLabel(6,"1000");

 // setup output observed plot
 TH2D * obslimits = (TH2D*) limits->Clone();

 // setup canvas
 //TCanvas * c = new TCanvas("c","",889,768);
 TCanvas * c = new TCanvas("c","",1);
 c->cd();
 gStyle->SetOptStat(0);
 //gStyle->SetPaintTextFormat("2.1f");
 // c->SetLeftMargin(0.1);
 c->SetRightMargin(0.11);

 Double_t limitval1[nMasses];
 Double_t limitval10[nMasses];
 Double_t limitval50[nMasses];
 Double_t limitval150[nMasses];
 Double_t limitval500[nMasses];
 Double_t limitval1000[nMasses];

 Double_t limitval1_obs[nMasses];
 Double_t limitval10_obs[nMasses];
 Double_t limitval50_obs[nMasses];
 Double_t limitval150_obs[nMasses];
 Double_t limitval500_obs[nMasses];
 Double_t limitval1000_obs[nMasses];

 Double_t xsec1[nMasses];
 Double_t xsec10[nMasses];
 Double_t xsec50[nMasses];
 Double_t xsec150[nMasses];
 Double_t xsec500[nMasses];
 Double_t xsec1000[nMasses]; 

 Double_t explimit1[nMasses];
 Double_t explimit10[nMasses];
 Double_t explimit50[nMasses];
 Double_t explimit150[nMasses];
 Double_t explimit500[nMasses];
 Double_t explimit1000[nMasses];
 
 Double_t obslimit1[nMasses];
 Double_t obslimit10[nMasses];
 Double_t obslimit50[nMasses];
 Double_t obslimit150[nMasses];
 Double_t obslimit500[nMasses];
 Double_t obslimit1000[nMasses];
 // std::cout<<"debug1"<<std::endl;
 for (Int_t n=0; n<nMasses; n++){
   getLimits(higgsCombineFiles_Mchi1[n],limitval1[n],0.5); 
   getLimits(higgsCombineFiles_Mchi10[n],limitval10[n],0.5); 
   getLimits(higgsCombineFiles_Mchi50[n],limitval50[n],0.5); 
   getLimits(higgsCombineFiles_Mchi150[n],limitval150[n],0.5); 
   getLimits(higgsCombineFiles_Mchi500[n],limitval500[n],0.5); 
   getLimits(higgsCombineFiles_Mchi1000[n],limitval1000[n],0.5); 
   //std::cout<<"---------------------------------------------------------"<<std::endl;
   //   std::cout<<mass[n]<<" "<<limitval1[n]<< " "<<limitval10[n]<<" "<<limitval50[n]<< " "<<limitval150[n] <<" "<<limitval500[n]<< " "<<limitval1000[n]<<std::endl;
   //std::cout<<"---------------------------------------------------------"<<std::endl;
   getLimits(higgsCombineFiles_Mchi1[n],limitval1_obs[n],-1.); 
   getLimits(higgsCombineFiles_Mchi10[n],limitval10_obs[n],-1.); 
   getLimits(higgsCombineFiles_Mchi50[n],limitval50_obs[n],-1.); 
   getLimits(higgsCombineFiles_Mchi150[n],limitval150_obs[n],-1.); 
   getLimits(higgsCombineFiles_Mchi500[n],limitval500_obs[n],-1.); 
   getLimits(higgsCombineFiles_Mchi1000[n],limitval1000_obs[n],-1.); 

   //   std::cout<<"debug1"<<std::endl;
   getXsecScalar(theory,1,(Int_t)mass[n],xsec1[n]);
   getXsecScalar(theory,10,(Int_t)mass[n],xsec10[n]);
   getXsecScalar(theory,50,(Int_t)mass[n],xsec50[n]);
   getXsecScalar(theory,150,(Int_t)mass[n],xsec150[n]);
   getXsecScalar(theory,500,(Int_t)mass[n],xsec500[n]);
   getXsecScalar(theory,1000,(Int_t)mass[n],xsec1000[n]);
   //   std::cout<<"debug1"<<std::endl;
   explimit1[n] = limitval1[n]/xsec1[n];
   explimit10[n] = limitval10[n]/xsec10[n];
   explimit50[n] = limitval50[n]/xsec50[n];
   explimit150[n] = limitval150[n]/xsec150[n];
   explimit500[n] = limitval500[n]/xsec500[n];
   explimit1000[n] = limitval1000[n]/xsec1000[n];

   obslimit1[n] = limitval1_obs[n]/xsec1[n];
   obslimit10[n] = limitval10_obs[n]/xsec10[n];
   obslimit50[n] = limitval50_obs[n]/xsec50[n];
   obslimit150[n] = limitval150_obs[n]/xsec150[n];
   obslimit500[n] = limitval500_obs[n]/xsec500[n];
   obslimit1000[n] = limitval1000_obs[n]/xsec1000[n];
   //   std::cout<<obslimit10[n]<<" ,------------------------------------- ,  "<<n<<std::endl;
     
   // fill limit plot
   limits->Fill(((Double_t)n+0.5),0.5,limitval1[n]/xsec1[n]);
   //   limits->Fill(((Double_t)n+0.5),1.5,1968); // --hardcode in some numbers
   limits->Fill(((Double_t)n+0.5),1.5,limitval10[n]/xsec10[n]);
   limits->Fill(((Double_t)n+0.5),2.5,limitval50[n]/xsec50[n]);
   //   limits->Fill(((Double_t)n+0.5),3.5,199.4);
   //limits->Fill(((Double_t)n+0.5),3.5,78.8);
   limits->Fill(((Double_t)n+0.5),3.5,limitval150[n]/xsec150[n]);
   limits->Fill(((Double_t)n+0.5),4.5,limitval500[n]/xsec500[n]);
   limits->Fill(((Double_t)n+0.5),5.5,limitval1000[n]/xsec1000[n]);

   obslimits->Fill(((Double_t)n+0.5),0.5,limitval1_obs[n]/xsec1[n]);
   //if(n==nMasses-1) obslimits->Fill(((Double_t)n+0.5),1.5,1941.1);
   obslimits->Fill(((Double_t)n+0.5),1.5,limitval10_obs[n]/xsec10[n]);
   obslimits->Fill(((Double_t)n+0.5),2.5,limitval50_obs[n]/xsec50[n]);
   //if (n==1)        obslimits->Fill(((Double_t)n+0.5),3.5,196.7);
   //else if (n==3)   obslimits->Fill(((Double_t)n+0.5),3.5,77.6);
   obslimits->Fill(((Double_t)n+0.5),3.5,limitval150_obs[n]/xsec150[n]);
   obslimits->Fill(((Double_t)n+0.5),4.5,limitval500_obs[n]/xsec500[n]);
   obslimits->Fill(((Double_t)n+0.5),5.5,limitval1000_obs[n]/xsec1000[n]);


 }


 // limits->Print("V ALL");
 // only pick up the limits that are non-zero
 Double_t mass_1[7] = {10,50,200,300,500,1000,10000};
 Double_t mass_10[4] = {10,50,100,10000};
 Double_t mass_50[3] = {200,300,10000};
 Double_t mass_150[4] = {200,500,1000,10000};
 Double_t mass_500[3] = {10,500,10000};
 Double_t mass_1000[2] = {10,10000};

 Double_t limitval_exp_1[7] = {explimit1[0],explimit1[1],explimit1[3],explimit1[4],explimit1[5],explimit1[6],explimit10[7]};
 Double_t limitval_obs_1[7] = {obslimit1[0],obslimit1[1],obslimit1[3],obslimit1[4],obslimit1[5],obslimit1[6],obslimit10[7]};

 Double_t limitval_exp_10[4] = {explimit10[0],explimit10[1],explimit10[2],explimit10[7]};
 Double_t limitval_obs_10[4] = {obslimit10[0],obslimit10[1],obslimit10[2],obslimit10[7]};
 

 Double_t limitval_exp_50[3] = {explimit50[3],explimit50[4],explimit50[7]};
 Double_t limitval_obs_50[3] = {obslimit50[3],obslimit50[4],obslimit50[7]};


 Double_t limitval_exp_150[4] = {explimit150[3],explimit150[5],explimit150[6],explimit150[7]};
 Double_t limitval_obs_150[4] = {obslimit150[3],obslimit150[5],obslimit150[6],obslimit150[7]};



 Double_t limitval_exp_500[3] = {explimit500[0],explimit500[5],explimit500[7]};
 Double_t limitval_obs_500[3] = {obslimit500[0],obslimit500[5],obslimit500[7]};


 Double_t limitval_exp_1000[2] = {explimit1000[0],explimit1000[7]};
 Double_t limitval_obs_1000[2] = {obslimit1000[0],obslimit1000[7]};

 limit1 = new TGraph(7,mass_1,limitval_exp_1);
 limit10 = new TGraph(4,mass_10,limitval_exp_10);
 limit50 = new TGraph(3,mass_50,limitval_exp_50);
 limit150 = new TGraph(4,mass_150,limitval_exp_150);
 limit500 = new TGraph(3,mass_500,limitval_exp_500);
 limit1000 = new TGraph(2,mass_1000,limitval_exp_1000);
 
 limit1_obs = new TGraph(7,mass_1,limitval_obs_1);
 limit10_obs = new TGraph(4,mass_10,limitval_obs_10);
 limit50_obs = new TGraph(3,mass_50,limitval_obs_50);
 limit150_obs = new TGraph(4,mass_150,limitval_obs_150);
 limit500_obs = new TGraph(3,mass_500,limitval_obs_500);
 limit1000_obs = new TGraph(2,mass_1000,limitval_obs_1000);
 
 //styling
 limit1->GetXaxis()->SetTitle("m_{Z'} [GeV]");
 limit1->GetYaxis()->SetTitle("#sigma_{95\% CL} / #sigma_{th}");
 limit1->SetTitle("");
 limit1->GetYaxis()->SetRangeUser(0.001,100000000000);
 limit1_obs->GetYaxis()->SetRangeUser(0.001,100000000000);
 //limit1->SetMinimum(0.9);
 limit1->SetLineWidth(2);
 limit10->SetLineWidth(2);
 limit50->SetLineWidth(2);
 limit150->SetLineWidth(2);
 limit500->SetLineWidth(2);
 limit1000->SetLineWidth(2);
 limit1->SetMarkerStyle(8);
 limit10->SetMarkerStyle(8);
 limit50->SetMarkerStyle(8);
 limit150->SetMarkerStyle(8);
 limit500->SetMarkerStyle(8);
 limit1000->SetMarkerStyle(8);
 // set up colors to match Hbb
 limit1->SetLineColor(kBlack);
 limit10->SetLineColor(kCyan);
 limit50->SetLineColor(kGreen);
 limit150->SetLineColor(kBlue);
 limit500->SetLineColor(kYellow);
 limit1000->SetLineColor(kMagenta);
 limit1->SetMarkerColor(kBlack);
 limit10->SetMarkerColor(kCyan);
 limit50->SetMarkerColor(kGreen);
 limit150->SetMarkerColor(kBlue);
 limit500->SetMarkerColor(kYellow);
 limit1000->SetMarkerColor(kMagenta);


 //styling
 limit1_obs->GetXaxis()->SetTitle("m_{Z'} [GeV]");
 limit1_obs->GetYaxis()->SetTitle("#sigma_{95\% CL} / #sigma_{th}");
 limit1_obs->SetTitle("");
 //limit1_obs->SetMaximum(10);
 //limit1_obs->SetMinimum(0.9);
 limit1_obs->SetLineWidth(2);
 limit10_obs->SetLineWidth(2);
 limit50_obs->SetLineWidth(2);
 limit150_obs->SetLineWidth(2);
 limit500_obs->SetLineWidth(2);
 limit1000_obs->SetLineWidth(2);
 limit1_obs->SetMarkerStyle(8);
 limit10_obs->SetMarkerStyle(8);
 limit50_obs->SetMarkerStyle(8);
 limit150_obs->SetMarkerStyle(8);
 limit500_obs->SetMarkerStyle(8);
 limit1000_obs->SetMarkerStyle(8);
 // set up colors to match Hbb
 limit1_obs->SetMarkerColor(kBlack);
 limit10_obs->SetMarkerColor(kCyan);
 limit50_obs->SetMarkerColor(kGreen);
 limit150_obs->SetMarkerColor(kBlue);
 limit500_obs->SetMarkerColor(kYellow);
 limit1000_obs->SetMarkerColor(kMagenta);
 limit1_obs->SetLineColor(kBlack);
 limit10_obs->SetLineColor(kCyan);
 limit50_obs->SetLineColor(kGreen);
 limit150_obs->SetLineColor(kBlue);
 limit500_obs->SetLineColor(kYellow);
 limit1000_obs->SetLineColor(kMagenta);
 


 // draw expected limits plot
 gStyle->SetPaintTextFormat("2.2e");
 gStyle->SetPalette(57);
 limits->Draw("COLZ TEXT"); 
 //latex.DrawLatex(0.08,5.7,thestring);
 // save plot
 //CMS_lumi(c,false,0);
 //l1->Draw("same");
 //l2->Draw("same");
 // l3->Draw("same");
 //l4->Draw("same");
 // c->cd();
  gStyle->SetPaintTextFormat("2.2e");
  c->SaveAs(Form("%s/limits2D_Scalar_exp.png",outDir.Data()));
  c->SaveAs(Form("%s/limits2D_Scalar_exp.pdf",outDir.Data()));

 // draw observed limits plot
  gStyle->SetPaintTextFormat("2.2e");
 gStyle->SetPalette(57);
 obslimits->Draw("COLZ TEXT"); 
 gStyle->SetPaintTextFormat("2.2e");
 //latex.DrawLatex(0.08,5.7,thestring);
 // save plot
 //CMS_lumi(c,false,0);
 //l1->Draw("same");
 //l2->Draw("same");
 // l3->Draw("same");
 //l4->Draw("same");
 c->cd();
 c->SaveAs(Form("%s/limits2D_Scalar_obs.png",outDir.Data()));
 c->SaveAs(Form("%s/limits2D_Scalar_obs.pdf",outDir.Data()));

 TLegend* leg = new TLegend(0.25,0.6,0.45,0.9,NULL,"brNDC"); // (x1,y1,x2,y2)
 leg->SetTextSize(0.046);
 leg->SetBorderSize(0);
 leg->SetLineColor(1);
 leg->SetLineWidth(1);
 leg->SetLineStyle(1);
 leg->SetFillColor(0);
 leg->SetFillStyle(0);
 leg->SetTextFont(62);
 leg->AddEntry(limit1,"m_{#chi} = 1 GeV","pl");
 leg->AddEntry(limit10,"m_{#chi} = 10 GeV","pl");
 leg->AddEntry(limit50,"m_{#chi} = 50 GeV","pl");
 leg->AddEntry(limit150,"m_{#chi} = 150 GeV","pl");
 leg->AddEntry(limit500,"m_{#chi} = 500 GeV","pl");
 leg->AddEntry(limit1000,"m_{#chi} = 1000 GeV","pl");
 leg->SetTextSize(0.03);

 TLine* line1 = new TLine();
 line1->SetX1(limit1->GetXaxis()->GetXmin());
 line1->SetY1(1.0);
 line1->SetX2(limit1->GetXaxis()->GetXmax());
 line1->SetY2(1.0);
 line1->SetLineColor(kRed);
 line1->SetLineWidth(2);

 // draw 1D comparisons --expected
 //c->Clear();
 //c->SetLogy(1);
 // c->SetLogx(1);
 // limit1->GetYaxis()->SetRangeUser(0.1, 1000);
 limit1->GetYaxis()->SetRangeUser(0.01, 100000000000);
 limit1->Draw("APL");

 limit10->Draw("PL SAME");
 limit50->Draw("PL SAME");
 limit150->Draw("PL SAME");
 limit500->Draw("PL SAME");
 limit1000->Draw("PL SAME");


 leg->Draw("SAME");
 line1->Draw("SAME");
 //latex.DrawLatex(50,2100,thestring);
 //CMS_lumi(c,false,0);
 //l1->Draw("same");
 //l2->Draw("same");
 // l3->Draw("same");
 //l4->Draw("same");
 // c->cd();
 c->SaveAs(Form("%s/limits_comparison_Scalar_exp.png",outDir.Data()));
 c->SaveAs(Form("%s/limits_comparison_Scalar_exp.pdf",outDir.Data()));

 // draw 1D comparisons --observed
 c->Clear();
 // c->SetLogy(1);
 // c->SetLogx(1);
 limit1_obs->GetYaxis()->SetRangeUser(0.01, 100000000000);
 // limit1_obs->GetYaxis()->SetRangeUser(0.1, 1000);
 limit1_obs->Draw("APL");
 //limit1_obs->GetYaxis()->SetRangeUser(0.1, 1000);
 limit10_obs->Draw("PL SAME");
 limit50_obs->Draw("PL SAME");
 limit150_obs->Draw("PL SAME");
 limit500_obs->Draw("PL SAME");
 limit1000_obs->Draw("PL SAME");
 leg->Draw("SAME");
 line1->Draw("SAME");
 //latex.DrawLatex(50,2100,thestring);
 //CMS_lumi(c,false,0);
 //l1->Draw("same");
 //l2->Draw("same");
 //l3->Draw("same");
 // l4->Draw("same");
 // c->cd();
 c->SaveAs(Form("%s/limits_comparison_Scalar_obs.png",outDir.Data()));
 c->SaveAs(Form("%s/limits_comparison_Scalar_obs.pdf",outDir.Data()));

 // draw 1D comparisons --expected & observed
 c->Clear();
 // c->SetLogy(1);
 limit1->GetYaxis()->SetRangeUser(0.001,10);
 limit1_obs->GetYaxis()->SetRangeUser(0.001,10);

 limit1->SetLineStyle(9);
 limit10->SetLineStyle(9);
 limit50->SetLineStyle(9);
 limit150->SetLineStyle(9);
 limit500->SetLineStyle(9);
 limit1000->SetLineStyle(9);
 limit1->Draw("APL");
 limit10->Draw("PL SAME");
 limit50->Draw("PL SAME");
 limit150->Draw("PL SAME");
 limit500->Draw("PL SAME");
 limit1000->Draw("PL SAME");
 limit1_obs->Draw("PL SAME");
 limit10_obs->Draw("PL SAME");
 limit50_obs->Draw("PL SAME");
 limit150_obs->Draw("PL SAME");
 limit500_obs->Draw("PL SAME");
 limit1000_obs->Draw("PL SAME");
 //leg->AddEntry(limit1,"Expected, m_{A0} = 1 GeV","pl");
 //leg->AddEntry(limit1_obs,"Observed, m_{A0} = 1 GeV","pl");
 //leg->AddEntry(limit10,"Expected, m_{A0} = 10 GeV","pl");
 //leg->AddEntry(limit10_obs,"Observed, m_{A0} = 1 GeV","pl");
 //leg->AddEntry(limit50,"Expected, m_{A0} = 50 GeV","pl");
 //leg->AddEntry(limit50_obs,"Observed, m_{A0} = 1 GeV","pl");
 //leg->AddEntry(limit150,"Expected, m_{A0} = 150 GeV","pl");
 //leg->AddEntry(limit150_obs,"Observed, m_{A0} = 1 GeV","pl");
 //leg->AddEntry(limit500,"Expected, m_{A0} = 500 GeV","pl");
 //leg->AddEntry(limit500_obs,"Observed, m_{A0} = 1 GeV","pl");
 //leg->AddEntry(limit1000,"Expected, m_{A0} = 1000 GeV","pl");
 //leg->AddEntry(limit1000_obs,"Observed, m_{A0} = 1 GeV","pl");
 leg->Draw("SAME");
 line1->Draw("SAME");
 //latex.DrawLatex(50,2100,thestring);
 //CMS_lumi(c,false,0);
 // l1->Draw("same");
 //l2->Draw("same");
 //l3->Draw("same");
 //l4->Draw("same");
 c->cd();
 c->SaveAs(Form("%s/limits_comparison_Scalar.png",outDir.Data()));
 c->SaveAs(Form("%s/limits_comparison_Scalar.pdf",outDir.Data()));


 delete c;


 // make plot with both expected and observed on same graph
 //TCanvas* cboth = new TCanvas("cboth","",889,768);
 TCanvas * cboth = new TCanvas("cboth","",1);
 cboth->cd();
 gStyle->SetOptStat(0);
 gStyle->SetPaintTextFormat("2.2e");
 gStyle->SetPalette(57);
 gStyle->SetFrameLineWidth(3);
 gStyle->SetPadRightMargin(0.109);
 gStyle->SetPadLeftMargin(0.13);
 cboth->SetLeftMargin(0.1);
 cboth->SetRightMargin(0.1);

 gStyle->SetHistMinimumZero(kFALSE);
 Double_t pad1_x1 = 0.01; 
 Double_t pad1_x2 = 0.98;
 Double_t pad1_y1 = 0.03;
 Double_t pad1_y2 = 1.00;
 Double_t pad2_x1 = pad1_x1;
 Double_t pad2_x2 = pad1_x2;
 Double_t pad2_y1 = pad1_y1-0.03;
 Double_t pad2_y2 = pad1_y2-0.03;


 //TPad* p1 = new TPad("p1","",0,0.12,0.95,0.98);
 //TPad* p1 = new TPad("p1","",0,0.09,0.95,0.89); //x1,y1,x2,y2
 //TPad* p1 = new TPad("p1","",0.02,0.03,0.98,0.95); //x1,y1,x2,y2
 TPad* p1 = new TPad("p1","",pad1_x1,pad1_y1,pad1_x2,pad1_y2); //x1,y1,x2,y2
 p1->Draw();
 p1->cd();
 p1->SetLogz();

 gStyle->SetPaintTextFormat("2.2e");
 obslimits->Draw("COLZ"); 
 // obslimits->Draw("TEXT SAME");

 float obscontent; 
 double obsbinX;
 double obsbinY;
 TString obsbincon;
 TLatex *obsbintxt; 
 for (int obsbinx = 1; obsbinx <= obslimits->GetXaxis()->GetNbins(); obsbinx++){
    for (int obsbiny = 1; obsbiny <= obslimits->GetYaxis()->GetNbins(); obsbiny++){
       obscontent   = obslimits->GetBinContent(obsbinx,obsbiny);
       //       if (obsbinx==1 && obsbiny>=3) continue;
       //if (obsbinx==2 && obsbiny>=5) continue;
       if (obscontent > 1000) obsbincon = TString::Format("%2.2e",obscontent);
       else obsbincon = TString::Format("%2.2e",obscontent);
       obsbinX      = obslimits->GetXaxis()->GetBinCenter(obsbinx);
       obsbinY      = obslimits->GetYaxis()->GetBinCenter(obsbiny);
       obsbintxt = new TLatex(obsbinX,obsbinY-0.08,obsbincon);
       obsbintxt->SetTextAlign(21);
       obsbintxt->SetTextSize(0.022);
       if(obscontent>0)obsbintxt->Draw("same");
    }
 }
 //latex.DrawLatex(0.08,5.7,thestring);
 p1->Update();

 Double_t x1,y1,x2,y2;
 p1->GetRange(x1,y1,x2,y2);

 cboth->cd();
 TPad* p2 = new TPad("p2","",pad2_x1,pad2_y1,pad2_x2,pad2_y2); //x1,y1,x2,y2
 p2->SetFillStyle(0);
 p2->SetFillColor(0);
 p2->Draw();
 p2->cd();
 p2->Range(x1,y1,x2,y2);
 gStyle->SetFrameLineWidth(3);
 TFrame *f = (TFrame*)cboth->FindObject("TFrame");
 Double_t px1 = f->GetX1();
 Double_t px2 = f->GetX2();
 Double_t py1 = f->GetY1()+0.23;
 Double_t py2 = f->GetY2()+0.23;

 limits->SetMarkerSize(1.4);
 limits->GetXaxis()->SetTitle("");
 limits->GetYaxis()->SetTitle("");
 limits->SetTitle("");
 limits->GetXaxis()->SetBinLabel(1,"");
 limits->GetXaxis()->SetBinLabel(2,"");
 limits->GetXaxis()->SetBinLabel(3,"");
 limits->GetXaxis()->SetBinLabel(4,"");
 limits->GetXaxis()->SetBinLabel(5,"");
 limits->GetXaxis()->SetBinLabel(6,"");
 limits->GetXaxis()->SetBinLabel(7,"");
 limits->GetXaxis()->SetBinLabel(8,"");
 limits->GetYaxis()->SetBinLabel(1,"");
 limits->GetYaxis()->SetBinLabel(2,"");
 limits->GetYaxis()->SetBinLabel(3,"");
 limits->GetYaxis()->SetBinLabel(4,"");
 limits->GetYaxis()->SetBinLabel(5,"");
 limits->GetYaxis()->SetBinLabel(6,"");

 float content; 
 double binX;
 double binY;
 TString bincon;
 TLatex *bintxt; 
 for (int binx = 1; binx <= limits->GetXaxis()->GetNbins(); binx++){
    for (int biny = 1; biny <= limits->GetYaxis()->GetNbins(); biny++){
       content   = limits->GetBinContent(binx,biny);
       //       if (binx==1 && biny>=3) continue;
       //if (binx==2 && biny>=5) continue;
       if (binx==limits->GetXaxis()->GetNbins()) bincon = TString::Format("(%2.2f)",content);
       else bincon = TString::Format("(%2.2e)",content);
       binX      = limits->GetXaxis()->GetBinCenter(binx);
       binY      = limits->GetYaxis()->GetBinCenter(biny);
       bintxt = new TLatex(binX,binY-0.08,bincon);
       bintxt->SetTextAlign(21);
       bintxt->SetTextSize(0.022);
       if(content>0)bintxt->Draw("same");
    }
 }

 //limits->Draw("TEXT SAME"); 
 p1->Update();

 // redraw the frame around the histogram
 TLine l;
 l.SetLineWidth(3);
 l.DrawLine(px1,py2,px2,py2);
 l.DrawLine(px2,py1,px2,py2);
 l.DrawLine(px1,py1,px2,py1);
 l.DrawLine(px1,py1,px1,py2);

 //CMS_lumi(cboth,false,0);
 // l1b->Draw("same");
 //l2b->Draw("same");
 //l3b->Draw("same");
 //l4b->Draw("same");
 cboth->cd();
 cboth->SaveAs(Form("%s/limits2D_Scalar_ExpAndObs.png",outDir.Data()));
 cboth->SaveAs(Form("%s/limits2D_Scalar_ExpAndObs.pdf",outDir.Data()));
 delete cboth;


}









void makePlotsBaryo(TString inDir, TString outDir){

  // lumi
  Float_t lumi = 2.3;

  // mZp masses 
  Double_t mass[8] = {10,50,100,200,500,1000,2000,10000};
  Int_t nMasses = 8;

  // pick up the higgsCombine files for each A0 mass
  std::vector<TFile* > higgsCombineFiles_Mchi1;
  std::vector<TFile* > higgsCombineFiles_Mchi10;
  std::vector<TFile* > higgsCombineFiles_Mchi50;
  std::vector<TFile* > higgsCombineFiles_Mchi150;
  std::vector<TFile* > higgsCombineFiles_Mchi500;
  std::vector<TFile* > higgsCombineFiles_Mchi1000;

  higgsCombineFiles_Mchi1.resize(nMasses);
  higgsCombineFiles_Mchi10.resize(nMasses);
  higgsCombineFiles_Mchi50.resize(nMasses);
  higgsCombineFiles_Mchi150.resize(nMasses);
  higgsCombineFiles_Mchi500.resize(nMasses);
  higgsCombineFiles_Mchi1000.resize(nMasses);

  for (int n=0; n<nMasses; n++){  
    higgsCombineFiles_Mchi1[n] = new TFile(Form("%s/higgsCombineMonoHgg_sig_ZpBaryonic_mZP%d_mChi1.Asymptotic.mH1.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
    higgsCombineFiles_Mchi10[n] = new TFile(Form("%s/higgsCombineMonoHgg_sig_ZpBaryonic_mZP%d_mChi10.Asymptotic.mH10.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
    higgsCombineFiles_Mchi50[n] = new TFile(Form("%s/higgsCombineMonoHgg_sig_ZpBaryonic_mZP%d_mChi50.Asymptotic.mH50.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
    higgsCombineFiles_Mchi150[n] = new TFile(Form("%s/higgsCombineMonoHgg_sig_ZpBaryonic_mZP%d_mChi150.Asymptotic.mH150.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
    higgsCombineFiles_Mchi500[n] = new TFile(Form("%s/higgsCombineMonoHgg_sig_ZpBaryonic_mZP%d_mChi500.Asymptotic.mH500.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
    higgsCombineFiles_Mchi1000[n] = new TFile(Form("%s/higgsCombineMonoHgg_sig_ZpBaryonic_mZP%d_mChi1000.Asymptotic.mH1000.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
  }

 // pick up theory xsec
 TFile* theory = new TFile("theory_ZPBaryonic.root");

 // pick up efficiencies 
 // std::cout<<"debug1"<<std::endl;

 // make TLatex label
 TString latexCMSname = "CMS";
 TLatex *l1 = new TLatex(0.12,0.99,latexCMSname);
 l1->SetTextSize(0.036);
 l1->SetTextAlign(12);
 l1->SetNDC(kTRUE);
 l1->SetTextFont(62);
 TLatex *l1b = new TLatex(0.11,0.99,latexCMSname);
 l1b->SetTextSize(0.040);
 l1b->SetTextAlign(12);
 l1b->SetNDC(kTRUE);
 l1b->SetTextFont(62);

 TString latexlumi = Form("%1.1f fb^{-1}",lumi);
 TString latexenergy = " (13 TeV)";
 TString latexname = latexlumi+latexenergy;
 TLatex *l2 = new TLatex(0.65,0.99,latexname);
 l2->SetTextSize(0.034);
 l2->SetTextAlign(12);
 l2->SetNDC(kTRUE);
 l2->SetTextFont(42);
 //TLatex *l2b = new TLatex(0.74,0.90,latexname);
 //TLatex *l2b = new TLatex(0.725,0.95,latexname);
 TLatex *l2b = new TLatex(0.68,0.95,latexname);
 l2b->SetTextSize(0.034);
 l2b->SetTextAlign(12);
 l2b->SetNDC(kTRUE);
 l2b->SetTextFont(42);

 TString thestring = "Z'#rightarrow DM+h(#gamma#gamma)";
 TString add2hdm   = "(Baryonic)";
 //latex.SetTextSize(0.036);
 //latex.SetTextAlign(12); // centered
 //const char *thestring = "Z'#rightarrow DM+H(#gamma#gamma)";
 TLatex *l3 = new TLatex(0.11,0.82,thestring);
 l3->SetTextSize(0.036);
 l3->SetTextAlign(12);
 l3->SetNDC(kTRUE);
 l3->SetTextFont(42);
 TLatex *l3b = new TLatex(0.11,0.86,thestring);
 l3b->SetTextSize(0.031);
 l3b->SetTextAlign(12);
 l3b->SetNDC(kTRUE);
 l3b->SetTextFont(42);
 TLatex *l4 = new TLatex(0.11,0.85,add2hdm);
 l4->SetTextSize(0.036);
 l4->SetTextAlign(12);
 l4->SetNDC(kTRUE);
 l4->SetTextFont(42);
 //TLatex *l4b = new TLatex(0.115,0.73,add2hdm);
 TLatex *l4b = new TLatex(0.11,0.83,add2hdm);
 l4b->SetTextSize(0.031);
 l4b->SetTextAlign(12);
 l4b->SetNDC(kTRUE);
 l4b->SetTextFont(42);


 // setup 1D plots - expected
 TGraph* limit1;
 TGraph* limit10;
 TGraph* limit50;
 TGraph* limit150;
 TGraph* limit500;
 TGraph* limit1000;

 // setup 1D plots - observed
 TGraph* limit1_obs;
 TGraph* limit10_obs;
 TGraph* limit50_obs;
 TGraph* limit150_obs;
 TGraph* limit500_obs;
 TGraph* limit1000_obs;

 // setup output plot
 
 TH2D * limits = new TH2D("limits","limits",8,0,8,6,0,6);
 limits->GetXaxis()->SetTitle("m_{Z'} [GeV]");
 limits->GetYaxis()->SetTitle("m_{#chi} [GeV]");
 limits->GetZaxis()->SetTitle("#sigma_{95\% CL} / #sigma_{th}");
 limits->SetTitle("");
 limits->GetZaxis()->SetRangeUser(0.1,1e04);
 limits->SetMarkerSize(1.6);
 //limits->GetZaxis()->SetLimits(0.,10000);
 gStyle->SetHistMinimumZero(kFALSE);

 // limits->GetYaxis()->SetTitleOffset(1.3);
 //limits->GetZaxis()->SetTitleOffset(0.75);
 // size of axis labels
 limits->GetXaxis()->SetTitleSize(0.04);
 limits->GetYaxis()->SetTitleSize(0.04);
 limits->GetZaxis()->SetTitleSize(0.035);
 limits->GetXaxis()->SetLabelSize(0.05);
 limits->GetYaxis()->SetLabelSize(0.05); 
 limits->GetZaxis()->SetLabelSize(0.025);

 // set the lables for the Xaxis (mZp)
 limits->GetXaxis()->SetBinLabel(1,"10");
 limits->GetXaxis()->SetBinLabel(2,"50");
 limits->GetXaxis()->SetBinLabel(3,"100");
 limits->GetXaxis()->SetBinLabel(4,"200");
 limits->GetXaxis()->SetBinLabel(5,"500");
 limits->GetXaxis()->SetBinLabel(6,"1000");
 limits->GetXaxis()->SetBinLabel(7,"2000");
 limits->GetXaxis()->SetBinLabel(8,"10000");

 // set the lables for the Yaxis (mA0)
 limits->GetYaxis()->SetBinLabel(1,"1");
 limits->GetYaxis()->SetBinLabel(2,"10");
 limits->GetYaxis()->SetBinLabel(3,"50");
 limits->GetYaxis()->SetBinLabel(4,"150");
 limits->GetYaxis()->SetBinLabel(5,"500");
 limits->GetYaxis()->SetBinLabel(6,"1000");

 // setup output observed plot
 TH2D * obslimits = (TH2D*) limits->Clone();

 // setup canvas
 //TCanvas * c = new TCanvas("c","",889,768);
 TCanvas * c = new TCanvas("c","",1);
 c->cd();
 gStyle->SetOptStat(0);
 //gStyle->SetPaintTextFormat("2.1f");
 // c->SetLeftMargin(0.1);
 c->SetRightMargin(0.11);

 Double_t limitval1[nMasses];
 Double_t limitval10[nMasses];
 Double_t limitval50[nMasses];
 Double_t limitval150[nMasses];
 Double_t limitval500[nMasses];
 Double_t limitval1000[nMasses];

 Double_t limitval1_obs[nMasses];
 Double_t limitval10_obs[nMasses];
 Double_t limitval50_obs[nMasses];
 Double_t limitval150_obs[nMasses];
 Double_t limitval500_obs[nMasses];
 Double_t limitval1000_obs[nMasses];

 Double_t xsec1[nMasses];
 Double_t xsec10[nMasses];
 Double_t xsec50[nMasses];
 Double_t xsec150[nMasses];
 Double_t xsec500[nMasses];
 Double_t xsec1000[nMasses]; 

 Double_t explimit1[nMasses];
 Double_t explimit10[nMasses];
 Double_t explimit50[nMasses];
 Double_t explimit150[nMasses];
 Double_t explimit500[nMasses];
 Double_t explimit1000[nMasses];
 
 Double_t obslimit1[nMasses];
 Double_t obslimit10[nMasses];
 Double_t obslimit50[nMasses];
 Double_t obslimit150[nMasses];
 Double_t obslimit500[nMasses];
 Double_t obslimit1000[nMasses];
 // std::cout<<"debug1"<<std::endl;
 for (Int_t n=0; n<nMasses; n++){
   getLimits(higgsCombineFiles_Mchi1[n],limitval1[n],0.5); 
   getLimits(higgsCombineFiles_Mchi10[n],limitval10[n],0.5); 
   getLimits(higgsCombineFiles_Mchi50[n],limitval50[n],0.5); 
   getLimits(higgsCombineFiles_Mchi150[n],limitval150[n],0.5); 
   getLimits(higgsCombineFiles_Mchi500[n],limitval500[n],0.5); 
   getLimits(higgsCombineFiles_Mchi1000[n],limitval1000[n],0.5); 
   //std::cout<<n<<" "<<limitval300[n]<< " "<<limitval400[n]<<limitval500[n]<< " "<<limitval600[n] <<limitval700[n]<< " "<<limitval800[n]<<std::endl;
   std::cout<<"---------------------Mchi 1---------------- "<<mass[n]<<std::endl;
   getLimits(higgsCombineFiles_Mchi1[n],limitval1_obs[n],-1.); 
   std::cout<<"---------------------Mchi 10---------------- "<<mass[n]<<std::endl;
   getLimits(higgsCombineFiles_Mchi10[n],limitval10_obs[n],-1.); 
   std::cout<<"---------------------Mchi 50---------------- "<<mass[n]<<std::endl;
   getLimits(higgsCombineFiles_Mchi50[n],limitval50_obs[n],-1.); 
   std::cout<<"---------------------Mchi 150---------------- "<<mass[n]<<std::endl;
   getLimits(higgsCombineFiles_Mchi150[n],limitval150_obs[n],-1.); 
   std::cout<<"---------------------Mchi 500---------------- "<<mass[n]<<std::endl;
   getLimits(higgsCombineFiles_Mchi500[n],limitval500_obs[n],-1.); 
   std::cout<<"---------------------Mchi 1000---------------- "<<mass[n]<<std::endl;
   getLimits(higgsCombineFiles_Mchi1000[n],limitval1000_obs[n],-1.); 
   std::cout<<"------------------------------------- "<<std::endl;
   //   std::cout<<"debug1"<<std::endl;
   getXsecBaryo(theory,1,(Int_t)mass[n],xsec1[n]);
   // std::cout<<"MAssCHI=1:  massZp: "<<mass[n]<< " xsec: "<<xsec1[n]<<std::endl;
   getXsecBaryo(theory,10,(Int_t)mass[n],xsec10[n]);
   getXsecBaryo(theory,50,(Int_t)mass[n],xsec50[n]);
   getXsecBaryo(theory,150,(Int_t)mass[n],xsec150[n]);
   getXsecBaryo(theory,500,(Int_t)mass[n],xsec500[n]);
   getXsecBaryo(theory,1000,(Int_t)mass[n],xsec1000[n]);
   // std::cout<<"debug1"<<std::endl;
   explimit1[n] = limitval1[n]/xsec1[n];
   explimit10[n] = limitval10[n]/xsec10[n];
   explimit50[n] = limitval50[n]/xsec50[n];
   explimit150[n] = limitval150[n]/xsec150[n];
   explimit500[n] = limitval500[n]/xsec500[n];
   explimit1000[n] = limitval1000[n]/xsec1000[n];

   obslimit1[n] = limitval1_obs[n]/xsec1[n];
   obslimit10[n] = limitval10_obs[n]/xsec10[n];
   obslimit50[n] = limitval50_obs[n]/xsec50[n];
   obslimit150[n] = limitval150_obs[n]/xsec150[n];
   obslimit500[n] = limitval500_obs[n]/xsec500[n];
   obslimit1000[n] = limitval1000_obs[n]/xsec1000[n];

   
   std::cout<<"Mass Chi =1: massZ: "<<mass[n]<<" explim: "<<limitval1[n]<< " obslim :"<<limitval1_obs[n]<< " theory: "<<xsec1[n]<<std::endl;
   //   std::cout<<obslimit10[n]<<" ,------------------------------------- ,  "<<n<<std::endl;
   if(limitval1[n]==0) limits->SetBinContent(((Double_t)n+0.5),0.5,0);
   if(limitval1[n]==0) limits->SetBinError(((Double_t)n+0.5),0.5,0);
   if(limitval10[n]==0) limits->SetBinContent(((Double_t)n+0.5),0.5,0);
   if(limitval10[n]==0) limits->SetBinError(((Double_t)n+0.5),0.5,0);
   if(limitval50[n]==0) limits->SetBinContent(((Double_t)n+0.5),0.5,0);
   if(limitval50[n]==0) limits->SetBinError(((Double_t)n+0.5),0.5,0);
   if(limitval150[n]==0) limits->SetBinContent(((Double_t)n+0.5),0.5,0);
   if(limitval150[n]==0) limits->SetBinError(((Double_t)n+0.5),0.5,0);
   if(limitval500[n]==0) limits->SetBinContent(((Double_t)n+0.5),0.5,0);
   if(limitval500[n]==0) limits->SetBinError(((Double_t)n+0.5),0.5,0);
   if(limitval1000[n]==0) limits->SetBinContent(((Double_t)n+0.5),0.5,0);
   if(limitval1000[n]==0) limits->SetBinError(((Double_t)n+0.5),0.5,0);
     
   // fill limit plot
   limits->Fill(((Double_t)n+0.5),0.5,limitval1[n]/xsec1[n]);
   //   limits->Fill(((Double_t)n+0.5),1.5,1968); // --hardcode in some numbers
   limits->Fill(((Double_t)n+0.5),1.5,limitval10[n]/xsec10[n]);
   limits->Fill(((Double_t)n+0.5),2.5,limitval50[n]/xsec50[n]);
   //   limits->Fill(((Double_t)n+0.5),3.5,199.4);
   //limits->Fill(((Double_t)n+0.5),3.5,78.8);
   limits->Fill(((Double_t)n+0.5),3.5,limitval150[n]/xsec150[n]);
   limits->Fill(((Double_t)n+0.5),4.5,limitval500[n]/xsec500[n]);
   limits->Fill(((Double_t)n+0.5),5.5,limitval1000[n]/xsec1000[n]);

   obslimits->Fill(((Double_t)n+0.5),0.5,limitval1_obs[n]/xsec1[n]);
   //if(n==nMasses-1) obslimits->Fill(((Double_t)n+0.5),1.5,1941.1);
   obslimits->Fill(((Double_t)n+0.5),1.5,limitval10_obs[n]/xsec10[n]);
   obslimits->Fill(((Double_t)n+0.5),2.5,limitval50_obs[n]/xsec50[n]);
   //if (n==1)        obslimits->Fill(((Double_t)n+0.5),3.5,196.7);
   //else if (n==3)   obslimits->Fill(((Double_t)n+0.5),3.5,77.6);
   obslimits->Fill(((Double_t)n+0.5),3.5,limitval150_obs[n]/xsec150[n]);
   obslimits->Fill(((Double_t)n+0.5),4.5,limitval500_obs[n]/xsec500[n]);
   obslimits->Fill(((Double_t)n+0.5),5.5,limitval1000_obs[n]/xsec1000[n]);


 }


 // limits->Print("V ALL");
 // only pick up the limits that are non-zero
 Double_t mass_10[4] = {10,50,100,10000};
 Double_t mass_50[3] = {10,200,10000};
 Double_t mass_150[4] = {200,500,1000,10000};
 Double_t mass_500[3] = {500,2000,10000};
 Double_t mass_1000[3] = {10,1000,10000};

 Double_t limitval_exp_10[4] = {explimit10[0],explimit10[1],explimit10[2],explimit10[7]};
 Double_t limitval_obs_10[4] = {obslimit10[0],obslimit10[1],obslimit10[2],obslimit10[7]};
 

 Double_t limitval_exp_50[3] = {explimit50[0],explimit50[3],explimit50[7]};
 Double_t limitval_obs_50[3] = {obslimit50[0],obslimit50[3],obslimit50[7]};


 Double_t limitval_exp_150[4] = {explimit150[3],explimit150[4],explimit150[5],explimit150[7]};
 Double_t limitval_obs_150[4] = {obslimit150[3],obslimit150[4],obslimit150[5],obslimit150[7]};



 Double_t limitval_exp_500[3] = {explimit500[4],explimit500[6],explimit500[7]};
 Double_t limitval_obs_500[3] = {obslimit500[4],obslimit500[6],obslimit500[7]};


 Double_t limitval_exp_1000[3] = {explimit1000[0],explimit1000[5],explimit1000[7]};
 Double_t limitval_obs_1000[3] = {obslimit1000[0],obslimit1000[5],obslimit1000[7]};

 limit1 = new TGraph(nMasses,mass,explimit1);
 limit10 = new TGraph(4,mass_10,limitval_exp_10);
 limit50 = new TGraph(3,mass_50,limitval_exp_50);
 limit150 = new TGraph(4,mass_150,limitval_exp_150);
 limit500 = new TGraph(3,mass_500,limitval_exp_500);
 limit1000 = new TGraph(3,mass_1000,limitval_exp_1000);
 
 limit1_obs = new TGraph(nMasses,mass,obslimit1);
 limit10_obs = new TGraph(4,mass_10,limitval_obs_10);
 limit50_obs = new TGraph(3,mass_50,limitval_obs_50);
 limit150_obs = new TGraph(4,mass_150,limitval_obs_150);
 limit500_obs = new TGraph(3,mass_500,limitval_obs_500);
 limit1000_obs = new TGraph(3,mass_1000,limitval_obs_1000);
 
 //styling
 limit1->GetXaxis()->SetTitle("m_{Z'} [GeV]");
 limit1->GetYaxis()->SetTitle("#sigma_{95\% CL} / #sigma_{th}");
 limit1->SetTitle("");
 limit1->GetYaxis()->SetRange(0.001,100000000);
 limit1_obs->GetYaxis()->SetRange(0.001,100000000);
 //limit1->SetMinimum(0.9);
 limit1->SetLineWidth(2);
 limit10->SetLineWidth(2);
 limit50->SetLineWidth(2);
 limit150->SetLineWidth(2);
 limit500->SetLineWidth(2);
 limit1000->SetLineWidth(2);
 limit1->SetMarkerStyle(8);
 limit10->SetMarkerStyle(8);
 limit50->SetMarkerStyle(8);
 limit150->SetMarkerStyle(8);
 limit500->SetMarkerStyle(8);
 limit1000->SetMarkerStyle(8);
 // set up colors to match Hbb
 limit1->SetLineColor(kBlack);
 limit10->SetLineColor(kCyan);
 limit50->SetLineColor(kGreen);
 limit150->SetLineColor(kBlue);
 limit500->SetLineColor(kYellow);
 limit1000->SetLineColor(kMagenta);
 limit1->SetMarkerColor(kBlack);
 limit10->SetMarkerColor(kCyan);
 limit50->SetMarkerColor(kGreen);
 limit150->SetMarkerColor(kBlue);
 limit500->SetMarkerColor(kYellow);
 limit1000->SetMarkerColor(kMagenta);


 //styling
 limit1_obs->GetXaxis()->SetTitle("m_{Z'} [GeV]");
 limit1_obs->GetYaxis()->SetTitle("#sigma_{95\% CL} / #sigma_{th}");
 limit1_obs->SetTitle("");
 //limit1_obs->SetMaximum(10);
 //limit1_obs->SetMinimum(0.9);
 limit1_obs->SetLineWidth(2);
 limit10_obs->SetLineWidth(2);
 limit50_obs->SetLineWidth(2);
 limit150_obs->SetLineWidth(2);
 limit500_obs->SetLineWidth(2);
 limit1000_obs->SetLineWidth(2);
 limit1_obs->SetMarkerStyle(8);
 limit10_obs->SetMarkerStyle(8);
 limit50_obs->SetMarkerStyle(8);
 limit150_obs->SetMarkerStyle(8);
 limit500_obs->SetMarkerStyle(8);
 limit1000_obs->SetMarkerStyle(8);
 // set up colors to match Hbb
 limit1_obs->SetMarkerColor(kBlack);
 limit10_obs->SetMarkerColor(kCyan);
 limit50_obs->SetMarkerColor(kGreen);
 limit150_obs->SetMarkerColor(kBlue);
 limit500_obs->SetMarkerColor(kYellow);
 limit1000_obs->SetMarkerColor(kMagenta);
 limit1_obs->SetLineColor(kBlack);
 limit10_obs->SetLineColor(kCyan);
 limit50_obs->SetLineColor(kGreen);
 limit150_obs->SetLineColor(kBlue);
 limit500_obs->SetLineColor(kYellow);
 limit1000_obs->SetLineColor(kMagenta);
 


 // draw expected limits plot
 gStyle->SetPaintTextFormat("2.1e");
 gStyle->SetPalette(57);
 limits->Draw("COLZ TEXT"); 
 //latex.DrawLatex(0.08,5.7,thestring);
 // save plot
 //CMS_lumi(c,false,0);
 //l1->Draw("same");
 //l2->Draw("same");
 // l3->Draw("same");
 //l4->Draw("same");
 // c->cd();
  gStyle->SetPaintTextFormat("2.1e");
  c->SaveAs(Form("%s/limits2D_Baryo_exp.png",outDir.Data()));
  c->SaveAs(Form("%s/limits2D_Baryo_exp.pdf",outDir.Data()));

 // draw observed limits plot
  gStyle->SetPaintTextFormat("2.1e");
 gStyle->SetPalette(57);
 obslimits->Draw("COLZ TEXT"); 
 gStyle->SetPaintTextFormat("2.1e");
 //latex.DrawLatex(0.08,5.7,thestring);
 // save plot
 //CMS_lumi(c,false,0);
 //l1->Draw("same");
 //l2->Draw("same");
 // l3->Draw("same");
 //l4->Draw("same");
 c->cd();
 c->SaveAs(Form("%s/limits2D_Baryo_obs.png",outDir.Data()));
 c->SaveAs(Form("%s/limits2D_Baryo_obs.pdf",outDir.Data()));

 TLegend* leg = new TLegend(0.25,0.6,0.45,0.9,NULL,"brNDC"); // (x1,y1,x2,y2)
 leg->SetTextSize(0.046);
 leg->SetBorderSize(0);
 leg->SetLineColor(1);
 leg->SetLineWidth(1);
 leg->SetLineStyle(1);
 leg->SetFillColor(0);
 leg->SetFillStyle(0);
 leg->SetTextFont(62);
 leg->AddEntry(limit1,"m_{#chi} = 1 GeV","pl");
 leg->AddEntry(limit10,"m_{#chi} = 10 GeV","pl");
 leg->AddEntry(limit50,"m_{#chi} = 50 GeV","pl");
 leg->AddEntry(limit150,"m_{#chi} = 150 GeV","pl");
 leg->AddEntry(limit500,"m_{#chi} = 500 GeV","pl");
 leg->AddEntry(limit1000,"m_{#chi} = 1000 GeV","pl");
 leg->SetTextSize(0.03);

 TLine* line1 = new TLine();
 line1->SetX1(limit1->GetXaxis()->GetXmin());
 line1->SetY1(1.0);
 line1->SetX2(1050);
 line1->SetY2(1.0);
 line1->SetLineColor(kRed);
 line1->SetLineWidth(2);

 // draw 1D comparisons --expected
 //c->Clear();
 // c->SetLogy(1);
 //c->SetLogx(1);
 limit1->GetYaxis()->SetRangeUser(0.01, 6);
 limit1->GetXaxis()->SetRangeUser(0.01, 1050);
 limit1->Draw("APL");
 // limit1->GetYaxis()->SetRangeUser(0.1, 1000);
 limit10->Draw("PL SAME");
 limit50->Draw("PL SAME");
 limit150->Draw("PL SAME");
 limit500->Draw("PL SAME");
 limit1000->Draw("PL SAME");


 leg->Draw("SAME");
 line1->Draw("SAME");
 //latex.DrawLatex(50,2100,thestring);
 //CMS_lumi(c,false,0);
 //l1->Draw("same");
 //l2->Draw("same");
 // l3->Draw("same");
 //l4->Draw("same");
 // c->cd();
 c->SaveAs(Form("%s/limits_comparison_Baryo_exp.png",outDir.Data()));
 c->SaveAs(Form("%s/limits_comparison_Baryo_exp.pdf",outDir.Data()));

 // draw 1D comparisons --observed
 c->Clear();

 // c->SetLogy(1);
 //c->SetLogx(1);
 limit1_obs->GetYaxis()->SetRangeUser(0.01, 6);
 limit1_obs->GetXaxis()->SetRangeUser(0.01, 1050);

 // limit1_obs->GetYaxis()->SetRangeUser(0.1, 1000);
 limit1_obs->Draw("APL");
 //limit1_obs->GetYaxis()->SetRangeUser(0.1, 1000);
 limit10_obs->Draw("PL SAME");
 limit50_obs->Draw("PL SAME");
 limit150_obs->Draw("PL SAME");
 limit500_obs->Draw("PL SAME");
 limit1000_obs->Draw("PL SAME");
 leg->Draw("SAME");
 line1->Draw("SAME");
 //latex.DrawLatex(50,2100,thestring);
 //CMS_lumi(c,false,0);
 //l1->Draw("same");
 //l2->Draw("same");
 //l3->Draw("same");
 // l4->Draw("same");
 // c->cd();
 c->SaveAs(Form("%s/limits_comparison_Baryo_obs.png",outDir.Data()));
 c->SaveAs(Form("%s/limits_comparison_Baryo_obs.pdf",outDir.Data()));

 // draw 1D comparisons --expected & observed
 c->Clear();
 // c->SetLogy(1);
 limit1->SetLineStyle(9);
 limit10->SetLineStyle(9);
 limit50->SetLineStyle(9);
 limit150->SetLineStyle(9);
 limit500->SetLineStyle(9);
 limit1000->SetLineStyle(9);
 //c->SetLogy(1);
 // c->SetLogx(1);
 limit1->Draw("APL");
 limit10->Draw("PL SAME");
 limit50->Draw("PL SAME");
 limit150->Draw("PL SAME");
 limit500->Draw("PL SAME");
 limit1000->Draw("PL SAME");
 limit1_obs->Draw("PL SAME");
 limit10_obs->Draw("PL SAME");
 limit50_obs->Draw("PL SAME");
 limit150_obs->Draw("PL SAME");
 limit500_obs->Draw("PL SAME");
 limit1000_obs->Draw("PL SAME");
 //leg->AddEntry(limit1,"Expected, m_{A0} = 1 GeV","pl");
 //leg->AddEntry(limit1_obs,"Observed, m_{A0} = 1 GeV","pl");
 //leg->AddEntry(limit10,"Expected, m_{A0} = 10 GeV","pl");
 //leg->AddEntry(limit10_obs,"Observed, m_{A0} = 1 GeV","pl");
 //leg->AddEntry(limit50,"Expected, m_{A0} = 50 GeV","pl");
 //leg->AddEntry(limit50_obs,"Observed, m_{A0} = 1 GeV","pl");
 //leg->AddEntry(limit150,"Expected, m_{A0} = 150 GeV","pl");
 //leg->AddEntry(limit150_obs,"Observed, m_{A0} = 1 GeV","pl");
 //leg->AddEntry(limit500,"Expected, m_{A0} = 500 GeV","pl");
 //leg->AddEntry(limit500_obs,"Observed, m_{A0} = 1 GeV","pl");
 //leg->AddEntry(limit1000,"Expected, m_{A0} = 1000 GeV","pl");
 //leg->AddEntry(limit1000_obs,"Observed, m_{A0} = 1 GeV","pl");
 leg->Draw("SAME");
 line1->Draw("SAME");
 //latex.DrawLatex(50,2100,thestring);
 //CMS_lumi(c,false,0);
 // l1->Draw("same");
 //l2->Draw("same");
 //l3->Draw("same");
 //l4->Draw("same");
 c->cd();

 c->SaveAs(Form("%s/limits_comparison_Baryo.png",outDir.Data()));
 c->SaveAs(Form("%s/limits_comparison_Baryo.pdf",outDir.Data()));


 delete c;


 // make plot with both expected and observed on same graph
 //TCanvas* cboth = new TCanvas("cboth","",889,768);
 TCanvas * cboth = new TCanvas("cboth","",1);
 cboth->cd();
 gStyle->SetOptStat(0);
 gStyle->SetPaintTextFormat("2.1e");
 gStyle->SetPalette(57);
 gStyle->SetFrameLineWidth(3);
 gStyle->SetPadRightMargin(0.109);
 gStyle->SetPadLeftMargin(0.13);
 cboth->SetLeftMargin(0.1);
 cboth->SetRightMargin(0.1);

 gStyle->SetHistMinimumZero(kFALSE);
 Double_t pad1_x1 = 0.01; 
 Double_t pad1_x2 = 0.98;
 Double_t pad1_y1 = 0.03;
 Double_t pad1_y2 = 1.00;
 Double_t pad2_x1 = pad1_x1;
 Double_t pad2_x2 = pad1_x2;
 Double_t pad2_y1 = pad1_y1-0.03;
 Double_t pad2_y2 = pad1_y2-0.03;


 //TPad* p1 = new TPad("p1","",0,0.12,0.95,0.98);
 //TPad* p1 = new TPad("p1","",0,0.09,0.95,0.89); //x1,y1,x2,y2
 //TPad* p1 = new TPad("p1","",0.02,0.03,0.98,0.95); //x1,y1,x2,y2
 TPad* p1 = new TPad("p1","",pad1_x1,pad1_y1,pad1_x2,pad1_y2); //x1,y1,x2,y2
 p1->Draw();
 p1->cd();
 p1->SetLogz();

 gStyle->SetPaintTextFormat("2.1e");
 obslimits->Draw("COLZ"); 
 // obslimits->Draw("TEXT SAME");

 float obscontent; 
 double obsbinX;
 double obsbinY;
 TString obsbincon;
 TLatex *obsbintxt; 
 for (int obsbinx = 1; obsbinx <= obslimits->GetXaxis()->GetNbins(); obsbinx++){
    for (int obsbiny = 1; obsbiny <= obslimits->GetYaxis()->GetNbins(); obsbiny++){
       obscontent   = obslimits->GetBinContent(obsbinx,obsbiny);
       //       if (obsbinx==1 && obsbiny>=3) continue;
       //if (obsbinx==2 && obsbiny>=5) continue;
       if (obscontent > 1000) obsbincon = TString::Format("%2.1e",obscontent);
       else obsbincon = TString::Format("%2.1e",obscontent);
       obsbinX      = obslimits->GetXaxis()->GetBinCenter(obsbinx);
       obsbinY      = obslimits->GetYaxis()->GetBinCenter(obsbiny);
       obsbintxt = new TLatex(obsbinX,obsbinY-0.08,obsbincon);
       obsbintxt->SetTextAlign(21);
       obsbintxt->SetTextSize(0.022);
       if(obscontent>0)obsbintxt->Draw("same");
    }
 }
 //latex.DrawLatex(0.08,5.7,thestring);
 p1->Update();

 Double_t x1,y1,x2,y2;
 p1->GetRange(x1,y1,x2,y2);

 cboth->cd();
 TPad* p2 = new TPad("p2","",pad2_x1,pad2_y1,pad2_x2,pad2_y2); //x1,y1,x2,y2
 p2->SetFillStyle(0);
 p2->SetFillColor(0);
 p2->Draw();
 p2->cd();
 p2->Range(x1,y1,x2,y2);
 gStyle->SetFrameLineWidth(3);
 TFrame *f = (TFrame*)cboth->FindObject("TFrame");
 Double_t px1 = f->GetX1();
 Double_t px2 = f->GetX2();
 Double_t py1 = f->GetY1()+0.23;
 Double_t py2 = f->GetY2()+0.23;

 limits->SetMarkerSize(1.4);
 limits->GetXaxis()->SetTitle("");
 limits->GetYaxis()->SetTitle("");
 limits->SetTitle("");
 limits->GetXaxis()->SetBinLabel(1,"");
 limits->GetXaxis()->SetBinLabel(2,"");
 limits->GetXaxis()->SetBinLabel(3,"");
 limits->GetXaxis()->SetBinLabel(4,"");
 limits->GetXaxis()->SetBinLabel(5,"");
 limits->GetXaxis()->SetBinLabel(6,"");
 limits->GetXaxis()->SetBinLabel(7,"");
 limits->GetXaxis()->SetBinLabel(8,"");
 limits->GetYaxis()->SetBinLabel(1,"");
 limits->GetYaxis()->SetBinLabel(2,"");
 limits->GetYaxis()->SetBinLabel(3,"");
 limits->GetYaxis()->SetBinLabel(4,"");
 limits->GetYaxis()->SetBinLabel(5,"");
 limits->GetYaxis()->SetBinLabel(6,"");

 float content; 
 double binX;
 double binY;
 TString bincon;
 TLatex *bintxt; 
 for (int binx = 1; binx <= limits->GetXaxis()->GetNbins(); binx++){
    for (int biny = 1; biny <= limits->GetYaxis()->GetNbins(); biny++){
       content   = limits->GetBinContent(binx,biny);
       //       if (binx==1 && biny>=3) continue;
       //if (binx==2 && biny>=5) continue;
       if (binx==limits->GetXaxis()->GetNbins()) bincon = TString::Format("(%2.0f)",content);
       else bincon = TString::Format("(%2.1e)",content);
       binX      = limits->GetXaxis()->GetBinCenter(binx);
       binY      = limits->GetYaxis()->GetBinCenter(biny);
       bintxt = new TLatex(binX,binY-0.08,bincon);
       bintxt->SetTextAlign(21);
       bintxt->SetTextSize(0.022);
       if(content>0)bintxt->Draw("same");
    }
 }

 //limits->Draw("TEXT SAME"); 
 p1->Update();

 // redraw the frame around the histogram
 TLine l;
 l.SetLineWidth(3);
 l.DrawLine(px1,py2,px2,py2);
 l.DrawLine(px2,py1,px2,py2);
 l.DrawLine(px1,py1,px2,py1);
 l.DrawLine(px1,py1,px1,py2);

 //CMS_lumi(cboth,false,0);
 // l1b->Draw("same");
 //l2b->Draw("same");
 //l3b->Draw("same");
 //l4b->Draw("same");
 cboth->cd();
 cboth->SaveAs(Form("%s/limits2D_Baryo_ExpAndObs.png",outDir.Data()));
 cboth->SaveAs(Form("%s/limits2D_Baryo_ExpAndObs.pdf",outDir.Data()));
 delete cboth;






 ///compute interpolating points
 double exp_mZP_mChi1=makeInterpolation(limit1, false);
 std::cout<<"Expected excluded MZP for MCHI=1: "<<exp_mZP_mChi1<<std::endl;
 double exp_mZP_mChi10_first=makeInterpolation(limit10,true);
 std::cout<<"Expected excluded MZP First for MCHI=10: "<<exp_mZP_mChi10_first<<std::endl;
 double exp_mZP_mChi10_2nd=makeInterpolation(limit10,false);
 std::cout<<"Expected excluded MZP First for MCHI=10: "<<exp_mZP_mChi10_2nd<<std::endl;
 double exp_mZP_mChi50_first=makeInterpolation(limit50,true);
 std::cout<<"Expected excluded MZP First for MCHI=50: "<<exp_mZP_mChi50_first<<std::endl;
 /// double exp_mZP_mChi50_2nd=makeInterpolation(limit50,false);
 // std::cout<<"Expected excluded MZP First for MCHI=50: "<<exp_mZP_mChi50_2nd<<std::endl;
 double exp_mZP_mChi150_first=makeInterpolation(limit150,true);
 std::cout<<"Expected excluded MZP First for MCHI=150: "<<exp_mZP_mChi150_first<<std::endl;
 double exp_mZP_mChi150_2nd=makeInterpolation(limit150,false);
 std::cout<<"Expected excluded MZP First for MCHI=150: "<<exp_mZP_mChi150_2nd<<std::endl;

 double mzp_excl_exp[7]={exp_mZP_mChi1,exp_mZP_mChi150_2nd,exp_mZP_mChi150_first,exp_mZP_mChi50_first,exp_mZP_mChi10_first,10,exp_mZP_mChi1};
 double mchi_excl_exp[7]={1,150,150,50,10,1,1};
 double obs_mZP_mChi1=makeInterpolation(limit1_obs, false);
 std::cout<<"Obsected excluded MZP for MCHI=1: "<<obs_mZP_mChi1<<std::endl;
 double obs_mZP_mChi10_first=makeInterpolation(limit10_obs,true);
 std::cout<<"Obsected excluded MZP First for MCHI=10: "<<obs_mZP_mChi10_first<<std::endl;
 double obs_mZP_mChi10_2nd=makeInterpolation(limit10_obs,false);
 std::cout<<"Obsected excluded MZP First for MCHI=10: "<<obs_mZP_mChi10_2nd<<std::endl;
 double obs_mZP_mChi50_first=makeInterpolation(limit50_obs,true);
 std::cout<<"Obsected excluded MZP First for MCHI=50: "<<obs_mZP_mChi50_first<<std::endl;
 /// double obs_mZP_mChi50_2nd=makeInterpolation(limit50,false);
 // std::cout<<"Obsected excluded MZP First for MCHI=50: "<<obs_mZP_mChi50_2nd<<std::endl;
 double obs_mZP_mChi150_first=makeInterpolation(limit150_obs,true);
 std::cout<<"Obsected excluded MZP First for MCHI=150: "<<obs_mZP_mChi150_first<<std::endl;
 double obs_mZP_mChi150_2nd=makeInterpolation(limit150_obs,false);
 std::cout<<"Obsected excluded MZP First for MCHI=150: "<<obs_mZP_mChi150_2nd<<std::endl;
 double  mzp_excl_obs[7]={obs_mZP_mChi1,obs_mZP_mChi150_2nd,obs_mZP_mChi150_first,obs_mZP_mChi50_first,obs_mZP_mChi10_first,10,obs_mZP_mChi1};
 double mchi_excl_obs[7]={1,150,150,50,10,1,1};


 double massZP[13]={10,50,100,200,500,1000,2000,50,100,200,500,1000,2000};
 double massChi[13]={ 1., 1., 1., 1., 1., 1.,1.,10.,10.,50.,150.,150.,500.};
 // double exp_vec[17]={explimit1[0],explimit1[1],explimit1[2],explimit1[3],explimit1[4],explimit1[5],explimit1[6],explimit1[7],explimit10[0],explimit10[1],explimit10[2],explimit10[7],explimit50[0],explimit50[3],explimit50[7],explimit150[3],explimit150[4],explimit150[5],explimit150[7],explimit500[4],explimit500[6],explimit500[7],explimit1000[0],explimit1000[5],explimit1000[7]};
 //double obs_vec[17]={obslimit1[0],obslimit1[1],obslimit1[2],obslimit1[3],obslimit1[4],obslimit1[5],obslimit1[6],obslimit1[7],obslimit10[0],obslimit10[1],obslimit10[2],obslimit10[7],obslimit50[0],obslimit50[3],obslimit50[7],obslimit150[3],obslimit150[4],obslimit150[5],obslimit150[7],obslimit500[4],obslimit500[6],obslimit500[7],obslimit1000[0],obslimit1000[5],obslimit1000[7]};
 double obs_vec[13]={obslimit1[0],obslimit1[1],obslimit1[2],obslimit1[3],obslimit1[4],obslimit1[5],obslimit1[6],obslimit10[1],obslimit10[2],obslimit50[3],obslimit150[4],obslimit150[5],obslimit500[6]};

 for(int i=0;i<13; i++) std::cout<<massZP[i]<<" "<<massChi[i]<<" "<<obs_vec[i]<<std::endl;
 TCanvas* cex = new TCanvas("cex","cex",400,400);
 cex->SetLeftMargin(0.1);
 cex->SetRightMargin(0.1);
 gStyle->SetFrameLineWidth(3);
 gStyle->SetPadRightMargin(0.109);
 gStyle->SetPadLeftMargin(0.13);
 cex->SetLeftMargin(0.1);
 cex->SetRightMargin(0.1);
 gStyle->SetHistMinimumZero(kFALSE);
 Double_t pad1x1 = 0.01;
 Double_t pad1x2 = 0.93;
 Double_t pad1y1 = 0.03;
 Double_t pad1y2 = 1.00;
 Double_t pad2x1 = pad1x1;
 Double_t pad2x2 = pad1x2;
 Double_t pad2y1 = pad1y1-0.03;
 Double_t pad2y2 = pad1y2-0.03;


 //TPad* p1 = new TPad("p1","",0,0.12,0.95,0.98);                                                                                 
 //TPad* p1 = new TPad("p1","",0,0.09,0.95,0.89); //x1,y1,x2,y2                                                                   
 //TPad* p1 = new TPad("p1","",0.02,0.03,0.98,0.95); //x1,y1,x2,y2                                                                
 TPad* pp1 = new TPad("pp1","",pad1x1,pad1y1,pad1x2,pad1y2); //x1,y1,x2,y2                                                      
 pp1->Draw();
 pp1->cd();
 pp1->SetLogz();


 TGraph* exp_g = new TGraph(7,mzp_excl_exp,mchi_excl_exp);
 TGraph* obs_g = new TGraph(7,mzp_excl_obs,mchi_excl_obs);

 // cex->cd();
 TGraph2D* mu = new TGraph2D(13,massZP,massChi,obs_vec);
 TH2F* h = new TH2F("h","h",200,10,2000,200,1,500);
 mu->SetHistogram(h);
 h->GetZaxis()->SetRangeUser(0.1,100);
 h->Draw("COLZ");
 h->GetZaxis()->SetLabelOffset(0.5*h->GetZaxis()->GetLabelOffset());
 gStyle->SetPalette(57,0);
 // mu->GetYaxis()->SetTitle("m_{#chi} [GeV]");
 //mu->GetXaxis()->SetTitle("m_{Z'} [GeV]");
 //mu->GetZaxis()->SetTitle("Expected #mu 95%CL");
 
 
 cex->SetLogz();                                                                                                                 // cex->SetLogx();                                                                                                                         
 // cex->SetLogy();                                                                                                                         
 mu->Print("V ALL");
 /* mu->GetZaxis()->SetRangeUser(0.001,1e03);
 mu->GetYaxis()->SetRangeUser(0.1,1e03);
 mu->GetYaxis()->SetRangeUser(1,3e03);
 */
 // mu->GetZaxis()->SetRangeUser(0.00001,1e03);
 // mu->GetXaxis()->SetRangeUser(0.01,2e03);
 //mu->GetYaxis()->SetRangeUser(0.01,500);

 mu->Draw("COLZsame");
 // mu->GetZaxis()->SetRangeUser(0.00001,1e03);
 // mu->GetXaxis()->SetRangeUser(0.01,2e03);
 //mu->GetYaxis()->SetRangeUser(0.01,500);
 mu->GetZaxis()->SetLabelOffset(0.5*mu->GetZaxis()->GetLabelOffset());
 exp_g->SetLineWidth(2);
 exp_g->SetLineStyle(kDashed);
 exp_g->SetLineColor(kRed);
 exp_g->Draw("LSAME");
 obs_g->SetLineWidth(2);
 obs_g->SetLineColor(kRed);
 obs_g->Draw("LSAME");
 TLegend* legex = new TLegend(0.2,0.7, 0.5, 0.9);
 legex->SetFillColor(0);
 legex->AddEntry(exp_g,"#splitline{Median}{Expected 95% CL}", "L");
 legex->AddEntry(obs_g,"#splitline{Median}{Observed 95% CL}", "L");
 TString latexCMSnameex = "#bf{CMS} #it{Preliminary}";
 cex->cd();
 TLatex *l1ex = new TLatex(0.13,0.97,latexCMSname);
 l1ex->SetTextSize(0.036);
 l1ex->SetTextAlign(12);
 l1ex->SetNDC(kTRUE);
 l1ex->SetTextFont(42);

 TString latexlumiex = Form("34.7 fb^{-1}");
 TString latexenergyex = " (13 TeV)";
 TString latexnameex = latexlumiex+latexenergyex;
 TLatex *l2ex = new TLatex(0.59,0.97,latexnameex);
 l2ex->SetTextSize(0.034);
 l2ex->SetTextAlign(12);
 l2ex->SetNDC(kTRUE);
 l2ex->SetTextFont(42);


 TString thestringex = "Z'_{B}#rightarrow DM+H(#gamma#gamma)";
 TLatex *l3ex = new TLatex(0.21,0.65,thestringex);
 l3ex->SetTextSize(0.036);
 l3ex->SetTextAlign(12);
 l3ex->SetNDC(kTRUE);
 l3ex->SetTextFont(42);
 l1ex->Draw("same");
 l2ex->Draw("same");


 TLatex *lxex = new TLatex(0.7,0.05,"m_{Z'_{B}} [GeV]");
 lxex->SetTextSize(0.054);
 lxex->SetTextAlign(12);
 lxex->SetNDC(kTRUE);
 lxex->SetTextFont(42);
 lxex->Draw("same");
 TLatex *lyex = new TLatex(0.03,0.77,"m_{#chi} [GeV]");
 lyex->SetTextSize(0.054);
 lyex->SetTextAlign(12);
 lyex->SetTextAngle(90);
 lyex->SetNDC(kTRUE);
 lyex->SetTextFont(42);
 lyex->Draw("same");
 TLatex *lzex = new TLatex(0.97,0.95,"Observed #mu 95% CL");
 lzex->SetTextSize(0.054);
 lzex->SetTextAlign(12);
 lzex->SetTextAngle(270);
 lzex->SetNDC(kTRUE);
 lzex->SetTextFont(42);
 lzex->Draw("same");
 legex->Draw("same");
 l3ex->Draw("same");
 cex->SaveAs("~/www/plotsMonoH/FitLimits/2HDM-2017/2D_ZPBaryonic.png");
 cex->SaveAs("~/www/plotsMonoH/FitLimits/2HDM-2017/2D_ZPBaryonic.pdf");

}





void makePlots2HDM(TString inDir, TString outDir){

  // lumi
  Float_t lumi = 2.3;

  // mZp masses 
  Double_t mass[8] = {600,800,1000,1200,1400,1700,2000,2500};
  Int_t nMasses = 8;

  // pick up the higgsCombine files for each A0 mass
  std::vector<TFile* > higgsCombineFiles_MA0300;
  std::vector<TFile* > higgsCombineFiles_MA0400;
  std::vector<TFile* > higgsCombineFiles_MA0500;
  std::vector<TFile* > higgsCombineFiles_MA0600;
  std::vector<TFile* > higgsCombineFiles_MA0700;
  std::vector<TFile* > higgsCombineFiles_MA0800;

  higgsCombineFiles_MA0300.resize(nMasses);
  higgsCombineFiles_MA0400.resize(nMasses);
  higgsCombineFiles_MA0500.resize(nMasses);
  higgsCombineFiles_MA0600.resize(nMasses);
  higgsCombineFiles_MA0700.resize(nMasses);
  higgsCombineFiles_MA0800.resize(nMasses);

  for (int n=0; n<nMasses; n++){
    higgsCombineFiles_MA0300[n] = new TFile(Form("%s/higgsCombineMonoHgg_mZP%d_mA0300.Asymptotic.mH%d.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
    higgsCombineFiles_MA0400[n] = new TFile(Form("%s/higgsCombineMonoHgg_mZP%d_mA0400.Asymptotic.mH%d.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
    higgsCombineFiles_MA0500[n] = new TFile(Form("%s/higgsCombineMonoHgg_mZP%d_mA0500.Asymptotic.mH%d.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
    higgsCombineFiles_MA0600[n] = new TFile(Form("%s/higgsCombineMonoHgg_mZP%d_mA0600.Asymptotic.mH%d.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
    higgsCombineFiles_MA0700[n] = new TFile(Form("%s/higgsCombineMonoHgg_mZP%d_mA0700.Asymptotic.mH%d.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
    higgsCombineFiles_MA0800[n] = new TFile(Form("%s/higgsCombineMonoHgg_mZP%d_mA0800.Asymptotic.mH%d.root",inDir.Data(),(Int_t)mass[n],(Int_t)mass[n]));
  }

 // pick up theory xsec
 TFile* theory_gz08 = new TFile("ScanPlot_gz08.root");

 // pick up efficiencies 


 // make TLatex label
 TString latexCMSname = "CMS";
 TLatex *l1 = new TLatex(0.12,0.99,latexCMSname);
 l1->SetTextSize(0.036);
 l1->SetTextAlign(12);
 l1->SetNDC(kTRUE);
 l1->SetTextFont(62);
 TLatex *l1b = new TLatex(0.11,0.99,latexCMSname);
 l1b->SetTextSize(0.040);
 l1b->SetTextAlign(12);
 l1b->SetNDC(kTRUE);
 l1b->SetTextFont(62);

 TString latexlumi = Form("%1.1f fb^{-1}",lumi);
 TString latexenergy = " (13 TeV)";
 TString latexname = latexlumi+latexenergy;
 TLatex *l2 = new TLatex(0.65,0.99,latexname);
 l2->SetTextSize(0.034);
 l2->SetTextAlign(12);
 l2->SetNDC(kTRUE);
 l2->SetTextFont(42);
 //TLatex *l2b = new TLatex(0.74,0.90,latexname);
 //TLatex *l2b = new TLatex(0.725,0.95,latexname);
 TLatex *l2b = new TLatex(0.68,0.99,latexname);
 l2b->SetTextSize(0.034);
 l2b->SetTextAlign(12);
 l2b->SetNDC(kTRUE);
 l2b->SetTextFont(42);

 TString thestring = "Z'#rightarrow DM+h(#gamma#gamma)";
 TString add2hdm   = "(2HDM)";
 //latex.SetTextSize(0.036);
 //latex.SetTextAlign(12); // centered
 //const char *thestring = "Z'#rightarrow DM+H(#gamma#gamma)";
 TLatex *l3 = new TLatex(0.11,0.82,thestring);
 l3->SetTextSize(0.036);
 l3->SetTextAlign(12);
 l3->SetNDC(kTRUE);
 l3->SetTextFont(42);
 TLatex *l3b = new TLatex(0.11,0.86,thestring);
 l3b->SetTextSize(0.031);
 l3b->SetTextAlign(12);
 l3b->SetNDC(kTRUE);
 l3b->SetTextFont(42);
 TLatex *l4 = new TLatex(0.11,0.85,add2hdm);
 l4->SetTextSize(0.036);
 l4->SetTextAlign(12);
 l4->SetNDC(kTRUE);
 l4->SetTextFont(42);
 //TLatex *l4b = new TLatex(0.115,0.73,add2hdm);
 TLatex *l4b = new TLatex(0.11,0.83,add2hdm);
 l4b->SetTextSize(0.031);
 l4b->SetTextAlign(12);
 l4b->SetNDC(kTRUE);
 l4b->SetTextFont(42);


 // setup 1D plots - expected
 TGraph* limit300;
 TGraph* limit400;
 TGraph* limit500;
 TGraph* limit600;
 TGraph* limit700;
 TGraph* limit800;

 // setup 1D plots - observed
 TGraph* limit300_obs;
 TGraph* limit400_obs;
 TGraph* limit500_obs;
 TGraph* limit600_obs;
 TGraph* limit700_obs;
 TGraph* limit800_obs;

 // setup output plot
 TH2D * limits = new TH2D("limits","limits",8,0,8,6,0,6);
 limits->GetXaxis()->SetTitle("m_{Z'} [GeV]");
 limits->GetYaxis()->SetTitle("m_{A} [GeV]");
 limits->GetZaxis()->SetTitle("#sigma_{95\% CL} / #sigma_{th}");
 limits->SetTitle("");
 // limits->SetMaximum(3000);
 limits->SetMarkerSize(1.6);

 // limits->GetYaxis()->SetTitleOffset(1.3);
 //limits->GetZaxis()->SetTitleOffset(0.75);
 // size of axis labels
 limits->GetXaxis()->SetTitleSize(0.04);
 limits->GetYaxis()->SetTitleSize(0.04);
 limits->GetZaxis()->SetTitleSize(0.035);
 limits->GetXaxis()->SetLabelSize(0.05);
 limits->GetYaxis()->SetLabelSize(0.05); 
 limits->GetZaxis()->SetLabelSize(0.025);

 // set the lables for the Xaxis (mZp)
 limits->GetXaxis()->SetBinLabel(1,"600");
 limits->GetXaxis()->SetBinLabel(2,"800");
 limits->GetXaxis()->SetBinLabel(3,"1000");
 limits->GetXaxis()->SetBinLabel(4,"1200");
 limits->GetXaxis()->SetBinLabel(5,"1400");
 limits->GetXaxis()->SetBinLabel(6,"1700");
 limits->GetXaxis()->SetBinLabel(7,"2000");
 limits->GetXaxis()->SetBinLabel(8,"2500");

 // set the lables for the Yaxis (mA0)
 limits->GetYaxis()->SetBinLabel(1,"300");
 limits->GetYaxis()->SetBinLabel(2,"400");
 limits->GetYaxis()->SetBinLabel(3,"500");
 limits->GetYaxis()->SetBinLabel(4,"600");
 limits->GetYaxis()->SetBinLabel(5,"700");
 limits->GetYaxis()->SetBinLabel(6,"800");

 // setup output observed plot
 TH2D * obslimits = (TH2D*) limits->Clone();

 // setup canvas
 //TCanvas * c = new TCanvas("c","",889,768);
 TCanvas * c = new TCanvas("c","",1);
 c->cd();
 gStyle->SetOptStat(0);
 gStyle->SetPaintTextFormat("2.1e");
 // c->SetLeftMargin(0.1);
 //c->SetRightMargin(0.11);

 Double_t limitval300[nMasses];
 Double_t limitval400[nMasses];
 Double_t limitval500[nMasses];
 Double_t limitval600[nMasses];
 Double_t limitval700[nMasses];
 Double_t limitval800[nMasses];

 Double_t limitval300_obs[nMasses];
 Double_t limitval400_obs[nMasses];
 Double_t limitval500_obs[nMasses];
 Double_t limitval600_obs[nMasses];
 Double_t limitval700_obs[nMasses];
 Double_t limitval800_obs[nMasses];

 Double_t xsecA0300[nMasses];
 Double_t xsecA0400[nMasses];
 Double_t xsecA0500[nMasses];
 Double_t xsecA0600[nMasses];
 Double_t xsecA0700[nMasses];
 Double_t xsecA0800[nMasses];

 Double_t explimit300[nMasses];
 Double_t explimit400[nMasses];
 Double_t explimit500[nMasses];
 Double_t explimit600[nMasses];
 Double_t explimit700[nMasses];
 Double_t explimit800[nMasses];
 
 Double_t obslimit300[nMasses];
 Double_t obslimit400[nMasses];
 Double_t obslimit500[nMasses];
 Double_t obslimit600[nMasses];
 Double_t obslimit700[nMasses];
 Double_t obslimit800[nMasses];

 Double_t eff_A0300[nMasses];
 Double_t eff_A0400[nMasses];
 Double_t eff_A0500[nMasses];
 Double_t eff_A0600[nMasses];
 Double_t eff_A0700[nMasses];
 Double_t eff_A0800[nMasses];
 Double_t efferr_A0300[nMasses];
 Double_t efferr_A0400[nMasses];
 Double_t efferr_A0500[nMasses];
 Double_t efferr_A0600[nMasses];
 Double_t efferr_A0700[nMasses];
 Double_t efferr_A0800[nMasses];

 for (Int_t n=0; n<nMasses; n++){
   getLimits(higgsCombineFiles_MA0300[n],limitval300[n],0.5); 
   getLimits(higgsCombineFiles_MA0400[n],limitval400[n],0.5); 
   getLimits(higgsCombineFiles_MA0500[n],limitval500[n],0.5); 
   getLimits(higgsCombineFiles_MA0600[n],limitval600[n],0.5); 
   getLimits(higgsCombineFiles_MA0700[n],limitval700[n],0.5); 
   getLimits(higgsCombineFiles_MA0800[n],limitval800[n],0.5); 
   //std::cout<<n<<" "<<limitval300[n]<< " "<<limitval400[n]<<limitval500[n]<< " "<<limitval600[n] <<limitval700[n]<< " "<<limitval800[n]<<std::endl;

   getLimits(higgsCombineFiles_MA0300[n],limitval300_obs[n],-1.); 
   getLimits(higgsCombineFiles_MA0400[n],limitval400_obs[n],-1.); 
   getLimits(higgsCombineFiles_MA0500[n],limitval500_obs[n],-1.); 
   getLimits(higgsCombineFiles_MA0600[n],limitval600_obs[n],-1.); 
   getLimits(higgsCombineFiles_MA0700[n],limitval700_obs[n],-1.); 
   getLimits(higgsCombineFiles_MA0800[n],limitval800_obs[n],-1.); 

   getXsec(theory_gz08,300,(Int_t)mass[n],xsecA0300[n]);
   getXsec(theory_gz08,400,(Int_t)mass[n],xsecA0400[n]);
   getXsec(theory_gz08,500,(Int_t)mass[n],xsecA0500[n]);
   getXsec(theory_gz08,600,(Int_t)mass[n],xsecA0600[n]);
   getXsec(theory_gz08,700,(Int_t)mass[n],xsecA0700[n]);
   getXsec(theory_gz08,800,(Int_t)mass[n],xsecA0800[n]);

   explimit300[n] = limitval300[n]/xsecA0300[n];
   explimit400[n] = limitval400[n]/xsecA0400[n];
   explimit500[n] = limitval500[n]/xsecA0500[n];
   explimit600[n] = limitval600[n]/xsecA0600[n];
   explimit700[n] = limitval700[n]/xsecA0700[n];
   explimit800[n] = limitval800[n]/xsecA0800[n];

   obslimit300[n] = limitval300_obs[n]/xsecA0300[n];
   obslimit400[n] = limitval400_obs[n]/xsecA0400[n];
   obslimit500[n] = limitval500_obs[n]/xsecA0500[n];
   obslimit600[n] = limitval600_obs[n]/xsecA0600[n];
   obslimit700[n] = limitval700_obs[n]/xsecA0700[n];
   obslimit800[n] = limitval800_obs[n]/xsecA0800[n];

   // fill limit plot
   limits->Fill(((Double_t)n+0.5),0.5,limitval300[n]/xsecA0300[n]);
   //   if(n==nMasses-1) limits->Fill(((Double_t)n+0.5),1.5,1968); // --hardcode in some numbers
   limits->Fill(((Double_t)n+0.5),1.5,limitval400[n]/xsecA0400[n]);
   limits->Fill(((Double_t)n+0.5),2.5,limitval500[n]/xsecA0500[n]);
   //if (n==1)        limits->Fill(((Double_t)n+0.5),3.5,199.4);
   //else if (n==3)   limits->Fill(((Double_t)n+0.5),3.5,78.8);
   limits->Fill(((Double_t)n+0.5),3.5,limitval600[n]/xsecA0600[n]);
   limits->Fill(((Double_t)n+0.5),4.5,limitval700[n]/xsecA0700[n]);
   limits->Fill(((Double_t)n+0.5),5.5,limitval800[n]/xsecA0800[n]);

   obslimits->Fill(((Double_t)n+0.5),0.5,limitval300_obs[n]/xsecA0300[n]);
   //if(n==nMasses-1) obslimits->Fill(((Double_t)n+0.5),1.5,1941.1);
   obslimits->Fill(((Double_t)n+0.5),1.5,limitval400_obs[n]/xsecA0400[n]);
   obslimits->Fill(((Double_t)n+0.5),2.5,limitval500_obs[n]/xsecA0500[n]);
   //if (n==1)        obslimits->Fill(((Double_t)n+0.5),3.5,196.7);
   //else if (n==3)   obslimits->Fill(((Double_t)n+0.5),3.5,77.6);
   obslimits->Fill(((Double_t)n+0.5),3.5,limitval600_obs[n]/xsecA0600[n]);
   obslimits->Fill(((Double_t)n+0.5),4.5,limitval700_obs[n]/xsecA0700[n]);
   obslimits->Fill(((Double_t)n+0.5),5.5,limitval800_obs[n]/xsecA0800[n]);


 }

 // only pick up the limits that are non-zero
 Double_t mass_400[8] = {600,800,1000,1200,1400,1700,2000,2500};
 Double_t mass_500[6] = {800,1200,1400,1700,2000,2500};
 Double_t mass_600[6] = {800,1200,1400,1700,2000,2500};
 Double_t mass_700[5] = {1200,1400,1700,2000,2500};
 Double_t mass_800[6] = {1000,1200,1400,1700,2000,2500};
 // mA0 400
 Double_t limitval_exp_400[8] = {explimit400[0],explimit400[1],explimit400[2],explimit400[3],explimit400[4],explimit400[5],explimit400[6],explimit400[7]};
 Double_t limitval_obs_400[8] = {obslimit400[0],obslimit400[1],obslimit400[2],obslimit400[3],obslimit400[4],obslimit400[5],obslimit400[6],obslimit400[7]};

 // mA0 500
 Double_t limitval_exp_500[6] = {explimit500[1],explimit500[3],explimit500[4],explimit500[5],explimit500[6],explimit500[7]};
 Double_t limitval_obs_500[6] = {obslimit500[1],obslimit500[3],obslimit500[4],obslimit500[5],obslimit500[6],obslimit500[7]};

 // mA0 600
 Double_t limitval_exp_600[6] = {explimit600[1],explimit600[3],explimit600[4],explimit600[5],explimit600[6],explimit600[7]};
 Double_t limitval_obs_600[6] = {obslimit600[1],obslimit600[3],obslimit600[4],obslimit600[5],obslimit600[6],obslimit600[7]};
 // for(int i=0; i<6;i++) std::cout<<i<<" "<<mass_600[i]<<" "<<limitval_exp_600[i]<<" "<<limitval_obs_600[i]<<std::endl;
 // mA0 700
 Double_t limitval_exp_700[5] = {explimit700[3],explimit700[4],explimit700[5],explimit700[6],explimit700[7]};
 Double_t limitval_obs_700[5] = {obslimit700[3],obslimit700[4],obslimit700[5],obslimit700[6],obslimit700[7]};
 // for(int i=0; i<5;i++) std::cout<<i<<" "<<mass_700[i]<<" "<<limitval_exp_700[i]<<" "<<limitval_obs_700[i]<<std::endl;

 // mA0 800
 Double_t limitval_exp_800[6] = {explimit800[2],explimit800[3],explimit800[4],explimit800[5],explimit800[6],explimit800[7]};
 Double_t limitval_obs_800[6] = {obslimit800[2],obslimit800[3],obslimit800[4],obslimit800[5],obslimit800[6],obslimit800[7]};

 limit300 = new TGraph(nMasses,mass,explimit300);
 limit400 = new TGraph(nMasses,mass_400,limitval_exp_400);
 limit500 = new TGraph(6,mass_500,limitval_exp_500);
 limit600 = new TGraph(6,mass_600,limitval_exp_600);
 limit700 = new TGraph(5,mass_700,limitval_exp_700);
 limit800 = new TGraph(6,mass_800,limitval_exp_800);
 
 limit300_obs = new TGraph(nMasses,mass,obslimit300);
 limit400_obs = new TGraph(nMasses,mass_400,limitval_obs_400);
 limit500_obs = new TGraph(6,mass_500,limitval_obs_500);
 limit600_obs = new TGraph(6,mass_600,limitval_obs_600);
 limit700_obs = new TGraph(5,mass_700,limitval_obs_700);
 limit800_obs = new TGraph(6,mass_800,limitval_obs_800);

 //styling
 limit300->GetXaxis()->SetTitle("m_{Z'} [GeV]");
 limit300->GetYaxis()->SetTitle("#sigma_{95\% CL} / #sigma_{th}");
 limit300->SetTitle("");
 // limit300->SetMaximum(3000);
 //limit300->SetMinimum(0.9);
 limit300->SetLineWidth(2);
 limit400->SetLineWidth(2);
 limit500->SetLineWidth(2);
 limit600->SetLineWidth(2);
 limit700->SetLineWidth(2);
 limit800->SetLineWidth(2);
 limit300->SetMarkerStyle(8);
 limit400->SetMarkerStyle(8);
 limit500->SetMarkerStyle(8);
 limit600->SetMarkerStyle(8);
 limit700->SetMarkerStyle(8);
 limit800->SetMarkerStyle(8);
 // set up colors to match Hbb
 limit300->SetLineColor(kBlack);
 limit400->SetLineColor(kCyan);
 limit500->SetLineColor(kGreen);
 limit600->SetLineColor(kBlue);
 limit700->SetLineColor(kYellow);
 limit800->SetLineColor(kMagenta);
 limit300->SetMarkerColor(kBlack);
 limit400->SetMarkerColor(kCyan);
 limit500->SetMarkerColor(kGreen);
 limit600->SetMarkerColor(kBlue);
 limit700->SetMarkerColor(kYellow);
 limit800->SetMarkerColor(kMagenta);


 //styling
 limit300_obs->GetXaxis()->SetTitle("m_{Z'} [GeV]");
 limit300_obs->GetYaxis()->SetTitle("#sigma_{95\% CL} / #sigma_{th}");
 limit300_obs->SetTitle("");
 //limit300_obs->SetMaximum(3000);
 //limit300_obs->SetMinimum(0.9);
 limit300_obs->SetLineWidth(2);
 limit400_obs->SetLineWidth(2);
 limit500_obs->SetLineWidth(2);
 limit600_obs->SetLineWidth(2);
 limit700_obs->SetLineWidth(2);
 limit800_obs->SetLineWidth(2);
 limit300_obs->SetMarkerStyle(8);
 limit400_obs->SetMarkerStyle(8);
 limit500_obs->SetMarkerStyle(8);
 limit600_obs->SetMarkerStyle(8);
 limit700_obs->SetMarkerStyle(8);
 limit800_obs->SetMarkerStyle(8);
 // set up colors to match Hbb
 limit300_obs->SetMarkerColor(kBlack);
 limit400_obs->SetMarkerColor(kCyan);
 limit500_obs->SetMarkerColor(kGreen);
 limit600_obs->SetMarkerColor(kBlue);
 limit700_obs->SetMarkerColor(kYellow);
 limit800_obs->SetMarkerColor(kMagenta);
 limit300_obs->SetLineColor(kBlack);
 limit400_obs->SetLineColor(kCyan);
 limit500_obs->SetLineColor(kGreen);
 limit600_obs->SetLineColor(kBlue);
 limit700_obs->SetLineColor(kYellow);
 limit800_obs->SetLineColor(kMagenta);
 


 // draw expected limits plot

 gStyle->SetHistMinimumZero(kFALSE);
limits->Draw("COLZ TEXT"); 

 //latex.DrawLatex(0.08,5.7,thestring);
 // save plot
 //CMS_lumi(c,false,0);
 // l1->Draw("same");
 //l2->Draw("same");
 // l3->Draw("same");
 //l4->Draw("same");
 c->cd();
 c->SaveAs(Form("%s/limits2D_2HDM_exp.png",outDir.Data()));
 c->SaveAs(Form("%s/limits2D_2HDM_exp.pdf",outDir.Data()));

 // draw observed limits plot
 gStyle->SetHistMinimumZero(kFALSE);
 obslimits->Draw("COLZ TEXT"); 
 //latex.DrawLatex(0.08,5.7,thestring);
 // save plot
 //CMS_lumi(c,false,0);
 // l1->Draw("same");
 //l2->Draw("same");
 //l3->Draw("same");
 //l4->Draw("same");
 c->cd();
 c->SaveAs(Form("%s/limits2D_2HDM_obs.png",outDir.Data()));
 c->SaveAs(Form("%s/limits2D_2HDM_obs.pdf",outDir.Data()));

 TLegend* leg = new TLegend(0.65,0.2,0.85,0.4,NULL,"brNDC"); // (x1,y1,x2,y2)
 leg->SetTextSize(0.046);
 leg->SetBorderSize(0);
 leg->SetLineColor(1);
 leg->SetLineWidth(1);
 leg->SetLineStyle(1);
 leg->SetFillColor(0);
 leg->SetFillStyle(0);
 leg->SetTextFont(62);
 leg->AddEntry(limit300,"m_{A} = 300 GeV","pl");
 leg->AddEntry(limit400,"m_{A} = 400 GeV","pl");
 leg->AddEntry(limit500,"m_{A} = 500 GeV","pl");
 leg->AddEntry(limit600,"m_{A} = 600 GeV","pl");
 leg->AddEntry(limit700,"m_{A} = 700 GeV","pl");
 leg->AddEntry(limit800,"m_{A} = 800 GeV","pl");
 leg->SetTextSize(0.03);

 TLine* line1 = new TLine();
 line1->SetX1(limit300->GetXaxis()->GetXmin());
 line1->SetY1(1.0);
 line1->SetX2(limit300->GetXaxis()->GetXmax());
 line1->SetY2(1.0);
 line1->SetLineColor(kRed);
 line1->SetLineWidth(2);

 // draw 1D comparisons --expected
 //c->Clear();
 // c->SetLogy(1);
 limit300->GetYaxis()->SetRangeUser(0.1, 1000);
 limit300->Draw("APL");
 limit300->GetYaxis()->SetRangeUser(0.1, 1000);
 limit400->Draw("PL SAME");
 limit500->Draw("PL SAME");
 limit600->Draw("PL SAME");
 limit700->Draw("PL SAME");
 limit800->Draw("PL SAME");
 leg->Draw("SAME");
 line1->Draw("SAME");
 //latex.DrawLatex(500,2100,thestring);
 //CMS_lumi(c,false,0);
 //l1->Draw("same");
 //l2->Draw("same");
 //l3->Draw("same");
 //l4->Draw("same");
 // c->cd();
 c->SaveAs(Form("%s/limits_comparison_2HDM_exp.png",outDir.Data()));
 c->SaveAs(Form("%s/limits_comparison_2HDM_exp.pdf",outDir.Data()));

 // draw 1D comparisons --observed
 c->Clear();
 // c->SetLogy(1);
 limit300_obs->GetYaxis()->SetRangeUser(0.1, 1000);
 limit300_obs->Draw("APL");
 limit300_obs->GetYaxis()->SetRangeUser(0.1, 1000);
 limit400_obs->Draw("PL SAME");
 limit500_obs->Draw("PL SAME");
 limit600_obs->Draw("PL SAME");
 limit700_obs->Draw("PL SAME");
 limit800_obs->Draw("PL SAME");
 leg->Draw("SAME");
 line1->Draw("SAME");
 //latex.DrawLatex(500,2100,thestring);
 //CMS_lumi(c,false,0);
 //l1->Draw("same");
 //l2->Draw("same");
 // /l3->Draw("same");
 //l4->Draw("same");
 // c->cd();
 c->SaveAs(Form("%s/limits_comparison_2HDM_obs.png",outDir.Data()));
 c->SaveAs(Form("%s/limits_comparison_2HDM_obs.pdf",outDir.Data()));

 // draw 1D comparisons --expected & observed
 // c->Clear();
 // c->SetLogy(1);
 limit300->SetLineStyle(9);
 limit400->SetLineStyle(9);
 limit500->SetLineStyle(9);
 limit600->SetLineStyle(9);
 limit700->SetLineStyle(9);
 limit800->SetLineStyle(9);
 limit300->GetYaxis()->SetRangeUser(0., 50);
 limit300->Draw("APL");
 limit400->Draw("PL SAME");
 limit500->Draw("PL SAME");
 limit600->Draw("PL SAME");
 limit700->Draw("PL SAME");
 limit800->Draw("PL SAME");
 limit300_obs->Draw("PL SAME");
 limit400_obs->Draw("PL SAME");
 limit500_obs->Draw("PL SAME");
 limit600_obs->Draw("PL SAME");
 limit700_obs->Draw("PL SAME");
 limit800_obs->Draw("PL SAME");
 //leg->AddEntry(limit300,"Expected, m_{A0} = 300 GeV","pl");
 //leg->AddEntry(limit300_obs,"Observed, m_{A0} = 300 GeV","pl");
 //leg->AddEntry(limit400,"Expected, m_{A0} = 400 GeV","pl");
 //leg->AddEntry(limit400_obs,"Observed, m_{A0} = 300 GeV","pl");
 //leg->AddEntry(limit500,"Expected, m_{A0} = 500 GeV","pl");
 //leg->AddEntry(limit500_obs,"Observed, m_{A0} = 300 GeV","pl");
 //leg->AddEntry(limit600,"Expected, m_{A0} = 600 GeV","pl");
 //leg->AddEntry(limit600_obs,"Observed, m_{A0} = 300 GeV","pl");
 //leg->AddEntry(limit700,"Expected, m_{A0} = 700 GeV","pl");
 //leg->AddEntry(limit700_obs,"Observed, m_{A0} = 300 GeV","pl");
 //leg->AddEntry(limit800,"Expected, m_{A0} = 800 GeV","pl");
 //leg->AddEntry(limit800_obs,"Observed, m_{A0} = 300 GeV","pl");
 leg->Draw("SAME");
 line1->Draw("SAME");
 //latex.DrawLatex(500,2100,thestring);
 //CMS_lumi(c,false,0);
 //l1->Draw("same");
 //l2->Draw("same");
 //l3->Draw("same");
 //l4->Draw("same");
 c->cd();
 c->SaveAs(Form("%s/limits_comparison_2HDM.png",outDir.Data()));
 c->SaveAs(Form("%s/limits_comparison_2HDM.pdf",outDir.Data()));


 delete c;


 // make plot with both expected and observed on same graph
 //TCanvas* cboth = new TCanvas("cboth","",889,768);
 TCanvas * cboth = new TCanvas("cboth","",1);
 cboth->cd();
 gStyle->SetHistMinimumZero(kFALSE);
 gStyle->SetOptStat(0);
 gStyle->SetPaintTextFormat("2.1e");
 gStyle->SetPalette(57);
 gStyle->SetFrameLineWidth(3);
 gStyle->SetPadRightMargin(0.109);
 gStyle->SetPadLeftMargin(0.13);
 cboth->SetLeftMargin(0.1);
 cboth->SetRightMargin(0.1);


 Double_t pad1_x1 = 0.01; 
 Double_t pad1_x2 = 0.98;
 Double_t pad1_y1 = 0.03;
 Double_t pad1_y2 = 1.00;
 Double_t pad2_x1 = pad1_x1;
 Double_t pad2_x2 = pad1_x2;
 Double_t pad2_y1 = pad1_y1-0.03;
 Double_t pad2_y2 = pad1_y2-0.03;

 gStyle->SetPaintTextFormat("2.1e");
 //TPad* p1 = new TPad("p1","",0,0.12,0.95,0.98);
 //TPad* p1 = new TPad("p1","",0,0.09,0.95,0.89); //x1,y1,x2,y2
 //TPad* p1 = new TPad("p1","",0.02,0.03,0.98,0.95); //x1,y1,x2,y2
 TPad* p1 = new TPad("p1","",pad1_x1,pad1_y1,pad1_x2,pad1_y2); //x1,y1,x2,y2
 p1->Draw();
 p1->cd();
 p1->SetLogz();
 gStyle->SetHistMinimumZero(kFALSE);
 gStyle->SetPaintTextFormat("2.1e");
 obslimits->Draw("COLZ"); 
 //obslimits->Draw("TEXT SAME");

 float obscontent; 
 double obsbinX;
 double obsbinY;
 TString obsbincon;
 TLatex *obsbintxt; 
 for (int obsbinx = 1; obsbinx <= obslimits->GetXaxis()->GetNbins(); obsbinx++){
    for (int obsbiny = 1; obsbiny <= obslimits->GetYaxis()->GetNbins(); obsbiny++){
       obscontent   = obslimits->GetBinContent(obsbinx,obsbiny);
       if (obsbinx==1 && obsbiny>=3) continue;
       if (obsbinx==2 && obsbiny>=5) continue;
       if (obscontent > 1000) obsbincon = TString::Format("%2.1e",obscontent);
       else obsbincon = TString::Format("%2.1e",obscontent);
       obsbinX      = obslimits->GetXaxis()->GetBinCenter(obsbinx);
       obsbinY      = obslimits->GetYaxis()->GetBinCenter(obsbiny);
       obsbintxt = new TLatex(obsbinX,obsbinY-0.08,obsbincon);
       obsbintxt->SetTextAlign(21);
       obsbintxt->SetTextSize(0.022);
       if(obscontent>0)obsbintxt->Draw("same");
    }
 }
 //latex.DrawLatex(0.08,5.7,thestring);
 p1->Update();
 gStyle->SetPaintTextFormat("2.1e");
 Double_t x1,y1,x2,y2;
 p1->GetRange(x1,y1,x2,y2);

 cboth->cd();
 TPad* p2 = new TPad("p2","",pad2_x1,pad2_y1,pad2_x2,pad2_y2); //x1,y1,x2,y2
 p2->SetFillStyle(0);
 p2->SetFillColor(0);
 p2->Draw();
 gStyle->SetHistMinimumZero(kFALSE);
 gStyle->SetPaintTextFormat("2.1e");
 p2->cd();
 p2->Range(x1,y1,x2,y2);
 gStyle->SetFrameLineWidth(3);
 TFrame *f = (TFrame*)cboth->FindObject("TFrame");
 Double_t px1 = f->GetX1();
 Double_t px2 = f->GetX2();
 Double_t py1 = f->GetY1()+0.23;
 Double_t py2 = f->GetY2()+0.23;

 limits->SetMarkerSize(1.4);
 limits->GetXaxis()->SetTitle("");
 limits->GetYaxis()->SetTitle("");
 limits->SetTitle("");
 limits->GetXaxis()->SetBinLabel(1,"");
 limits->GetXaxis()->SetBinLabel(2,"");
 limits->GetXaxis()->SetBinLabel(3,"");
 limits->GetXaxis()->SetBinLabel(4,"");
 limits->GetXaxis()->SetBinLabel(5,"");
 limits->GetXaxis()->SetBinLabel(6,"");
 limits->GetXaxis()->SetBinLabel(7,"");
 limits->GetXaxis()->SetBinLabel(8,"");
 limits->GetYaxis()->SetBinLabel(1,"");
 limits->GetYaxis()->SetBinLabel(2,"");
 limits->GetYaxis()->SetBinLabel(3,"");
 limits->GetYaxis()->SetBinLabel(4,"");
 limits->GetYaxis()->SetBinLabel(5,"");
 limits->GetYaxis()->SetBinLabel(6,"");

 float content; 
 double binX;
 double binY;
 TString bincon;
 TLatex *bintxt; 
 gStyle->SetPaintTextFormat("2.1e");
 for (int binx = 1; binx <= limits->GetXaxis()->GetNbins(); binx++){
    for (int biny = 1; biny <= limits->GetYaxis()->GetNbins(); biny++){
       content   = limits->GetBinContent(binx,biny);
       if (binx==1 && biny>=3) continue;
       if (binx==2 && biny>=5) continue;
       if (binx==limits->GetXaxis()->GetNbins()) bincon = TString::Format("(%2.1e)",content);
       else bincon = TString::Format("(%2.1e)",content);
       binX      = limits->GetXaxis()->GetBinCenter(binx);
       binY      = limits->GetYaxis()->GetBinCenter(biny);
       bintxt = new TLatex(binX,binY-0.08,bincon);
       bintxt->SetTextAlign(21);
       bintxt->SetTextSize(0.022);
       if(content>0)bintxt->Draw("same");
    }
 }

 //limits->Draw("TEXT SAME"); 
 p1->Update();
 gStyle->SetPaintTextFormat("2.1e");
 // redraw the frame around the histogram
 TLine l;
 l.SetLineWidth(3);
 l.DrawLine(px1,py2,px2,py2);
 l.DrawLine(px2,py1,px2,py2);
 l.DrawLine(px1,py1,px2,py1);
 l.DrawLine(px1,py1,px1,py2);

 //CMS_lumi(cboth,false,0);
 l1b->Draw("same");
 l2b->Draw("same");
 // l3b->Draw("same");
 //l4b->Draw("same");
 cboth->cd();
 cboth->SaveAs(Form("%s/limits2D_2HDM_ExpAndObs.png",outDir.Data()));
 cboth->SaveAs(Form("%s/limits2D_2HDM_ExpAndObs.pdf",outDir.Data()));
 delete cboth;


}

void getEff(TFile* file, int mA0, int mZp, Double_t & eff, Double_t & eff_err){
  TH2F* sigeff = (TH2F*)file->Get("effhisto");
  if (sigeff!=(TH2F*)NULL){
    Int_t binX = sigeff->GetXaxis()->FindBin(mZp);
    Int_t binY = sigeff->GetYaxis()->FindBin(mA0);  
    eff = sigeff->GetBinContent(binX,binY);
    eff_err = sigeff->GetBinError(binX,binY);
  }
  else{
    eff = 0;
    eff_err = 0;
    std::cout << " No eff for mA0 = " << mA0 << " and mZp = " << mZp << std::endl;
  }
}


void getXsec(TFile* file, int mA0, int mZp, Double_t & xsec){

  TH2F* xsecs = (TH2F*)file->Get("xsec1"); 
  if (xsecs!=(TH2F*)NULL){
     Int_t binX = xsecs->GetXaxis()->FindBin(mZp);
     Int_t binY = xsecs->GetYaxis()->FindBin(mA0);
     xsec = xsecs->GetBinContent(binX,binY); 
  }
  else{
   xsec = 1;
   std::cout << "Couldn't find xsec histogram" << std::endl;
  }
} 


void getXsecBaryo(TFile* file, int mchi, int mZp, Double_t & xsec){
  //  file->Print("V");
  TGraph* xsecs = (TGraph*)file->Get(Form("zp%d",mZp)); 
  if (xsecs!=(TGraph*)NULL){
     xsec = xsecs->Eval(mchi);
     //     xsec = xsecs->GetBinContent(binX);
     //     std::cout<<mZp<< " "<<mchi<< " "<<xsec<<std::endl;
  }
  else{
   xsec = 1;
   std::cout << "Couldn't find xsec histogram" << std::endl;
  }
} 


void getXsecScalar(TFile* file, int mchi, int mZp, Double_t & xsec){
  //  file->Print("V");
  TGraph* xsecs = (TGraph*)file->Get(Form("zp%d",mZp)); 
  if (xsecs!=(TGraph*)NULL){
     xsec = xsecs->Eval(mchi);
     //     xsec = xsecs->GetBinContent(binX);
     std::cout<<mZp<< " "<<mchi<< " "<<xsec<<std::endl;
  }
  else{
   xsec = 1;
   //   std::cout << "Couldn't find xsec histogram" << std::endl;
  }
} 


void getLimits(TFile* file, Double_t & Limit, Double_t quantile){

  Double_t limit;
  Float_t quantileExpected;

  TBranch *b_limit;
  TBranch *b_quantileExpected;

  TTree* tree = (TTree*)file->Get("limit");
  if (tree!=(TTree*)NULL){
 
    tree->SetBranchAddress("limit", &limit, &b_limit);
    tree->SetBranchAddress("quantileExpected", &quantileExpected, &b_quantileExpected);
 
    Limit = 0;
    int nentries = tree->GetEntries();
    for (int entry = 0; entry < nentries; entry++){
      tree->GetEntry(entry);
      std::cout << "Quantile = " << quantileExpected << std::endl;
      std::cout << "Limit    = " << limit << std::endl;
      if (quantileExpected==quantile) Limit=limit;
    }

  }// end valid tree
  else Limit = 0;
  //std::cout << "Limit    = " << Limit << std::endl;
  
  delete tree;

}

