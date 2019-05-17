#include "Riostream.h"
#include <string>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map> 
#include <algorithm>
#include <functional>
#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>
#include <limits>
#include <getopt.h>
using namespace std;
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>


#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TRandom.h"
#include <TMath.h>
#include "TSpline.h"


void smearing2(double muonSmear = 0.04, double pipmSmear = 0.0, double neutronSmear = 0., double protonSmear = 0., double pi0Smear = 0.,int fileLimit = 100,TString filepath = "/dune/app/users/gyang/genie-0.0-0.101/") {
   double neutronMass = 0.939565;
double protonMass = 0.938272;
double pipmMass = 0.13957;
double pi0Mass = 0.13498;
double muonMass = 0.105658;
double electronMass = 0.0005;

  // Getting the shifted spectra
     TFile fluxS("/dune/app/users/gyang/lblpwgtools/code/CAFAna/CAFAna/Systs/flux_shifts.root");
     TH1F* shift_numu[10];
     TSpline3* spline[10];
     for(Int_t i=0;i<10;i++){
     shift_numu[i] = (TH1F*)fluxS.Get(Form("syst%d/ND_numu_FHC",i));
     spline[i] = new TSpline3(shift_numu[i]);
     }
                     

 TFile *rootFile1 = new TFile(Form("Smearing-mu%0.2f-pipm%0.2f-n%0.2f-p%0.2f-pi0%0.2f-noKEneutron.root",muonSmear,pipmSmear,neutronSmear,protonSmear,pi0Smear),"RECREATE");
 
   int num=0;
   int num1=0;
   string line;
   ifstream myfile(filepath+"test1.txt");
 while (getline(myfile, line)){
        ++num;
        ++num1;
        }
    num = fileLimit;
    std::cout << "Number of lines in text file: " << num1 <<endl;;
   if(num>num1) exit(1);
   const int n=num;
   const int n1=num1;

  TH1F* histo[n];
  TH1F* histo1[n];
  TH1F* histo1s[n];
   TH1F* xpos[40];
   TH1F* xposm[40];
   TH1F* xposms[40]; 
  TH2F* histo2;

   double EE[n];
   double bins[n];
   double binsm[n];
   double binsms[n];
  int c=0;
   double max=0;
   double maxm=0;
   double min=0;
   double minm=0;
   double maxms=0;
   double minms=0;
   double mean[40];
   double rms[40];
    int totNFile = 0;
   string input1;
   ifstream infile;
   infile.open (filepath+"test1.txt");
   
  histo2 =  new TH2F("Energy vs Days","Energy vs Days",40,0,10,num,0,num);
  while(infile >> input1) // To get you all the lines.
       {       
       
       totNFile ++;
    if(totNFile> fileLimit) break;
 
         
  
  cout << input1 <<endl;
  TFile input(Form("%s",input1.c_str()));
  TTree* h1 = (TTree*) input.Get("gRooTracker");
           
  Int_t StdHepN;
  Int_t StdHepPdg[1000];
  Int_t StdHepStatus[1000];
  Double_t vtx[4]={};
  Double_t ivtx[1000][4]={};
  Double_t imom[1000][4]={};
  Bool_t flag=false;
  Bool_t wflag=false;
  Int_t rightSign=0,wrongSign=0;
    Int_t StdHepFm[50000];     // stdhep-like particle array: first mother
  
  h1->SetBranchAddress("StdHepPdg", &StdHepPdg);
  h1->SetBranchAddress("StdHepStatus", &StdHepStatus);
  h1->SetBranchAddress("EvtVtx", &vtx);
  h1->SetBranchAddress("StdHepX4", &ivtx);
  h1->SetBranchAddress("StdHepP4", &imom); 
  h1->SetBranchAddress("StdHepN", &StdHepN);
    h1->SetBranchAddress("StdHepFm",&StdHepFm);
  
  double inFV=0;
  double inFV1=0;
  double inFV2=0;
  int allEvt=0;
  
  double vetoY=0.10;
  double vetoX=0.10;
  
  rootFile1->cd();
  Int_t nentries = (Int_t) h1->GetEntries();
  //cout<< nentries <<endl;
  histo[c] =  new TH1F(Form("Energy Spectrum Day %i",c+1),Form("Energy Spectrum Day %i",c+1),40,0,10);
  histo1[c] =  new TH1F(Form("E_reco Spectrum Day %i",c+1),Form("E_reco Spectrum Day %i",c+1),40,0,10);
  //histo1s[c] =  new TH1F(Form("Muon Energy Spectrum Day %i - Smearing 4%%",c+1),Form("Muon Energy Spectrum Day %i - Smearing 4%%",c+1),40,0,10);

  for(Int_t i=0; i<nentries ; i++) 
  {

      
      
      h1->GetEntry(i);
      
      double muonMom = 0;
      double protonMom = 0;
      double pipmMom = 0;
      double pi0Mom = 0;
      double neutronMom = 0;
      
      double nuE = 0;
      double Efsmuon = 0;
      double Efsproton = 0;
      double Efspipm = 0;
      double Efspi0 = 0;
      double Efsneutron = 0;
      
      if(vtx[0]>-1.1 && vtx[0]<1.1 && vtx[1]>-1.1 && vtx[1]<1.1 && vtx[2]<0.9 && vtx[2]>-0.9 ){

      inFV1++;
      for(int ip=0; ip<StdHepN; ip++) {
      inFV2++;
      
      
      
      
      if(StdHepPdg[ip]==14 && StdHepStatus[ip]== 0 && StdHepFm[ip]==-1){
      histo[c]->Fill(imom[ip][3]);
      histo2->Fill(imom[ip][3],c); 
      nuE=imom[ip][3];        
      inFV++;
      }
      
      
      

   if(StdHepPdg[ip]==13 && StdHepFm[ip]==0 && StdHepStatus[ip]==1){
     muonMom = TMath::Sqrt(imom[ip][0]*imom[ip][0]+imom[ip][1]*imom[ip][1]+imom[ip][2]*imom[ip][2]);
     Efsmuon=TMath::Sqrt(TMath::Power(gRandom->Gaus(1,muonSmear)* muonMom,2) + TMath::Power(muonMass,2));
   }

    // take out neutron kinetic energy
    if(StdHepPdg[ip]==2112 && StdHepStatus[ip]==1){
      //neutronMom = TMath::Sqrt(imom[ip][0]*imom[ip][0]+imom[ip][1]*imom[ip][1]+imom[ip][2]*imom[ip][2]);
      //Efsneutron=TMath::Sqrt(TMath::Power(gRandom->Gaus(1,neutronSmear)*neutronMom,2) + TMath::Power(neutronMass,2));
      // you lose all this energy
      Efsneutron = imom[ip][3] - neutronMass;
    }
    
    if(StdHepPdg[ip]==2212 && StdHepStatus[ip]==1){
      protonMom = TMath::Sqrt(imom[ip][0]*imom[ip][0]+imom[ip][1]*imom[ip][1]+imom[ip][2]*imom[ip][2]);
      Efsproton=TMath::Sqrt(TMath::Power(gRandom->Gaus(1,protonSmear)*protonMom,2) + TMath::Power(protonMass,2));
    }

    if(abs(StdHepPdg[ip])==211 && StdHepStatus[ip]==1){
      pipmMom = TMath::Sqrt(imom[ip][0]*imom[ip][0]+imom[ip][1]*imom[ip][1]+imom[ip][2]*imom[ip][2]);
      Efspipm=TMath::Sqrt(TMath::Power(gRandom->Gaus(1,pipmSmear)*pipmMom,2) + TMath::Power(pipmMass,2));
    }
    
    if(abs(StdHepPdg[ip])==111 && StdHepStatus[ip]==1){
      pi0Mom = TMath::Sqrt(imom[ip][0]*imom[ip][0]+imom[ip][1]*imom[ip][1]+imom[ip][2]*imom[ip][2]);
      Efspi0=TMath::Sqrt(TMath::Power(gRandom->Gaus(1,pi0Smear)*pi0Mom,2) + TMath::Power(pi0Mass,2));
    }
    
      
      //histo1s[c]->Fill(gRandom->Gaus(1,0.04)*imom[ip][3]);
      //cout<<"smearing: "<<gRandom->Gaus(1,0.04)*imom[ip][3]<<endl;
      //cout<<"energy: "<<imom[ip][3]<<endl;
      
      }
      if(nuE>0){
      histo1[c]->Fill(Efsmuon + Efsproton + Efspipm + Efspi0 -neutronMass );
      }
      
      }
     
      
      
      
      
    
      
     allEvt ++;  
      
  }
           
  histo[c]->Write();
  histo1[c]->Write();
  //histo1s[c]->Write();

  c++;
        }
        
  histo2->Write();
  
      TH1F* xp =  new TH1F("Mean Energy Spectrum per Day","Mean Energy Spectrum per Day",40,0,10);
      TH1F* xpm =  new TH1F("Mean E_reco Energy Spectrum per Day","Mean E_reco Energy Spectrum per Day",40,0,10);
      //TH1F* xpms =  new TH1F("Mean Muon Energy Spectrum per Day - Smearing 4%","Mean Muon Energy Spectrum per Day - Smearing 4%",40,0,10);
      TH1F* xp1 =  new TH1F("Mean RMS per Day","Mean RMS per Day",40,0,10);
      TH1F* xp1m =  new TH1F("Mean E_reco RMS per Day","Mean E_reco RMS per Day",40,0,10);
      //TH1F* xp1ms =  new TH1F("Mean Muon RMS per Day - Smearing 4%","Mean Muon RMS per Day - Smearing 4%",40,0,10);
  
  
  for(int b=0;b<40;b++) 
  {
  
  //cout<<"N: "<<b<<endl;
             
  for(int p=0;p<num;p++)      
  {
  
  bins[p]=histo[p]->GetBinContent(b+1);
  binsm[p]=histo1[p]->GetBinContent(b+1);
  //binsms[p]=histo1s[p]->GetBinContent(b+1);
  // cout<<"Bin: "<<bins[p]<<endl;
   //cout<<"Bin 1: "<<binsm[p]<<endl;
  
  }
  
   std::vector<double> v(bins, bins + num);
   max = *max_element(v.begin(), v.end());
   min = *min_element(v.begin(), v.end());
   
   std::vector<double> v1(binsm, binsm + num);
   maxm = *max_element(v1.begin(), v1.end());
   minm = *min_element(v1.begin(), v1.end());
   
   //std::vector<double> v1s(binsms, binsms + num);
   //maxms = *max_element(v1s.begin(), v1s.end());
   //minms = *min_element(v1s.begin()+1, v1s.end());
   
//   cout<<"Max value: "<<max<<endl;
//   cout<<"Min value: "<<min<<endl;
  cout<<"________________________"<<endl;
   xpos[b] =  new TH1F(Form("Events vs Days (%0.2f GeV)",0.25*(b+1)),Form("Events vs Days (%0.2f GeV)",0.25*(b+1)),20,min,max+(max-min)/20);
   xposm[b] =  new TH1F(Form("E_reco Events vs Days (%0.2f GeV)",0.25*(b+1)),Form("E_reco Events vs Days (%0.2f GeV)",0.25*(b+1)),20,minm,maxm+(max-min)/20);
   //xposms[b] =  new TH1F(Form("Muon Events vs Days (%0.2f GeV) - Smearing 4%%",0.25*(b+1)),Form("Muon Events vs Days (%0.2f GeV) - Smearing 4%%",0.25*(b+1)),20,minms,maxms);
   for(int r=0;r<num;r++) {xpos[b]->Fill(bins[r]); xposm[b]->Fill(binsm[r]);}
   xpos[b]->Write();
   xposm[b]->Write();
   //xposms[b]->Write();
   //mean[b]=xpos[b]->GetMean();
   //rms[b]=xpos[b]->GetRMS();
   
   
      xp->SetBinContent(b+1,xpos[b]->GetMean());
      xp->SetBinError(b+1,xpos[b]->GetRMS());
      xpm->SetBinContent(b+1,xposm[b]->GetMean());
      xpm->SetBinError(b+1,xposm[b]->GetRMS());
      //xpms->SetBinContent(b+1,xposms[b]->GetMean());
      //xpms->SetBinError(b+1,xposms[b]->GetRMS());
      
      xp1->SetBinContent(b+1,xpos[b]->GetRMS());
      xp1m->SetBinContent(b+1,xposm[b]->GetRMS());
      //xp1ms->SetBinContent(b+1,xposms[b]->GetRMS());
  }
  

  TH1F* xp_shift[10];
  for (Int_t i=0; i<10; i++){
    xp_shift[i] =  new TH1F(Form("shfited_%d",i),Form("shifted_%d",i),40,0,10);
  }

  // use spline to get the shifted spectra
     for(Int_t ii=0;ii<10;ii++){
         for(Int_t i=0;i<xpm->GetNbinsX();i++){
               xp_shift[ii] -> SetBinContent(i+1, xpm->GetBinContent(i+1) * ( 1 + spline[ii] -> Eval(xpm->GetBinCenter(i+1) ) ) );
                   }
                     }
                     

      //gPad->SaveAs("plot.pdf");
      xp->Write();
      xpm->Write();
      //xpms->Write();
            for(Int_t i=0;i<10;i++) xp_shift[i]->Write();

      //gPad->SaveAs("plot1.pdf");
      xp1->Write();
      xp1m->Write();
      //xp1ms->Write();

   rootFile1->Close();
   
   }
   

