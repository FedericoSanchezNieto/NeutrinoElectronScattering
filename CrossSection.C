#include "NeutrinoElectronScattering.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TRandom.h"
#include <iostream> 

void CrossSection(double EnuIn = 10.,double Emin = 0 ){

  TH1F *histang = new TH1F("histang","",1000,0.,3.141592/2);
  TH1F *histanganu = new TH1F("histanganu","",1000,0.,3.141592/2.);
  TH1F *histangl = new TH1F("histangl","",1000,0.,3.141592/2);
  TH1F *histanganul = new TH1F("histanganul","",1000,0.,3.141592/2.);
  
  TH1F *histEang = new TH1F("histEang","",1000,0.,1.);
  TH1F *histEanganu = new TH1F("histEanganu","",1000,0.,1.); 
  TH1F *histEangl = new TH1F("histEangl","",1000,0.,1.);
  TH1F *histEanganul = new TH1F("histEanganul","",1000,0.,1.); 

  TH1F *histE = new TH1F("histE","",500,0.,TMath::Abs(EnuIn));
  TH1F *histEanu = new TH1F("histEanu","",500,0.,TMath::Abs(EnuIn)); 
  TH1F *histEl = new TH1F("histEl","",500,0.,TMath::Abs(EnuIn));
  TH1F *histEanul = new TH1F("histEanul","",500,0.,TMath::Abs(EnuIn)); 

  TH1F *histEMC = new TH1F("histEMC","",500,0.,TMath::Abs(EnuIn));
  TH1F *histEMCanu = new TH1F("histEMCanu","",500,0.,TMath::Abs(EnuIn));
  TH1F *histEMCl = new TH1F("histEMCl","",500,0.,TMath::Abs(EnuIn));
  TH1F *histEMCanul = new TH1F("histEMCanul","",500,0.,TMath::Abs(EnuIn));

  TH2F *histEang2vsang = new TH2F("histEang2vsang"," ",100,0.,3.141592/2,500,0.,1.2);
  TH2F *histEang2vsangl = new TH2F("histEang2vsangl"," ",100,0.,3.141592/2,500,0.,1.2); 

  TH2F *histEvsang = new TH2F("histEvsang"," ",100,0.,3.141592/2,100,0.,TMath::Abs(EnuIn));
  TH2F *histEvsangl = new TH2F("histEvsangl"," ",100,0.,3.141592/2,100,0.,TMath::Abs(EnuIn)); 

  TH2F *ErecovsE = new TH2F("ErecovsE","",100,0.,TMath::Abs(EnuIn),100,0.,TMath::Abs(EnuIn)); 
  
  TRandom r;

  NeutrinoElectronScattering nul; 

  int neutrinotype = neutrinoelectron; 
  
  double Enu;
  
  for(int i = 0; i < 1e+6; i++ ) {

    double a,b;

    if( EnuIn < 0 )
      Enu = r.Uniform(Emin,-EnuIn);
    else
      Enu = EnuIn;
	
    double xsecttot = nul.Totalcrosssection(Enu,neutrinotype);
    double anuxsecttot = nul.Totalcrosssection(Enu,-neutrinotype);

    
    nul.GENcrosssection(Enu,a,b,neutrinotype); 

    histEMC->Fill(a,1./1.e+8); 

    nul.GENcrosssection(Enu,a,b,-neutrinotype); 

    histEMCanu->Fill(a,1./1.e+8);
    
    double Ee  = r.Uniform(nul.electronmass(),Enu);
    
    double xs = nul.diffcrosssection(Enu,Ee,neutrinotype)/xsecttot/1.e+8;
    double axs = nul.diffcrosssection(Enu,Ee,-neutrinotype)/anuxsecttot/1.e+8;

    double cos = nul.GetCosine(Enu,Ee,neutrinotype);
    double theta = acos(cos); 
    double Etheta2 = Ee*theta*theta;


    ErecovsE->Fill(Enu,nul.GetEnu(cos,Ee,neutrinotype)); 
    
    //  if( abs(cos) > 1 ) continue;
    
    if(xs > 0 ) {
      histEang->Fill(Etheta2,xs);
      histE->Fill(Ee,xs);
      histang->Fill(theta,xs);
      histEang2vsang->Fill(theta,Etheta2,xs);
      histEvsang->Fill(theta,Ee,xs); 
    }

    if(axs > 0 ) {
      histEanganu->Fill(Etheta2,axs);
      histEanu->Fill(Ee,axs);
      histanganu->Fill(theta,axs); 
    }

    xsecttot = nul.Totalcrosssection(Enu,neutrinomuon);
    anuxsecttot = nul.Totalcrosssection(Enu,-neutrinomuon);

    
    nul.GENcrosssection(Enu,a,b,neutrinomuon); 

    histEMCl->Fill(a,1./1.e+8); 

    nul.GENcrosssection(Enu,a,b,-neutrinomuon); 

    histEMCanul->Fill(a,1./1.e+8); 
    
    xs = nul.diffcrosssection(Enu,Ee,neutrinomuon)/xsecttot/1.e+8;
    axs = nul.diffcrosssection(Enu,Ee,-neutrinomuon)/anuxsecttot/1.e+8;

    if(xs > 0 ) {
      histEangl->Fill(Etheta2,xs);
      histEl->Fill(Ee,xs);
      histangl->Fill(theta,xs);
      histEang2vsangl->Fill(theta,Etheta2,xs);
      histEvsangl->Fill(theta,Ee,xs); 
    }

    if(axs > 0 ) {
      histEanganul->Fill(Etheta2,axs);
      histEanul->Fill(Ee,axs);
      histanganul->Fill(theta,axs); 
    }
 
    
  }


  TH1F *histEnu = new TH1F("histEnu","",1000,0.,5000.);  

  for(int i = 0;i < 1000;i++) {
    double E = histEnu->GetBinCenter(i+1);
    double xs = nul.Totalcrosssection(E,neutrinotype)*1.e-27;
    histEnu->SetBinContent(i+1,xs);
  }
  
  TCanvas *c1 = new TCanvas("c1","",800,800);
  histEnu->Draw();
  
  c1->Update();

  TCanvas *c2 = new TCanvas("c2","",800,800);
  histEang2vsang->Draw("colz");
  histEang2vsangl->Draw("cont2 same");
  c2->Update();
  
  TCanvas *c = new TCanvas("c","",800,800);
  c->Divide(2,2);
  c->cd(1); 
  histang->Draw("hist L");
  histangl->Draw("hist L same");
  histangl->SetLineStyle(10);
  histanganu->Draw("hist L same");
  histanganul->Draw("hist L same");
  histanganu->SetLineColor(2);
  histanganul->SetLineColor(2);
  histanganul->SetLineStyle(10);
  c->cd(2);
  histEanu->Draw("hist L");
  histEanul->Draw("hist L same");
  histEanu->SetLineColor(2);
  histEanul->SetLineColor(2);
  histEanul->SetLineStyle(10);
  histE->Draw("hist L same");
  histEl->Draw("hist L same");
  histEl->SetLineStyle(10);

  c->cd(3);
  histEang->Draw("hist L");
  histEangl->Draw("hist L same");
  histEangl->SetLineStyle(10);
  histEanganu->Draw("hist L same");
  histEanganul->Draw("hist L same");
  histEanganu->SetLineColor(2);
  histEanganul->SetLineColor(2);
  histEanganul->SetLineStyle(10);
  c->cd(4);
  histEMCanu->Draw("hist");
  histEMC->Draw("hist same");
  c->Update();

  TCanvas *c3 = new TCanvas("c3","",800,800);
  histEvsang->Draw("colz");
  histEvsang->GetXaxis()->SetTitle("#theta_{e} (rads)");
  histEvsang->GetYaxis()->SetTitle("E_{e} (MeV)"); 
  c3->Update(); 

  TCanvas *c4 = new TCanvas("c4","",800,800);
  ErecovsE->Draw("colz");
  ErecovsE->GetXaxis()->SetTitle("E_{#nu} (MeV)");
  ErecovsE->GetYaxis()->SetTitle("E^{Rec}_{#nu} (MeV)"); 
  c4->Update(); 


  
}

