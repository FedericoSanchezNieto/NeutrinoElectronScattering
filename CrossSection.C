#include "NeutrinoElectronScattering.h" 


void CrossSection(double EnuIn = 10.){

  TH1F *histang = new TH1F("histang","",1000,0.,3.141592/2);
  TH1F *histanganu = new TH1F("histanganu","",1000,0.,3.141592/2.); 

  TH1F *histEang = new TH1F("histEang","",1000,0.,1.);
  TH1F *histEanganu = new TH1F("histEanganu","",1000,0.,1.); 

  TH1F *histE = new TH1F("histE","",1000,0.,TMath::Abs(EnuIn));
  TH1F *histEanu = new TH1F("histEanu","",1000,0.,TMath::Abs(EnuIn)); 

  TH1F *histEMC = new TH1F("histEMC","",1000,0.,TMath::Abs(EnuIn));
  TH1F *histEMCanu = new TH1F("histEMCanu","",1000,0.,TMath::Abs(EnuIn));

  TH2F *histEang2vsang = new TH2F("histEang2vsang"," ",100,0.,3.141592/2,500,0.,1.2); 

  TH2F *histEvsang = new TH2F("histEvsang"," ",100,0.,3.141592/2,100,0.,TMath::Abs(EnuIn)); 
  
  TRandom r;

  NeutrinoElectronScattering nul; 

  double Enu;
  
  for(int i = 0; i < 1e+7; i++ ) {

    double a,b;

    if( EnuIn < 0 )
      Enu = r.Uniform(0.,-EnuIn);
    else
      Enu = EnuIn;
	
    double xsecttot = nul.Totalcrosssection(Enu,neutrinoelectron);
    double anuxsecttot = nul.Totalcrosssection(Enu,-neutrinoelectron);

    
    nul.GENcrosssection(Enu,a,b,neutrinoelectron); 

    histEMC->Fill(a,1./1.e+8); 

    nul.GENcrosssection(Enu,a,b,-neutrinoelectron); 

    histEMCanu->Fill(a,1./1.e+8); 
    
    double Ee  = r.Uniform(nul.electronmass(),Enu);
    
    double xs = nul.diffcrosssection(Enu,Ee,neutrinoelectron)/xsecttot/1.e+8;
    double axs = nul.diffcrosssection(Enu,Ee,-neutrinoelectron)/anuxsecttot/1.e+8;

    double cos = nul.GetCosine(Enu,Ee,neutrinoelectron);
    double theta = acos(cos); 
    double Etheta2 = Ee*theta*theta; 
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
    
  }


  TH1F *histEnu = new TH1F("histEnu","",1000,0.,5000.);  

  for(int i = 0;i < 1000;i++) {
    double E = histEnu->GetBinCenter(i+1);
    double xs = nul.Totalcrosssection(E,neutrinoelectron)*1.e-27;
    histEnu->SetBinContent(i+1,xs);
  }
  
  TCanvas *c1 = new TCanvas("c1","",800,800);
  histEnu->Draw();
  
  c1->Update();

  TCanvas *c2 = new TCanvas("c2","",800,800);
  histEang2vsang->Draw("colz");
  c2->Update();
  
  TCanvas *c = new TCanvas("c","",800,800);
  c->Divide(2,2);
  c->cd(1); 
  histang->Draw("hist L");
  histanganu->Draw("hist L same");
  histanganu->SetLineColor(2);
  c->cd(2);
  histEanu->Draw("hist L");
  histE->Draw("hist L same");
  histEanu->SetLineColor(2);
  histEMC->Draw("hist L same");
  histEMC->SetLineColor(3);
  histEMCanu->Draw("hist L same");
  histEMCanu->SetLineColor(4);
  c->cd(3);
  histEang->Draw("hist L");
  histEanganu->Draw("hist L same");
  histEanganu->SetLineColor(2);
  c->cd(4);
  histEMCanu->Draw("hist");
  histEMC->Draw("hist same");
  c->Update();

  TCanvas *c3 = new TCanvas("c3","",800,800);
  histEvsang->Draw("colz");
  histEvsang->GetXaxis()->SetTitle("#theta_{e} (rads)");
  histEvsang->GetYaxis()->SetTitle("E_{e} (MeV)"); 
  c3->Update(); 
  
}

