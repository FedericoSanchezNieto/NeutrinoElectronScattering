#define neutrinoelectron  12
#define neutrinomuon      14
#define neutrinotau       16

#include <stdlib.h>
#include <math.h>
#include <iostream>

class NeutrinoElectronScattering {
  
 private:
  double sin2w ; 
  double gv; 
  double ga;
  double me;
  double mmu;
  double mtau;
  
  double Pi;
  double GF;
  double GF2;
  double MeV2mbarn;

  
 public:             // Access specifier
  NeutrinoElectronScattering () {

    sin2w = 0.23121; 
    gv = -1./2.+2.*sin2w;
    ga = -1./2.;
    me = 0.510998950; // MeV
    mmu = 105.6583755; // MeV
    mtau = 1776.86; // MeV
  
    Pi = 3.1415927;
    GF = 1.166378e-5; //  GeV^-2 
    GF2 =  GF*GF;
    MeV2mbarn = 1./2.56819e+6;   

    std::cout<< " Neutrino-Electron scattering Initialised " << std::endl;
  }

  double electronmass(void) { return me;}
  double muonmass(void) { return mmu;}
  double taumass(void) { return mtau;}
  double S(double Enu) { return (Enu+electronmass())*(Enu+electronmass())-Enu*Enu;}
  double Y(double Enu,double Ee) { return (Ee-electronmass())/Enu; }
  double Normalization(void) { return GF2/Pi*MeV2mbarn;}
  
  double diffcrosssection(double Enu, double El, int neutrino){
    bool aneut = false; 

    if( neutrino < 0 ) aneut = true;  // It is anti-neutrino 

    double xs = 0;

    if( El < electronmass() ) return 0; 
    
    if( abs(neutrino) ==  neutrinoelectron )
      xs = nuee(Enu,El,aneut);
    else
      xs = nule(Enu,El,aneut); 

    if( xs < 0 ) xs = 0;
    
    return xs; 
  }

  double GENcrosssection(double Enu, double &El,double &cos, int neutrino){
    bool aneut = false; 
    
    if( neutrino < 0 ) aneut = true;  // It is anti-neutrino
    
    double xs = 0;

    if( abs(neutrino) ==  neutrinoelectron )
      xs = GENnuee(Enu,El,cos,aneut);
    else
      xs = GENnule(Enu,El,cos,aneut); 

    //  if( aneut ) std::cout << xs << std::endl;
    
    if( xs < 0 ) xs = 0;

    return xs;
  }


     
  double Totalcrosssection(double Enu, int neutrino){
    bool aneut = false; 

    if( neutrino < 0 ) aneut = true;  // It is anti-neutrino 

    double xs = 0;

    if( abs(neutrino) ==  neutrinoelectron )
      xs =  nueeInt(Enu,aneut);
    else
      xs = nuleInt(Enu,aneut); 
    
    if( xs < 0 ) xs = 0;

    return xs;
  }

  
  double GetCosine(double Enu,double El,int neutrino){

    double mass = me; // This is only for electron in final state.

    double Pl = sqrt(El*El-mass*mass); 
   
    return (El*mass+Enu*El-mass*mass-Enu*mass)/Enu/Pl;

  }

  double GetEnu(double Cosine,double El,int neutrino){

    double mass = me; // This is only for electron in final state.

    double Pl = sqrt(El*El-mass*mass); 
   
    return (El*mass-mass*mass)/(-El+mass+Pl*Cosine);

  }
  
  
 private:
  
  
  double GenerateY(double CLL,double CLR,double N) {
    double r = (double)rand()/(double)RAND_MAX;
    
    double A = CLL*CLL/N;
    double B = CLR*CLR/N;
    double C = r; 

    // x = 1-y 
    // B/3 x^3 + A x + ( C - B/3 - A ) = 0
    //  
    //

    double a2 = 3.*A/B;
    double a3 = 3.*(C-B/3.-A)/B;

    double Q = a2/3.;
    double R = -a3/2.;
    
    double discriminant = Q*Q*Q+R*R;

    double S1 = pow(R+sqrt(discriminant),1./3.);

    double S2 = -pow(sqrt(discriminant)-R,1./3.); 


    double x1 = S1+S2;

    return 1-x1;

 }




  
  double nule(double Enu, double Ee, bool aneut = false ){
    
    double y = Y(Enu,Ee);
    
    double CLL;
    double CLR; 
    
    if( !aneut ) {
      CLL = -1./2.+sin2w;
      CLR = sin2w;
    }
    else {      
      CLL = sin2w;
      CLR = -1./2.+sin2w;
    }
    
    return Normalization()*S(Enu)*(CLL*CLL+CLR*CLR*(1-y)*(1-y)); 
    
  }

  double nuleInt(double Enu, bool aneut = false ){
    
    double norm = GF2/Pi*MeV2mbarn;
    
    double s = S(Enu); 
    
    double CLL;
    double CLR; 
     
   if( !aneut ) {
      CLL = -1./2.+sin2w;
      CLR = sin2w;
    }
    else {      
      CLL = sin2w;
      CLR = -1./2.+sin2w;
    }

    double ymax = (Enu-me)/Enu;
    
    return norm*s*(CLL*CLL*ymax+CLR*CLR*1./3.*(-(1-ymax)*(1.-ymax)*(1-ymax)+1.));
  }

  double GENnule(double Enu, double &El,double &cos, bool aneut = false ){
    
    double norm = GF2/Pi*MeV2mbarn;

    double s = S(Enu); 
    
    double CLL;
    double CLR; 
    
     if( !aneut ) {
      CLL = -1./2.+sin2w;
      CLR = sin2w;
    }
    else {      
      CLL = sin2w;
      CLR = -1./2.+sin2w;
    }
    
    double ymax = (Enu-me)/Enu;
    
    double N = nuleInt(Enu,aneut)/norm/s;

    double y1 = GenerateY(CLL,CLR,N); 
    
    El = y1*Enu+me;

    cos = GetCosine(Enu,El,neutrinoelectron);
    
    
    return  nule( Enu,  El, aneut );
    
  }

  double nueeInt(double Enu, bool aneut = false ){
        
    double CLL;
    double CLR; 
    
    if( !aneut ) {
      CLL = 1./2.+sin2w;
      CLR = sin2w;
    }
    else {
      CLL = -1./2.+sin2w;
      CLR = 1./2.+sin2w;  
    }

    double ymax = (Enu-me)/Enu;
    
    return Normalization()*S(Enu)*(CLL*CLL*ymax+CLR*CLR*1./3.*(-(1-ymax)*(1.-ymax)*(1-ymax)+1.)); 
    
  }

  
  double nuee(double Enu, double Ee, bool aneut = false ){
    
    double norm = GF2/Pi*MeV2mbarn;

    double s = S(Enu); 
    
    double y = Y(Enu,Ee); 
    
    double CLL;
    double CLR; 
    
    if( !aneut ) {
      CLL = 1./2.+sin2w;
      CLR = sin2w;
    }
    else {
      CLL = -1./2.+sin2w;
      CLR = 1./2.+sin2w;  
    }
    
    return norm*s*(CLL*CLL+CLR*CLR*(1-y)*(1-y)); 
    
  }

 double GENnuee(double Enu, double &El,double &cos, bool aneut = false ){
            
    double CLL;
    double CLR; 
    
    if( !aneut ) {
      CLL = 1./2.+sin2w;
      CLR = sin2w;
    }
    else {
      CLL = -1./2.+sin2w;
      CLR = 1./2.+sin2w;  
    }

    double ymax = (Enu-me)/Enu;
    
    double N = nueeInt(Enu,aneut)/Normalization()/S(Enu); 
    
    double y1 = GenerateY(CLL,CLR,N);  

    El = y1*Enu+me;

    cos = GetCosine(Enu,El,neutrinoelectron);
        
    return   nuee( Enu,  El, aneut );
    
  }
  
};




