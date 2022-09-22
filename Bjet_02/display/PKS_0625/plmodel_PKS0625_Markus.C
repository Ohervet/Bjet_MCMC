#include "math.h"
#include "Riostream.h"


//************************************************************************
//             Interpolation function
//         Interpolate sum of two emission fields f1 & f2
//
//	y1temp:		number of points f1
//	y2temp:		number of points f2
//	x1temp:		number of points f1 for discretization (free param)
//	x2temp:		number of points f2 for discretization (free param)
//	ftemp:		log10 flux density f1
//	f2temp:		log10 flux density f2
//	nutemp:		freq f1
//	nu2temp:	freq f2
//	ft:		output log10 flux density
//	nut:		output freq
//
//************************************************************************

int interpol(const int y1temp,const int y2temp,const int x1temp, const int x2temp, double ftemp[x1temp],double f2temp[x2temp],double nutemp[x1temp],double nu2temp[x2temp], double ft[x1temp+x2temp],double nut[x1temp+x2temp]){ 
  int i,j,a,k;
  double b;
  double y1, y2, x1, x2;
  double f[x1temp+x2temp],f2[x1temp+x2temp],nu[x1temp+x2temp],nu2[x1temp+x2temp];

  if( (y1temp==0) && (y2temp==0) ){
    return 0;
  }else if(y1temp == 0){
     for(int ii=0; ii<y2temp; ii++){
       ft[ii]=f2temp[ii];
       nut[ii]=nu2temp[ii];
    }
     return y2temp;
  }else if(y2temp == 0){
     for(int ii=0; ii<y1temp; ii++){
       ft[ii]=ftemp[ii];
       nut[ii]=nutemp[ii];
    }
     return y1temp;
  }
  // make sure that nu2[0] comes before nu[0]

  if(nu2temp[0]<=nutemp[0]){
    for(int ii=0; ii<x1temp; ii++){
      f[ii]=ftemp[ii];
      nu[ii]=nutemp[ii];
    }
    for(int ii=0; ii<x2temp; ii++){
      f2[ii]=f2temp[ii];
      nu2[ii]=nu2temp[ii];
    }
    x1=x1temp;
    y1=y1temp;
    x2=x2temp;
    y2=y2temp;
  }else{
    for(int ii=0; ii<x1temp; ii++){
      f2[ii]=ftemp[ii];
      nu2[ii]=nutemp[ii];
    }
    for(int ii=0; ii<x2temp; ii++){
      f[ii]=f2temp[ii];
      nu[ii]=nu2temp[ii];
    }
    x2=x1temp;
    y2=y1temp;
    x1=x2temp;
    y1=y2temp;
    
  }

  a=0;
       while(nu2[a]<nu[0] && a<x1+x2)
        { nut[a]=nu2[a];
          ft[a]=f2[a];
          a+=1; }
    j=a;
    i=0; 


    while(a<y2 && i<y1)
      { while(nu[i]<nu2[a] && i<y1)
	  { nut[j+i]=nu[i];
	    b= f2[a-1]+(nu[i]-nu2[a-1])*(f2[a]-f2[a-1])/(nu2[a]-nu2[a-1]);
	    ft[j+i]=log10(10**f[i]+10**b);
	  i+=1;}
	if(i<y1)a+=1;}   



    if(a == y2 && i<y1) 
      { while(i<y1)
      { nut[j+i]=nu[i];
	ft[j+i]=f[i];
	i+=1;}}

    else if(i == y1 && a<y2) 
      { while(a<y2)
      { nut[j+i]=nu2[a];
	ft[j+i]=f2[a];
	a+=1;
        j+=1;}}
        

    k=i+j;

    return k;
}




//*******************************************
//             MAIN ROUTINE
//*******************************************

 void plmodel()
{
  ifstream file0h,file0l,file1,file2,fileG;
  char   name;
  char   stmp[256];
  const int  nb0 =600, nb1 = 151, nb2 = 315, nbG = 1304, nbint = nb0+nbG;
  int i, k0 =0, k1 =0, k2 =0, k3=0,;
  Double_t freqh[nb0],ffluxh[nb0],freql[nb0],ffluxl[nb0],freq1[nb1],fflux1[nb1],freq2[nb2],fflux2[nb2];
  Double_t freqGa[nbG],fluxG[nbG],ffluxGa[nbG],freqG[nbG],ffluxG[nbG];
  Double_t freqtot1[nbint], fluxtot1[nbint];
  
  const double MeV=log10(1.6e-13),h=log10(6.626e-34),erg=-7,Hz2erg=6.6261e-27;

  //Initialisation

      gROOT->Reset();
      
      file0h.open("data/PKS_0625-354/Markus_hadron_4_high.dat");
      file0l.open("data/PKS_0625-354/Markus_hadron_4_low.dat");
      //file1.open("data/PKS_0625-354/Markus_e_synchrotron2.dat");
      //file2.open("data/PKS_0625-354/Markus_p_synchrotron2.dat");
      fileG.open("data/PKS_0625-354/PKS_0625_gal.dat"); // host galaxy


//Donnees simulees SED_Markus
//    nu [Hz],nu Fnu absorbed [Jy]
  if(!file0h){
      cout<<"SED file not found !"<<endl;     
    } 
    else{
      for(i=0;i<nb0;i++){ 
	file0h>>freqh[i]>>ffluxh[i];
	freqh[i] = log10(freqh[i]);
	ffluxh[i] = log10(ffluxh[i]/1.0e23);

	if (!file0h.good()) break;
      }
      k0=i;
    }
  if(!file0l){
      cout<<"SED file not found !"<<endl;     
    } 
    else{
      for(i=0;i<nb0;i++){ 
	file0l>>freql[i]>>ffluxl[i];
	freql[i] = log10(freql[i]);
	ffluxl[i] = log10(ffluxl[i]/1.0e23);

	if (!file0l.good()) break;
      }
      k0=i;
    }
  
  if(!file1){
      cout<<"SED file not found !"<<endl;     
    } 
    else{
      for(i=0;i<nb1;i++){ 
	file1>>freq1[i]>>fflux1[i];
	freq1[i] = log10(freq1[i]);
	fflux1[i] = log10(fflux1[i]/1.0e23);

	if (!file1.good()) break;
      }
      k1=i;
    }
    
  if(!file2){
      cout<<"SED file not found !"<<endl;     
    } 
    else{
      for(i=0;i<nb2;i++){ 
	file2>>freq2[i]>>fflux2[i];
	freq2[i] = log10(freq2[i]);
	fflux2[i] = log10(fflux2[i]/1.0e23);

	if (!file2.good()) break;
      }
      k2=i;
    }

//******************************************************************
// Host galaxy spectrum
//******************************************************************

      if(!fileG){
         cout<<"Host galaxy file not found !"<<endl;     
       } 
       else{
         for(i=nbG;i>0;i--){ 
            fileG>>freqGa[i]>>fluxG[i]>>ffluxGa[i]>>freqG[i]>>ffluxG[i];
            if (!fileG.good()) break;
         }
       }

//*****************************************************************
  //Data: read with the readdata.C script into TMultiGraph grdata
//******************************************************************

  std::cout<<"read data...."<<std::endl;
  
       wehavedata = gROOT->ProcessLine(".x readdata.C");
       
       
   std::cout<<"interpolate"<<std::endl;


    k3 = interpol(nb0,nbG,nbint,nbG,ffluxl,ffluxG,freql,freqG,fluxtot1,freqtot1); //host gal
      
  
    // 2. SED

   TMultiGraph *sed = new TMultiGraph();
   sed->SetTitle("");
   
   //SED_Markus
     if(k0>0){
      gr0l = new TGraph(k0,freql,ffluxl);
      gr0l->SetLineColor(kBlue-9);
      gr0l->SetLineStyle(3);
      gr0l->SetLineWidth(3);
      sed->Add(gr0l);    
      
      gr0h = new TGraph(k0,freqh,ffluxh);
      gr0h->SetLineColor(kRed);
      gr0h->SetLineStyle(3);
      gr0h->SetLineWidth(3);
          
    }
    
     if(k1>0){
      gr1 = new TGraph(k1,freq1,fflux1);
      gr1->SetLineColor(kBlue-9);
      gr1->SetLineStyle(1);
      gr1->SetLineWidth(3);
      sed->Add(gr1);    
    }
    
     if(k2>0){
      gr2 = new TGraph(k2,freq2,fflux2);
      gr2->SetLineColor(kRed-9);
      gr2->SetLineStyle(1);
      gr2->SetLineWidth(3);
      sed->Add(gr2);    
    }

    // Host galaxy
 
    if(nbG>0){
          grG = new TGraph(nbG,freqG,ffluxG);
	  grG->SetLineColor(1);
          grG->SetLineStyle(1);
          //sed->Add(grG);
    }
    
    //interpolation
    if(k3>0){
      //gr = new TGraph(k0,freq,fflux);
      gr = new TGraph(k3,freqtot1,fluxtot1);
      gr->SetLineColor(kBlue-9);
      gr->SetLineStyle(1);
      gr->SetLineWidth(3);
      sed->Add(gr);   
      sed->Add(gr0h);
    }

  
    TCanvas *csed = new TCanvas("SED","SED");
    csed->cd();
    // csed->SetBorderSize(0);   
    csed->GetFrame()->SetBorderSize(12);
    gPad->SetTopMargin(0.1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetTicky();


    Double_t minHz, maxHz,minF,maxF;
    minHz=7.5;//8
    maxHz=27.8;    
    minF=-14.;
    maxF=-9.5; //9.
    //synch zoom
    /*
    minHz=12;//8
    maxHz=19;    
    minF=-12.;
    maxF=-10; //9.
 */
    //  TH2D *hframe = new TH2D("hframe","",500,minHz,maxHz,500,minF,maxF);
    sed->Draw("AC");
    sed->GetXaxis()->SetRangeUser(minHz,maxHz);
    sed->GetXaxis()->SetLabelSize(0);
    sed->GetYaxis()->SetTitle("log (#nu F_{#nu} [erg cm^{-2} s^{-1} ] )");
    //sed->GetYaxis()->SetTitleFont(62);
    sed->GetYaxis()->SetTitleSize(0.04);
    sed->GetYaxis()->SetTitleOffset(1.5);
    sed->GetYaxis()->SetLabelSize(0.04);
    //sed->GetYaxis()->SetLabelFont(42);
    sed->SetMinimum(minF);
    sed->SetMaximum(maxF);

    Double_t mineV, maxeV;
    const Double_t hplanck=4.13566733E-15; // in eV s
    Double_t logh=log10(hplanck);
    mineV=minHz + logh;
    maxeV=maxHz + logh;

    // add energy axis
    TGaxis *eVaxis = new TGaxis(minHz,maxF,maxHz, maxF,mineV,maxeV,510,"-");
    eVaxis->SetTitle("log(E [eV])");
    eVaxis->SetLabelSize(0.04);
    eVaxis->SetTitleOffset(1.2);
    eVaxis->SetTitleFont(42);
    eVaxis->SetLabelFont(42);
    eVaxis->Draw();
    
    // add DATA
    if(wehavedata) grdata->Draw("P"); 
    sed->Draw("C"); // to have the fit above the data points
  

    
    //superpose the x axis over the legend
    TGaxis *Nuaxis = new TGaxis(minHz,minF,maxHz, minF,minHz,maxHz,510,"");
    Nuaxis->SetTitle("log (#nu [Hz])");
    Nuaxis->SetTitleOffset(1.2);
    Nuaxis->SetTitleFont(42);
    Nuaxis->SetLabelFont(42);
    Nuaxis->Draw();
}

