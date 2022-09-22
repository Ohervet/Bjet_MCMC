#include "math.h"
#include "Riostream.h"


//*******************************************
//             Interpolation functions
//*******************************************


int Quad_interpol(const int xlength, double f[xlength], double nu[xlength], const int xoutlength, double fout[xoutlength], double nuout[xoutlength]){
// transforme fction de taille xlength en xoutlength par interpolation quadratique
    int i, j, ii;
double step, a2, a1, b2, b1,c;
double AA, BB, CC;
Double_t fout1[xoutlength], fout2[xoutlength];

step = (nu[xlength-1] - nu[1])/xoutlength;
nuout[1] = nu[1];

for(int i=2; i<xoutlength; i++){
   nuout[i] =  nuout[i-1] + step;
}
ii = 0;
for(int j=2; j<xlength-1; j++){
    c= 1.0e5;
    //c = (nu[j+1] - nu[j+1]*nu[j+1]/nu[j-1]) / (nu[j] - nu[j]*nu[j]/nu[j-1]); //a reflechir pourquoi une grande valeure marche mieux...
   a2 = f[j+1] - (nu[j]/nu[j-1])*(nu[j]/nu[j-1])*f[j-1] - c*(f[j] - (nu[j]/nu[j-1])*(nu[j]/nu[j-1])*f[j-1]);
   b2 = 1.0 - (nu[j+1]/nu[j-1])*(nu[j+1]/nu[j-1]) - (1.0-(nu[j]/nu[j-1])*(nu[j]/nu[j-1])) * c;
   CC = a2/b2;
   //std::cout<<"a2 "<<a2<<std::endl;
   a1 = f[j] - (nu[j]/nu[j-1])*(nu[j]/nu[j-1])*f[j-1] - CC * (1.0 - (nu[j]/nu[j-1])*(nu[j]/nu[j-1]));
   b1 = nu[j] - nu[j]*nu[j]/nu[j-1];
   BB = a1/b1;
   
   
   AA = (f[j-1] - BB*nu[j-1] - CC) / (nu[j-1]*nu[j-1]);
   
    for(int ii=1; ii<xoutlength; ii++){
        if( (nuout[ii] >= nu[j-1]) && (nuout[ii] < nu[j]) ){
            fout1[ii] = AA*nuout[ii]*nuout[ii] + BB*nuout[ii] + CC;
            
        }
        if( (nuout[ii] >= nu[j]) && (nuout[ii] < nu[j+1]) ){
            fout2[ii] = AA*nuout[ii]*nuout[ii] + BB*nuout[ii] + CC;
        }
    }
}
for(int ii=1; ii<xoutlength; ii++){
    //fout[ii] = (fout1[ii] + 0.5 *fout2[ii])/1.5;
    fout[ii] = fout1[ii];
}
return 0;
}



int interpol(const int y1temp,const int y2temp,const int x1temp, const int x2temp, double ftemp[x1temp],double f2temp[x2temp],double nutemp[x1temp],double nu2temp[x2temp], double ft[x1temp+x2temp],double nut[x1temp+x2temp]){ 
  //y1temp: taille du domaine ou faire l'interpolation fction 1 -> from [0] to [y1temp]
  //y2temp: taille du domaine ou faire l'interpolation fction 2 -> from [0] to [y2temp]
  //x1temp: length fction 1
  //x1temp: length fction 2
  //nutemp: freq fction 1
  //nu2temp: freq fction 2
  int i,j,a,k;
  double b;
  double y1, y2, x1, x2;
  double f[x1temp+x2temp],f2[x1temp+x2temp],nu[x1temp+x2temp],nu2[x1temp+x2temp];

  //AZ 
  if( (y1temp==0) && (y2temp==0) ){ // 2 fonctions vides
    return 0;
  }else if(y1temp == 0){ // 1ere fonctoin vide
     for(int ii=0; ii<y2temp; ii++){
       ft[ii]=f2temp[ii];
       nut[ii]=nu2temp[ii];
    }
     return y2temp;
  }else if(y2temp == 0){ // 2eme fonction vide
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
      //std::cout<<"nu[ii] "<<nu[ii]<<std::endl;
    }
    x2=x1temp;
    y2=y1temp;
    x1=x2temp;
    y1=y2temp;
    
  }
  // end AZ
  

  a=0;
       while(nu2[a]<nu[0] && a<x1+x2)
        { nut[a]=nu2[a];
          ft[a]=f2[a];
       	//std::cout<<a<<" "<<nut[a]<<" "<<ft[a]<<" "<<std::endl;
          a+=1; }
    j=a;
    i=0; 

    //std::cout<<"etape 1 "<<y1<<" "<<y2<<std::endl;

    while(a<y2 && i<y1)
      { while(nu[i]<nu2[a] && i<y1)
	  { nut[j+i]=nu[i];
	    b= f2[a-1]+(nu[i]-nu2[a-1])*(f2[a]-f2[a-1])/(nu2[a]-nu2[a-1]);
            ft[j+i]=log10(pow(10,f[i])+pow(10,b));
	    //ft[j+i]=log10(10**f[i]+10**b);
	    //std::cout<<j+i<<" "<<a<<" "<<nut[j+i]<<" "<<ft[j+i]<<" "<<std::endl;
	  i+=1;}
	if(i<y1)a+=1;}   

    //std::cout<<"etape 2 "<<y1<<std::endl;


    if(a == y2 && i<y1) 
      { while(i<y1)
      { nut[j+i]=nu[i];
	ft[j+i]=f[i];
	//std::cout<<j+i<<" "<<nut[j+i]<<" "<<ft[j+i]<<" "<<std::endl;	
	i+=1;}}

    else if(i == y1 && a<y2) 
      { while(a<y2)
      { nut[j+i]=nu2[a];
	ft[j+i]=f2[a];
	//std::cout<<i+j<<" "<<a<<" "<<nut[j+i]<<" "<<ft[j+i]<<" "<<std::endl;	
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
  ifstream file1,file1p,file2,file2p,file3,file3p,file4,file5,file6,file7,file7p,file8,file101,file9,file9p,file9d,file9dp,file9a,file9ap,file10,file11,file11p,file12,file12p,file13,file13p,file14,file14p, FILE1;
  char   name;
  char   stmp[256];
  const int nb=500, ng = 1304, nda = 775,ndb = 129,n1 = nb, n3 = nb*3,n4 = nb*4,n5 = nb*5,n6 = nb*6,n7= nb*7,n2= n7+ng;
  int i,j,a,k1=0,k2=0,k3=0,k4=0,k4d=0,k5=0,k6=0,k7=0,k8=0,k9=0,k10=0,k11=0,k12=0,k13=0,k0=0,nb2=0,nb3=0,nb4=0,nb5=0,nb6=0,nb7=0,nb8=0,nb9=0,nb9a=0,nb9d=0,nb10=0,nb11=0,nb12=0,nb12a=0,nb12b=0,nb13=0,nb14=0,lx,l,wehavedata;
  int pk0=0,pnb2=0,pnb3=0,pnb7=0,pnb9=0,pnb9a=0,pnb9d=0,k101=0,pnb11=0,pnb12=0,pnb13=0,pnb14=0;
  Double_t freq[nb],flux[nb],ffluxc[nb],fluxc[nb],fflux[nb],
  freq2[nb],flux2[nb],fflux2[nb],freq3[nb],flux3[nb],fflux3[nb],
  freq101[nb],fluxc101[nb],ffluxc101[nb],flux101[nb],fflux101[nb],
  freq9[nb],fluxc9[nb],ffluxc9[nb],flux9[nb],fflux9[nb],
  freq9a[nb],fluxc9a[nb],ffluxc9a[nb],flux9a[nb],fflux9a[nb],
  freq9d[nb],fluxc9d[nb],ffluxc9d[nb],flux9d[nb],fflux9d[nb],
  freq10a[ng],freq10[ng],flux10[ng],fflux10a[ng],fflux10[ng],
  freq11[nb],flux11[nb],fflux11[nb],
  freq12[nb],flux12[nb],fflux12[nb],freq12a[nda],fflux12a[nda],freq12b[ndb],fflux12b[ndb];
  Double_t freq_highe[nb],fflux_highe[nb];
  Double_t pfreq[nb],pflux[nb],pffluxc[nb],pfluxc[nb],pfflux[nb],pfreq2[nb],pflux2[nb],pfflux2[nb],
  pfreq3[nb],pflux3[nb],pfflux3[nb],
  pfreq9[nb],pfluxc9[nb],pffluxc9[nb],pflux9[nb],pfflux9[nb],pfreq9a[nb],pfluxc9a[nb],pffluxc9a[nb],pflux9a[nb],pfflux9a[nb],pfreq9d[nb],pfluxc9d[nb],pffluxc9d[nb],pflux9d[nb],pfflux9d[nb],
  pfreq11[nb],pflux11[nb],pfflux11[nb],pfreq12[nb],pflux12[nb],pfflux12[nb];
  Double_t freqti[n1],fluxti[n1],freqt[n2],fluxt[n2],ge[nb],ne[nb],gse[nb],nse[nb],gp[nb],np[nb],eqm[3],freqinter[n2],ffluxinter[n2],tot[n1],freqtot[n1],tot2[n3],freqtot2[n3],
  tot3[n4],freqtot3[n4],tot3d[nb],freqtot3d[nb],tot4[n5],freqtot4[n5],tot5[n6],freqtot5[n6],
           tot6[n7],freqtot6[n7],tot7[n7],freqtot7[n7],tot8[n7],freqtot8[n7],tot9[n7],freqtot9[n7],tot10[n7],freqtot10[n7],tot11[n7],freqtot11[n7],tot12[n7],freqtot12[n7];
  Double_t A,B,b,k,s,km;
  Double_t freqext[n1],fluxext[n1],ffluxext[n1],logfreqext[n1],logfluxext[n1],logffluxext[n1],pfreqext[n1],pfluxext[n1],pffluxext[n1],plogfreqext[n1],plogfluxext[n1],plogffluxext[n1],
           freqjssc[n1],fluxjssc[n1],ffluxjssc[n1],logfreqjssc[n1],logfluxjssc[n1],logffluxjssc[n1],pfreqjssc[n1],pfluxjssc[n1],pffluxjssc[n1],plogfreqjssc[n1],plogfluxjssc[n1],plogffluxjssc[n1],
	   freqjeic[n1],fluxcjeic[n1],ffluxcjeic[n1],fluxjeic[n1],ffluxjeic[n1],pfreqjeic[n1],pfluxcjeic[n1],pffluxcjeic[n1],pfluxjeic[n1],pffluxjeic[n1];
  const double MeV=log10(1.6e-13),h=log10(6.626e-34),erg=-7,Hz2erg=6.6261e-27;

  //Initialisation

      gROOT->Reset();

      file1.open("data/test_bj_cs.dat"); // e ssc 
      file1p.open("data/test_bj_prev_cs.dat"); // previous e ssc file
      file2.open("data/test_bj_ss.dat"); // e synch
      file2p.open("data/test_bj_prev_ss.dat"); // previous e synch file
      file3.open("data/test_bj_pss.dat");  // p synch
      file3p.open("data/test_bj_prev_pss.dat"); // previous p synch file
      file4.open("data/test_bj_es.dat"); // e spectrum
      file5.open("data/test_bj_ps.dat");  // p spectrum
      file6.open("data/test_bj_ses.dat"); // secondary e spectrum
      file7.open("data/F_jet_syn.dat"); // e synch from the extended jet
      
      //file7p.open("data/F_jet_syn_prev.dat"); // previous extended jet
      file7p.open("data/OJ_287/F_jet_syn_core.dat");
      
      
      file9.open("data/test_bj_ecs.dat"); // EIC disc-BLR file
      file9p.open("data/test_bj_prev_ecs.dat"); // previous EIC disc-BLR file
      file9d.open("data/test_bj_ecdisc.dat"); // EIC discfile
      //file9dp.open("data/test_bj_prev_ecdisc.dat"); // previous EIC disc-BLR file
      file9a.open("data/test_bj_ecs1.dat"); // EIC torus file
      file9ap.open("data/test_bj_prev_ecs1.dat"); // previous EIC torus file
      file101.open("data/test_bj_cs2.dat"); // 2n order SSC
      //file10.open("data/HESS_J1943+213/HGS_15_cut_z02.dat"); // host galaxy
      file10.open("data/HESS_J1943+213/HGS_15_cut_z016.dat"); // host galaxy
      file11.open("data/test_bj_nuc.dat"); // disk spectrum
      file11p.open("data/test_bj_prev_nuc.dat"); // previous disk spectrum
      file12.open("data/test_bj_tor.dat"); // torus spectrum
      file12p.open("data/test_bj_prev_tor.dat"); // previous torus spectrum
      file13.open("data/F_jet_com.dat"); // e ssc from the extended jet
      
      //file13p.open("data/F_jet_com_prev.dat"); // previous e ssc from the extended jet
      file13p.open("data/OJ_287/F_jet_com_core.dat");

      
      file14.open("data/test_bj_ecs_jet.dat"); // eic from the extended jet
      file14p.open("data/test_bj_prev_ecs_jet.dat"); // previous eic from the extended jet

//*******************************************************************************
// SBLOB CURVES
//*******************************************************************************

  //Donnees simulees Compton inverse
  //    nu [log Hz], F absorbed, nu Fnu absorbed [J cm^-2 s^-1 ?], F, nu Fnu

      if(!file1){
         cout<<"SSC file not found !"<<endl;     
       } 
       else{
         for(i=0;i<nb;i++){ 
            file1>>freq[i]>>fluxc[i]>>ffluxc[i]>>flux[i]>>fflux[i];

	    //        std::cout<<i<<" "<<freq[i]<<" "<<flux[i]<<" "<<fflux[i]<<" "<<fluxc[i]<<" "<<ffluxc[i]<<" "<<std::endl;

            if (!file1.good()) break;
         }
         k0=i;
       }

      if(!file1p){
         cout<<"previous SSC file not found !"<<endl;     
       } 
       else{
         for(i=0;i<nb;i++){ 
	    file1p>>pfreq[i]>>pfluxc[i]>>pffluxc[i]>>pflux[i]>>pfflux[i];
            if (!file1p.good()) break;
         }
         pk0=i;
       }

  //Donnees simulees 2nd order SSC
  //    nu [log Hz], F absorbed, nu Fnu absorbed [J cm^-2 s^-1 ?], F, nu Fnu

      if(!file101){
         cout<<"2nd order SSC file not found !"<<endl;     
       } 
       else{
         for(i=0;i<nb;i++){ 
            file101>>freq101[i]>>fluxc101[i]>>ffluxc101[i]>>flux101[i]>>fflux101[i];

	    //        std::cout<<i<<" "<<freq[i]<<" "<<flux[i]<<" "<<fflux[i]<<" "<<fluxc[i]<<" "<<ffluxc[i]<<" "<<std::endl;

            if (!file101.good()) break;
         }
         k101=i;
       }



   // Donnees simulees synchrotron des electrons
	    
       if(!file2){
             cout<<"e synch file not found !"<<endl;     
       }
       else{
         for(i=0;i<nb;i++){ 
            file2>>freq2[i]>>flux2[i]>>fflux2[i];
            if (!file2.good()) break;
            //std::cout<<i<<" "<<freq2[i]<<" "<<flux2[i]<<" "<<fflux2[i]<<std::endl;
            nb2=i;
         }
       }

       if(!file2p){
             cout<<"previous e synch file not found !"<<endl;     
       }
       else{
         for(i=0;i<nb;i++){ 
            file2p>>pfreq2[i]>>pflux2[i]>>pfflux2[i];
 
            if (!file2p.good()) break;
         }
	pnb2=i;
       }

  // Donnees simulees synchrotron des protons
	    
       if(!file3){
            cout<<" p synch file not found !"<<endl;     
       }
       else{
         for(i=0;i<nb;i++){ 
           file3>>freq3[i]>>flux3[i]>>fflux3[i];
           if (!file3.good()) break;
         }
	nb3=i;
       }

       if(!file3p){
            cout<<"previous p synch file not found !"<<endl;     
       }
       else{
         for(i=0;i<nb;i++){ 
           file3p>>pfreq3[i]>>pflux3[i]>>pfflux3[i];
           if (!file3p.good()) break;
         }
	pnb3=i;
       }


   // Spectre des electrons
	    
       if(!file4){
          cout<<"e spectrum file not found !"<<endl;     
       }
       else{
         for(i=0;i<nb;i++){ 
            file4>>ge[i]>>ne[i]>>A>>B;
            if (!file4.good()) break;
         }
	 nb4=i;
       }

 // Spectre des protons
	    
       if(!file5){
          cout<<"proton spectrum not found !"<<endl; 
       }
       else{    
         for(i=0;i<nb;i++){ 
            file5>>gp[i]>>np[i]>>A>>B;
            if (!file5.good()) break;
         }
	 nb5=i;
       }

  // Spectre des electrons secondaires
	    
       if(!file6){
          cout<<"secondary e file not found !"<<endl;     
       }
       else{
         for(i=0;i<nb;i++){ 
            file6>>gse[i]>>nse[i]>>A>>B;
            if (!file6.good()) break;
         }
	 nb6=i;
       }

 // Sum of e-synch and SSC (absorbed)

       if(!file8){
	 cout<<"interpolated e-synch / SSC file not found !"<<endl;
       }
       else{
         for(i=0;i<n2;i++){
	   file8>>freqinter[i]>>ffluxinter[i];
           if(!file8.good()) break;
         }
         nb8=i;
       }
       
       
//******************************************************************
// EIC disc & disc-BLR spectra 
//******************************************************************

      if(!file9){
         cout<<"EIC disc-BLR file not found !"<<endl;     
       } 
       else{
         for(i=0;i<nb;i++){ 
            file9>>freq9[i]>>fluxc9[i]>>ffluxc9[i]>>flux9[i]>>fflux9[i];

            if (!file9.good()) break;
         }
         nb9=i;
       }

      if(!file9p){
         cout<<"previous EIC disc-BLR file not found !"<<endl;     
       } 
       else{
         for(i=0;i<nb;i++){ 
	    file9p>>pfreq9[i]>>pfluxc9[i]>>pffluxc9[i]>>pflux9[i]>>pfflux9[i];
            if (!file9p.good()) break;
         }
         pnb9=i;
       }
       
      if(!file9d){
         cout<<"EIC disc file not found !"<<endl;     
       } 
       else{
         for(i=0;i<nb;i++){ 
            file9d>>freq9d[i]>>fluxc9d[i]>>ffluxc9d[i]>>flux9d[i]>>fflux9d[i];

            if (!file9d.good()) break;
         }
         nb9d=i;
       }

      if(!file9dp){
         cout<<"previous EIC disc file not found !"<<endl;     
       } 
       else{
         for(i=0;i<nb;i++){ 
	    file9dp>>pfreq9d[i]>>pfluxc9d[i]>>pffluxc9d[i]>>pflux9d[i]>>pfflux9d[i];
            if (!file9dp.good()) break;
         }
         pnb9d=i;
       }
       
// EIC torus spectrum 

      if(!file9a){
         cout<<"EIC torus file not found !"<<endl;     
       } 
       else{
         for(i=0;i<nb;i++){ 
            file9a>>freq9a[i]>>fluxc9a[i]>>ffluxc9a[i]>>flux9a[i]>>fflux9a[i];

            if (!file9a.good()) break;
         }
         nb9a=i;
       }

      if(!file9ap){
         cout<<"previous EIC torusfile not found !"<<endl;     
       } 
       else{
         for(i=0;i<nb;i++){ 
	    file9ap>>pfreq9a[i]>>pfluxc9a[i]>>pffluxc9a[i]>>pflux9a[i]>>pfflux9a[i];
            if (!file9ap.good()) break;
         }
         pnb9a=i;
       }
       
// Disk blackbody spectrum 

      if(!file11){
         cout<<"disk file not found !"<<endl;     
       } 
       else{
         for(i=0;i<nb;i++){ 
            file11>>freq11[i]>>flux11[i]>>fflux11[i];

            if (!file11.good()) break;
         }
         nb11=i;
       }

      if(!file11p){
         cout<<"previous disk file not found !"<<endl;     
       } 
       else{
         for(i=0;i<nb;i++){ 
	    file11p>>pfreq11[i]>>pflux11[i]>>pfflux11[i];
            if (!file11p.good()) break;
         }
         pnb11=i;
       }

// Dust torus blackbody spectrum 

      if(!file12){
         cout<<"torus file not found !"<<endl;     
       } 
       else{
         for(i=0;i<nb;i++){ 
            file12>>freq12[i]>>flux12[i]>>fflux12[i];

            if (!file12.good()) break;
         }
         nb12=i;
       }

      if(!file12p){
         cout<<"previous torus file not found !"<<endl;     
       } 
       else{
         for(i=0;i<nb;i++){ 
	    file12p>>pfreq12[i]>>pflux12[i]>>pfflux12[i];
            if (!file12p.good()) break;
         }
         pnb12=i;
       }
//******************************************************************
// SimpleJet curve
//******************************************************************

       if(!file7){
	 cout<<"no extended synchrotron jet file found !"<<endl;
       }
       else{
         for(i=0;i<n1;i++){ 
	    file7>>freqext[i]>>fluxext[i]>>ffluxext[i]>>logfreqext[i]>>logfluxext[i]>>logffluxext[i];
            if (!file7.good()) break;
            nb7=i;
            //std::cout<<logffluxext[i]<<std::endl;
          }
       } 
       
       if(!file7p){
         cout<<"previous extended synchrotron jet file not found !"<<endl;     
       } 
       else{
         for(i=0;i<n1;i++){ 
	    file7p>>pfreqext[i]>>pfluxext[i]>>pffluxext[i]>>plogfreqext[i]>>plogfluxext[i]>>plogffluxext[i];
            if (!file7p.good()) break;
         }
         pnb7=i;
       }
       
       if(!file13){
	 cout<<"no extended ssc jet file found !"<<endl;
       }
       else{
         for(i=0;i<n1;i++){ 
	    file13>>freqjssc[i]>>fluxjssc[i]>>ffluxjssc[i]>>logfreqjssc[i]>>logfluxjssc[i]>>logffluxjssc[i];
          
            if (!file13.good()) break;
          }
         nb13=i;
       } 
       
       if(!file13p){
         cout<<"previous extended ssc jet file not found !"<<endl;     
       } 
       else{
         for(i=0;i<n1;i++){ 
	    file13p>>pfreqjssc[i]>>pfluxjssc[i]>>pffluxjssc[i]>>plogfreqjssc[i]>>plogfluxjssc[i]>>plogffluxjssc[i];
            if (!file13p.good()) break;
         }
         pnb13=i;
       }
       
              if(!file14){
	 cout<<"no eic jet file found !"<<endl;
       }
       else{
         for(i=0;i<n1;i++){ 
	    file14>>freqjeic[i]>>fluxcjeic[i]>>ffluxcjeic[i]>>fluxjeic[i]>>ffluxjeic[i];
          
            if (!file14.good()) break;
          }
         nb14=i;
       } 
       
       if(!file14p){
         cout<<"previous eic jet file not found !"<<endl;     
       } 
       else{
         for(i=0;i<n1;i++){ 
	    file14p>>pfreqjeic[i]>>pfluxcjeic[i]>>pffluxcjeic[i]>>pfluxjeic[i]>>pffluxjeic[i];
            if (!file14p.good()) break;
         }
         pnb14=i;
       }

//******************************************************************
// Host galaxy spectrum
//******************************************************************

      if(!file10){
         cout<<"Host galaxy file not found !"<<endl;     
       } 
       else{
         for(i=ng-1;i>0;i--){ 
            file10>>freq10a[i]>>flux10[i]>>fflux10a[i]>>freq10[i]>>fflux10[i];
            if (!file10.good()) break;
         }
         nb10=ng;
       }


//*****************************************************************
  //Data: read with the readdata.C script into TMultiGraph grdata
//******************************************************************

  std::cout<<"read data...."<<std::endl;
  
       wehavedata = gROOT->ProcessLine(".x readdata.C");
       


   
   std::cout<<"interpolate"<<std::endl;

    k2 = interpol(nb7,nb2,500,nb, logffluxext,fflux2,logfreqext,freq2,tot,freqtot); //Jet syn + blob syn
    k3 = interpol(k2,k0,500,nb, tot,ffluxc,freqtot,freq,tot2,freqtot2); //+ blob ssc
    k4d = interpol(k3,nb9d,500,nb,tot2,fflux9d,freqtot2,freq9d,tot3d,freqtot3d); // +EIC disc
    k4 = interpol(k4d,nb9,n4,nb,tot3d,fflux9,freqtot3d,freq9,tot3,freqtot3); // +EIC disc-BLR
    k5 = interpol(k4,nb6,n5,nb,tot3,gse,freqtot3,nse,tot4,freqtot4);
    k6 = interpol(k5,nb11,n6,nb,tot4,fflux11,freqtot4,freq11,tot5,freqtot5);//disk
    k7 = interpol(k6,k101,n7,nb, tot5,fflux101,freqtot5,freq101,tot6,freqtot6); //2nd order SSC
    //Double_t tot7b[500], freqtot7b[500];
    //Quad_interpol(k8, tot7, freqtot7, 500, tot7b, freqtot7b);
    k9 = interpol(k7,nb13,n7,nb,tot6,logffluxjssc,freqtot6,logfreqjssc,tot8,freqtot8); // ssc jet etendu
    //k9 = interpol(500,nb13,498,nb,tot7b,logffluxjssc,freqtot7b,logfreqjssc,tot8,freqtot8); // ssc jet etendu
    k10 = interpol(k9,nb14,n7,nb,tot8,ffluxjeic,freqtot8,freqjeic,tot9,freqtot9); // eic jet etendu
    /*
    k11 = interpol(k10,nb12,n7,nb,tot9,fflux12,freqtot9,freq12,tot10,freqtot10); //ajout tore
    k12 = interpol(k11,nb9a,n7,nb,tot10,fflux9a,freqtot10,freq9a,tot11,freqtot11); //ajout eic du tore
    */
    //for Chi2 calculation
    Double_t tot11b[500], freqtot11b[500];
    Quad_interpol(k10, tot9, freqtot9, 500, tot11b, freqtot11b);
    /*
    k13 = interpol(498,nb10,498,nb10,tot11b,fflux10,freqtot11b,freq10,tot12,freqtot12); //host gal
    
    */

    
    
// Chi2
    
   // Chi2 of the fit on Fermi and VERITAS

   double Chi2_Fermi = 0.;
   Double_t ymodel_Fermi[n_F];
   double Chi2_VERITAS = 0.;
   Double_t ymodel_VERITAS[n_V];
   double F;
   double deltaFmin;
   double deltaFmax;
   double ymodel;
   for (i=0;i<(500);i++){
       
       //std::cout<<"freqtot11b[i] = "<<freqtot11b[i]<<std::endl;
       /*
        for(j=0;j<nx;j++){
           if((freqtot11b[i]<=E_xrt[j]) && (freqtot11b[i+1]>E_xrt[j])){
               //linear interpolation
               ymodel_XRT[j] = (tot11b[i+1] - (tot11b[i+1]-tot11b[i])*(freqtot11b[i+1]-E_xrt[j])/(freqtot11b[i+1]-freqtot11b[i]));
               if (F_xrt[j] >= ymodel_XRT[j]){
                  Chi2_XRT += pow((F_xrt[j]-ymodel_XRT[j]) / deltaFmin_xrt[j],2);
               }else{
                  Chi2_XRT += pow((F_xrt[j]-ymodel_XRT[j]) / deltaFmax_xrt[j],2);
               }
           }
       }
        for(j=0;j<7;j++){
           if((freqtot11b[i]<=E_bat[j]) && (freqtot11b[i+1]>E_bat[j])){
               //linear interpolation
               ymodel_BAT[j] = (tot11b[i+1] - (tot11b[i+1]-tot11b[i])*(freqtot11b[i+1]-E_bat[j])/(freqtot11b[i+1]-freqtot11b[i]));
               if (F_bat[j] >= ymodel_BAT[j]){
                  Chi2_BAT += pow((F_bat[j]-ymodel_BAT[j]) / deltaFmin_bat[j],2);
               }else{
                  Chi2_BAT += pow((F_bat[j]-ymodel_BAT[j]) / deltaFmax_bat[j],2);
               }
           }
       }
       */
        for(j=0;j<(n_F-2);j++){
            F = pow(10,F_F[j]);
            deltaFmin = F - pow(10,F_F[j]-deltaFmin_F[j]);
            deltaFmax = pow(10,F_F[j]+deltaFmax_F[j]) - F;
           if((freqtot11b[i]<=E_F[j]) && (freqtot11b[i+1]>E_F[j])){
               //linear interpolation
               ymodel_Fermi[j] = (tot11b[i+1] - (tot11b[i+1]-tot11b[i])*(freqtot11b[i+1]-E_F[j])/(freqtot11b[i+1]-freqtot11b[i]));
               ymodel = pow(10,ymodel_Fermi[j]);
               if (F_F[j] >= ymodel_Fermi[j]){
                  Chi2_Fermi += pow((F-ymodel) / deltaFmin,2);
               }else{
                  Chi2_Fermi += pow((F-ymodel) / deltaFmax,2);
               }
           }
       }
        for(j=1;j<(n_V-3);j++){
            F = pow(10,F_V[j]);
            deltaFmin = F - pow(10,F_V[j]-deltaFmin_V[j]);
            deltaFmax = pow(10,F_V[j]+deltaFmax_V[j]) - F;
           if((freqtot11b[i]<=E_V[j]) && (freqtot11b[i+1]>E_V[j])){
               //linear interpolation
               ymodel_VERITAS[j] = (tot11b[i+1] - (tot11b[i+1]-tot11b[i])*(freqtot11b[i+1]-E_V[j])/(freqtot11b[i+1]-freqtot11b[i]));
               ymodel = pow(10,ymodel_VERITAS[j]);
               if (F_V[j] >= ymodel_VERITAS[j]){
                  Chi2_VERITAS += pow((F-ymodel) / deltaFmin,2);
               }else{
                  Chi2_VERITAS += pow((F-ymodel) / deltaFmax,2);
               }
               std::cout<<"pow((F-ymodel) / deltaFmin,2) = "<<pow((F-ymodel) / deltaFmin,2)<<std::endl;
           }
       }
   }
   std::cout<<"Chi2_Fermi/DOF = "<<Chi2_Fermi<<"/"<<5<<std::endl;
   std::cout<<"Chi2_VERITAS/DOF = "<<Chi2_VERITAS<<"/"<<3<<std::endl;
   std::cout<<"Chi2_tot = "<<Chi2_VERITAS+Chi2_Fermi<<std::endl;



    
// Make Plots

 // 1. Particle spectra

   TMultiGraph *graph1 = new TMultiGraph();  
   graph1->SetTitle("Particle spectra");
  
   if(nb4>0){
          TGraph *s1 = new TGraph(nb4-1,ge,ne);
	  s1->SetLineColor(4);
          graph1->Add(s1);
    }
    
    if(nb5>0){
          TGraph *s2 = new TGraph(nb5-1,gp,np);
	  s2->SetLineColor(2);
          graph1->Add(s2);
    }

    if(nb6>0){
          TGraph *s3 = new TGraph(nb6,gse,nse);
	  s3->SetLineColor(3);
	  s3->SetLineStyle(2);
          graph1->Add(s3);
    }

    TCanvas *inputspec = new TCanvas("inputspec","inputspec");
    inputspec->cd();
    graph1->Draw("AC");
   

    // 2. SED

   TMultiGraph *sed = new TMultiGraph();
   sed->SetTitle("");

   // SSC
    if(pk0>0){
      TGraph *pgr = new TGraph(pk0,pfreq,pffluxc);  // previous SSC absorbed
      pgr->SetLineColor(1);
      pgr->SetLineStyle(3);
      sed->Add(pgr);
    }
    
    if(k0>0){
        TGraph *gr = new TGraph(k0,freq,ffluxc); // SSC absorbed
	gr->SetLineColor(kBlue-9);
        gr->SetLineStyle(1);
	gr->SetLineWidth(2);
        sed->Add(gr);

	//**** SSC unabsorbed: remove low frequency points for plot
        int i_highe=0;
        for(int i=0; i<k0; i++){
	   if(freq[i]>22.){
              freq_highe[i_highe]=freq[i];
              fflux_highe[i_highe]=fflux[i];
              i_highe++;
           }
        }
	//****
	//gr2 = new TGraph(k0,freq,fflux); // SSC unabsorbed  (plot all points)
        TGraph *gr2 = new TGraph(i_highe, freq_highe, fflux_highe); //plot only high energy points
        gr2->SetLineColor(kBlue-9);
        gr2->SetLineStyle(3);
	gr2->SetLineWidth(2);
        sed->Add(gr2);
    }

    if(k101>0){
        TGraph *gr101 = new TGraph(k101,freq101,ffluxc101); // SSC 2nd absorbed
	gr101->SetLineColor(kBlue-9);
        gr101->SetLineStyle(2);
	gr101->SetLineWidth(2);
        sed->Add(gr101);
    }

    // sum of e-synch and absorbed SSC

    if(nb8>0){
      TGraph *gr_inter = new TGraph(nb8,freqinter,ffluxinter);
      gr_inter->SetLineColor(4);
      gr_inter->SetLineStyle(1);
      //sed->Add(gr_inter);
    }

    //   e synch 
    if(pnb2>0){
          TGraph *pgr3 = new TGraph(pnb2,pfreq2,pfflux2);
	  pgr3->SetLineColor(1);
          pgr3->SetLineStyle(3);
	  pgr3->SetLineWidth(2);
          //sed->Add(pgr3);
    }
    
    if(nb2>0){
          TGraph *gr3 = new TGraph(nb2,freq2,fflux2);
	  gr3->SetLineColor(kBlue-9);
          gr3->SetLineStyle(1);
	  gr3->SetLineWidth(2);
          sed->Add(gr3);
    }

    // p synch

    if(nb3>0){
          TGraph *gr4 = new TGraph(nb3,freq3,fflux3);
	  gr4->SetLineColor(2);
          gr4->SetLineStyle(1);
          //sed->Add(gr4);
    }
 
    if(pnb3>0){
          TGraph *pgr4 = new TGraph(nb3,freq3,fflux3);
	  pgr4->SetLineColor(6);
          pgr4->SetLineStyle(3);
          //sed->Add(pgr4);
    }

    //  extended jet synch

    if(nb7>0){
        TGraph *gr7 = new TGraph(nb7,logfreqext,logffluxext);
	gr7->SetLineColor(kRed-9);
        gr7->SetLineStyle(5);
	gr7->SetLineWidth(4);
        sed->Add(gr7);
    }
    
    if(pnb7>0){
        TGraph *pgr7 = new TGraph(pnb7,plogfreqext,plogffluxext);
	pgr7->SetLineColor(kRed-9);
        //pgr7->SetLineStyle(8);
        sed->Add(pgr7);
    }
    
    //  extended jet ssc

    if(nb13>0){
        TGraph *gr13 = new TGraph(nb13,logfreqjssc,logffluxjssc);
	gr13->SetLineColor(kRed-9);
        gr13->SetLineStyle(5);
	gr13->SetLineWidth(4);
        sed->Add(gr13);
    }
    
    if(pnb13>0){
        TGraph *pgr13 = new TGraph(pnb13,plogfreqjssc,plogffluxjssc);
	pgr13->SetLineColor(kRed-9);
        //pgr13->SetLineStyle(8);
        sed->Add(pgr13);
    }
    
    //  EIC from extended jet

    if(nb14>0){
        TGraph *gr14 = new TGraph(nb14,freqjeic,ffluxjeic);
	gr14->SetLineColor(kMagenta-9);
        gr14->SetLineStyle(7);
	gr14->SetLineWidth(2);
        sed->Add(gr14);
    }
    
    if(pnb14>0){
        TGraph *pgr14 = new TGraph(pnb14,pfreqjeic,pffluxjeic);
	pgr14->SetLineColor(kMagenta-6);
        pgr14->SetLineStyle(3);
        //sed->Add(pgr14);
    }
    
    // EIC from disc & disc-BLR
    
    if(nb9>0){
          TGraph *gr9 = new TGraph(nb9,freq9,fflux9);
	  gr9->SetLineColor(kGreen-6);
	  //gr9->SetLineColorAlpha(2, 0.35)
          gr9->SetLineStyle(2);
	  gr9->SetLineWidth(2);
          sed->Add(gr9);
    }
    
    if(pnb9>0){
          TGraph *pgr9 = new TGraph(pnb9,pfreq9,pfflux9);
	  pgr9->SetLineColor(kGreen-9);
          pgr9->SetLineStyle(3);
          //sed->Add(pgr9);
    }
    
    if(nb9d>0){
          TGraph *gr9d = new TGraph(nb9d,freq9d,fflux9d);
	  gr9d->SetLineColor(kGreen-6);
	  //gr9->SetLineColorAlpha(2, 0.35)
          gr9d->SetLineStyle(2);
	  gr9d->SetLineWidth(2);
          sed->Add(gr9d);
    }
    
    if(pnb9d>0){
          TGraph *pgr9d = new TGraph(pnb9d,pfreq9d,pfflux9d);
	  pgr9d->SetLineColor(kGreen-9);
          pgr9d->SetLineStyle(3);
          //sed->Add(pgr9d);
    }   
   // EIC from dust torus
    
    if(nb9a>0){
          TGraph *gr9a = new TGraph(nb9a,freq9a,fflux9a);
	  gr9a->SetLineColor(2);
          gr9a->SetLineStyle(10);
	  gr9a->SetLineWidth(2);
          //sed->Add(gr9a);
    }
    
    if(pnb9a>0){
          TGraph *pgr9a = new TGraph(nb9a,freq9a,fflux9a);
	  pgr9a->SetLineColor(2);
          pgr9a->SetLineStyle(3);
          //sed->Add(pgr9a);
    }
    
    // disk blackbody 
    
    if(nb11>0){
          TGraph *gr11 = new TGraph(nb11,freq11,fflux11);
	  //Int_t trans_orange = GetColorTransparent(kOrange+1, 0.3);
	  //gr11->SetLineColorAlpha(2, 0.35)
	  gr11->SetLineColor(kGreen-6);
          gr11->SetLineStyle(8);
	  gr11->SetLineWidth(4);
          sed->Add(gr11);
    }
    
    if(pnb11>0){
          TGraph *pgr11 = new TGraph(nb11,pfreq11,pfflux11);
	  pgr11->SetLineColor(kGreen-6);
          pgr11->SetLineStyle(6);
          //sed->Add(pgr11);
    }
    
    // dust torus blackbody 
    
    if(nb12>0){
          TGraph *gr12 = new TGraph(nb12,freq12,fflux12);
	  gr12->SetLineColor(1);
          gr12->SetLineStyle(9);
	  gr12->SetLineWidth(2);
          sed->Add(gr12);
    }
    
    if(pnb12>0){
          TGraph *pgr12 = new TGraph(pnb12,pfreq12,pfflux12);
	  pgr12->SetLineColor(1);
          pgr12->SetLineStyle(3);
          //sed->Add(pgr12);
    }

    // Host galaxy
    
    if(nb10>0){/*
          gr10 = new TGraph(nb10,freq10,fflux10);
	  gr10->SetLineColor(2);
          gr10->SetLineStyle(1);*/
          //sed->Add(gr10);
    }
    
    // Accretion disk template
    
    if(nb12a>0){
          TGraph *gr12a = new TGraph(nb12a,freq12a,fflux12a);
	  gr12a->SetLineColor(8);
          gr12a->SetLineStyle(1);
          //sed->Add(gr12a);
    }
    if(nb12b>0){
          TGraph *gr12b = new TGraph(nb12b,freq12b,fflux12b);
	  gr12b->SetLineColor(8);
          gr12b->SetLineStyle(1);
          //sed->Add(gr12b);
    }
    //interpolation
	  TGraph *grtot = new TGraph(k10-1,freqtot9,tot9);
	  grtot->SetLineColor(kOrange+4);//kOrange+1
          grtot->SetLineStyle(1);
	  grtot->SetLineWidth(2);
          sed->Add(grtot);
	  
	  
    // ecriture de l'enveloppe dans un fichier
    FILE*  stream_dat;
    
    stream_dat = fopen("data/model_bj.dat", "w+");
  
    for (i = 2; i <= k10-1; i++) {
      fprintf(stream_dat, "%f %f \n", freqtot9[i],  tot9[i]);	
    }
/*
    for (i = 2; i <= 498; i++) {
      //fprintf(stream_dat, "%f %f \n", freqtot13[i],  tot13[i]);	
      std::cout<<freqtotb[i]<<totb[i]<<std::endl; 
    }
    */
    fclose(stream_dat);


    TCanvas *csed = new TCanvas("SED","SED");
    csed->cd();
    // csed->SetBorderSize(0);   
    csed->GetFrame()->SetBorderSize(12);
    gPad->SetTopMargin(0.1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetTicky();


    Double_t minHz, maxHz,minF,maxF;
    minHz=8;
    maxHz=28;    
    minF=-14;
    maxF=-9; 
 
    //TH2D *hframe = new TH2D("hframe","",500,minHz,maxHz,500,minF,maxF);
    sed->Draw("AC");
    //sed->GetXaxis()->SetRangeUser(minHz,maxHz);
    sed->GetXaxis()->SetLimits(minHz,maxHz);
    sed->GetXaxis()->SetLabelSize(0);
    sed->GetYaxis()->SetTitle("log (#nu F#nu [erg cm^{-2} s^{-1} ] )");
    sed->GetYaxis()->SetTitleSize(0.04);
    sed->GetYaxis()->SetTitleOffset(1.4);
    sed->GetYaxis()->SetLabelSize(0.04);
   // sed->GetXaxis()->SetTitle("log (#nu [Hz])");
    //sed->GetXaxis()->SetTitleSize(0.04);
    //sed->GetXaxis()->SetTitleOffset(1.3);
    sed->SetMinimum(minF);
    sed->SetMaximum(maxF);

    Double_t mineV, maxeV;
    const Double_t hplanck=4.13566733E-15; // in eV s
    Double_t logh=log10(hplanck);
    mineV=minHz + logh;
    maxeV=maxHz + logh;

    


    
    // add DATA
    
    if(wehavedata) grdata->Draw("P"); 
    sed->Draw("C"); // to have the fit above the data points
    
    
    
    TLegend *leg = new TLegend(0.15,0.78,0.92,0.87);
    leg->SetFillColor(kWhite);
    leg->SetNColumns(4);
    leg->SetMargin(0.15);

    leg->AddEntry(gr,"Blob Syn. & SSC","l");
    leg->AddEntry(pgr13,"Jet Syn. & SSC, 16pc","l");
    leg->AddEntry(gr14,"EIC blob - jet","l");
    leg->AddEntry(grtot,"All components","l");
    leg->AddEntry(gr11,"Disk thermal emission","l");
    leg->AddEntry(gr13,"Jet Syn. & SSC, 50pc","l");
    leg->AddEntry(gr9,"EIC blob - BLR","l");
    
    
    leg->SetEntrySeparation(0.0);
    leg->SetTextSize(0.03);
    leg->Draw();
    

    
    
    
    // add energy axis
    TGaxis *eVaxis = new TGaxis(minHz,maxF,maxHz, maxF,mineV,maxeV,510,"-");
    eVaxis->SetTitle("log(E [eV])");
    eVaxis->SetTitleOffset(1.15);
    eVaxis->SetTitleFont(42);
    eVaxis->SetLabelFont(42);
    eVaxis->Draw();
    
    
    
    //superpose the x axis over the legend
    TGaxis *Nuaxis = new TGaxis(minHz,minF,maxHz, minF,minHz,maxHz,510,"");
    Nuaxis->SetTitle("log (#nu [Hz])");
    Nuaxis->SetTitleOffset(1.15);
    Nuaxis->SetTitleFont(42);
    Nuaxis->SetLabelFont(42);
    Nuaxis->Draw();

    


    /*    int npoints= 10;
    double xvec[10]={5.6,4.7,... };
    double yvec[10]={ };
    TGraph *mygraph = new TGraph(npoints,xvec,yvec);  
    mygraph->SetMarkerStyle(20);
    mygraph->Draw("P");  */

    // add zooms 
   
//    TPad *xpad= new TPad("xpad", "", 0.2,0.2,0.45,0.45); 
//    xpad->SetRightMargin(0.01);
//    xpad->SetTopMargin(0.); 
//    xpad->Draw(); 
//    xpad->cd(); 
//    TMultiGraph *xsed= sed->Clone("xsed");  
//    xsed->GetXaxis()->SetRangeUser(17.5,19.);
//    xsed->GetXaxis()->SetTitle("");
//    xsed->GetYaxis()->SetTitle("");
//    xsed->SetMinimum(-12.);
//    xsed->SetMaximum(-10.);
//   xsed->Draw("AC");
//    if(wehavedata){ 
//        TMultiGraph *xgrdata = grdata->Clone();
//        xgrdata->Draw("P"); 
//	} 
//      xsed->Draw("C");

//    csed->cd();
//    TPad *gpad= new TPad("gpad", "", 0.6,0.2,0.85,0.45);
//    gpad->SetRightMargin(0.01);
//    gpad->SetTopMargin(0.); 
//    gpad->Draw(); 
//    gpad->cd(); 
//    TMultiGraph *gsed= sed->Clone("gsed");  
//    gsed->GetXaxis()->SetRangeUser(25.5,27.);
//    gsed->GetXaxis()->SetTitle("");
//    gsed->GetYaxis()->SetTitle("");
//    gsed->Draw("AC");
//    if(wehavedata){ 
//        TMultiGraph *ggrdata = grdata->Clone();
//        ggrdata->Draw("P");
//	}
//    gsed->Draw("C"); // to have the fit above the data points


    // Calculate the Chi^2 between data points and model
 
    grdata->SaveAs("test_data.C");
   

}

