{
  ifstream file0,file_gal;
  char str[20],str1[20],str2[40];
  int i, low_state; 
  
    low_state = 0;
  
  if (low_state == 1){
    int n0 = 106, n0_1 =106; //nb de pts SED
   file0.open("data/VER_J0521+211/VERJ0521_low_state.dat");
  }else{
        int n0 = 92, n0_1 =92; //nb de pts SED
   file0.open("data/VER_J0521+211/VERJ0521_high_state.dat");
  }
  
  int n0_0=78,l_1,l_2,nl =20,nb0,nb1,nb2,nb3;
  int n1 = n0_1-n0_0;
  Double_t A,B,E[n0],F[n0],C[n0],D[n0],deltaEmin[n0],deltaEmax[n0],deltaF[n0],deltaFmax[n0],E_0[n0_0],F_0[n0_0],deltaEmin_0[n0_0],deltaEmax_0[n0_0],deltaF_0[n0_0],deltaFmax_0[n0_0]
  ,E_1[n1],F_1[n1],deltaEmin_1[n1],deltaEmax_1[n1],deltaF_1[n1],deltaFmax_1[n1],Euplim_1[nl],Fuplim_1[nl],dEminuplim_1[nl],dEmaxuplim_1[nl],dFuplim_1[nl];
  const double planck = 4.135667517e-15; // h in [eV s]



   
   
   double     nb_g=0;
   const int   ng = 1298;
   Double_t   freq_g1[ng],flux_g1[ng],freq_g[ng],flux_g[ng],fflux_g1[ng],fflux_g[ng];
   
   file_gal.open("data/HGS_13.dat"); // host galaxy
   
//******************************************************************
// Host galaxy spectrum
//******************************************************************

      if(!file_gal){
         cout<<"Host galaxy file not found !"<<endl;     
       } 
       else{
         for(i=0;i<ng;i++){ 
            file_gal>>freq_g1[i]>>flux_g[i]>>fflux_g1[i]>>freq_g[i]>>fflux_g[i];
            if (!file_gal.good()) break;
         }
         nb_g=i;
       }       
       
       
       


    if(!file0){
       cout<<"data file not found"<<endl;
       return 0;
    }
 
     int j=0, k=0, l=0, l1=0;
     
     for(int i=1;i<=11;i++){
        file0>>str;
     }

      std::cout<<"i\tE\tF\t\tdE-\tdE+\tdF-\tdF+\tObservatory"<<std::endl;
      
      for(int i=0;i<n0;i++)
      { file0>>A>>B>>deltaEmin[i]>>deltaEmax[i]>>deltaF[i]>>deltaFmax[i]>>str2;


            //A = A / planck; 
            E[i]=log10(A/planck);
 	    F[i]=log10(B);
	    if (deltaEmin[i] != 0){
	      deltaEmin[i]= E[i] - log10(A/planck - deltaEmin[i]/planck);
	      deltaEmax[i]= -E[i] + log10(A/planck + deltaEmax[i]/planck);
	    }
	      
	    
	
 	    deltaF[i]=F[i]-log10(B-deltaF[i]);
	    deltaFmax[i]=-F[i]+log10(B+deltaFmax[i]);
	    /*
	    if (deltaFmax[i] == 0){
	      Euplim[l] = E[i];
	      Fuplim[l] = F[i];
	      dEminuplim[l]= deltaEmin[i];
	      dEmaxuplim[l]= deltaEmax[i];
	      dFuplim[l]= 0.2;
	      l+=1;
	    }
*/
	  std::cout<<i<<"\t"<<E[i]<<"\t"<<F[i]<<" \t"<<deltaEmin[i]<<"\t"<<deltaEmax[i]<<"\t"<<deltaF[i]<<"\t"<<deltaFmax[i]<<"\t"<<str2<<std::endl;

            if (!file0.good()) break;

        } 
        
                if(nb0>0){
          TGraph *gls = new TGraph(nb0,C,D);
 	  gls->SetMarkerStyle(1);
          gls->SetMarkerColor(15);
          gls->Draw("*");
        }
        
        
        //archival data
         memcpy(E_0, E, n0_0*sizeof(double)); 
	 memcpy(F_0, F, n0_0*sizeof(double)); 
	 memcpy(deltaEmin_0, deltaEmin, n0_0*sizeof(double)); 
	 memcpy(deltaEmax_0, deltaEmax, n0_0*sizeof(double)); 
	 memcpy(deltaF_0, deltaF, n0_0*sizeof(double)); 
	 memcpy(deltaFmax_0, deltaFmax, n0_0*sizeof(double)); 
	 


        

	

	Double_t minHz, maxHz,minF,maxF;
	minHz=5.;
	maxHz=28;    
	minF=-15.;
	maxF=-10.;

	//display errorbars
	TGraphAsymmErrors *gr0 = new TGraphAsymmErrors(n0_0,E_0,F_0,deltaEmin_0,deltaEmax_0,deltaF_0,deltaFmax_0);
	gr0->SetLineColor(16);
	gr0->SetMarkerStyle(20);
	gr0->SetMarkerColor(16);
	gr0->Draw("AP");
	gr0->SetTitle("VER_J0521+211 Data");
	gr0->GetYaxis()->SetTitle("log (#nu F#nu [erg cm^{-2} s^{-1} ] )");
	gr0->GetYaxis()->SetTitleSize(0.04);
	gr0->GetYaxis()->SetTitleOffset(1.3);
	gr0->GetYaxis()->SetLabelSize(0.04);
	gr0->GetXaxis()->SetTitle("log (#nu [Hz])");
	gr0->GetXaxis()->SetTitleSize(0.04);
	gr0->GetXaxis()->SetTitleOffset(1.3);
	gr0->GetXaxis()->SetLabelSize(0.04);
	gr0->GetXaxis()->SetRangeUser(minHz,maxHz);
	gr0->SetMinimum(minF);
	gr0->SetMaximum(maxF);
	

	//Low state
         memcpy(E_1, E+n0_0, (n0_1-n0_0)*sizeof(double)); 
	 memcpy(F_1, F+n0_0, (n0_1-n0_0)*sizeof(double)); 
	 memcpy(deltaEmin_1, deltaEmin+n0_0, (n0_1-n0_0)*sizeof(double)); 
	 memcpy(deltaEmax_1, deltaEmax+n0_0, (n0_1-n0_0)*sizeof(double)); 
	 memcpy(deltaF_1, deltaF+n0_0, (n0_1-n0_0)*sizeof(double)); 
	 memcpy(deltaFmax_1, deltaFmax+n0_0, (n0_1-n0_0)*sizeof(double)); 
	 
	for(int i=0;i<(n0_1-n0_0);i++){
	  if (deltaFmax_1[i] == 0){
	      Euplim_1[l_1] = E_1[i];
	      Fuplim_1[l_1] = F_1[i];
	      dEminuplim_1[l_1]= deltaEmin_1[i];
	      dEmaxuplim_1[l_1]= deltaEmax_1[i];
	      dFuplim_1[l_1]= 0.2;
	      l_1+=1;
	    }
	}
	 
	TGraphAsymmErrors *gr1 = new TGraphAsymmErrors((n0_1-n0_0),E_1,F_1,deltaEmin_1,deltaEmax_1,deltaF_1,deltaFmax_1);
	gr1->SetMarkerStyle(20);
	if (low_state == 1){
	gr1->SetMarkerColor(4);
	gr1->SetLineColor(4);
	}else{
	gr1->SetMarkerColor(2);
	gr1->SetLineColor(2);
	}
	gr1->Draw("P");
	
	//display upper limits
	if (l_1>0){
	  TGraphAsymmErrors *gr1l = new TGraphAsymmErrors(l_1,Euplim_1,Fuplim_1,0,0,dFuplim_1,0);
	  gr1l->SetMarkerStyle(20);
	  gr1l->SetMarkerColor(4);
	  gr1l->SetLineColor(4);
	  gr1l->SetFillColor(4);
	  gr1l->Draw("P|>");
	}
	
	
	//High state
	/*
         memcpy(E_2, E+n0_1, (n0-n0_1)*sizeof(double)); 
	 memcpy(F_2, F+n0_1, (n0-n0_1)*sizeof(double)); 
	 memcpy(deltaEmin_2, deltaEmin+n0_1, (n0-n0_1)*sizeof(double)); 
	 memcpy(deltaEmax_2, deltaEmax+n0_1, (n0-n0_1)*sizeof(double)); 
	 memcpy(deltaF_2, deltaF+n0_1, (n0-n0_1)*sizeof(double)); 
	 memcpy(deltaFmax_2, deltaFmax+n0_1, (n0-n0_1)*sizeof(double)); 
	 

	for(int i=0;i<(n0-n0_1);i++){
	  if (deltaFmax_2[i] == 0){
	      Euplim_2[l_2] = E_2[i];
	      Fuplim_2[l_2] = F_2[i];
	      dEminuplim_2[l_2]= deltaEmin_2[i];
	      dEmaxuplim_2[l_2]= deltaEmax_2[i];
	      dFuplim_2[l_2]= 0.2;
	      l_2+=1;
	    }
	}


	TGraphAsymmErrors *gr2a = new TGraphAsymmErrors(n2,E_2,F_2,deltaEmin_2,deltaEmax_2,deltaF_2,deltaFmax_2);
	gr2a->SetMarkerStyle(20);
	gr2a->SetMarkerColor(2);
	gr2a->SetLineColor(2);
	gr2a->Draw("P");
	

	//display upper limits
	if (l_2>0){
	  TGraphAsymmErrors *gr2l = new TGraphAsymmErrors(l_2,Euplim_2,Fuplim_2,0,0,dFuplim_2,0);
	  gr2l->SetMarkerStyle(20);
	  gr2l->SetMarkerColor(2);
	  gr2l->SetLineColor(2);
	  gr2l->SetFillColor(2);
	  gr2l->Draw("P|>");
	}
	*/
	
	
	    // Host galaxy
    
        if(nb_g>0){
          gr10 = new TGraph(nb_g,freq_g,fflux_g);
	  gr10->SetLineColor(2);
          gr10->SetLineStyle(1);
	  //gr10->Draw("");
	  //sed->Add(gr10);
        }
	
	
	TMultiGraph *grdata = new TMultiGraph();
	grdata->SetTitle("SED Data");
	grdata->Add(gr0);
	//grdata->Add(gr2a);
	grdata->Add(gr1);
	if (l_1>0){
	  grdata->Add(gr1l,"P|>");
	  /*
	if (l_2>0){
	  grdata->Add(gr2l,"P|>");
	}
	*/


	
	return 1;
}
