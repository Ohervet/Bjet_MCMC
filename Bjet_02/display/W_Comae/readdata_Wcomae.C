{
  ifstream file0,file_gal,file_VERITAS;
  char str[20],str1[20],str2[50],info[10];
  int i, n0=79, n0_0=24,n0_V=8, n1a_2=43,n1a_3=8,n2=32,n2_1=9,n2_2=4,n2_3=6,n2_4=8,n2_5=5,nb0,nb1,nb2,nb3,n_Vb = 50;
  Double_t A,B,E[n0],F[n0],C[n0],D[n0],deltaEmin[n0],deltaEmax[n0],deltaF[n0],deltaFmax[n0],E1[n2],F1[n2],C1[n2],D1[n2],deltaEmin1a[n2],deltaEmax1[n2],deltaF1[n2],deltaF1max[n2],E1uplim[n2],F1uplim[n2],dEmin1auplim[n2],dEmax1uplim[n2],dF1uplim[n2],
  ,Euplim[n2],Fuplim[n2],dEminuplim[n2],dEmaxuplim[n2],dFuplim[n2];
  Double_t E1_1[n2_1],F1_1[n2_1],E1_2[n2_2],F1_2[n2_2],E1_3[n2_3],F1_3[n2_3],E1_4[n2_4],F1_4[n2_4],E1_5[n2_5],F1_5[n2_5];
  Double_t E_V[n0_V],F_V[n0_V],deltaEmin_V[n0_V],deltaEmax_V[n0_V],deltaF_V[n0_V],deltaFmax_V[n0_V];
  Double_t Freq_Vb[n_Vb], Flux_Vb[n_Vb],DFmin_Vb[n_Vb],DFmax_Vb[n_Vb];
  
  const double planck = 4.135667517e-15; // h in [eV s]
   
   //low state
   //file0.open("data/W_Comae/Wcomae_low_state.dat");
   //file_VERITAS.open("data/W_Comae/VERITAS_butterfly_low_allPeriods.dat");
   
   //high state
   file0.open("data/W_Comae/Wcomae_high_state.dat");
   file_VERITAS.open("data/W_Comae/VERITAS_butterfly_high.dat");
   
   
   
   double     nb_g=0;
   const int   ng = 1298;
   double    Freq, FreqLow, FreqHigh, Flux, FluxLow, FluxHigh, DFlux;
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

      std::cout<<"i\tE\tF\t\tdE-\tdE+\tdF-\t\tdF+\t\tObservatory"<<std::endl;
      
      for(int i=0;i<n0;i++)
      {file0>>A>>B>>deltaEmin[i]>>deltaEmax[i]>>deltaF[i]>>deltaFmax[i]>>str2;


            A = A / planck; 
            E[i]=log10(A);
 	    F[i]=log10(B);
	    deltaEmin[i]=E[i]-log10(A-deltaEmin[i]);
	    deltaEmax[i]=-E[i]+log10(A+deltaEmax[i]);
	
 	    deltaF[i]=F[i]-log10(B-deltaF[i]);
	    deltaFmax[i]=-F[i]+log10(B+deltaFmax[i]);
	    if (deltaFmax[i] == 0){
	      Euplim[l] = E[i];
	      Fuplim[l] = F[i];
	      dEminuplim[l]= deltaEmin[i];
	      dEmaxuplim[l]= deltaEmax[i];
	      dFuplim[l]= 0.2;
	      l+=1;
	    }
	  std::cout<<i<<"\t"<<E[i]<<"\t"<<F[i]<<" \t"<<deltaEmin[i]<<"\t"<<deltaEmax[i]<<"\t"<<deltaF[i]<<"\t"<<deltaFmax[i]<<"\t"<<str2<<std::endl;

            if (!file0.good()) break;

        } 
         
        
	//VERITAS Data
	for (int i=0;i<n0_V;i++){
	  E_V[i] = E[i+n0_0];
	  F_V[i] = F[i+n0_0];
	  deltaEmin_V[i] = deltaEmin[i+n0_0];
	  deltaEmax_V[i] = deltaEmax[i+n0_0];
	  deltaF_V[i] = deltaF[i+n0_0];
	  deltaFmax_V[i] = deltaFmax[i+n0_0];
	}   

	
	
        if(nb0>0){
          TGraph *gls = new TGraph(nb0,C,D);
 	  gls->SetMarkerStyle(1);
          gls->SetMarkerColor(15);
          gls->Draw("*");
        }
        
        
       //VERITAS butterfly 
      for(int i=0;i<n_Vb;i++)
      {
	file_VERITAS>>Freq>>Flux>>FluxLow>>FluxHigh;

            Freq_Vb[i]=log10(Freq);
 	    Flux_Vb[i]=log10(Flux);
	    DFmin_Vb[i] = Flux_Vb[i] - log10(Flux-FluxLow);
	    DFmax_Vb[i] = log10(FluxHigh+Flux)-Flux_Vb[i];
	    
	    std::cout<<i<<"\t"<<Freq_Vb[i]<<"\t"<<Flux_Vb[i]<<" \t"<<DFmin_Vb[i]<<"\t"<<DFmax_Vb[i]<<std::endl;
   
            if (!file_VERITAS.good()) break;
      }
        
        
        
        
        
	

	Double_t minHz, maxHz,minF,maxF;
	minHz=22.;//5.;
	maxHz=28;    
	minF=-14.;
	maxF=-10.5;//-9.8;

 	nb1=k;
	//display errorbars
	TGraphAsymmErrors *gr0 = new TGraphAsymmErrors(n0,E,F,deltaEmin,deltaEmax,deltaF,deltaFmax);
	//gr1->SetMarkerStyle(20);
	gr0->SetMarkerColor(1);
	gr0->Draw("AP");
	gr0->SetTitle("PKS_0625-354");
	gr0->GetYaxis()->SetTitle("log (#nu F#nu [erg cm^{-2} s^{-1} ] )");
	gr0->GetYaxis()->SetTitleSize(0.04);
	gr0->GetYaxis()->SetTitleOffset(1.3);
	gr0->GetYaxis()->SetLabelSize(0.04);
	gr0->GetXaxis()->SetTitle("log (#nu [Hz])");
	gr0->GetXaxis()->SetTitleSize(0.04);
	gr0->GetXaxis()->SetTitleOffset(1.1);
	gr0->GetXaxis()->SetLabelSize(0.04);
	gr0->GetXaxis()->SetRangeUser(minHz,maxHz);
	gr0->SetMinimum(minF);
	gr0->SetMaximum(maxF);

	
	
	//display dots
	TGraph *gr0_0 = new TGraph(n0,E,F);
	gr0_0->SetMarkerStyle(20);
	gr0_0->SetMarkerColor(1);
	gr0_0->Draw("P");
	
	//display upper limits
	if (l>0){
	  TGraphAsymmErrors *gr0a = new TGraphAsymmErrors(l,Euplim,Fuplim,0,0,dFuplim,0);
	  gr0a->SetMarkerStyle(20);
	  gr0a->SetMarkerColor(1);
	  gr0a->Draw("|>");
	}
	
	
	//VERITAS butterfly 
	gr_Vb = new TGraphAsymmErrors(n_Vb,Freq_Vb,Flux_Vb,0,0,DFmin_Vb,DFmax_Vb);
	gr_Vb->SetFillColor(kRed-9);
	gr_Vb->Draw("3");
	
        
        
        //VERITAS data
        gr_V = new TGraphAsymmErrors(n0_V,E_V,F_V,deltaEmin_V,deltaEmax_V,deltaF_V,deltaFmax_V);
        gr_V->SetLineColor(kRed+2);
        gr_V->Draw("P");
        //display dots
	grV_0 = new TGraph(n0_V,E_V,F_V);
	grV_0->SetMarkerStyle(20);
	grV_0->SetMarkerColor(kRed+2);
        grV_0->Draw("P");
        
        
        
        
        
	    // Host galaxy
    
        if(nb_g>0){
          gr10 = new TGraph(nb_g,freq_g,fflux_g);
	  gr10->SetLineColor(2);
          gr10->SetLineStyle(1);
	  //gr10->Draw("");
	  //sed->Add(gr10);
        }
	
	
	TMultiGraph *grdata = new TMultiGraph();
	grdata->Add(gr_Vb,"3");
	grdata->SetTitle("SED Data");
	grdata->Add(gr0);
	grdata->Add(gr0_0);
        grdata->Add(gr_V);
        grdata->Add(grV_0);
	if (l>0){
	  grdata->Add(gr0a,"|>");
	}
	
	


	
	return 1;
}
