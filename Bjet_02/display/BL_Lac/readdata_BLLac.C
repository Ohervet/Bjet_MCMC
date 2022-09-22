{
  ifstream file0,file_gal,file_MAGIC;
  char str[20],str1[20],str2[50],info[10];
  int i, n0=79, n0_0=58,n0_1=7, n1a_2=43,n1a_3=8,n2=32,n2_1=9,n2_2=4,n2_3=6,n2_4=8,n2_5=5,nb0,nb1,nb2,nb3,n_Mb = 50;
  Double_t A,B,E[n0],F[n0],C[n0],D[n0],deltaEmin[n0],deltaEmax[n0],deltaF[n0],deltaFmax[n0],E1[n2],F1[n2],C1[n2],D1[n2],deltaEmin1a[n2],deltaEmax1[n2],deltaF1[n2],deltaF1max[n2],E1uplim[n2],F1uplim[n2],dEmin1auplim[n2],dEmax1uplim[n2],dF1uplim[n2],
  ,Euplim[n2],Fuplim[n2],dEminuplim[n2],dEmaxuplim[n2],dFuplim[n2];
  Double_t E1_1[n2_1],F1_1[n2_1],E1_2[n2_2],F1_2[n2_2],E1_3[n2_3],F1_3[n2_3],E1_4[n2_4],F1_4[n2_4],E1_5[n2_5],F1_5[n2_5],E_1[n0_1],F_1[n0_1],deltaEmin_1[n0_1],deltaEmax_1[n0_1],deltaF_1[n0_1],deltaFmax_1[n0_1],E_2[n1a_3],F_2[n1a_3];
  Double_t Freq_Mb[n_Mb], Flux_Mb[n_Mb],DFmin_Mb[n_Mb],DFmax_Mb[n_Mb];
  
  const double planck = 4.135667517e-15; // h in [eV s]

   file0.open("data/BL_Lac/BLLac_data.dat");
   file_MAGIC.open("data/BL_Lac/MAGIC_butterfly.dat");
   
   
   
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
	    /*
	    if(strcmp(info,"yes")>=0) 
	    { Ul[i]=1;	  
	      C[j]=E[i];
 	      D[j]=F[i];
	      j+=1;}
 	    else
	    { Ul[i]=0;}
   */
	  std::cout<<i<<"\t"<<E[i]<<"\t"<<F[i]<<" \t"<<deltaEmin[i]<<"\t"<<deltaEmax[i]<<"\t"<<deltaF[i]<<"\t"<<deltaFmax[i]<<"\t"<<str2<<std::endl;

            if (!file0.good()) break;

        } 
         
        
	//HESS Data
	for (int i=0;i<n0_1;i++){
	  E_1[i] = E[i+n0_0];
	  F_1[i] = F[i+n0_0];
	  deltaEmin_1[i] = deltaEmin[i+n0_0];
	  deltaEmax_1[i] = deltaEmax[i+n0_0];
	  deltaF_1[i] = deltaF[i+n0_0];
	  deltaFmax_1[i] = deltaFmax[i+n0_0];
	}   

	
	
        if(nb0>0){
          TGraph *gls = new TGraph(nb0,C,D);
 	  gls->SetMarkerStyle(1);
          gls->SetMarkerColor(15);
          gls->Draw("*");
        }
        

        
        //Papillon MAGIC low state 2007
        /*
        const int nt = 1000;
	Double_t EM[nt], EM[nt], FM[nt], FMmin[nt], FMmax[nt], DFM[nt];
	EM[0] = log10(1.5e11/planck);
	EM[nt] = log10(9.0e11/planck); 
	double step = (EM[nt]-EM[0])/nt;
	
	FM[0] = log10(2.0e-12 * pow(pow(10,EM[0])/3.0e11*planck,-1.1));
        for(i=1;i<nt;i++){
	  EM[i] = EM[i-1] + step;
	  FM[i] = log10(2.0e-12 * pow(pow(10,EM[i])/3.0e11*planck,-1.1));
	  //std::cout<<FM[i]<<std::endl;
	  FMmax[i] = FM[i];
	  FMmin[i] = FM[i];
	}
	//FM[0] = log10(2.671e-7 * pow(pow(10,EH[0])*planck,-0.46));
	FM[0] = log10(2.0e-12 * pow(pow(10,EM[0])/3.0e11*planck,-2.1));
        for(i=1;i<nt;i++){
	  EM[i] = EM[i-1] + step;
	  FM[i] = log10(2.0e-12 * pow(pow(10,EM[i])/3.0e11*planck,-2.1));
	  if(FM[i] > FMmax[i]){
	    FMmax[i] = FM[i];
	  }
	  if(FM[i] < FMmin[i]){
	    FMmin[i] = FM[i];
	  }
	}
	//FM[0] = log10(9.216e-3 * pow(pow(10,EM[0])*planck,-0.84));
	FM[0] = log10(3.4e-12 * pow(pow(10,EM[0])/3.0e11*planck,-1.1));
        for(i=1;i<nt;i++){
	  EM[i] = EM[i-1] + step;
	  FM[i] =log10(3.4e-12 * pow(pow(10,EM[i])/3.0e11*planck,-1.1));
	  if(FM[i] > FMmax[i]){
	    FMmax[i] = FM[i];
	  }
	  if(FM[i] < FMmin[i]){
	    FMmin[i] = FM[i];
	  }
	}
	//FM[0] = log10(3.487e-7 * pow(pow(10,EM[0])*planck,-0.46));
	FM[0] = log10(3.4e-12 * pow(pow(10,EM[0])/3.0e11*planck,-2.1));
        for(i=1;i<nt;i++){
	  EM[i] = EM[i-1] + step;
	  FM[i] = log10(3.4e-12 * pow(pow(10,EM[i])/3.0e11*planck,-2.1));
	  if(FM[i] > FMmax[i]){
	    FMmax[i] = FM[i];
	  }
	  if(FM[i] < FMmin[i]){
	    FMmin[i] = FM[i];
	  }
	  DFM[i] = FMmax[i] - FMmin[i];
	}*/

       //MAGIC butterfly 
      for(int i=0;i<n_Mb;i++)
      {
	file_MAGIC>>Freq>>Flux>>FluxLow>>FluxHigh;

            Freq_Mb[i]=log10(Freq);
 	    Flux_Mb[i]=log10(Flux);
	    DFmin_Mb[i] = Flux_Mb[i] - log10(Flux-FluxLow);
	    DFmax_Mb[i] = log10(FluxHigh+Flux)-Flux_Mb[i];
	    
	    //std::cout<<i<<"\t"<<Freq_Mb[i]<<"\t"<<Flux_Mb[i]<<" \t"<<DFmin_Mb[i]<<"\t"<<DFmax_Mb[i]<<std::endl;
   
            if (!file_MAGIC.good()) break;
      }
        
        
        
        
	

	Double_t minHz, maxHz,minF,maxF;
	minHz=8.;
	maxHz=28;    
	minF=-13.5;
	maxF=-9.5;

 	nb1=k;
	//display errorbars
	TGraphAsymmErrors *gr0 = new TGraphAsymmErrors(n0,E,F,deltaEmin,deltaEmax,deltaF,deltaFmax);
	//gr1->SetMarkerStyle(20);
	gr0->SetMarkerColor(1);
	gr0->Draw("AP");
	gr0->SetTitle("");
	gr0->GetYaxis()->SetTitle("log (#nu F#nu [erg cm^{-2} s^{-1} ] )");
	gr0->GetYaxis()->SetTitleSize(0.04);
	gr0->GetYaxis()->SetTitleOffset(1.3);
	gr0->GetYaxis()->SetLabelSize(0.04);
	gr0->GetXaxis()->SetTitle("log (#nu [Hz])");
	gr0->GetXaxis()->SetTitleSize(0.04);
	gr0->GetXaxis()->SetTitleOffset(1.2);
	gr0->GetXaxis()->SetLabelSize(0.04);
	gr0->GetXaxis()->SetRangeUser(minHz,maxHz);
	gr0->SetMinimum(minF);
	gr0->SetMaximum(maxF);

	Double_t mineV, maxeV;
	Double_t logh=log10(planck);
	mineV=minHz + logh;
	maxeV=maxHz + logh;

	// add energy axis
	TGaxis *eVaxis = new TGaxis(minHz,maxF,maxHz, maxF,mineV,maxeV,510,"-");
	eVaxis->SetTitle("log(E [eV])");
	eVaxis->SetTitleOffset(1.2);
	eVaxis->Draw();
	

	
	
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

	//MAGIC butterfly 
	gr_Mb = new TGraphAsymmErrors(n_Mb,Freq_Mb,Flux_Mb,0,0,DFmin_Mb,DFmax_Mb);
	gr_Mb->SetFillColor(17);
	gr_Mb->Draw("3");
	
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
	grdata->Add(gr0_0);
	grdata->Add(gr_Mb,"3");
	if (l>0){
	  grdata->Add(gr0a,"|>");
	}
	
	


	
	return 1;
}
