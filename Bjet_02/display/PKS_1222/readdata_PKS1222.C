{
  ifstream file0,file1,file_gal,file_Fermi, file_VERITAS;
  char str[20],str1[20],str2[50],info[10];
  int i, count, countUL, n0=79,n1=28,nb0,nb1, gr0_0, nfs = 50, nVs = 50;
  Double_t A,B,E[n0],F[n0],C[n0],D[n0],deltaEmin[n0],deltaEmax[n0],deltaFmin[n0],deltaFmax[n0],Euplim[n0],Fuplim[n0],dEminuplim[n0],dEmaxuplim[n0],dFuplim[n0],
  E1[n1],F1[n1],deltaEmin1[n1],deltaEmax1[n1],deltaFmin1[n1],deltaFmax1[n1],Euplim1[n1],Fuplim1[n1],dEminuplim1[n1],dEmaxuplim1[n1],dFuplim1[n1];

  int nO = 0;
  Double_t E_O[nO], F_O[nO],deltaEmin_O[nO],deltaEmax_O[nO],deltaFmin_O[nO],deltaFmax_O[nO];
  int nS = 2;
  Double_t E_S[nS], F_S[nS],deltaEmin_S[nS],deltaEmax_S[nS],deltaFmin_S[nS],deltaFmax_S[nS];
  int nU = 6;
  Double_t E_uvot[nU], F_uvot[nU],deltaEmin_uvot[nU],deltaEmax_uvot[nU],deltaFmin_uvot[nU],deltaFmax_uvot[nU];
  int nx = 62;
  Double_t E_xrt[nx], F_xrt[nx],deltaEmin_xrt[nx],deltaEmax_xrt[nx],deltaFmin_xrt[nx],deltaFmax_xrt[nx];
  //Double_t Freq_F[nfs], Flux_F[nfs], DFmin_F[nfs], DFmax_F[nfs];
  int n_Fl = 0;
  Double_t E_Fl[n_Fl], F_Fl[n_Fl],deltaF_Fl[n_Fl];
  int n_F = 6;
  Double_t E_F[n_F], F_F[n_F],deltaEmin_F[n_F],deltaEmax_F[n_F],deltaFmin_F[n_F],deltaFmax_F[n_F];
  int n_Vl = 1;
  Double_t E_Vl[n_Vl], F_Vl[n_Vl],deltaF_Vl[n_Vl];
  int n_V = 3;
  Double_t E_V[n_V], F_V[n_V],deltaEmin_V[n_V],deltaEmax_V[n_V],deltaFmin_V[n_V],deltaFmax_V[n_V];

  
  const double planck = 4.135667517e-15; // h in [eV s]

   file0.open("data/PKS_1222+216/PKS1222+216_SED.dat");
   file1.open("data/NGC_1275/NGC1275_Archival_SED.dat");
   
   
   double     nb_g=0;
   const int   ng = 1108;
   double    Freq, FreqLow, FreqHigh, Flux, FluxLow, FluxHigh, DFlux, DFluxLow, DFluxHigh;
   Double_t   freq_g1[ng],flux_g1[ng],freq_g[ng],flux_g[ng],fflux_g1[ng],fflux_g[ng];
   
   //file_gal.open("data/HESS_J1943+213/HGS_15_cut_z02.dat"); // host galaxy
   file_gal.open("data/HESS_J1943+213/HGS_15_cut_z016.dat"); // host galaxy
   
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
    if(!file_Fermi){
       cout<<"Fermi spectrum file not found"<<endl;
    }
 
     int j=0, k=0, l=0, l1=0;
     
     for(int i=1;i<=11;i++){
        file0>>str;
     }

      std::cout<<"i\tE\tF\t\tdE-\tdE+\tdF-\t\tdF+\t\tObservatory"<<std::endl;
      
      for(int i=0;i<n0;i++)
      {file0>>A>>B>>deltaEmin[i]>>deltaEmax[i]>>deltaFmin[i]>>deltaFmax[i]>>str2;


            A = A / planck; 
            E[i]=log10(A);
 	    F[i]=log10(B);
	    deltaEmin[i]=E[i]-log10(A-deltaEmin[i]/ planck);
	    deltaEmax[i]=-E[i]+log10(A+deltaEmax[i]/ planck);
	
 	    deltaFmin[i]=F[i]-log10(B-deltaFmin[i]);
	    deltaFmax[i]=-F[i]+log10(B+deltaFmax[i]);
	    if (deltaFmax[i] == 0){
	      Euplim[l] = E[i];
	      Fuplim[l] = F[i];
	      dEminuplim[l]= deltaEmin[i];
	      dEmaxuplim[l]= deltaEmax[i];
	      dFuplim[l]= 0.2;
	      l+=1;
	    }
	  std::cout<<i<<"\t"<<E[i]<<"\t"<<F[i]<<" \t"<<deltaEmin[i]<<"\t"<<deltaEmax[i]<<"\t"<<deltaFmin[i]<<"\t"<<deltaFmax[i]<<"\t"<<str2<<std::endl;

            if (!file0.good()) break;

        } 
        
        
     for(int i=1;i<=11;i++){
        file1>>str;
     }
      
      for(int i=0;i<n1;i++)
      {file1>>A>>B>>deltaEmin1[i]>>deltaEmax1[i]>>deltaFmin1[i]>>deltaFmax1[i]>>str2;


            A = A / planck; 
            E1[i]=log10(A);
 	    F1[i]=log10(B);
	    deltaEmin1[i]=E1[i]-log10(A-deltaEmin1[i]/ planck);
	    deltaEmax1[i]=-E1[i]+log10(A+deltaEmax1[i]/ planck);
	
 	    deltaFmin1[i]=F1[i]-log10(B-deltaFmin1[i]);
	    deltaFmax1[i]=-F1[i]+log10(B+deltaFmax1[i]);
	    if (deltaFmax1[i] == 0){
	      Euplim1[l] = E1[i];
	      Fuplim1[l] = F1[i];
	      dEminuplim1[l]= deltaEmin1[i];
	      dEmaxuplim1[l]= deltaEmax1[i];
	      dFuplim1[l]= 0.2;
	      l+=1;
	    }

            if (!file1.good()) break;

        }         
        
        
        
        
        
        
 /*       
         
       //Fermi spectrum
      for(int i=0;i<nfs;i++)
      {
	file_Fermi>>Freq>>Flux>>DFluxLow>>DFluxHigh;

            Freq_F[i]=log10(Freq);
 	    Flux_F[i]=log10(Flux);
	    DFmin_F[i] = Flux_F[i] - log10((Flux-DFluxLow));
	    DFmax_F[i] = log10((Flux+DFluxHigh)) - Flux_F[i];
	    
	    //std::cout<<i<<"\t"<<Freq_F[i]<<"\t"<<Flux_F[i]<<" \t"<<DFmin_F[i]<<"\t"<<DFmax_F[i]<<std::endl;
   
            if (!file_Fermi.good()) break;

        }
        
       //VERITAS spectrum
      for(int i=0;i<nVs;i++)
      {
	file_VERITAS>>Freq>>Flux>>DFluxLow>>DFluxHigh;

            Freq_Vs[i]=log10(Freq);
 	    Flux_Vs[i]=log10(Flux);
	    DFmin_Vs[i] = Flux_Vs[i] - log10((Flux-DFluxLow));
	    DFmax_Vs[i] = log10((Flux+DFluxHigh)) - Flux_Vs[i];
   
            if (!file_VERITAS.good()) break;

        }

	
	*/
        if(nb0>0){
          TGraph *gls = new TGraph(nb0,C,D);
 	  gls->SetMarkerStyle(1);
          gls->SetMarkerColor(15);
          gls->Draw("*");
        }
        
        TCanvas *csed = new TCanvas("SED1","SED1");
        csed->cd();
        csed->GetFrame()->SetBorderSize(12);
        gPad->SetTopMargin(0.1);
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.05);
        gPad->SetTicky();
    
	Double_t minHz, maxHz,minF,maxF;
        minHz=9.0;
        maxHz=28.2;    
        minF=-14.;
        maxF=-7.5; 

 	nb1=k;
	//display errorbars
	TGraphAsymmErrors *gr0 = new TGraphAsymmErrors(n1,E1,F1,deltaEmin1,deltaEmax1,deltaFmin1,deltaFmax1);
	gr0->SetMarkerStyle(20);
	gr0->SetMarkerColor(0);//15
        gr0->SetLineColor(0);
	gr0->Draw("AP");
	gr0->SetTitle("");
	gr0->GetYaxis()->SetTitle("log (#nu F#nu [erg cm^{-2} s^{-1} ] )");
	gr0->GetYaxis()->SetTitleSize(0.04);
	gr0->GetYaxis()->SetTitleOffset(1.4);
	gr0->GetYaxis()->SetLabelSize(0.04);
	gr0->GetXaxis()->SetTitle("log (#nu [Hz])");
	gr0->GetXaxis()->SetTitleSize(0.04);
	gr0->GetXaxis()->SetTitleOffset(1.2);
	gr0->GetXaxis()->SetLabelSize(0.04);
        gr0->GetXaxis()->SetLimits(minHz,maxHz);
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
        eVaxis->SetTitleFont(42);
        eVaxis->SetLabelFont(42);
        eVaxis->Draw();
	

	
	
	//display dots
        /*
	TGraph *gr0_0 = new TGraph(5,E,F);
	gr0_0->SetMarkerStyle(20);
	gr0_0->SetMarkerColor(1);
	gr0_0->Draw("P");
	
	//display upper limits
	if (l>0){
	  TGraphAsymmErrors *gr0a = new TGraphAsymmErrors(l,Euplim,Fuplim,0,0,dFuplim,0);
	  gr0a->SetMarkerStyle(20);
	  gr0a->SetMarkerColor(1);
	  //gr0a->Draw("|>");
	}*/
/*
	// Host galaxy
        if(nb_g>0){
          gr10 = new TGraph(nb_g,freq_g,fflux_g);
	  gr10->SetLineColor(2);
          gr10->SetLineStyle(1);
	  //gr10->Draw("C");
        }
*/
	//OVRO
        count = 0;
	for(i=0;i<nO;i++){
	  E_O[i] = E[count+i];
	  F_O[i] = F[count+i];
	  deltaEmin_O[i] = deltaEmin[count+i];
	  deltaEmax_O[i] = deltaEmax[count+i];
	  deltaFmin_O[i] = deltaFmin[count+i];
	  deltaFmax_O[i] = deltaFmax[count+i];
	  }

	gr_OVRO = new TGraphAsymmErrors(nO,E_O,F_O,deltaEmin_O,deltaEmax_O,deltaFmin_O,deltaFmax_O);
	gr_OVRO->SetMarkerColor(1);
	gr_OVRO->SetLineColor(1);
	gr_OVRO->SetMarkerStyle(20);
	gr_OVRO->SetMarkerColor(1);
	gr_OVRO->Draw("P");  
        

	//Steward Observatory
        count += nO;
	for(i=0;i<nS;i++){
	  E_S[i] = E[count+i];
	  F_S[i] = F[count+i];
	  deltaEmin_S[i] = deltaEmin[count+i];
	  deltaEmax_S[i] = deltaEmax[count+i];
	  deltaFmin_S[i] = deltaFmin[count+i];
	  deltaFmax_S[i] = deltaFmax[count+i];
	  }

	gr_Steward = new TGraphAsymmErrors(nS,E_S,F_S,deltaEmin_S,deltaEmax_S,deltaFmin_S,deltaFmax_S);
	gr_Steward->SetMarkerColor(2);
	gr_Steward->SetLineColor(2);
	gr_Steward->SetMarkerStyle(21);
	gr_Steward->SetMarkerColor(2);
	gr_Steward->Draw("P");        

        
	//UVOT
        count += nS;
	for(i=0;i<nU;i++){
	  E_uvot[i] = E[count+i];
	  F_uvot[i] = F[count+i];
	  deltaEmin_uvot[i] = deltaEmin[count+i];
	  deltaEmax_uvot[i] = deltaEmax[count+i];
	  deltaFmin_uvot[i] = deltaFmin[count+i];
	  deltaFmax_uvot[i] = deltaFmax[count+i];
	  }

	gr_uvot = new TGraphAsymmErrors(nU,E_uvot,F_uvot,deltaEmin_uvot,deltaEmax_uvot,deltaFmin_uvot,deltaFmax_uvot);
        gr_uvot->SetLineColor(kOrange+2);
	gr_uvot->SetMarkerStyle(33);
        gr_uvot->SetMarkerSize(1.5);
	gr_uvot->SetMarkerColor(kOrange+2);
	gr_uvot->Draw("P");      

        
	//XRT
        count += nU;
	for(i=0;i<nx;i++){
	  E_xrt[i] = E[count+i];
	  F_xrt[i] = F[count+i];
	  deltaEmin_xrt[i] = deltaEmin[count+i];
	  deltaEmax_xrt[i] = deltaEmax[count+i];
	  deltaFmin_xrt[i] = deltaFmin[count+i];
	  deltaFmax_xrt[i] = deltaFmax[count+i];
	  }

	gr_xrt = new TGraphAsymmErrors(nx,E_xrt,F_xrt,deltaEmin_xrt,deltaEmax_xrt,deltaFmin_xrt,deltaFmax_xrt);
        gr_xrt->SetLineColor(kGreen+2);
	gr_xrt->SetMarkerStyle(29);
	gr_xrt->SetMarkerColor(kGreen+2);
	gr_xrt->Draw("P");          
/*
	//Fermi spectrum
	gr_FERMI = new TGraphAsymmErrors(nfs,Freq_F,Flux_F,0,0,DFmin_F,DFmax_F);
	gr_FERMI->SetFillColor(kBlue-7);
	//gr_FERMI->Draw("3");
*/

	//Fermi data
	count += nx;
	for(i=0;i<n_F;i++){
	  E_F[i] = E[count+i];
	  F_F[i] = F[count+i];
	  deltaEmin_F[i] = deltaEmin[count+i];
	  deltaEmax_F[i] = deltaEmax[count+i];
	  deltaFmin_F[i] = deltaFmin[count+i];
	  deltaFmax_F[i] = deltaFmax[count+i];
	}
	gr_F = new TGraphAsymmErrors(n_F,E_F,F_F,deltaEmin_F,deltaEmax_F,deltaFmin_F,deltaFmax_F);
	gr_F->SetMarkerStyle(24);
        gr_F->SetMarkerSize(0.8);
        gr_F->SetMarkerColor(kBlue+2);
	gr_F->SetLineColor(kBlue+2);
	gr_F->Draw("P");

	if (l>0){
          countUL = 0;
	  for(i=0;i<n_Fl;i++){
	  E_Fl[i] = Euplim[countUL+i];
	  F_Fl[i] = Fuplim[countUL+i];
	  deltaF_Fl[i] = dFuplim[countUL+i];
	  std::cout<<i<<"\t"<<E_Fl[i]<<"\t"<<F_Fl[i]<<" \t"<<deltaF_Fl[i]<<"\t"<<std::endl;
	}
	  TGraphAsymmErrors *gr_Fl = new TGraphAsymmErrors(n_Fl,E_Fl,F_Fl,0,0,deltaF_Fl,0);
	  gr_Fl->SetMarkerStyle(24);
          gr_Fl->SetMarkerSize(0.8);
	  gr_Fl->SetLineColor(kBlue+2);
	  gr_Fl->SetFillColor(kBlue+2);
          gr_Fl->SetMarkerColor(kBlue+2);
	  gr_Fl->Draw("|>");
	}

	//VERITAS data
	count += n_F;
	for(i=0;i<n_V;i++){
	  E_V[i] = E[count+i];
	  F_V[i] = F[count+i];
	  deltaEmin_V[i] = deltaEmin[count+i];
	  deltaEmax_V[i] = deltaEmax[count+i];
	  deltaFmin_V[i] = deltaFmin[count+i];
	  deltaFmax_V[i] = deltaFmax[count+i];
	}
	gr_V = new TGraphAsymmErrors(n_V,E_V,F_V,deltaEmin_V,deltaEmax_V,deltaFmin_V,deltaFmax_V);
	gr_V->SetLineColor(1);
        gr_V->SetMarkerStyle(20);
        //gr_V->SetMarkerSize(0.8);
	gr_V->SetMarkerColor(1);
	gr_V->Draw("P");
        
	if (l>0){
          countUL += n_Fl;
	  for(i=0;i<n_Vl;i++){
	  E_Vl[i] = Euplim[countUL+i];
	  F_Vl[i] = Fuplim[countUL+i];
	  deltaF_Vl[i] = dFuplim[countUL+i];
	  std::cout<<i<<"\t"<<E_Vl[i]<<"\t"<<F_Vl[i]<<" \t"<<deltaF_Vl[i]<<"\t"<<std::endl;
	}
	  TGraphAsymmErrors *gr_Vl = new TGraphAsymmErrors(n_Vl,E_Vl,F_Vl,0,0,deltaF_Vl,0);
	  gr_Vl->SetMarkerStyle(24);
          gr_Vl->SetMarkerSize(0.8);
	  gr_Vl->SetLineColor(1);
	  gr_Vl->SetFillColor(1);
          gr_Vl->SetMarkerColor(1);
	  gr_Vl->Draw("|>");
	}
	
	TMultiGraph *grdata = new TMultiGraph();
	grdata->SetTitle("SED Data");/*
	//grdata->Add(gr_FERMI,"3");*/
        //grdata->Add(gr_VERITAS,"3");
        
	//grdata->Add(gr0);
	//grdata->Add(gr0_0);
        
        grdata->Add(gr_OVRO);
        grdata->Add(gr_Steward);
        grdata->Add(gr_uvot);
        grdata->Add(gr_xrt);
        //grdata->Add(gr_Vd);
        
	if (l>0){
	  //grdata->Add(gr0a,"|>");
          grdata->Add(gr_Fl,"|>");
	  grdata->Add(gr_Vl,"|>");
	}
	grdata->Add(gr_F);
        grdata->Add(gr_V);
	//host gal
	//grdata->Add(gr10,"C");
        
        
        //add Legend
 	leg = new TLegend(0.13,0.80,0.95,0.87);
        //lowX, lowY, highX, highY
	leg->SetTextSize(0.03);
	leg->SetFillColor(kWhite);
	leg->SetNColumns(6);
	leg->SetMargin(0.2);
        //leg->AddEntry(gr_OVRO,"OVRO","lep");
        leg->AddEntry(gr_Steward,"Steward","lep");
        leg->AddEntry(gr_uvot,"Swift-UVOT","lep");
	leg->AddEntry(gr_xrt,"Swift-XRT","lep");
        leg->AddEntry(gr_F,"Fermi-LAT","lep");
        leg->AddEntry(gr_V,"VERITAS","lep");

	leg->SetEntrySeparation(0.0);
	leg->Draw();
	
	


	
	return 1;
}
