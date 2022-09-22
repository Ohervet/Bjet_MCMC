{
  ifstream file0,file_gal,file_Fermi, file_FermiC, file_VERITAS;
  char str[20],str1[20],str2[50],info[10];
  int i, count, n0=46,n0_1=7, n1a_2=43,n1a_3=8,n2=32,n2_1=9,n2_2=4,n2_3=6,n2_4=8,n2_5=5,nb0,nb1,nb2,nb3, nfs = 50, nVs = 50, nS = 50;
  Double_t A,B,E[n0],F[n0],C[n0],D[n0],deltaEmin[n0],deltaEmax[n0],deltaFmin[n0],deltaFmax[n0],E1[n2],F1[n2],C1[n2],D1[n2],deltaEmin1a[n2],deltaEmax1[n2],deltaF1[n2],deltaF1max[n2],E1uplim[n2],F1uplim[n2],dEmin1auplim[n2],dEmax1uplim[n2],dF1uplim[n2],
  ,Euplim[n2],Fuplim[n2],dEminuplim[n2],dEmaxuplim[n2],dFuplim[n2];
  Double_t E1_1[n2_1],F1_1[n2_1],E1_2[n2_2],F1_2[n2_2],E1_3[n2_3],F1_3[n2_3],E1_4[n2_4],F1_4[n2_4],E1_5[n2_5],F1_5[n2_5],E_1[n0_1],F_1[n0_1],deltaEmin_1[n0_1],deltaEmax_1[n0_1],deltaF_1[n0_1],deltaFmax_1[n0_1],E_2[n1a_3],F_2[n1a_3];
  int nw = 4;
  Double_t E_wise[nw], F_wise[nw],deltaEmin_wise[nw],deltaEmax_wise[nw],deltaFmin_wise[nw],deltaFmax_wise[nw];
  Double_t E_wisel[1], F_wisel[1],deltaF_wisel[1];
  int ns = 3;
  Double_t E_caha[ns], F_caha[ns],deltaEmin_caha[ns],deltaEmax_caha[ns],deltaFmin_caha[ns],deltaFmax_caha[ns];
  int nf = 3;
  Double_t E_flwo[nf], F_flwo[nf],deltaEmin_flwo[nf],deltaEmax_flwo[nf],deltaFmin_flwo[nf],deltaFmax_flwo[nf];
  int nx = 10;
  Double_t E_xrt[nx], F_xrt[nx],deltaEmin_xrt[nx],deltaEmax_xrt[nx],deltaFmin_xrt[nx],deltaFmax_xrt[nx];
  int nb = 8;
  Double_t E_bat[nb], F_bat[nb],deltaEmin_bat[nb],deltaEmax_bat[nb],deltaFmin_bat[nb],deltaFmax_bat[nb];
  Double_t E_batl[1], F_batl[1],deltaF_batl[1];
  Double_t Freq_F[nfs], Flux_F[nfs], DFmin_F[nfs], DFmax_F[nfs];
  Double_t Freq_Vs[nVs], Flux_Vs[nVs], DFmin_Vs[nVs], DFmax_Vs[nVs];
  //Fermi spectrum Caitlin
  Double_t Freq_FCs[nS], Flux_FCs[nS], DFmin_FCs[nS], DFmax_FCs[nS];
  //Fermi data Caitlin
  int n_FC = 6;
  Double_t E_FC[n_FC], F_FC[n_FC],deltaEmin_FC[n_FC],deltaEmax_FC[n_FC],deltaFmin_FC[n_FC],deltaFmax_FC[n_FC];
  int n_F = 6, s_F = 32;
  Double_t E_F[n_F], F_F[n_F],deltaEmin_F[n_F],deltaEmax_F[n_F],deltaFmin_F[n_F],deltaFmax_F[n_F];
  int n_Fl = 1, s_Fl = 2;
  Double_t E_Fl[n_Fl], F_Fl[n_Fl],deltaF_Fl[n_Fl];
  int n_V = 11, s_V = 38;

  
  const double planck = 4.135667517e-15; // h in [eV s]

   file0.open("data/3C66A/SED_3C66A.dat");
   file_Fermi.open("data/HESS_J1943+213/Fermi_spectrum.dat");
   file_FermiC.open("data/3C66A/Fermi_Caitlin_spectrum.dat");
   file_VERITAS.open("data/HESS_J1943+213/VERITAS_spectrum.dat");
   
   
   double     nb_g=0;
   const int   ng = 1298;
   double    Freq, FreqLow, FreqHigh, Flux, FluxLow, FluxHigh, DFlux, DFluxLow, DFluxHigh;
   Double_t   freq_g1[ng],flux_g1[ng],freq_g[ng],flux_g[ng],fflux_g1[ng],fflux_g[ng];
   
   file_gal.open("data/HGS_15.dat"); // host galaxy
   
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
        
      //Fermi Caitlin spectrum
      for(int i=0;i<nS;i++)
      {
	file_FermiC>>Freq>>Flux>>DFluxLow>>DFluxHigh;

            Freq_FCs[i]=log10(Freq*1.0e6/planck);;
 	    Flux_FCs[i]=log10(Flux);
	    DFmin_FCs[i] = Flux_FCs[i] - log10((Flux-DFluxLow));
	    DFmax_FCs[i] = log10((Flux+DFluxHigh)) - Flux_FCs[i];
   
            if (!file_FermiC.good()) break;

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

	
	
        if(nb0>0){
          TGraph *gls = new TGraph(nb0,C,D);
 	  gls->SetMarkerStyle(1);
          gls->SetMarkerColor(15);
          gls->Draw("*");
        }
        
        
        
	

	Double_t minHz, maxHz,minF,maxF;
	minHz=22.;
	maxHz=26.;    
	minF=-11.5.;
	maxF=-9.5;

 	nb1=k;
	//display errorbars
	TGraphAsymmErrors *gr0 = new TGraphAsymmErrors(4,E,F,deltaEmin,deltaEmax,deltaFmin,deltaFmax);
	//gr1->SetMarkerStyle(20);
	gr0->SetMarkerColor(1);
	gr0->SetTitle("");
	gr0->GetYaxis()->SetTitle("log (#nu F#nu [erg cm^{-2} s^{-1} ] )");
	gr0->GetYaxis()->SetTitleSize(0.04);
	gr0->GetYaxis()->SetTitleOffset(1.3);
	gr0->GetYaxis()->SetLabelSize(0.04);
	gr0->GetXaxis()->SetTitle("log (#nu [Hz])");
	gr0->GetXaxis()->SetTitleSize(0.04);
	gr0->GetXaxis()->SetTitleOffset(1.2);
	gr0->GetXaxis()->SetLabelSize(0.04);
	gr0->GetXaxis()->SetLimits(minHz,maxHz);
	gr0->SetMinimum(minF);
	gr0->SetMaximum(maxF);
        gr0->Draw("AP");

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
	TGraph *gr0_0 = new TGraph(5,E,F);
	gr0_0->SetMarkerStyle(20);
	gr0_0->SetMarkerColor(1);
	gr0_0->Draw("P");
	
	//display upper limits
        /*
	if (l>0){
	  TGraphAsymmErrors *gr0a = new TGraphAsymmErrors(l,Euplim,Fuplim,0,0,dFuplim,0);
	  gr0a->SetMarkerStyle(20);
	  gr0a->SetMarkerColor(1);
	  //gr0a->Draw("|>");
	}*/

	// Host galaxy
        if(nb_g>0){
          gr10 = new TGraph(nb_g,freq_g,fflux_g);
	  gr10->SetLineColor(2);
          gr10->SetLineStyle(1);
	  gr10->Draw("C");
        }

        
	//WISE
        count = 4;
	for(i=0;i<nw;i++){
	  E_wise[i] = E[count+i];
	  F_wise[i] = F[count+i];
	  deltaEmin_wise[i] = deltaEmin[count+i];
	  deltaEmax_wise[i] = deltaEmax[count+i];
	  deltaFmin_wise[i] = deltaFmin[count+i];
	  deltaFmax_wise[i] = deltaFmax[count+i];
	  }

	gr_wise = new TGraphAsymmErrors(nw,E_wise,F_wise,deltaEmin_wise,deltaEmax_wise,deltaFmin_wise,deltaFmax_wise);
	gr_wise->SetMarkerColor(2);
	gr_wise->SetLineColor(2);
	gr_wise->SetMarkerStyle(21);
	gr_wise->SetMarkerColor(2);
	gr_wise->Draw("P");        
        
	if (l>0){
	  for(i=0;i<1;i++){
	  E_wisel[i] = Euplim[0];
	  F_wisel[i] = Fuplim[0];
	  deltaF_wisel[i] = dFuplim[0];
	  std::cout<<i<<"\t"<<E_wisel[i]<<"\t"<<F_wisel[i]<<" \t"<<deltaF_wisel[i]<<"\t"<<std::endl;
	}
	  TGraphAsymmErrors *gr_wise_l = new TGraphAsymmErrors(1,E_wisel,F_wisel,0,0,deltaF_wisel,0);
	  gr_wise_l->SetMarkerStyle(20);
	  gr_wise_l->SetLineColor(2);
	  gr_wise_l->SetFillColor(2);
          gr_wise_l->SetMarkerColor(2);
	  gr_wise_l->Draw("|>");
	}
        
        
	//CAHA
        count += nw;
	for(i=0;i<ns;i++){
	  E_caha[i] = E[count+i];
	  F_caha[i] = F[count+i];
	  deltaEmin_caha[i] = deltaEmin[count+i];
	  deltaEmax_caha[i] = deltaEmax[count+i];
	  deltaFmin_caha[i] = deltaFmin[count+i];
	  deltaFmax_caha[i] = deltaFmax[count+i];
	  }

	gr_caha = new TGraphAsymmErrors(ns,E_caha,F_caha,deltaEmin_caha,deltaEmax_caha,deltaFmin_caha,deltaFmax_caha);
        gr_caha->SetLineColor(kOrange+2);
	gr_caha->SetMarkerStyle(33);
        gr_caha->SetMarkerSize(1.5);
	gr_caha->SetMarkerColor(kOrange+2);
	gr_caha->Draw("P");      

        
	//FLWO
        count += ns;
	for(i=0;i<nf;i++){
	  E_flwo[i] = E[count+i];
	  F_flwo[i] = F[count+i];
	  deltaEmin_flwo[i] = deltaEmin[count+i];
	  deltaEmax_flwo[i] = deltaEmax[count+i];
	  deltaFmin_flwo[i] = deltaFmin[count+i];
	  deltaFmax_flwo[i] = deltaFmax[count+i];
	  }

	gr_flwo = new TGraphAsymmErrors(nf,E_flwo,F_flwo,deltaEmin_flwo,deltaEmax_flwo,deltaFmin_flwo,deltaFmax_flwo);
        gr_flwo->SetLineColor(kYellow+2);
	gr_flwo->SetMarkerStyle(34);
	gr_flwo->SetMarkerColor(kYellow+2);
	gr_flwo->Draw("P");    
        
        
	//XRT
        count += nf;
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
        
        
 	//BAT
        count += nx;
	for(i=0;i<nb;i++){
	  E_bat[i] = E[count+i];
	  F_bat[i] = F[count+i];
	  deltaEmin_bat[i] = deltaEmin[count+i];
	  deltaEmax_bat[i] = deltaEmax[count+i];
	  deltaFmin_bat[i] = deltaFmin[count+i];
	  deltaFmax_bat[i] = deltaFmax[count+i];
	  }

	gr_bat = new TGraphAsymmErrors(7,E_bat,F_bat,deltaEmin_bat,deltaEmax_bat,deltaFmin_bat,deltaFmax_bat);
        gr_bat->SetLineColor(kCyan+2);
	gr_bat->SetMarkerStyle(22);
	gr_bat->SetMarkerColor(kCyan+2);
	gr_bat->Draw("P");     
        
        E_batl[0] = Euplim[1];
        F_batl[0] = Fuplim[1];
        deltaF_batl[0] = dFuplim[1];

        TGraphAsymmErrors *gr_bat_l = new TGraphAsymmErrors(1,E_batl,F_batl,0,0,deltaF_batl,0);
        gr_bat_l->SetMarkerStyle(20);
        gr_bat_l->SetLineColor(kCyan+2);
        gr_bat_l->SetFillColor(kCyan+2);
        gr_bat_l->SetMarkerColor(kCyan+2);
        gr_bat_l->Draw("|>");
     
        
        
        
	//Fermi spectrum
	gr_FERMI = new TGraphAsymmErrors(nfs,Freq_F,Flux_F,0,0,DFmin_F,DFmax_F);
	gr_FERMI->SetFillColor(kBlue-7);
	//gr_FERMI->Draw("3");
        
        //Fermi Caitlin spectrum
	gr_FERMIC = new TGraphAsymmErrors(nS,Freq_FCs,Flux_FCs,0,0,DFmin_FCs,DFmax_FCs);
	gr_FERMIC->SetFillColor(kBlue-7);
	gr_FERMIC->Draw("3");
        
	//VERITAS spectrum
	gr_VERITAS = new TGraphAsymmErrors(nVs,Freq_Vs,Flux_Vs,0,0,DFmin_Vs,DFmax_Vs);
	gr_VERITAS->SetFillColor(kMagenta-7);
	//gr_VERITAS->Draw("3");
	
	//Fermi data
        count += nb;
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
	gr_F->SetLineColor(kBlue+2);
        gr_F->SetMarkerColor(kBlue+2);
	//gr_F->Draw("P");
        
	if (l>0){
	  for(i=0;i<n_Fl;i++){
	  E_Fl[i] = Euplim[s_Fl+i];
	  F_Fl[i] = Fuplim[s_Fl+i];
	  deltaF_Fl[i] = dFuplim[s_Fl+i];
	  std::cout<<i<<"\t"<<E_Fl[i]<<"\t"<<F_Fl[i]<<" \t"<<deltaF_Fl[i]<<"\t"<<std::endl;
	}
	  TGraphAsymmErrors *gr_F_l = new TGraphAsymmErrors(n_Fl,E_Fl,F_Fl,0,0,deltaF_Fl,0);
	  gr_F_l->SetMarkerStyle(24);
          gr_F_l->SetMarkerSize(0.8);
	  gr_F_l->SetLineColor(kBlue+2);
	  gr_F_l->SetFillColor(kBlue+2);
          gr_F_l->SetMarkerColor(kBlue+2);
	  //gr_F_l->Draw("|>");
	}
        
	//Fermi data Caitlin
        count += n_F;
	for(i=0;i<n_FC;i++){
	  E_FC[i] = E[count+i];
	  F_FC[i] = F[count+i];
	  deltaEmin_FC[i] = deltaEmin[count+i];
	  deltaEmax_FC[i] = deltaEmax[count+i];
	  deltaFmin_FC[i] = deltaFmin[count+i];
	  deltaFmax_FC[i] = deltaFmax[count+i];
	}
	gr_FC = new TGraphAsymmErrors(n_FC,E_FC,F_FC,deltaEmin_FC,deltaEmax_FC,deltaFmin_FC,deltaFmax_FC);
	gr_FC->SetMarkerStyle(24);
        gr_FC->SetMarkerSize(0.8);
        gr_FC->SetMarkerColor(kMagenta+2);
	gr_FC->SetLineColor(kMagenta+2);
	gr_FC->Draw("P");
        
	 
	  
	
	TMultiGraph *grdata = new TMultiGraph();
	grdata->SetTitle("SED Data");
	grdata->Add(gr_FERMI,"3");
        grdata->Add(gr_VERITAS,"3");
	grdata->Add(gr0);
	grdata->Add(gr0_0);
        grdata->Add(gr_wise);
        grdata->Add(gr_caha);
        grdata->Add(gr_flwo);
        grdata->Add(gr_xrt);
        grdata->Add(gr_bat);
	if (l>0){
	  //grdata->Add(gr0a,"|>");
          grdata->Add(gr_wise_l,"|>");
          //grdata->Add(gr_bat_l,"|>");
	  grdata->Add(gr_F_l,"|>");
	}
	grdata->Add(gr_F);
	//host gal
	grdata->Add(gr10,"C");
	
	


	
	return 1;
}
