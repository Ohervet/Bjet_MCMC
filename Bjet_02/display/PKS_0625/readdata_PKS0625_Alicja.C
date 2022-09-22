{
  ifstream file0,file1,file_gal,file_hess,file_fermi, file_XRT, file_XRT2, file_XRT3, file_XRT4, file_XMM, file_Suzaku, file15;
  char str[20],str1[20],str2[40],info[10];
  int i, count, n0=38, n1=163, nw = 4, na =1,nu =6, nf =6,nh = 8, nhb = 46, nfb = 50, n_Xb = 50, nbFu = 59;
  int n0_0=58,n0_1=7, n1a_2=43,n1a_3=8,n2=32,n2_1=9,n2_2=4,n2_3=6,n2_4=8,n2_5=5,nb0,nb1,nb2,nb3,,nb15 = 0;
  Double_t E[n0],F[n0],deltaEmin[n0],deltaEmax[n0],deltaFmin[n0],deltaFmax[n0],Euplim[n0],Fuplim[n0],dEminuplim[n0],dEmaxuplim[n0],dFuplim[n0];
  Double_t E1[n1],F1[n1],deltaEmin1[n1],deltaEmax1[n1],deltaFmin1[n1],deltaFmax1[n1],Euplim1[n1],Fuplim1[n1],dEminuplim1[n1],dEmaxuplim1[n1],dFuplim1[n1];
  //Double_t E1_1[n2_1],F1_1[n2_1],E1_2[n2_2],F1_2[n2_2],E1_3[n2_3],F1_3[n2_3],E1_4[n2_4],F1_4[n2_4],E1_5[n2_5],F1_5[n2_5],E_1[n0_1],F_1[n0_1],deltaEmin_1[n0_1],deltaEmax_1[n0_1],deltaF_1[n0_1],deltaFmax_1[n0_1],E_2[n1a_3],F_2[n1a_3];
  Double_t E_Rcore[na], F_Rcore[na];
  Double_t E_wise[nw], F_wise[nw],deltaEmin_wise[nw],deltaEmax_wise[nw],deltaFmin_wise[nw],deltaFmax_wise[nw];
  Double_t E_atom[na], F_atom[na],deltaEmin_atom[na],deltaEmax_atom[na],deltaFmin_atom[na],deltaFmax_atom[na];
  Double_t E_uvot[nu], F_uvot[nu],deltaEmin_uvot[nu],deltaEmax_uvot[nu],deltaFmin_uvot[nu],deltaFmax_uvot[nu];
  int nx1 = 8;
  Double_t E_xrt1[nx1], F_xrt1[nx1],deltaEmin_xrt1[nx1],deltaEmax_xrt1[nx1],deltaFmin_xrt1[nx1],deltaFmax_xrt1[nx1];
  int nx2 = 56;
  Double_t E_xrt2[nx2], F_xrt2[nx2],deltaEmin_xrt2[nx2],deltaEmax_xrt2[nx2],deltaFmin_xrt2[nx2],deltaFmax_xrt2[nx2];
  int nxl = 26;
  Double_t E_xl[nxl], F_xl[nxl],deltaEmin_xl[nxl],deltaEmax_xl[nxl],deltaFmin_xl[nxl],deltaFmax_xl[nxl];
  int nxh = 21;
  Double_t E_xh[nxh], F_xh[nxh],deltaEmin_xh[nxh],deltaEmax_xh[nxh],deltaFmin_xh[nxh],deltaFmax_xh[nxh];
  Double_t E_fermi[nf], F_fermi[nf],deltaEmin_fermi[nf],deltaEmax_fermi[nf],deltaFmin_fermi[nf],deltaFmax_fermi[nf];
  Double_t E_hess[nh], F_hess[nh],deltaEmin_hess[nh],deltaEmax_hess[nh],deltaFmin_hess[nh],deltaFmax_hess[nh];
  Double_t Freq_H[nhb], Flux_H[nhb], DFmin_H[nhb], DFmax_H[nhb];
  Double_t Freq_F[nfb], Flux_F[nfb], DFmin_F[nfb], DFmax_F[nfb];
  Double_t Freq_Xb[n_Xb], Flux_Xb[n_Xb],DFmin_Xb[n_Xb],DFmax_Xb[n_Xb];
  Double_t freq15[nbFu], fflux15[nbFu];
  Double_t Euplimh[2],Fuplimh[2];
  
  const double planck = 4.135667517e-15; // h in [eV s]

   file0.open("data/PKS_0625-354/PKS_0625_Alicja_archive.dat");
   file1.open("data/PKS_0625-354/PKS_0625_Alicja_sed.dat");
   file_hess.open("data/PKS_0625-354/hess_butterfly.dat");
   //file_fermi.open("data/PKS_0625-354/fermi_butterfly_2.dat");
   file_fermi.open("data/PKS_0625-354/fermi_PowerLaw_01_2018.dat");
   file_XRT.open("data/PKS_0625-354/XRT_butterfly_1e.dat");
   file_XRT2.open("data/PKS_0625-354/XRT_butterfly_2e.dat");
   file_XRT3.open("data/PKS_0625-354/XRT_butterfly_3e.dat");
   file_XRT4.open("data/PKS_0625-354/XRT_butterfly_4e.dat");
   file_XMM.open("data/PKS_0625-354/XMM_PowerLaw.dat");
   file_Suzaku.open("data/PKS_0625-354/Suzaku_PowerLaw.dat");
   //file_XRT.open("data/PKS_0625-354/XRT_butterfly_highstate.dat");
   //file_gal.open("data/PKS_0625-354/PKS_0625_gal.dat"); // host galaxy
   file_gal.open("data/PKS_0625-354/host_gal_PKS_0625_2011.dat");
   file15.open("data/PKS_0625-354/Fukazawa_model.dat"); // model from Fukazawa 2015
   
   
   
   double     nb_g=0;
   const int   ng = 1304;
   double    Freq, FreqLow, FreqHigh, Flux, FluxLow, FluxHigh, DFlux;
   Double_t   freq_g1[ng],flux_g1[ng],freq_g[ng],flux_g[ng],fflux_g1[ng],fflux_g[ng];
   

   
   
   
   
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
       
       
       


    if(!file0 || !file1){
       cout<<"data file not found"<<endl;
       return 0;
    }
 
     int j=0, k=0, l=0, l1=0;
     
     for(int i=0;i<=14;i++){
        file0>>str;
	file1>>str1;
     }

      std::cout<<"i\tFreq\tFreqLow\tFreqHigh\tFlux\tFluxLow\tFluxHigh\tData"<<std::endl;
      
      //archival data
      for(int i=3;i<n0;i++)
      {
	file0>>Freq>>FreqLow>>FreqHigh>>Flux>>FluxLow>>FluxHigh>>str2;


            E[i]=log10(Freq);
 	    F[i]=log10(Flux);
	    deltaEmin[i]=E[i]-log10(Freq-FreqLow);
	    deltaEmax[i]=-E[i]+log10(Freq+FreqHigh);
	
 	    deltaFmin[i]=F[i]-log10(Flux-FluxLow);
	    deltaFmax[i]=-F[i]+log10(Flux+FluxHigh);
	    if (deltaFmax[i] == 0){
	      Euplim[l] = E[i];
	      Fuplim[l] = F[i];
	      dEminuplim[l]= deltaEmin[i];
	      dEmaxuplim[l]= deltaEmax[i];
	      dFuplim[l]= 0.2;
	      l+=1;
	    }
   
	  std::cout<<i<<"\t"<<E[i]<<"\t"<<deltaEmin[i]<<" \t"<<deltaEmax[i]<<"\t"<<F[i]<<"\t"<<deltaFmin[i]<<"\t"<<deltaFmax[i]<<"\t"<<str2<<std::endl;

            if (!file0.good()) break;

        } 
         
      //sed data
      for(int i=3;i<n1;i++)
      {
	file1>>Freq>>FreqLow>>FreqHigh>>Flux>>FluxLow>>FluxHigh>>str2;


            E1[i]=log10(Freq);
 	    F1[i]=log10(Flux);
	    deltaEmin1[i]=E1[i]-log10(Freq-FreqLow);
	    deltaEmax1[i]=-E1[i]+log10(Freq+FreqHigh);
	
 	    deltaFmin1[i]=F1[i]-log10(Flux-FluxLow);
	    deltaFmax1[i]=-F1[i]+log10(Flux+FluxHigh);
	    if (deltaFmax1[i] == 0){
	      Euplim1[l1] = E1[i];
	      Fuplim1[l1] = F1[i];
	      dEminuplim1[l1]= deltaEmin1[i];
	      dEmaxuplim1[l1]= deltaEmax1[i];
	      dFuplim1[l1]= 0.2;
	      l1+=1;
	      //std::cout<<Euplim1[l]<<std::endl;
	    }
   
	  std::cout<<i<<"\t"<<E1[i]<<"\t"<<deltaEmin1[i]<<" \t"<<deltaEmax1[i]<<"\t"<<F1[i]<<"\t"<<deltaFmin1[i]<<"\t"<<deltaFmax1[i]<<"\t"<<str2<<std::endl;

            if (!file1.good()) break;

        }
        
        
       //HESS butterfly
      for(int i=0;i<nhb;i++)
      {
	file_hess>>Freq>>FluxLow>>FluxHigh;


            Freq_H[i]=log10(Freq*1.0e12/planck);
 	    Flux_H[i]=log10((FluxHigh+FluxLow)/2. *Freq*1e-4* Freq*1.0e12 * 1.60217e-12 );
	    DFmin_H[i] = Flux_H[i] - log10(FluxLow *Freq*1e-4* Freq*1.0e12 * 1.60217e-12);//1.60217657
	    DFmax_H[i] = log10(FluxHigh *Freq*1e-4* Freq*1.0e12 * 1.60217e-12) - Flux_H[i];
	    
	    //std::cout<<i<<"\t"<<Freq_H[i]<<"\t"<<Flux_H[i]<<" \t"<<DFmin_H[i]<<"\t"<<DFmax_H[i]<<std::endl;
   
            if (!file_hess.good()) break;

        }
        
       //Fermi butterfly
      for(int i=0;i<nfb;i++)
      {
	//file_fermi>>Freq>>Flux>>DFlux;
        file_fermi>>Freq>>Flux>>FluxLow>>FluxHigh;

            Freq_F[i]=log10(Freq*1.0e6/planck);
 	    Flux_F[i]=log10(Flux);
	    //DFmin_F[i] = Flux_F[i] - log10(Flux-DFlux);
	    //DFmax_F[i] = log10(Flux+DFlux) - Flux_F[i];
            DFmin_F[i] = Flux_F[i] - log10(Flux-FluxLow);
	    DFmax_F[i] = log10(FluxHigh+Flux)-Flux_F[i];
	    
	    //std::cout<<i<<"\t"<<Freq_F[i]<<"\t"<<Flux_F[i]<<" \t"<<DFmin_F[i]<<"\t"<<DFmax_F[i]<<std::endl;
   
            if (!file_fermi.good()) break;

        }


	
        if(nb0>0){
          TGraph *gls = new TGraph(nb0,C,D);
 	  gls->SetMarkerStyle(1);
          gls->SetMarkerColor(15);
          gls->Draw("*");
        }
        

	

	Double_t minHz, maxHz,minF,maxF;
	minHz=7.5;
	maxHz=28;    
	minF=-14.5;
	maxF=-9.5;

	
//----------archival data----------//
	
	TCanvas *c1 = new TCanvas("c1","c1",700,500);
//	gStyle->SetLegendTextSize(15);
	
 	nb1=k;
	//display errorbars
	TGraphAsymmErrors *gr0 = new TGraphAsymmErrors(n0,E,F,deltaEmin,deltaEmax,deltaFmin,deltaFmax);
	//gr1->SetMarkerStyle(20);
	gr0->SetLineColor(16);
	gr0->SetMarkerStyle(20);
	gr0->SetMarkerColor(16);
	gr0->SetTitle("");
	gr0->GetYaxis()->SetTitle("log (#nu F_{#nu} [erg cm^{-2} s^{-1} ] )");
	gr0->GetYaxis()->SetTitleSize(0.04);
	gr0->GetYaxis()->SetTitleOffset(1.4);
	gr0->GetYaxis()->SetLabelSize(0.04);
	gr0->GetXaxis()->SetTitle("log (#nu [Hz])");
	gr0->GetXaxis()->SetTitleSize(0.04);
	gr0->GetXaxis()->SetTitleOffset(1.3);
	gr0->GetXaxis()->SetLabelSize(0.04);
	gr0->SetMinimum(minF);
	gr0->SetMaximum(maxF);
	gr0->GetXaxis()->SetLimits(minHz,maxHz);
	gr0->Draw("AP");
	
	//display dots
	/*
	TGraph *gr0_0 = new TGraph(n0,E,F);
	gr0_0->SetMarkerStyle(20);
	gr0_0->SetMarkerColor(16);
	gr0_0->Draw("P");*/
	
	//display upper limits
	if (l>0){
	  gr0a = new TGraphAsymmErrors(l,Euplim,Fuplim,0,0,dFuplim,0);
	  gr0a->SetMarkerStyle(20);
	  gr0a->SetMarkerColor(1);
	  gr0a->Draw("|>");
	}
	
	
	
//----------SED data----------//


//******************************************************************
// Fukazawa 2015 model
//******************************************************************
       if(!file15){
            cout<<" Fukazawa 2015 model not found !"<<endl;     
       }
       else{
         for(i=0;i<nbFu;i++){ 
           file15>>freq15[i]>>fflux15[i];
           freq15[i] = log10(freq15[i]);
           fflux15[i] = log10(fflux15[i]);
           if (!file15.good()) break;
         }
	nb15=i;
        
        grFu = new TGraph(nb15,freq15,fflux15);  // previous SSC absorbed
        grFu->SetLineColor(16);
        grFu->SetLineStyle(2);
        grFu->SetLineWidth(2);
        grFu->Draw("C");
       }   

/*
	gr1 = new TGraphAsymmErrors(n1,E1,F1,deltaEmin1,deltaEmax1,deltaFmin1,deltaFmax1);
	gr1->SetMarkerColor(1);
	gr1->SetLineColor(1);
	gr1->Draw("P");
	
	//display dots
	TGraph *gr1_0 = new TGraph(n1,E1,F1);
	gr1_0->SetMarkerSize(0.5);
	gr1_0->SetMarkerStyle(20);
	gr1_0->SetMarkerColor(1);
	gr1_0->Draw("P");
*/	
	//display upper limits
	if (l1>0){
	  gr1_l = new TGraphAsymmErrors(l1,Euplim1,Fuplim1,0,0,dFuplim1,0);
	  gr1_l->SetMarkerStyle(20);
	  gr1_l->SetMarkerColor(1);
	  gr1_l->Draw("|>");
	}	
	
	//radio_core
	E_Rcore[0] = E1[3];
	F_Rcore[0] = F1[3];
	gr_Rcore = new TGraph(na,E_Rcore,F_Rcore);
	gr_Rcore->SetMarkerStyle(20);
	gr_Rcore->SetMarkerColor(2);
	//gr_Rcore->Draw("P");

	//WISE
        count = 4;
	for(i=0;i<nw;i++){
	  E_wise[i] = E1[count+i];
	  F_wise[i] = F1[count+i];
	  deltaEmin_wise[i] = deltaEmin1[count+i];
	  deltaEmax_wise[i] = deltaEmax1[count+i];
	  deltaFmin_wise[i] = deltaFmin1[count+i];
	  deltaFmax_wise[i] = deltaFmax1[count+i];
	  }

	gr_wise = new TGraphAsymmErrors(nw,E_wise,F_wise,deltaEmin_wise,deltaEmax_wise,deltaFmin_wise,deltaFmax_wise);
	gr_wise->SetMarkerColor(1);
        gr_wise->SetMarkerStyle(33);
        gr_wise->SetMarkerSize(1.3);
	gr_wise->SetLineColor(1);
	gr_wise->SetMarkerColor(1);
	gr_wise->Draw("P");
	
	//ATOM
        count += nw;
	E_atom[0] = E1[count];
	F_atom[0] = F1[count];
	deltaEmin_atom[0] = deltaEmin1[count];
	deltaEmax_atom[0] = deltaEmax1[count];
	deltaFmin_atom[0] = deltaFmin1[count];
	deltaFmax_atom[0] = deltaFmax1[count];

	gr_atom = new TGraphAsymmErrors(na,E_atom,F_atom,deltaEmin_atom,deltaEmax_atom,deltaFmin_atom,deltaFmax_atom);
	gr_atom->SetMarkerColor(2);
	gr_atom->SetLineColor(2);
	gr_atom->SetMarkerStyle(21);
	gr_atom->SetMarkerColor(2);
	gr_atom->Draw("P");
	
	
	//UVOT
        count += 1;
	for(i=0;i<nu;i++){
	  E_uvot[i] = E1[count+i];
	  F_uvot[i] = F1[count+i];
	  deltaEmin_uvot[i] = deltaEmin1[count+i];
	  deltaEmax_uvot[i] = deltaEmax1[count+i];
	  deltaFmin_uvot[i] = deltaFmin1[count+i];
	  deltaFmax_uvot[i] = deltaFmax1[count+i];
	}

	gr_uvot = new TGraphAsymmErrors(nu,E_uvot,F_uvot,deltaEmin_uvot,deltaEmax_uvot,deltaFmin_uvot,deltaFmax_uvot);
	//gr_uvot->SetMarkerColor(7);
	gr_uvot->SetLineColor(kCyan+2);
	gr_uvot->SetMarkerStyle(20);
	gr_uvot->SetMarkerColor(kCyan+2);
	gr_uvot->Draw("P");

       //XMM butterfly
        for(int i=0;i<n_Xb;i++)
        {
            file_XMM>>Freq>>Flux>>FluxLow>>FluxHigh;

                Freq_Xb[i]=log10(Freq);
                Flux_Xb[i]=log10(Flux);
                DFmin_Xb[i] = Flux_Xb[i] - log10(Flux-FluxLow);
                DFmax_Xb[i] = log10(FluxHigh+Flux)-Flux_Xb[i];
    
                if (!file_XMM.good()) break;
        }
 	gr_XMM = new TGraphAsymmErrors(n_Xb,Freq_Xb,Flux_Xb,0,0,DFmin_Xb,DFmax_Xb);
	gr_XMM->SetFillColor(14);
        gr_XMM->SetLineColor(14);
        gr_XMM->SetLineWidth(2);
	gr_XMM->Draw("3");  
        gr_XMMa = new TGraphAsymmErrors(n_Xb,Freq_Xb,Flux_Xb,0,0,DFmin_Xb,DFmax_Xb);
        gr_XMMa->SetFillColor(14);
        gr_XMMa->SetFillStyle(0);
        gr_XMMa->Draw("3");

       //Suzaku butterfly
        for(int i=0;i<n_Xb;i++)
        {
            file_Suzaku>>Freq>>Flux>>FluxLow>>FluxHigh;

                Freq_Xb[i]=log10(Freq);
                Flux_Xb[i]=log10(Flux);
                DFmin_Xb[i] = Flux_Xb[i] - log10(Flux-FluxLow);
                DFmax_Xb[i] = log10(FluxHigh+Flux)-Flux_Xb[i];
    
                if (!file_Suzaku.good()) break;
        }
 	gr_Suz = new TGraphAsymmErrors(n_Xb,Freq_Xb,Flux_Xb,0,0,DFmin_Xb,DFmax_Xb);
	gr_Suz->SetFillColor(14);
        gr_Suz->SetLineColor(14);
        gr_Suz->SetLineWidth(2);
	gr_Suz->Draw("3");  
        gr_Suza = new TGraphAsymmErrors(n_Xb,Freq_Xb,Flux_Xb,0,0,DFmin_Xb,DFmax_Xb);
        gr_Suza->SetFillColor(14);
        gr_Suza->SetFillStyle(0);
        gr_Suza->Draw("3");
	
	//XRT
        count += nu;
	for(i=0;i<nx1;i++){
	  E_xrt1[i] = E1[count+i];
	  F_xrt1[i] = F1[count+i];
	  deltaEmin_xrt1[i] = deltaEmin1[count+i];
	  deltaEmax_xrt1[i] = deltaEmax1[count+i];
	  deltaFmin_xrt1[i] = deltaFmin1[count+i];
	  deltaFmax_xrt1[i] = deltaFmax1[count+i];
	}

	gr_xrt1 = new TGraphAsymmErrors(nx1,E_xrt1,F_xrt1,deltaEmin_xrt1,deltaEmax_xrt1,deltaFmin_xrt1,deltaFmax_xrt1);
	gr_xrt1->SetLineColor(kRed+1);
	gr_xrt1->SetMarkerStyle(20);
	gr_xrt1->SetMarkerSize(0.5);
	gr_xrt1->SetMarkerColor(kRed+1);
	//gr_xrt1->Draw("P");

        count += nx1;
	for(i=0;i<nx2;i++){
	  E_xrt2[i] = E1[count+i];
	  F_xrt2[i] = F1[count+i];
	  deltaEmin_xrt2[i] = deltaEmin1[count+i];
	  deltaEmax_xrt2[i] = deltaEmax1[count+i];
	  deltaFmin_xrt2[i] = deltaFmin1[count+i];
	  deltaFmax_xrt2[i] = deltaFmax1[count+i];
	}

	gr_xrt2 = new TGraphAsymmErrors(nx2,E_xrt2,F_xrt2,deltaEmin_xrt2,deltaEmax_xrt2,deltaFmin_xrt2,deltaFmax_xrt2);
	gr_xrt2->SetLineColor(kMagenta-7);
	gr_xrt2->SetMarkerStyle(20);
	gr_xrt2->SetMarkerSize(0.5);
	gr_xrt2->SetMarkerColor(kMagenta-7);
	//gr_xrt2->Draw("P");


        //XRT high state
        count += nx2;
	for(i=0;i<nxh;i++){
	  E_xh[i] = E1[count+i];
	  F_xh[i] = F1[count+i];
	  deltaEmin_xh[i] = deltaEmin1[count+i];
	  deltaEmax_xh[i] = deltaEmax1[count+i];
	  deltaFmin_xh[i] = deltaFmin1[count+i];
	  deltaFmax_xh[i] = deltaFmax1[count+i];
	}

	gr_xh = new TGraphAsymmErrors(nxh,E_xh,F_xh,deltaEmin_xh,deltaEmax_xh,deltaFmin_xh,deltaFmax_xh);
	gr_xh->SetLineColor(kRed+1);
        gr_xh->SetMarkerStyle(20);
        gr_xh->SetMarkerSize(0.5);
	gr_xh->SetMarkerColor(kRed+1);
	//gr_xh->Draw("P");
        
	//XRTlow states
        count += nxh;
	for(i=0;i<nxl;i++){
	  E_xl[i] = E1[count+i];
	  F_xl[i] = F1[count+i];
	  deltaEmin_xl[i] = deltaEmin1[count+i];
	  deltaEmax_xl[i] = deltaEmax1[count+i];
	  deltaFmin_xl[i] = deltaFmin1[count+i];
	  deltaFmax_xl[i] = deltaFmax1[count+i];
	}

	gr_xl = new TGraphAsymmErrors(nxl,E_xl,F_xl,deltaEmin_xl,deltaEmax_xl,deltaFmin_xl,deltaFmax_xl);
	gr_xl->SetLineColor(kMagenta-7);
        gr_xl->SetMarkerStyle(20);
        gr_xl->SetMarkerSize(0.5);
	gr_xl->SetMarkerColor(kMagenta-7);
	//gr_xl->Draw("P");
	
	//XRT butterfly
      for(int i=0;i<n_Xb;i++)
      {
	file_XRT>>Freq>>Flux>>FluxLow>>FluxHigh;

            Freq_Xb[i]=log10(Freq);
 	    Flux_Xb[i]=log10(Flux);
	    DFmin_Xb[i] = Flux_Xb[i] - log10(Flux-FluxLow);
	    DFmax_Xb[i] = log10(FluxHigh+Flux)-Flux_Xb[i];
   
            if (!file_XRT.good()) break;
      }
        gr_Xba = new TGraphAsymmErrors(n_Xb,Freq_Xb,Flux_Xb,0,0,DFmin_Xb,DFmax_Xb);
        gr_Xba->SetFillColor(kRed+1);
        gr_Xba->SetLineColor(kRed+1);
        gr_Xba->SetLineWidth(2);
        gr_Xba->SetFillStyle(0); 
        gr_Xba->Draw("3");
        
       //XRT butterfly 2
        for(int i=0;i<n_Xb;i++)
        {
            file_XRT2>>Freq>>Flux>>FluxLow>>FluxHigh;

                Freq_Xb[i]=log10(Freq);
                Flux_Xb[i]=log10(Flux);
                DFmin_Xb[i] = Flux_Xb[i] - log10(Flux-FluxLow);
                DFmax_Xb[i] = log10(FluxHigh+Flux)-Flux_Xb[i];
    
                if (!file_XRT2.good()) break;
        }
        gr_Xb2a = new TGraphAsymmErrors(n_Xb,Freq_Xb,Flux_Xb,0,0,DFmin_Xb,DFmax_Xb);
        gr_Xb2a->SetFillColor(kMagenta-7);
        gr_Xb2a->SetLineColor(kMagenta-7);
        gr_Xb2a->SetLineWidth(2);
        gr_Xb2a->SetFillStyle(0);
        gr_Xb2a->Draw("3");

       //XRT butterfly 3
        for(int i=0;i<n_Xb;i++)
        {
            file_XRT3>>Freq>>Flux>>FluxLow>>FluxHigh;

                Freq_Xb[i]=log10(Freq);
                Flux_Xb[i]=log10(Flux);
                DFmin_Xb[i] = Flux_Xb[i] - log10(Flux-FluxLow);
                DFmax_Xb[i] = log10(FluxHigh+Flux)-Flux_Xb[i];
    
                if (!file_XRT3.good()) break;
        }
 	gr_Xb3 = new TGraphAsymmErrors(n_Xb,Freq_Xb,Flux_Xb,0,0,DFmin_Xb,DFmax_Xb);
	gr_Xb3->SetFillColor(kMagenta-7);
        //gr_Xb3->SetFillStyle(0);
        gr_Xb3->SetLineColor(kMagenta+1);
        gr_Xb3->SetLineWidth(2);
	//gr_Xb3->Draw("3");
        gr_Xb3a = new TGraphAsymmErrors(n_Xb,Freq_Xb,Flux_Xb,0,0,DFmin_Xb,DFmax_Xb);
        gr_Xb3a->SetFillColor(kMagenta-7);
        gr_Xb3a->SetFillStyle(0); 

       //XRT butterfly 4
        for(int i=0;i<n_Xb;i++)
        {
            file_XRT4>>Freq>>Flux>>FluxLow>>FluxHigh;

                Freq_Xb[i]=log10(Freq);
                Flux_Xb[i]=log10(Flux);
                DFmin_Xb[i] = Flux_Xb[i] - log10(Flux-FluxLow);
                DFmax_Xb[i] = log10(FluxHigh+Flux)-Flux_Xb[i];
    
                if (!file_XRT4.good()) break;
        }
 	gr_Xb4 = new TGraphAsymmErrors(n_Xb,Freq_Xb,Flux_Xb,0,0,DFmin_Xb,DFmax_Xb);
	gr_Xb4->SetFillColor(kMagenta-7);
        gr_Xb4->SetLineColor(kMagenta+1);
        gr_Xb4->SetLineWidth(2);
	//gr_Xb4->Draw("3");    
        gr_Xb4a = new TGraphAsymmErrors(n_Xb,Freq_Xb,Flux_Xb,0,0,DFmin_Xb,DFmax_Xb);
        gr_Xb4a->SetFillColor(kMagenta-7);
        gr_Xb4a->SetFillStyle(0);
        gr_Xb3a->Draw("3");
        gr_Xb4a->Draw("3");  
   
        
	
	//Fermi butterfly
	gr_FERMI = new TGraphAsymmErrors(nfb,Freq_F,Flux_F,0,0,DFmin_F,DFmax_F);
	gr_FERMI->SetFillColor(kBlue-7);
	gr_FERMI->Draw("3");
	
	//Fermi data
        count += nxl;
	for(i=0;i<nf;i++){
	  E_fermi[i] = E1[count+i];
	  F_fermi[i] = F1[count+i];
	  deltaEmin_fermi[i] = deltaEmin1[count+i];
	  deltaEmax_fermi[i] = deltaEmax1[count+i];
	  deltaFmin_fermi[i] = deltaFmin1[count+i];
	  deltaFmax_fermi[i] = deltaFmax1[count+i];
	}

	gr_fermi = new TGraphAsymmErrors(nf,E_fermi,F_fermi,deltaEmin_fermi,deltaEmax_fermi,deltaFmin_fermi,deltaFmax_fermi);
	gr_fermi->SetLineColor(kBlue+1);
	gr_fermi->Draw("P");

	gr_fermi_0 = new TGraph(nf,E_fermi,F_fermi);
	gr_fermi_0->SetMarkerStyle(20);
	gr_fermi_0->SetMarkerColor(kBlue+1);
	gr_fermi_0->Draw("P");	
	/*
	//display upper limits
	if (l1>0){
	  gr1_l = new TGraphAsymmErrors(1,Euplim1,Fuplim1,0,0,dFuplim1,0);
	  gr1_l->SetMarkerStyle(20);
	  gr1_l->SetMarkerColor(kBlue+1);
	  gr1_l->SetLineColor(kBlue+2);
	  gr1_l->SetFillColor(kBlue+1);
	  gr1_l->Draw("|>");
	}	
	*/

	
	
	
	//HESS butterfly
	gr_HESS = new TGraphAsymmErrors(nhb,Freq_H,Flux_H,0,0,DFmin_H,DFmax_H);
	gr_HESS->SetFillColor(kGreen-7);
	gr_HESS->Draw("3");
	
	//HESS data
        count += nf;
	for(i=0;i<nh;i++){
	  E_hess[i] = E1[count+i];
	  F_hess[i] = F1[count+i];
	  deltaEmin_hess[i] = deltaEmin1[count+i];
	  deltaEmax_hess[i] = deltaEmax1[count+i];
	  deltaFmin_hess[i] = deltaFmin1[count+i];
	  deltaFmax_hess[i] = deltaFmax1[count+i];
	}

	gr_hess = new TGraphAsymmErrors(nh,E_hess,F_hess,deltaEmin_hess,deltaEmax_hess,deltaFmin_hess,deltaFmax_hess);
	gr_hess->SetLineColor(kGreen+2);
	gr_hess->Draw("P");

	gr_hess_0 = new TGraph(nh,E_hess,F_hess);
	gr_hess_0->SetMarkerStyle(20);
	gr_hess_0->SetMarkerColor(kGreen+2);
	gr_hess_0->Draw("P");	
	
	//display upper limits
	if (l1>0){
	  for(i=0;i<2;i++){
	    Euplimh[i] = Euplim1[i];
	    Fuplimh[i] = Fuplim1[i];
	  }
	  grh_l = new TGraphAsymmErrors(2,Euplimh,Fuplimh,0,0,dFuplim1,0);
	  grh_l->SetMarkerStyle(20);
	  grh_l->SetMarkerColor(kGreen+2);
	  grh_l->SetLineColor(kGreen+2);
	  grh_l->SetFillColor(kGreen+2);
	  grh_l->Draw("|>");
	}	
	
	
	
	// Host galaxy
    
        if(nb_g>0){
          gr10 = new TGraph(nb_g,freq_g,fflux_g);
	  gr10->SetLineColor(1);
          gr10->SetLineStyle(1);
	  gr10->Draw("C");
        }
        
       
       
       
	// add legend
	leg = new TLegend(0.15,0.76,0.9,0.87);
	leg->SetTextSize(0.03);
	leg->SetFillColor(kWhite);
	leg->SetNColumns(3);
	leg->SetMargin(0.2);
	leg->AddEntry(gr_wise,"WISE","lep");
	
	leg->AddEntry(gr_Xba,"Swift-XRT high state","f");
        leg->AddEntry(gr_FERMI,"Fermi-LAT","f");
        leg->AddEntry(gr_atom,"ATOM","lep");
        
	
	leg->AddEntry(gr_Xb2a,"Swift-XRT low states","f");
        leg->AddEntry(gr_HESS,"H.E.S.S.","f");
	leg->AddEntry(gr_uvot,"Swift-UVOT","lep");
        
	leg->AddEntry(gr0,"Archival data","lep");
	leg->AddEntry(gr10,"Host galaxy","l");
	leg->SetEntrySeparation(0.0);
	leg->Draw();
	
	//superpose the x axis over the legend
	/*
	TGaxis *Nuaxis = new TGaxis(minHz,minF,maxHz, minF,minHz,maxHz,510,"");
	Nuaxis->SetTitle("log (#nu [Hz])");
	Nuaxis->SetTitleOffset(1.2);
	Nuaxis->Draw();
        */
        
        //add text on the graph
        TText *tt = new TText(15,-10.95,"XMM-Newton (2005)");
        tt->SetTextColor(14); 
        tt->SetTextSize(0.018);
        tt->Draw();
        TText *tt = new TText(18.7,-11.5,"Suzaku (2011)");
        tt->SetTextColor(14); 
        tt->SetTextSize(0.018);
        tt->Draw();	
        

	TMultiGraph *grdata = new TMultiGraph();
	gPad->SetBottomMargin(0.11);
	gPad->SetLeftMargin(0.12);
	grdata->Add(gr0);
	//grdata->Add(gr0_0);
	//grdata->Add(gr1);
	//grdata->Add(gr1_0);
	//grdata->Add(gr_Rcore);
	grdata->Add(gr_wise);
	grdata->Add(gr_atom);
	grdata->Add(gr_uvot);
        //grdata->Add(gr_Xb2,"3");
        
        //grdata->Add(gr_Xb3,"3");
        //grdata->Add(gr_Xb4,"3");
	grdata->Add(gr_XMM,"3");
	grdata->Add(gr_XMMa,"3");
        grdata->Add(gr_Suz,"3");
        grdata->Add(gr_Suza,"3");
        grdata->Add(gr_Xba,"3");
        grdata->Add(gr_Xb2a,"3");
        grdata->Add(gr_Xb3a,"3");
        grdata->Add(gr_Xb4a,"3");
        //grdata->Add(gr_xl);
        //grdata->Add(gr_xh);

	
	if (l>0){
	  grdata->Add(gr0a,"|>");
	}
	grdata->Add(gr_FERMI,"3");
	grdata->Add(gr_HESS,"3");
	if (l1>0){
	  grdata->Add(gr1_l,"|>");
	  grdata->Add(grh_l,"|>");
	}
	
	grdata->Add(gr_fermi);
	grdata->Add(gr_fermi_0);
	grdata->Add(gr_hess);
	grdata->Add(gr_hess_0);	
	//host gal
	grdata->Add(gr10,"C");

        // add energy axis
        Double_t mineV, maxeV;
        const Double_t hplanck=4.13566733E-15; // in eV s
        Double_t logh=log10(hplanck);
        mineV=minHz + logh;
        maxeV=maxHz + logh;
        
        TGaxis *eVaxis = new TGaxis(minHz,maxF,maxHz, maxF,mineV,maxeV,510,"-");
        eVaxis->SetTitle("log(E [eV])");
        eVaxis->SetTitleOffset(1.2);
        eVaxis->SetTitleFont(42);
        eVaxis->SetLabelFont(42);
        eVaxis->Draw();
        
        // add zoom
        //TPad *xpad= new TPad("xpad", "", 0.3,0.15,0.53,0.45); 
        //TPad *xpad= new TPad("xpad", "", 0.12,0.45,0.35,0.76); 
        TPad *xpad= new TPad("xpad", "", 0.37,0.15,0.59,0.45); 
        xpad->SetRightMargin(0.1);
        xpad->SetLeftMargin(0.2);
        xpad->SetTopMargin(0.05); 
        xpad->Draw(); 
        xpad->cd(); 
        
        TMultiGraph *xgrdata = grdata->Clone("xgrdata");
        xgrdata->SetMinimum(-11.6);
	xgrdata->SetMaximum(-10.5);
        xgrdata->Draw("AP");
        xgrdata->GetXaxis()->SetLimits(14.5,15.4);
        xgrdata->GetXaxis()->SetLabelSize(0.08);
        xgrdata->GetYaxis()->SetLabelSize(0.08);
        xgrdata->GetXaxis()->SetNdivisions(7);
        xgrdata->GetXaxis()->SetTickLength(0.08);
        xgrdata->GetYaxis()->SetTickLength(0.08);
        xgrdata->Draw("P");
        


	
	return 1;
}
