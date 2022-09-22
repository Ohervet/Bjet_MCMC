{
  ifstream file0, file_gal, file_3FGL, file_1FHL, file_2FHL, file_3FHL, file_HESS;
  char str[20],str1[20],str2[40],info[10];
  int i, count=0, n0=651, ns =50;
  Double_t A, B, E[n0],F[n0],deltaEmin[n0],deltaEmax[n0],deltaFmin[n0],deltaFmax[n0],Euplim[n0],Fuplim[n0],dEminuplim[n0],dEmaxuplim[n0],dFuplim[n0];
  int nR = 14;
  Double_t E_R[nR], F_R[nR],deltaEmin_R[nR],deltaEmax_R[nR],deltaFmin_R[nR],deltaFmax_R[nR];
  int nRc = 2;
  Double_t E_Rc[nRc], F_Rc[nRc],deltaEmin_Rc[nRc],deltaEmax_Rc[nRc],deltaFmin_Rc[nRc],deltaFmax_Rc[nRc];
  int nRvlbi = 2;
  Double_t E_Rvlbi[nRvlbi], F_Rvlbi[nRvlbi],deltaEmin_Rvlbi[nRvlbi],deltaEmax_Rvlbi[nRvlbi],deltaFmin_Rvlbi[nRvlbi],deltaFmax_Rvlbi[nRvlbi];
  int nW = 4;
  Double_t E_W[nW], F_W[nW],deltaEmin_W[nW],deltaEmax_W[nW],deltaFmin_W[nW],deltaFmax_W[nW];
  int n2M = 3;
  Double_t E_2M[n2M], F_2M[n2M],deltaEmin_2M[n2M],deltaEmax_2M[n2M],deltaFmin_2M[n2M],deltaFmax_2M[n2M];
  int nXS = 592;
  Double_t E_XS[nXS], F_XS[nXS],deltaEmin_XS[nXS],deltaEmax_XS[nXS],deltaFmin_XS[nXS],deltaFmax_XS[nXS];
  int nG = 2;
  Double_t E_G[nG], F_G[nG],deltaEmin_G[nG],deltaEmax_G[nG],deltaFmin_G[nG],deltaFmax_G[nG];
  //Double_t E_N[nG], F_N[nG],deltaEmin_N[nG],deltaEmax_N[nG],deltaFmin_N[nG],deltaFmax_N[nG];
  int nX = 14;
  Double_t E_X[nX], F_X[nX],deltaEmin_X[nX],deltaEmax_X[nX],deltaFmin_X[nX],deltaFmax_X[nX];
  int n3F = 5;
  Double_t E_3F[n3F], F_3F[n3F],deltaEmin_3F[n3F],deltaEmax_3F[n3F],deltaFmin_3F[n3F],deltaFmax_3F[n3F];
  int nF = 5;
  Double_t E_F[nF], F_F[nF],deltaEmin_F[nF],deltaEmax_F[nF],deltaFmin_F[nF],deltaFmax_F[nF];
  int nHC = 6;
  Double_t E_HC[nHC], F_HC[nHC],deltaEmin_HC[nHC],deltaEmax_HC[nHC],deltaFmin_HC[nHC],deltaFmax_HC[nHC];
  
  Double_t Freq_3Fs[ns], Flux_3Fs[ns], DFmin_3Fs[ns], DFmax_3Fs[ns];
  Double_t Freq_1FHs[ns], Flux_1FHs[ns], DFmin_1FHs[ns], DFmax_1FHs[ns];
  Double_t Freq_2FHs[ns], Flux_2FHs[ns], DFmin_2FHs[ns], DFmax_2FHs[ns];
  Double_t Freq_3FHs[ns], Flux_3FHs[ns], DFmin_3FHs[ns], DFmax_3FHs[ns];
  Double_t Freq_Hs[ns], Flux_Hs[ns], DFmin_Hs[ns], DFmax_Hs[ns];
  
  double    Freq, FreqLow, FreqHigh, Flux, FluxLow, FluxHigh, DFlux, DFluxLow, DFluxHigh;
  
  const double planck = 4.135667517e-15; // h in [eV s]

  file0.open("data/PMNJ_1603-4904/PMN_J1603_SED_3.dat");
  file_3FGL.open("data/PMNJ_1603-4904/PMNJ_1603_3FGL_spectrum.dat");
  file_1FHL.open("data/PMNJ_1603-4904/PMNJ_1603_1FHL_spectrum.dat");
  file_2FHL.open("data/PMNJ_1603-4904/PMNJ_1603_2FHL_spectrum.dat");
  file_3FHL.open("data/PMNJ_1603-4904/PMNJ_1603_3FHL_spectrum.dat");
  //file_HESS.open("data/PMNJ_1603-4904/PMNJ_1603_HESSMono_spectrum.dat");
  file_HESS.open("data/PMNJ_1603_HESS1Loose_UL.dat");
   

    
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
 
     int j=0, k=0, l=0, l1=0;
     
     for(int i=0;i<=8;i++){
        file0>>str;
     }

      std::cout<<"i\tFreq\tFlux\tFreqLow\tFreqHight\tFluxLow\tFluxHigh\tData"<<std::endl;
      
        
      //sed data
      for(int i=0;i<n0;i++)
      {file0>>A>>B>>deltaEmin[i]>>deltaEmax[i]>>deltaFmin[i]>>deltaFmax[i]>>str2;

            A = A / planck; 
            E[i]=log10(A);
 	    F[i]=log10(B);
	    deltaEmin[i]=E[i]-log10(deltaEmin[i]/ planck);
            //deltaEmin[i]= log10(deltaEmin[i]/ planck);
            //std::cout<<deltaEmin[i]<<std::endl;
	    deltaEmax[i]=-E[i]+log10(deltaEmax[i]/ planck);
	
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
  
        
        
	// Host galaxy
        if(nb_g>0){
          gr10 = new TGraph(nb_g,freq_g,fflux_g);
	  gr10->SetLineColor(2);
          gr10->SetLineStyle(1);
	  gr10->Draw("C");
        }
        
    //Radio data kpc
    for(i=0;i<nR;i++){
        E_R[i] = E[count+i];
        F_R[i] = F[count+i];
        deltaEmin_R[i] = deltaEmin[count+i];
        deltaEmax_R[i] = deltaEmax[count+i];
        deltaFmin_R[i] = deltaFmin[count+i];
        deltaFmax_R[i] = deltaFmax[count+i];
    }   
    
    //Radio core vlbi
    count += nR;
    for(i=0;i<nRc;i++){
        E_Rc[i] = E[count+i];
        F_Rc[i] = F[count+i];
        deltaEmin_Rc[i] = deltaEmin[count+i];
        deltaEmax_Rc[i] = deltaEmax[count+i];
        deltaFmin_Rc[i] = deltaFmin[count+i];
        deltaFmax_Rc[i] = deltaFmax[count+i];
    }   
    
    //Radio total vlbi
    count += nRc;
    for(i=0;i<nRvlbi;i++){
        E_Rvlbi[i] = E[count+i];
        F_Rvlbi[i] = F[count+i];
        deltaEmin_Rvlbi[i] = deltaEmin[count+i];
        deltaEmax_Rvlbi[i] = deltaEmax[count+i];
        deltaFmin_Rvlbi[i] = deltaFmin[count+i];
        deltaFmax_Rvlbi[i] = deltaFmax[count+i];
    }   

    //Wise data
    count += nRvlbi;
    for(i=0;i<nW;i++){
        E_W[i] = E[count+i];
        F_W[i] = F[count+i];
        deltaEmin_W[i] = deltaEmin[count+i];
        deltaEmax_W[i] = deltaEmax[count+i];
        deltaFmin_W[i] = deltaFmin[count+i];
        deltaFmax_W[i] = deltaFmax[count+i];
    }  
    
    //2MASS data
    count += nW;
    for(i=0;i<n2M;i++){
        E_2M[i] = E[count+i];
        F_2M[i] = F[count+i];
        deltaEmin_2M[i] = deltaEmin[count+i];
        deltaEmax_2M[i] = deltaEmax[count+i];
        deltaFmin_2M[i] = deltaFmin[count+i];
        deltaFmax_2M[i] = deltaFmax[count+i];
    }  

    //X-Shooter data
    count += n2M;
    for(i=0;i<nXS;i++){
        E_XS[i] = E[count+i];
        F_XS[i] = F[count+i];
        deltaEmin_XS[i] = deltaEmin[count+i];
        deltaEmax_XS[i] = deltaEmax[count+i];
        deltaFmin_XS[i] = deltaFmin[count+i];
        deltaFmax_XS[i] = deltaFmax[count+i];
    }  
    
    //GMOS data
    count += nXS;
    for(i=0;i<nG;i++){
        E_G[i] = E[count+i];
        F_G[i] = F[count+i];
        deltaEmin_G[i] = deltaEmin[count+i];
        deltaEmax_G[i] = deltaEmax[count+i];
        deltaFmin_G[i] = deltaFmin[count+i];
        deltaFmax_G[i] = deltaFmax[count+i];
    }  
/*
    //NTT data
    count += nG;
    for(i=0;i<nG;i++){
        E_N[i] = E[count+i];
        F_N[i] = F[count+i];
        deltaEmin_N[i] = deltaEmin[count+i];
        deltaEmax_N[i] = deltaEmax[count+i];
        deltaFmin_N[i] = deltaFmin[count+i];
        deltaFmax_N[i] = deltaFmax[count+i];
    }  
*/
    //X-ray data
    count += nG;
    for(i=0;i<nX;i++){
        E_X[i] = E[count+i];
        F_X[i] = F[count+i];
        deltaEmin_X[i] = deltaEmin[count+i];
        deltaEmax_X[i] = deltaEmax[count+i];
        deltaFmin_X[i] = deltaFmin[count+i];
        deltaFmax_X[i] = deltaFmax[count+i];
    }  

    //Fermi 3FGL data
    count += nX;
    for(i=0;i<n3F;i++){
        E_3F[i] = E[count+i];
        F_3F[i] = F[count+i];
        deltaEmin_3F[i] = deltaEmin[count+i];
        deltaEmax_3F[i] = deltaEmax[count+i];
        deltaFmin_3F[i] = deltaFmin[count+i];
        deltaFmax_3F[i] = deltaFmax[count+i];
    }  
    //Fermi Michael data
    count += n3F;
    for(i=0;i<nF;i++){
        E_F[i] = E[count+i];
        F_F[i] = F[count+i];
        deltaEmin_F[i] = deltaEmin[count+i];
        deltaEmax_F[i] = deltaEmax[count+i];
        deltaFmin_F[i] = deltaFmin[count+i];
        deltaFmax_F[i] = deltaFmax[count+i];
    }  

      //Fermi 1FHL spectrum
      for(int i=0;i<ns;i++)
      {
	file_1FHL>>Freq>>Flux>>DFluxLow>>DFluxHigh;

            Freq_1FHs[i]=log10(Freq);
 	    Flux_1FHs[i]=log10(Flux);
	    DFmin_1FHs[i] = Flux_1FHs[i] - log10((Flux-DFluxLow));
	    DFmax_1FHs[i] = log10((Flux+DFluxHigh)) - Flux_1FHs[i];
            
            if (!file_1FHL.good()) break;
        }
        
      //Fermi 2FHL spectrum
      for(int i=0;i<ns;i++)
      {
	file_2FHL>>Freq>>Flux>>DFluxLow>>DFluxHigh;

            Freq_2FHs[i]=log10(Freq);
 	    Flux_2FHs[i]=log10(Flux);
	    DFmin_2FHs[i] = Flux_2FHs[i] - log10((Flux-DFluxLow));
	    DFmax_2FHs[i] = log10((Flux+DFluxHigh)) - Flux_2FHs[i];
            
            if (!file_2FHL.good()) break;
        }
        
      //Fermi 3FHL spectrum
      for(int i=0;i<ns;i++)
      {
	file_3FHL>>Freq>>Flux>>DFluxLow>>DFluxHigh;

            Freq_3FHs[i]=log10(Freq);
 	    Flux_3FHs[i]=log10(Flux);
	    DFmin_3FHs[i] = Flux_3FHs[i] - log10((Flux-DFluxLow));
	    DFmax_3FHs[i] = log10((Flux+DFluxHigh)) - Flux_3FHs[i];
            
            if (!file_3FHL.good()) break;
        }
        
      //Fermi 3FGL spectrum
      for(int i=0;i<ns;i++)
      {
	file_3FGL>>Freq>>Flux>>DFluxLow>>DFluxHigh;

            Freq_3Fs[i]=log10(Freq);
 	    Flux_3Fs[i]=log10(Flux);
	    DFmin_3Fs[i] = Flux_3Fs[i] - log10((Flux-DFluxLow));
	    DFmax_3Fs[i] = log10((Flux+DFluxHigh)) - Flux_3Fs[i];
            
            if (!file_3FGL.good()) break;
        }
        
      //HESS spectrum
      for(int i=0;i<ns;i++)
      {
	file_HESS>>Freq>>Flux>>DFluxLow>>DFluxHigh;

            Freq_Hs[i]=log10(Freq);
 	    Flux_Hs[i]=log10(Flux);
	    DFmin_Hs[i] = Flux_Hs[i] - log10((Flux-DFluxLow));
	    DFmax_Hs[i] = log10((Flux+DFluxHigh)) - Flux_Hs[i];
            
            if (!file_HESS.good()) break;
        }
        
    //HESS combined data
    count += nF;
    for(i=0;i<nHC;i++){
        E_HC[i] = E[count+i];
        F_HC[i] = F[count+i];
        deltaEmin_HC[i] = deltaEmin[count+i];
        deltaEmax_HC[i] = deltaEmax[count+i];
        deltaFmin_HC[i] = deltaFmin[count+i];
        deltaFmax_HC[i] = deltaFmax[count+i];
    }  
    

//----------Display----------//	
    gPad->SetTopMargin(0.1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetTicky();

    Double_t minHz, maxHz,minF,maxF;
    /*
    minHz=7.5;
    maxHz=27.5;    
    minF=-15.5;
    maxF=-9.0;
*/
    //Zoom

    minHz=22;
    maxHz=26.5; 
    minF=-12;
    maxF=-10;

    //nb1=k;
    //display errorbars
    TGraphAsymmErrors *gr0 = new TGraphAsymmErrors(n0,E,F,deltaEmin,deltaEmax,deltaFmin,deltaFmax);
    gr0->SetMarkerColor(0);
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
    gr0->GetXaxis()->SetRangeUser(minHz,maxHz);
    gr0->SetMinimum(minF);
    gr0->SetMaximum(maxF);
    
 
    //Radio data 
    gr_R = new TGraphAsymmErrors(nR,E_R,F_R,deltaEmin_R,deltaEmax_R,deltaFmin_R,deltaFmax_R);
    gr_R->SetMarkerStyle(20);
    gr_R->SetMarkerSize(0.8);/*
    gr_R->SetMarkerColor(kRed+2);
    gr_R->SetLineColor(kRed+2);
    gr_R->SetFillColor(kRed+2);*/
    gr_R->Draw("P");
    
    gr_Rc = new TGraphAsymmErrors(nRc,E_Rc,F_Rc,deltaEmin_Rc,deltaEmax_Rc,deltaFmin_Rc,deltaFmax_Rc);
    gr_Rc->SetMarkerStyle(20);
    gr_Rc->SetMarkerSize(0.8);
    gr_Rc->SetMarkerColor(kRed);
    gr_Rc->SetLineColor(kRed);
    gr_Rc->SetFillColor(kRed);
    gr_Rc->Draw("P");
    
    gr_Rvlbi = new TGraphAsymmErrors(nRvlbi,E_Rvlbi,F_Rvlbi,deltaEmin_Rvlbi,deltaEmax_Rvlbi,deltaFmin_Rvlbi,deltaFmax_Rvlbi);
    gr_Rvlbi->SetMarkerStyle(20);
    gr_Rvlbi->SetMarkerSize(0.8);
    gr_Rvlbi->SetMarkerColor(kRed+2);
    gr_Rvlbi->SetLineColor(kRed+2);
    gr_Rvlbi->SetFillColor(kRed+2);
    gr_Rvlbi->Draw("P");

    //Wise data 
    gr_W = new TGraphAsymmErrors(nW,E_W,F_W,deltaEmin_W,deltaEmax_W,deltaFmin_W,deltaFmax_W);
    gr_W->SetMarkerStyle(20);
    gr_W->SetMarkerSize(0.8);
    gr_W->SetMarkerColor(kOrange+2);
    gr_W->SetLineColor(kOrange+2);
    gr_W->SetFillColor(kOrange+2);
    gr_W->Draw("P");
    
    //X-Shooter data
    gr_XS = new TGraphAsymmErrors(nXS,E_XS,F_XS,deltaEmin_XS,deltaEmax_XS,deltaFmin_XS,deltaFmax_XS);
    gr_XS->SetMarkerStyle(20);
    gr_XS->SetMarkerSize(0.5);
    gr_XS->SetMarkerColor(kSpring+2);
    gr_XS->SetLineColor(kSpring+2);
    gr_XS->SetFillColor(kSpring+2);
    gr_XS->Draw("P");
    
    //2MASS data
    gr_2M = new TGraphAsymmErrors(n2M,E_2M,F_2M,deltaEmin_2M,deltaEmax_2M,deltaFmin_2M,deltaFmax_2M);
    gr_2M->SetMarkerStyle(20);
    gr_2M->SetMarkerSize(0.8);
    gr_2M->SetMarkerColor(kYellow+2);
    gr_2M->SetLineColor(kYellow+2);
    gr_2M->SetFillColor(kYellow+2);
    gr_2M->Draw("P");

    //GMOS data
    gr_G = new TGraphAsymmErrors(nG,E_G,F_G,deltaEmin_G,deltaEmax_G,deltaFmin_G,deltaFmax_G);
    gr_G->SetMarkerStyle(20);
    gr_G->SetMarkerSize(0.8);
    gr_G->SetMarkerColor(kTeal+2);
    gr_G->SetLineColor(kTeal+2);
    gr_G->SetFillColor(kTeal+2);
    gr_G->Draw("P");
/*
    //NTT data
    gr_N = new TGraphAsymmErrors(nG,E_N,F_N,deltaEmin_N,deltaEmax_N,deltaFmin_N,deltaFmax_N);
    gr_N->SetMarkerStyle(20);
    gr_N->SetMarkerSize(0.8);
    gr_N->SetMarkerColor(kTeal+2);
    gr_N->SetLineColor(kTeal+2);
    gr_N->SetFillColor(kTeal+2);
    gr_N->Draw("P");
*/
    //X-ray data
    gr_X = new TGraphAsymmErrors(nX,E_X,F_X,deltaEmin_X,deltaEmax_X,deltaFmin_X,deltaFmax_X);
    gr_X->SetMarkerStyle(20);
    gr_X->SetMarkerSize(0.8);
    gr_X->SetMarkerColor(kCyan+2);
    gr_X->SetLineColor(kCyan+2);
    gr_X->SetFillColor(kCyan+2);
    gr_X->Draw("P");
    
    //Fermi 1FHL spectrum 
    /*
    gr_1FHLs = new TGraphAsymmErrors(ns,Freq_1FHs,Flux_1FHs,0,0,DFmin_1FHs,DFmax_1FHs);
    gr_1FHLs->SetFillColor(kBlue-10);
    gr_1FHLs->Draw("3");
    */

    //Fermi 2FHL spectrum 
    gr_2FHLs = new TGraphAsymmErrors(ns,Freq_2FHs,Flux_2FHs,0,0,DFmin_2FHs,DFmax_2FHs);
    gr_2FHLs->SetFillColor(kBlue-10);
    //gr_2FHLs->Draw("3");
    
    //Fermi 3FHL spectrum 
    gr_3FHLs = new TGraphAsymmErrors(ns,Freq_3FHs,Flux_3FHs,0,0,DFmin_3FHs,DFmax_3FHs);
    gr_3FHLs->SetLineColor(kViolet+2);
    gr_3FHLs->SetLineWidth(2);
    gr_3FHLs->SetFillColor(kViolet+2);
    gr_3FHLs->SetFillStyle(0);
    //gr_3FHLs->Draw("3");

    //Fermi 3FGL spectrum 
    gr_3FGLs = new TGraphAsymmErrors(ns,Freq_3Fs,Flux_3Fs,0,0,DFmin_3Fs,DFmax_3Fs);
    gr_3FGLs->SetFillColor(18);
    gr_3FGLs->Draw("3");
    
    //Fermi 3FGL data
    gr_3F = new TGraphAsymmErrors(n3F,E_3F,F_3F,deltaEmin_3F,deltaEmax_3F,deltaFmin_3F,deltaFmax_3F);
    gr_3F->SetMarkerStyle(20);
    gr_3F->SetMarkerSize(0.8);
    gr_3F->SetMarkerColor(16);
    gr_3F->SetLineColor(16);
    gr_3F->SetFillColor(16);
    gr_3F->Draw("P");    

    //Fermi Felicia data
    gr_F = new TGraphAsymmErrors(nF,E_F,F_F,deltaEmin_F,deltaEmax_F,deltaFmin_F,deltaFmax_F);
    gr_F->SetMarkerStyle(20);
    gr_F->SetMarkerSize(0.8);
    gr_F->SetMarkerColor(kRed+2);
    gr_F->SetLineColor(kRed+2);
    gr_F->SetFillColor(kRed+2);
    gr_F->Draw("P");   
    
    //HESS Mono spectrum 
    gr_Hs = new TGraphAsymmErrors(ns,Freq_Hs,Flux_Hs,0,0,DFmin_Hs,DFmax_Hs);
    gr_Hs->SetFillColor(kRed-7);
    //gr_Hs->Draw("3");
    
    
    //Upper limits
    if (l>0){

        TGraphAsymmErrors *gr0u = new TGraphAsymmErrors(l,Euplim,Fuplim,0,0,dFuplim,0);
        gr0u->SetLineColor(kMagenta+2);
        gr0u->SetFillColor(kMagenta+2);
        //gr0u->Draw("|>");
        
        TGraphAsymmErrors *gr0uF = new TGraphAsymmErrors(2,Euplim,Fuplim,0,0,dFuplim,0);
        gr0uF->SetLineColor(kRed+2);
        gr0uF->SetFillColor(kRed+2);
        gr0uF->Draw("|>");
    }

    //HESS combined data
    gr_HC = new TGraphAsymmErrors(nHC,E_HC,F_HC,deltaEmin_HC,deltaEmax_HC,deltaFmin_HC,deltaFmax_HC);
    gr_HC->SetMarkerStyle(20);
    gr_HC->SetMarkerSize(0.8);
    gr_HC->SetMarkerColor(kMagenta+2);
    gr_HC->SetLineColor(kMagenta+2);
    //gr_HC->Draw("P");
    
    TMultiGraph *grdata = new TMultiGraph();
    grdata->Add(gr10,"C");
        
//----------Legende----------//        
        
	leg = new TLegend(0.12,0.78,0.95,0.87);
	leg->SetTextSize(0.03);
	leg->SetFillColor(kWhite);
	leg->SetMargin(0.2);
	leg->SetNColumns(7);/*
        leg->AddEntry(gr_R,"Radio kpc","lep");
        leg->AddEntry(gr_Rc,"Radio core","lep");
        leg->AddEntry(gr_2M,"2MASS","lep");
        leg->AddEntry(gr_G,"GMOS","lep");*/
        leg->AddEntry(gr_3FGLs,"3FGL","f");/*
        leg->AddEntry(gr_3FHLs,"3FHL","f");
        leg->AddEntry(gr_HC,"H.E.S.S.","lep");
        leg->AddEntry(gr_Rvlbi,"Radio vlbi","lep");
        leg->AddEntry(gr_W,"Wise","lep");
        leg->AddEntry(gr_XS,"X-Shooter","lep");
        leg->AddEntry(gr_X,"Suzaku,XMM","lep");
        //leg->AddEntry(gr_1FHLs,"1FHL","f");
        leg->AddEntry(gr_2FHLs,"2FHL","f");*/
        leg->AddEntry(gr_F,"Our analysis Fermi-LAT","lep");
        
        
        /*
        leg->SetNColumns(3);
        leg->AddEntry(gr_2M,"2MASS","lep");
        leg->AddEntry(gr_XS,"X-Shooter","lep");
        leg->AddEntry(gr_G,"GMOS","lep");*/
	leg->SetEntrySeparation(0.0);
	leg->Draw();
        
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
        

//----------Add data to the model----------//
        
	TMultiGraph *grdata = new TMultiGraph();
	grdata->SetTitle("SED Data");
        grdata->Add(gr_R);
        grdata->Add(gr_W);
        grdata->Add(gr_XS);
        grdata->Add(gr_2M);
        grdata->Add(gr_G);
        grdata->Add(gr_X);
        grdata->Add(gr_2FHLs,"3");
        grdata->Add(gr_3FHLs,"3");
        grdata->Add(gr_3FGLs,"3");
        grdata->Add(gr_3F);
        grdata->Add(gr_F);
	if (l>0){
          grdata->Add(gr_HC);
	  grdata->Add(gr0u,"|>");
          grdata->Add(gr0uF,"|>");
	}
        
        
        
        
        
        
        
        /*
	grdata->Add(gr_FERMI,"3");
        grdata->Add(gr_VERITAS,"3");
	grdata->Add(gr0);
	grdata->Add(gr0_0);
	if (l>0){
	  grdata->Add(gr0a,"|>");
	  grdata->Add(gr_F_l,"|>");
          grdata->Add(gr_V_l,"|>");
	}
	grdata->Add(gr_R);
	grdata->Add(gr_F_0);
	grdata->Add(gr_V);
	grdata->Add(gr_V_0);
	//host gal
	grdata->Add(gr10,"C");
*/
	
	return 1;
}
