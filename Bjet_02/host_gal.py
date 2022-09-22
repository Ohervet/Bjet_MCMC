#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

from pylab import*
import csv
from scipy.interpolate import UnivariateSpline

#-----Flags-----#
Abraham = 0
Rc = 0
R_Johnson = 1
Ks_SOFI = 0

Aplib = 0
PKS0625 = 0
Bllac = 0
Wcomae = 0
HESSJ1943 = 0
NGC1275 = 1

UVOT = 1
ATOM = 0
WISE = 0
TMASS = 0

c = 2.998e8

#----------Fonctions usuelles----------#

def convfreq(lambd):
    #angström
    return 2.99792458e8/(lambd*1e-10)

def moy(liste):
    somme=0.
    for i in range (len(liste)):
        somme= somme+liste[i]
    if len(liste) != 0:
        return somme/len(liste)
    else:
        return 0

def factorielle(arg):
  if arg==0:
    return 1
  return arg*factorielle(arg-1)

def YOUNG(a,ae):
    #luminosité relative intégrée dans un cercle (YOUNG 1976)
    #Loi de Vaucouleur en r 1/4
    a = float(a)
    ae = float(ae)
    b = 7.66924944
    S = 0
    T = a/ae
    for n in range (0,8):
        S += b**n * T**(n/4.)/(factorielle(n))
        #print n
    return 1-exp(-b*T**(1./4.))*S
    
def YOUNG_Sersic(a,ae,N):
    #luminosité relative intégrée dans un cercle (YOUNG 1976)
    #Profil de Sersic, generalisation de Vaucouleur en r 1/N
    #égal a de Vaucouleur pour N = 4
    a = float(a)
    ae = float(ae)
    b = 7.66924944
    S = 0
    T = a/ae
    for i in range (0,8):
        S += b**i * T**(i/N)/(factorielle(i))
    return 1-exp(-b*T**(1./N))*S
    
def ConvMagnFlux_Ks2MASS(magn):
    #conv en erg.cm-2.s-2
    return 1.387e14*10**((magn +1.85 + 48.6)/(-2.5))


#----------lecture des données PEGASE----------#

fichier = csv.reader(open ("data/HGS_13.dat", "r"), delimiter=" ")
#fichier = csv.reader(open ("data/HGS_13.dat","r"), delimiter=" ")

matrix = []
flux_pond= []
flux_pond_uvot= []
flux_moy= []
flux_pond1= []
flux_moy1= []
flux_pond2= []
flux_moy2= []
	
for line in fichier:
        matrix.append(line)




#-----------------------------#
#          Ap Lib
#-----------------------------#
if Aplib == 1:

    #rayon effectif d'Ap Lib en U ["]
    reffU = 2.55
    
    #rayon effectif d'Ap Lib en B ["]
    reffB = 2.9
    
    #rayon effectif d'Ap Lib en V ["]
    reffV = 5.3 # Valeur assez improbable (Visvanatanh 1977)
    
    #rayon effectif d'Ap Lib en R ["]
    reffRkp = 5.66 #old Abraham 1991 (R kitt peak)
    reffRc = 3.7 #HST scarpa 2000
    reffR = 6.72 #Pursimo 2002
    
    reff = [reffRc, reffB, reffU]
    
    #Flux UVOT non corrigé:
    uvot = [2.6093E-11, 1.9526E-11, 1.3217E-11, 1.0630E-11, 1.1317E-11, 1.1070E-11]
    
    #Flux wise non corrigé
    wise = [2.12e-11, 2.13e-11, 2.51e-11, 2.76e-11] #Vizier
    
    #Flux 2MASS non corrigé
    mass = [2.708e-11, 2.652e-11, 2.372e-11] #5 arcsec Vizier


#-----------------------------#
#          PKS 0625-354
#-----------------------------#

if PKS0625 == 1:
    print "PKS 0625-354"

    #rayon effectif en Rc ["] (Govoni et al. 2000)
    #Reff = 13.08 #avec table 3
        
    #rayon effectif Ks ["] (Inskip et al. 2010)
    Reff = 18.77
    
    #flux de la galaxie hôte en Rc [erg.cm-2.s-1], déduit de (Govoni et al. 2000)
    #F_Rc = 7.366e-11 #table 3
    #F_Rc = 4.332e-11 #table 2
    #F_Rc = 6.354e-11 #table 3 Mhost tot en corrigeant l'absorption par une diff de 0.2 magn (entre Schlafly 2011 et leur article 2000)
    
    #flux de la galaxie hôte en Ks [erg.cm-2.s-1]
    #F_Ks = 5.428e-11   (Inskip et al. 2010) 
    
    #Flux UVOT non corrigés V,B,U,UVW1,UVM2,UVW2 :
    uvot = [1.404722e-11, 8.13053906379e-12, 5.05785149231e-12, 3.99987733815e-12, 4.95892865495e-12, 4.75420812174e-12]
    
    
    #flux high uvot
#    uvot = [1.88751485e-11 + 9.67956332e-13, 1.14083968e-11 + 5.45856305e-13, 6.78772635e-12 + 4.53522919e-13,
#            4.39498354e-12 + 4.04559155e-13, 4.82402005e-12+ 2.96862772e-13, 4.08855756e-12 + 2.66645058e-13]
    
    #Flux ATOM non corrigé R:
    atom = 2.00021599e-11
    

#-----------------------------#
#          Bl Lac
#-----------------------------#

if Bllac == 1:

    #rayon effectif en Rc ["] (HST, Scarpa et al. 2000)
    Reff = 4.8
    
    #densité de flux de la galaxie hôte en Rc [erg.cm-2.s-1], déduit de (Scarpa et al. 2000)
    F_Rc = 9.09e-12
    
    #Flux UVOT non corrigés V,B,U,UVW1,UVM2,UVW2 (from ASDC) :
    uvot = [3.4277E-11, 2.5942E-11, 1.9953E-11, 1.2823E-11, 1.5668E-11, 1.2388E-11]
    

#-----------------------------#
#          W Comae
#-----------------------------#

if Wcomae == 1:

    #rayon effectif en Rc ["] (HST, Scarpa et al. 2000)
    Reff = 2.1 
    DReff = 0.4
    
    Reff = Reff - DReff
    
    #densité de flux de la galaxie hôte en Rc [erg.cm-2.s-1], déduit de (Nilsson et al. 2003)
    F_Rc = 3.319e-12
    DF_Rc = 0.305e-12
    
    F_Rc =F_Rc + DF_Rc
    
    #Flux UVOT non corrigés V,B,U,UVW1,UVM2,UVW2 (from ASDC) :  (UVM2 non dispo)
    uvot = [2.06063E-011, 2.07014E-011, 1.93197E-011, 1.77828E-011, 0.0, 1.65577E-011]

 
#-----------------------------#
#          HESS J1943
#-----------------------------#

if HESSJ1943 == 1:

    #rayon effectif max en 2MASS K["] (HST, Scarpa et al. 2000)
    Reff = 2.5 
    N = 8.0
    ouv = 4.0
    
    print "Fraction de luminosite de la galaxy hôte pour une ouverture de 5 arcsec:", YOUNG_Sersic(ouv,Reff, N) , "\n"
    

#-----------------------------#
#          NGC 1275
#-----------------------------#

if NGC1275 == 1:

    #rayon effectif de Spitzer (3.6um) [kpc] (Sani et al. 2018)
    Reff_kpc = 42 #+-15
    
    #Mathews 2006 state an effective radius of 6.41 kpc
    Reff_kpc = 6.41
    
    #using the cosmology of Sani 2018 (H0 =70), we have 0.357kpc/"
    Reff = Reff_kpc/0.357 #["]    
    
    #most reliable source "Third reference catalogue of bright galaxies, De Vaucouleur 1991"
    Reff = 16.9#+2.5 -2.2 ["]
    
    N = 4.0 # mean it follows a de Vaucouleur profile (as stated by Sani 2018)
    ouv = 5.0
    print "Fraction de luminosite de la galaxy hôte pour une ouverture de 5 arcsec:", YOUNG_Sersic(ouv,Reff, N) , "\n"
    
    #Flux UVOT non corrigés V,B,U,UVW1,UVM2,UVW2 [erg cm-2 s-1]:
    #from NGC1275_Flare2017_SED.dat
    F_uvot = [1.2233E-10, 1.0898E-10, 7.3274E-11, 5.8938E-11, 6.4023E-11, 5.1152E-11]
    print('Flux UVOT non corriges:',F_uvot)


#----------calcul du flux de la galaxy hote dans les filtres d'UVOT----------#

if UVOT == 1:

    #ouverture UVOT ["]
    ouv0 = 5.0
    #PouvB = ouv0/reffB
    #PouvU = ouv0/reffU
    
    
    freq_min= []
    freq_max= []
    #ponderationV= YOUNG(ouv0,reffV)
    #ponderationB= YOUNG(ouv0,reffB)    #Luminosite dans le champ de UVOT/ luminosite totale de la galaxie (Young)
    #ponderationU= YOUNG(ouv0,reffU)
    
    #----------Filtes de UVOT----------#
    
    #longueurs d'onde centrales des filtres V,B,U,UVW1,UVM2,UVW2 [A]
    lambd_centr=[5468, 4392, 3465, 2600, 2246, 1928]
    #nu_centr = c/(array(lambd_centr)*1e-10)
    nu_centr = convfreq(array(lambd_centr))
    ponderation = zeros(len(lambd_centr))
    for i in range(len(lambd_centr)):
    #    if i == 0:
    #        reffmin_UVOT = reffB
    #        reffmax_UVOT = reffR
    #        ponderationVmin= YOUNG(ouv0,reffmax_UVOT)
    #        ponderationVmax= YOUNG(ouv0,reffmin_UVOT)
    #    if i == 1:
    #        reff_UVOT.append(reffB)
    #    elif i >= 2:
    #        reff_UVOT.append(reffU)
    #    else:
    #        reff_UVOT.append(s(nu_centr[i]))
        ponderation[i]= YOUNG(ouv0,Reff)
    print ponderation
    
    
    
    #largeur a mi-hauteur des filtres V,B,U,UVW1,UVM2,UVW2 [A]
    D_lambd= [769, 975, 785, 693, 498, 657]
    
    
    for i in range (len(lambd_centr)):
        freq_min.append(convfreq(lambd_centr[i] + D_lambd[i]/2.))
        freq_max.append(convfreq(lambd_centr[i] - D_lambd[i]/2.))
        
    #lecture des donnees PEGASE
    for i in range(len(freq_min)):
        freq=[]
        flux=[]
        for j in range (len(matrix)):
            if (matrix [j][0][0])!= "#" and (float(matrix [j][0]) >= freq_min[i]) and (float(matrix [j][0]) <= freq_max[i]) :
                freq.append(float(matrix [j][0]))
                flux.append(float(matrix [j][2]))
    #    if i == 0:      
    #        flux_Vmin = moy(flux)*ponderationVmin
    #        flux_Vmax = moy(flux)*ponderationVmax
    #        flux_pond.append(moy(flux)*ponderation[i])
    #    elif i == 1:
    #        #flux_pond.append(moy(flux)*ponderationB)
    #        flux_pond.append(moy(flux)*ponderation[i])
    #    elif i > 1:
    #        flux_pond.append(moy(flux)*ponderationU)
    #        #flux_pond.append(moy(flux)*ponderation[i])
        flux_pond_uvot.append(moy(flux)*ponderation[i])
        flux_moy.append(moy(flux))
        
    flux_moy_uvot = flux_moy
        
   
    #print "flux de la galaxie hote dans les filtres V, valeurs min et max: \n",flux_Vmin, flux_Vmax
    print "flux moyen de la galaxie hote sur les bandes UVOT sans prendre en compte l'ouverture\n", flux_moy,"\n"
    print "flux de la galaxie hote dans les filtres V,B,U,UVW1,UVM2,UVW2 de UVOT [erg.cm-2.s-1]: \n",flux_pond_uvot ,"\n"
    print "flux corrige UVOT \n", array(F_uvot)- array(flux_pond_uvot), "\n"
    #print "flux corrigé UVOT en V, valeurs min et max \n", uvot[0]- flux_Vmax, uvot[0]- flux_Vmin,"\n"
    print "percentage Host contamination UVOT \n", array(flux_pond_uvot) /(array(F_uvot)- array(flux_pond_uvot))*100
    
    print "B-V UVOT=",-2.5*log10(flux_moy[1]/flux_moy[0])
    

#----------calcul du flux de la galaxy hote dans les filtres de WISE----------#

if WISE == 1:

    #ouverture wise ["]
    ouv = array([12.0, 6.5, 6.4, 6.1])
    #Pouv = ouv/reffR
    #ponderation1 = [8.51e-1, 7.14e-1, 7.14e-1, 7.03e-1] #tables de young #Luminosite dans ls champ de WISE/ luminosite totale de la galaxie (reffB)
    #ponderation1 = [7.14e-1, 5.51e-1, 5.27e-1, 5.27e-1] #tables de young #Luminosite dans ls champ de WISE/ luminosite totale de la galaxie (reffR)
    #ponderation1 = [6.215e-1, 4.55e-1, 4.4e-1, 4.3e-1] #extrapolation linéaire de reff (voir plot d'en bas pour la freq max de 2Mass) Reff = 7.73"
    
    
    freq_min1= []
    freq_max1= []
    freq_centr =[]
       
    #----------Filtes de WISE----------#
    
    #longueurs d'onde centrales des filtres de WISE [A]
    lambd_centr1=[2.2e5, 1.2e5, 4.6e4, 3.4e4]
    nu_centr1 = c/(array(lambd_centr1)*1e-10)
    reff_W= []
    ponderation1 = zeros(len(ouv))
    for i in range(len(lambd_centr1)):
        #reff_W.append(s(nu_centr1[i])) #extrapolation des rayons effectifs
        reff_W.append(reffR) #limite inférieure du rayon effectif
        ponderation1[i]= YOUNG(ouv[i],reff_W[i]) #fraction de luminosité ds le champ de wise (limite sup)
    
    #largeur a mi-hauteur des filtres W4, W3, W2, W1 [A]
    D_lambd1= [3.5e4, 8.8e4, 3.1e4, 3.5e4]
    
    
    for i in range (len(lambd_centr1)):
        freq_min1.append(convfreq(lambd_centr1[i] + D_lambd1[i]/2.))
        freq_max1.append(convfreq(lambd_centr1[i] - D_lambd1[i]/2.))
        freq_centr.append(convfreq(lambd_centr1[i]))
    
    #lecture des donnees PEGASE
    for i in range(len(freq_min1)):
        freq1=[]
        flux1=[]
        for j in range (len(matrix)):          
            if (matrix [j][0][0])!= "#" and (float(matrix [j][0]) >= freq_min1[i]) and (float(matrix [j][0]) <= freq_max1[i]) :
                freq1.append(float(matrix [j][0]))
                flux1.append(float(matrix [j][2]))
        
        flux_pond1.append(moy(flux1)*ponderation1[i])
        flux_moy1.append(moy(flux1))
    
    print "flux de la galaxie hote dans les filtres W4, W3, W2, W1 de WISE [erg.cm-2.s-1]: \n",flux_pond1
    print "flux moyen \n", flux_moy1
    print "flux corrigé WISE \n", array(wise)- array(flux_pond1), "\n"
    #print "flux  WISE intermediare (entre valeur sans et avec corr)\n", array(wise)- array(flux_pond1)/2, "\n"

   

#----------calcul du flux de la galaxy hote dans les filtres de 2MASS----------#

if TMASS == 1:
    #ouverture 2MASS ["]
    ouv1 = 5.0
    Pouv1 = ouv1/reffRc
    #ponderation2 = 8.67e-1 #reffB 12.9sec
    #ponderation2 = 7.26e-1 #reffR 12.9sec
    #ponderation2 = 4.71e-1 #reffR 5 sec
    ponderation2 = 3.84e-1 #extrapolation linéaire de reff 5 sec (voir plot d'en bas pour la freq max de 2Mass) Reff = 7.73"
    
    freq_min2= []
    freq_max2= []
       
    #----------Filtes de 2MASS----------#
    
    #longueurs d'onde centrales des filtres de 2MASS [A]
    lambd_centr2=[2.16e4, 1.65e4, 1.25e4]
    
    #largeur a mi-hauteur des filtres Ks, H, J [A]
    D_lambd2= [2.2e3, 2.8e3, 3.0e3]
    
    
    
    for i in range (len(lambd_centr2)):
        freq_min2.append(convfreq(lambd_centr2[i] + D_lambd2[i]/2.))
        freq_max2.append(convfreq(lambd_centr2[i] - D_lambd2[i]/2.))


    for i in range(len(freq_min2)):
        freq2=[]
        flux2=[]
        for j in range (len(matrix)):          
            if (matrix [j][0][0])!= "#" and (float(matrix [j][0]) >= freq_min2[i]) and (float(matrix [j][0]) <= freq_max2[i]) :
                freq2.append(float(matrix [j][0]))
                flux2.append(float(matrix [j][2]))
        
        flux_pond2.append(moy(flux2)*ponderation2)
        flux_moy2.append(moy(flux2))
        
    print "flux de la galaxie hote dans les filtres Ks, H, J de 2MASS[erg.cm-2.s-1]: \n",flux_pond2
    print "flux moyen \n", flux_moy2
    print "flux corrigé 2MASS \n", array(mass)- array(flux_pond2), "\n"     
    

#----------calcul du flux de la galaxy hote dans les filtres de 2MASS----------#

if ATOM == 1:
    #ouverture ATOM ["]
    ouvA = 4.0
    #Pouv1 = ouv1/reffRc
    
#    freq_min2= []
#    freq_max2= []
       
    #----------Filtes de ATOM----------#
    
    #longueurs d'onde centrales du filtre Rc d'ATOM (these marcus hauser) [A]
    lambd_centr =  6250
    #nu_centr = c/(array(lambd_centr)*1e-10)
    nu_centr = convfreq(lambd_centr)
    ponderation= YOUNG(ouvA,Reff)
    
    
    #largeur a mi-hauteur du filtre R (these marcus hauser)[A]
    D_lambd= 1150
    
    
    freq_min = convfreq(lambd_centr + D_lambd/2.)
    freq_max = convfreq(lambd_centr - D_lambd/2.)
        
    #lecture des donnees PEGASE

    freq=[]
    flux=[]
    for j in range (len(matrix)):
        if (matrix [j][0][0])!= "#" and (float(matrix [j][0]) >= freq_min) and (float(matrix [j][0]) <= freq_max) :
            freq.append(float(matrix [j][0]))
            flux.append(float(matrix [j][2]))
            
    flux_pond = moy(flux)*ponderation
    flux_moy = moy(flux)
        
   
    print "flux moyen sur la bande R de ATOM sans prendre en compte l'ouverture\n", flux_moy,"\n"
    print "flux de la galaxie hote dans le filtre R de ATOM [erg.cm-2.s-1]: \n",flux_pond ,"\n"
    print "flux corrigé ATOM \n", atom - flux_pond, "\n" 
        
        
#----------Etalonnage de la galaxie hôte sur un filtre----------#      

if Abraham == 1:
    # etalonnage de Pegase sur la bande R kitt peak (abraham 1991)
    numin = 4.212e14
    numax = 5.102e14
    freqR = []
    fluxR = []
    for j in range (len(matrix)):
        if (matrix [j][0][0])!= "#" and (float(matrix [j][0]) >= numin) and (float(matrix [j][0]) <= numax) :
            freqR.append(float(matrix [j][0]))
            fluxR.append(float(matrix [j][2]))
    flux_moyR = (moy(fluxR))
    print "flux moyen sur la bande R-kitt-peak\n", flux_moyR

elif Rc == 1:
    #étalonnage de Pegase sur la bande R Cousin (scarpa 2000)
    numin = 3.997e14
    numax = 5.451e14
    freqR = []
    fluxR = []
    for j in range (len(matrix)):
        if (matrix [j][0][0])!= "#" and (float(matrix [j][0]) >= numin) and (float(matrix [j][0]) <= numax) :
            freqR.append(float(matrix [j][0]))
            fluxR.append(float(matrix [j][2]))
    flux_moyR = (moy(fluxR))
    print "flux moyen de Pegase sur la bande R Cousin \n", flux_moyR , "\n"
    print "flux moyen de Pegase sur la bande R Cousin pour une ouverture de 5 arcsec\n", flux_moyR*YOUNG(ouv0,Reff) , "\n"
    
elif R_Johnson == 1:
    #étalonnage de Pegase sur la bande R Johnson (Pursimo 2002)
    numin = convfreq(9200)
    numax = convfreq(4800)
    freqR = []
    fluxR = []
    for j in range (len(matrix)):
        if (matrix [j][0][0])!= "#" and (float(matrix [j][0]) >= numin) and (float(matrix [j][0]) <= numax) :
            freqR.append(float(matrix [j][0]))
            fluxR.append(float(matrix [j][2]))
    flux_moyR = (moy(fluxR))
    print "flux moyen de Pegase sur la bande R Johnson \n", flux_moyR , "\n"
    print "flux moyen de Pegase sur la bande R Johnson pour une ouverture de 5 arcsec\n", flux_moyR*YOUNG(ouv0,Reff) , "\n"
    

#print "B-V: \n", -2.5*log10(flux_moy[1])+2.5*log10(flux_moy[0]), "\n"
#print "V-R: \n", -2.5*log10(flux_moy_uvot[0])+2.5*log10(flux_moyR), "\n"  


elif Ks_SOFI == 1:
    #étalonnage de Pegase sur la bande Ks de SOFI
    numin = convfreq(22995)
    numax = convfreq(20245) #https://www.eso.org/sci/facilities/lasilla/instruments/sofi/inst/Imaging.html
    freqKs = []
    fluxKs = []
    for j in range (len(matrix)):
        if (matrix [j][0][0])!= "#" and (float(matrix [j][0]) >= numin) and (float(matrix [j][0]) <= numax) :
            freqKs.append(float(matrix [j][0]))
            fluxKs.append(float(matrix [j][2]))
    flux_moyKs = (moy(fluxKs))
    print "flux moyen de Pegase sur la bande Ks \n", flux_moyKs , "\n"
    print "flux moyen de Pegase sur la bande Ks pour une ouverture de 5 arcsec\n", flux_moyKs*YOUNG(ouv0,Reff) , "\n"




#flux corrigé UVOT 
f = array([  3.17346716e-12,1.74857449e-12,3.95466226e-12,3.35478668e-12,3.92352096e-12,3.09039856e-12])

fmax = array([  4.14142350e-12,2.29443079e-12,4.40818518e-12,3.75934584e-12,4.22038374e-12,3.35704361e-12])

dmax_fcorr = fmax-f





