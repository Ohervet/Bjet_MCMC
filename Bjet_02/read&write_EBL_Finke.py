# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 14:02:07 2016

@author: olivier
"""

from pylab import *
from os import path, access, R_OK
import csv

n = 0.00
font = []
for i in range(500):
    #print str(n)+',',
    font.append('tau_modelC_total/tau_modelC_total_z%3.2f.dat'% n)
    n += 0.01

E = [] #[TeV]
MTau = []
for i in range(500):
    Tau = []
    if path.isfile(font[i]) and access(font[i], R_OK):
        fich = open(font[i], "r")
        fich1= csv.reader(fich, delimiter = ' ', skipinitialspace = True)
        for line in fich1:
            if i == 0:
                E.append(float(line[0]))
                print str(E[-1])+',',
            Tau.append(float(line[1]))
        MTau.append(Tau)
        fich.close()


#
#if path.isfile(font2) and access(font, R_OK):
#    fich = open(font2, "r")
#    fich1= csv.reader(fich)
#    for line in fich1:
#        line_results.append(line)
#        line_results[-1] = line_results[-1][0].replace('.root','')
#        line_results[-1] = int(line_results[-1])
#    fich.close()
#    
with open('TauEBL.txt', 'wb') as csvfile:
    fich1 = csv.writer(csvfile,lineterminator="\n", delimiter ='a', quoting=csv.QUOTE_NONE)
    for i in range(len(E)):
        row = ''
        for j in range(500):
            row += str(MTau[j][i])+', '
        fich1.writerow([row])