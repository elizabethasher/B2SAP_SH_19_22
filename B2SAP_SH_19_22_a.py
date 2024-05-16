#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 14:44:27 2021

@author: easher
"""
#calculate climate-relevant forward scattering for particle size distribution
import PyMieScatt
import sys
import xarray
import os
import csv
import math
#from datetime import datetime as dt
#from datetime import timedelta
import datetime as dt
from datetime import datetime
from datetime import timedelta
from datetime import date
import time
import pandas as pd
import numpy as np
import numpy.matlib 
import netCDF4 as nc
import glob
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
import seaborn as seaborn
from scipy import stats
from scipy.optimize import fsolve
#select one size distribution...
 

#Use the look-up table to find the particle diameter and dLogDp calculated from POPS signal 
def DiaLookUpTable1(Composition, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15):
    file = '/Users/asher/Documents/NOAAweb/NZ_Paper_IOR1.csv' #'/Users/easher/Documents/NOAAweb/POPS_Sizes.csv'
    df1=pd.read_csv(file,sep=',', dtype=None, engine='python')

    if Composition == 'Smoky':
        dp_in = df1[(df1.Composition == 'Smoky') & (df1.Binning == 'Manual')]
        
    elif Composition == 'Smoky Auto':
        dp_in = df1[(df1.Composition == 'Smoky') & (df1.Binning == 'Automatic')]

        
    elif Composition == 'Sulfate Auto':
        dp_in = df1[(df1.Composition == 'Sulfate') & (df1.Binning == 'Automatic')]
       

    else:
        dp_in = df1[(df1.Composition == 'Sulfate') & (df1.Binning == 'Manual')]
        
    ndp = np.array([B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15])
    dp = np.array(dp_in.Diameter*1000).tolist()
    
    dp_w = WeightedBinMean(ndp, dp) #returns 15 bin weighted means instead of 16 bin edges (input)

    #print(dp_w)
    
    #del dp_in, dp
    #return dfp_w
    return dp_w


    
def WeightedBinMean(ndp, dp):
    dp_f = np.ones(15, dtype=float)
    dp_f[0] = (dp[0] + dp[1])/2
    #print(ndp)
    #print(dp)
    for num in range (1,15): 
        #print (num)
           
        if (dp[num]) == 0:
            dp_f[num] = (dp[num] + dp[num+1])/2
        elif (ndp[num - 1]) == 0:
            dp_f[num] = (dp[num] + dp[num+1])/2
        else:
            dp_f[num] = ((dp[num] * ndp[num - 1] + dp[num+1] * ndp[num])/ (ndp[num - 1] + ndp[num]))
        #except:
        
    #print(dp_f)

    return dp_f


def SimpleSigAndBinMean(dp, Sig):
    Sig_f = np.ones(15, dtype=float)
    dp_f = np.ones(15, dtype=float)
    
    for num in range (0,15): 
        #print (num)
        dp_f[num] = (dp[num] + dp[num+1])/2
        Sig_f[num] = (Sig[num] + Sig[num+1])/2

    return Sig_f, dp_f

#Use the look-up table to plot POPS signal vs diameter at two different IORs...
def plotIOR(df):#take one profile and calculate extinction with two different indices of refraction
    #test the outcome look at variable
    file = '/Users/asher/Documents/NOAAweb/NZ_Paper_IOR1.csv' #'/Users/easher/Documents/NOAAweb/POPS_Sizes.csv'
    df1=pd.read_csv(file,sep=',', dtype=None, engine='python')
    
    
    dp_in_A = df1[(df1.Composition == 'Smoky') & (df1.Binning == 'Manual')] 
    dp_A =  np.array(dp_in_A.Diameter*1000).tolist()
    signal_A =  np.array(dp_in_A.Signal).tolist()
    
    sig_A_f = np.ones(15, dtype=float)
    dp_A_f = np.ones(15, dtype=float)
    
    for num in range (0,15): 
        #print (num)
        dp_A_f[num] = (dp_A[num] + dp_A[num+1])/2
        sig_A_f[num] = (signal_A[num] + signal_A[num+1])/2
   
    dp_A = dp_A_f
    signal_A = sig_A_f
    
    dp_in_B = df1[(df1.Composition == 'Sulfate') & (df1.Binning == 'Manual')]
    dp_B = np.array(dp_in_B.Diameter*1000).tolist()
    signal_B =  np.array(dp_in_B.Signal).tolist()
    
    
    sig_B_f = np.ones(15, dtype=float)
    dp_B_f = np.ones(15, dtype=float)
    
    for num in range (0,15): 
        #print (num)
        dp_B_f[num] = (dp_B[num] + dp_B[num+1])/2
        sig_B_f[num] = (signal_B[num] + signal_B[num+1])/2
   
    dp_B = dp_B_f
    signal_B = sig_B_f
    
    
    
    plots = 1
    if plots == 1:
        #figure 1a
        fig = plt.figure()
        plt.plot(signal_A, dp_A, 'rs-')
        plt.plot(signal_B, dp_B, 'ks-')
        plt.xlabel('POPS signal (ADC)')
        plt.ylabel('Mean POPS Bin Diameter (nm)')
        plt.legend(['Smoke', 'Baseline (Sulfate)'])
        plt.yscale('log')
        plt.xscale('log')
        plt.grid()
        
        #figure 1b & 1c
        #select one profile
        df_select = df_AscDesc[df_AscDesc['Launch'] == '2020-07-23'] #this could also be checked out with POPS from 08/31/2021
        Title = '2020-07-23'
        df_select['Composition_A'] = 'Smoky'
        df_select['Composition_B'] = 'Sulfate'
        
        df_select[['Sim_SA_A', 'Sim_SA_B']] = df_select.apply(lambda x: pd.Series(CalcSA(x['Composition_A'], x['Composition_B'], x['B1'], x['B2'], x['B3'], x['B4'], x['B5'], x['B6'], x['B7'], x['B8'], x['B9'], x['B10'], x['B11'], x['B12'], x['B13'], x['B14'], x['B15'])), axis=1)
        df_select[['Sim Extinction_532 A', 'Sim Ambient Extinction_532 A']] = df_select.apply(lambda x: pd.Series(VPsca(x['Wavelength1'], x['Composition_A'], x['RH'], x['Air Temp (K)'], x['H2SO4_%wt_pops'], x['Density_pops'], x['H2SO4_%wt_amb'], x['Density_amb'],x['B1'], x['B2'], x['B3'], x['B4'], x['B5'], x['B6'], x['B7'], x['B8'], x['B9'], x['B10'], x['B11'], x['B12'], x['B13'], x['B14'], x['B15'])), axis=1)
        df_select[['Sim Extinction_532 B', 'Sim Ambient Extinction_532 B']] = df_select.apply(lambda x: pd.Series(VPsca(x['Wavelength1'], x['Composition_B'], x['RH'], x['Air Temp (K)'], x['H2SO4_%wt_pops'], x['Density_pops'], x['H2SO4_%wt_amb'], x['Density_amb'], x['B1'], x['B2'], x['B3'], x['B4'], x['B5'], x['B6'], x['B7'], x['B8'], x['B9'], x['B10'], x['B11'], x['B12'], x['B13'], x['B14'], x['B15'])), axis=1)
        
        
        df_select = df_select.filter(['Launch', 'Altitude (km)', 'Composition_A', 'Composition_B', 'Sim_SA_A', 'Sim_SA_B', 'Sim Extinction_532 A', 'Sim Ambient Extinction_532 A', 'Sim Extinction_532 B', 'Sim Ambient Extinction_532 B', ])
        df_select['% Difference relative to Sulfate Case'] = (df_select['Sim Extinction_532 A'] - df_select['Sim Extinction_532 B'])/ df_select['Sim Extinction_532 B']*100
    
        
        #Figure 1b
        fig = plt.figure()
        plt.plot(df_select['Sim_SA_A'], df_select['Altitude (km)'], 'r-')
        plt.plot(df_select['Sim_SA_B'], df_select['Altitude (km)'], 'k--')
        plt.xlabel('Surface Area')
        plt.ylabel('Altitude (km)')
        plt.legend(['Smoke', 'Sulfate'])
        plt.ylim([12.0, 28.0])
        plt.xlim([-1, 25])
        plt.title(Title)
        
        [slope_a, intercept_a, r_value_a, p_value_a, std_err_a] = stats.linregress(df_select['Sim_SA_B'], df_select['Sim_SA_A'])
        fig = plt.figure()
        plt.plot(df_select['Sim_SA_B'], df_select['Sim_SA_A'], 'ko')
        line_x = np.linspace(0.5, 15, num=100)
        plt.plot(line_x, line_x*slope_a+intercept_a, 'b-')
        plt.xlabel('Sulfate (Baseline) Surface Area')
        plt.ylabel('Smoke Surface Area')
        plt.legend([str(slope_a)]) #['y = 1.38x + -0.41'])
        plt.title(Title)
        
        
        #Figure 1c
        fig = plt.figure()
        plt.plot(df_select['Sim Extinction_532 A'], df_select['Altitude (km)'], 'r-')
        plt.plot(df_select['Sim Extinction_532 B'], df_select['Altitude (km)'], 'k--')
        plt.xlabel('Extinction 532 nm')
        plt.ylabel('Altitude (km)')
        plt.legend(['Smoke', 'Sulfate'])
        plt.ylim([12.0, 28.0])
        plt.xlim([-0.0001, 0.010])
        plt.title(Title)
        
        
        #df_select = df_select[df_select['Sim Extinction_532 B'] < 0.01]
        
        [slope_b, intercept_b, r_value_b, p_value_b, std_err_b] = stats.linregress(df_select['Sim Extinction_532 B'], df_select['Sim Extinction_532 A'])
        fig = plt.figure()
        plt.plot(df_select['Sim Extinction_532 B'], df_select['Sim Extinction_532 A'], 'ko')
        line_x = np.linspace(0.00001, 0.010, num=100)
        plt.plot(line_x, line_x*slope_b+intercept_b, 'b-')
        plt.xlabel('Sulfate (Baseline) Extinction')
        plt.ylabel('Smoke Extinction')
        plt.legend([str(slope_b)]) #['y = 1.20x + -1.6E-5'])
        plt.title(Title)
    
    return 

 
  
def DnDlogDp(Launch, Composition, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15):
    
    file_ior = '/Users/asher/Documents/NOAAweb/NZ_Paper_IOR1.csv'
    dfIOR = pd.read_csv(file_ior,sep=',', dtype=None, engine='python')
    #dLogDp = DLogDpLookUpTable(Composition)
    
    if Composition == 'Smoky':
        dp_in = dfIOR[(dfIOR.Composition == 'Smoky') & (dfIOR.Binning == 'Manual')]
        dLogDp = np.array(dp_in.dLogDp)

        ndp = np.array([B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15])


        #dNdLogDp = np.divide(np.array([B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15, B16, B17, B18, B19, B20]), np.array(dLogDp)).tolist()
    elif Composition == 'Smoky Auto':
        dp_in = dfIOR[(dfIOR.Composition == 'Smoky') & (dfIOR.Binning == 'Automatic')]
        dLogDp = np.array(dp_in.dLogDp)

        ndp = np.array([B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15])
        #dNdLogDpA = np.divide(ndp, dLogDp)
            
    elif Composition == 'Sulfate Auto':
        dp_in = dfIOR[(dfIOR.Composition == 'Sulfate') & (dfIOR.Binning == 'Automatic')]
        dLogDp = np.array(dp_in.dLogDp)

        ndp = np.array([B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15])
        #dNdLogDpA = np.divide(ndp, dLogDp)
            
    else:
        dp_in = dfIOR[(dfIOR.Composition == 'Sulfate') & (dfIOR.Binning == 'Manual')]
        dLogDp = np.array(dp_in.dLogDp)

        ndp = np.array([B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15])
        #dNdLogDpA = np.divide(ndp, dLogDp)
        
        # for x in range(5):
        #     dNdLogDpA = np.append(dNdLogDpA, [np.nan])
        #     #print(x)
    #print(dNdLogDpA)  
    #print(dLogDp)
    dLogDp_f = dLogDp[0:15]
    dNdLogDpA = np.divide(ndp, dLogDp_f)
    dNdLogDp = dNdLogDpA.tolist()   
    return dNdLogDp

def CalcAerConc(Choice, T, P, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15):
    #constants
    Rgas = 287.058 #J K-1 kg-1 or Pa m3 K-1 kg-1
    R_unv = 8.3144 #m3 Pa K-1 mol -1
    convFac = 1E6 #g to ug convert mass of aerosols
    Pres = P * 100  #convert Pressure hPa to Pa
    Rho_gas_STP = 1.225 #kg/m
    
    #calc STP choice
    AtSTP = 0 #if 1, calculate #/cm3 at STP, if 0, calculate #/cm3 ambient
        
    if Choice == 1: #do for all bins
        
        if AtSTP == 0:
            ndp = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15] 
            Conc = pd.Series(ndp)
            AerConc_f = Conc.sum()
        else:
            ndp = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15] 
            Conc = pd.Series(ndp)
            AerConc_i = Conc.sum()
            AerConc_i2 = AerConc_i * Rgas * T / (Pres)
            AerConc_f = AerConc_i2 * Rho_gas_STP 
    
    if Choice == 2: #do for all bins

            if AtSTP == 0:
                ndp = [B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15] 
                Conc = pd.Series(ndp)
                AerConc_f = Conc.sum()
            else:
                ndp = [B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15] 
                Conc = pd.Series(ndp)
                AerConc_i = Conc.sum()
                AerConc_i2 = AerConc_i * Rgas * T / (Pres)
                AerConc_f = AerConc_i2 * Rho_gas_STP 
        
    return AerConc_f



def Calc_STP_AB(T, P, PartAmb):
    #constants
    Temp = T + 273.15
    Rgas = 287.058 #J K-1 kg-1 or Pa m3 K-1 kg-1
    R_unv = 8.3144 #m3 Pa K-1 mol -1
    convFac = 1E6 #g to ug convert mass of aerosols
    Pres = P * 100  #convert Pressure hPa to Pa
    Rho_gas_STP = 1.225 #kg/m
    
    #calc Particle concentration OR surface area at STP using ambient particle concentration or surface area

    Part_i = PartAmb * Rgas * Temp / (Pres)
    Part_STP = Part_i * Rho_gas_STP 

    return Part_STP


    
   
def VPscaCARMA(B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15, B16, B17, B18, B19, B20):
   # 385, 448, 521, 596, 676, 754, 868, 1019, and ~1550 nm
    #wavelength =521 #nm
    Wavelength = 532
    Composition = 'Smoky'
    dp = [100, 131.2, 172.2, 226, 296, 388, 510, 670, 878, 1152,
                      1512, 1984, 2603.22, 3415.62, 4481.56, 5880.16, 7715.22, 10123, 13282.1, 17427.1]
    ndp = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15, B16, B17, B18, B19, B20] 
    
    IR = 1.54+0.03j
    #dpAmbnp_a = KKcalcDamb(RH, T, Composition)

    #IR = 1.45+0.00j
        
    try:
        [Bext_m, Bsca_m, Babs_m, bigG_m, Bpr_m, Bback_m, Bratio_m] = PyMieScatt.Mie_SD(IR, Wavelength, dp, ndp, nMedium=1.0, interpolate=True, SMPS = True,  asDict=False) 
    except:
        #print(dp)
        Bext_m = np.NaN
    # try:
    #     [Bext_K, Bsca_K, Babs_K, bigG_K, Bpr_K, Bback_K, Bratio_K] = PyMieScatt.Mie_SD(IR, Wavelength, dpAmbnp_a, ndp, nMedium=1.0, interpolate=True, SMPS = True,  asDict=False) 
    # except:
    #     print(dpAmbnp_a)
    #     Bext_K = np.NaN

    #convert Mm-1 to km -1
    B_m = Bext_m/1000
    #B_K = Bext_K/1000

    
    return B_m

def VPsca(Wavelength, Composition, RH, T,  Fpops, rho_pops, Famb, rho_amb, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15):
   # 385, 448, 521, 596, 676, 754, 868, 1019, and ~1550 nm
    #wavelength =521 #nm
    dp_m = DiaLookUpTable1(Composition, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15)  
    ndp = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15] 
    
    

    
    #return backscattering if Bback = 1 and Extinction if Bback = 0
    Bback = 0
    
    # if LaunchNo == 'run004':
    #     ndp = [B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15, B16, B17, B18, B19, B20]
    #     dpAmbnp_K = dpAmbnp_K[2:]
    #     dp_m = dp_m[2:]
    
    if (Composition == 'Smoky'):
        #Kchem = 0.17 # for H2SO4
        IR = 1.54+0.18j
        #IR = 1.40+0.01j
        dpAmbnp_a = KKcalcDamb(RH, T, Composition,  B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15)

    elif (Composition == 'Smoky Auto'):
        #Kchem = 0.17 # for H2SO4
        IR = 1.54+0.18j
        #IR = 1.40+0.01j
        dpAmbnp_a = KKcalcDamb(RH, T, Composition, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15)

    elif (Composition == 'Sulfate'): 
        #Kchem = 0.87 # for H2SO4
        IR = 1.45+0.00j
        dpAmbnp_a = Steele_Hamill_1981_J1995_calcDamb(Composition, Fpops, rho_pops, Famb, rho_amb, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15) #Steele and Hamill, 1981 method, also used by Jonsson et al. 1995. F is % wt., rho is density.
        

    else:
        #Kchem = 0.87 # for H2SO4
        IR = 1.45+0.00j
        dpAmbnp_a = Steele_Hamill_1981_J1995_calcDamb(Composition, Fpops, rho_pops, Famb, rho_amb, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15) #Steele and Hamill, 1981 method, also used by Jonsson et al. 1995. F is % wt., rho is density.

    #m = 1.54+0.01j #all smoky assumed
    #dp = [133, 145, 160, 179, 203,  235, 324, 456, 552, 684, 1348, 1663, 2541, 3000, 3500]
    ndp = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15]
    #[Bext, Bsca, Babs, bigG, Bpr, Bback, Bratio] = PyMieScatt.Mie_SD(m, Wavelength, dp, ndp, nMedium=1.0, interpolate=False, SMPS = True,  asDict=False)
    try:
        [Bext_m, Bsca_m, Babs_m, bigG_m, Bpr_m, Bback_m, Bratio_m] = PyMieScatt.Mie_SD(IR, Wavelength, dp_m, ndp, nMedium=1.0, interpolate=True, SMPS = True,  asDict=False) 
    except:
        print(dp_m)
        Bext_m = np.NaN
    try:
        [Bext_K, Bsca_K, Babs_K, bigG_K, Bpr_K, Bback_K, Bratio_K] = PyMieScatt.Mie_SD(IR, Wavelength, dpAmbnp_a, ndp, nMedium=1.0, interpolate=True, SMPS = True,  asDict=False) 
    except:
        print(dpAmbnp_a)
        Bext_K = np.NaN
    #print(Bext)
    #print(Bsca)
    if Bback == 0:
        #convert Mm-1 to km -1
        B_m = Bext_m/1000
        B_K = Bext_K/1000
    else:
        #convert Mm-1 to km -1
        B_m = Bback_m/1000
        B_K = Bback_K/1000
    
    return B_m, B_K

def Steele_Hamill_1981_J1995_calcDamb(Composition, Fpops, rho_pops, Famb, rho_amb, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15):
#calculate ambient aerosol diameter for H2SO4 specifically 
#This calculation is identical to method using Steele & Hamill 1981 (above). It assumes that particles are at equilibrium both at the time of measurement as well as in ambient air,
#and that only H2O, not H2SO4 is lost from particles during sampling. Also used by Jonsson et al., 1995 J. of Atmos. & Oceanic Tech. Vol. 12 (1), as discussed. in Asher et al.
    
    dp = DiaLookUpTable1(Composition, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15)
    dpAmb = []
        
    for i in dp:
        Dpops = i
        Damb = Dpops*(Fpops*rho_pops/(Famb*rho_amb))**(1/3)
        dpAmb.append(Damb)
    dpAmbn= np.array(dpAmb).tolist()
       
    return dpAmbn

def PwLookuptable(T1, Pwa):
    #Lookup table for H2SO4 aerosol %wt given partial pressure of H2O and T (K) 
    #Tabazadeh et al., 1997 GRL vol. 24 (15)
    
    Pw = np.round_(Pwa, decimals = 5)
    #print('Lookup table Pw : ' + str(Pw))
    
    T = T1 + 273.15 #can interpolate between table values..
    #Pwa must be in hPa/mb
    Wi = np.linspace(0.1, 0.8, 15) #Tabazadeh et al., 1997 GRL vol. 24 (15) Table 1.
    W_0 = np.linspace(0.1, 0.8, 701) #71 to the 10's place (or nearest % of % wt)
    ai = [19.726, 19.747, 19.761, 19.794, 19.883, 20.078, 20.379, 20.637, 20.682, 20.55, 20.405, 20.383, 20.585, 21.169, 21.808 ]
    bi = [-4364.8, -4390.9, -4414.7,  -4451.1, -4519.2, -4644.0, -4828.5, -5011.5, -5121.3, -5177.6, -5252.1, -5422.4, -5743.8, -6310.6, -6985.9 ]
    ci = [-147620, -144690, -142940, -140870, -136500, -127240, -112550, -98811, -94033, -96984, -100840, -97966, -83701, -48396, -12170]
    a = pd.Series(np.interp(W_0, Wi, ai))
    b = pd.Series(np.interp(W_0, Wi, bi))
    c = pd.Series(np.interp(W_0, Wi, ci))
    W_1 = pd.Series(W_0)
    
    PwTable = pd.DataFrame({'a': a, 'b': b, 'c': c, 'W_1': W_1})#(a, b, c, W_1)
    PwTable['T'] = T
    
    PwTinit = pd.DataFrame({'a': ai, 'b': bi, 'c': ci, 'W': Wi})
    PwTinit['T'] = T
    
    PwTinit['Pw'] = PwTinit.apply(lambda x: PwCalc(x['a'], x['b'], x['c'], x['T']), axis=1)
    PwTable['Pw_calc'] = PwTable.apply(lambda x: PwCalc(x['a'], x['b'], x['c'], x['T']), axis=1)
    #create dictionary and look up H2SO4 weight fraction (apply to mass)
    PwTable.drop(['a', 'b', 'c', 'T'], axis = 1, inplace=True)
    PwDict = pd.Series(PwTable.Pw_calc.values,index=PwTable.W_1).to_dict()
    #find the closest dictionary look-up value
    W_T, Pw_t = min(PwDict.items(), key=lambda x: abs(Pw - x[1]))
    W = W_T
    
    #if the wt% is greater than or equal to 0.8, use the Gmitro and Vermeulen  1964 parameterization eq. 16; partial molar properties are listed in Table 3.
    if W_T >= 0.8:
        W_GV = pd.Series([0.80, 0.82, 0.84, 0.86, 0.88, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98]) #maybe go up to 0.98 W?
        R = 1.98726
        Cp_298 = pd.Series([9.77, 10.36, 13.78, 18.96, 22.13, 22.76, 22.30, 21.48, 20.44, 19.32, 18.06, 16.64, 15.05, 13.25])
        L_298 = pd.Series([-3475, -4015, -4656, -5319, -5938, -6419, -6627, -6816, -6983, -7139, -7286, -7433, -7574, -7712])
        F_F0 = pd.Series([-3090, -3427, -3789, -4167, -4557, -4960, -5165, -5375, -5595, -5830, -6090, -6390, -6741, -7204])
        alpha = pd.Series([0.0114, 0.0568, 0.1233, 0.0666, -0.0120, -0.0346, -0.0398, -0.0427, -0.0436, -0.0428, -0.0405, -0.0368, -0.0314, -0.024])
        
        A1 = -3.67340
        B1 = -4143.5
        C1 =  10.24353
        D1 = 0.618943E-3
        E = 0
        #calculate the constants based on constants above and Temperature   
        A = A1 + 1/R * (Cp_298 - 298.15 * alpha)
        B = B1 + 1/R * (L_298 - 298.15 * Cp_298 + 298.15**2/2 * alpha)
        C = C1 + 1/R * (Cp_298 + (F_F0 - L_298) * 1/298.15)
        D = D1 - alpha/(2 * R)
        #calculate Pw
        
        PwTable_GV = pd.DataFrame({'A': A, 'B': B, 'C': C, 'D': D, 'E': E, 'W': W_GV})#
        PwTable_GV['T'] = T
        
        PwTable_GV['pH2O'] = PwTable_GV.apply(lambda x: PwCalc_GV1964(x['A'], x['B'], x['C'], x['D'], x['E'], x['T']), axis = 1)                          
        
        PwTable_GV.drop(['A', 'B', 'C', 'D', 'E', 'T'], axis = 1, inplace=True)
        PwDict_GV = pd.Series(PwTable_GV.pH2O.values,index=PwTable_GV.W).to_dict()
        W_GV, Pw_t = min(PwDict_GV.items(), key=lambda x: abs(Pw - x[1]))
        W = W_GV

    return W


def PwCalc(a, b, c, T):
    #use polynomial to calc. eq. partial pressure of H2O given T
    #Tabazadeh et al., 1997 GRL vol. 24 (15) Table 1. eq
        Pw_0 = math.exp(a + b/T + c/(T**2))
        Pw = np.round_(Pw_0, decimals = 5)
        return Pw
    
def PwCalc_GV1964(A,B,C,D,E,T):
    #Gmitro and Vermeulen 1964 AIChE Journal vol. 10 (5) pg. 740 746 eq. 16
        Pwlog = A * np.log(298.15/T) + B/T + C + D*T + E * np.power(T,2)
        Pw = np.exp(Pwlog) #units of partial pressure in atmospheres
        Pwa = Pw * 1013.25#conversion factor from atmospheres to mb
        return Pwa


def DensityCalc(perWt, Tc):
#density parameterization Oca et al., 2018 J. Chem. & Eng. Data vol. 63 (9)
        T = Tc + 273.15 #convert C to K
        PolyCoeff = [1022, -0.5076, 2.484E-4, 976.4, -1.015, 1, 237.8, 1, 1]
        Density = PolyCoeff[0] + PolyCoeff[1] * T + PolyCoeff[2] * T**2 + PolyCoeff[3] * perWt + PolyCoeff[4] + PolyCoeff[5] * perWt**2 
        Density_c = Density / 1000 #convert from kg/m3 to g/cm3
        Density_unc = 0.02 #uncertainty reported 
        
        if perWt >= 0.68:
            #use Washburn analytical tables of chemistry 1928 (0 C) similar to POPS internal temperatures ~ -5 C to +5 C when percent weight 
            Wi = np.linspace(0.68, 0.98, 31)
            W = pd.Series(Wi)
            rho = pd.Series([1.6059, 1.6176, 1.6293, 1.6411, 1.6529, 1.6648, 1.6768, 1.6888, 1.7008, 1.7128, 1.7247, 1.7365, 1.7482, 1.7597,
                1.7709, 1.7815, 1.7916, 1.8009, 1.8095, 1.8173, 1.8243, 1.8306, 1.8361, 1.8410, 1.8453, 1.8490, 1.8520, 1.8544, 1.8560, 1.8569, 1.8567])
          
            WDensityTable = pd.DataFrame({'rho': rho, 'W': W})#

            WDensityDict = pd.Series(WDensityTable.W.values,index=WDensityTable.rho).to_dict()
            
            Density_W, W = min(WDensityDict.items(), key=lambda x: abs(perWt - x[1]))
            
            #if W >0.6
            p1 = 473.52 + 4903.99 * W -11916.50 * np.power(W,2) + 15057.60 * np.power(W,3) - 6668.37 * np.power(W,4)
            p2 =  250.52 +  5733.14* W - 13138.14 * np.power(W,2) + 15565.78 * np.power(W,3) - 6618.70* np.power(W,4)
            
            Density = p1 + (p2-p1) * ((T - 273.15)/69) #density in kg/m3 want it in g cm-3
            Density_c = Density/1000
        
            Density_unc = (Density_c - Density_W)/Density_c
        return Density_c, Density_unc
    
def GoffGratch(TFpHyg, Temp):
    #GoffGratch 1984 formulation
    TFpHyg = TFpHyg + 273.15
    Temp = Temp + 273.15
    #RH = VaporPressLiquid(hygroFrostpoint) /VaporPressLiquid(temperature) * 100.0;
    #RH over ice uses Goff Gratch saturation vapor pressure from http://cires.colorado.edu/~voemel/vp.html
    #RH over liquid water uses Hyland and Wexler saturation vapor pressure over liquid water from http://cires.colorado.edu/~voemel/vp.html
    
    
    if (Temp < 273.15):
        e =  10**(-9.09718*(273.16/TFpHyg - 1) - 3.56654*math.log10(273.16/TFpHyg) + 0.876793*(1-TFpHyg/273.16) + math.log10(6.1071))
        es = 10**(-9.09718*(273.16/Temp - 1) -   3.56654*math.log10(273.16/Temp) +   0.876793*(1-Temp/273.16) +   math.log10(6.1071))
        RH = e/es*100

    else:
        e = 10**(-9.09718*(273.16/TFpHyg - 1) - 3.56654*math.log10(273.16/TFpHyg) + 0.876793*(1-TFpHyg/273.16) + math.log10(6.1173))
        es = (math.exp(-0.58002206*10**4/Temp + 0.13914993*10**1 - 0.48640239*10**(-1)*Temp + 0.41764768*10**(-4)*Temp**2 -0.14452093*10**(-7)*Temp**3 + 0.65459673*10**1*math.log(Temp)))/100
        RH = e/es*100
        
    return (RH, e)
    
    

def KKcalcDamb(RH, T, Composition, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15):
    #assume something else for smoky aerosol - differnt Kchem??
    dp = DiaLookUpTable1(Composition, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15) #assuming Diameter with RH = 0% (at equilibirum)
    dpAmb = []
    
    if (Composition == 'Smoky'):
        Kchem = 0.18 # for organic aerosol
        #Kchem = 0.0
    elif (Composition == 'Smoky Auto'):
        Kchem = 0.18 # for organic aerosol
        #Kchem = 0.0
    elif (Composition == 'Sulfate'): 
        Kchem = 0.87 # for H2SO4
    else:
        Kchem = 0.87 # for H2SO4


    for i in dp:
        Dd = i
        func = lambda Drh : (Drh**3 - Dd**3)/ (Drh**3 - Dd**3*(1-Kchem))*np.exp(4*0.072*0.018/(8.314*T*1000*Drh))-RH/100
        Drh_initial_guess = Dd # this is true unless RH is very high and the particle is very small... then its size may change a lot
        Drh_solution = fsolve(func, Drh_initial_guess)
        dpAmb.append(Drh_solution[0])
    dpAmbn= np.array(dpAmb).tolist() #calculated aerosol Diameter with ambient RH and T)
    
    
    
    #Diameter difference (not plotted here)
    DpDiff = pd.DataFrame({'Dp':dp, 'dpAmbn':dpAmbn})
    DpDiff['Diff'] =  DpDiff['dpAmbn'] - DpDiff['Dp']
    DpDiff['Ratio'] = DpDiff['dpAmbn'] / DpDiff['Dp']

    DpDiff['RH'] = RH

    return dpAmbn

#calculates dry surface area
def CalcVolume(Composition, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15):
    dp = DiaLookUpTable1(Composition, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15)
    #print(dp)
    ndp = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15] 
            
    Conc = pd.Series(ndp)
    Rad =  pd.Series(np.divide(dp,2))/1000
    Vol = np.pi*4/3*np.power(Rad,3)
    DryVol = (Conc.mul(Vol)).sum()
    
    return DryVol #units are um3/cm3

def CalcDrySA(T, P, Composition, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15):
    #constants
    Rgas = 287.058 #J K-1 kg-1 or Pa m3 K-1 kg-1
    R_unv = 8.3144 #m3 Pa K-1 mol -1
    convFac = 1E6 #g to ug convert mass of aerosols
    Pres = P * 100  #convert Pressure hPa to Pa
    Rho_gas_STP = 1.225 #kg/m
    
    dp = DiaLookUpTable1(Composition, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15)
    ndp = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15] 
    Conc = pd.Series(ndp)
    
    #calc choice
    AtSTP = 0 #if 1, calculate #/cm3 at STP, if 0, calculate #/cm3 ambient
    
    if AtSTP == 0:
        Rad_A =  pd.Series(np.divide(dp,2))/1000
        SA_A = np.pi*4*np.power(Rad_A,2)
        AerSA_f = (Conc.mul(SA_A)).sum()
    else:
        Rad_A =  pd.Series(np.divide(dp,2))/1000
        SA_A = np.pi*4*np.power(Rad_A,2)
        DrySA_A = (Conc.mul(SA_A)).sum()

        AerSA_i2 = DrySA_A * Rgas * T / (Pres)
        AerSA_f = AerSA_i2 * Rho_gas_STP 
    return AerSA_f


def CalcSA(Composition_A, Composition_B, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15):
    dp_A = DiaLookUpTable1(Composition_A, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15)
    dp_B = DiaLookUpTable1(Composition_B, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15)
    
    ndp = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15] 
    Conc = pd.Series(ndp)
    
    Rad_A =  pd.Series(np.divide(dp_A,2))/1000
    SA_A = np.pi*4*np.power(Rad_A,2)
    DrySA_A = (Conc.mul(SA_A)).sum()
    
    Rad_B =  pd.Series(np.divide(dp_B,2))/1000
    SA_B = np.pi*4*np.power(Rad_B,2)
    DrySA_B = (Conc.mul(SA_B)).sum()
    
    return DrySA_A, DrySA_B #units are um3/cm3

#calculate effective radius. computationally intensive
def eRadCalc(Composition, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15):

    dp = DiaLookUpTable1(Composition, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15)
    ndp = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15] 
        
    Conc = pd.Series(ndp)
    Rad =  pd.Series(np.divide(dp,2))
    r3 = np.power(Rad,3)
    r2 = np.power(Rad,2)
    m3= (Conc.mul(r3)).sum()
    m2= (Conc.mul(r2)).sum()
    Rad = m3/m2

    return Rad

def Plume(Launch, Flight, Altitude):
    if (Launch == '2020-01-27'):
        # if ((Flight == 'Ascent') and (Altitude == 16.5)):
        #     composition = 'Smoky Auto'
        # elif ((Flight == 'Ascent') and (Altitude == 16.0)):
        #     composition = 'Smoky Auto'
        # elif ((Flight == 'Ascent') and (Altitude == 13.0)):
        #     composition = 'Sulfate Auto'
        # elif ((Flight == 'Ascent') and (Altitude == 12.5)):
        #     composition = 'Sulfate Auto'
        # elif ((Flight == 'Ascent') and (Altitude == 12.0)):
        #     composition = 'Sulfate Auto'
        # elif ((Flight == 'Descent') and (Altitude == 15.5)):
        #    composition = 'Smoky Auto'
        if ((Altitude >= 15.0) and (Altitude <= 17.5)):
            composition = 'Smoky Auto'
        #elif ((Altitude >= 11.75) and (Altitude <= 14.0)):
        #            composition = 'Smoky Auto'
        else:
            composition = 'Sulfate Auto'
    elif (Launch == '2019-04-30'):
        composition = 'Sulfate Auto'
    elif (Launch == '2019-09-03'):
        composition = 'Sulfate Auto'
    elif (Launch == '2020-06-19'):
        composition = 'Smoky Auto'
    elif (Launch == '2020-07-23'):
        composition = 'Smoky'
    elif (Launch == '2020-11-17'):
        composition = 'Sulfate'
    elif (Launch == '2021-08-31'):
        composition = 'Sulfate'
    elif (Launch == '2021-11-03'):
        composition = 'Sulfate'
#    elif (Launch == '2020-06-19'):
#        composition = 'Sulfate Auto'
    else:
        composition = 'Sulfate'
    return composition


def remove_char(lst, char):
    if not lst:
        return []
    else:
        return [lst[0].replace(char, '')] + remove_char(lst[1:], char)
    
def FLE(path):
    
    FLE_colnames = ['Time [min]', 'Time [sec]', 'Time GMT [hh:mm:ss GMT]','Press [hPa]', 'Praw [hPa]', 'RS P Flag', 'Alt  [km]', 
                    'Tcorr [deg C]', 'Temp [deg C]', 'Traw [deg C]', 'RS T Flag', 'Theta [K]', 'RH [%]', 'RHraw [%]', 'RS U Flag', 
                    'TFp V [deg C]', 'RS H2O Mr [ppmv]', 'IPW V [mm]', 'TVaisI [deg C]', 'RS Bat [V]', 'RS P Temp [deg C]', 
                    'RS U Temp [deg C]', 'Rise Rate [m/s]', 'O3 Cell I [uA]', 'O3 P [mPa]', 'O3_ppmv_r', 'O3 Flag', 'T Pump [deg C]',
                    'T Pump Raw [deg C]','O3Bat [V]', 'I Pump [mA]', 'Total Column O3 [DU]', 'Total w/ Extrap O3 [DU]', 'Total w/ SBUV O3 [DU]',
                    'Current Uncert', 'O3 Conversion Uncert', 'O3 Flowrate Uncert', 'O3 Pump Temp Uncert', 'O3 Total Uncert', 
                    'GPS lat  [deg]','Lon_deg','GPS alt [km]', 'GPS Flag', 'GPS sats', 'GPS Time [hh:mm:ss GMT]', 'GPS Press [mb]', 'GPS Heading [deg]', 
                    'GPS Elev Angle [deg]', 'Wind Speed [m/s]','Wind Dir [deg]','TFp Hyg [deg C]', 'H2O_ppmv_r','Fl', 'Co', 'RH FP [%]','RHI FP [%]',
                    'IPW FP [mm]','Frost [ADC]', 'Filt Frost [ADC]','Hyg Sun [ADC]','Hyg Heat [PWM]','Hyg FPT ADC [ADC]', 'Hyg Opt Temp [deg C]',
                    'Hyg Opt ADC [ADC]','Hyg Opt Heat [PWM]','Hyg Press [mb]', 'Hyg Geo Alt[km]', 'Hyg P Temp [deg C]', 'Hyg Avg FPT [deg C]',
                    'Hyg Bat [V]', 'TP2 [V]']
    
    VertRes = 100 #250 m bins like lev files   
    if VertRes == 250: #data binned to to every 250 m
        AltLabels = ["0.25", "0.5","0.75", "1", "1.25", "1.5", "1.75", "2", "2.25", "2.5", "2.75", "3", "3.25", "3.5", "3.75", "4","4.25", "4.5", "4.75", "5", "5.25", "5.5", "5.75", "6",\
                  "6.25", "6.5", "6.75", "7", "7.25", "7.5", "7.75", "8", "8.25", "8.5", "8.75", "9", "9.25", "9.5", "9.75", "10",\
                  "10.25", "10.5", "10.75", "11", "11.25", "11.5", "11.75", "12", "12.25", "12.5", "12.75", "13", "13.25", "13.5", "13.75", "14",\
                      "14.25", "14.5", "14.75", "15", "15.25", "15.5", "15.75", "16", "16.25", "16.5", "16.75", "17", "17.25", "17.5", "17.75", "18",\
                          "18.25", "18.5", "18.75", "19", "19.25", "19.5", "19.75", "20", "20.25", "20.5", "20.75", "21", "21.25", "21.5", "21.75", "22",\
                              "22.25", "22.5", "22.75", "23", "23.25", "23.5", "23.75", "24", "24.25", "24.5", "24.75", "25", "25.25", "25.5", "25.75", "26",\
                                  "26.25", "26.5", "26.75", "27", "27.25", "27.5", "27.75", "28", "28.25", "28.5", "28.75", "29",\
                          "29.25", "29.5", "29.75", "30.00"]
        AvertBins = np.linspace(0.125, 30.125, 121)
        
    elif VertRes == 100: #data binned to to every 100 m
        AltLabels = ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0", "1.1","1.2","1.3","1.4", "1.5", "1.6", "1.7", "1.8", "1.9", "2.0", "2.1", "2.2","2.3","2.4","2.5","2.6", "2.7", "2.8", "2.9", "3.0",
                "3.1", "3.2", "3.3", "3.4", "3.5", "3.6", "3.7", "3.8", "3.9", "4.0", "4.1","4.2","4.3","4.4", "4.5", "4.6", "4.7", "4.8", "4.9", "5.0", "5.1", "5.2","5.3","5.4","5.5","5.6", "5.7", "5.8", "5.9", "6.0",
                  "6.1", "6.2", "6.3", "6.4", "6.5", "6.6", "6.7", "6.8", "6.9", "7.0", "7.1","7.2","7.3","7.4", "7.5", "7.6", "7.7", "7.8", "7.9", "8.0", "8.1", "8.2","8.3","8.4","8.5","8.6", "8.7", "8.8", "8.9", "9.0",
                  "9.1", "9.2", "9.3", "9.4", "9.5", "9.6", "9.7", "9.8", "9.9", "10.0", "10.1","10.2","10.3","10.4", "10.5", "10.6", "10.7", "10.8", "10.9", "11.0", "11.1", "11.2", "11.3", "11.4", "11.5", "11.6", "11.7", "11.8", "11.9", "12.0",
                  "12.1", "12.2", "12.3", "12.4", "12.5", "12.6", "12.7", "12.8", "12.9", "13.0", "13.1","13.2","13.3","13.4", "13.5", "13.6", "13.7", "13.8", "13.9", "14.0", "14.1", "14.2", "14.3", "14.4", "14.5", "14.6", "14.7", "14.8", "14.9", "15.0",
                  "15.1", "15.2", "15.3", "15.4", "15.5", "15.6", "15.7", "15.8", "15.9", "16.0", "16.1","16.2","16.3","16.4", "16.5", "16.6", "16.7", "16.8", "16.9", "17.0", "17.1", "17.2", "17.3", "17.4", "17.5", "17.6", "17.7", "17.8", "17.9", "18.0",
                  "18.1", "18.2", "18.3", "18.4", "18.5", "18.6", "18.7", "18.8", "18.9", "19.0", "19.1","19.2","19.3","19.4", "19.5", "19.6", "19.7", "19.8", "19.9", "20.0", "20.1", "20.2", "20.3", "20.4", "20.5", "20.6", "20.7", "20.8", "20.9", "21.0",
                  "21.1","21.2","21.3","21.4", "21.5", "21.6", "21.7", "21.8", "21.9", "22.0", "22.1", "22.2", "22.3", "22.4", "22.5", "22.6", "22.7", "22.8", "22.9", "23.0", "23.1", "23.2", "23.3", "23.4", "23.5", "23.6", "23.7", "23.8", "23.9", "24.0", 
                  "24.1","24.2","24.3","24.4", "24.5", "24.6", "24.7", "24.8", "24.9", "25.0", "25.1", "25.2", "25.3", "25.4", "25.5", "25.6", "25.7", "25.8", "25.9", "26.0", "26.1", "26.2", "26.3", "26.4", "26.5", "26.6", "26.7", "26.8", "26.9", "27.0",
                  "27.1","27.2","27.3","27.4", "27.5", "27.6", "27.7", "27.8", "27.9", "28.0", "28.1", "28.2", "28.3", "28.4", "28.5", "28.6", "28.7", "28.8", "28.9", "29.0", "29.1", "29.2", "29.3", "29.4", "29.5", "29.6", "29.7", "29.8", "29.9", "30.0"]
        AvertBins = np.linspace(0.05, 30.05, 301)
        
        AltLabels2 = ["-19.9", "-19.8", "-19.7", "-19.6", "-19.5", "-19.4", "-19.3", "-19.2", "-19.1", "-19.0", "-18.9","-18.8","-18.7","-18.6", "-18.5", "-18.4", "-18.3", "-18.2", "-18.1", "-18.0", "-17.9", "-17.8","-17.7","-17.6","-17.5","-17.4", "-17.3", "-17.2", "-17.1", "-17.0",
            "-16.9", "-16.8", "-16.7", "-16.6", "-16.5", "-16.4", "-16.3", "-16.2", "-16.1", "-16.0", "-15.9","-15.8","-14.7","-15.6", "-15.5", "-15.4", "-15.3", "-15.2", "-15.2", "-15.0", "-14.9", "-14.8","-14.7","-14.6","-14.5","-14.4", "-14.3", "-14.2", "-14.1", "-14.0",
              "-13.9", "-13.8", "-13.7", "-13.6", "-13.5", "-13.4", "-13.3", "-13.2", "-13.1", "-13.0", "-12.9","-12.8","-12.7","-12.6", "-12.5", "-12.4", "-12.3", "-12.2", "-12.1", "-12.0", "-11.9", "-11.8","-11.7","-11.6","-11.5","-11.4", "-11.3", "-11.2", "-11.1", "-11.0",
              "-10.9", "-10.8", "-10.7", "-10.6", "-10.5", "-10.4", "-10.3", "-10.2", "-10.1", "-10.0", "-9.9","-9.8","-9.7","-9.6", "-9.5", "-9.4", "-9.3", "-9.2", "-9.1", "-9.0", "-8.9", "-8.8","-8.7","-8.6","-8.5","-8.4", "-8.3", "-8.2", "-8.1", "-8.0",
                      "-7.9", "-7.8", "-7.7", "-7.6", "-7.5", "-7.4", "-7.3", "-7.2", "-7.1", "-7.0", "-6.9","-6.8","-6.7","-6.6", "-6.5", "-6.4", "-6.3", "-6.2", "-6.1", "-6.0", "-5.9", "-5.8","-5.7","-5.6","-5.5","-5.4", "-5.3", "-5.2", "-5.1", "-5.0",
                        "-4.9", "-4.8", "-4.7", "-4.6", "-4.5", "-4.4", "-4.3", "-4.2", "-4.1", "-4.0", "-3.9","-3.8","-3.7","-3.6", "-3.5", "-3.4", "-3.3", "-3.2", "-3.1", "-3.0", "-2.9", "-2.8","-2.7","-2.6","-2.5","-2.4", "-2.3", "-2.2", "-2.1", "-2.0",
              "-1.9", "-1.8", "-1.7", "-1.6", "-1.5", "-1.4", "-1.3", "-1.2", "-1.1", "-1.0", "-0.9","-0.8","-0.7","-0.6", "-0.5", "-0.4", "-0.3", "-0.2", "-0.1", "0.0",
              "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0", "1.1","1.2","1.3","1.4", "1.5", "1.6", "1.7", "1.8", "1.9", "2.0", "2.1", "2.2","2.3","2.4","2.5","2.6", "2.7", "2.8", "2.9", "3.0",
                      "3.1", "3.2", "3.3", "3.4", "3.5", "3.6", "3.7", "3.8", "3.9", "4.0", "4.1","4.2","4.3","4.4", "4.5", "4.6", "4.7", "4.8", "4.9", "5.0", "5.1", "5.2","5.3","5.4","5.5","5.6", "5.7", "5.8", "5.9", "6.0",
                        "6.1", "6.2", "6.3", "6.4", "6.5", "6.6", "6.7", "6.8", "6.9", "7.0", "7.1","7.2","7.3","7.4", "7.5", "7.6", "7.7", "7.8", "7.9", "8.0", "8.1", "8.2","8.3","8.4","8.5","8.6", "8.7", "8.8", "8.9", "9.0",
              "9.1", "9.2", "9.3", "9.4", "9.5", "9.6", "9.7", "9.8", "9.9", "10.0", "10.1","10.2","10.3","10.4", "10.5", "10.6", "10.7", "10.8", "10.9", "11.0", "11.1", "11.2", "11.3", "11.4", "11.5", "11.6", "11.7", "11.8", "11.9", "12.0",
              "12.1", "12.2", "12.3", "12.4", "12.5", "12.6", "12.7", "12.8", "12.9", "13.0", "13.1","13.2","13.3","13.4", "13.5", "13.6", "13.7", "13.8", "13.9", "14.0", "14.1", "14.2", "14.3", "14.4", "14.5", "14.6", "14.7", "14.8", "14.9", "15.0",
              "15.1", "15.2", "15.3", "15.4", "15.5", "15.6", "15.7", "15.8", "15.9", "16.0", "16.1","16.2","16.3","16.4", "16.5", "16.6", "16.7", "16.8", "16.9", "17.0", "17.1", "17.2", "17.3", "17.4", "17.5", "17.6", "17.7", "17.8", "17.9", "18.0",
              "18.1", "18.2", "18.3", "18.4", "18.5", "18.6", "18.7", "18.8", "18.9", "19.0", "19.1","19.2","19.3","19.4", "19.5", "19.6", "19.7", "19.8", "19.9", "20.0"]
        AvertBins2 = np.linspace(-19.95, 20.05, 401)
        
        
    AltLabels3 = ["1", "3", "5", "7", "9", "11",  "13", "15", "17", "19", "21", "23", "25",\
                                  "27",  "29"]
    vertBins3 = np.linspace(0., 30., 16)
    
    PLabels = ["15.0", "26.0", "28.0", "32.0", "35.0", "38.0", "42.0", "46.0", "50.0", "56.0", "62.0", "68.0", "74.0", "83.0", "91.0", "100.0", "450.0"]
    vertBins = [10.0, 24.7, 27.3, 30.4, 33.76, 36.1, 39.9, 43.7, 48.3, 53.2, 58.8, 64.6, 71.4, 78.85, 87.15, 95.0, 105.0, 1000.0]
        

    os.chdir(path)
    
    filelist = glob.glob('*.dat')
    #df0 = pd.DataFrame()
    df = pd.DataFrame()

    for file in filelist:
            
        print(file)
            
        folder_list = path.split("/")
        location = folder_list[5]
    #folderlist = glob.glob('*') 
    #filelist = glob.glob('*') 
    #df = pd.DataFrame()

    #loop over different locaiton folders
    # for folder in folderlist:
    #     location = folder
    #     print(location)
    #     os.chdir(path + folder)
    #     filelist = glob.glob('*.dat')


    #     for file in filelist:
        
        #folder_list = path.split("/")
        
        #determine the number of header lines
        count = 0
        FLE_header = 0 #initialize this variable - it will get overwritten
        with open(file) as csv_file:
                csv_reader = csv.reader(csv_file)
                #in the file header find the date (also can find lat lon for a map)
                for row in csv_reader:
                    if count == 1: #specific to FLE file format
                        row_list_ = str(row).split("=")
                        header_str = row_list_[1].split("'")
                        FLE_header = int(header_str[0])
                    elif count == (FLE_header-2):
                            FLE_colnames_i = str(row).split(",")
                    elif count == 4: #specific to FLE file format
                            row_list_ = str(row).split("=")
                            date_str = row_list_[1].split("'")
                            FLE_date_ = (date_str[0])
                            FLE_date_euro = FLE_date_[1:]
                            
                            #Reformat the date from FLE file
                            yr = FLE_date_euro[6:10]
                            mon = FLE_date_euro[3:5]
                            day = FLE_date_euro[0:2]
                            
                            FLE_date = yr + '-' + mon + '-' + day
                            FLE_date_num =  datetime(int(yr),int(mon),int(day))
                    elif count == 8: #latitude
                            row_list_ = str(row).split("=")
                            latitude_str = row_list_[1].split("'")
                            Launch_latitude = float(latitude_str[0])
                                                
                    elif count == 13: #troppause 1
                            row_list_ = str(row).split("=")
                            trop1_str = row_list_[1].split("'")
                            tropopause1 = float(trop1_str [0])
                                    
                    elif count == 14: #tropopause 2
                            row_list_ = str(row).split("=")
                            trop2_str = row_list_[1].split("'")
                            tropopause2 = float(trop2_str[0])

                    count = count+1
                    #dynamic colname list depending on the number of MLS profiles found/ number of columns
                    
        #determine the number of columns/column names! 
        #(& concatenate to units automatically)?
        # remove [] and extra single quotations ' and whitespace from list
        char = '['
        FLE_colnames_i =remove_char(FLE_colnames_i,char)
        char = ']'
        FLE_colnames_i =remove_char(FLE_colnames_i,char)
        char = "'"
        FLE_colnames_i =remove_char(FLE_colnames_i,char)
        char = " "
        FLE_colnames =remove_char(FLE_colnames_i,char)
 

        try:
            df1=pd.read_csv(file, delimiter = ',', header = FLE_header, names = FLE_colnames, dtype=None, low_memory=False)
        except:
            print('Warning: There was an issue reading in this file')
        
            
        df2 = df1.filter(['Press', 'RSPFlag', 'RiseRate', 'GPSalt', 'Theta', 'Tcorr', 'RSTFlag' , 'H2OMr', 'O3Mr','Fl', 'O3Flag', 'RHFP', 'RHIFP', 'TFpHyg',])
        
        df2['Launch'] = FLE_date
        df2['Date'] = FLE_date_num
        df2['Launch_Latitude'] = Launch_latitude
        df2['Location'] = location
        df2['Tropopause1'] = tropopause1
        df2['Tropopause2'] = tropopause2
        df2['Tropopause'] = np.where(df2.Tropopause2 == 99999, df2.Tropopause1,  df2.Tropopause2)
        #df2['Tropopause'] = df2['Tropopause1']
         
        df2['Altitude'] = pd.to_numeric(df2['GPSalt'], errors='coerce').fillna(0) 
        df2['Altitude (km)'] = pd.cut(df2['GPSalt'], bins=AvertBins,labels = AltLabels, include_lowest=True) 
        df2['Altitude (km)'] = pd.to_numeric(df2['Altitude (km)'], errors='coerce').fillna(0)
        df2['Tropopause_Altitude'] = df2['Altitude'] - df2['Tropopause']
 
        df2['Tropopause_Altitude (km)'] = pd.cut(df2['Tropopause_Altitude'], bins=AvertBins2,labels = AltLabels2, include_lowest=True,ordered=False) 
        df2['Tropopause_Altitude (km)'] = pd.to_numeric(df2['Tropopause_Altitude (km)'], errors='coerce').fillna(0)
        
        df2['Air Temperature (K)'] = df1['Tcorr'] + 273.15
        df2['Flight'] = df1['RiseRate'].apply(lambda x: 'Ascent' if x > 0 else 'Descent')
        #df1['Longitude'] = pd.cut(df1['GPS lon [deg]'], bins=lonBins,labels = LonLabels, include_lowest=True) 
        #df1['Longitude'] = pd.to_numeric(df1['Longitude'], errors='coerce').fillna(0) 
        df2['Date'] = pd.to_numeric(df2['Date'], errors='coerce').fillna(0) 
        
        df3 = df2

        del df1, df2
        
        if file == 'la190fle2.dat':
            df3['H2O_ppmv'] = np.NaN
            df3['RH'] = np.NaN
            df3['RHi'] = np.NaN
            df3['TFpHyg']= np.NaN
        else:
            df3['H2O_ppmv'] = df3.apply(lambda x: QAQC(x['H2OMr'], x['Fl']), axis=1)
            df3['RHi'] = df3.apply(lambda x: QAQC(x['RHIFP'], x['Fl']), axis=1)
            df3['RH'] = df3.apply(lambda x: QAQC(x['RHFP'], x['Fl']), axis=1)
            df3['FTFpHyg'] =df3.apply(lambda x: QAQC(x['TFpHyg'], x['Fl']), axis=1)

        df3['O3_ppmv'] = df3.apply(lambda x: QAQC(x['O3Mr'], x['O3Flag']), axis=1)
        df3['Tcorr_degC'] = df3.apply(lambda x: QAQC(x['Tcorr'], x['RSTFlag']), axis=1) #to compare with boulder use Tcorr (C)
        df3['Air Temperature (K)'] = df3.apply(lambda x: QAQC(x['Air Temperature (K)'], x['RSTFlag']), axis=1) #to compare with boulder use Tcorr (C)
        df3['Theta [K]'] = df3.apply(lambda x: QAQC(x['Theta'], x['RSTFlag']), axis=1)
        df3['Pressure_hPa'] = df3.apply(lambda x: QAQC(x['Press'], x['RSPFlag']), axis=1)
        df3['Air Temp (K)'] = df3['Tcorr_degC'] + 273.15 #change to Kelvin
        #df3 = df3[df3['H2O_ppmv'] > 0.0] #this shouldn't be necessary - no actual values 
        
        df = pd.concat([df, df3], ignore_index=True)
        #Calculate median and IQR
    
    df['H2O_q1'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['H2O_ppmv'].transform(lambda x: x.quantile(0.25))
    df['H2O_q2'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['H2O_ppmv'].transform(lambda x: x.quantile(0.5))
    df['H2O_q3'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['H2O_ppmv'].transform(lambda x: x.quantile(0.75))
    
    df['O3_q1'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['O3_ppmv'].transform(lambda x: x.quantile(0.25))
    df['O3_q2'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['O3_ppmv'].transform(lambda x: x.quantile(0.5))
    df['O3_q3'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['O3_ppmv'].transform(lambda x: x.quantile(0.75))
    
    df['T_q1'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Air Temperature (K)'].transform(lambda x: x.quantile(0.25))
    df['T_q2'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Air Temperature (K)'].transform(lambda x: x.quantile(0.5))
    df['T_q3'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Air Temperature (K)'].transform(lambda x: x.quantile(0.75))
    
    df['H2O_q1_TA'] = df.groupby(['Launch', 'Tropopause_Altitude (km)', 'Flight'])['H2O_ppmv'].transform(lambda x: x.quantile(0.25))
    df['H2O_q2_TA'] = df.groupby(['Launch', 'Tropopause_Altitude (km)', 'Flight'])['H2O_ppmv'].transform(lambda x: x.quantile(0.5))
    df['H2O_q3_TA'] = df.groupby(['Launch', 'Tropopause_Altitude (km)', 'Flight'])['H2O_ppmv'].transform(lambda x: x.quantile(0.75))
    
    df['O3_q1_TA'] = df.groupby(['Launch', 'Tropopause_Altitude (km)', 'Flight'])['O3_ppmv'].transform(lambda x: x.quantile(0.25))
    df['O3_q2_TA'] = df.groupby(['Launch', 'Tropopause_Altitude (km)', 'Flight'])['O3_ppmv'].transform(lambda x: x.quantile(0.5))
    df['O3_q3_TA'] = df.groupby(['Launch', 'Tropopause_Altitude (km)', 'Flight'])['O3_ppmv'].transform(lambda x: x.quantile(0.75))
    

    df_test = df.groupby(['Launch','Altitude (km)', 'Flight'],as_index=False).mean() #if smps try mean for all variables except size distribution
    df_test2 = df.groupby(['Launch','Tropopause_Altitude (km)', 'Flight'],as_index=False).mean() #if smps try mean for all variables except size distribution

    df_final = df_test.filter(['Launch','Altitude (km)', 'Flight', 'Pressure_hPa', 'Press [hPa]', 'Theta [K]', 'Tcorr_degC', 'Air Temp (K)', 'O3_ppmv', 'O3Mr', 'O3Flag',
                          'H2O_ppmv', 'H2OMr', 'Fl', 'TFpHyg', 'RH', 'RHi', 'H2O_q1', 'H2O_q2', 'H2O_q3',
                          'O3_q1', 'O3_q2', 'O3_q3', 'T_q1', 'T_q2', 'T_q3', 'H2O_q1_TA', 'H2O_q2_TA', 'H2O_q3_TA',
                          'O3_q1_TA', 'O3_q2_TA', 'O3_q3_TA','Tropopause', 'Tropopause1', 'Tropopause2', 'Tropopause_Altitude (km)'])
    
    df_final2 = df_test2.filter(['Launch','Altitude (km)', 'Flight', 'Pressure_hPa', 'Press [hPa]', 'Theta [K]', 'Tcorr_degC', 'Air Temp (K)', 'O3_ppmv', 'O3Mr', 'O3Flag',
                          'H2O_ppmv', 'H2OMr', 'Fl', 'TFpHyg', 'RH', 'RHi', 'H2O_q1', 'H2O_q2', 'H2O_q3',
                          'O3_q1', 'O3_q2', 'O3_q3', 'T_q1', 'T_q2', 'T_q3', 'H2O_q1_TA', 'H2O_q2_TA', 'H2O_q3_TA',
                          'O3_q1_TA', 'O3_q2_TA', 'O3_q3_TA','Tropopause', 'Tropopause1', 'Tropopause2', 'Tropopause_Altitude (km)'])
 
    
    header_WV = ['\n' + 'Processed data file' + '\n' + 'File Type:' + ',' + '0.1 km altitude binned' + '\n'
              + 'Data Table:' + ',' + 'H2O for Fig.2 and Fig.3' + '\n' + 'Instruments Used:' + ',' + 'Radiosonde' + ',' + 'ECC Ozonesonde,' + 'NOAA FPH,'
              + '\n' + 'Station Name:' +',' + 'Lauder, NZ NIWA' + '\n' + 'Station Latitude:' + ',' + '-45.29' +'\n' + 'Station Longitude:'
              + ','+ '169.67' + '\n' + 'Station Altitude:' + ',' + '370.00 [m]' + '\n'  + 'Pressure/Altitude Source:' 
              + ',' + 'iMet Radiosonde' +'\n' + 'Solar Radiation Correction:' + ',' + 'Yes' +'\n' + 'Variable Descriptions:' + ',' + 'GMT Launch Date [YYYY-mm-DD]' + ',' + 'Binned Altitude [km]'
              + ',' + 'Pressure [mb]' + ',' +'O3 Mixing Ratio [ppmv]' + ',' + 'H2O Mixing Ratio [ppmv]' + ',' + 'Frost Point Hygrometer Temperature [deg C]' + 
              ',' + 'Air Temperature [deg C]' + ',' + 'Potential Temperature [K]' + ',' + 'Tropopause 1 Height [m]' + ',' + 'Tropopause 2 Height [m]' + '\n']
    
    with open('/Users/asher/Documents/PapersInProgress/JGR_ANYSO/Submission/POPSDataForPublication/O3_H2O_table.csv', 'w') as ict:
        for line in header_WV:
            ict.write(line)

    dfWV1 = df_final.filter(['Launch', 'Altitude (km)', 'Pressure_hPa', 'O3_q2', 'H2O_q2', 'TFpHyg', 'Tcorr_degC', 'Theta [K]', 'Tropopause1', 'Tropopause2'])
    dfWV1.to_csv('/Users/asher/Documents/PapersInProgress/JGR_ANYSO/Submission/POPSDataForPublication/O3_H2O_table.csv', mode='a', float_format='%g')


    return df_final, df_final2

def QAQC(Quant, Fl):
    if (Quant == 99999):
        Quant = np.nan
    elif Fl != 1:
        Quant = np.nan
    return Quant


def ChooseRH(Co, RHw, RHi):
    if (Co == 1):
        RH = RHw
    elif (Co == 2):
        RH = RHi
    else:
        RH = np.nan
    return RH 

def StratTrop(PotT, Altitude, Tropopause):
    if (PotT < 380):
        AtmosLayer = 'Troposphere'
    elif (np.isnan(PotT)):
            AtmosLayer = 'Troposphere'
    else:
        AtmosLayer = 'Stratosphere' #maybe should group lowermost strat with lower stratosphere (as discussed in Figure 1)...
        

    return AtmosLayer

def StratTropLidar(Launch, Altitude):
    if Launch == '2019-04-30':
        if Altitude < 14.75:
            AtmosLayer = 'Troposphere'
        else:
            AtmosLayer = 'Stratosphere' 
    elif Launch == '2019-09-03':
        if Altitude < 13.25:
                AtmosLayer = 'Troposphere'
        else:
                AtmosLayer = 'Stratosphere' 
    elif Launch == '2020-01-27':
        if Altitude < 14.25:
                AtmosLayer = 'Troposphere'
        else:
                AtmosLayer = 'Stratosphere' 
    elif Launch == '2020-06-19':
        if Altitude < 13.75:
                AtmosLayer = 'Troposphere'
        else:
                AtmosLayer = 'Stratosphere' 
    elif Launch == '2020-07-23':
        if Altitude < 13.0:
                AtmosLayer = 'Troposphere'
        else:
                AtmosLayer = 'Stratosphere' 
    elif Launch == '2020-11-17':
            if Altitude < 13.25:
                    AtmosLayer = 'Troposphere'
            else:
                    AtmosLayer = 'Stratosphere' 
    elif Launch == '2021-08-31':
                if Altitude < 14.25:
                        AtmosLayer = 'Troposphere'
                else:
                        AtmosLayer = 'Stratosphere' 
    elif Launch == '2021-11-03':
                if Altitude < 14.5:
                        AtmosLayer = 'Troposphere'
                else:
                        AtmosLayer = 'Stratosphere' 
    elif Launch == '2022-01-24':
                if Altitude < 14.5:
                        AtmosLayer = 'Troposphere'
                else:
                        AtmosLayer = 'Stratosphere' 
    else:
        AtmosLayer = 'Troposphere'
    return AtmosLayer

def CARMA_Trop(Launch):
    if Date == '2019-09-03':
            Tropopause = 10.13
    elif Date == '2020-01-27':
            Tropopause = 16.13
    elif Date == '2020-06-19':
            Tropopause = 10.9
    elif Date == '2020-07-23':
                Tropopause = 16.9
    else:
            Tropopause = 9.4
        
    
    return Tropopause #for consistency not using the second WMO tropopause...

def calcTsUnc(Launch, Property, VertIntTot):
    
    #print(Launch)
    #print(Property)
    #print(VertIntTot)
    POPS_Unc = 0
    if (Property == 'sAOD'):
        if (Launch == '2019-04-30'):
            POPS_Unc = VertIntTot * 0.23
        elif (Launch == '2019-09-03'):
                POPS_Unc = VertIntTot * 0.23
        elif (Launch == '2020-01-27'):
                POPS_Unc = VertIntTot * 0.25        
        elif (Launch == '2020-06-19'):
                POPS_Unc = VertIntTot * 0.25  
        elif (Launch == '2020-07-23'):
                POPS_Unc = VertIntTot * 0.29  
        elif (Launch == '2020-11-17'):
                POPS_Unc = VertIntTot * 0.33 
        elif (Launch == '2021-08-31'):
                POPS_Unc = VertIntTot * 0.57
        elif (Launch == '2021-11-03'):
                POPS_Unc = VertIntTot * 0.61
        elif (Launch == '2020-01-24'):
                POPS_Unc = VertIntTot * 0.60  
    else:
        POPS_Unc = VertIntTot * 0.05
    
    return POPS_Unc


def conDS2D(Launch):
    Date = dt.datetime.strptime(Launch, '%Y-%m-%d').toordinal()
    DaysSince2022 = Date - 738156
    return DaysSince2022

def extDY(Launch):
    Date = dt.datetime.strptime(Launch, '%Y-%m-%d')
    Year = Date.year
    return Year

#from netCDF4 import Dataset and create one dataframe seperated by launch...
os.chdir('/Users/asher/Documents/PyroCBPaper/Aerosol_AB')
#os.chdir('/Users/asher/Documents/AnalysesInProgress/MLSdailyV5subset/')
filelist = glob.glob('*')
path = '/Users/asher/Documents/PyroCBPaper/Aerosol_AB'
#path = '/Users/asher/Documents/AnalysesInProgress/MLSdailyV5subset'


for file in filelist:
#for x in range(1): #testing with one file
    #file = 'MLS-Aura_L3DB-H2O_v05-02-c01_2022d045.nc'

    fn = file
    
    if fn == 'B2SAP_Boulder_CO.nc':
        print(file)
        # #netcdf method
        ds = nc.Dataset(fn, "r")
        #print(ds.variables)
        # PV = ds['pv'][:]
        
        Bpops15bins = ds['Dp15'][:]
        #Bpops20bins = ds['Dp20'][:]
        
        Blat_asc = ds['latA'][:]
        Blat_dsc = ds['latD'][:]
        Btime_asc = ds['TimeA'][:] #average UTC time of ascent Unit: seconds after 01/01/2000 midnight UTC
        Btime_dsc = ds['TimeD'][:] #average UTC time of descent Unit: seconds after 01/01/2000 midnight UTC
        Balt = ds['Altitude'][:] #unit is km
        Beff_rad_asc = ds['EffRA'][:] #unit is um
        Beff_rad_dsc = ds['EffRD'][:] #unit is um
        Bsur_area_asc = ds['SACA'][:] #unit is um2 cm3 SACA
        Bpart_conc_asc = ds['PCA'][:] #unit is # cm3 PCA
        BPress = ds['PA'][:] #unit hPa
        BTemp = ds['TA'][:] #unit C
        ds.close()
    
    if fn == 'B2SAP_Lauder_NZ.nc':
    
        print(file)
        # #netcdf method
        ds = nc.Dataset(fn, "r")
        #print(ds.variables)
        # PV = ds['pv'][:]
        
        Lpops15bins = ds['Dp15'][:]
        #Lpops20bins = ds['Dp20'][:]
        LpopsSD_asc = ds['dCdlogDpD'][:]
        Llat_asc = ds['latA'][:]
        Llat_dsc = ds['latD'][:]
        Ltime_asc = ds['TimeA'][:] #average UTC time of ascent Unit: seconds after 01/01/2000 midnight UTC
        Ltime_dsc = ds['TimeD'][:] #average UTC time of descent Unit: seconds after 01/01/2000 midnight UTC
        Lalt = ds['Altitude'][:] #unit is km
        Leff_rad_asc = ds['EffRA'][:] #unit is um
        Leff_rad_dsc = ds['EffRD'][:] #unit is um
        Lsur_area_asc = ds['SACA'][:] #unit is um2 cm3 SACA
        Lpart_conc_asc = ds['PCA'][:] #unit is # cm3 PCA
        LPress = ds['PA'][:] #unit hPa
        LTemp = ds['TA'][:] #unit C
        
        file = '/Users/asher/Documents/NOAAweb/NZ_Paper_IOR1.csv' #'/Users/easher/Documents/NOAAweb/POPS_Sizes.csv'
        df1=pd.read_csv(file,sep=',', dtype=None, engine='python')
        Lpops16binEdges_autoSmoky = np.array( df1[(df1.Composition == 'Smoky') & (df1.Binning == 'Automatic')].Diameter)
        Lpops16binEdges_autoSulfate = np.array( df1[(df1.Composition == 'Sulfate') & (df1.Binning == 'Automatic')].Diameter)
        Lpops16binEdges_Smoky = np.array( df1[(df1.Composition == 'Sulfate') & (df1.Binning == 'Manual')].Diameter)
        
        ds.close()

Btime_asc_int = Btime_asc.astype(int).tolist()
Launches_dt = [dt.datetime(2000,1,1) + dt.timedelta(seconds=each) for each in Btime_asc_int]
Launch = [date_obj.strftime('%Y-%m-%d') for date_obj in Launches_dt]

T_B =  pd.DataFrame(BTemp.T)
T_B.columns = Launch
T_B['Altitude'] = Balt
T_B['Site'] = 'Boulder'
T_B_f = pd.melt(T_B, id_vars=['Altitude'], value_vars= Launch, var_name='Launch', value_name='Temperature_C')

P_B =  pd.DataFrame(BPress.T)
P_B.columns = Launch
P_B['Altitude'] = Balt
P_B['Site'] = 'Boulder'
P_B_f = pd.melt(P_B, id_vars=['Altitude'], value_vars= Launch, var_name='Launch', value_name='Pressure_hPa')


Bpops =  pd.DataFrame(Beff_rad_asc.T)
Bpops.columns = Launch
Bpops['Altitude'] = Balt
Bpops['Site'] = 'Boulder'
Bpops_f = pd.melt(Bpops, id_vars=['Altitude'], value_vars= Launch, var_name='Launch', value_name='EffectiveRadius')
Bpops_f['Days since 2022'] = Bpops_f.apply(lambda x: conDS2D(x['Launch']), axis=1)
Bpops_f['Year'] = Bpops_f.apply(lambda x: extDY(x['Launch']), axis=1)

Bpops_sa =  pd.DataFrame(Bsur_area_asc.T)
Bpops_sa.columns = Launch
Bpops_sa['Altitude'] = Balt
Bpops_sa['Site'] = 'Boulder'
B_sa_pops_f = pd.melt(Bpops_sa, id_vars=['Altitude'], value_vars= Launch, var_name='Launch', value_name='SurfaceArea')
B_sa_pops_f['Days since 2022'] = B_sa_pops_f.apply(lambda x: conDS2D(x['Launch']), axis=1)
B_sa_pops_f['Year'] = B_sa_pops_f.apply(lambda x: extDY(x['Launch']), axis=1)
B_sa_pops_f = pd.merge(B_sa_pops_f,T_B_f, how='left', on=['Launch', 'Altitude'])
B_sa_pops_f = pd.merge(B_sa_pops_f,P_B_f, how='left', on=['Launch', 'Altitude'])


Bpops_pc =  pd.DataFrame(Bpart_conc_asc.T)
Bpops_pc.columns = Launch
Bpops_pc['Altitude'] = Balt
Bpops_pc['Site'] = 'Boulder'
B_pc_pops_f = pd.melt(Bpops_pc, id_vars=['Altitude'], value_vars= Launch, var_name='Launch', value_name='ParticleConcentration')
B_pc_pops_f['Days since 2022'] = B_pc_pops_f.apply(lambda x: conDS2D(x['Launch']), axis=1)
B_pc_pops_f['Year'] = B_pc_pops_f.apply(lambda x: extDY(x['Launch']), axis=1)
B_pc_pops_f = pd.merge(B_pc_pops_f,T_B_f, how='left', on=['Launch', 'Altitude'])
B_pc_pops_f = pd.merge(B_pc_pops_f,P_B_f, how='left', on=['Launch', 'Altitude'])

#only look at 2019 - 2022
Bpops_f = Bpops_f[Bpops_f['Days since 2022'] <= 30]
B_pc_pops_f = B_pc_pops_f[B_pc_pops_f['Days since 2022'] <= 30]
B_sa_pops_f = B_sa_pops_f[B_sa_pops_f['Days since 2022'] <= 30]

#remove data from all three where particle concentration is zero.
IndexRows = B_pc_pops_f.index[B_pc_pops_f['ParticleConcentration'] > 0].tolist()
B_pc_pops_f = B_pc_pops_f.iloc[IndexRows]
B_sa_pops_f = B_sa_pops_f.iloc[IndexRows]
Bpops_f = Bpops_f.iloc[IndexRows]

#define a baseline - removing obviously/freshly pertubed profiles
Bpops_f = Bpops_f[Bpops_f['Launch'] != '2019-08-07']
Bpops_f = Bpops_f[Bpops_f['Launch'] != '2019-08-27']
Bpops_f = Bpops_f[Bpops_f['Launch'] != '2021-04-30']
Bpops_f = Bpops_f[Bpops_f['Launch'] != '2021-05-14']

B_sa_pops_f = B_sa_pops_f[B_sa_pops_f['Launch'] != '2019-08-07']
B_sa_pops_f = B_sa_pops_f[B_sa_pops_f['Launch'] != '2019-08-27']
B_sa_pops_f = B_sa_pops_f[B_sa_pops_f['Launch'] != '2021-04-30']
B_sa_pops_f = B_sa_pops_f[B_sa_pops_f['Launch'] != '2021-05-14']

B_pc_pops_f = B_pc_pops_f[B_pc_pops_f['Launch'] != '2019-08-07']
B_pc_pops_f = B_pc_pops_f[B_pc_pops_f['Launch'] != '2019-08-27']
B_pc_pops_f = B_pc_pops_f[B_pc_pops_f['Launch'] != '2021-04-30']
B_pc_pops_f = B_pc_pops_f[B_pc_pops_f['Launch'] != '2021-05-14']

B_pc_pops_f['ParticleConcentrationAtSTP'] = B_pc_pops_f.apply(lambda x: Calc_STP_AB(x['Temperature_C'], x['Pressure_hPa'], x['ParticleConcentration']), axis=1)
B_sa_pops_f['SurfaceAreaAtSTP'] = B_sa_pops_f.apply(lambda x: Calc_STP_AB(x['Temperature_C'], x['Pressure_hPa'], x['SurfaceArea']), axis=1)

df_Bbase_ab = Bpops_f
df_Bbase_ab['EF_base_q1'] = df_Bbase_ab.groupby(['Altitude'])['EffectiveRadius'].transform(lambda x: x.quantile(0.25))
df_Bbase_ab['EF_base_q2'] = df_Bbase_ab.groupby(['Altitude'])['EffectiveRadius'].transform(lambda x: x.quantile(0.5))
df_Bbase_ab['EF_base_q3'] = df_Bbase_ab.groupby(['Altitude'])['EffectiveRadius'].transform(lambda x: x.quantile(0.75))
df_Bbase_ab = df_Bbase_ab[df_Bbase_ab['Altitude'] >= 1.8]

df_Bbase_pc_ab = B_pc_pops_f
df_Bbase_pc_ab['PC_base_q1'] = df_Bbase_pc_ab.groupby(['Altitude'])['ParticleConcentration'].transform(lambda x: x.quantile(0.25))
df_Bbase_pc_ab['PC_base_q2'] = df_Bbase_pc_ab.groupby(['Altitude'])['ParticleConcentration'].transform(lambda x: x.quantile(0.5))
df_Bbase_pc_ab['PC_base_q3'] = df_Bbase_pc_ab.groupby(['Altitude'])['ParticleConcentration'].transform(lambda x: x.quantile(0.75))

df_Bbase_pc_ab['PC_STP_base_q1'] = df_Bbase_pc_ab.groupby(['Altitude'])['ParticleConcentrationAtSTP'].transform(lambda x: x.quantile(0.25))
df_Bbase_pc_ab['PC_STP_base_q2'] = df_Bbase_pc_ab.groupby(['Altitude'])['ParticleConcentrationAtSTP'].transform(lambda x: x.quantile(0.5))
df_Bbase_pc_ab['PC_STP_base_q3'] = df_Bbase_pc_ab.groupby(['Altitude'])['ParticleConcentrationAtSTP'].transform(lambda x: x.quantile(0.75))

df_Bbase_pc_ab = df_Bbase_pc_ab[df_Bbase_pc_ab['Altitude'] >= 1.8]

df_Bbase_sa_ab = B_sa_pops_f
df_Bbase_sa_ab['SA_base_q1'] = df_Bbase_sa_ab.groupby(['Altitude'])['SurfaceArea'].transform(lambda x: x.quantile(0.25))
df_Bbase_sa_ab['SA_base_q2'] = df_Bbase_sa_ab.groupby(['Altitude'])['SurfaceArea'].transform(lambda x: x.quantile(0.5))
df_Bbase_sa_ab['SA_base_q3'] = df_Bbase_sa_ab.groupby(['Altitude'])['SurfaceArea'].transform(lambda x: x.quantile(0.75))

df_Bbase_sa_ab['SA_STP_base_q1'] = df_Bbase_sa_ab.groupby(['Altitude'])['SurfaceAreaAtSTP'].transform(lambda x: x.quantile(0.25))
df_Bbase_sa_ab['SA_STP_base_q2'] = df_Bbase_sa_ab.groupby(['Altitude'])['SurfaceAreaAtSTP'].transform(lambda x: x.quantile(0.5))
df_Bbase_sa_ab['SA_STP_base_q3'] = df_Bbase_sa_ab.groupby(['Altitude'])['SurfaceAreaAtSTP'].transform(lambda x: x.quantile(0.75))

df_Bbase_sa_ab = df_Bbase_sa_ab[df_Bbase_sa_ab['Altitude'] >= 1.8]

df_Bbase_pc_ab = df_Bbase_pc_ab[df_Bbase_pc_ab['Launch'] == '2019-03-20']
df_Bbase_sa_ab = df_Bbase_sa_ab[df_Bbase_sa_ab['Launch'] == '2019-03-20']
df_Bbase_ab = df_Bbase_ab[df_Bbase_ab['Launch'] == '2019-03-20']


Ltime_asc_int = Ltime_asc.astype(int).tolist()
Launches_dt = [dt.datetime(2000,1,1) + dt.timedelta(seconds=each) for each in Ltime_asc_int]
Launch = [date_obj.strftime('%Y-%m-%d') for date_obj in Launches_dt]

T_L =  pd.DataFrame(LTemp.T)
T_L.columns = Launch
T_L['Altitude'] = Lalt
T_L['Site'] = 'Lauder'
T_L_f = pd.melt(T_L, id_vars=['Altitude'], value_vars= Launch, var_name='Launch', value_name='Temperature_C')

P_L =  pd.DataFrame(LPress.T)
P_L.columns = Launch
P_L['Altitude'] = Lalt
P_L['Site'] = 'Lauder'
P_L_f = pd.melt(P_L, id_vars=['Altitude'], value_vars= Launch, var_name='Launch', value_name='Pressure_hPa')

Lpops =  pd.DataFrame(Leff_rad_asc.T) 
Lpops.columns = Launch
Lpops['Altitude'] = Lalt
Lpops['Site'] = 'Lauder'
Lpops_f = pd.melt(Lpops, id_vars=['Altitude'], value_vars= Launch, var_name='Launch', value_name='EffectiveRadius')
Lpops_f['Days since 2022'] = Lpops_f.apply(lambda x: conDS2D(x['Launch']), axis=1)
Lpops_f['Year'] = Lpops_f.apply(lambda x: extDY(x['Launch']), axis=1)

Lpops_sa =  pd.DataFrame(Lsur_area_asc.T)
Lpops_sa.columns = Launch
Lpops_sa['Altitude'] = Lalt
Lpops_sa['Site'] = 'Lauder'
L_sa_pops_f = pd.melt(Lpops_sa, id_vars=['Altitude'], value_vars= Launch, var_name='Launch', value_name='SurfaceArea')
L_sa_pops_f['Days since 2022'] = L_sa_pops_f.apply(lambda x: conDS2D(x['Launch']), axis=1)
L_sa_pops_f['Year'] = L_sa_pops_f.apply(lambda x: extDY(x['Launch']), axis=1)
L_sa_pops_f = pd.merge(L_sa_pops_f,T_L_f, how='left', on=['Launch', 'Altitude'])
L_sa_pops_f = pd.merge(L_sa_pops_f,P_L_f, how='left', on=['Launch', 'Altitude'])

Lpops_pc =  pd.DataFrame(Lpart_conc_asc.T)
Lpops_pc.columns = Launch
Lpops_pc['Altitude'] = Lalt
Lpops_pc['Site'] = 'Lauder'
L_pc_pops_f = pd.melt(Lpops_pc, id_vars=['Altitude'], value_vars= Launch, var_name='Launch', value_name='ParticleConcentration')
L_pc_pops_f['Days since 2022'] = L_pc_pops_f.apply(lambda x: conDS2D(x['Launch']), axis=1)
L_pc_pops_f['Year'] =L_pc_pops_f.apply(lambda x: extDY(x['Launch']), axis=1)
L_pc_pops_f = pd.merge(L_pc_pops_f,T_L_f, how='left', on=['Launch', 'Altitude'])
L_pc_pops_f = pd.merge(L_pc_pops_f,P_L_f, how='left', on=['Launch', 'Altitude'])

##remove data from all three where particle concentration is zero.
IndexRows = L_pc_pops_f.index[L_pc_pops_f['ParticleConcentration'] > 0].tolist()
L_pc_pops_f = L_pc_pops_f.iloc[IndexRows]
L_sa_pops_f = L_sa_pops_f.iloc[IndexRows]
Lpops_f = Lpops_f.iloc[IndexRows]

Lpops_f = Lpops_f[Lpops_f['Altitude'] >= 0.4]
Lpops_f = Lpops_f[Lpops_f['Launch'] != '2019-04-25']
#Lpops_f = Lpops_f[Lpops_f['Launch'] != '2024-02-22']
#Lpops_f.drop(Lpops_f[(Lpops_f['Altitude'] > 15) & (Lpops_f['Launch'] == '2023-03-02')].index, inplace = True)
#Lpops_f.drop(Lpops_f[(Lpops_f['Altitude'] < 5) & (Lpops_f['Launch'] == '2022-10-11')].index, inplace = True)

#df_Lbase_sa_ab = L_sa_pops_f
L_sa_pops_f = L_sa_pops_f[L_sa_pops_f['Altitude'] >= 0.4]
#L_sa_pops_f = L_sa_pops_f[L_sa_pops_f['SurfaceArea'] != 0]
L_sa_pops_f = L_sa_pops_f[L_sa_pops_f['Launch'] != '2019-04-25']
#L_sa_pops_f = L_sa_pops_f[L_sa_pops_f['Launch'] != '2024-02-22']
#L_sa_pops_f.drop(L_sa_pops_f[(L_sa_pops_f['Altitude'] > 15) & (L_sa_pops_f['Launch'] == '2023-03-02')].index, inplace = True)
#L_sa_pops_f.drop(L_sa_pops_f[(L_sa_pops_f['Altitude'] < 5) & (L_sa_pops_f['Launch'] == '2022-10-11')].index, inplace = True)

L_pc_pops_f = L_pc_pops_f[L_pc_pops_f['Altitude'] >= 0.4]
#L_pc_pops_f = L_pc_pops_f[L_pc_pops_f['ParticleConcentration'] != 0]
L_pc_pops_f = L_pc_pops_f[L_pc_pops_f['Launch'] != '2019-04-25']
#L_pc_pops_f = L_pc_pops_f[L_pc_pops_f['Launch'] != '2024-02-22']
#L_pc_pops_f.drop(L_pc_pops_f[(L_pc_pops_f['Altitude'] > 15) & (L_pc_pops_f['Launch'] == '2023-03-02')].index, inplace = True)
#L_pc_pops_f.drop(L_pc_pops_f[(L_pc_pops_f['Altitude'] < 5) & (L_pc_pops_f['Launch'] == '2022-10-11')].index, inplace = True)
 
#only look at 2019 - 2022
Lpops_f = Lpops_f[Lpops_f['Days since 2022'] <= 30]
L_pc_pops_f = L_pc_pops_f[L_pc_pops_f['Days since 2022'] <= 30]
L_sa_pops_f = L_sa_pops_f[L_sa_pops_f['Days since 2022'] <= 30]

#calculate particle conentration and surface area at STP
L_pc_pops_f['ParticleConcentrationAtSTP'] = L_pc_pops_f.apply(lambda x: Calc_STP_AB(x['Temperature_C'], x['Pressure_hPa'], x['ParticleConcentration']), axis=1)
L_sa_pops_f['SurfaceAreaAtSTP'] = L_sa_pops_f.apply(lambda x: Calc_STP_AB(x['Temperature_C'], x['Pressure_hPa'], x['SurfaceArea']), axis=1)

#define a baseline - removing obviously/freshly pertubed profiles
df_Lbase_ab = Lpops_f
#df_Lbase_ab = df_Lbase_ab[df_Lbase_ab['Launch'] != '2019-09-03']
#df_Lbase_ab = df_Lbase_ab[df_Lbase_ab['Launch'] != '2020-01-27']
df_Lbase_ab = df_Lbase_ab[df_Lbase_ab['Launch'] != '2022-05-25']
df_Lbase_ab['EF_base_q1'] = df_Lbase_ab.groupby(['Altitude'])['EffectiveRadius'].transform(lambda x: x.quantile(0.25))
df_Lbase_ab['EF_base_q2'] = df_Lbase_ab.groupby(['Altitude'])['EffectiveRadius'].transform(lambda x: x.quantile(0.5))
df_Lbase_ab['EF_base_q3'] = df_Lbase_ab.groupby(['Altitude'])['EffectiveRadius'].transform(lambda x: x.quantile(0.75))
df_Lbase_ab = df_Lbase_ab[df_Lbase_ab['Altitude'] > 0.0]

df_Lbase_pc_ab = L_pc_pops_f
#df_Lbase_pc_ab = df_Lbase_pc_ab[df_Lbase_pc_ab['Launch'] != '2019-09-03']
#df_Lbase_pc_ab = df_Lbase_pc_ab[df_Lbase_pc_ab['Launch'] != '2020-01-27']
df_Lbase_pc_ab = df_Lbase_pc_ab[df_Lbase_pc_ab['Launch'] != '2022-05-25']
df_Lbase_pc_ab['PC_base_q1'] = df_Lbase_pc_ab.groupby(['Altitude'])['ParticleConcentration'].transform(lambda x: x.quantile(0.25))
df_Lbase_pc_ab['PC_base_q2'] = df_Lbase_pc_ab.groupby(['Altitude'])['ParticleConcentration'].transform(lambda x: x.quantile(0.5))
df_Lbase_pc_ab['PC_base_q3'] = df_Lbase_pc_ab.groupby(['Altitude'])['ParticleConcentration'].transform(lambda x: x.quantile(0.75))

df_Lbase_pc_ab['PC_STP_base_q1'] = df_Lbase_pc_ab.groupby(['Altitude'])['ParticleConcentrationAtSTP'].transform(lambda x: x.quantile(0.25))
df_Lbase_pc_ab['PC_STP_base_q2'] = df_Lbase_pc_ab.groupby(['Altitude'])['ParticleConcentrationAtSTP'].transform(lambda x: x.quantile(0.5))
df_Lbase_pc_ab['PC_STP_base_q3'] = df_Lbase_pc_ab.groupby(['Altitude'])['ParticleConcentrationAtSTP'].transform(lambda x: x.quantile(0.75))
df_Lbase_pc_ab = df_Lbase_pc_ab[df_Lbase_pc_ab['Altitude'] > 0.0]

df_Lbase_sa_ab = L_sa_pops_f
#df_Lbase_sa_ab = df_Lbase_sa_ab[df_Lbase_sa_ab['Launch'] != '2019-09-03']
#df_Lbase_sa_ab = df_Lbase_sa_ab[df_Lbase_sa_ab['Launch'] != '2020-01-27']
df_Lbase_sa_ab = df_Lbase_sa_ab[df_Lbase_sa_ab['Launch'] != '2022-05-25']
df_Lbase_sa_ab['SA_base_q1'] = df_Lbase_sa_ab.groupby(['Altitude'])['SurfaceArea'].transform(lambda x: x.quantile(0.25))
df_Lbase_sa_ab['SA_base_q2'] = df_Lbase_sa_ab.groupby(['Altitude'])['SurfaceArea'].transform(lambda x: x.quantile(0.5))
df_Lbase_sa_ab['SA_base_q3'] = df_Lbase_sa_ab.groupby(['Altitude'])['SurfaceArea'].transform(lambda x: x.quantile(0.75))

df_Lbase_sa_ab['SA_STP_base_q1'] = df_Lbase_sa_ab.groupby(['Altitude'])['SurfaceAreaAtSTP'].transform(lambda x: x.quantile(0.25))
df_Lbase_sa_ab['SA_STP_base_q2'] = df_Lbase_sa_ab.groupby(['Altitude'])['SurfaceAreaAtSTP'].transform(lambda x: x.quantile(0.5))
df_Lbase_sa_ab['SA_STP_base_q3'] = df_Lbase_sa_ab.groupby(['Altitude'])['SurfaceAreaAtSTP'].transform(lambda x: x.quantile(0.75))

df_Lbase_pc_ab = df_Lbase_pc_ab[df_Lbase_pc_ab['Altitude'] > 0.0]

df_Lbase_pc_ab = df_Lbase_pc_ab[df_Lbase_pc_ab['Launch'] == '2020-11-17']
df_Lbase_sa_ab = df_Lbase_sa_ab[df_Lbase_sa_ab['Launch'] == '2020-11-17']
df_Lbase_ab = df_Lbase_ab[df_Lbase_ab['Launch'] == '2020-11-17']

df_effrad = pd.concat([Bpops_f, Lpops_f], ignore_index=True)




testing = 0
if testing == 0:
    # dp = [100, 200, 300, 400, 500, 1000, 15000] #sp size in nm
    # ndp = [1000.5, 800.3, 600.1, 400, 300, 100, 5] #single particle first?
    
    #create a function to read in the model data, starting with bulk propertiefor s
    CESM_CARMA_dia = np.array([5.00E-07,	6.56E-07, 8.61E-07, 1.13E-06, 1.48E-06, 1.94E-06, 2.55E-06, 3.35E-06, 4.39E-06, 5.76E-06, 7.56E-06, 9.92E-06,
                      1.30E-05, 1.71E-05, 2.24E-05, 2.94E-05, 3.86E-05, 5.06E-05, 6.64E-05, 8.71E-05])*2.0
    
    
    #os.chdir('/Users/asher/Documents/PyroCBPaper/CESM_CARMA/10dayaverages/Bulk/')
    #os.chdir('/Users/asher/Documents/PyroCBPaper/CESM_CARMA/dailyNewest/Bulk/')
    os.chdir('/Users/asher/Documents/PyroCBPaper/CESM_CARMA/neardatesNewest/Bulk/')
    filelist = glob.glob('*.txt')
    #path = '/Users/asher/Documents/PyroCBPaper/CESM_CARMA/Bulk'
    df_CARMA = pd.DataFrame()
            

    
    for file in filelist:
                #print(file)
                
                yr = file[8:12]
                mon = file[12:14]
                day = file[14:16]
                Date = yr + '-' + mon + '-' + day
                
                colNames_CESM_CARMA = ['GeoPotential Height [km]', 'Potential T [K]', 'Aerosol Concentration [#/cm3]', 'Mass Mixing Raio [ug/kg]', 'Extinction [Mm-1]']
    
                df1=pd.read_csv(file, delimiter=r"\s+", names = colNames_CESM_CARMA, dtype=None, low_memory=False)
                df1['Launch'] = Date
                df1['Extinction [km-1]'] = df1['Extinction [Mm-1]'] * 1E-3
                df1.drop(['Mass Mixing Raio [ug/kg]', 'Extinction [Mm-1]'], axis=1, inplace = True)
                df1['diffs_alt'] = df1['GeoPotential Height [km]'].diff() *-1.0
                df1['Level'] = df1.index #LEVEL IS NOT LISTED CORRECTLY (FOR EACH LAUNCH) IT IS MESSING UP THE MERGE
                #df_CARMA = df_CARMA.append(df1)
                df_CARMA = pd.concat([df_CARMA, df1],  ignore_index=True, sort=False) 
                #del df1

    df_CARMA['Number Concentration Column(# m^-2)'] = df_CARMA['Aerosol Concentration [#/cm3]'] * 1E6 * df_CARMA['diffs_alt'] *1E3 # not at STP
    df_CARMA['sAOD'] = df_CARMA['Extinction [km-1]'] * df_CARMA['diffs_alt']
    df_CARMA['Tropopause'] = df_CARMA.apply(lambda x: CARMA_Trop(x['Launch']), axis=1)
    df_CARMA['AtmosLayer'] = df_CARMA.apply(lambda x: StratTrop(x['Potential T [K]'], x['GeoPotential Height [km]'], x['Tropopause']), axis=1)
    df_CARMAvp = df_CARMA #to merge with size distribution data
    
    df_CARMA = df_CARMA[df_CARMA['GeoPotential Height [km]'] > 2.0]
    df_CARMA= df_CARMA[df_CARMA['GeoPotential Height [km]'] <= 28.0]
    df_CARMA.drop(['Potential T [K]', 'GeoPotential Height [km]', 'Tropopause', 'Aerosol Concentration [#/cm3]', 'Extinction [km-1]'], axis = 1, inplace = True)
    

    
    dfc_s = df_CARMA.groupby(['Launch', 'AtmosLayer'], as_index=False).sum()
    dfc_s = dfc_s[dfc_s['AtmosLayer'] != 'Troposphere']  
    Column =['Number Concentration Column(# m^-2)', 'sAOD']
    dfc_int = pd.melt(dfc_s, id_vars=['Launch', 'AtmosLayer'], value_vars= Column, var_name = 'Property', value_name = 'Vertically Integrated Total') 
    dfc_int['Date'] = pd.to_datetime(dfc_int['Launch']).apply(lambda date: date.toordinal())
    dfc_int['Source'] = 'CARMA'
    
    
    
    #pull in FLE2 files & Merge
    path  = '/Users/asher/Documents/PyroCBPaper/Lauder/FLE'
    
    df_w, df_w_TA = FLE(path)
    df_w['Year'] = df_w.apply(lambda x: extDY(x['Launch']), axis=1)
    
    Launches = df_w.Launch.unique().tolist()
    #calculate WMO tropopoause for all launches
    df_trop1 = df_w[df_w['Flight'] == 'Ascent']
    df_trop1.dropna(subset=['Air Temp (K)'], inplace = True)
    df_trop1 = df_trop1.filter(['Launch', 'Altitude (km)', 'Air Temp (K)'])
    
    #defining the WMO tropopause (lapserate > 2 K / km instantaneous and in the next 2 )
    df_trop1['LapseRate_Inst'] = (df_trop1['Air Temp (K)'].shift(-1) - df_trop1['Air Temp (K)'])/0.25   
    df_trop1['LapseRate_2km'] = df_trop1['LapseRate_Inst'].rolling(8).mean().shift(-7)
    
    df_trop1['LRinst_cond_met'] = df_trop1['LapseRate_Inst'].apply(lambda x: 'Y' if x >= -2 else 'N')
    df_trop1['LR2km_cond_met'] = df_trop1['LapseRate_2km'].apply(lambda x: 'Y' if x >= -2 else 'N')
    
    #plot the IOR and an extinction profile calculated for different IORs


df = pd.DataFrame()

#download, plot or write our csv of reer conc for Lauder NZ or Boulder, CO...
choice = 2
if choice == 2:
        Dia = np.array([187.63, 227.39, 275.56, 333.94, 404.69, 490.43, 594.34,
                720.26, 872.86, 1057.79, 1281.89, 1553.48, 1882.61, 2281.47])
        #Dia = np.divide(Dia,1000)
        Rad = np.divide(Dia,2)
        Area = np.pi*4*np.power(Rad,2)
    #lauder
        os.chdir('/Users/asher/Documents/PyroCBPaper/Lauder/CSV/')
        filelist = glob.glob('*.csv')
        path = '/Users/sher/Documents/PyroCBPaper/Lauder/CSV'

        for file in filelist:
            print(file)
            
            colNames = ['GMTdateYmD', 'GMTtime'	, 'milliseconds', 'GMTSSM', 'elapsedMin', 
            'AltitudeIMETkm','frostADC','filteredFrostADC', 'sunlightADC', 'lowSunADC',
            'MirorheatPWM', 'FrostPtTempC','AvgFrostPtTempAVRcalc','H2OmrPPMV', 
            'opticsTempC',	 'opticsHeatADC', 'PressureHygrMB', 'PressSensorTempC', 
            'AltitudeCalchygrometerPressKm',	 'hygrometerBatV',
            'TotColH2Ohygrmm','hygrometerRH%',	 'iMetPressMB', 'iMetTempCorrC',
            'iMetAirTempCraw', 'iMetRH%', 'iMetFrostPtC','iMetinternalTempC',
            'iMetBatV', 'iMetThetaK','iMetPressSensorTempC', 'iMetRHsensorTempC',
            'iMetASCENTms', 'iMetWVPPMV','iMetTotColH2Omm', 'O3mPa', 'O3ppmv',	 'TotColO3DU', 
            'TotColO3DUextrap', 'O3currentUA', 'O3pumpTempC', 'O3pumpIMA', 'O3BatV',
            'GPSlat', 'GPSlon', 'GPSaltkm',	 'GPSsat', 'GPSpressMB', 'GPSwsms',
            'GPSwdDeg',	 'GPSASCENTms',	 'GPSxeastVelms', 'GPSxNorthvelms', 'GPSupvelms',
            'GPSTimeGMT', 'GPSHeadDeg','GPSelevAngleDeg', 'GPSDistkm',	 'predictedlandlat',
            'predictedlandlon',	 'predtimetolandmin','MirrorNum','Cal0degC','Caln45degC',
            'Caln79degC', 'POPSflowCCpS', 'POPStempC',	'POPSbaseline',	 'POPSbaselineStdDev',
            'POPSavgWdithUS', 'B1', 'B2', 'B3', 'B4',
            'B5',	'B6', 'B7',	'B8', 'B9','B10',
            'B11', 'B12', 'B13', 'B14','B15']
            
            
            if file == 'la185_20200127.csv':
                
                colNames = ['GMTdateYmD', 'GMTtime'	, 'milliseconds', 'GMTSSM', 'elapsedMin', 
                    'AltitudeIMETkm','frostADC','filteredFrostADC', 'sunlightADC', 'lowSunADC',
                    'MirorheatPWM', 'FrostPtTempC','AvgFrostPtTempAVRcalc','H2OmrPPMV', 
                    'opticsTempC',	 'opticsHeatADC', 'PressureHygrMB', 'PressSensorTempC', 
                    'AltitudeCalchygrometerPressKm',	 'hygrometerBatV',
                    'TotColH2Ohygrmm','hygrometerRH%',	 'iMetPressMB', 'iMetTempCorrC',
                    'iMetAirTempCraw', 'iMetRH%', 'iMetFrostPtC','iMetinternalTempC',
                    'iMetBatV', 'iMetThetaK','iMetPressSensorTempC', 'iMetRHsensorTempC',
                    'iMetASCENTms', 'iMetWVPPMV','iMetTotColH2Omm', 'O3mPa', 'O3ppmv',	 'TotColO3DU', 
                    'TotColO3DUextrap', 'O3currentUA', 'O3pumpTempC', 'O3pumpIMA', 'O3BatV',
                    'GPSlat', 'GPSlon', 'GPSaltkm',	 'GPSsat', 'GPSpressMB', 'GPSwsms',
                    'GPSwdDeg',	 'GPSASCENTms',	 'GPSxeastVelms', 'GPSxNorthvelms', 'GPSupvelms',
                    'GPSTimeGMT', 'GPSHeadDeg','GPSelevAngleDeg', 'GPSDistkm',	 'predictedlandlat',
                    'predictedlandlon',	 'predtimetolandmin','MirrorNum','Cal0degC','Caln45degC',
                    'Caln79degC', 'POPSflowCCpS', 'POPStempC',	'POPS Data Status', 'POPSbaseline',	 'POPSbaselineStdDev',
                    'POPSavgWdithUS', 'B1', 'B2', 'B3', 'B4',
                    'B5',	'B6', 'B7',	'B8',	'B9','B10',
                    'B11', 'B12', 'B13', 'B14','B15']
            
            if file == 'la190_20200619.csv':
                
                colNames = ['GMTdateYmD', 'GMTtime'	, 'milliseconds', 'GMTSSM', 'elapsedMin', 
                            'AltitudeIMETkm','iMetPressMB', 'iMetTempCorrC',
                            'iMetAirTempCraw', 'iMetRH%', 'iMetFrostPtC','iMetinternalTempC',
                            'iMetBatV', 'iMetThetaK','iMetPressSensorTempC', 'iMetRHsensorTempC',
                            'iMetASCENTms', 'iMetWVPPMV','iMetTotColH2Omm', 'O3mPa', 'O3ppmv', 'TotColO3DU', 
                            'TotColO3DUextrap', 'O3currentUA', 'O3pumpTempC', 'O3pumpIMA', 'O3BatV',
                            'GPSlat', 'GPSlon', 'GPSaltkm','GPSsat', 'GPSpressMB', 'GPSwsms',
                            'GPSwdDeg',	 'GPSASCENTms',	 'GPSxeastVelms', 'GPSxNorthvelms', 'GPSupvelms',
                            'GPSTimeGMT', 'GPSHeadDeg','GPSelevAngleDeg', 'GPSDistkm',	 'predictedlandlat',
                            'predictedlandlon',	 'predtimetolandmin','POPSflowCCpS', 'POPStempC','POPSmaxStd',
                            'PSbaseline',	 'POPSbaselineStdDev',
                            'POPSavgWdithUS', 'B1', 'B2', 'B3', 'B4','B5',	'B6', 'B7',	'B8',	'B9','B10',
                            'B11', 'B12', 'B13', 'B14','B15']
                    
            if file == 'la191_20200723.csv':
                
                    colNames = ['GMTdateYmD', 'GMTtime'	, 'milliseconds', 'GMTSSM', 'elapsedMin', 
            'AltitudeIMETkm','frostADC','filteredFrostADC', 'sunlightADC', 'lowSunADC',
            'MirorheatPWM', 'FrostPtTempC','AvgFrostPtTempAVRcalc','H2OmrPPMV', 
            'opticsTempC',	 'opticsHeatADC', 'PressureHygrMB', 'PressSensorTempC', 
            'AltitudeCalchygrometerPressKm',	 'hygrometerBatV',
            'TotColH2Ohygrmm','hygrometerRH%',	 'iMetPressMB', 'iMetTempCorrC',
            'iMetAirTempCraw', 'iMetRH%', 'iMetFrostPtC','iMetinternalTempC',
            'iMetBatV', 'iMetThetaK','iMetPressSensorTempC', 'iMetRHsensorTempC',
            'iMetASCENTms', 'iMetWVPPMV','iMetTotColH2Omm', 'O3mPa', 'O3ppmv',	 'TotColO3DU', 
            'TotColO3DUextrap', 'O3currentUA', 'O3pumpTempC', 'O3pumpIMA', 'O3BatV',
            'GPSlat', 'GPSlon', 'GPSaltkm',	 'GPSsat', 'GPSpressMB', 'GPSwsms',
            'GPSwdDeg',	 'GPSASCENTms',	 'GPSxeastVelms', 'GPSxNorthvelms', 'GPSupvelms',
            'GPSTimeGMT', 'GPSHeadDeg','GPSelevAngleDeg', 'GPSDistkm',	 'predictedlandlat',
            'predictedlandlon',	 'predtimetolandmin','MirrorNum','Cal0degC','Caln45degC',
            'Caln79degC', 'POPSflowCCpS', 'POPStempC',	'POPSMaxStdDev', 'POPSbaseline',	 'POPSbaselineStdDev',
            'POPSavgWdithUS', 'B1', 'B2', 'B3', 'B4',
            'B5',	'B6', 'B7',	'B8',	'B9','B10',
            'B11', 'B12', 'B13', 'B14','B15']
                    
            if file == 'la194_20201117.csv':
                
                    colNames = ['GMTdateYmD', 'GMTtime'	, 'milliseconds', 'GMTSSM', 'elapsedMin', 
            'AltitudeIMETkm','frostADC','filteredFrostADC', 'sunlightADC', 'lowSunADC',
            'MirorheatPWM', 'FrostPtTempC','AvgFrostPtTempAVRcalc','H2OmrPPMV', 
            'opticsTempC',	 'opticsHeatADC', 'PressureHygrMB', 'PressSensorTempC', 
            'AltitudeCalchygrometerPressKm',	 'hygrometerBatV',
            'TotColH2Ohygrmm','hygrometerRH%',	 'iMetPressMB', 'iMetTempCorrC',
            'iMetAirTempCraw', 'iMetRH%', 'iMetFrostPtC','iMetinternalTempC',
            'iMetBatV', 'iMetThetaK','iMetPressSensorTempC', 'iMetRHsensorTempC',
            'iMetASCENTms', 'iMetWVPPMV','iMetTotColH2Omm', 'O3mPa', 'O3ppmv',	 'TotColO3DU', 
            'TotColO3DUextrap', 'O3currentUA', 'O3pumpTempC', 'O3pumpIMA', 'O3BatV',
            'GPSlat', 'GPSlon', 'GPSaltkm',	 'GPSsat', 'GPSpressMB', 'GPSwsms',
            'GPSwdDeg',	 'GPSASCENTms',	 'GPSxeastVelms', 'GPSxNorthvelms', 'GPSupvelms',
            'GPSTimeGMT', 'GPSHeadDeg','GPSelevAngleDeg', 'GPSDistkm',	 'predictedlandlat',
            'predictedlandlon',	 'predtimetolandmin','MirrorNum','Cal0degC','Caln45degC',
            'Caln79degC', 'POPSflowCCpS', 'POPStempC', 'POPSMaxStdDev',	'POPSbaseline',	 'POPSbaselineStdDev',
            'POPSavgWdithUS', 'B1', 'B2', 'B3', 'B4',
            'B5',	'B6', 'B7',	'B8',	'B9','B10',
            'B11', 'B12', 'B13', 'B14','B15']
                    
            if file == 'la203_20210831.csv':
                
                    colNames = ['GMTdateYmD', 'GMTtime'	, 'milliseconds', 'GMTSSM', 'elapsedMin', 
            'AltitudeIMETkm','frostADC','filteredFrostADC', 'sunlightADC', 'lowSunADC',
            'MirorheatPWM', 'FrostPtTempC','AvgFrostPtTempAVRcalc','H2OmrPPMV', 
            'opticsTempC',	 'opticsHeatADC', 'PressureHygrMB', 'PressSensorTempC', 
            'AltitudeCalchygrometerPressKm',	 'hygrometerBatV',
            'TotColH2Ohygrmm','hygrometerRH%',	 'iMetPressMB', 'iMetTempCorrC',
            'iMetAirTempCraw', 'iMetRH%', 'iMetFrostPtC','iMetinternalTempC',
            'iMetBatV', 'iMetThetaK','iMetPressSensorTempC', 'iMetRHsensorTempC',
            'iMetASCENTms', 'iMetWVPPMV','iMetTotColH2Omm', 'O3mPa', 'O3ppmv',	 'TotColO3DU', 
            'TotColO3DUextrap', 'O3currentUA', 'O3pumpTempC', 'O3pumpIMA', 'O3BatV',
            'GPSlat', 'GPSlon', 'GPSaltkm',	 'GPSsat', 'GPSpressMB', 'GPSwsms',
            'GPSwdDeg',	 'GPSASCENTms',	 'GPSxeastVelms', 'GPSxNorthvelms', 'GPSupvelms',
            'GPSTimeGMT', 'GPSHeadDeg','GPSelevAngleDeg', 'GPSDistkm',	 'predictedlandlat',
            'predictedlandlon',	 'predtimetolandmin','MirrorNum','Cal0degC','Caln45degC',
            'Caln79degC', 'POPSflowCCpS', 'POPStempC', 'POPSMaxStdDev',	'POPSbaseline',	 'POPSbaselineStdDev',
            'POPSavgWdithUS', 'B1', 'B2', 'B3', 'B4',
            'B5',	'B6', 'B7',	'B8',	'B9','B10',
            'B11', 'B12', 'B13', 'B14','B15']
                
            if file == 'la205_20211103.csv':
                
                colNames = ['GMTdateYmD', 'GMTtime'	, 'milliseconds', 'GMTSSM', 'elapsedMin', 
        'AltitudeIMETkm','frostADC','filteredFrostADC', 'sunlightADC', 'lowSunADC',
        'MirorheatPWM', 'FrostPtTempC','AvgFrostPtTempAVRcalc','H2OmrPPMV', 
        'opticsTempC',	 'opticsHeatADC', 'PressureHygrMB', 'PressSensorTempC', 
        'AltitudeCalchygrometerPressKm',	 'hygrometerBatV',
        'TotColH2Ohygrmm','hygrometerRH%',	 'iMetPressMB', 'iMetTempCorrC',
        'iMetAirTempCraw', 'iMetRH%', 'iMetFrostPtC','iMetinternalTempC',
        'iMetBatV', 'iMetThetaK','iMetPressSensorTempC', 'iMetRHsensorTempC',
        'iMetASCENTms', 'iMetWVPPMV','iMetTotColH2Omm', 'O3mPa', 'O3ppmv',	 'TotColO3DU', 
        'TotColO3DUextrap', 'O3currentUA', 'O3pumpTempC', 'O3pumpIMA', 'O3BatV',
        'GPSlat', 'GPSlon', 'GPSaltkm',	 'GPSsat', 'GPSpressMB', 'GPSwsms',
        'GPSwdDeg',	 'GPSASCENTms',	 'GPSxeastVelms', 'GPSxNorthvelms', 'GPSupvelms',
        'GPSTimeGMT', 'GPSHeadDeg','GPSelevAngleDeg', 'GPSDistkm',	 'predictedlandlat',
        'predictedlandlon',	 'predtimetolandmin','MirrorNum','Cal0degC','Caln45degC',
        'Caln79degC', 'POPSflowCCpS', 'POPStempC', 'POPSMaxStdDev',	'POPSbaseline',	 'POPSbaselineStdDev',
        'POPSavgWdithUS', 'B1', 'B2', 'B3', 'B4',
        'B5',	'B6', 'B7',	'B8',	'B9','B10',
        'B11', 'B12', 'B13', 'B14','B15']
            
            
            if file == 'la206_20220124.csv':
                
                colNames = ['GMTdateYmD', 'GMTtime'	, 'milliseconds', 'GMTSSM', 'elapsedMin', 
        'AltitudeIMETkm','frostADC','filteredFrostADC', 'sunlightADC', 'lowSunADC',
        'MirorheatPWM', 'FrostPtTempC','AvgFrostPtTempAVRcalc','H2OmrPPMV', 
        'opticsTempC',	 'opticsHeatADC', 'PressureHygrMB', 'PressSensorTempC', 
        'AltitudeCalchygrometerPressKm',	 'hygrometerBatV',
        'TotColH2Ohygrmm','hygrometerRH%',	 'iMetPressMB', 'iMetTempCorrC',
        'iMetAirTempCraw', 'iMetRH%', 'iMetFrostPtC','iMetinternalTempC',
        'iMetBatV', 'iMetThetaK','iMetPressSensorTempC', 'iMetRHsensorTempC',
        'iMetASCENTms', 'iMetWVPPMV','iMetTotColH2Omm', 'O3mPa', 'O3ppmv',	 'TotColO3DU', 
        'TotColO3DUextrap', 'O3currentUA', 'O3pumpTempC', 'O3pumpIMA', 'O3BatV',
        'GPSlat', 'GPSlon', 'GPSaltkm',	 'GPSsat', 'GPSpressMB', 'GPSwsms',
        'GPSwdDeg',	 'GPSASCENTms',	 'GPSxeastVelms', 'GPSxNorthvelms', 'GPSupvelms',
        'GPSTimeGMT', 'GPSHeadDeg','GPSelevAngleDeg', 'GPSDistkm',	 'predictedlandlat',
        'predictedlandlon',	 'predtimetolandmin','MirrorNum','Cal0degC','Caln45degC',
        'Caln79degC', 'POPSflowCCpS', 'POPStempC', 'POPSMaxStdDev',	'POPSbaseline',	 'POPSbaselineStdDev',
        'POPSavgWdithUS', 'B1', 'B2', 'B3', 'B4',
        'B5',	'B6', 'B7',	'B8',	'B9','B10',
        'B11', 'B12', 'B13', 'B14','B15']
            # if file == 'bu811_20200124.csv':
            #            colNames = ['GMTdateYmD', 'GMTtime'	, 'milliseconds', 'GMTSSM', 'elapsedMin', 
            # 'AltitudeIMETkm','frostADC','filteredFrostADC', 'sunlightADC', 'lowSunADC',
            # 'MirorheatPWM', 'FrostPtTempC','AvgFrostPtTempAVRcalc','H2OmrPPMV', 
            # 'opticsTempC',	 'opticsHeatADC', 'PressureHygrMB', 'PressSensorTempC', 
            # 'AltitudeCalchygrometerPressKm',	 'hygrometerBatV',
            # 'TotColH2Ohygrmm','hygrometerRH%',	 'iMetPressMB', 'iMetTempCorrC',
            # 'iMetAirTempCraw', 'iMetRH%', 'iMetFrostPtC','iMetinternalTempC',
            # 'iMetBatV', 'iMetThetaK','iMetPressSensorTempC', 'iMetRHsensorTempC',
            # 'iMetASCENTms', 'iMetWVPPMV','iMetTotColH2Omm', 'O3mPa', 'O3ppmv',	 'TotColO3DU', 
            # 'TotColO3DUextrap', 'O3currentUA', 'O3pumpTempC', 'O3pumpIMA', 'O3BatV',
            # 'GPSlat', 'GPSlon', 'GPSaltkm',	 'GPSsat', 'GPSpressMB', 'GPSwsms',
            # 'GPSwdDeg',	 'GPSASCENTms',	 'GPSxeastVelms', 'GPSxNorthvelms', 'GPSupvelms',
            # 'GPSTimeGMT', 'GPSHeadDeg','GPSelevAngleDeg', 'GPSDistkm',	 'predictedlandlat',
            # 'predictedlandlon',	 'predtimetolandmin','MirrorNum','Cal0degC','Caln45degC',
            # 'Caln79degC', 'POPSflowCCpS', 'POPStempC',	'MaxSTD', 'POPSbaseline',	 'POPSbaselineStdDev',
            # 'POPSavgWdithUS', 'B1', 'B2', 'B3', 'B4',
            # 'B5',	'B6', 'B7',	'B8',	'B9','B10',
            # 'B11', 'B12', 'B13', 'B14','B15']

            df1=pd.read_csv(file, delimiter = ',', header=26,names = colNames, dtype=None, low_memory=False)
            
            #Filter out data that you aren't using for plotting etc        

            df1 = df1.filter(['GMTdateYmD','iMetTempCorrC','iMetPressMB',
                    'iMetASCENTms', 'GPSlat', 'GPSlon','GPSaltkm', 'POPSflowCCpS', 
                    'POPStempC', 'POPSbaselineStdDev','POPSavgWdithUS', 
                    'B1', 'B2', 'B3','B4','B5','B6','B7','B8','B9','B10','B11', 'B12', 
                    'B13', 'B14','B15'], axis=1)
            #df = df.append(df1)
            df = pd.concat([df,df1], ignore_index=True)
            del df1

            
df['POPSflowCCpS'] = pd.to_numeric(df['POPSflowCCpS'], errors='coerce').fillna(0)            
#keep only ascent data when POPS is turned on 
df = df[df.POPSflowCCpS >= 3.0]
df = df[df.POPSavgWdithUS > 0]
df = df[df.B1 != 99999]
df = df[df.GPSaltkm != 99999]
df = df[df.POPSbaselineStdDev <= 13]
#df = df[df.GMTdateYmD != '1/28/20']
df = df[df.GMTdateYmD != '2019-04-25']

df["GMTdateYmD"]= df["GMTdateYmD"].str.replace("1/27/20", "2020-01-27", case = False) 
df["GMTdateYmD"]= df["GMTdateYmD"].str.replace("1/28/20", "2020-01-28", case = False) 
df["GMTdateYmD"]= df["GMTdateYmD"].str.replace("6/19/20", "2020-06-19", case = False) 
df["GMTdateYmD"]= df["GMTdateYmD"].str.replace("11/3/21", "2021-11-03", case = False) 
df["GMTdateYmD"]= df["GMTdateYmD"].str.replace("11/4/21", "2021-11-04", case = False) 
#df["GMTdateYmD"]= df["GMTdateYmD"].str.replace("7/23/23", "2020-07-23", case = False) 
#df["GMTdateYmD"]= df["GMTdateYmD"].str.replace("1/24/20", "2020-01-24", case = False) 
df["GMTdateYmD"]= df["GMTdateYmD"].str.replace("6/7/19", "2019-06-19", case = False) 
#df["GMTdateYmD"]= df["GMTdateYmD"].str.replace("3/20/19", "2019-03-20", case = False) 
#df["GMTdateYmD"]= df["GMTdateYmD"].str.replace("11/17/20", "2020-11-17", case = False) 
#df["GMTdateYmD"]= df["GMTdateYmD"].str.replace("11/18/20", "2020-11-18", case = False) 

#create a launch date so that if a flight covers two days, all gets included as one "launch"
conditions = [
    (df['GMTdateYmD'] == '2019-04-30') ,
    (df['GMTdateYmD'] == '2019-05-01') ,
    (df['GMTdateYmD'] == '2019-09-03') ,
    (df['GMTdateYmD'] == '2019-09-04') ,
    (df['GMTdateYmD'] == '2020-01-27') ,
    (df['GMTdateYmD'] == '2020-01-28') ,
    (df['GMTdateYmD'] == '2020-06-19') ,
    (df['GMTdateYmD'] == '2020-06-20') ,
    (df['GMTdateYmD'] == '2020-07-23') ,
    (df['GMTdateYmD'] == '2020-07-24') ,
    (df['GMTdateYmD'] == '2020-11-17') ,
    (df['GMTdateYmD'] == '2020-11-18') ,
    (df['GMTdateYmD'] == '2021-08-31') ,
    (df['GMTdateYmD'] == '2021-09-01') ,
    (df['GMTdateYmD'] == '2021-11-03') ,
    (df['GMTdateYmD'] == '2021-11-04') ,
    (df['GMTdateYmD'] == '2022-01-24') ,
    (df['GMTdateYmD'] == '2022-01-25') 
    ]

choices = ['2019-04-30', '2019-04-30', '2019-09-03', '2019-09-03', '2020-01-27', '2020-01-27', '2020-06-19', '2020-06-19', '2020-07-23','2020-07-23', '2020-11-17', '2020-11-17', '2021-08-31', '2021-08-31','2021-11-03', '2021-11-03', '2022-01-24', '2022-01-24']

df['Launch']=np.select(conditions, choices, default = 'Unassigned')

#seperate out ascent and descent

dfa = df[df.iMetASCENTms >= 0]
dfa['Flight'] = 'Ascent'
# else:
dfb = df[df.iMetASCENTms <= 0]
dfb['Flight']= 'Descent'

df = pd.concat([dfa,dfb], axis=0, join='outer', ignore_index=False)


# AltLabels = ["0.25", "0.5","0.75", "1", "1.25", "1.5", "1.75", "2", "2.25", "2.5", "2.75", "3", "3.25", "3.5", "3.75", "4","4.25", "4.5", "4.75", "5", "5.25", "5.5", "5.75", "6",\
#               "6.25", "6.5", "6.75", "7", "7.25", "7.5", "7.75", "8", "8.25", "8.5", "8.75", "9", "9.25", "9.5", "9.75", "10",\
#               "10.25", "10.5", "10.75", "11", "11.25", "11.5", "11.75", "12", "12.25", "12.5", "12.75", "13", "13.25", "13.5", "13.75", "14",\
#                   "14.25", "14.5", "14.75", "15", "15.25", "15.5", "15.75", "16", "16.25", "16.5", "16.75", "17", "17.25", "17.5", "17.75", "18",\
#                       "18.25", "18.5", "18.75", "19", "19.25", "19.5", "19.75", "20", "20.25", "20.5", "20.75", "21", "21.25", "21.5", "21.75", "22",\
#                           "22.25", "22.5", "22.75", "23", "23.25", "23.5", "23.75", "24", "24.25", "24.5", "24.75", "25", "25.25", "25.5", "25.75", "26",\
#                               "26.25", "26.5", "26.75", "27", "27.25", "27.5", "27.75", "28", "28.25", "28.5", "28.75", "29",\
#                       "29.25", "29.5", "29.75", "30.00"]
# AvertBins = np.linspace(0.125, 30.125, 121)

VertRes = 100 #250 m bins like lev files   
if VertRes == 250: #data binned to to every 250 m
    AltLabels = ["0.25", "0.5","0.75", "1", "1.25", "1.5", "1.75", "2", "2.25", "2.5", "2.75", "3", "3.25", "3.5", "3.75", "4","4.25", "4.5", "4.75", "5", "5.25", "5.5", "5.75", "6",\
              "6.25", "6.5", "6.75", "7", "7.25", "7.5", "7.75", "8", "8.25", "8.5", "8.75", "9", "9.25", "9.5", "9.75", "10",\
              "10.25", "10.5", "10.75", "11", "11.25", "11.5", "11.75", "12", "12.25", "12.5", "12.75", "13", "13.25", "13.5", "13.75", "14",\
                  "14.25", "14.5", "14.75", "15", "15.25", "15.5", "15.75", "16", "16.25", "16.5", "16.75", "17", "17.25", "17.5", "17.75", "18",\
                      "18.25", "18.5", "18.75", "19", "19.25", "19.5", "19.75", "20", "20.25", "20.5", "20.75", "21", "21.25", "21.5", "21.75", "22",\
                          "22.25", "22.5", "22.75", "23", "23.25", "23.5", "23.75", "24", "24.25", "24.5", "24.75", "25", "25.25", "25.5", "25.75", "26",\
                              "26.25", "26.5", "26.75", "27", "27.25", "27.5", "27.75", "28", "28.25", "28.5", "28.75", "29",\
                      "29.25", "29.5", "29.75", "30.00"]
    AvertBins = np.linspace(0.125, 30.125, 121)
    
elif VertRes == 100: #data binned to to every 100 m
    AltLabels = ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0", "1.1","1.2","1.3","1.4", "1.5", "1.6", "1.7", "1.8", "1.9", "2.0", "2.1", "2.2","2.3","2.4","2.5","2.6", "2.7", "2.8", "2.9", "3.0",
            "3.1", "3.2", "3.3", "3.4", "3.5", "3.6", "3.7", "3.8", "3.9", "4.0", "4.1","4.2","4.3","4.4", "4.5", "4.6", "4.7", "4.8", "4.9", "5.0", "5.1", "5.2","5.3","5.4","5.5","5.6", "5.7", "5.8", "5.9", "6.0",
              "6.1", "6.2", "6.3", "6.4", "6.5", "6.6", "6.7", "6.8", "6.9", "7.0", "7.1","7.2","7.3","7.4", "7.5", "7.6", "7.7", "7.8", "7.9", "8.0", "8.1", "8.2","8.3","8.4","8.5","8.6", "8.7", "8.8", "8.9", "9.0",
              "9.1", "9.2", "9.3", "9.4", "9.5", "9.6", "9.7", "9.8", "9.9", "10.0", "10.1","10.2","10.3","10.4", "10.5", "10.6", "10.7", "10.8", "10.9", "11.0", "11.1", "11.2", "11.3", "11.4", "11.5", "11.6", "11.7", "11.8", "11.9", "12.0",
              "12.1", "12.2", "12.3", "12.4", "12.5", "12.6", "12.7", "12.8", "12.9", "13.0", "13.1","13.2","13.3","13.4", "13.5", "13.6", "13.7", "13.8", "13.9", "14.0", "14.1", "14.2", "14.3", "14.4", "14.5", "14.6", "14.7", "14.8", "14.9", "15.0",
              "15.1", "15.2", "15.3", "15.4", "15.5", "15.6", "15.7", "15.8", "15.9", "16.0", "16.1","16.2","16.3","16.4", "16.5", "16.6", "16.7", "16.8", "16.9", "17.0", "17.1", "17.2", "17.3", "17.4", "17.5", "17.6", "17.7", "17.8", "17.9", "18.0",
              "18.1", "18.2", "18.3", "18.4", "18.5", "18.6", "18.7", "18.8", "18.9", "19.0", "19.1","19.2","19.3","19.4", "19.5", "19.6", "19.7", "19.8", "19.9", "20.0", "20.1", "20.2", "20.3", "20.4", "20.5", "20.6", "20.7", "20.8", "20.9", "21.0",
              "21.1","21.2","21.3","21.4", "21.5", "21.6", "21.7", "21.8", "21.9", "22.0", "22.1", "22.2", "22.3", "22.4", "22.5", "22.6", "22.7", "22.8", "22.9", "23.0", "23.1", "23.2", "23.3", "23.4", "23.5", "23.6", "23.7", "23.8", "23.9", "24.0", 
              "24.1","24.2","24.3","24.4", "24.5", "24.6", "24.7", "24.8", "24.9", "25.0", "25.1", "25.2", "25.3", "25.4", "25.5", "25.6", "25.7", "25.8", "25.9", "26.0", "26.1", "26.2", "26.3", "26.4", "26.5", "26.6", "26.7", "26.8", "26.9", "27.0",
              "27.1","27.2","27.3","27.4", "27.5", "27.6", "27.7", "27.8", "27.9", "28.0", "28.1", "28.2", "28.3", "28.4", "28.5", "28.6", "28.7", "28.8", "28.9", "29.0", "29.1", "29.2", "29.3", "29.4", "29.5", "29.6", "29.7", "29.8", "29.9", "30.0"]
    AvertBins = np.linspace(0.05, 30.05, 301)
        

        
#bins=np.linspace(0, 30, 61)
# AltLabels = ["0", "0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5","5.5","6","6.5","7","7.5","8","8.5","9","9.5","10",\
#               "10.5","11","11.5","12","12.5","13","13.5","14","14.5","15","15.5","16","16.5","17","17.5","18","18.5","19","19.5",\
#                   "20","20.5", "21","21.5","22","22.5","23","23.5","24","24.5","25","25.5","26","26.5","27","27.5","28","28.5","29",\
#                       "29.5"]

#create more altitude bins
df['Altitude (km)'] = pd.cut(df['GPSaltkm'], bins = AvertBins,labels = AltLabels, include_lowest=True) # bin every 100 m;  if 61 bin every 500m 
df['Altitude (km)'] = pd.to_numeric(df['Altitude (km)'], errors='coerce').fillna(0) 

#flow corrected for differences in POPS temperature and outside air temperature (in Kelvin)
df['POPSTempCorrection'] = (df['iMetTempCorrC']+273.15)/(df['POPStempC']+273.15)
df['Air Temp K'] = df['iMetTempCorrC']+273.15
df['POPSflowCorrCCpS'] = df['POPSflowCCpS']*df['POPSTempCorrection']

#divide counts in all the bins by flow...
df[['B1', 'B2','B3','B4', 'B5', 'B6', 'B7','B8', 'B9', 'B10', 'B11', 'B12', 'B13', 'B14','B15']]\
= df[['B1', 'B2','B3','B4', 'B5', 'B6', 'B7','B8', 'B9', 'B10', 'B11', 'B12', 'B13','B14','B15']].div(df['POPSflowCorrCCpS'].values,axis=0)

df['Composition'] = df.apply(lambda x: Plume(x['Launch'], x['Flight'], x['GPSaltkm']), axis=1)
choice1 = 1
choice2 = 2
df['Aer_Conc'] = df.apply(lambda x: CalcAerConc(choice1, x['Air Temp K'], x['iMetPressMB'], x['B1'], x['B2'], x['B3'], x['B4'], x['B5'], x['B6'], x['B7'], x['B8'], x['B9'], x['B10'], x['B11'], x['B12'], x['B13'], x['B14'], x['B15']), axis=1)
df['Aer_ConcAtSTP'] = df.apply(lambda x: Calc_STP_AB(x['iMetTempCorrC'], x['iMetPressMB'], x['Aer_Conc']), axis=1)
df['Aer_Conc_205nm'] = df.apply(lambda x: CalcAerConc( choice2, x['Air Temp K'], x['iMetPressMB'], x['B1'], x['B2'], x['B3'], x['B4'], x['B5'], x['B6'], x['B7'], x['B8'], x['B9'], x['B10'], x['B11'], x['B12'], x['B13'], x['B14'], x['B15']), axis=1)
df['Surf_Area'] = df.apply(lambda x: CalcDrySA(x['Air Temp K'], x['iMetPressMB'], x['Composition'], x['B1'], x['B2'], x['B3'], x['B4'], x['B5'], x['B6'], x['B7'], x['B8'], x['B9'], x['B10'], x['B11'], x['B12'], x['B13'], x['B14'], x['B15']), axis=1)
df['Surf_AreaAtSTP'] = df.apply(lambda x: Calc_STP_AB(x['iMetTempCorrC'], x['iMetPressMB'], x['Surf_Area']), axis=1)
df['Dry Volume'] = df.apply(lambda x: CalcVolume(x['Composition'], x['B1'], x['B2'], x['B3'], x['B4'], x['B5'], x['B6'], x['B7'], x['B8'], x['B9'], x['B10'], x['B11'], x['B12'], x['B13'], x['B14'], x['B15']), axis=1)
df['Effective_Radius'] = df.apply(lambda x: eRadCalc(x['Composition'], x['B1'], x['B2'], x['B3'], x['B4'], x['B5'], x['B6'], x['B7'], x['B8'], x['B9'], x['B10'], x['B11'], x['B12'], x['B13'], x['B14'], x['B15']), axis=1)

df['AV_q1'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Dry Volume'].transform(lambda x: x.quantile(0.25))
df['AV_q2'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Dry Volume'].transform(lambda x: x.quantile(0.5))
df['AV_q3'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Dry Volume'].transform(lambda x: x.quantile(0.75))

df['Aconc_q1'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Aer_Conc'].transform(lambda x: x.quantile(0.25))
df['Aconc_q2'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Aer_Conc'].transform(lambda x: x.quantile(0.5))
df['Aconc_q3'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Aer_Conc'].transform(lambda x: x.quantile(0.75))

df['AconcSTP_q1'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Aer_ConcAtSTP'].transform(lambda x: x.quantile(0.25))
df['AconcSTP_q2'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Aer_ConcAtSTP'].transform(lambda x: x.quantile(0.5))
df['AconcSTP_q3'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Aer_ConcAtSTP'].transform(lambda x: x.quantile(0.75))

df['Asa_q1'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Surf_AreaAtSTP'].transform(lambda x: x.quantile(0.25))
df['Asa_q2'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Surf_AreaAtSTP'].transform(lambda x: x.quantile(0.5))
df['Asa_q3'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Surf_AreaAtSTP'].transform(lambda x: x.quantile(0.75))

df['Erad_q1'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Effective_Radius'].transform(lambda x: x.quantile(0.25))
df['Erad_q2'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Effective_Radius'].transform(lambda x: x.quantile(0.5))
df['Erad_q3'] = df.groupby(['Launch', 'Altitude (km)', 'Flight'])['Effective_Radius'].transform(lambda x: x.quantile(0.75))

df = df.groupby(['Launch','Altitude (km)', 'Flight'],as_index=False).mean()

df = df.filter(['Launch', 'Altitude (km)', 'Flight', 'Air Temp K', 'GPSlat', 'GPSlon','GPSaltkm', 'POPSflowCCpS', 
                    'POPStempC', 'POPSbaselineStdDev','POPSavgWdithUS', 'POPSTempCorrection', 'POPSflowCorrCCpS',
                    'B1', 'B2', 'B3','B4','B5','B6','B7','B8','B9','B10','B11', 'B12', 
                    'B13', 'B14','B15', 'Dry Extinction', 'Wavelength', 'Aer_Conc', 'Aer_Conc_205nm', 'SurfArea', 'Dry Volume', 
                    'Aconc_q1', 'Aconc_q2', 'Aconc_q3', 'AconcSTP_q1', 'AconcSTP_q2', 'AconcSTP_q3', 'Asa_q1', 'Asa_q2', 'Asa_q3', 'AV_q1', 'AV_q2', 'AV_q3', 'Erad_q1', 'Erad_q2', 'Erad_q3'], axis=1)
            
df_p = df
del df

header_B2SAP_POPS = ['\n' + 'NOAA/GMD balloon flight ascent data file' + '\n' + 'SkySonde Client Version:' + ',' + '1.3.1.3' + '\n'
          + 'Balloon Flight Name:' + ',' + file[0:24]+ '\n' + 'Attached Instruments:' + ',' + 'Radiosonde' + ',' + 'EN-SCI ECC Ozonesonde' + ',' + 'NOAA Frostpoint Hygrometer' + ',' 
          + 'POPS' + '\n' + 'Station Name:' + ',' + 'Maido Observatory' + '\n' + 'Station Latitude:' + ',' + '-21.07972' + '\n' + 'Station Longitude:'
          + ',' +  '55.38278' + '\n' + 'Station Altitude:' + '2182.00 [m]' + '\n' + 'Station Pressure:' + ',' + '-99999 [mb]' + '\n' + 'Station Temperature'
          + ',' +  '99999 [deg C]' + '\n' + 'Station Humidity' + ',' + '-99999 [% RH]' + '\n' + 'Radiosonde Serial Number / ID:' + ',' + 'Pressure/Altitude Source:' 
          + ',' +  'iMet Radiosonde' + '\n' + 'Solar Radiation Correction:' + ',' + 'Yes' + '\n' + 'POPS Number of Bins:' + ',' + '15' + ',' + 'Units: Counts' +'\n' 
          + 'Mean Dp [um] (rounded to three significant digits assuming IOR = 1.45 + 0j):,0.150, 0.182, 0.221, 0.267, 0.326, 0.406, 0.493, 0.590, 0.722, 0.905, 1.20, 1.23, 1.32, 1.86, 2.24,' 
          + '\n' + 'dLogDp (rounded to four significant digits assuming IOR = 1.45 + 0j):, 0.0835, 0.0835, 0.0835, 0.0835, 0.0835, 0.0835, 0.0835, 0.0835, 0.0835, 0.0835, 0.0835, 0.0835, 0.0835, 0.0835, 0.0835' 
          + '\n' + 'Variable Descriptions:, Date [YYYY-mm-DD GMT], Elapsed minute since launch, iMET Ascent Rate [m/s], iMET Potential Temperature [K], iMET Pressure [MB],' 
          + 'iMet Air Temperature [deg C], iMet Altitude [km], FPH H2O Mixing Ratio [PPMV], FPH RH [%], FPH (GG) H2O partial press. [hPa], GPS Latitiude, GPS Longitude, GPS Altitude [km], GPS Wind Speed [m/s],' 
          + 'GPS Wind Direction [Degrees], POPS Flow [cc/s], POPS Laminar Flow Element Temperature [deg C], POPS Standard Deviation, POPS Mean Particle Width [uS],' + '\n']
        
# newpath = '/Users/asher/Documents/PapersInProgress/JGR_ANYSO/Submission/POPSDataForPublication/POPS/' + file
        

# with open(newpath, 'w') as ict:
#     for line in header_B2SAP_POPS:
#         ict.write(line)            
# df_p.to_csv(newpath,  mode ='a', float_format='%g') #write a new file for each flight (same header)


#read in Boulder Unperturbed median +/- IQR
#aerosol
Boulder_file = '/Users/asher/Documents/PyroCBPaper/Boulder/Background_IQR.xlsx'
df_Boulder=pd.read_excel(Boulder_file, index_col=0)  
#water vapor, ozone, Temperature
df_wv_Boulder=pd.read_csv('/Users/asher/Documents/PyroCBPaper/Boulder/Background_IQR_H2O.csv', delimiter = ',', header=0, dtype=None, low_memory=False)
df_o3_Boulder=pd.read_csv('/Users/asher/Documents/PyroCBPaper/Boulder/Background_IQR_O3.csv', delimiter = ',', header=0, dtype=None, low_memory=False)
df_t_Boulder=pd.read_csv('/Users/asher/Documents/PyroCBPaper/Boulder/Background_IQR_T.csv', delimiter = ',', header=0, dtype=None, low_memory=False)

df_wv_Boulder2=pd.read_csv('/Users/asher/Documents/PyroCBPaper/Boulder/Background_IQR_H2O_TA.csv', delimiter = ',', header=0, dtype=None, low_memory=False)
df_o3_Boulder2=pd.read_csv('/Users/asher/Documents/PyroCBPaper/Boulder/Background_IQR_O3_TA.csv', delimiter = ',', header=0, dtype=None, low_memory=False)


#merge fle2 data and csv data
df_AscDesc = pd.merge(df_p, df_w, how='left', on=['Launch', 'Altitude (km)', 'Flight'])
#df_AscDesc = pd.merge(df_w, df_p, how='left', on=['Launch', 'Altitude (km)', 'Flight'])



Wavelength1 = 532
Wavelength2 = 1064
df_AscDesc['Wavelength1'] = 532
df_AscDesc['Wavelength2'] = 1064
df_AscDesc['Composition'] = df_AscDesc.apply(lambda x: Plume(x['Launch'], x['Flight'], x['Altitude (km)']), axis=1)

df_AscDesc[['RH_FPHc', 'Pw_frostpoint']] = df_AscDesc.apply(lambda x: pd.Series(GoffGratch(x['TFpHyg'], x['Tcorr_degC'])), axis=1)
df_AscDesc['H2SO4_%wt_amb'] = df_AscDesc.apply(lambda x: PwLookuptable(x['Tcorr_degC'], x['Pw_frostpoint']), axis = 1)
df_AscDesc['H2SO4_%wt_pops'] = df_AscDesc.apply(lambda x: PwLookuptable(x['POPStempC'], x['Pw_frostpoint']), axis = 1)
#density calculation
df_AscDesc[['Density_amb', 'Density_amb_unc']] = df_AscDesc.apply(lambda x: pd.Series(DensityCalc(x['H2SO4_%wt_amb'], x['Tcorr_degC'])), axis = 1)
df_AscDesc[['Density_pops', 'Density_pops_unc']] = df_AscDesc.apply(lambda x: pd.Series(DensityCalc(x['H2SO4_%wt_pops'], x['POPStempC'])), axis = 1)

df_AscDesc['BinDp_S&H1981'] = df_AscDesc.apply(lambda x: Steele_Hamill_1981_J1995_calcDamb(x['Composition'], x['H2SO4_%wt_pops'], x['Density_pops'], x['H2SO4_%wt_amb'], x['Density_amb'], x['B1'], x['B2'], x['B3'], x['B4'], x['B5'], x['B6'], x['B7'], x['B8'], x['B9'], x['B10'], x['B11'], x['B12'], x['B13'], x['B14'], x['B15']), axis = 1)

plotIOR(df_AscDesc)

df_AscDesc[['Measured Extinction_532', 'Ambient Extinction_532']] = df_AscDesc.apply(lambda x: pd.Series(VPsca(x['Wavelength1'], x['Composition'], x['RH'], x['Air Temp (K)'], x['H2SO4_%wt_pops'], x['Density_pops'], x['H2SO4_%wt_amb'], x['Density_amb'], x['B1'], x['B2'], x['B3'], x['B4'], x['B5'], x['B6'], x['B7'], x['B8'], x['B9'], x['B10'], x['B11'], x['B12'], x['B13'], x['B14'], x['B15'])), axis=1)
df_AscDesc[['Measured Extinction_1064', 'Ambient Extinction_1064']] = df_AscDesc.apply(lambda x: pd.Series(VPsca(x['Wavelength2'], x['Composition'], x['RH'], x['Air Temp (K)'], x['H2SO4_%wt_pops'], x['Density_pops'], x['H2SO4_%wt_amb'], x['Density_amb'], x['B1'], x['B2'], x['B3'], x['B4'], x['B5'], x['B6'], x['B7'], x['B8'], x['B9'], x['B10'], x['B11'], x['B12'], x['B13'], x['B14'], x['B15'])), axis=1)



# #calculating extinction using the measured (dehydrated size distributions, as well as using KK theory and Steele and Hamill 1981 methods; discussed at length in Methods section Asher et al.). Done using the two different wavelengths defined above (532 and 997)
# new_df[['Extinction_m', 'Extinction_SH']] = new_df.apply(lambda x: pd.Series(VPsca(x['Wavelength'],  x['Composition'],x['RH_frostpoint'], x['Air Temp (K)'], x['H2SO4_%wt_pops'], x['Density_pops'], x['H2SO4_%wt_amb'], x['Density_amb'],
#                                                     x['B1'], x['B2'], x['B3'], x['B4'], x['B5'], x['B6'], x['B7'], x['B8'], x['B9'], x['B10'], x['B11'], 
#                                               x['B12'], x['B13'], x['B14'], x['B15'], x['B16'], x['B17'], x['B18'], x['B19'], x['B20'])), axis=1)


df_AscDesc['dNdLogDp'] = df_AscDesc.apply(lambda x: DnDlogDp(x['Launch'], x['Composition'], x['B1'], x['B2'], x['B3'], x['B4'], x['B5'], x['B6'], x['B7'], x['B8'], x['B9'], x['B10'], x['B11'], x['B12'], x['B13'], x['B14'], x['B15']), axis=1)
df_AscDesc['Dp'] = df_AscDesc.apply(lambda x: DiaLookUpTable1(x['Composition'], x['B1'], x['B2'], x['B3'], x['B4'], x['B5'], x['B6'], x['B7'], x['B8'], x['B9'], x['B10'], x['B11'], x['B12'], x['B13'], x['B14'], x['B15']), axis=1)


df_POPSsizedist = df_AscDesc.filter(['Launch', 'Altitude (km)', 'Flight', 'Dp', 'dNdLogDp', 'Theta [K]'])



#import lidar data
os.chdir('/Users/asher/Documents/PyroCBPaper/Lidar/')
file = 'Lauder_aerosol_lidar532_POPSdates_2019_22.csv'

Lidar_colnames = ['Altitude (km)', '30-Apr-2019', '03-Sep-2019',
       '27-Jan-2020', '19-Jun-2020','23-Jul-2020', '11-Nov-2020',
       '02-Sep-2021', '04-Nov-2021', '24-Jan-2022', 'Unc 30-Apr-2019', 'Unc 03-Sep-2019',
              'Unc 27-Jan-2020', 'Unc 19-Jun-2020','Unc 23-Jul-2020', 'Unc 11-Nov-2020',
              'Unc 02-Sep-2021', 'Unc 04-Nov-2021', 'Unc 24-Jan-2022']

df_Lidar=pd.read_csv(file, delimiter=",",  names = Lidar_colnames, header = 1, dtype=None, low_memory=False) #names = colNames_CESM_CARMA,

Lidar_Obs = ['30-Apr-2019', '03-Sep-2019',
       '27-Jan-2020', '19-Jun-2020','23-Jul-2020', '11-Nov-2020',
       '02-Sep-2021', '04-Nov-2021', '24-Jan-2022']

df_Lidar_unc = df_Lidar.filter(['Altitude (km)', 'Unc 30-Apr-2019', 'Unc 03-Sep-2019',
       'Unc 27-Jan-2020', 'Unc 19-Jun-2020','Unc 23-Jul-2020', 'Unc 11-Nov-2020',
       'Unc 02-Sep-2021', 'Unc 04-Nov-2021', 'Unc 24-Jan-2022'])

df_Lidar = df_Lidar.filter(['Altitude (km)', '30-Apr-2019', '03-Sep-2019',
       '27-Jan-2020', '19-Jun-2020','23-Jul-2020', '11-Nov-2020',
       '02-Sep-2021', '04-Nov-2021', '24-Jan-2022'])

df_Lidar_unc.rename(columns={"Unc 30-Apr-2019": "30-Apr-2019", "Unc 03-Sep-2019": "03-Sep-2019", "Unc 27-Jan-2020": "27-Jan-2020", "Unc 19-Jun-2020": "19-Jun-2020", "Unc 23-Jul-2020": "23-Jul-2020", "Unc 11-Nov-2020": "11-Nov-2020", "Unc 02-Sep-2021": "02-Sep-2021","Unc 04-Nov-2021": "04-Nov-2021",  "Unc 24-Jan-2022": "24-Jan-2022"}, inplace = True)

dfl_int_1 = pd.melt(df_Lidar, id_vars='Altitude (km)', value_vars= Lidar_Obs, var_name = 'GMTdateDmY', value_name = 'Extinction_per_quarter_km') 
dfl_int_1u = pd.melt(df_Lidar_unc, id_vars='Altitude (km)', value_vars= Lidar_Obs, var_name = 'GMTdateDmY', value_name = 'Unc Extinction_per_quarter_km') 
dfl_int = pd.merge(dfl_int_1, dfl_int_1u, how='left', on=['GMTdateDmY', 'Altitude (km)'])

dfl_int[ dfl_int['Extinction_per_quarter_km'] < 0] = 0 #assign extrintion less than zero equal to zero...
dfl_int[ dfl_int['Unc Extinction_per_quarter_km'] < 0] =  dfl_int[ dfl_int['Unc Extinction_per_quarter_km'] < 0] * -1 #assign extrintion less than zero equal to zero...
dfl_int['sAOD'] = dfl_int['Extinction_per_quarter_km']* 0.25 #this will get summed momentarily (units are km-1)
dfl_int['sAOD unc'] = dfl_int['Unc Extinction_per_quarter_km'] * 0.25 #this will get summed momentarily (units are km-1)
#
#dfl_int = dfl_int[dfl_int['Altitude (km)'] <= 30.0] #so mmiddle stratosphere is comparable to POPS...
#dfl_int['Tropopause'] = dfl_int.apply(lambda x: Lidar_Trop(x['Launch']), axis=1)
#dfl_int['Potential T [K]'] = dfl_int.apply(lambda x: Lidar_PotT(x['Launch']), axis=1)


#assign closest POPS launch date as a column
#create a launch date so that if a flight covers two days, all gets included as one "launch"
conditions = [(dfl_int['GMTdateDmY'] == '30-Apr-2019'), (dfl_int['GMTdateDmY'] == '31-Aug-2019'), (dfl_int['GMTdateDmY'] == '01-Sep-2019'),
       (dfl_int['GMTdateDmY'] == '02-Sep-2019'), (dfl_int['GMTdateDmY'] == '03-Sep-2019'), (dfl_int['GMTdateDmY'] == '09-Sep-2019'), (dfl_int['GMTdateDmY'] == '20-Jan-2020'),
       (dfl_int['GMTdateDmY'] == '22-Jan-2020'), (dfl_int['GMTdateDmY'] == '23-Jan-2020'), (dfl_int['GMTdateDmY'] == '25-Jan-2020'), (dfl_int['GMTdateDmY'] == '26-Jan-2020'),
       (dfl_int['GMTdateDmY'] == '27-Jan-2020'), (dfl_int['GMTdateDmY'] == '29-Jan-2020'), (dfl_int['GMTdateDmY'] == '30-Jan-2020'), (dfl_int['GMTdateDmY'] == '31-Jan-2020'),
       (dfl_int['GMTdateDmY'] == '01-Feb-2020'), (dfl_int['GMTdateDmY'] == '02-Feb-2020'), (dfl_int['GMTdateDmY'] == '03-Feb-2020'), (dfl_int['GMTdateDmY'] == '19-Jun-2020'),
       (dfl_int['GMTdateDmY'] == '23-Jun-2020'), (dfl_int['GMTdateDmY'] == '23-Jul-2020'), (dfl_int['GMTdateDmY'] == '27-Jul-2020'), (dfl_int['GMTdateDmY'] == '11-Nov-2020'),
       (dfl_int['GMTdateDmY'] == '02-Sep-2021'), (dfl_int['GMTdateDmY'] == '04-Nov-2021'), (dfl_int['GMTdateDmY'] == '17-Jan-2022'), (dfl_int['GMTdateDmY'] == '21-Jan-2022'), (dfl_int['GMTdateDmY'] == '22-Jan-2022'),
       (dfl_int['GMTdateDmY'] == '23-Jan-2022'), (dfl_int['GMTdateDmY'] == '24-Jan-2022'), (dfl_int['GMTdateDmY'] == '25-Jan-2022'), (dfl_int['GMTdateDmY'] == '26-Jan-2022'),
       (dfl_int['GMTdateDmY'] == '27-Jan-2022'), (dfl_int['GMTdateDmY'] == '28-Jan-2022'), (dfl_int['GMTdateDmY'] == '30-Jan-2022'), (dfl_int['GMTdateDmY'] == '31-Jan-2022')]


choices = ['2019-04-30', '2019-09-03', '2019-09-03', '2019-09-03', '2019-09-03', '2019-09-03',
            '2020-01-27', '2020-01-27', '2020-01-27', '2020-01-27', '2020-01-27', '2020-01-27', '2020-01-27', '2020-01-27', '2020-01-27', '2020-01-27', '2020-01-27', '2020-01-27',
            '2020-06-19', '2020-06-19', '2020-07-23','2020-07-23', '2020-11-17',  '2021-08-31', '2021-11-03', 
            '2022-01-24', '2022-01-24', '2022-01-24', '2022-01-24', '2022-01-24', '2022-01-24', '2022-01-24', '2022-01-24', '2022-01-24', '2022-01-24', '2022-01-24']

dfl_int['Launch']=np.select(conditions, choices, default = 'Unassigned')
dfl_int['AtmosLayer'] = dfl_int.apply(lambda x: StratTropLidar(x['Launch'], x['Altitude (km)']), axis = 1)
df_psub = df_AscDesc.filter(['Launch', 'Flight', 'Altitude (km)', 'Theta [K]', 'Tropopause'], axis=1)
df_psub = df_psub[df_psub['Flight'] == 'Ascent']
df_psub['AtmosLayer'] = df_psub.apply(lambda x: StratTrop(x['Theta [K]'], x['Altitude (km)'], x['Tropopause']), axis=1)
df_psub = df_psub.loc[[df_psub.loc[df_psub.AtmosLayer == 'Stratosphere', 'Altitude (km)'].idxmin()]]
#dfl_f = pd.merge(dfl_int,df_psub, how='left', on=['Launch', 'Altitude (km)'])
dfl_f = dfl_int

dfl_s = dfl_f.groupby(['GMTdateDmY', 'AtmosLayer'], as_index=False).sum()
dfl_s = dfl_s[dfl_s['AtmosLayer'] != 'Troposphere']  
Column =['sAOD']
Column2 =['sAOD unc']
#rename datestr
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("30-Apr-2019", "2019-04-30", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("31-Aug-2019", "2019-08-31", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("01-Sep-2019", "2019-09-01", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("02-Sep-2019", "2019-09-02", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("03-Sep-2019", "2019-09-03", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("09-Sep-2019", "2019-09-09", case = False) 

dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("20-Jan-2020", "2020-01-20", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("22-Jan-2020", "2020-01-22", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("23-Jan-2020", "2020-01-23", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("25-Jan-2020", "2020-01-25", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("26-Jan-2020", "2020-01-26", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("27-Jan-2020", "2020-01-27", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("29-Jan-2020", "2020-01-29", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("30-Jan-2020", "2020-01-30", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("31-Jan-2020", "2020-01-31", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("01-Feb-2020", "2020-02-01", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("02-Feb-2020", "2020-02-02", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("03-Feb-2020", "2020-02-03", case = False) 

dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("19-Jun-2020", "2020-06-19", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("23-Jun-2020", "2020-06-23", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("23-Jul-2020", "2020-07-23", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("27-Jul-2020", "2020-07-27", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("11-Nov-2020", "2020-11-11", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("02-Sep-2021", "2021-09-02", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("03-Nov-2021", "2021-11-03", case = False) 

dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("17-Jan-2022", "2022-01-17", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("21-Jan-2022", "2022-01-21", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("22-Jan-2022", "2022-01-22", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("23-Jan-2022", "2022-01-23", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("24-Jan-2022", "2022-01-24", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("25-Jan-2022", "2022-01-25", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("26-Jan-2022", "2022-01-26", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("27-Jan-2022", "2022-01-27", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("28-Jan-2022", "2022-01-28", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("30-Jan-2022", "2022-01-30", case = False) 
dfl_s["GMTdateDmY"]= dfl_s["GMTdateDmY"].str.replace("31-Jan-2022", "2022-01-31", case = False) 
dfl_s["GMTdateYmD"] = dfl_s["GMTdateDmY"]

dfl_sf = pd.melt(dfl_s, id_vars=['GMTdateYmD', 'AtmosLayer'], value_vars= Column, var_name = 'Property', value_name = 'Vertically Integrated Total') 
dfl_sf['Date'] = pd.to_datetime(dfl_sf['GMTdateYmD']).apply(lambda date: date.toordinal())
dfl_sf['Launch'] = dfl_sf['Date']
dfl_sf['Source'] = 'Lidar'

dfl_sf_u = pd.melt(dfl_s, id_vars=['GMTdateYmD', 'AtmosLayer'], value_vars= Column2, var_name = 'Property', value_name = 'Vertically Integrated Total') 
dfl_sf_u['Date'] = pd.to_datetime(dfl_sf_u['GMTdateYmD']).apply(lambda date: date.toordinal())
dfl_sf_u['Source'] = 'Lidar'
    
    
            
   
df = df_AscDesc.groupby(['Launch','Altitude (km)'],as_index=False).median(numeric_only= True)
df = df_AscDesc[df_AscDesc['Flight'] == 'Ascent']
df2 = df
#df['AtmosLayer'] = df.apply(lambda x: StratTrop(x['Theta [K]'], x['Altitude (km)'], x['Tropopause']), axis=1)
df2 = df[df['Altitude (km)'] >= 10]
#for travis show a single example 
#df = df[df.Launch == '2020-01-27']
#for line plots (filter, make sure altitude for y-axis is numeric)

#df['Altitude (km)'] = pd.to_numeric(df['Altitude (km)'], errors='coerce').fillna(0)

#define a baseline - removing obviously/freshly pertubed profiles
df_base = df
df_base = df[df.Launch != "2019-09-03"]
df_base = df[df.Launch != "2020-01-27"]
df_base['AC_base_q1'] = df_base.groupby(['Altitude (km)'])['Aconc_q2'].transform(lambda x: x.quantile(0.25))
df_base['AC_base_q2'] = df_base.groupby(['Altitude (km)'])['Aconc_q2'].transform(lambda x: x.quantile(0.5))
df_base['AC_base_q3'] = df_base.groupby(['Altitude (km)'])['Aconc_q2'].transform(lambda x: x.quantile(0.75))

df_base['ACSTP_base_q1'] = df_base.groupby(['Altitude (km)'])['AconcSTP_q2'].transform(lambda x: x.quantile(0.25))
df_base['ACSTP_base_q2'] = df_base.groupby(['Altitude (km)'])['AconcSTP_q2'].transform(lambda x: x.quantile(0.5))
df_base['ACSTP_base_q3'] = df_base.groupby(['Altitude (km)'])['AconcSTP_q2'].transform(lambda x: x.quantile(0.75))

df_base['ASA_base_q1'] = df_base.groupby(['Altitude (km)'])['Asa_q2'].transform(lambda x: x.quantile(0.25))
df_base['ASA_base_q2'] = df_base.groupby(['Altitude (km)'])['Asa_q2'].transform(lambda x: x.quantile(0.5))
df_base['ASA_base_q3'] = df_base.groupby(['Altitude (km)'])['Asa_q2'].transform(lambda x: x.quantile(0.75))

df_base['ER_base_q1'] = df_base.groupby(['Altitude (km)'])['Erad_q2'].transform(lambda x: x.quantile(0.25))
df_base['ER_base_q2'] = df_base.groupby(['Altitude (km)'])['Erad_q2'].transform(lambda x: x.quantile(0.5))
df_base['ER_base_q3'] = df_base.groupby(['Altitude (km)'])['Erad_q2'].transform(lambda x: x.quantile(0.75))


df_base['H2O_base_q1'] = df_base.groupby(['Altitude (km)'])['H2O_q2'].transform(lambda x: x.quantile(0.25))
df_base['H2O_base_q2'] = df_base.groupby(['Altitude (km)'])['H2O_q2'].transform(lambda x: x.quantile(0.5))
df_base['H2O_base_q3'] = df_base.groupby(['Altitude (km)'])['H2O_q2'].transform(lambda x: x.quantile(0.75))

df_base['O3_base_q1'] = df_base.groupby(['Altitude (km)'])['O3_q2'].transform(lambda x: x.quantile(0.25))
df_base['O3_base_q2'] = df_base.groupby(['Altitude (km)'])['O3_q2'].transform(lambda x: x.quantile(0.5))
df_base['O3_base_q3'] = df_base.groupby(['Altitude (km)'])['O3_q2'].transform(lambda x: x.quantile(0.75))

df_base['T_base_q1'] = df_base.groupby(['Altitude (km)'])['T_q2'].transform(lambda x: x.quantile(0.25))
df_base['T_base_q2'] = df_base.groupby(['Altitude (km)'])['T_q2'].transform(lambda x: x.quantile(0.5))
df_base['T_base_q3'] = df_base.groupby(['Altitude (km)'])['T_q2'].transform(lambda x: x.quantile(0.75))

#tropopause relevant coordinates
df_base['H2O_base_q1_TA'] = df_base.groupby(['Tropopause_Altitude (km)'])['H2O_q2_TA'].transform(lambda x: x.quantile(0.25))
df_base['H2O_base_q2_TA'] = df_base.groupby(['Tropopause_Altitude (km)'])['H2O_q2_TA'].transform(lambda x: x.quantile(0.5))
df_base['H2O_base_q3_TA'] = df_base.groupby(['Tropopause_Altitude (km)'])['H2O_q2_TA'].transform(lambda x: x.quantile(0.75))

df_base['O3_base_q1_TA'] = df_base.groupby(['Tropopause_Altitude (km)'])['O3_q2_TA'].transform(lambda x: x.quantile(0.25))
df_base['O3_base_q2_TA'] = df_base.groupby(['Tropopause_Altitude (km)'])['O3_q2_TA'].transform(lambda x: x.quantile(0.5))
df_base['O3_base_q3_TA'] = df_base.groupby(['Tropopause_Altitude (km)'])['O3_q2_TA'].transform(lambda x: x.quantile(0.75))

df_base2 = df_w_TA.groupby(['Tropopause_Altitude (km)'], as_index=False).median(numeric_only= True)
#df_base = df_w.groupby(['Altitude (km)'], as_index=False).median(numeric_only= True)
#only show above the tropopause


#df_extinction = df2.filter(['Launch','Altitude (km)', 'GPSlat', 'GPSlon', 'GPSaltkkm','Measured Extinction_532', 'Measured Extinction_1064', 'Ambient Extinction_532', 'Ambient Extinction_1064', 'Tropopause'])
df2['Measured Backscattering_532'] = df2['Measured Extinction_532']
df2['Ambient Backscattering_532'] = df2['Ambient Extinction_532']
df_extinction = df2.filter(['Launch','Altitude (km)', 'GPSlat', 'GPSlon', 'GPSaltkkm','Measured Extinction_532', 'Ambient Extinction_532', 'Tropopause'])
#df_extinction.to_csv('/Users/asher/Desktop/LauderExtinction_UTLS.csv', mode='w', float_format='%g')

#df_extinction_fullprofiles = df_AscDesc.filter(['Launch',  'Flight', 'Altitude (km)', 'GPSlat', 'GPSlon', 'GPSaltkkm','Measured Extinction_532', 'Measured Extinction_1064', 'Ambient Extinction_532', 'Ambient Extinction_1064', 'Tropopause'])
df_AscDesc['Measured Backscattering_532'] = df_AscDesc['Measured Extinction_532']
df_AscDesc['Ambient Backscattering_532'] = df_AscDesc['Ambient Extinction_532']
df_extinction_fullprofiles = df_AscDesc.filter(['Launch',  'Flight', 'Altitude (km)', 'GPSlat', 'GPSlon', 'GPSaltkkm','Measured Extinction_532', 'Ambient Extinction_532', 'Tropopause'])
#df_extinction_fullprofiles.to_csv('/Users/asher/Desktop/LauderExtinction_fullprofiles.csv', mode='w', float_format='%g')


df5 = df_AscDesc[df_AscDesc['Altitude (km)'] >= 0.5]
#these are automatically sized (log signal space bin means)
x_axis_labels = np.array(['154','190','230','270','330','400', '490','590','720',\
                          '870','1050','1280','1550', '1870','2270'])

x_axis_smoke_labels = np.array(['133', '145', '160', '179', '203', '235', '324',\
                                '456', '552', '684', '1348', '1663', '2541', '3263'])


    
    
#assume sulfate based sizes...
#ndp = [0.407, 0.459, 0.404, 1.205, 1.195, 1.730, 1.795, 0.771, 0.270, 0.117, 0.0390, 0.0, 0.0, 0.0, 0.0]


#df4['Composition'] = df4.apply(lambda x: Plume(x['Flight'], x['Altitude (km)']), axis=1)
df5['Composition'] = df5.apply(lambda x: Plume(x['Launch'], x['Flight'], x['Altitude (km)']), axis=1)
#df6['Composition'] = df6.apply(lambda x: Plume(x['Flight'], x['Altitude (km)']), axis=1)
#df7['Composition'] = df7.apply(lambda x: Plume(x['Flight'], x['Altitude (km)']), axis=1)



#Compare Size distributions

ThetaLabels = ['310', '330', '350', '370', '390', '410', '430', '450', '470', '490', '510', '530', '550', '570', '590', '610', '630', '650', '670', '690', '710', '730', '750', '770', '790']

theta_a =   300 #290 #295 # 307.5 
theta_b =   800 #810 #795 # 782.5 
theta_bins =  11 #27 # # 27 

theta_a =   300 #290 #295 # 307.5 
theta_b =   800 #810 #795 # 782.5 
theta_bins =  26 #27 # # 27 

os.chdir('/Users/asher/Documents/PyroCBPaper/CESM_CARMA/neardatesNewest/SD/')
filelist = glob.glob('*.txt')
path = '/Users/asher/Documents/PyroCBPaper/CESM_CARMA/neardatesNewest/SD/'
CARMAsizedist= pd.DataFrame()
# SD_colnames = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10',
#                 'B11', 'B12', 'B13', 'B14', 'B15', 'B16', 'B17', 'B18', 'B19', 'B20']
SD_colnames = [100, 131.2, 172.2, 226, 296, 388, 510, 670, 878, 1152,
                  1512, 1984, 2603.22, 3415.62, 4481.56, 5880.16, 7715.22, 10123, 13282.1, 17427.1]
# SD_colnam_df = pd.read_csv('/Users/asher/Documents/PyroCBPaper/CESM_CARMA/10dayaverages/radius_unit_cm.csv') * 1E4 * 2
# SD_colnames = SD_colnam_df.values.tolist()

for file in filelist:    
    df1=pd.read_csv(file,sep=r'(?:,|\s+)',names = SD_colnames, dtype=None, engine='python')
    #add date
    date = (file[11:15]+ '-'+ file[15:17]+ '-'+ file[17:19])
    df1['Launch'] = date
    df1['Level'] = df1.index
    #CARMAsizedist = CARMAsizedist.append(df1) #change to pd.concat (ignore index)
    CARMAsizedist = pd.concat([CARMAsizedist, df1], ignore_index = True, axis=0)
    del df1
    
#join CARMA datasets
CARMAsizedist.reset_index(drop=True, inplace=True)
df_CARMAvp.reset_index(drop=True, inplace=True)
CARMA = df_CARMAvp.merge(CARMAsizedist, on=['Launch', 'Level'], how='left', sort=False)

#add header manually
#CARMA.to_csv('/Users/asher/Documents/PapersInProgress/JGR_ANYSO/Submission/POPSDataForPublication/CARMA_VP_table.csv', mode='a', float_format='%g')

    

CARMA_ac_test = CARMA
CARMA_ac_test['Aer Conc calc. SD >= 200 nm'] = CARMA_ac_test[[131.2, 172.2, 226.0, 296.0, 388.0, 510.0, 670.0, 878.0, 1152.0,
                  1512.0, 1984.0, 2603.22]].sum(axis=1) #100.0, 131.2, 3415.62, 4481.56, 5880.16, 7715.22, 10123.0, 13282.1, 17427.1
CARMA_ac_test['Aer Conc calc. SD > 140 nm'] = CARMA_ac_test[[172.2, 226.0, 296.0, 388.0, 510.0, 670.0, 878.0, 1152.0,
                  1512.0, 1984.0, 2603.22]].sum(axis=1)

CARMA_ac_test['Extinction calc. SD'] =  CARMA_ac_test.apply(lambda x: VPscaCARMA(x[100.0], x[131.2], x[172.2], x[226.0], x[296.0], x[388.0], x[510.0], x[670.0], x[878.0], x[1152.0], x[1512.0], x[1984.0], x[2603.22], x[3415.62], x[4481.56], x[5880.16], x[7715.22], x[10123.0], x[13282.1], x[17427.1]), axis=1)

CARMA_ac_test = CARMA_ac_test.filter(['Launch', 'Potential T [K]', 'GeoPotential Height [km]', 'diffs_alt', 'Aerosol Concentration [#/cm3]', 'Aer Conc calc. SD >= 200 nm', 'Aer Conc calc. SD > 140 nm', 'Extinction calc. SD', 'sAOD'], axis=1)

CARMA_ac_test['Number Concentration Column(# m^-2)'] = CARMA_ac_test['Aer Conc calc. SD >= 200 nm'] * 1E6 * CARMA_ac_test['diffs_alt'] *1E3 # not at STP
CARMA_ac_test['sAOD']  = CARMA_ac_test['Extinction calc. SD']* CARMA_ac_test['diffs_alt']
CARMA_ac_test['Tropopause'] = CARMA_ac_test.apply(lambda x: CARMA_Trop(x['Launch']), axis=1)
CARMA_ac_test['AtmosLayer'] = CARMA_ac_test.apply(lambda x: StratTrop(x['Potential T [K]'], x['GeoPotential Height [km]'], x['Tropopause']), axis=1)
#df_CARMAvp = df_CARMA #to merge with size distribution data

CARMA_ac_test = CARMA_ac_test[CARMA_ac_test['GeoPotential Height [km]'] > 2.0]
CARMA_ac_test = CARMA_ac_test[CARMA_ac_test['GeoPotential Height [km]'] <= 28.0]
CARMA_ac_test.drop(['Potential T [K]', 'GeoPotential Height [km]', 'Tropopause', 'Aerosol Concentration [#/cm3]', 'Aer Conc calc. SD > 140 nm', 'Extinction calc. SD'], axis = 1, inplace = True)
dfc_s2 = CARMA_ac_test.groupby(['Launch', 'AtmosLayer'], as_index=False).sum()
dfc_s2 = dfc_s2[dfc_s2['AtmosLayer'] != 'Troposphere']  
Column2 =['Number Concentration Column(# m^-2)', 'sAOD']
dfc_int2 = pd.melt(dfc_s2, id_vars=['Launch', 'AtmosLayer'], value_vars= Column2, var_name = 'Property', value_name = 'Vertically Integrated Total') 
dfc_int2['Date'] = pd.to_datetime(dfc_int2['Launch']).apply(lambda date: date.toordinal())
dfc_int2['Source'] = 'CARMA calc.'

dfc_int.rename(columns={"Vertically Integrated Total": "Vertically Integrated Total Min"}, inplace = True)
dfc_int["Vertically Integrated Total Max"] = dfc_int2["Vertically Integrated Total"]
dfc_int["Vertically Integrated Total"] = (dfc_int["Vertically Integrated Total Min"] + dfc_int["Vertically Integrated Total Max"])/2


CARMA.drop(['GeoPotential Height [km]', 'Aerosol Concentration [#/cm3]', 'Extinction [km-1]', 'diffs_alt', 'Number Concentration Column(# m^-2)', 'sAOD', 'Tropopause', 'AtmosLayer', 'Level' , 'Potential T [K]'], axis=1)
#bin by altitude of pot T
CARMA['Theta'] = pd.cut(CARMA['Potential T [K]'], bins=np.linspace(theta_a,theta_b,theta_bins),labels = ThetaLabels, include_lowest=True) #(290, 810, 27)
CARMA['Theta'] = pd.to_numeric(CARMA ['Theta'], errors='coerce').fillna(0)
#CARMA.to_csv('/Users/asher/Documents/PyroCBPaper/CESM_CARMA/10dayaverages/CARMAmerge.csv', index=True)
 
CARMA['Altitude (km)'] = CARMA['GeoPotential Height [km]']
CARMA = CARMA.drop(['GeoPotential Height [km]', 'Aerosol Concentration [#/cm3]', 'Extinction [km-1]', 'diffs_alt', 'Number Concentration Column(# m^-2)', 'sAOD', 'Tropopause', 'AtmosLayer', 'Level' , 'Potential T [K]'], axis=1)


CARMA = CARMA[((CARMA['Theta'] == 410) | (CARMA['Theta'] == 510) | (CARMA['Theta'] == 610))]
CARMA = CARMA[(CARMA['Launch'] != '2020-06-19')]
#CARMA = CARMA[((((CARMA['Launch'] != '2019-09-03') & (CARMA['Launch'] != '2020-06-19')) & (CARMA['Launch'] != '2020-01-27')) & (CARMA['Launch'] != '2021-11-02')) & (CARMA['Launch'] != '2021-08-30')]
#plot number concentraiton etc...compared to specific profiles...
Carma_test = pd.melt(CARMA,id_vars=['Launch','Theta'], value_vars=SD_colnames, var_name = 'Dp', value_name = 'Concentration') #, 'Altitude (km)'
Carma_test['Dp'] = pd.to_numeric(Carma_test['Dp'])
Carma_test = Carma_test.groupby(['Launch', 'Theta', 'Dp'],as_index=False).mean()
Carma_test['Source'] = 'CARMA'
conditions = [(Carma_test['Launch'] == '2019-09-03'), (Carma_test['Launch'] == '2020-01-27'), (Carma_test['Launch'] == '2020-06-19'),
       (Carma_test['Launch'] == '2020-07-23'), (Carma_test['Launch'] == '2020-11-17')]   
choices = ['Volcanic', 'BB', 'Aged BB', 'Aged BB', 'Approaching Baseline']
Carma_test['Type']=np.select(conditions, choices, default = 'Unassigned')
Carma_test['Concentration'] = Carma_test['Concentration']/1.18E-1
Carma_test = Carma_test.filter(['Launch', 'Theta', 'Altitude (km)', 'Dp', 'Concentration', 'Source', 'Type'])


df_POPSsizedist['Theta'] = pd.cut(df_POPSsizedist['Theta [K]'], bins=np.linspace(theta_a,theta_b,theta_bins),labels = ThetaLabels, include_lowest=True) #(290, 810, 27)
df_POPSsizedist['Theta'] = pd.to_numeric(df_POPSsizedist['Theta'], errors='coerce').fillna(0)

df_POPSsizedist = df_POPSsizedist[((df_POPSsizedist['Theta'] == 410) | (df_POPSsizedist['Theta'] == 510) | (df_POPSsizedist['Theta'] == 610))]
#df_POPSsizedist = df_POPSsizedist[(((df_POPSsizedist['Launch'] != '2019-04-30') & (df_POPSsizedist['Launch'] != '2020-06-19')) & (df_POPSsizedist['Launch'] != '2021-08-31')) & (df_POPSsizedist['Launch'] != '2022-01-24') ]

df_POPSsizedist[['dNdLogDp1','dNdLogDp2', 'dNdLogDp3', 'dNdLogDp4', 'dNdLogDp5', 'dNdLogDp6', 'dNdLogDp7', 'dNdLogDp8', 'dNdLogDp9', 'dNdLogDp10', 'dNdLogDp11', 'dNdLogDp12', 'dNdLogDp13', 'dNdLogDp14', 'dNdLogDp15']] = pd.DataFrame(df_POPSsizedist.dNdLogDp.tolist(), index= df_POPSsizedist.index)
df_POPSsizedist[['Dp1','Dp2', 'Dp3', 'Dp4', 'Dp5', 'Dp6', 'Dp7', 'Dp8', 'Dp9', 'Dp10', 'Dp11', 'Dp12', 'Dp13', 'Dp14', 'Dp15']] = pd.DataFrame(df_POPSsizedist.Dp.tolist(), index= df_POPSsizedist.index)

POPSsd_m = df_POPSsizedist.groupby(['Launch', 'Theta'],as_index=False).mean()


POPSsd_m_conc = POPSsd_m.filter(['Launch', 'Theta', 'Altitude (km)', 'Theta [K]', 'dNdLogDp1',
       'dNdLogDp2', 'dNdLogDp3', 'dNdLogDp4', 'dNdLogDp5', 'dNdLogDp6',
       'dNdLogDp7', 'dNdLogDp8', 'dNdLogDp9', 'dNdLogDp10', 'dNdLogDp11',
       'dNdLogDp12', 'dNdLogDp13', 'dNdLogDp14', 'dNdLogDp15', 'Source', 'Type'])

POPSsd_m_dp = POPSsd_m.filter(['Launch', 'Theta', 'Altitude (km)', 'Theta [K]', 'Dp1', 'Dp2',
       'Dp3', 'Dp4', 'Dp5', 'Dp6', 'Dp7', 'Dp8', 'Dp9', 'Dp10', 'Dp11', 'Dp12',
       'Dp13', 'Dp14', 'Dp15', 'Source', 'Type'])

POPS_SD_colnames = ['dNdLogDp1','dNdLogDp2', 'dNdLogDp3', 'dNdLogDp4', 'dNdLogDp5', 'dNdLogDp6',
       'dNdLogDp7', 'dNdLogDp8', 'dNdLogDp9', 'dNdLogDp10', 'dNdLogDp11',
       'dNdLogDp12', 'dNdLogDp13', 'dNdLogDp14', 'dNdLogDp15']

POPS_DP_colnames = ['Dp1', 'Dp2',
       'Dp3', 'Dp4', 'Dp5', 'Dp6', 'Dp7', 'Dp8', 'Dp9', 'Dp10', 'Dp11', 'Dp12',
       'Dp13', 'Dp14', 'Dp15']

POPSsd_test_conc = pd.melt(POPSsd_m_conc,id_vars=['Launch','Theta','Altitude (km)'], value_vars=POPS_SD_colnames, var_name = 'Bins', value_name = 'Concentration')
POPSsd_test_dp = pd.melt(POPSsd_m_dp,id_vars=['Launch','Theta', 'Altitude (km)'], value_vars=POPS_DP_colnames, var_name = 'Diameter', value_name = 'Dp')
POPSsd_test_conc = POPSsd_test_conc.drop(['Launch','Theta', 'Altitude (km)'], axis = 1)

POPSsd_mf = pd.concat([POPSsd_test_dp, POPSsd_test_conc], axis = 1)
POPSsd_mf = POPSsd_mf.drop(['Bins', 'Diameter'], axis=1)

POPSsd_mf['Source'] = 'POPS'

#POPSsd_mf = POPSsd_mf [POPSsd_mf['Launch'] != '2020-06-19']
conditions = [(POPSsd_mf['Launch'] == '2019-09-03'), (POPSsd_mf['Launch'] == '2020-01-27'), (POPSsd_mf['Launch'] == '2020-06-19'),
       (POPSsd_mf['Launch'] == '2020-07-23'), (POPSsd_mf['Launch'] == '2020-11-17')]   
choices = ['Volcanic', 'BB', 'Other', 'Aged BB', 'Approaching Baseline']
POPSsd_mf['Type']=np.select(conditions, choices, default = 'Unassigned')


##vertcat model and measu size dist dataframes together
SizeDist = pd.concat([Carma_test, POPSsd_mf], ignore_index = True, axis=0)


palette1 = {"2019-04-30":"cadetblue",
            "2019-09-03":"black",
            "2020-01-27":"firebrick", 
            "2020-06-19":"peru", 
            "2020-07-23":"peru",
            "2020-11-17":"steelblue"
                }

palette1 = {"POPS":"black",
            "CARMA":"firebrick"
                }

#seaborn.set(rc={"xtick.bottom" : True, "ytick.left" : True})
#seaborn.set_style("whitegrid")
seaborn.set(font_scale=1, style="ticks")
g = seaborn.relplot(x="Dp", y="Concentration", style = 'Source', style_order=(['POPS','CARMA']),
                    hue = 'Source', hue_order = (['POPS','CARMA']),  #hue = "Launch", hue_order = ([ '2019-09-03', '2020-01-27', '2020-07-23', '2020-11-17']),
                    col = 'Theta', col_order = ([410, 510]), row = 'Type', row_order = ([ 'Volcanic', 'BB', 'Aged BB', 'Approaching Baseline']), kind="line",  marker = 's', palette=palette1, data= SizeDist, facet_kws={'sharey': True, 'sharex': True}) #drawstyle='steps-post', 
g.set_axis_labels("Dp (nm)", "dN/dlogDp")    
#g.set_axis_labels("Dp (nm)", "Aerosol Concentration (# $\mathregular{cm^{-3}}$)") 
g.set(xlim=(90, 3000))
g.set(ylim=(0.08, 200))
g.set(xscale="log")
g.set(yscale="log")
uticks = [0.1, 1, 10, 100]
g.set(yticks=uticks)
g.set(yticklabels=uticks)



g.axes[0,0].set_title('Theta: 400-420 K; Volcanic Plume')
g.axes[0,1].set_title('Theta: 500-520 K; Baseline - Above Volcanic Plume')
g.axes[1,0].set_title('Theta: 400-420 K; Biomass Burning Plume')
g.axes[1,1].set_title('Theta: 500-520 K; Above Biomass Burning Plume')
g.axes[2,0].set_title('Theta: 400-420 K; Aged Biomass Burning Plume')
g.axes[2,1].set_title('Theta: 500-520 K; Aged Biomass Burning Plume')
g.axes[3,0].set_title('Theta: 400-420 K; Approaching Baseline')
g.axes[3,1].set_title('Theta: 500-520 K; Approaching Baseline')
# SizeDist = b_final.append(Carma_test)
# SizeDist= SizeDist.groupby(['Dp', 'Launch', 'Theta', 'Source', 'Type'],as_index=False).mean()

# seaborn.set_style("whitegrid")
# seaborn.set(rc={"xtick.bottom" : True, "ytick.left" : True})
# g = seaborn.relplot(x="Dp", y="Concentration", style = 'Source', style_order=(['POPS','Model']),
#                     hue = "Launch", hue_order = (['2020-11-17', '2019-09-03', '2020-01-27', '2020-06-19']),
#                     col = 'Theta', col_order = ([420, 520]), row = 'Type', row_order = (['Approaching Baseline', 'Volcanic', 'BB', 'Aged BB']), kind="line", drawstyle='steps-pre', palette="colorblind", data=SizeDist, facet_kws={'sharey': True, 'sharex': True})
# g.set_axis_labels("Dp (nm)", "dN/dlogDp")    
# #g.set_axis_labels("Dp (nm)", "Aerosol Concentration (# $\mathregular{cm^{-3}}$)") 
# g.set(xlim=(90, 3000))
# g.set(ylim=(0.001, 100))
# g.set(xscale="log")
# g.set(yscale="log")













#figures below

custom_lines = [Line2D([0], [0], color="gray", lw=4),
                Line2D([0], [0], color="peru", lw=4),
                Line2D([0], [0], color="silver", lw=4),
                Line2D([0], [0], color="blue", lw=4)]
    
palette0 = {"2019-04-30":"cadetblue",
            "2019-09-03":"cadetblue",
            "2020-01-27":"peru", 
            "2020-06-19":"peru", 
            "2020-07-23":"peru",
            "2020-11-17":"peru",
            "2021-08-31":"gray",
            "2021-11-03":"gray",
            "2022-01-24":"gray",
            "Unassigned": "gray"
                }

palette4 = {"2019-04-30":"gray",
            "2019-09-03":"gray",
            "2020-01-27":"peru", 
            "2020-06-19":"peru", 
            "2020-07-23":"peru",
            "2020-11-17":"peru",
            "2021-08-31":"silver",
            "2021-11-03":"silver",
            "2022-01-24":"steelblue",
            "2022-05-25":"steelblue",
            "2022-08-22":"steelblue",
            "2022-10-11":"steelblue",
            "2022-12-30":"steelblue",
            "2023-01-21":"cadetblue",
            "2023-03-02":"cadetblue",
            "2023-06-30":"cadetblue",
            "2023-09-04":"cadetblue",
            "2023-11-25":"cadetblue",
            "Unassigned": "gray"
                }

palette5 = {2019:"gray",
            2020:"peri",
            2021:"silver", 
            2022:"cadetblue"
            }



palette0 = {"2019-04-30":"lightgray",
            "2019-09-03":"orangered",
            "2020-01-27":"firebrick", 
            "2020-06-19":"lightgray", 
            "2020-07-23":"peru",
            "2020-11-17":"lightgray",
            "2021-08-31":"lightgray",
            "2021-11-03":"lightgray",
            "2022-01-24":"lightgray",
            "Unassigned": "lightgray"
                }


palette2 = {"2019-01-23": "gray",
            "2019-02-27": "gray",
            "2019-04-02": "gray", 
            "2019-04-25": "gray",
            "2019-04-30": "gray",
            "2019-06-04": "gray",
            "2019-07-04": "gray",
            "2019-09-03": "gray",
            "2019-10-07": "gray", 
            "2019-11-11": "gray",
            "2019-12-30": "gray",     
            "2020-01-27": "peru",
            "2020-03-04": "peru",
            "2020-03-29": "peru",
            "2020-05-04": "peru",
            "2020-05-27": "peru",
            "2020-06-19": "peru", 
            "2020-07-23": "peru",
            "2020-09-01": "peru",
            "2020-10-08": "peru",
            "2020-11-17": "peru",
            "2020-12-14": "peru",
            "2021-01-11": "silver",
            "2021-02-04": "silver",
            "2021-03-01": "silver",
            "2021-04-08": "silver",
            "2021-05-26": "silver",
            "2021-06-30": "silver",
            "2021-08-03": "silver",  
            "2021-08-31": "silver",
            "2021-10-15": "silver",
            "2021-11-03": "silver",
            "2022-01-24": "steelblue",
            "2022-03-01": "steelblue",
            "2022-04-06": "steelblue",
            "2022-05-25": "steelblue",
            "2022-08-22": "steelblue",
            "2022-10-11": "steelblue",
            "2023-01-21": "cadetblue",
            "2023-03-02": "cadetblue",
            "2023-04-19": "cadetblue",
            "2023-05-31": "cadetblue",
            "2023-09-04": "cadetblue",
            }


palette2 = {"2019-01-23": "lightgray",
            "2019-02-27": "lightgray",
            "2019-04-02": "lightgray", 
            "2019-04-25": "lightgray",
            "2019-06-04": "lightgray",
            "2019-07-04": "lightgray",
            "2019-10-07": "lightgray", 
            "2019-11-11": "lightgray",
            "2019-12-30": "lightgray",     
            "2020-03-04": "lightgray",
            "2020-03-29": "lightgray",
            "2020-05-04": "lightgray",
            "2020-05-27": "lightgray",
            "2020-09-01": "lightgray",
            "2020-10-08": "lightgray",
            "2020-12-14": "lightgray",
            "2021-01-11": "steelblue",
            "2021-02-04": "lightgray",
            "2021-03-01": "lightgray",
            "2021-04-08": "lightgray",
            "2021-05-26": "lightgray",
            "2021-06-30": "lightgray",
            "2021-08-03": "lightgray",  
            "2021-10-15": "lightgray",
            "2019-04-30": "lightgray",
            "2019-09-03": "orangered",
            "2020-01-27": "firebrick", 
            "2020-06-19": "lightgray", 
            "2020-07-23": "peru",
            "2020-11-17": "lightgray",
            "2021-08-31": "lightgray",
            "2021-11-03": "lightgray",
            "2022-01-24": "lightgray",
            "Unassigned": "lightgray"
            }
           

palette1 = seaborn.color_palette("cubehelix", 8) #6
#palette1 = seaborn.color_palette("rocket", 1) #6
markers1 = [".", "."]

plots = 2
if plots == 1:
    #check extinction with plots before sending to Alex B.
    df = df.sort_values(by=['Launch', 'Altitude (km)', 'Measured Extinction_532'])
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    g = seaborn.scatterplot(x="Measured Extinction_532", y="Altitude (km)", data=df, hue="Launch", 
                              palette= palette0)
    # g.set(xscale='log', yscale='linear')
    # df = df.sort_values(by=['Launch',  'Altitude (km)', 'Measured Extinction_1064'])
    # fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    # g.set(xscale='log', yscale='linear')
    # g = seaborn.scatterplot(x="Measured Extinction_1064", y="Altitude (km)", data=df,  hue="Launch", 
    #                           palette= palette0)
    df = df.sort_values(by=['Launch', 'Altitude (km)', 'Ambient Extinction_532'])
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    g.set(xscale='log', yscale='linear')
    g = seaborn.scatterplot(x="Ambient Extinction_532", y="Altitude (km)", data=df,  hue="Launch", 
                              palette= palette0)
    g.set(xscale='log', yscale='linear')
    # df = df.sort_values(by=['Launch', 'Altitude (km)', 'Ambient Extinction_1064'])
    # fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    # g = seaborn.scatterplot(x="Ambient Extinction_1064", y="Altitude (km)", data=df, hue="Launch", 
    #                           palette= palette0)
    # g.set(xscale='log', yscale='linear')


 
if plots == 2:
    
    lauder_line = 'k-'
    lauder_shade = 'black'
    boulder_line = 'r-'
    boulder_shade = 'red'
    
    #effective radius
    Lpops_f = Lpops_f[Lpops_f['EffectiveRadius'] != np.nan]
    Lpops_f['EffectiveRadiusMean'] = Lpops_f['EffectiveRadius'].rolling(3, center=True).mean()
    Lpops_f = Lpops_f.sort_values(by=['Launch', 'Altitude','EffectiveRadius'])
    seaborn.set(font_scale=1, style="ticks")
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    
    g = seaborn.lineplot(x="EffectiveRadiusMean", y="Altitude", data=Lpops_f, orient='y', sort = False, hue="Launch", 
                            linestyle='dashed', palette = palette0, legend = False, zorder = 0)


    #g.set(xlim=[0.1, 1000])
    g.set(ylim=[0, 30])
    plt.xlabel(u"Effective Radius \u03bc$\mathregular{m}$) ") 
    plt.ylabel("Altitude (km) ") 
    g.set(xscale='log', yscale='linear')
    
    plt.axhline(y = 13.9, color = 'k', linestyle = ':')
    plt.axhline(y = 15.75, color = 'r', linestyle = ':')
    
    df_Bbase_ab = df_Bbase_ab.sort_values(by=['Launch', 'Altitude'])
    ax.plot(df_Bbase_ab['EF_base_q2'].rolling(3, center=True).mean(), df_Bbase_ab['Altitude'], boulder_line, linewidth = 3)
    ax.fill_betweenx(y=df_Bbase_ab['Altitude'], x1= df_Bbase_ab['EF_base_q1'], x2 = df_Bbase_ab['EF_base_q3'], color= boulder_shade, alpha=0.3)
    
    df_Lbase_ab = df_Lbase_ab.sort_values(by=['Launch', 'Altitude'])
    ax.plot(df_Lbase_ab['EF_base_q2'].rolling(3, center=True).mean(), df_Lbase_ab['Altitude'], lauder_line, linewidth = 3)
    ax.fill_betweenx(y=df_Lbase_ab['Altitude'], x1= df_Lbase_ab['EF_base_q1'], x2 = df_Lbase_ab['EF_base_q3'], color=lauder_shade, alpha=0.3)
    
    
    #aerosol concentration
    L_pc_pops_f = L_pc_pops_f[L_pc_pops_f['ParticleConcentration'] != np.nan]
    L_pc_pops_f = L_pc_pops_f.sort_values(by=['Launch', 'Altitude','ParticleConcentration'])
    seaborn.set(font_scale=1, style="ticks")
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    g = seaborn.lineplot(x="ParticleConcentration", y="Altitude", data=L_pc_pops_f, orient='y', sort = False, hue="Launch", 
                            linestyle='dashed', palette = palette0, legend = False, zorder = 0)
    
    #ax.legend(custom_lines, ['2019', '2020', '2021', '2022'])
    #g.set(xlim=[0.1, 1000])
    g.set(ylim=[0, 30])
    plt.xlabel(u"Particle Concentration (# $\mathregular{cm^{-3}}$)") #at S.T.P.
    plt.ylabel("Altitude (km) ") 
    g.set(xscale='log', yscale='linear')
    
    plt.axhline(y = 13.9, color = 'k', linestyle = ':')
    plt.axhline(y = 15.75, color = 'r', linestyle = ':')
    
    df_Bbase_pc_ab = df_Bbase_pc_ab.sort_values(by=['Launch', 'Altitude','PC_base_q2'])
    ax.plot(df_Bbase_pc_ab['PC_base_q2'], df_Bbase_pc_ab['Altitude'], boulder_line, linewidth = 3)
    ax.fill_betweenx(y=df_Bbase_pc_ab['Altitude'], x1= df_Bbase_pc_ab['PC_base_q1'], x2 = df_Bbase_pc_ab['PC_base_q3'], color=boulder_shade, alpha=0.3)
    
    df_Lbase_pc_ab = df_Lbase_pc_ab.sort_values(by=['Launch', 'Altitude','PC_base_q2'])
    ax.plot(df_Lbase_pc_ab['PC_base_q2'], df_Lbase_pc_ab['Altitude'], lauder_line, linewidth = 3)
    ax.fill_betweenx(y=df_Lbase_pc_ab['Altitude'], x1= df_Lbase_pc_ab['PC_base_q1'], x2 = df_Lbase_pc_ab['PC_base_q3'], color=lauder_shade, alpha=0.3)
    
    #aerosol concentration at STP
    L_pc_pops_f = L_pc_pops_f[L_pc_pops_f['ParticleConcentrationAtSTP'] != np.nan]
    L_pc_pops_f = L_pc_pops_f.sort_values(by=['Launch', 'Altitude','ParticleConcentrationAtSTP'])
    seaborn.set(font_scale=1, style="ticks")
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    g = seaborn.lineplot(x="ParticleConcentrationAtSTP", y="Altitude", data=L_pc_pops_f, orient='y', sort = False, hue="Launch", 
                            linestyle='dashed', palette = palette0, legend = True, zorder = 0)
    
    #ax.legend(custom_lines, ['2019', '2020', '2021', '2022'])
    #g.set(xlim=[0.1, 1000])
    g.set(ylim=[0, 30])
    plt.xlabel(u"Particle Concentration at S.T.P. (# $\mathregular{cm^{-3}}$)") #at S.T.P.
    plt.ylabel("Altitude (km) ") 
    g.set(xscale='log', yscale='linear')
    
    plt.axhline(y = 13.9, color = 'k', linestyle = ':')
    plt.axhline(y = 15.75, color = 'r', linestyle = ':')
    
    df_Bbase_pc_ab = df_Bbase_pc_ab.sort_values(by=['Launch', 'Altitude','PC_STP_base_q2'])
    ax.plot(df_Bbase_pc_ab['PC_STP_base_q2'], df_Bbase_pc_ab['Altitude'], boulder_line, linewidth = 3)
    ax.fill_betweenx(y=df_Bbase_pc_ab['Altitude'], x1= df_Bbase_pc_ab['PC_STP_base_q1'], x2 = df_Bbase_pc_ab['PC_STP_base_q3'], color=boulder_shade, alpha=0.3)
    
    df_Lbase_pc_ab = df_Lbase_pc_ab.sort_values(by=['Launch', 'Altitude','PC_STP_base_q2'])
    ax.plot(df_Lbase_pc_ab['PC_STP_base_q2'], df_Lbase_pc_ab['Altitude'], lauder_line, linewidth = 3)
    ax.fill_betweenx(y=df_Lbase_pc_ab['Altitude'], x1= df_Lbase_pc_ab['PC_STP_base_q1'], x2 = df_Lbase_pc_ab['PC_STP_base_q3'], color=lauder_shade, alpha=0.3)
    
    
    #surface area
    L_sa_pops_f = L_sa_pops_f[L_sa_pops_f['SurfaceArea'] != np.nan]
    L_sa_pops_f = L_sa_pops_f.sort_values(by=['Launch', 'Altitude','SurfaceArea'])
    seaborn.set(font_scale=1, style="ticks")
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    g = seaborn.lineplot(x="SurfaceArea", y="Altitude", data=L_sa_pops_f, orient='y', sort = False, hue="Launch", 
                            linestyle='dashed', legend = False, palette = palette0, zorder = 0)

    #g.set(xlim=[0.1, 1000])
    g.set(ylim=[0, 30])
    plt.xlabel(u"Particle Surface Area (\u03bc$\mathregular{m^{2}}$ $\mathregular{cm^{-3}}$)")
    plt.ylabel(" ") 
    g.set(xscale='log', yscale='linear')
    
    plt.axhline(y = 13.9, color = 'k', linestyle = ':')
    plt.axhline(y = 15.75, color = 'r', linestyle = ':')
    
    df_Bbase_sa_ab = df_Bbase_sa_ab.sort_values(by=['Launch', 'Altitude','SA_base_q2'])
    ax.plot(df_Bbase_sa_ab['SA_base_q2'], df_Bbase_sa_ab['Altitude'], boulder_line, linewidth = 3)
    ax.fill_betweenx(y=df_Bbase_sa_ab['Altitude'], x1= df_Bbase_sa_ab['SA_base_q1'], x2 = df_Bbase_sa_ab['SA_base_q3'], color=boulder_shade, alpha=0.3)
    
    df_Lbase_sa_ab = df_Lbase_sa_ab.sort_values(by=['Launch', 'Altitude','SA_base_q2'])
    ax.plot(df_Lbase_sa_ab['SA_base_q2'], df_Lbase_sa_ab['Altitude'], lauder_line, linewidth = 3)
    ax.fill_betweenx(y=df_Lbase_sa_ab['Altitude'], x1= df_Lbase_sa_ab['SA_base_q1'], x2 = df_Lbase_sa_ab['SA_base_q3'], color=lauder_shade, alpha=0.3)
    
    
    
    #surface area at STP
    L_sa_pops_f = L_sa_pops_f[L_sa_pops_f['SurfaceAreaAtSTP'] != np.nan]
    L_sa_pops_f = L_sa_pops_f.sort_values(by=['Launch', 'Altitude','SurfaceAreaAtSTP'])
    seaborn.set(font_scale=1, style="ticks")
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    g = seaborn.lineplot(x="SurfaceAreaAtSTP", y="Altitude", data=L_sa_pops_f, orient='y', sort = False, hue="Launch", 
                            linestyle='dashed', legend = False, palette = palette0, zorder = 0)

    #g.set(xlim=[0.1, 1000])
    g.set(ylim=[0, 30])
    plt.xlabel(u"Particle Surface Area at S.T.P. (\u03bc$\mathregular{m^{2}}$ $\mathregular{cm^{-3}}$ at S.T.P.)")
    plt.ylabel("") 
    g.set(xscale='log', yscale='linear')
    
    plt.axhline(y = 13.9, color = 'k', linestyle = ':')
    plt.axhline(y = 15.75, color = 'r', linestyle = ':')
    
    df_Bbase_sa_ab = df_Bbase_sa_ab.sort_values(by=['Launch', 'Altitude','SA_STP_base_q2'])
    ax.plot(df_Bbase_sa_ab['SA_STP_base_q2'], df_Bbase_sa_ab['Altitude'], boulder_line, linewidth = 3)
    ax.fill_betweenx(y=df_Bbase_sa_ab['Altitude'], x1= df_Bbase_sa_ab['SA_STP_base_q1'], x2 = df_Bbase_sa_ab['SA_STP_base_q3'], color=boulder_shade, alpha=0.3)
    
    df_Lbase_sa_ab = df_Lbase_sa_ab.sort_values(by=['Launch', 'Altitude','SA_STP_base_q2'])
    ax.plot(df_Lbase_sa_ab['SA_STP_base_q2'], df_Lbase_sa_ab['Altitude'], lauder_line, linewidth = 3)
    ax.fill_betweenx(y=df_Lbase_sa_ab['Altitude'], x1= df_Lbase_sa_ab['SA_STP_base_q1'], x2 = df_Lbase_sa_ab['SA_STP_base_q3'], color=lauder_shade, alpha=0.3)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #Figure 1 or 2 in pyroCB paper
    #AEROSOL CONCENTRATION AMBIENT
    df_AC = df[df['Aconc_q2'] != np.nan]
    df_AC = df[df['Aconc_q2'] > 0.3]
    #seaborn.set(style="whitegrid")
    seaborn.set(font_scale=1, style="ticks")
    
    df_AC = df
    df_AC = df_AC.sort_values(by=['Launch', 'Altitude (km)','Aconc_q2'])
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    #plotorder = df5.sort_values(by=['Launch','Flight','Altitude (km)'])
    #seaborn.set(style="darkgrid")
    g = seaborn.lineplot(x="Aconc_q2", y="Altitude (km)", data=df_AC, orient='y', sort = False, hue="Launch", 
                            linestyle='dashed', palette= palette0,  zorder = 0, legend=False)

    g.set(xlim=[0.1, 1000])
    g.set(ylim=[0, 30])
    plt.xlabel(u"Particle Concentration (# $\mathregular{cm^{-3}}$)") #at S.T.P.
    g.set(xscale='log', yscale='linear')
    #seaborn.set(rc={"xtick.bottom" : True, "ytick.left" : True})
    
    # specifying horizontal line type
    plt.axhline(y = 13.9, color = 'k', linestyle = ':')
    plt.axhline(y = 15.75, color = 'r', linestyle = ':')
    
    # #plot Boulder median +/- IQR
    df_Bbase_pc_ab = df_Bbase_pc_ab.sort_values(by=['Launch', 'Altitude','PC_base_q2'])
    ax.plot(df_Bbase_pc_ab['PC_base_q2'].rolling(3, center=True).mean(), df_Bbase_pc_ab['Altitude'], boulder_line, linewidth = 3)
    ax.fill_betweenx(y=df_Bbase_pc_ab['Altitude'], x1= df_Bbase_pc_ab['PC_base_q1'], x2 = df_Bbase_pc_ab['PC_base_q3'], color=boulder_shade, alpha=0.3)
    
    #ax.plot(df_Boulder['50th_cm3'], df_Boulder['Altitude (km)'], 'k-', linewidth = 3)
    #ax.fill_betweenx(y=df_Boulder['Altitude (km)'], x1=df_Boulder['25th'], x2 = df_Boulder['75th'], color="black", alpha=0.3)
    
    #plot median +/- IQR
    df_base1 = df_base.groupby(['Altitude (km)' ], as_index=False).mean()
    df_base1= df_base1[df_base1['Altitude (km)'] > 0.3]
    df_base1= df_base1[df_base1['Altitude (km)'] != 0]
    df_base1 = df_base1.sort_values(by=['Altitude (km)'])
    ax.plot(df_base1['AC_base_q2'].rolling(3, center=True).mean(), df_base1['Altitude (km)'], lauder_line, linewidth = 3)
    ax.fill_betweenx(y=df_base1['Altitude (km)'], x1=df_base1['AC_base_q1'], x2 = df_base1['AC_base_q3'], color=lauder_shade, alpha=0.3)
    
    
    #AEROSOL CONCENTRATION AT STP
    df_AC_STP = df
    df_AC_STP = df_AC_STP.sort_values(by=['Launch', 'Altitude (km)','AconcSTP_q2'])
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    #plotorder = df5.sort_values(by=['Launch','Flight','Altitude (km)'])
    #seaborn.set(style="darkgrid")
    g = seaborn.lineplot(x="AconcSTP_q2", y="Altitude (km)", data=df_AC_STP, orient='y', sort = False, hue="Launch", 
                            linestyle='dashed', palette= palette0,  zorder = 0, legend = False)

    g.set(xlim=[0.1, 1000])
    g.set(ylim=[0, 30])
    plt.xlabel(u"Particle Concentration (# $\mathregular{cm^{-3}}$) at S.T.P.") #at S.T.P.
    plt.ylabel('')
    g.set(xscale='log', yscale='linear')
    #seaborn.set(rc={"xtick.bottom" : True, "ytick.left" : True})
    
    # specifying horizontal line type
    plt.axhline(y = 13.9, color = 'k', linestyle = ':')
    plt.axhline(y = 15.75, color = 'r', linestyle = ':')
    
    # #plot Boulder median +/- IQR
    df_Bbase_pc_ab = df_Bbase_pc_ab.sort_values(by=['Launch', 'Altitude','PC_STP_base_q2'])
    ax.plot(df_Bbase_pc_ab['PC_STP_base_q2'].rolling(3, center=True).mean(), df_Bbase_pc_ab['Altitude'], boulder_line, linewidth = 3)
    ax.fill_betweenx(y=df_Bbase_pc_ab['Altitude'], x1= df_Bbase_pc_ab['PC_STP_base_q1'], x2 = df_Bbase_pc_ab['PC_STP_base_q3'], color=boulder_shade, alpha=0.3)
    
    #ax.plot(df_Boulder['50th_cm3'], df_Boulder['Altitude (km)'], 'k-', linewidth = 3)
    #ax.fill_betweenx(y=df_Boulder['Altitude (km)'], x1=df_Boulder['25th'], x2 = df_Boulder['75th'], color="black", alpha=0.3)
    
    #plot median +/- IQR
    df_base1 = df_base1.sort_values(by=['Altitude (km)'])
    ax.plot(df_base1['ACSTP_base_q2'].rolling(3, center=True).mean(), df_base1['Altitude (km)'], lauder_line, linewidth = 3)
    ax.fill_betweenx(y=df_base1['Altitude (km)'], x1=df_base1['ACSTP_base_q1'], x2 = df_base1['ACSTP_base_q3'], color=lauder_shade, alpha=0.3)
    
    
    #EFFECTIVE RADIUS
    df_ER = df
    df_ER[df_ER['Erad_q2']==0]=np.NaN
    df_ER = df_ER.sort_values(by=['Launch', 'Altitude (km)','Erad_q2'])
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    g = seaborn.lineplot(x="Erad_q2", y="Altitude (km)", data=df_ER, orient='y', sort = False, hue="Launch", 
                              palette= palette0, linestyle='dashed', zorder = 0, legend = False) #
    #g.set(xlim=[0.1, 1000])
    g.set(ylim=[0, 30])
    plt.xlabel(u"Effective Radius (nm)")
    plt.ylabel('Altitude (km)')
    g.set(xscale='log', yscale='linear')
    #seaborn.set(rc={"xtick.bottom" : True, "ytick.left" : True})
    
    # specifying horizontal line type
    plt.axhline(y = 13.9, color = 'k', linestyle = ':')
    plt.axhline(y = 15.75, color = 'r', linestyle = ':')
    
    # #plot Boulder median +/- IQR
    df_Bbase_ab = df_Bbase_ab.sort_values(by=['Launch', 'Altitude'])
    ax.plot(df_Bbase_ab['EF_base_q2'].rolling(3, center=True).mean()*1000, df_Bbase_ab['Altitude'], boulder_line, linewidth = 3)
    ax.fill_betweenx(y=df_Bbase_ab['Altitude'], x1= df_Bbase_ab['EF_base_q1']*1000, x2 = df_Bbase_ab['EF_base_q3']*1000, color= boulder_shade, alpha=0.3)
    
    #plot median +/- IQR
    df_base1 = df_base1.sort_values(by=['Altitude (km)'])
    ax.plot(df_base1['ER_base_q2'].rolling(3, center=True).mean(), df_base1['Altitude (km)'], lauder_line, linewidth = 3)
    ax.fill_betweenx(y=df_base1['Altitude (km)'], x1=df_base1['ER_base_q1'], x2 = df_base1['ER_base_q3'], color=lauder_shade, alpha=0.3)
    
    #SURFACE AREA AT STP
    df_SA = df
    df_SA[df_SA['Asa_q2'] == 0] = np.NaN
    df_SA = df_AC.sort_values(by=['Launch', 'Altitude (km)','Asa_q2'])
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    #plotorder = df5.sort_values(by=['Launch','Flight','Altitude (km)'])
    #seaborn.set(style="darkgrid")
    g = seaborn.lineplot(x="Asa_q2", y="Altitude (km)", data=df_SA, orient='y', sort = False, hue="Launch", 
                              palette= palette0, linestyle='dashed', zorder = 0, legend = False) #

    g.set(xlim=[0.1, 1000])
    g.set(ylim=[0, 30])
    plt.xlabel(u"Particle Surface Area (\u03bc$\mathregular{m^{2}}$ $\mathregular{cm^{-3}}$ at S.T.P.)")
    plt.ylabel('')
    g.set(xscale='log', yscale='linear')
    #seaborn.set(rc={"xtick.bottom" : True, "ytick.left" : True})
    
    # specifying horizontal line type
    plt.axhline(y = 13.9, color = 'k', linestyle = ':')
    plt.axhline(y = 15.75, color = 'r', linestyle = ':')
    
    # #plot Boulder median +/- IQR
    df_Bbase_sa_ab = df_Bbase_sa_ab.sort_values(by=['Launch', 'Altitude','SA_STP_base_q2'])
    ax.plot(df_Bbase_sa_ab['SA_STP_base_q2'].rolling(3, center=True).mean(), df_Bbase_sa_ab['Altitude'], boulder_line, linewidth = 3)
    ax.fill_betweenx(y=df_Bbase_sa_ab['Altitude'], x1= df_Bbase_sa_ab['SA_STP_base_q1'], x2 = df_Bbase_sa_ab['SA_STP_base_q3'], color=boulder_shade, alpha=0.3)
    
    #plot median +/- IQR
    ax.plot(df_base1['ASA_base_q2'].rolling(3, center=True).mean(), df_base1['Altitude (km)'], lauder_line, linewidth = 3)
    ax.fill_betweenx(y=df_base1['Altitude (km)'], x1=df_base1['ASA_base_q1'], x2 = df_base1['ASA_base_q3'], color=lauder_shade, alpha=0.3)
    
    
    
    df_Boulder['Aer_Decr'] = (df_Boulder['50th_cm3'].shift(-1) - df_Boulder['50th_cm3'])/0.25 
    df_base['Aer_Decr'] = (df_base['AC_base_q2'].shift(-1) - df_base['AC_base_q2'])/0.25     
    
    
    
    
    
    #remove nan rows so it isn't difficult to plot a line vertical profile...
    df_w = df_w[df_w['Flight'] == 'Ascent']
    df_w = df_w[df_w['Altitude (km)'] <= 28.0]
    df_w = df_w[df_w['Altitude (km)'] > 0.0]
    df_H2O = df_w.filter(['Launch', 'Altitude (km)', 'H2O_q2'])
    #df_H2O = df_H2O[df_H2O['H2O_q2'] != np.nan]
    df_H2O = df_H2O.groupby(['Launch','Altitude (km)' ], as_index=False).mean()
    df_H2O = df_H2O.sort_values(by=['Altitude (km)', 'Launch' ])
    df_H2O = df_H2O.reset_index(drop=True)
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    # g = seaborn.lineplot(data=df_H2O, x="H2O_q2", y="Altitude (km)", orient='y', sort = False, palette= palette2, linestyle='dashed', legend = False, hue="Launch", hue_order = ["2019-01-23", "2019-02-27", "2019-04-02", "2019-04-25",  "2019-06-04", "2019-07-04", 
    #             "2019-10-07", "2019-11-11", "2019-12-30","2020-03-04","2020-03-29","2020-05-04", "2020-05-27",  "2020-09-01", "2020-10-08",  "2020-12-14", "2021-01-11", "2021-02-04",
    #             "2021-03-01", "2021-04-08", "2021-05-26", "2021-06-30", "2021-08-03", "2021-10-15","2019-04-30", "2019-09-03", "2020-01-27", "2020-07-23", "2020-11-17", "2021-08-31", "2021-11-03", "2022-01-24"],  zorder = 0) #
 
    g = seaborn.lineplot(data=df_H2O, x="H2O_q2", y="Altitude (km)", orient='y', sort = False, palette= palette2, linestyle='dashed', legend = False, hue="Launch", zorder = 0) #
 
    plt.xlabel("H2O (ppmv)")
    plt.ylabel('')
    g.set(xscale='log', yscale='linear')
    g.set(xlim=[2, 2000])
    g.set(ylim=[0, 30])

    
    # specifying horizontal line type
    plt.axhline(y = 13.9, color = 'k', linestyle = ':')
    plt.axhline(y = 15.75, color = 'r', linestyle = ':')
    
    #plot Boulder water IQR +/- median
    df.rolling(3, center=True).mean()

    ax.plot(df_wv_Boulder['H2O_q2'], df_wv_Boulder['Altitude (km)'], boulder_line, linewidth = 3)
    ax.fill_betweenx(y=df_wv_Boulder['Altitude (km)'], x1=df_wv_Boulder['H2O_q1'], x2 = df_wv_Boulder['H2O_q3'], color=boulder_shade, alpha=0.3)
    
    ax.plot( df_base1['H2O_base_q2'], df_base1['Altitude (km)'], lauder_line, linewidth = 3)
    ax.fill_betweenx(y=df_base1['Altitude (km)'], x1=df_base1['H2O_base_q1'], x2 = df_base1['H2O_base_q3'], color=lauder_shade, alpha=0.3)
    
    
    #tropopause relevant coordinates
    df_w_TA2 = df_w_TA[df_w_TA['Flight'] == 'Ascent']
    #df_w_TA = df_w_TA[df_w_TA['Tropopause_Altitude (km)'] <= 28.0]
    #df_w_TA = df_w_TA[df_w_TA['Tropopause_Altitude (km)'] > 0.0]
    
    df_H2O = df_w_TA2.filter(['Launch', 'Tropopause_Altitude (km)', 'H2O_q2_TA'])
    df_H2O = df_H2O.groupby(['Launch','Tropopause_Altitude (km)' ], as_index=False).mean()
    df_H2O = df_H2O.sort_values(by=['Tropopause_Altitude (km)', 'Launch' ])
    df_H2O = df_H2O.reset_index(drop=True)
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    # g = seaborn.lineplot(data=df_H2O, x="H2O_q2", y="Altitude (km)", orient='y', sort = False, palette= palette2, linestyle='dashed', legend = False, hue="Launch", hue_order = ["2019-01-23", "2019-02-27", "2019-04-02", "2019-04-25",  "2019-06-04", "2019-07-04", 
    #             "2019-10-07", "2019-11-11", "2019-12-30","2020-03-04","2020-03-29","2020-05-04", "2020-05-27",  "2020-09-01", "2020-10-08",  "2020-12-14", "2021-01-11", "2021-02-04",
    #             "2021-03-01", "2021-04-08", "2021-05-26", "2021-06-30", "2021-08-03", "2021-10-15","2019-04-30", "2019-09-03", "2020-01-27", "2020-07-23", "2020-11-17", "2021-08-31", "2021-11-03", "2022-01-24"],  zorder = 0) #
 
    g = seaborn.lineplot(data=df_H2O, x="H2O_q2_TA", y="Tropopause_Altitude (km)", orient='y', sort = False, palette= palette2, linestyle='dashed', legend = False, hue="Launch", zorder = 0) #
 
    plt.xlabel("H2O (ppmv)")
    plt.ylabel('Tropopause Relevant Altitude (km)')
    g.set(xscale='linear', yscale='linear')
    g.set(xlim=[2, 10])
    g.set(ylim=[-2, 20])
    
    #plot Boulder water IQR +/- median
        
    ax.plot(df_base2['H2O_q2_TA'], df_base2['Tropopause_Altitude (km)'], boulder_line, linewidth = 3)
    ax.fill_betweenx(y=df_base2['Tropopause_Altitude (km)'], x1=df_base2['H2O_q1_TA'], x2 = df_base2['H2O_q3_TA'], color=boulder_shade, alpha=0.3)
    
    ax.plot(df_wv_Boulder2['H2O_q2_TA'], df_wv_Boulder2['Tropopause_Altitude (km)'], lauder_line, linewidth = 3)
    ax.fill_betweenx(y=df_wv_Boulder2['Tropopause_Altitude (km)'], x1=df_wv_Boulder2['H2O_q1_TA'], x2 = df_wv_Boulder2['H2O_q3_TA'], color=lauder_shade, alpha=0.3)

    
    plt.xlabel("H2O (ppmv)")
    plt.ylabel('Tropopause Relative Altitude (km)')
    
    #df_O3 = df_w[df_w['O3_q2'] != np.nan]
    df_O3 = df_w.filter(['Launch', 'Altitude (km)', 'O3_q2'])
    df_O3 = df_O3[df_O3['O3_q2'] > 0.0]
    df_O3 = df_O3[df_O3['Altitude (km)'] >= 0.4]
    df_O3 = df_O3.groupby(['Launch','Altitude (km)' ], as_index=False).mean()
    df_O3 = df_O3.sort_values(by=['Altitude (km)','Launch'])
    df_O3 = df_O3.reset_index(drop=True)
    #seaborn.set(rc={"xtick.bottom" : True, "ytick.left" : True})
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    
    
    # g = seaborn.lineplot(x="O3_q2", y="Altitude (km)", data=df_O3, hue="Launch", hue_order = ["2019-01-23", "2019-02-27", "2019-04-02", "2019-04-25",  "2019-06-04", "2019-07-04", 
    #             "2019-10-07", "2019-11-11", "2019-12-30","2020-03-04","2020-03-29","2020-05-04", "2020-05-27",  "2020-09-01", "2020-10-08",  "2020-12-14", "2021-01-11", "2021-02-04",
    #             "2021-03-01", "2021-04-08", "2021-05-26", "2021-06-30", "2021-08-03", "2021-10-15","2019-04-30", "2019-09-03", "2020-01-27", "2020-07-23", "2020-11-17", "2021-08-31", "2021-11-03", "2022-01-24"], orient='y', sort = False,
    #                       palette= palette2, legend = False, linestyle='dashed', zorder = 0) 

    g = seaborn.lineplot(x="O3_q2", y="Altitude (km)", data=df_O3, hue="Launch", orient='y', sort = False,
                      palette= palette2, legend = False, linestyle='dashed', zorder = 0) 

    plt.xlabel("O3 (ppmv)")
    plt.ylabel('Altitude (km)')
    g.set(ylim=[0, 30])
    g.set(xlim=[0.001, 10])
    g.set(xscale='log', yscale='linear')
    #seaborn.set(rc={"xtick.bottom" : True, "ytick.left" : True})
    
    # specifying horizontal line type
    plt.axhline(y = 13.9, color = 'k', linestyle = ':')
    plt.axhline(y = 15.75, color = 'r', linestyle = ':')
    
    #plot Boulder water IQR +/- median
    df_base = df_base[df_base['Altitude (km)'] >= 0.4]
    ax.plot(df_o3_Boulder['O3_q2'], df_o3_Boulder['Altitude (km)'], boulder_line, linewidth = 3)
    ax.fill_betweenx(y=df_o3_Boulder['Altitude (km)'], x1=df_o3_Boulder['O3_q1'], x2 = df_o3_Boulder['O3_q3'], color=boulder_shade, alpha=0.3)
    
    ax.plot(df_base1['O3_base_q2'], df_base1['Altitude (km)'], lauder_line, linewidth = 3)
    ax.fill_betweenx(y=df_base1['Altitude (km)'], x1=df_base1['O3_base_q1'], x2 = df_base1['O3_base_q3'], color=lauder_shade, alpha=0.3)
    
    
        
    #tropopause relevant coordinates
    
    df_O3 = df_w_TA2.filter(['Launch', 'Tropopause_Altitude (km)', 'O3_q2_TA'])
    df_O3 = df_O3.groupby(['Launch','Tropopause_Altitude (km)' ], as_index=False).mean()
    df_O3 = df_O3.sort_values(by=['Tropopause_Altitude (km)', 'Launch' ])
    df_O3 = df_O3.reset_index(drop=True)
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    # g = seaborn.lineplot(data=df_H2O, x="H2O_q2", y="Altitude (km)", orient='y', sort = False, palette= palette2, linestyle='dashed', legend = False, hue="Launch", hue_order = ["2019-01-23", "2019-02-27", "2019-04-02", "2019-04-25",  "2019-06-04", "2019-07-04", 
    #             "2019-10-07", "2019-11-11", "2019-12-30","2020-03-04","2020-03-29","2020-05-04", "2020-05-27",  "2020-09-01", "2020-10-08",  "2020-12-14", "2021-01-11", "2021-02-04",
    #             "2021-03-01", "2021-04-08", "2021-05-26", "2021-06-30", "2021-08-03", "2021-10-15","2019-04-30", "2019-09-03", "2020-01-27", "2020-07-23", "2020-11-17", "2021-08-31", "2021-11-03", "2022-01-24"],  zorder = 0) #
 
    g = seaborn.lineplot(data=df_O3, x="O3_q2_TA", y="Tropopause_Altitude (km)", orient='y', sort = False, palette= palette2, linestyle='dashed', legend = False, hue="Launch", zorder = 0) #
 
    plt.xlabel("O3 (ppmv)")
    plt.ylabel('Tropopause Relative Altitude (km)')
    g.set(xscale='linear', yscale='linear')
    g.set(xlim=[0, 12])
    g.set(ylim=[-2, 20])
    
    #plot Boulder water IQR +/- median
        
    ax.plot(df_base2['O3_q2_TA'], df_base2['Tropopause_Altitude (km)'], boulder_line, linewidth = 3)
    ax.fill_betweenx(y=df_base2['Tropopause_Altitude (km)'], x1=df_base2['O3_q1_TA'], x2 = df_base2['O3_q3_TA'], color=boulder_shade, alpha=0.3)
    
    ax.plot(df_o3_Boulder2['O3_q2_TA'], df_o3_Boulder2['Tropopause_Altitude (km)'], lauder_line, linewidth = 3)
    ax.fill_betweenx(y=df_o3_Boulder2['Tropopause_Altitude (km)'], x1=df_o3_Boulder2['O3_q1_TA'], x2 = df_o3_Boulder2['O3_q3_TA'], color=lauder_shade, alpha=0.3)


#OTHER FIGUERS. UNCOMMENT ALL NEXT LINES TO PLOT MORE THAN FIGURE 2  
    
    
 #    #add a figure for temperature - does the WMO tropopause seem reasonable?
 #    #correct Boulder to show K.,,
 #    df_T = df_w[df_w['T_q2'] != np.nan]
 #    df_T = df_T.sort_values(by=['Launch','Altitude (km)'])
 #    #seaborn.set(rc={"xtick.bottom" : True, "ytick.left" : True})
 #    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
 #    g = seaborn.lineplot(x="T_q2", y="Altitude (km)", data=df_T, hue="Launch", sort = False,
 #                          palette= palette2, linestyle='dashed') #

 #    plt.xlabel("Temperatue (C)")
 #    plt.ylabel('')
 #    g.set(ylim=[0, 30])
 #    g.set(xscale='linear', yscale='linear')
 #    #seaborn.set(rc={"xtick.bottom" : True, "ytick.left" : True})
    
 #    # specifying horizontal line type
 #    plt.axhline(y = 13.9, color = 'r', linestyle = ':')
 #    plt.axhline(y = 15.75, color = 'k', linestyle = ':')
    
 #    #plot Boulder Temperature IQR +/- median
 #    ax.plot(df_t_Boulder['T_q2'], df_t_Boulder['Altitude (km)'], 'k-', linewidth = 3)
 #    ax.fill_betweenx(y=df_t_Boulder['Altitude (km)'], x1=df_t_Boulder['T_q1'], x2 = df_t_Boulder['T_q3'], color="black", alpha=0.3)
    
 #    g.plot(df_base['T_base_q2'], df_base['Altitude (km)'], 'r-', linewidth = 3)
 #    g.fill_betweenx(y=df_base['Altitude (km)'], x1=df_base['T_base_q1'], x2 = df_base['T_base_q3'], color="red", alpha=0.3)
    
   
 #    #Figure 3 tracking the plume
 #    fig = plt.figure()
 #    ax1 = fig.add_subplot(111)
 #    ax2 = ax1.twiny()
    
 #    df_AC = df_AC[df_AC.Launch == '2020-01-27']
 #    #remove nan rows so it isn't difficult to plot a line vertical profile...
 #    df_tp = df_AC.groupby(['Launch'],as_index=False).median()
 #    tropopause1_plot = df_tp['Tropopause1'].values
 #    tropopause2_plot = df_tp['Tropopause2'].values
    
 #    #ax1.plot(df_base['AV_q2'], df_base['Altitude (km)'], '-', linewidth = 1, color = 'darkred')
 #    ax1.plot(df_AC['AV_q2'], df_AC['Altitude (km)'], 'r-', linewidth = 1)
 #    ax1.fill_betweenx(y= df_AC['Altitude (km)'], x1= df_AC['AV_q1'], x2 = df_AC['AV_q3'], color="red", alpha=0.3)
 #    ax1.set_xlabel(u"Aerosol Volume (\u03bc$\mathregular{m^{3}}$ $\mathregular{cm^{-3}}$)")
 #    ax1.set(ylim=[10, 30])
 #    ax1.set_xscale('log')
 #    ax1.axhline(y = tropopause1_plot , color = 'k', linestyle = '--')
 #    ax1.axhline(y = tropopause2_plot , color = 'k', linestyle = '--')
 #    ax1.xaxis.label.set_color('darkred') #0.11 median vs 5.8
    
 #    df_H2O = df_w[df_w.Launch == '2020-01-27']
 #    #ax2.plot( df_base['H2O_base_q2'], df_base['Altitude (km)'], '-', linewidth = 1, color = 'darkblue')
 #    ax2.plot(df_H2O['H2O_q2'], df_H2O['Altitude (km)'], 'b-', linewidth = 1)
 #    ax2.fill_betweenx(y=df_H2O['Altitude (km)'], x1=df_H2O['H2O_q1'], x2 = df_H2O['H2O_q3'], color="blue", alpha=0.3)
 #    ax2.set_xlabel("H2O (ppmv)")
 #    ax2.set(xlim=[3, 10])
 #    ax2.set(ylim=[10, 30])
 #    ax2.xaxis.label.set_color('darkblue') #3.7 median vs 4.6

    
 #    fig = plt.figure()
 #    ax1 = fig.add_subplot(111)
 #    ax2 = ax1.twiny()
    
 #    #df_AC = df_AC[df_AC.Launch == '2020-01-27']
 #    df_O3 = df_w
 #    df_O3 = df_O3[df_O3.Launch == '2020-01-27']
 #    ax1.plot(df_O3['O3_q2'], df_O3['Altitude (km)'], 'r-', linewidth = 1)
 #    ax1.fill_betweenx(y= df_O3['Altitude (km)'], x1= df_O3['O3_q1'], x2 = df_O3['O3_q3'], color="red", alpha=0.3)
 #    ax1.axhline(y = tropopause1_plot , color = 'k', linestyle = '--')
 #    ax1.axhline(y = tropopause2_plot , color = 'k', linestyle = '--')
 #    ax1.set_xlabel("O3 (ppmv)")
 #    ax1.set(ylim=[10, 30])
 #    ax1.set_xscale('log')
 #    ax1.set_ylabel('Altitude (km) ')
 #    ax1.xaxis.label.set_color('darkred')
    
 #    ax2.plot(df_AC['T_q2'], df_AC['Altitude (km)'], 'b-', linewidth = 1)
 #    ax2.fill_betweenx(y= df_AC['Altitude (km)'], x1= df_AC['T_q1'], x2 = df_AC['T_q3'], color="blue", alpha=0.3)
 #    ax2.set_xlabel(u"Air Temperature (K)")
 # #  ax2.set(ylim=[10, 30])
 #    ax2.set(xlim=[200, 240])
 #    #ax2.set(xlim=[-70, 30])
 #    ax2.set_ylabel('Altitude (km)')
 #    ax2.xaxis.label.set_color('darkblue')
    
    
 #    #figure 4 (timeseries
    df3  = df2.filter([ 'Launch', 'Altitude (km)', 'Tropopause', 'Theta [K]', 'Aer_Conc' ,  'Ambient Extinction_532']) #vs measured extinction....
    df3['AtmosLayer'] = df3.apply(lambda x: StratTrop(x['Theta [K]'], x['Altitude (km)'], x['Tropopause']), axis=1)
    df3 = df3 [df3['AtmosLayer'] != 'Troposphere']   
    df3['Number Concentration Column(# m^-2)'] = df3['Aer_Conc']*1E6*100# this is number concentration (#/m2) which can be summed to give a column
    df3['sAOD']  = df3['Ambient Extinction_532']*0.1 #units already K-1 which can be summed...
    df3.drop(['Altitude (km)', 'Tropopause', 'Theta [K]', 'Aer_Conc', 'Ambient Extinction_532'], axis = 1, inplace = True)
    df3_s = df3.groupby(['Launch', 'AtmosLayer'], as_index=False).sum()
    Column =['Number Concentration Column(# m^-2)',  'sAOD']
    df_int = pd.melt(df3_s, id_vars=['Launch', 'AtmosLayer'], value_vars= Column, var_name = 'Property', value_name = 'Vertically Integrated Total')
    df_int['Date'] = pd.to_datetime(df_int['Launch']).apply(lambda date: date.toordinal())
    df_int['Source'] = 'POPS'
    
    
    df4  = df2.filter([ 'Launch', 'Altitude (km)', 'Tropopause', 'Theta [K]', 'Aer_Conc_205nm',  'Measured Extinction_532']) #vs measured extinction....
    df4['AtmosLayer'] = df4.apply(lambda x: StratTrop(x['Theta [K]'], x['Altitude (km)'], x['Tropopause']), axis=1)
    df4 = df4 [df4['AtmosLayer'] != 'Troposphere']   
    df4['Number Concentration Column(# m^-2)'] = df4['Aer_Conc_205nm']*0.1# this is number concentration (#/m2) which can be summed to give a column
    df4['sAOD']  = df4['Measured Extinction_532']*100/1000 #units already K-1 which can be summed...
    df4.drop(['Altitude (km)', 'Tropopause', 'Theta [K]', 'Aer_Conc_205nm',  'Measured Extinction_532'], axis = 1, inplace = True)
    df4_s = df4.groupby(['Launch', 'AtmosLayer'], as_index=False).sum()
    Column4 =['Number Concentration Column(# m^-2)', 'sAOD']
    df_int4 = pd.melt(df4_s, id_vars=['Launch', 'AtmosLayer'], value_vars= Column4, var_name = 'Property', value_name = 'Vertically Integrated Total')
    df_int4['Date'] = pd.to_datetime(df_int4['Launch']).apply(lambda date: date.toordinal())
    df_int4['Source'] = 'POPS > 205 nm'
    
 #    #add a column for uncertainty 
    df_int['Unc'] = df_int.apply(lambda x: calcTsUnc(x['Launch'], x['Property'], x['Vertically Integrated Total']), axis=1)
    dfc_int['Unc'] = np.nan
    dfc_int2['Unc'] = np.nan
    dfl_sf['Unc'] = dfl_sf_u['Vertically Integrated Total']
    
 #    # cont = 0
 #    # if cont == 1: 
    df_int_all = pd.concat([df_int, dfc_int, dfc_int2, dfl_sf], ignore_index=True, axis=0)
    ord_labels = [ '2018-11-2', '2019-05-21', '2019-12-7', '2020-06-24', '2021-01-10', '2021-07-29', '2022-02-14', '2022-09-02'] # 
    
    
 #    #timeseries for shaded error region
    df_int_nc_ls = df_int[df_int['Property'] == 'Number Concentration Column(# m^-2)']
    df_int_nc_ms = df_int_nc_ls[df_int_nc_ls['AtmosLayer'] == 'Stratosphere']
    
    dfc_int_nc_ls = dfc_int[dfc_int['Property'] == 'Number Concentration Column(# m^-2)']
    
     #df_int_nc_ls = df_int_nc_ls[df_int_nc_ls['AtmosLayer'] != 'Middle Stratosphere']
    
    df_int_saod_ls = df_int[df_int['Property'] == 'sAOD']
    df_int_saod_ms = df_int_saod_ls[df_int_saod_ls['AtmosLayer'] == 'Stratosphere']
     #df_int_saod_ls = df_int_saod_ls[df_int_saod_ls['AtmosLayer'] != 'Middle Stratosphere']
    dfl_sf = dfl_sf.sort_values(by=['Date','Vertically Integrated Total'])
    
    #Figure 6
   
    #plotting palette
    palette1 = {"POPS":"darkblue",
                "POPS > 205 nm": "black",
                "CARMA":"maroon",
                "CARMA calc.":"brown",
                "Lidar":"teal"
                    }
   

   
    #dates of stratospheric perturbations
    d1= datetime.strptime('2019-12-29', '%Y-%m-%d').date() #'2019-12-29') #pyroCB
    d1_vl = d1.toordinal()
    d2= datetime.strptime('2019-06-26', '%Y-%m-%d').date()  #('2019-06-26')
    d2_vl = d2.toordinal()
    d3= datetime.strptime('2019-08-03', '%Y-%m-%d').date()  #('2019-08-03')
    d3_vl = d3.toordinal()
    d4= datetime.strptime('2021-04-09', '%Y-%m-%d').date() #('2021-04-09')
    d4_vl = d4.toordinal()
   


   
    #ax = g.facet_axis(0, 1)
    ax.fill_between(x=df_int_nc_ms['Date'], y1=df_int_nc_ms['Vertically Integrated Total']*0.95, y2 = df_int_nc_ms['Vertically Integrated Total']*1.05, color="blue", alpha=0.3)
    ax.set_title("Middle Stratosphere")
    plt.axvline(x=d1_vl, color = 'm', linestyle = '--')
    plt.axvline(x=d2_vl, color = 'm', linestyle = '--')
    plt.axvline(x=d3_vl, color = 'm', linestyle = '--')
    plt.axvline(x=d4_vl, color = 'm', linestyle = '--')
   
    #g = ax.axes.flatten()
    #g[0].set_title("Lower/Lowermost Stratosphere")
    #g[1].set_title("Middle Stratosphere")
   

       
   
   
   
#    df_int_all["Source"]= df_int_all["Source"].str.replace("CARMA >= 200 nm", "CARMA > 198 nm", case = False)  
#    df_int_all["Source"]= df_int_all["Source"].str.replace("POPS", "POPS > 140 nm", case = False) 
#    df_int_all["Source"]= df_int_all["Source"].str.replace("CARMA", "CARMA > 150 nm", case = False) 

    seaborn.set(style="whitegrid")
    g = seaborn.relplot(x="Date", y="Vertically Integrated Total", data=df_int_all, hue="Source", hue_order = ['POPS', 'CARMA'], row="Property", row_order = ["Number Concentration Column(# m^-2)"], col = "AtmosLayer", 
                          sort=False, legend = True, palette = palette1, zorder = 0, kind = "line", marker = 's')
    g.set(xlabel="Date", ylabel = u"Number Concentration Column (# $\mathregular{m^{-2}}$)") 
    g.set_xticklabels(ord_labels, rotation=45)
    g.set(xlim = (737100, 738300))
   
    ax = g.facet_axis(0, 0)
    ax.fill_between(x=df_int_nc_ls['Date'], y1=(df_int_nc_ls['Vertically Integrated Total']-df_int_nc_ls['Unc']), y2 = (df_int_nc_ls['Vertically Integrated Total']+df_int_nc_ls['Unc']), color="blue", alpha=0.3)
    ax.fill_between(x=dfc_int_nc_ls['Date'], y1=(dfc_int_nc_ls['Vertically Integrated Total Min']), y2 = (dfc_int_nc_ls['Vertically Integrated Total Max']), color="maroon", alpha=0.3)
    ax.set_title("")
    plt.axvline(x=d1_vl, color = 'k', linestyle = '--')
    plt.axvline(x=d2_vl, color = 'k', linestyle = '--')
    plt.axvline(x=d3_vl, color = 'k', linestyle = '--')
    plt.axvline(x=d4_vl, color = 'k', linestyle = '--')
   
    df_int_all["Source"]= df_int_all["Source"].str.replace("CARMA", "CARMA - 521 nm", case = False)  
    df_int_all["Source"]= df_int_all["Source"].str.replace("POPS", "POPS - 532 nm", case = False)  
    df_int_all["Source"]= df_int_all["Source"].str.replace("Lidar", "Lidar - 532 nm", case = False)  
    
    
    
    palette0 = {"POPS - 532 nm":"darkblue",
                 "CARMA - 521 nm":"maroon",
                 "CARMA calc.":"brown",
                 "Lidar - 532 nm":"teal"}
                
     #fig, ax = plt.subplots(1, 1, figsize=(8, 4))
     #seaborn.set(rc={"xtick.bottom" : True, "ytick.left" : True})
    seaborn.set(style="whitegrid")
    g = seaborn.relplot(x="Date", y="Vertically Integrated Total", data=df_int_all, hue="Source", hue_order = ['Lidar - 532 nm', 'POPS - 532 nm', 'CARMA - 521 nm'], row="Property", row_order = ["sAOD"], col = "AtmosLayer", 
                           sort=True, legend = True, palette = palette0, zorder = 0, kind = "line", marker = 's')
    g.set(xlabel='Date', ylabel = 'sAOD')
    g.set_xticklabels(ord_labels, rotation=45)
    g.set(xlim = (737100, 738300))
    
    POPS_Unc_lb = [0.85, 0.87, 0.8, 0.6, 0.8, 0.83, 0.9, 0.75]
    POPS_Unc_ub = [1.15, 1.13, 1.2, 1.4, 1.2, 1.17, 1.1, 1.25]
    ax = g.facet_axis(0, 0)
     #ax.fill_between(x=df_int_saod_ls['Date'], y1=df_int_saod_ls['Vertically Integrated Total']*POPS_Unc_lb, y2 = df_int_saod_ls['Vertically Integrated Total']*POPS_Unc_ub, color="blue", alpha=0.3)
    ax.fill_between(x=df_int_saod_ls['Date'], y1=(df_int_saod_ls['Vertically Integrated Total']-df_int_saod_ls['Unc']), y2 = (df_int_saod_ls['Vertically Integrated Total']+df_int_saod_ls['Unc']), color="blue", alpha=0.3)
     #ax.fill_between(x=dfl_sf['Date'], y1=(dfl_sf['Vertically Integrated Total']- dfl_sf['Unc']), y2 = (dfl_sf['Vertically Integrated Total']+ dfl_sf['Unc']), color="green", alpha=0.3)
    ax.fill_between(x=dfl_sf['Date'], y1=(dfl_sf['Vertically Integrated Total']*0.75), y2 = (dfl_sf['Vertically Integrated Total']*1.25), color="green", alpha=0.3)
    ax.set_title("")
    plt.axvline(x=d1_vl, color = 'k', linestyle = '--')
    plt.axvline(x=d2_vl, color = 'k', linestyle = '--')
    plt.axvline(x=d3_vl, color = 'k', linestyle = '--')
    plt.axvline(x=d4_vl, color = 'k', linestyle = '--')
    
        
