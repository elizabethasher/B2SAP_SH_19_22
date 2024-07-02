#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 20:57:09 2024

@author: asher
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 17:40:30 2021

@author: easher
"""
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import sys
import os
from datetime import datetime as dt
from datetime import timedelta
import time
import pandas as pd
import numpy as np
import glob
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from matplotlib import ticker
import matplotlib.colors as colors
import seaborn as seaborn
from scipy import stats
from scipy.interpolate import griddata

#Code to plot JGR Figure 4
def Plume(Launch, Altitude):
    if (Launch == '2020-01-27'):
        if (Altitude == '17.5'):
            composition = 'Smoky Auto'
        elif (Altitude == '17.0'):
                composition = 'Smoky Auto'
        elif (Altitude == '16.5'):
                    composition = 'Smoky Auto'
        elif (Altitude == '16.0'):
            composition = 'Smoky Auto'
     
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
    else:
        composition = 'Sulfate'
    return composition
    #composition = 'Sulfate' #comment this out to change...
    return composition

def DnDlogDp(Launch, Composition, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15):
    file_ior = '/Users/asher/Documents/NOAAweb/NZ_Paper_IOR.csv'
    dfIOR = pd.read_csv(file_ior,sep=',', dtype=None, engine='python')
    #dLogDp = DLogDpLookUpTable(Composition)
    
    if Composition == 'SulfateAuto':
        dp_in = dfIOR[(dfIOR.Composition == 'Sulfate') & (dfIOR.Binning == 'Automatic')]
        dLogDp = np.array(dp_in.dLogDp)

        #dNdLogDp = np.divide(np.array([B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15, B16, B17, B18, B19, B20]), np.array(dLogDp)).tolist()
    elif Composition == 'Sulfate':
        dp_in = dfIOR[(dfIOR.Composition == 'Sulfate') & (dfIOR.Binning == 'Manual')]
        dLogDp = np.array(dp_in.dLogDp)
        
    elif Composition == 'SmokyAuto':
        dp_in = dfIOR[(dfIOR.Composition == 'Smoky') & (dfIOR.Binning == 'Automatic')]
        dLogDp = np.array(dp_in.dLogDp)
            
    else:
        dp_in = dfIOR[(dfIOR.Composition == 'Smoky') & (dfIOR.Binning == 'Manual')]
        dLogDp = np.array(dp_in.dLogDp)

    ndp = np.array([B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15])
    dNdLogDpA = np.divide(ndp, dLogDp)
        
    dNdLogDp = dNdLogDpA.tolist()   
    return dNdLogDp


#set up 
def BinDia(Bin, composition):
    if Bin == 'B0' and composition == 'Smoky':
        Dia = 140
    elif Bin == 'B1' and composition == 'Smoky':
        Dia = 169
    elif Bin == 'B2' and composition == 'Smoky':
        Dia = 203 #204
    elif Bin == 'B3' and composition == 'Smoky':
        Dia = 242 #243
    elif Bin == 'B4' and composition == 'Smoky':
        Dia = 300 #301
    elif Bin == 'B5' and composition == 'Smoky':
        Dia = 391
    elif Bin == 'B6' and composition == 'Smoky':
        Dia = 482 #484
    elif Bin == 'B7' and composition == 'Smoky':
        Dia = 688 #641
    elif Bin == 'B8' and composition == 'Smoky':
        Dia = 1028 #1048
    elif Bin == 'B9' and composition == 'Smoky':
        Dia = 1202 #1316
    elif Bin == 'B10' and composition == 'Smoky':
        Dia = 1352 #1449
    elif Bin == 'B11' and composition == 'Smoky':
        Dia = 1727 #1816
    elif Bin == 'B12' and composition == 'Smoky':
        Dia = 2268 #2367
    elif Bin == 'B13' and composition == 'Smoky':
        Dia = 2930 #3005
    elif Bin == 'B14' and composition == 'Smoky':
        Dia = 3638 #3652
    elif Bin == 'B0' and composition == 'Smoky Auto':
        Dia = 133
    elif Bin == 'B1' and composition == 'Smoky Auto':
        Dia = 146
    elif Bin == 'B2' and composition == 'Smoky Auto':
        Dia = 162
    elif Bin == 'B3' and composition == 'Smoky Auto':
        Dia = 181
    elif Bin == 'B4' and composition == 'Smoky Auto':
        Dia = 207
    elif Bin == 'B5' and composition == 'Smoky Auto':
        Dia = 247
    elif Bin == 'B6' and composition == 'Smoky Auto':
        Dia = 351 #353
    elif Bin == 'B7' and composition == 'Smoky Auto':
        Dia = 480 #483
    elif Bin == 'B8' and composition == 'Smoky Auto':
        Dia = 593 #595
    elif Bin == 'B9' and composition == 'Smoky Auto':
        Dia = 962 #830
    elif Bin == 'B10' and composition == 'Smoky Auto':
        Dia = 1484 #1380
    elif Bin == 'B11' and composition == 'Smoky Auto':
        Dia = 2052 #2128
    elif Bin == 'B12' and composition == 'Smoky Auto':
        Dia = 2946 #3118
    elif Bin == 'B13' and composition == 'Smoky Auto':
        Dia = 3744 #3868
    elif Bin == 'B14' and composition == 'Smoky Auto':
        Dia = 4000
    elif Bin == 'B0' and composition == 'Sulfate':
        Dia = 150
    elif Bin == 'B1' and composition == 'Sulfate':
        Dia = 182
    elif Bin == 'B2' and composition == 'Sulfate':
        Dia = 222 #221
    elif Bin == 'B3' and composition == 'Sulfate':
        Dia = 269 #269
    elif Bin == 'B4' and composition == 'Sulfate':
        Dia = 328 #326
    elif Bin == 'B5' and composition == 'Sulfate':
        Dia = 401 #406
    elif Bin == 'B6' and composition == 'Sulfate':
        Dia = 496 #493
    elif Bin == 'B7' and composition == 'Sulfate':
        Dia = 593 #590
    elif Bin == 'B8' and composition == 'Sulfate':
        Dia = 693 #722
    elif Bin == 'B9' and composition == 'Sulfate':
        Dia = 791 #905
    elif Bin == 'B10' and composition == 'Sulfate':
        Dia = 940 #1205
    elif Bin == 'B11' and composition == 'Sulfate':
        Dia = 1285 #1231
    elif Bin == 'B12' and composition == 'Sulfate':
        Dia = 1717 #1320
    elif Bin == 'B13' and composition == 'Sulfate':
        Dia = 2021 #1858
    elif Bin == 'B14' and composition == 'Sulfate':
        Dia = 2375 #2244
    elif Bin == 'B0' and composition == 'Sulfate Auto':
        Dia = 142
    elif Bin == 'B1' and composition == 'Sulfate Auto':
        Dia = 157 #156
    elif Bin == 'B2' and composition == 'Sulfate Auto':
        Dia = 175 #174
    elif Bin == 'B3' and composition == 'Sulfate Auto':
        Dia = 197 #196
    elif Bin == 'B4' and composition == 'Sulfate Auto':
        Dia = 226 #224
    elif Bin == 'B5' and composition == 'Sulfate Auto':
        Dia = 275 #272
    elif Bin == 'B6' and composition == 'Sulfate Auto':
        Dia = 376 #371
    elif Bin == 'B7' and composition == 'Sulfate Auto':
        Dia = 494 #490
    elif Bin == 'B8' and composition == 'Sulfate Auto':
        Dia = 593 #593
    elif Bin == 'B9' and composition == 'Sulfate Auto':
        Dia = 763 #775
    elif Bin == 'B10' and composition == 'Sulfate Auto':
        Dia = 1057 #1022
    elif Bin == 'B11' and composition == 'Sulfate Auto':
        Dia = 1246 #1341
    elif Bin == 'B12' and composition == 'Sulfate Auto':
        Dia = 1735 #1921
    elif Bin == 'B13' and composition == 'Sulfate Auto':
        Dia = 2587 #2601
    elif Bin == 'B14' and composition == 'Sulfate Auto':
        Dia = 3382 #2423
    
    else:
        Dia = -99
    return Dia



df = pd.DataFrame()



#plot dNdLogDp with altitude for Lauder NZ
Dia = np.array([154.83, 187.63, 227.39, 275.56, 333.94, 404.69, 490.43, 594.34,
        720.26, 872.86, 1057.79, 1281.89, 1553.48, 1882.61, 2281.47])

#Dia = np.divide(Dia,1000)
Rad = np.divide(Dia,2)
Area = np.pi*4*np.power(Rad,2)

file = '/Users/asher/Documents/PapersInProgress/JGR_ANYSO/Submission/POPSDataForPublication/POPS/POPS_table.csv'
df_all = pd.read_csv(file, delimiter = ',', header = 14,  dtype=None, low_memory=False)
fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(ncols=5,sharex=True, sharey=True)

for x in range(0,6):
    count = x
    if count == 0:
        df1 = df_all[df_all['Launch'] == '2020-06-19']
    if count == 1:
        df1 = df_all[df_all['Launch'] == '2020-07-23']
    if count == 2:
        df1 = df_all[df_all['Launch'] == '2020-11-17']
    if count == 3:
        df1 = df_all[df_all['Launch'] == '2019-09-03']
    if count == 4:
        df1 = df_all[df_all['Launch'] == '2020-01-27']
    if count == 5:
        df1 = df_all[df_all['Launch'] == '2021-08-31']

    
    df = df1
    #calculate altitude levels
    AltLabels = ["0", "0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5","5.5","6","6.5","7","7.5","8","8.5","9","9.5","10",\
                  "10.5","11","11.5","12","12.5","13","13.5","14","14.5","15","15.5","16","16.5","17","17.5","18","18.5","19","19.5",\
                      "20","20.5", "21","21.5","22","22.5","23","23.5","24","24.5","25","25.5","26","26.5","27","27.5","28","28.5","29",\
                          "29.5"]
    #create more altitude bins
    df['Altitudekm'] = pd.cut(df['Altitude (km)'], bins=np.linspace(0, 30, 61),labels = AltLabels, include_lowest=True) 
    
    
    dfad = df
    del df
    
    #mean of ascent and descent
    df = dfad
    df = df.sort_values(by=['Flight','Altitudekm'])
  
    
    #ID composition
    df['Composition'] = df.apply(lambda x: Plume(x['Launch'], x['Altitudekm']), axis=1)
    #calculate dNdLogDp instead of aerosol concentration in each bin...
    df['B'] = df.apply(lambda x: DnDlogDp(x['Launch'], x['Composition'], 
                                x['B1'], x['B2'], x['B3'], x['B4'], x['B5'], x['B6'], x['B7'], x['B8'], x['B9'], x['B10'], x['B11'], 
                                                  x['B12'], x['B13'], x['B14'], x['B15']), axis=1)

    #drop aero conc columns
    df.drop(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'B12', 'B13', 'B14', 'B15'], axis = 1, inplace = True)
 
    a = pd.DataFrame(df['B'])
    b = pd.concat([pd.DataFrame(a['B'].values.tolist()) for c in a.columns], 
                      axis=1, 
                      keys=a.columns)

    b.columns = ['{}{}'.format(i, j) for i, j in b.columns]
    df.reset_index(drop=False, inplace = True)
    b.reset_index(drop=False, inplace = True)
    df = pd.concat([df, b], axis = 1)
    df.reset_index(drop=True, inplace = True)

    
    #df1= df.groupby(['Launch','Altitudekm', 'Flight'],as_index=False).mean()
    df1= df.groupby(['Launch','Altitudekm', 'Composition'],as_index=False).mean()
    df2 = df1[df1['Altitude (km)'] <= 28.5]
    #df2 = df2.filter(['Launch', 'Altitudekm', 'Flight', 'B1', 'B2','B3','B4', 'B5', 'B6', 'B7','B8', 'B9', 'B10', 'B11', 'B12', 'B13', 'B14','B15'])
    #seperate out ascent and descent
    df2['Altitude_km'] = pd.to_numeric(df2['Altitudekm'], errors='coerce').fillna(0)

    
    #create long format dataframe 
    #long format with dNdlogdp
    dfNew = pd.melt(df2, id_vars=['Altitudekm', 'Launch', 'Composition'], value_vars=['B0', 'B1', 'B2','B3','B4', 'B5', 'B6', 'B7','B8', 'B9', 'B10', 'B11', 'B12', 'B13', 'B14'],var_name = 'Bin',value_name ='dNdLogDp')
    dfNew['Diameter'] = dfNew.apply(lambda x: BinDia(x['Bin'], x['Composition']), axis=1)
    
    #select out particular altitude from ascent and descent datafrmes
    
    ngridx = 350#350
    ngridy = 690 #690
    ngridy2 = 800 #800
    
    # xx = np.random.uniform(-2, 2, 100)
    # yy = np.random.uniform(-2, 2, 100)
    # zz = xx * np.exp(-xx**2 - yy**2)
    dfNewA = dfNew
    
    x = dfNewA['Diameter'].tolist()
    y = dfNewA['Altitudekm'].tolist()
    z = dfNewA['dNdLogDp'].tolist()
    
    # x2 = dfNewD['Diameter'].tolist()
    # y2= dfNewD['Altitudekm'].tolist()
    # z2 = dfNewD['aerosol concentration'].tolist()
    
    x3 = [255, 7005]
    y3 = [8.5, 11.5]
    
    # -----------------------
    # Interpolation on a grid
    # -----------------------
    # A contour plot of irregularly spaced data coordinates
    # via interpolation on a grid.
    
    # Create grid values first.
    xi = np.linspace(100, 2500, num=ngridx) #change this to smaller than bin increments
    yi = np.linspace(7.9, 28.1, num=ngridy) #change this to altitudinal bins 100 m??
    yi2 = np.linspace(7.9, 28.1,num=ngridy2)
    
    # Linearly interpolate the data (x, y) on a grid defined by (xi, yi).
    triang = tri.Triangulation(x, y)
    interpolator = tri.CubicTriInterpolator(triang, z)
    Xi, Yi = np.meshgrid(xi, yi)
    Zi = interpolator(Xi, Yi)

    
    # # Note that scipy.interpolate provides means to interpolate data on a grid
    # # as well. The following would be an alternative to the four lines above:
    # #from scipy.interpolate import griddata
    zi= griddata((x, y), z, (xi[None, :], yi2[:, None]), method='linear', rescale='True')
    seaborn.set(rc={"xtick.bottom" : True, "ytick.left" : True})
    clevels = [1, 3, 10, 30, 60]
    xlims = [155, 2200]
    ylims = [9, 24]

    
    
    #plotting each of these launches a seperate figure panel
    if count == 0:
        #ax3.plot(x, y, 'k.', ms=1)
        ax3.contour(xi, yi2, zi, linewidths=1, levels=clevels, colors='k') #plot interpolated contours 
        cntr3 = ax3.contourf(xi, yi2, zi,  levels=clevels, extend = 'both', cmap="RdBu_r") #vmin = 0, vmax = 30,
        #fig.colorbar(cntr1, ax=ax1)#plot actual median altide bin etc...
        # 
        ax3.set_title('2020-06-19')
        ax3.set_xlabel('Dp (nm)', fontsize=13)
        #ax1.set_title('CO 07-19-2021')
        ax3.set_xlim(xlims)
        ax3.set_xscale('log')
        ax3.set_ylim(ylims)

    
    if count == 1:
        ax4.contour(xi, yi2, zi,  linewidths=1, levels=clevels, colors='k') #plot interpolated contours           locator=ticker.LogLocator(),
        cntr4 = ax4.contourf(xi, yi2, zi,  levels=clevels, extend = 'both', cmap="RdBu_r") #locator=ticker.LogLocator(),
        #fig.colorbar(cntr2, ax=ax2)#plot actual median altide bin etc...
        #ax4.plot(x, y, 'k.', ms=2) 
        ax4.set_title('2020-07-23')
        ax4.set_xlim(xlims)
        ax4.set_xscale('log')
        ax4.set_ylim(ylims)
    

    
    if count == 2:
        ax5.contour(xi, yi2, zi, linewidths=1, levels=clevels, colors='k') #plot interpolated contours  locator=ticker.LogLocator(), 
        cntr5 = ax5.contourf(xi, yi2, zi, levels=clevels, extend = 'both', cmap="RdBu_r") #locator=ticker.LogLocator(), 
        fig.colorbar(cntr5, ax=ax5, label = 'dNdLogDp' )
        #fig.colorbar(cntr5, ax=ax5, label = "Aerosol Concentration (# $\mathregular{cm^{-3}}$)")
        #ax5.plot(x, y, 'k.', ms=2) 
        ax5.set_title('2020-11-17')
        ax5.set_xlim(xlims)
        ax5.set_xscale('log')
        ax5.set_ylim(ylims)
        
    if count == 3:
        ax1.contour(xi, yi2, zi,  linewidths=1, levels=clevels, colors='k') #plot interpolated contours locator=ticker.LogLocator(),
        cntr1 = ax1.contourf(xi, yi2, zi, levels=clevels, extend = 'both', cmap="RdBu_r") #locator=ticker.LogLocator(),
        #fig.colorbar(cntr4, ax=ax3, label = 'Aerosol Concentration (#/cm3)' )
        #ax1.plot(x, y, 'k.', ms=2) 
        ax1.set_title('2019-09-03')
        ax1.set_xlim(xlims)
        ax1.set_xscale('log')
        ax1.set_ylim(ylims)
        ax1.set_ylabel('Altitude (km) ', fontsize=13)
    
    if count == 4:
        ax2.contour(xi, yi2, zi, linewidths=1, levels=clevels, colors='k') #plot interpolated contours locator=ticker.LogLocator(),
        cntr2 = ax2.contourf(xi, yi2, zi,  levels=clevels, extend = 'both', cmap="RdBu_r") #locator=ticker.LogLocator(),
        #fig.colorbar(cntr5, ax=ax5, label = 'Aerosol Concentration (#/cm3)' )
        #ax2.plot(x, y, 'k.', ms=2) 
        ax2.set_title('2020-01-27')
        ax2.set_xlim(xlims)
        ax2.set_xscale('log')
        ax2.set_ylim(ylims)


plt.show()



