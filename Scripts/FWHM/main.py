from Common import config, Spectra
import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import CubicSpline , interp1d
from heapq import nsmallest
import argparse
import os




def Read_Data(Route):
    Data = pd.read_csv(Route, delimiter = ' ', header = None)
    Data.columns = ['L','I', 'C']
    Data['L'] = Data['L']*10 #units Amstrong
    return Data

def Cut_Data(configuration,Data):
    #Recortamos en el rango de la línea K
    New_Data = Data[(Data['L'] >= configuration['LMIN']) & (Data['L'] <= configuration['LMAX'])] #Se encuentran los valores cercanos a la config inicial
    New_Data.index = list(range(len(New_Data)))
    L_min = New_Data['L'].min()
    L_max = New_Data['L'].max()
    return [New_Data, L_min, L_max]

def step_continuum(Lmin,Lmax,splines):
    return (Lmax-Lmin)/splines

def Points_continuum(_lambda, Intensity,start,end,step):
    
    steps_list = np.arange(start,end, step)
    grouped_data = {}
    grouped_data['L'] = []
    grouped_data['I'] = []
    for i in range(1,len(steps_list)):
        grouped_lambda = _lambda[(_lambda >= steps_list[i-1]) & (_lambda < steps_list[i])]
        
        grouped_data['L'].append(grouped_lambda)
        grouped_data['I'].append(Intensity[  grouped_lambda.index  ]) # Agrup_lamnda.index  obtiene la llave compartida
    return grouped_data

def fit_continuum(grouped_data, L_min, L_max):  
    newdataL = [] 
    newdataI = []
    for i in range(len(grouped_data['I'])):
        newdataL.append(grouped_data['L'][i].median())
        newdataI.append(grouped_data['I'][i].median())
        data_fitL = [x for x in newdataL if str(x) != 'nan']
        data_fitI = [x for x in newdataI if str(x) != 'nan']
    cs = InterpolatedUnivariateSpline(data_fitL,data_fitI)
    xs =  np.linspace(L_min,  L_max, 100)
    fit =  pd.DataFrame({'L': xs, 'I': cs(xs)})
    return fit

def parameters_double(fit, LMIN, LMAX, LLIN, RANG):
#Creo dos dataFrame desde el mínimo K3 hacía la izquiera y la derecha, teniendo en cuenta un rango de ese mínimo
    b = fit[     (fit['L'] >= LMIN) & (fit['L'] <= (LLIN - RANG))  ]
    r = fit[     (fit['L'] >= (LLIN + RANG)) & (fit['L'] <= LMAX)  ]
    b.index = list(range(len(b)))
    r.index = list(range(len(r)))

    #Calculo los máximos y mínimos para cada uno de los picos
    k1b=b[b['I'] == b['I'].min()]
    k2b=b[b['I'] == b['I'].max()]
    k1r=r[r['I'] == r['I'].min()]
    k2r=r[r['I'] == r['I'].max()]
    #Se usa float() para poder hacer operaciones entre datos de un dataFrame
    #Hago un DataFrame que va desde K1 a K2, para B y R
    peak1 = b[(b['L'] >= float(k1b['L'])) & (b['L'] <= float(k2b['L']))]
    peak2 = r[(r['L'] >= float(k2r['L'])) & (r['L'] <= float(k1r['L']))]
    peak1.index = list(range(len(peak1)))
    peak2.index = list(range(len(peak2)))
    #Calculo la intensidad media entre K1 y K2 para B y R
    Ibm_teo = float(k1b['I']) + (float(k2b['I']) - float(k1b['I']))/2 #Intensidad teórica media
    Ibm_ajus = nsmallest(1, peak1['I'],  key = lambda x: abs(x - Ibm_teo))[0] #El más cercano en el ajuste a la intensidad media
    B_m = peak1[peak1['I'] == Ibm_ajus] 

    Irm_teo = float(k1r['I']) + (float(k2r['I']) - float(k1r['I']))/2
    Irm_ajus = nsmallest(1, peak2['I'],  key = lambda x: abs(x - Irm_teo))[0]
    R_m = peak2[peak2['I'] == Irm_ajus]

    plt.figure(figsize = [8,6])
    plt.plot(fit['L'],fit['I'],marker='', ls='-')
    plt.scatter(k1b['L'], k1b['I'], color = 'blue')
    plt.scatter(k1r['L'], k1r['I'], color = 'r')
    plt.scatter(k2b['L'], k2b['I'], color = 'blue')
    plt.scatter(k2r['L'], k2r['I'], color = 'r')
    plt.scatter(B_m['L'], B_m['I'], color = 'green')
    plt.scatter(R_m['L'], R_m['I'], color = 'green')
    plt.grid()
    plt.ylabel("Flujo",fontsize=18)
    plt.xlabel("$\lambda [A]$",fontsize=18)
    plt.show()
    return ""

def Analysis(Analysis_Types):    
    for Data_route in ["zetaAur-eclipse_B_2019_11_15_23_52_17.dat"]:
        [New_Data,L_min, L_max] = Cut_Data(  config()[Analysis_Types]  ,  Read_Data(Data_route)  )

        grouped_data = Points_continuum(New_Data['L'],
                                        New_Data['I'],
                                        L_min,
                                        L_max,
                                        step_continuum(L_min,L_max,config()[Analysis_Types]['NSPL']) )
        fit = fit_continuum(grouped_data, L_min, L_max)
        parameters = parameters_double( fit, 
                                        config()[Analysis_Types]['LMIN'], 
                                        config()[Analysis_Types]['LMAX'],
                                        config()[Analysis_Types]['LLIN'],
                                        config()[Analysis_Types]['RANG'] )

def Weight_equivalent(parameter_list):
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    Analysis_Types = list(config().keys()) #Obtiene lista de opciones para analizar

    parser.add_argument('Analysis_Types', 
                        help = 'Enter the name of line to analyzed',
                        type= str,
                        choices=Analysis_Types) #Agrega argumentos para ejecucion del usuario
    args = parser.parse_args()
    Analysis(args.Analysis_Types)


