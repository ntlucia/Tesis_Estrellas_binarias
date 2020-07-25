#// ============================================================================================= //
#//                                                                                               //
#//       Filename:  main.py                                                                      //
#//    Description:  Equivalent width and half height measurement program for spectroscopic lines //
#//                                                                                               //
#//        Version:  1                                                                            //
#//        Created:  25/07/2020                                                                   //
#//       Compiler:  Python                                                                       //
#//                                                                                               //
#//         Author:  Natalia Lucía Oliveros Gómez                                                 //
#//          Email:  onatalialucia@gmail.com                                                      //
#//        Company:  Grupo Halley UIS - Grupo Física Estelar U. Gto                               //
#//                                                                                               //
#// ============================================================================================= //
#
#// ============================================================================================= //
#//                                                                                               // 
#//  Compile with:                                                                                //
#//                                                                                               //
#//  main.py "type of width to measure" "type of line"                                            //
#//                                                                                               //
#// ============================================================================================= //


# Necessary libraries
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
from scipy import integrate
import datetime


#Import Data to analyze
def Read_Data(Route):
    Data = pd.read_csv(Route, delimiter = ' ', header = None)
    Data.columns = ['L','I', 'C']
    Data['L'] = Data['L']*10 #units Amstrong
    return Data

#Cut only the line you are interested in
def Cut_Data(configuration,Data):
    New_Data = Data[(Data['L'] >= configuration['LMIN']) & (Data['L'] <= configuration['LMAX'])]
    New_Data.index = list(range(len(New_Data)))
    L_min = New_Data['L'].min()
    L_max = New_Data['L'].max()
    return [New_Data, L_min, L_max]


def julian_day(Data_route):
    name_data = Data_route[ len(Data_route)-41 : ]
    fecha_dt = datetime.datetime.strptime(name_data, 'zetaAur-eclipse_B_%Y_%m_%d_%H_%M_%S.dat')
    julianday = fecha_dt.toordinal() + 1721425
    return julianday
#"--------------------------------------------------------------------------"#
#"                     STAGE 1: Half Intensity Emission Width               "#
#"--------------------------------------------------------------------------"#

#Number of splines
def step_continuum(Lmin,Lmax,splines):
    return (Lmax-Lmin)/splines

#New dots to found continuous
def Points_continuum(_lambda, Intensity,start,end,step):
    
    steps_list = np.arange(start,end, step)
    grouped_data = {}
    grouped_data['L'] = []
    grouped_data['I'] = []
    for i in range(1,len(steps_list)):
        grouped_lambda = _lambda[(_lambda >= steps_list[i-1]) & (_lambda < steps_list[i])]
        
        grouped_data['L'].append(grouped_lambda)
        grouped_data['I'].append(Intensity[  grouped_lambda.index  ])
    return grouped_data

#Fit contunuous to line data
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

#Width calculation
def FWHM_parameters_double(New_Data, fit, LMIN, LMAX, LLIN, RANG):
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


    #Calculo de la incertidumbre de lambda para el punto medio

    valor_medioB = New_Data[New_Data['L'] == nsmallest(1, New_Data['L'],  key = lambda x: abs(x - float(B_m['L'])))[0]]
    indb = valor_medioB.index[0]

    if float(valor_medioB['L']) < float(B_m['L']):
        Lmbinf = valor_medioB
        Lmbsup = pd.DataFrame( [[New_Data['L'][indb+1], New_Data['I'][indb+1]]], columns = ('L','I'))
    else:
        Lmbsup = valor_medioB
        Lmbinf = pd.DataFrame( [[New_Data['L'][indb-1], New_Data['I'][indb-1]]], columns = ('L','I') )

    valor_medioR = New_Data[New_Data['L'] == nsmallest(1, New_Data['L'],  key = lambda x: abs(x - float(R_m['L'])))[0]]
    indr = valor_medioR.index[0]

    if float(valor_medioR['L']) < float(R_m['L']):
        Lmrinf = valor_medioR
        Lmrsup = pd.DataFrame( [[New_Data['L'][indr+1], New_Data['I'][indr+1]]], columns = ('L','I') )
    else:
        Lmrsup = valor_medioR
        Lmrinf = pd.DataFrame( [[New_Data['L'][indr-1], New_Data['I'][indr-1]]], columns = ('L','I') )

    #Calculo de W_0
    W_0 = float(R_m['L']) -float(B_m['L'])
    deltabm = ( float(Lmbsup['L']) - float(B_m['L']) + float(B_m['L']) - float(Lmbinf['L']))/2
    deltarm = ( float(Lmrsup['L']) - float(R_m['L']) + float(R_m['L']) - float(Lmrinf['L']))/2
    deltaW = np.sqrt(deltabm**2 + deltarm**2)

    return [W_0, deltaW]

#"--------------------------------------------------------------------------"#
#"                           STAGE 2: Equivalent Width                      "#
#"--------------------------------------------------------------------------"#

def EW_parameters_double(Data, New_Data, L_min, L_max, LMIN, LMAX):
    pseudo_continuous1 = Data[(Data['L'] >= LMIN -4) & (Data['L'] <= LMIN)]
    pseudo_continuous2 = Data[(Data['L'] >= LMAX) & (Data['L'] <= LMAX + 4)]
    pseudo_continuous1.index = list(range(len(pseudo_continuous1)))
    pseudo_continuous2.index = list(range(len(pseudo_continuous2)))

    mean1 = np.mean(pseudo_continuous1)
    mean2 = np.mean(pseudo_continuous2)
    mean = (mean1['I'] + mean2['I'])/2

    New_Data['I'] = New_Data['I'] + 1 - mean

    Aprom = 1*(L_max - L_min)
    Abajo_abs = integrate.simps(New_Data['I'], New_Data['L'])
    Area = Aprom - Abajo_abs
    EW = Area
    return EW

def Analysis_WE(Analysis_Types):
    WE = pd.DataFrame(columns=("Julian_day","WE"))
    Direction_base = Spectra()
    with open(Direction_base + "archives.txt",mode='r') as f:
        Records =  f.readlines()
    for i in range(len(Records)):
        Records[i] = Direction_base + Records[i][:-1]
        
    for Data_route in Records:
        print("Processing {}".format(Data_route))
        Data = Read_Data(Data_route)
        [New_Data,L_min, L_max] = Cut_Data(  config()[Analysis_Types]  ,  Read_Data(Data_route)  )
        
        try:
            Equivalent_Weigth = EW_parameters_double(Data, 
                                            New_Data,
                                            L_min, 
                                            L_max,
                                            config()[Analysis_Types]['LMIN'], 
                                            config()[Analysis_Types]['LMAX'] )
            line = [julian_day(Data_route),Equivalent_Weigth]
            WE.loc[ len(WE["WE"]) ] =  line
            #WE.loc[ len(WE["Julian_day"])] = julian_day(Data_route)
        except IndexError:
            print("Index Error")
    return WE

def Analysis_FWHM(Analysis_Types):
    WE = pd.DataFrame(columns=("Julian_day","W_0","deltaW"))
    Direction_base = Spectra()
    with open(Direction_base + "archives.txt",mode='r') as f:
        Records =  f.readlines()
    for i in range(len(Records)):
        Records[i] = Direction_base + Records[i][:-1]
        
    for Data_route in Records:
        print("Processing {}".format(Data_route))
        [New_Data,L_min, L_max] = Cut_Data(  config()[Analysis_Types]  ,  Read_Data(Data_route)  )

        grouped_data = Points_continuum(New_Data['L'],
                                        New_Data['I'],
                                        L_min,
                                        L_max,
                                        step_continuum(L_min,L_max,config()[Analysis_Types]['NSPL']) )
        fit = fit_continuum(grouped_data, L_min, L_max)
        try:
            parameters = FWHM_parameters_double(New_Data,
                                            fit, 
                                            config()[Analysis_Types]['LMIN'], 
                                            config()[Analysis_Types]['LMAX'],
                                            config()[Analysis_Types]['LLIN'],
                                            config()[Analysis_Types]['RANG'] )
            
            line = [julian_day(Data_route),parameters[0], parameters[1]]
            WE.loc[ len(WE["W_0"]) ] =  line

        except IndexError:
            print("Index Error")
    return WE
    
#"--------------------------------------------------------------------------"#
#"                         STAGE 3: Decisions to calculate                  "#
#"--------------------------------------------------------------------------"#

def Analysis(Type_width, Analysis_Types):    
    if(Type_width == "FWHM"):
        WE = Analysis_FWHM(Analysis_Types)
    elif(Type_width == "WE"):
        WE = Analysis_WE(Analysis_Types)
    else:
        return
    name_data = "Results/{}_{}.csv".format(Type_width,Analysis_Types)
    WE.to_csv(name_data, index = False, header=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    Analysis_Types = list(config().keys()) #Obtiene lista de opciones para analizar
    Type_width = ["FWHM","EW"]
    parser.add_argument('Type_width', 
                        help = 'Enter the type width to measure WE/FWHM',
                        type= str) #Agrega argumentos para ejecucion del usuario

    parser.add_argument('Analysis_Types', 
                        help = 'Enter the name of line to analyzed',
                        type= str,
                        choices=Analysis_Types) #Agrega argumentos para ejecucion del usuario
    
    args = parser.parse_args()
    Analysis(args.Type_width, args.Analysis_Types)


