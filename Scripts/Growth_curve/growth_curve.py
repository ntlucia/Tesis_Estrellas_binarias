#// ============================================================================================= //
#//                                                                                               //
#//       Filename:  growth_curve.py                                                              //
#//    Description:  Equivalent width and loggf program for cromospheric lines                    //
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
#//                                                                                               //
#//                                                                                               //
#// ============================================================================================= //


import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt
from heapq import nsmallest
from scipy.signal import butter, lfilter, freqz
from scipy.fftpack import fft, ifft
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import CubicSpline , interp1d
from common import config
import argparse  #Para crear argumentos en el ejutable


#---------- FOUND LOGGF AND ADD TO LIST ----------#

def elements_analyze(lines_s, Atomic_lines, element):
    data_element = lines_s[lines_s['element'] == element]
    data_element.index = list(range(len(data_element)))

    #Separate VALD lines by elements for decrease the amount of data
    val = Atomic_lines[Atomic_lines['element']== element]
    
    
    #Used just workin data range
    atomic_data_element = val[(val['wave_A'] >= 376) & (val['wave_A'] <= 460)]
    atomic_data_element.index = list(range(len(atomic_data_element)))

    #Working in nm
    return [data_element, atomic_data_element]

def list_loggf(element, data_element, atomic_data_element,spectrum, theoric_abs_lines):
    element_ = theoric_abs_lines[theoric_abs_lines['element']== element]
    element_.index = list(range(len(element_)))

    new_listgf = []
    loggf = []
    for i in range(len(element_['waveobs'])):
        x = abs(nsmallest(1,atomic_data_element['wave_A'], key = lambda x: abs(x-element_['waveobs'][i]))[0] - element_['waveobs'][i])
        if x < 0.1:
            new_listgf.append(nsmallest(1,atomic_data_element['wave_A'], key = lambda x: abs(x-element_['waveobs'][i]))[0])
            loggf.append(atomic_data_element[atomic_data_element["wave_A"] == new_listgf[i]]["loggf"].tolist()[0])
    Ti2_gf = pd.DataFrame(columns = ['waveobs','element', 'flux', 'loggf', 'wave_base', 'wave_top', 'error_f'])
    Ti2_gf['waveobs'] = data_element['waveobs']
    Ti2_gf['element'] = data_element['element']
    Ti2_gf['flux'] = data_element['flux']
    Ti2_gf['loggf'] = loggf
    Ti2_gf['wave_base'] = data_element['wave_base']
    Ti2_gf['wave_top'] = data_element['wave_top']
    Ti2_gf['error_f'] = spectrum['err']
    return Ti2_gf


#---------- CALCULATE EQUIVALENT WIDTH AND ADD TO LIST ----------#

#FIT --------- HACERLO CON ISPEC
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


#Find the local continuum for each line
def pseudocontinuou(Spectrum,LMIN, LMAX):
    pseudo_continuous1 = Spectrum[(Spectrum['waveobs'] <= LMIN)  & (Spectrum['waveobs'] > (LMIN-0.05))]
    pseudo_continuous2 = Spectrum[(Spectrum['waveobs'] >= LMAX )  & (Spectrum['waveobs'] < (LMAX+0.05))]
    pseudo_continuous1.index = list(range(len(pseudo_continuous1)))
    pseudo_continuous2.index = list(range(len(pseudo_continuous2)))

    mean1 = np.mean(pseudo_continuous1)
    mean2 = np.mean(pseudo_continuous2)
    mean = (mean1['flux'] + mean2['flux'])/2
    return mean

#Calculate equivalent width
def Equivalent_width(fit, spectrum, mean, base, top):
    Area_rec = mean*(top - base)
    #Area_rec = mean*( fit['L'][top] - fit['L'][base])

    b = spectrum['waveobs']>base
    a = spectrum['waveobs']<top
    c = a&b

    Area_fit = integrate.simps(spectrum['flux'][c],spectrum["waveobs"][c])
    #Area_fit = integrate.simps(fit['I'], fit['L'])
    Area_real = Area_rec - Area_fit
    EW = Area_real/mean
    return EW

def Equivalent_width_comp(spectrum, base,top):
    equivalent_widths = []
    for i in range(len(base)):
        step = step_continuum( base[i] ,top[i],20)
        grouped_data = Points_continuum(spectrum["waveobs"], spectrum["flux"], base[i], top[i], step)
        fit= fit_continuum(grouped_data, base[i], top[i])
        mean = pseudocontinuou(spectrum,base[i], top[i])
        equivalent_width = Equivalent_width(fit, spectrum, mean, base[i], top[i])
        equivalent_widths.append(equivalent_width)
    print(len(equivalent_widths))
    return equivalent_widths


def data_growth_curve(_data_lines, equivalent_widths):
    element_growth_c = pd.DataFrame(columns = ['waveobs','element', 'flux', 'loggf', 'EW','wave_base', 'wave_top', 'error_f'])
    element_growth_c['waveobs'] = _data_lines['waveobs']
    element_growth_c['element'] = _data_lines['element']
    element_growth_c['flux'] = _data_lines['flux']
    element_growth_c['loggf'] = _data_lines['loggf']
    element_growth_c['EW'] = equivalent_widths
    element_growth_c['wave_base'] = _data_lines['wave_base']
    element_growth_c['wave_top'] = _data_lines['wave_top']
    element_growth_c['error_f'] = _data_lines['error_f']
    return element_growth_c


def growth_curve(n_CROMOSP_LINES, n_ATOMIC_LINES, n_SPECTRUM, n_THEORIC_CROMOSP_LINES,  element):
    
    #Lecture Data: klaus' lines, cromospheric spectrum, VALD
    spectrum =  pd.read_csv(n_SPECTRUM, delimiter = '\t', header = 0)
    lines = pd.read_csv(n_CROMOSP_LINES, delimiter = ',', header = 0)
    Atomic_lines = pd.read_csv(n_ATOMIC_LINES, delimiter = '\t', usecols = ['element', 'wave_A','loggf'], header = 0, low_memory=False, keep_default_na= False)
    Atomic_lines['wave_A'] = Atomic_lines['wave_A']/10
    lines_ = pd.read_excel("lines.xlsx",sheet_name="cromospheric_lines", columns = ['waveobs', 'element','wave_base', 'wave_top'] )
    lines_['waveobs'] = lines_['waveobs']/10
    lines_['wave_base'] = lines_['wave_base']/10
    lines_['wave_top'] = lines_['wave_top']/10


    [data_element,atomic_data_element] = elements_analyze(lines, Atomic_lines, element)
    _data_lines = list_loggf(element, data_element, atomic_data_element,spectrum, lines_)
    

    equivalent_widths = Equivalent_width_comp(spectrum, _data_lines['wave_base'], _data_lines['wave_top'])
    data_lines = data_growth_curve(_data_lines, equivalent_widths)
    data_lines.to_csv("{}_growth_curve.csv".format(element), index = False, header=True)


if __name__ == '__main__':
    _config = list(config().keys()) 

    parser = argparse.ArgumentParser()
    parser.add_argument('element', 
                        help = 'Elements to Analyzer',
                        type= str) #Agrega argumentos para ejecucion del usuario
    
    args = parser.parse_args()
    growth_curve( config()["CROMOSP_LINES"],
                  config()["ATOMIC_LINES"] ,
                  config()["SPECTRUM"]     ,
                  config()["THEORIC_CROMOSP_LINES"],
                  args.element             )
                  