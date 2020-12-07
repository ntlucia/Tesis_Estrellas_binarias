#// ============================================================================================= //
#//                                                                                               //
#//       Filename:  growth_curve.py                                                              //
#//    Description:  Equivalent width and loggf program for cromospheric lines                    //
#//                                                                                               //
#//        Version:  7.1                                                                          //
#//        Created:  25/07/2020                                                                   //
#//       Compiler:  Python                                                                       //
#//                                                                                               //
#//         Author:  Natalia Lucia Oliveros Gomez                                                 //
#//          Email:  onatalialucia@gmail.com                                                      //
#//        Company:  Grupo Halley UIS - Grupo Fisica Estelar U. Gto                               //
#//                                                                                               //
#// ============================================================================================= //
#
#// ============================================================================================= //
#//                                                                                               // 
#//  Compile with:  python growth_curve.py 'element'                                              //
#//                                                                                               //
#//                                                                                               //
#// ============================================================================================= //


import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt
from heapq import nsmallest
from scipy import integrate
from common import config
import argparse  #Para crear argumentos en el ejutable
import logging



# --- Change LOG level --- #

#LOG_LEVEL = "warning"
LOG_LEVEL = "info"
logger = logging.getLogger() # root logger, common for all
logger.setLevel(logging.getLevelName(LOG_LEVEL.upper()))



#---------- ATOMIC DATA, THEORETICAL LINE AND SPECTRUM MATCHES ----------#

def match_theoric_atomic_data(atomic_lines, theoric_lines, element):

    logging.info("Restricting to the specified element")
    #Filter according to the element to be analyzed
    data_element = theoric_lines[theoric_lines['note'] == element]
    data_element.index = list(range(len(data_element)))
    
    atomic_data_element = atomic_lines[atomic_lines['element'] == element]
    atomic_data_element.index = list(range(len(atomic_data_element)))


    logging.info("Atomic data and theoretical line matches")

    #Find the nearest value within the list of atomic data
    peak = []
    for i in range(len(data_element['wave_peak'])):
        peak.append(nsmallest(1,atomic_data_element['wave_nm'], key = lambda x: abs(x-data_element['wave_peak'][i]))[0])

    #Filter only the lines that interest me (theorical lines)
    Atomic_lines = []
    for i in range(len(peak)):
        Atomic_lines.append(atomic_data_element[(atomic_data_element['wave_nm'] == peak[i])])

    logging.info("Creating a new match list")

    #Filter only the lines that interest me (theorical lines) 
    #columns = ['wave_peak', 'wave_base', 'wave_top','note', 'loggf', 'lower_state_eV', 'upper_state_eV']
    match_atomic_theoric = pd.DataFrame(columns = ['wave_peak', 'wave_base', 'wave_top','note', 'loggf', 'lower_state_eV', 'upper_state_eV'])
    match_atomic_theoric['wave_base'] = data_element['wave_base']
    match_atomic_theoric['wave_top'] = data_element['wave_top']

    for i in range(len(Atomic_lines)):
        match_atomic_theoric['wave_peak'].loc[i] = Atomic_lines[i]['wave_nm'].values[0]
        match_atomic_theoric['note'].loc[i] = Atomic_lines[i]['element'].values[0]
        match_atomic_theoric['loggf'].loc[i] = Atomic_lines[i]['loggf'].values[0]
        match_atomic_theoric['lower_state_eV'].loc[i] = Atomic_lines[i]['lower_state_eV'].values[0]
        match_atomic_theoric['upper_state_eV'].loc[i] = Atomic_lines[i]['upper_state_eV'].values[0]

    return match_atomic_theoric


def match_spectrum(spectrum, match_atomic_theoric):

    logging.info("Matches with star spectrum data")

    #Find the closest values within the data spectrum for the ends of the line
    top = []
    for i in range(len(match_atomic_theoric['wave_top'])):
        top.append(nsmallest(1,spectrum['waveobs'], key = lambda x: abs(x-match_atomic_theoric['wave_top'][i]))[0])

    base = []
    for i in range(len(match_atomic_theoric['wave_base'])):
        base.append(nsmallest(1,spectrum['waveobs'], key = lambda x: abs(x-match_atomic_theoric['wave_base'][i]))[0])

    logging.info("Creating a new match list")

    c_lines_spectrum = pd.DataFrame(columns = ['wave_peak', 'wave_base', 'wave_top','note', 'loggf', 'lower_state_eV', 'upper_state_eV','flux', 'error_f'])
    L = []
    I = []
    mins = []
       
    for i in range(len(base)):
        #Find the actual peak, the minimum between the ends of the line
        b = spectrum['waveobs']>base[i]
        a = spectrum['waveobs']<top[i]
        c = a&b 
        min_ = min(spectrum['flux'][c])
        L.append(spectrum["waveobs"][c][spectrum["flux"][c] == min_].tolist()[0])
        I.append(spectrum['flux'][c][spectrum["flux"][c] == min_].tolist()[0])
        mins.append(min_)
    
    c_lines_spectrum['wave_peak'] = L
    c_lines_spectrum['wave_base'] = base
    c_lines_spectrum['wave_top'] = top
    c_lines_spectrum['note'] = match_atomic_theoric['note']
    c_lines_spectrum['loggf'] = match_atomic_theoric['loggf']
    c_lines_spectrum['lower_state_eV'] = match_atomic_theoric['lower_state_eV']
    c_lines_spectrum['upper_state_eV'] = match_atomic_theoric['upper_state_eV']
    c_lines_spectrum['flux'] = I
    c_lines_spectrum['error_f'] = spectrum['err']

    return c_lines_spectrum




#---------- CALCULATE EQUIVALENT WIDTH AND ADD TO LIST ----------#


#Calculate equivalent width
def Equivalent_width(spectrum, peak, base, top,i):
    logging.info("Calculating the local continuum")
    pseudo_continuous1 = spectrum[(spectrum['waveobs'] <= base)  & (spectrum['waveobs'] > (base-0.06))] #Left
    pseudo_continuous2 = spectrum[(spectrum['waveobs'] >= top )  & (spectrum['waveobs'] < (top+0.06))] #Right
    pseudo_continuous1.index = list(range(len(pseudo_continuous1)))
    pseudo_continuous2.index = list(range(len(pseudo_continuous2)))

    mean1 = np.mean(pseudo_continuous1)
    mean2 = np.mean(pseudo_continuous2)
    mean = (mean1['flux'] + mean2['flux'])/2

    logging.info("Calculating areas under the curve")
    #Area under the local continuum
    Area_rec = mean*(top - base)

    b = spectrum['waveobs']>base
    a = spectrum['waveobs']<top
    c = a&b
    
    #Area under the absortion line 
    Area_fit = integrate.simps(spectrum['flux'][c],spectrum["waveobs"][c])
    #Area within the absorption line
    Area_real = Area_rec - Area_fit

    logging.info("Calculating equivalent width")

    #Equivalent width
    EW = Area_real/mean
    EWR = np.log10(abs(EW/peak))

    logging.info("Calculating error propagation")

    errlambda = 0.03 #Instrumental error
    errflux = (abs(mean - mean1['flux']) + abs(mean - mean2['flux']))/2 #flux error
    errorA1 = Area_rec*(errflux/mean)*((errlambda + errlambda)/(top - base))
    #errorA2 = np.sqrt(np.sum(spectrum['err'][c]**2))
    #errorAreal = errorA1 + errorA2 #Area error
    errEW = EW*(errorA1/Area_real + errflux/mean) #Equivalent width error
    errEWR = errEW/(EW*np.log(10))
    return EW,EWR, errEWR



    


def Equivalent_width_comp(spectrum,_data_lines, base,top):

    logging.info("Calculating equivalent widths, for all lines of the specified element")

    equivalent_widths = []
    equivalent_widths_r = []
    errorsEWR = []
    for i in range(len(base)):
        try:
            EW,EWR, errEWR = Equivalent_width(spectrum, _data_lines[i], base[i], top[i],i)
            equivalent_widths.append(EW)
            equivalent_widths_r.append(EWR)
            errorsEWR.append(errEWR)
        except IndexError:
            print("IndexError, i={}".format(i))

        
    return [equivalent_widths,equivalent_widths_r, errorsEWR]

#Create final list, with which the growth curves can be plotted
def data_growth_curve(_data_lines, equivalent_widths,equivalent_widths_r, errorsEWR ):

    logging.info("Creating final list")
    
    element_growth_c = pd.DataFrame(columns = ['wave_peak','wave_base', 'wave_top','note', 'flux', 'loggf', 'lower_state_eV', 'upper_state_eV', 'EW', 'EWR', 'error_f'])
    element_growth_c['wave_peak'] = _data_lines['wave_peak']
    element_growth_c['wave_base'] = _data_lines['wave_base']
    element_growth_c['wave_top'] = _data_lines['wave_top']
    element_growth_c['note'] = _data_lines['note']
    element_growth_c['flux'] = _data_lines['flux']
    element_growth_c['loggf'] = _data_lines['loggf']
    element_growth_c['lower_state_eV'] = _data_lines['lower_state_eV']
    element_growth_c['upper_state_eV'] = _data_lines['upper_state_eV']
    element_growth_c['EW'] = equivalent_widths
    element_growth_c['EWR'] = equivalent_widths_r
    element_growth_c['errEWR'] = errorsEWR
    element_growth_c['error_f'] = _data_lines['error_f']
    return element_growth_c



def growth_curve(n_CROMOSP_LINES, n_ATOMIC_LINES, n_SPECTRUM, n_THEORIC_CROMOSP_LINES,  element):
    
    #Lecture Data: cromospheric spectrum, VALD and theoretical  lines
    logging.info("Reading spectrum")
    spectrum =  pd.read_csv(n_SPECTRUM, delimiter = '\t', header = 0)

    logging.info("Reading atomic data")
    Atomic_lines = pd.read_csv(n_ATOMIC_LINES, delimiter = '\t', usecols = ['element', 'wave_nm','loggf', 'lower_state_eV', 'upper_state_eV'], header = 0, low_memory=False, keep_default_na= False)
    
    logging.info("Reading theoretical data of chromospheric lines")
    lines_ = pd.read_csv(n_THEORIC_CROMOSP_LINES, delimiter = '\t', header = 0)


    #Execute all functions
    match_theoric_atomic = match_theoric_atomic_data(Atomic_lines, lines_, element)

    c_lines_spectrum = match_spectrum(spectrum, match_theoric_atomic)

    
    equivalent_widths,equivalent_widths_r, error_ew = Equivalent_width_comp(spectrum,c_lines_spectrum['wave_peak'], c_lines_spectrum['wave_base'], c_lines_spectrum['wave_top'])
    data_lines = data_growth_curve(c_lines_spectrum, equivalent_widths, equivalent_widths_r, error_ew )

    data_lines.to_csv("DataSet/Outputs/{}_growth_curve.dat".format(element), sep = '\t', index = False, header=True)


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
                  