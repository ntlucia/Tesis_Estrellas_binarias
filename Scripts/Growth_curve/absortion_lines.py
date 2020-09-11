import pandas as pd
import numpy as np
import scipy

def absortion_lines(lines_, spectrum):
    #found point that correspond to absortion lines in cromospheric spectrum
    new_list = []
    for i in range(len(lines_['waveobs'])):
        new_list.append(nsmallest(1,spectrum['waveobs'], key = lambda x: abs(x-lines_['waveobs'][i]))[0])

    top = []
    for i in range(len(lines_['wave_top'])):
        top.append(nsmallest(1,spectrum['waveobs'], key = lambda x: abs(x-lines_['wave_top'][i]))[0])

    base = []
    for i in range(len(lines_['wave_base'])):
        base.append(nsmallest(1,spectrum['waveobs'], key = lambda x: abs(x-lines_['wave_base'][i]))[0])

    #Create new list with contain the absortion lines in cromospheric spectrum with its respective element
    #Create new list with contain the absortion lines in cromospheric spectrum with its respective element
    lines_s = pd.DataFrame(columns = ['waveobs','element', 'flux', 'wave_base', 'wave_top', 'error_f'])
    L = []
    I = []
    mins = []
                
    for i in range(len(base)):
        b = spectrum['waveobs']>base[i]
        a = spectrum['waveobs']<top[i]
        c = a & b 
        min_ = min(spectrum['flux'][c])
        L.append(spectrum["waveobs"][c][spectrum["flux"][c] == min_].tolist()[0])
        I.append(spectrum['flux'][c][spectrum["flux"][c] == min_].tolist()[0])
        mins.append(min_)
        
    lines_s['waveobs'] = L
    lines_s['flux'] = I
    lines_s['element'] = lines_['element']
    lines_s['wave_base'] = base[0:len(base)]
    lines_s['wave_top'] = top[0:len(top)]
    lines_s['error_f'] = spectrum['err']
    return lines_s



    def growth_curve(n_CROMOSP_LINES, n_SPECTRUM):
        #Lecture Data: klaus' lines, cromospheric spectrum, VALD
        spectrum = pd.read_csv(n_SPECTRUM, delimiter = '\t', header = 0)
        #Atomic_lines = pd.read_csv(n_ATOMIC_LINES, delimiter = '\t', usecols = ['element', 'wave_A','loggf'], header = 0, low_memory=False, keep_default_na= False)
        lines_ =  pd.read_excel(n_CROMOSP_LINES, sheet_name="cromospheric_lines", columns = ['waveobs', 'element','wave_base', 'wave_top'] )
        lines_['waveobs'] = lines_['waveobs']/10
        lines_['wave_base'] = lines_['wave_base']/10
        lines_['wave_top'] = lines_['wave_top']/10
    
        lines = absortion_lines(lines_,spectrum)
    
        lines.to_csv("{}_growth_curve.csv".format(lines), index = False, header=True)


if __name__ == '__main__':
    _config = list(config().keys()) 

    args = parser.parse_args()
    growth_curve( config()["CROMOSP_LINES"],
                  config()["ATOMIC_LINES"] ,
                  config()["SPECTRUM"]            )