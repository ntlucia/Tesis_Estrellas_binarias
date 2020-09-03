#!/usr/bin/python

import numpy
import readData_v2
import astropy.io.fits as pyfits
import glob

print("***************************************************")
print("*                                                 *")
print("*  Convert the TIGRE spectrum into an ASCII file  *")
print("*                                                 *")
print("***************************************************")
print("")

Fits=glob.glob('*.fits') 	

for file in Fits:
    print("TIGRE Spektrum: "),
    print(file)
        
    print("Output file: ")
    out_file = file.replace('Sci_', '').replace('.sp_ech_main.fits', '') + '.dat'
        
    spec = readData_v2.readheros(file)

    him = spec[0]
    hs = spec[1]
    fields = hs['TFIELDS']
    if fields == 12:
            wave = spec[2]
            spec_blz = spec[3]
            blaze = spec[4]
            sig = spec[5]
            spec_norm = spec[6]
            wave_merg = spec[7]
            spec_merg = spec[8]
            f_spec_merg = spec[9]
            spec_merg_norm = spec[10]
            f_spec_merg_norm = spec[11]
            irf = spec[12]
            day_var = spec[13]
                
    if fields == 13:
            wave = spec[2]
            spec_blz = spec[3]
            blaze = spec[4]
            sig = spec[5]
            spec_norm = spec[6]
            single_spec = spec[7]
            wave_merg = spec[8]
            spec_merg = spec[9]
            f_spec_merg = spec[10]
            spec_merg_norm = spec[11]
            f_spec_merg_norm = spec[12]
            irf = spec[13]
            day_var = spec[14]
                
    # correct and save spectrum
    nwave = len(wave_merg)
    tab = numpy.zeros((nwave,3), float)
    tab[::,0] = wave_merg/10.0
    tab[::,1] = spec_merg_norm
    tab[::,2] = f_spec_merg_norm

    numpy.savetxt(out_file, tab)
    print("spectrum written to",out_file)
