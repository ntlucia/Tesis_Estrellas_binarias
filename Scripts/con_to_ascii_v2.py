#!/usr/bin/python


import numpy
import readData_v2
import pyfits

print "***************************************************"
print "*                                                 *"
print "*  Convert the TIGRE spectrum into an ASCII file  *"
print "*                                                 *"
print "***************************************************"
print ""

print "TIGRE Spektrum: ",
file = raw_input()
print "Output file: ",
out_file = raw_input()
 
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

print ""
count = 0
while (count < 3):
    
    print ""
    print "If you want a single order (y/n): ",
    question_order = raw_input()
    
    
    if question_order == 'y':
        print "Select order: ",
        order_num = raw_input()
        
        nwave = len(wave[::,order_num])
        tab = numpy.zeros((nwave,5), float)
        tab[::,0] = wave[::,order_num]
        tab[::,1] = spec_blz[::,order_num]
        tab[::,2] = blaze[::,order_num]
        tab[::,3] = sig[::,order_num]
        tab[::,4] = spec_norm[::,order_num]
            
        numpy.savetxt(out_file, tab)

        count=4
        
    elif question_order == 'n':
        nwave = len(wave_merg)
        tab = numpy.zeros((nwave,7), float)
        tab[::,0] = wave_merg
        tab[::,1] = spec_merg
        tab[::,2] = f_spec_merg
        tab[::,3] = spec_merg_norm
        tab[::,4] = f_spec_merg_norm
        tab[::,5] = irf
        tab[::,6] = day_var

        numpy.savetxt(out_file, tab)
        
        count=4
            
    else:
        print ""
        print "You have to use y or n!"
        print ""
        
        count = count+1


