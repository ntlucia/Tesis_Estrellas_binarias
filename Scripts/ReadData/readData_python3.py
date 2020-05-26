

# This is a collection of functions to use the HRT spectra.



# written by M. Mittag


import gzip
import numpy
import scipy
import matplotlib.pylab as mpl
import math
import pyfits
import scipy.optimize


def readheros(filename, radvel=None, barycorr=None, small=None):

# Function to read the HRT main and quicklook spectum

	if small is None:
	#	print radvel	
		spec = pyfits.open(filename)
		headim = spec[0].header
		headspec = spec[1].header

		if barycorr == 'no':
			barycorr = 0.
		else:
			barycorr = headim['BARYCORR']

		if radvel is None:
			radvel = 0.

		fields = headspec['TFIELDS']
		if fields == 12:
			dim = headspec['TDIM1']
			dim_n1 = dim.replace( '(' ,'')
			dim_n2 = dim_n1.replace( ')' ,'') 
			dim_n3 = dim_n2.replace( ',' ,'') 
			dim_f = dim_n3.split()
			npix = int(dim_f[0])
			nord = int(dim_f[1])
			wave = numpy.zeros((npix,nord),float)
			spec_blz = numpy.zeros((npix,nord),float)
			blaze = numpy.zeros((npix,nord),float)
			sig = numpy.zeros((npix,nord),float)
			spec_norm = numpy.zeros((npix,nord),float)
			d = spec[1].data.field('WAVE')
			d = numpy.reshape(d, (nord,npix))
			for i in range(nord):
				wave[::,i] = d[i,::]
			d = spec[1].data.field('SPEC_BLZ')
			d = numpy.reshape(d, (nord,npix))
			for i in range(nord):
				spec_blz[::,i] = d[i,::]
			d = spec[1].data.field('BLAZE')
			d= numpy.reshape(d, (nord,npix))
			for i in range(nord):
				blaze[::,i] = d[i,::]
			d = spec[1].data.field('SIG')
			d = numpy.reshape(d, (nord,npix))
			for i in range(nord):
				sig[::,i] = d[i,::]
			d = spec[1].data.field('SPEC_NORM')
			d = numpy.reshape(d, (nord,npix))
			for i in range(nord):
				spec_norm[::,i] = d[i,::]
			wave_merg = spec[1].data.field('WAVE_MERG')[0,::]
			spec_merg = spec[1].data.field('SPEC_MERG')[0,::]
			f_spec_merg = spec[1].data.field('F_SPEC_MERG')[0,::]
			spec_merg_norm = spec[1].data.field('SPEC_MERG_NORM')[0,::]
			f_spec_merg_norm = spec[1].data.field('F_SPEC_MERG_NORM')[0,::]
			irf = spec[1].data.field('INST_RESP_FUNC')[0,::]
			day_var = spec[1].data.field('RESPONSE_CURVE_DAY_VAR')[0,::]

			if barycorr == 'no':	
				return headim, headspec, wave, spec_blz, blaze, sig, spec_norm, wave_merg, spec_merg, f_spec_merg, spec_merg_norm, f_spec_merg_norm, irf, day_var
			else:
				wave = wave * (1.0 + barycorr / (2.9979246e5))
				wave_merg = wave_merg * (1.0 + barycorr / (2.9979246e5))
				
				wave = wave * (1.0 - radvel / (2.9979246e5))
				wave_merg = wave_merg * (1.0 - radvel / (2.9979246e5))

				return headim,headspec , wave, spec_blz, blaze, sig, spec_norm, wave_merg, spec_merg, f_spec_merg, spec_merg_norm, f_spec_merg_norm, irf, day_var

		if fields == 13:
			dim = headspec['TDIM6']
			dim_n1 = dim.replace( '(' ,'')
			dim_n2 = dim_n1.replace( ')' ,'') 
			dim_n3 = dim_n2.replace( ',' ,'') 
			dim_f = dim_n3.split()
			nspec = int(dim_f[0])
			npix = int(dim_f[1])
			nord = int(dim_f[2])
			print(dim)
			print(nord)
			print(nspec)
			print(npix)
			wave = numpy.zeros((npix,nord),float)
			spec_blz = numpy.zeros((npix,nord),float)
			blaze = numpy.zeros((npix,nord),float)
			sig = numpy.zeros((npix,nord),float)
			spec_norm = numpy.zeros((npix,nord),float)
			single_spec = numpy.zeros((nspec,npix,nord),float)
			d = spec[1].data.field('WAVE')
			d = numpy.reshape(d, (nord,npix))
			for i in range(nord):
				wave[::,i] = d[i,::]
			d = spec[1].data.field('SPEC_BLZ')
			d = numpy.reshape(d, (nord,npix))
			for i in range(nord):
				spec_blz[::,i] = d[i,::]
			d = spec[1].data.field('BLAZE')
			d= numpy.reshape(d, (nord,npix))
			for i in range(nord):
				blaze[::,i] = d[i,::]
			d = spec[1].data.field('SIG')
			d = numpy.reshape(d, (nord,npix))
			for i in range(nord):
				sig[::,i] = d[i,::]
			d = spec[1].data.field('SPEC_NORM')
			d = numpy.reshape(d, (nord,npix))
			for i in range(nord):
				spec_norm[::,i] = d[i,::]
			d = spec[1].data.field('SINGLE_OBS_SPEC')
			d = numpy.reshape(d, (nord,npix,nspec))
			for j in range(nspec):
				for i in range(nord):
					single_spec[j,::,i] = d[i,::,j]
			wave_merg = spec[1].data.field('WAVE_MERG')[0,::]
			spec_merg = spec[1].data.field('SPEC_MERG')[0,::]
			f_spec_merg = spec[1].data.field('F_SPEC_MERG')[0,::]
			spec_merg_norm = spec[1].data.field('SPEC_MERG_NORM')[0,::]
			f_spec_merg_norm = spec[1].data.field('F_SPEC_MERG_NORM')[0,::]
			irf = spec[1].data.field('INST_RESP_FUNC')[0,::]
			day_var = spec[1].data.field('RESPONSE_CURVE_DAY_VAR')[0,::]


			if barycorr == 'no':	
				return headim, headspec, wave, spec_blz, blaze, sig, spec_norm, single_spec,wave_merg, spec_merg, f_spec_merg, spec_merg_norm, f_spec_merg_norm, irf, day_var
			else:
				wave = wave * (1.0 + barycorr / (2.9979246e5))
				wave_merg = wave_merg * (1.0 + barycorr / (2.9979246e5))
				
				wave = wave * (1.0 - radvel / (2.9979246e5))
				wave_merg = wave_merg * (1.0 - radvel / (2.9979246e5))

				return headim, headspec, wave, spec_blz, blaze, sig, spec_norm, single_spec,wave_merg, spec_merg, f_spec_merg, spec_merg_norm, f_spec_merg_norm, irf, day_var

	if small == 'yes':
		spec = pyfits.open(filename)
		headspec = spec[0].header

		if barycorr == 'no':
			barycorr = 0.
		else:
			barycorr = headspec['BARYCORR']

		if radvel is None:
			radvel = 0.
	#	print radvel	

		npix = headspec['NAXIS1']
		w0 = headspec['CRVAL1']
		dw = headspec['CDELT1']

		wave = w0 + numpy.arange(npix) * dw

		spec = spec[0].data
		
		if barycorr == 'no':	
			return headspec, wave, spec
		else:
			wave = wave * (1.0 + barycorr / (2.9979246e5))
			
			wave = wave * (1.0 - radvel / (2.9979246e5))

			return headspec, wave, spec




def ReadSpecHeaderOnly(filename, small=None):

# Function to read only the header of the HRT sprectra

        if small is None:
            #	print radvel	
                spec = pyfits.open(filename)
                headim = spec[0].header
                headspec = spec[1].header

                return headspec, headim

        if small == 'yes':
            spec = pyfits.open(filename)
            headspec = spec[0].header

            return headspec



