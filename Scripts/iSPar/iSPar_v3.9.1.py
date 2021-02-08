#!/usr/bin/env python

#v3.9.1

import os
import sys
import math
import numpy as np
import logging
import os.path
from os import path
import pandas as pd
from iSPar_Config import * # Load configuration file

################################################################################

# --- iSpec directory ---

#ispec_dir = '/home/executus/Software/iSpec/'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec

# --- Change LOG level ---

#LOG_LEVEL = "warning"
LOG_LEVEL = "info"
logger = logging.getLogger() # root logger, common for all
logger.setLevel(logging.getLevelName(LOG_LEVEL.upper()))

################################################################################

def frange(start, stop, step):
    x = start
    while x < stop+step:
        yield x
        x += step

# End of frange

################################################################################
        
def spectrum_reduction(spectrum_filename, results_dirname, resolution, m_spectrum_filename):
    
    # --- Read spectrum ---
    
    logging.info("Reading spectrum")
    
    star_spectrum = ispec.read_spectrum(spectrum_filename)
    
    # --- Radial Velocity determination with linelist mask ---
    
    logging.info("Radial velocity determination with linelist mask...")
    
    ccf_mask = ispec.read_cross_correlation_mask(mask_file)
    
    models, ccf = ispec.cross_correlate_with_mask(star_spectrum, ccf_mask, \
                                                  lower_velocity_limit=-200, upper_velocity_limit=200, \
                                                  velocity_step=0.5, mask_depth=0.01, \
                                                  fourier=False)
    
    # Number of models represent the number of components
    
    components = len(models)
    
    # First component:
    
    rv = np.round(models[0].mu(), 2) # km/s
    rv_err = np.round(models[0].emu(), 2) # km/s
    
    # --- Radial Velocity correction ---
    
    logging.info("Radial velocity correction: %.2f +/- %.2f" % (rv, rv_err))
    m_star_spectrum = ispec.correct_velocity(star_spectrum, rv)
    
    # --- Cut ---
    
    logging.info("Cutting...")
    
    # Keep points between two given wavelengths
    
    wfilter = ispec.create_wavelength_filter(m_star_spectrum, wave_base, wave_top)
    m_star_spectrum = m_star_spectrum[wfilter]

    if resolution_degradation == True:

        # --- Resolution degradation ---
    
        logging.info("Resolution degradation...")
        
        # NOTE: The line selection was built based on a solar and stellar spectrum with R ~ 20,000 and VALD atomic linelist.

        from_resolution = resolution
        to_resolution = final_R
        m_star_spectrum = ispec.convolve_spectrum(m_star_spectrum, to_resolution, from_resolution)

    # --- Resample spectrum ---
    
    logging.info("Resampling spectrum...")

    wave_step = np.median(np.abs(m_star_spectrum['waveobs'][1:] - m_star_spectrum['waveobs'][:-1]))
    wavelengths = np.arange(wave_base, wave_top, wave_step)
    m_star_spectrum = ispec.resample_spectrum(m_star_spectrum, wavelengths, method="bessel", zero_edges=False)

    logging.info("Resampled to wave step: %.5f" % wave_step)
    
    # Cut bad edges
    
    wfilter = ispec.create_wavelength_filter(m_star_spectrum, wave_base + wave_step, wave_top - wave_step)
    m_star_spectrum = m_star_spectrum[wfilter]
        
    # --- Save modified (RV + cut) spectrum ---
    
    logging.info("Saving modified (RV + cut) spectrum...")
    ispec.write_spectrum(m_star_spectrum, "./%s/%s" % (results_dirname, m_spectrum_filename))
    
    return wave_step, m_star_spectrum

# End of spectrum_reduction

################################################################################

def precompute_synthetic_grid(precomputed_grid_dirname, modeled_layers_pack, atomic_linelist, isotopes, solar_abundances):
    
    precomputed_grid_dir = ispec_dir + "/input/minigrid/%s/" % precomputed_grid_dirname
    
    logging.info("The grid will be saved in: %s" % precomputed_grid_dir)
    
    # --- Read grid ranges from file ---
    
    from astropy.io import ascii
    ranges_file = ispec_dir + "/input/minigrid/initial_estimate_grid_ranges.tsv"
    
    logging.info("Reading grid ranges from file: %s" % ranges_file)
    
    ranges = ascii.read(ranges_file, delimiter="\t")
    
    # Wavelengths
    
    wave_base = 480
    wave_top = 680
    
    resolution = 20000 # Aprox. to TIGRE Resolution. Individual files will not be convolved but the grid will be (for fast comparison)
    
    wavelengths = np.arange(wave_base, wave_top, wave_step)
    
    number_of_processes = 4 # It can be parallelized for computers with multiple processors
    
    # --- Precompute synthetic grid for code ---
    
    logging.info("Precomputing synthetic grid for: %s" % code)
    
    ispec.precompute_synthetic_grid(precomputed_grid_dir, ranges, wavelengths, resolution, \
                                    modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, \
                                    segments=None, number_of_processes=number_of_processes, \
                                    code=code, steps=False)
    
    logging.info("The synthetic grid for %s was successfully created!" % code)

# End of precompute_synthetic_grid
    
################################################################################

def estimate_initial_parameters_with_precomputed_grid(precomputed_grid_dirname, m_star_spectrum, initial_R, vel_telluric, telluric_linelist, normalization=True):
    
    ############################################################################
    # WARNING!!!
    # This routine depends on the previous precomputation of the synthetic grid
    ############################################################################
    
    if os.path.exists(ispec_dir + "/input/minigrid/%s/" % precomputed_grid_dirname):
        
        precomputed_grid_dir = ispec_dir + "/input/minigrid/%s/" % precomputed_grid_dirname
        
    else:
        
        logging.info("ERROR: Please calculate a synthetic grid first!")
        sys.exit()
    
    # --- Cut spectrum acording to the synthetic grid ---
        
    logging.info("Cutting spectrum acording to the synthetic grid...")
    
    wave_base = 480.0
    wave_top = 680.0
    
    wfilter = ispec.create_wavelength_filter(m_star_spectrum, wave_base, wave_top)
    m_star_spectrum = m_star_spectrum[wfilter]    

    # --- Convolve spectrum ---
    
    from_resolution = initial_R
    to_resolution = 20000
    
    logging.info("Convolving spectrum to grid resolution...")
    
    m_star_spectrum = ispec.convolve_spectrum(m_star_spectrum, to_resolution, from_resolution)

    if normalization == True:
        
        # --- Continuum fit ---
        
        fit_model = "Splines" # "Polynomy"
        degree = 3 # 2, 3
        nknots = np.round((wave_top-wave_base)/2, 0) # None: 1 spline every 5 nm
        
        # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
        
        order='median+max'
        median_wave_range = 0.100
        max_wave_range = 0.50
        
        logging.info("Fiting continuum...")
        
        star_continuum_model = ispec.fit_continuum(m_star_spectrum, \
                                                   ignore=strong_lines, \
                                                   nknots=nknots, degree=degree, \
                                                   median_wave_range=median_wave_range, \
                                                   max_wave_range=max_wave_range, \
                                                   model=fit_model, order=order, \
                                                   automatic_strong_line_detection=True, \
                                                   strong_line_probability=0.70, \
                                                   use_errors_for_fitting=True)
        
        # --- Continuum normalization ---
        
        logging.info("Normalizing spectrum...")
        
        n_star_spectrum = ispec.normalize_spectrum(m_star_spectrum, star_continuum_model, consider_continuum_errors=False)    

    elif normalization == False:

        n_star_spectrum = m_star_spectrum
        
    n_star_continuum_model = ispec.fit_continuum(n_star_spectrum, fixed_value=1.0, model="Fixed value")
        
    # --- Read line regions ---
    
    logging.info("Reading line regions...")
    
    line_regions = ispec.read_line_regions("%s/input/regions/iSPar_Line_Regions_Filtered_Combined.dat" % ispec_dir)
    line_regions = ispec.adjust_linemasks(n_star_spectrum, line_regions, max_margin=0.05, check_derivatives=True)
    
    # --- Create line segments ---
    
    logging.info("Creating line segments...")
    
    segments = ispec.create_segments_around_lines(line_regions, margin=0.10)
    
    initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff = \
        ispec.estimate_initial_ap(n_star_spectrum, precomputed_grid_dir, to_resolution, line_regions)
    
    return initial_teff, initial_logg, initial_MH, initial_alpha

# End of estimate_initial_parameters_with_precomputed_grid

################################################################################

def synthesize_spectrum(initial_R, modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, results_dirname, synth_filename):
    
    # --- Prepare atmosphere model ---
    
    logging.info("Preparing atmosphere model...")
    
    atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, {'teff':initial_teff, 'logg':initial_logg, 'MH':initial_MH, 'alpha':initial_alpha}, code=code)
    
    # --- Spectrum synthesis ---
    
    logging.info("Synthesizing spectrum...")
    
    fixed_abundances = None
    regions = None
    
    synth_spectrum = ispec.create_spectrum_structure(np.arange(wave_base, wave_top, wave_step))
    
    synth_spectrum['flux'] = ispec.generate_spectrum(synth_spectrum['waveobs'], \
                                                     atmosphere_layers, initial_teff, initial_logg, initial_MH, initial_alpha, \
                                                     atomic_linelist, isotopes, solar_abundances, \
                                                     fixed_abundances, microturbulence_vel=initial_vmic, \
                                                     macroturbulence=initial_vmac, vsini=initial_vsini, \
                                                     limb_darkening_coeff=initial_limb_darkening_coeff, \
                                                     R=initial_R, regions=regions, verbose=0, \
                                                     code=code)
    
    # --- Save synthesized spectrum ---
    
    logging.info("Saving synthesized spectrum...")
    
    ispec.write_spectrum(synth_spectrum, "./%s/%s" % (results_dirname, synth_filename))
    
    return synth_spectrum

# End of synthesize_spectrum

################################################################################

def find_synth_continuum_regions(synth_spectrum, initial_R, SNR, results_dirname, synth_continuum_regions_filename):
    
    # Find synth continuum regions
    
    logging.info("Finding continuum regions...")
    
    synth_continuum_model = ispec.fit_continuum(synth_spectrum, fixed_value=1.0, model="Fixed value")
    
    max_std_continuum = 0.001 # Maximum standard deviation (0.1%)
    max_continuum_diff = np.round((1.0/float(SNR)), 5) # Maximum fitted continuum difference
    spectral_profile = 4*wave_step # Equivalent aprox. to spectral profile of TIGRE
    
    synth_continuum_regions = ispec.find_continuum(synth_spectrum, initial_R, \
                                                   max_std_continuum=max_std_continuum, \
                                                   continuum_model=synth_continuum_model, \
                                                   max_continuum_diff=max_continuum_diff, \
                                                   fixed_wave_step=spectral_profile)
    
    # --- Avoid wide telluric regions ---
    
    logging.info("Avoiding wide telluric regions...")
    
    telluric_regions_1 = np.logical_and(synth_continuum_regions['wave_base'] >= 687.35, synth_continuum_regions['wave_top'] <= 740.00)
    telluric_regions_2 = np.logical_and(synth_continuum_regions['wave_base'] >= 759.16, synth_continuum_regions['wave_top'] <= 771.50)
    telluric_regions_3 = np.logical_and(synth_continuum_regions['wave_base'] >= 786.00, synth_continuum_regions['wave_top'] <= 845.00)
    telluric_regions = np.logical_or(telluric_regions_1, telluric_regions_2)
    telluric_regions = np.logical_or(telluric_regions, telluric_regions_3)
    synth_continuum_regions = synth_continuum_regions[~telluric_regions]
    
    # --- Save continuum regions ---
    
    logging.info("Saving continuum regions...")
    
    ispec.write_continuum_regions(synth_continuum_regions, "./%s/%s" % (results_dirname, synth_continuum_regions_filename))
    
    return synth_continuum_regions

# End of find_synth_continuum_regions

################################################################################

def continuum_normalization(m_star_spectrum, initial_R, SNR, synth_spectrum, synth_continuum_regions, \
                            vel_telluric, telluric_linelist, atomic_linelist, results_dirname, normalizations_dirname, residual_results_filename):
    
    # --- Create normalizations directory ---
    
    logging.info("Creating normalizations directory...")
    
    if not os.path.exists("./%s/%s" % (results_dirname, normalizations_dirname)):
        
        os.mkdir("./%s/%s" % (results_dirname, normalizations_dirname))
        logging.info("Directory /%s/%s created" % (results_dirname, normalizations_dirname))
        
    else:
        
        logging.info("Directory /%s/%s already exists!" % (results_dirname, normalizations_dirname))    
        
    # --- Cut synth spectrum from continuum regions ---
    
    logging.info("Cutting synth spectrum from continuum regions...")
    
    # Keep only fluxes inside a list of continuum regions
    
    synth_wfilter = ispec.create_wavelength_filter(synth_spectrum, regions=synth_continuum_regions)
    cutted_synth_spectrum = synth_spectrum[synth_wfilter]
    
    # --- Resample cutted synth spectrum ---
    
    logging.info("Resampling cutted synth spectrum...")
    
    wavelengths = np.arange(wave_base, wave_top, wave_step)
    resampled_synth_spectrum = ispec.resample_spectrum(cutted_synth_spectrum, wavelengths, method="bessel", zero_edges=False)
        
    # --- Define noise level ---
    
    flux_base = np.round(1.0 - (1.0/float(SNR)), 5) # Noise Level
    flux_top = np.round(1.0 + (1.0/float(SNR)), 5)
    
    logging.info("The noise level is: %.6f" % (1.0/float(SNR)))    

    # --- Normalization ---
    
    fit_model = "Splines" # "Polynomy"
    degree = 3 # 2, 3
    nknots = int(np.round((wave_top-wave_base)/2, 0)) # None: 1 spline every 5 nm
    order='median+max' # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
    
    logging.info("Continuum model defined as: %s, with degree: %s, and nknots: %s" % (fit_model, degree, nknots))
    
    # median_wave_range = {0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 0.070, 0.080, 0.090, 0.100}
    
    median_base = 0.010
    median_top = 0.100
    median_step = 0.010
    
    # max_wave_range = {0.15, 0.30, 0.45, 0.60, 0.75, 0.90, 1.05, 1.20, 1.35, 1.50}
    
    max_base = 0.15
    max_top = 1.50
    max_step = 0.15
    
    residuals_columns = ['nfile', 'median_wave_range', 'max_wave_range', 'residuum', 'rms']
    residuals = pd.DataFrame(columns = residuals_columns)

    # Start the loop

    logging.info("Processing...")
    
    for median_wave_range in frange(median_base, median_top, median_step):
    #for median_wave_range in frange(0.100, 0.100, 0.010):
        for max_wave_range in frange(max_base, max_top, max_step):
        #for max_wave_range in frange(0.50, 0.50, 0.15):
            
            logging.info("Median Wavelenght Range: %.3f, Max Wavelenght Range: %.2f..." % (median_wave_range, max_wave_range))            
        
            star_continuum_model = ispec.fit_continuum(m_star_spectrum, \
                                                       ignore=strong_lines, \
                                                       nknots=nknots, degree=degree, \
                                                       median_wave_range=median_wave_range, \
                                                       max_wave_range=max_wave_range, \
                                                       model=fit_model, order=order, \
                                                       automatic_strong_line_detection=True, \
                                                       strong_line_probability=0.70, \
                                                       use_errors_for_fitting=True)        
            
            n_star_spectrum = ispec.normalize_spectrum(m_star_spectrum, star_continuum_model, consider_continuum_errors=False)
            n_spectrum_filename = 'N_S_%s_%.3f_%.2f_%s' % (degree, median_wave_range, max_wave_range, m_spectrum_filename)
            
            # --- Cut normalized spectrum from continuum regions ---
            
            logging.info("Cutting normalized spectrum from continuum regions...")    
            
            # Keep only fluxes inside a list of continuum regions 
            
            normalized_wfilter = ispec.create_wavelength_filter(n_star_spectrum, regions=synth_continuum_regions)
            cutted_normalized_spectrum = n_star_spectrum[normalized_wfilter]
            
            # --- Clean bad continuum fluxes ---
            
            logging.info("Cleaning bad continuum fluxes...")
            
            clean_ffilter = (cutted_normalized_spectrum['flux'] >= flux_base) & (cutted_normalized_spectrum['flux'] <= flux_top)
            cutted_normalized_spectrum = cutted_normalized_spectrum[clean_ffilter]
            
            # --- Resample cutted and normalized spectrum ---
            
            logging.info("Resampling cutted and normalized spectrum...")
            
            resampled_star_spectrum = ispec.resample_spectrum(cutted_normalized_spectrum, wavelengths, method="bessel", zero_edges=False)
            
            # --- Calculate continuum residuals ---
            
            logging.info("Calculating continuum residuals...")
            
            continuum_residuals = ispec.create_spectrum_structure(resampled_star_spectrum['waveobs'])
            continuum_residuals['flux'] = (resampled_star_spectrum['flux'] - resampled_synth_spectrum['flux'])
            continuum_residuals_number = len(continuum_residuals['flux'])
            mean_continuum_residuum = continuum_residuals['flux'].sum()/continuum_residuals_number
            
            logging.info("Mean of continuum residuals = %.10f" % mean_continuum_residuum)
            
            # --- Calculate continuum rms ---
            
            logging.info("Calculating continuum rms...")
            
            continuum_rms = np.sqrt(((continuum_residuals['flux'])**2).sum()/continuum_residuals_number)
            
            logging.info("RMS of continuum fit = %.10f" % continuum_rms)
        
            residuals_result = {'nfile':n_spectrum_filename, 'median_wave_range':median_wave_range, \
                                'max_wave_range':max_wave_range, 'residuum':mean_continuum_residuum, 'rms':continuum_rms}
            
            residuals = residuals.append(residuals_result, ignore_index=True)
            
            if save_all_residuals == True:
                
                # --- Save continuum residuals ---
                
                continuum_residuals_filename = "Residuals_%.3f_%.2f_%s" % (median_wave_range, max_wave_range, m_spectrum_filename)
                
                logging.info("Saving continuum residuals...")
                
                ispec.write_spectrum(continuum_residuals, "./%s/%s/%s" % (results_dirname, normalizations_dirname, continuum_residuals_filename))                                
                                
    # --- Look best normalized spectrum ---
    
    logging.info("Looking best normalized spectrum...")
    
    best_residual_results = residuals.sort_values(by='rms', ascending=True)    
    best_n_spectrum_filename = best_residual_results.iloc[0]['nfile']
    best_median_wave_range = best_residual_results.iloc[0]['median_wave_range']
    best_max_wave_range = best_residual_results.iloc[0]['max_wave_range']
    best_residuum = best_residual_results.iloc[0]['residuum']
    best_continuum_rms = best_residual_results.iloc[0]['rms']
    
    logging.info("The best normalized spectrum file is: %s" % best_n_spectrum_filename)
    logging.info("The RMS of best normalized spectrum file is: %.10f" % best_residual_results.iloc[0]['rms'])
    
    # Best continuum fit
    
    best_star_continuum_model = ispec.fit_continuum(m_star_spectrum, \
                                                    ignore=strong_lines, \
                                                    nknots=nknots, degree=degree, \
                                                    median_wave_range=best_median_wave_range, \
                                                    max_wave_range=best_max_wave_range, \
                                                    model=fit_model, order=order, \
                                                    automatic_strong_line_detection=True, \
                                                    strong_line_probability=0.70, \
                                                    use_errors_for_fitting=True)
    
    # Best continuum normalization
    
    best_n_star_spectrum = ispec.normalize_spectrum(m_star_spectrum, best_star_continuum_model, consider_continuum_errors=False)

    #final_fixed_value = np.round((1.0 + abs(best_residuum) - (1.0/float(SNR))), 5)        
    #final_star_continuum_model = ispec.fit_continuum(best_n_star_spectrum, fixed_value=final_fixed_value, model="Fixed value")
    #best_n_star_spectrum = ispec.normalize_spectrum(best_n_star_spectrum, final_star_continuum_model, consider_continuum_errors=False)

    # --- Second Normalization (only B channel) ---
    
    if channel == "B":
        
        logging.info("Applying second normalization...")
        
        fit_model = "Template"
        nknots = int(np.round((wave_top-wave_base)/2, 0)) # None: 1 spline every 5 nm

        # median_wave_range = {1.0, 2.0, 3.0, 4.0, 5.0}
    
        median_base = 1.0
        median_top = 5.0
        median_step = 1.0
        
        for median_wave_range in frange(median_base, median_top, median_step):
        
            star_continuum_model = ispec.fit_continuum(best_n_star_spectrum, \
                                                       ignore=strong_lines, \
                                                       nknots=nknots, \
                                                       median_wave_range=median_wave_range, \
                                                       model=fit_model, \
                                                       template=synth_spectrum)
        
            best_n_star_spectrum = ispec.normalize_spectrum(best_n_star_spectrum, star_continuum_model, consider_continuum_errors=False)
    
    # --- Radial Velocity verification with linelist mask ---
    
    logging.info("Radial velocity verification with linelist mask...")
    
    ccf_mask = ispec.read_cross_correlation_mask(mask_file)
    
    models, ccf = ispec.cross_correlate_with_mask(best_n_star_spectrum, ccf_mask, \
                                                  lower_velocity_limit=-200, upper_velocity_limit=200, \
                                                  velocity_step=0.5, mask_depth=0.01, \
                                                  fourier=False)
    
    # Number of models represent the number of components
    
    components = len(models)
    
    # First component:
    
    rv = np.round(models[0].mu(), 2) # km/s
    rv_err = np.round(models[0].emu(), 2) # km/s
    
    # --- Radial Velocity correction ---
    
    logging.info("Radial velocity correction: %.2f +/- %.2f" % (rv, rv_err))
    best_n_star_spectrum = ispec.correct_velocity(best_n_star_spectrum, rv)
    
    # --- Save best normalized spectrum ---                
    
    logging.info("Saving best normalized spectrum...")
    ispec.write_spectrum(best_n_star_spectrum, "./%s/%s/%s" % (results_dirname, normalizations_dirname, best_n_spectrum_filename))
    
    # --- Save residual results ---
    
    logging.info("Saving residual results...")
        
    residuals['median_wave_range'] = residuals['median_wave_range'].map(lambda x: '%.3f' % x)
    residuals['max_wave_range'] = residuals['max_wave_range'].map(lambda x: '%.2f' % x)
    residuals['residuum'] = residuals['residuum'].map(lambda x: '%.10f' % x)
    residuals['rms'] = residuals['rms'].map(lambda x: '%.10f' % x)
    
    residuals.to_csv ("./%s/%s" % (results_dirname, residual_results_filename), index=False, header=True, sep='\t')

    return best_n_spectrum_filename

# End of continuum_normalization

################################################################################

def find_line_regions_function(star_spectrum, star_continuum_model, initial_R, SNR, \
                               vel_telluric, telluric_linelist, atomic_linelist, \
                               results_dirname, star_line_regions_filename, star_abundances_filename):
    
    # --- Finding line regions for parameters ---
    
    logging.info("Finding line regions for parameters...")
    
    min_depth = 5.0*(1/float(SNR)) # Maximum fitted continuum difference (Rose criterion)
    #min_depth = 0.05
    max_depth = 1.0
    
    logging.info("minimum depth = %.3f, maximum depth = %.3f" % (min_depth, max_depth))
    
    # Fit the lines and cross-match with atomic linelist
    
    line_regions = ispec.find_linemasks(star_spectrum, star_continuum_model, \
                                        atomic_linelist, \
                                        max_atomic_wave_diff=wave_step, \
                                        telluric_linelist=telluric_linelist, \
                                        vel_telluric=vel_telluric, \
                                        minimum_depth=min_depth, maximum_depth=max_depth, \
                                        check_derivatives=True, \
                                        discard_gaussian=False, \
                                        discard_voigt=True, \
                                        closest_match=False)
    
    # --- Avoid wide telluric regions ---
    
    logging.info("Avoiding wide telluric regions...")
    
    telluric_regions_1 = np.logical_and(line_regions['wave_base'] >= 687.35, line_regions['wave_top'] <= 742.50)
    telluric_regions_2 = np.logical_and(line_regions['wave_base'] >= 759.16, line_regions['wave_top'] <= 772.50)
    telluric_regions_3 = np.logical_and(line_regions['wave_base'] >= 786.00, line_regions['wave_top'] <= 845.00)
    telluric_regions = np.logical_or(telluric_regions_1, telluric_regions_2)
    telluric_regions = np.logical_or(telluric_regions, telluric_regions_3)
    line_regions = line_regions[~telluric_regions]

    # --- Avoid wings of strong lines ---
    
    logging.info("Avoiding wings of strong lines...")
    
    wings_regions_1 = np.logical_and(line_regions['wave_base'] >= 652.50, line_regions['wave_top'] <= 659.50) # H-alpha
    wings_regions_2 = np.logical_and(line_regions['wave_base'] >= 588.00, line_regions['wave_top'] <= 590.50) # Na doublet
    wings_regions_3 = np.logical_and(line_regions['wave_base'] >= 853.00, line_regions['wave_top'] <= 855.50) # Ca doublet
    wings_regions_4 = np.logical_and(line_regions['wave_base'] >= 865.20, line_regions['wave_top'] <= 867.30) # Ca doublet
    wings_regions = np.logical_or(wings_regions_1, wings_regions_2)
    wings_regions = np.logical_or(wings_regions, wings_regions_3)
    wings_regions = np.logical_or(wings_regions, wings_regions_4)
    line_regions = line_regions[~wings_regions]
    
    # --- Discard bad masks ---
    
    logging.info("Discarting bad masks...")                                                      

    flux_peak = star_spectrum['flux'][line_regions['peak']]
    flux_base = star_spectrum['flux'][line_regions['base']]
    flux_top = star_spectrum['flux'][line_regions['top']]
    bad_mask = np.logical_or(line_regions['wave_peak'] <= line_regions['wave_base'], line_regions['wave_peak'] >= line_regions['wave_top'])
    bad_mask = np.logical_or(bad_mask, flux_peak >= flux_base)
    bad_mask = np.logical_or(bad_mask, flux_peak >= flux_top)
    line_regions = line_regions[~bad_mask]
    
    # Exclude lines that have not been successfully cross matched with the atomic data
    # because we cannot calculate the chemical abundance (it will crash the corresponding routines)
    
    rejected_by_atomic_line_not_found = (line_regions['wave_nm'] == 0)
    line_regions = line_regions[~rejected_by_atomic_line_not_found]
    
    # Exclude lines with EW equal to zero
    
    rejected_by_zero_ew = (line_regions['ew'] == 0)
    line_regions = line_regions[~rejected_by_zero_ew]
    
    # Exclude lines that may be affected by tellurics
    
    rejected_by_telluric_line = (line_regions['telluric_wave_peak'] != 0)
    line_regions = line_regions[~rejected_by_telluric_line]

    # Only iron peak elements (Fe, Cr, Ni) and alpha elements (Si, Ca, Ti)
    
    elements_1 = line_regions['element'] == "Fe 1"
    elements_2 = line_regions['element'] == "Fe 2"
    elements_3 = line_regions['element'] == "Cr 1"
    elements_4 = line_regions['element'] == "Cr 2"
    elements_5 = line_regions['element'] == "Ni 1"
    elements_6 = line_regions['element'] == "Ni 2"
    elements_7 = line_regions['element'] == "Si 1"
    elements_8 = line_regions['element'] == "Si 2"
    elements_9 = line_regions['element'] == "Ca 1"
    elements_10 = line_regions['element'] == "Ca 2"
    elements_11 = line_regions['element'] == "Ti 1"
    elements_12 = line_regions['element'] == "Ti 2"
    
    elements = np.logical_or(elements_1, elements_2)
    elements = np.logical_or(elements, elements_3)
    elements = np.logical_or(elements, elements_4)
    elements = np.logical_or(elements, elements_5)
    elements = np.logical_or(elements, elements_6)
    elements = np.logical_or(elements, elements_7)
    elements = np.logical_or(elements, elements_8)
    elements = np.logical_or(elements, elements_9)
    elements = np.logical_or(elements, elements_10)
    elements = np.logical_or(elements, elements_11)
    elements = np.logical_or(elements, elements_12)

    line_regions = line_regions[elements]
    
    # Save line_regions
    
    logging.info("Saving line regions...")

    # line masks + notes
    ispec.write_line_regions(line_regions, "./%s/%s" % (results_dirname, star_line_regions_filename))
    
    # line masks + atomic cross-matched information + fit information
    #ispec.write_line_regions(line_regions, "./%s/%s" % (results_dirname, star_line_regions_filename), extended=True)

    if calculate_individual_abundances == True:
        
        if path.exists("%s/%s" % (results_dirname, star_abundances_filename)) == False:
            
            # --- Calculate individual abundances ---
            
            free_params = ["vrad"]
            
            individual_abundance_columns = ['element', 'wave_peak', 'wave_base', 'wave_top', '[X/H]', 'e[X/H]']
            individual_abundance = pd.DataFrame(columns = individual_abundance_columns)
            
            for i, line in enumerate(line_regions):
                
                # Free individual element abundance (WARNING: it should be coherent with the selected line regions!)
                
                free_abundances = ispec.create_free_abundances_structure([line['element'].split()[0]], chemical_elements, solar_abundances)
                free_abundances['Abund'] += initial_MH # Scale to metallicity
                
                # Line by line
                individual_line_regions = line_regions[i:i+1] # Keep recarray structure
                
                # Segments
                segments = ispec.create_segments_around_lines(individual_line_regions, margin=0.10)
                wfilter = ispec.create_wavelength_filter(star_spectrum, regions=segments) # Only use the segment
                
                if len(star_spectrum[wfilter]) == 0 or np.any(star_spectrum['flux'][wfilter] == 0):
                    continue

                if estimate_vmic_and_vmac == True:

                    obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
                        ispec.model_spectrum(star_spectrum[wfilter], star_continuum_model, \
                                             modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, \
                                             initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, \
                                             initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, segments=segments, \
                                             linemasks=individual_line_regions, \
                                             enhance_abundances=enhance_abundances, \
                                             use_errors = use_errors, \
                                             vmic_from_empirical_relation = True, \
                                             vmac_from_empirical_relation = True, \
                                             max_iterations=max_iterations, \
                                             tmp_dir = None, \
                                             code=code)

                elif estimate_vmic_and_vmac == False:
                    
                    obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
                        ispec.model_spectrum(star_spectrum[wfilter], star_continuum_model, \
                                             modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, \
                                             initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, \
                                             initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, segments=segments, \
                                             linemasks=individual_line_regions, \
                                             enhance_abundances=enhance_abundances, \
                                             use_errors = use_errors, \
                                             vmic_from_empirical_relation = False, \
                                             vmac_from_empirical_relation = False, \
                                             max_iterations=max_iterations, \
                                             tmp_dir = None, \
                                             code=code)
                    
                individual_abundance_result = {'element':line['element'], 'wave_peak':line['wave_peak'], \
                                               'wave_base':line['wave_base'], 'wave_top':line['wave_top'], \
                                               '[X/H]':''.join(map(str, abundances_found['[X/H]'])), \
                                               'e[X/H]':''.join(map(str, abundances_found['e[X/H]']))}
                
                individual_abundance = individual_abundance.append(individual_abundance_result, ignore_index=True)
        
            logging.info("Saving results...")
        
            individual_abundance.to_csv ("./%s/%s" % (results_dirname, star_abundances_filename), index=False, header=True, sep='\t')

        # Read the individual abundances file
        
        individual_abundance = pd.read_csv("%s/%s" % (results_dirname, star_abundances_filename), sep='\t')
        
        # --- Filter lines ---
        
        logging.info("Filtering lines...")
        
        # Determine number of lines by element
        
        lines_by_element = individual_abundance.groupby(['element']).agg({'[X/H]':'size'})
        lines_by_element = lines_by_element.xs('[X/H]', axis=1, drop_level=True)
        lines_by_element = lines_by_element.reset_index('element')
        lines_by_element.rename(columns={'[X/H]':'lines'}, inplace=True)
        
        lines_by_element.to_csv("./%s/%s" % (results_dirname, lines_by_element_filename), index=False, header=True, sep='\t')
        
        filtered_individual_abundance_columns = ['element', 'wave_peak', 'wave_base', 'wave_top', '[X/H]', 'e[X/H]']
        filtered_individual_abundance = pd.DataFrame(columns = filtered_individual_abundance_columns)
        
        filtered_line_regions_columns = ['wave_peak', 'wave_base', 'wave_top', 'note']
        filtered_line_regions = pd.DataFrame(columns = filtered_line_regions_columns)
        
        for i in range(len(individual_abundance)):
            for j in range(len(lines_by_element)):
                
                element_IA = individual_abundance.loc[i, "element"]
                element_LBE = lines_by_element.loc[j, "element"]
                
                if element_IA == element_LBE:
                    
                    XH = individual_abundance.loc[i, "[X/H]"]
                    eXH = individual_abundance.loc[i, "e[X/H]"]
                    
                    if lines_by_element.loc[j, "lines"] > 20:
                        
                        U_XH = initial_MH + 0.25
                        L_XH = initial_MH - 0.25
                        
                    if lines_by_element.loc[j, "lines"] <= 20:
                    
                        U_XH = initial_MH + 0.50
                        L_XH = initial_MH - 0.50
                    
                    if XH <= U_XH and XH >= L_XH:
                        
                        if eXH <= 0.25:
                            
                            filtered_individual_abundance_result = {'element':individual_abundance.loc[i, "element"], \
                                                                    'wave_peak':individual_abundance.loc[i, "wave_peak"], \
                                                                    'wave_base':individual_abundance.loc[i, "wave_base"], \
                                                                    'wave_top':individual_abundance.loc[i, "wave_top"], \
                                                                    '[X/H]':individual_abundance.loc[i, "[X/H]"], \
                                                                    'e[X/H]':individual_abundance.loc[i, "e[X/H]"]}
                            
                            filtered_individual_abundance = filtered_individual_abundance.append(filtered_individual_abundance_result, ignore_index=True)
                            
                            filtered_line_regions_result = {'wave_peak':individual_abundance.loc[i, "wave_peak"], \
                                                            'wave_base':individual_abundance.loc[i, "wave_base"], \
                                                            'wave_top':individual_abundance.loc[i, "wave_top"], \
                                                            'note':individual_abundance.loc[i, "element"]}
                            
                            filtered_line_regions = filtered_line_regions.append(filtered_line_regions_result, ignore_index=True)

        logging.info("Saving results...")
                    
        filtered_individual_abundance.to_csv ("./%s/%s" % (results_dirname, filtered_individual_abundance_filename), index=False, header=True, sep='\t')
        
        filtered_line_regions.to_csv ("./%s/%s" % (results_dirname, filtered_line_regions_filename), index=False, header=True, sep='\t')
        
        line_regions = ispec.read_line_regions("%s/%s" % (results_dirname, filtered_line_regions_filename))

    return line_regions

# End of find_line_regions_function

################################################################################

def fit_lines_determine_ew_and_crossmatch_with_atomic_data(telluric_linelist, modeled_layers_pack, \
                                                           atomic_linelist, chemical_elements, isotopes, solar_abundances):
    
    # --- Read spectrum ---
    
    logging.info("Reading spectrum")

    c_star_spectrum = ispec.read_spectrum("Subtracted_spectrum.dat")
    c_star_continuum_model = ispec.fit_continuum(c_star_spectrum, fixed_value=1.0, model="Fixed value")
    
    # --- Fit lines ---
    
    logging.info("Fitting lines...")
    
    chromosphere_lines = ispec.read_line_regions(ispec_dir + "/input/regions/chromosphere_line_regions.txt")
    chromosphere_lines['wave_peak'] = chromosphere_lines['wave_peak']/10
    chromosphere_lines['wave_base'] = chromosphere_lines['wave_base']/10
    chromosphere_lines['wave_top']  = chromosphere_lines['wave_peak']/10
    chromosphere_lines = ispec.adjust_linemasks(c_star_spectrum, chromosphere_lines, max_margin=0.05, check_derivatives=True)
    
    chromosphere_line_regions = ispec.fit_lines(chromosphere_lines, c_star_spectrum, c_star_continuum_model, \
                                                atomic_linelist = atomic_linelist, \
                                                max_atomic_wave_diff = wave_step, \
                                                smoothed_spectrum = None, \
                                                telluric_linelist = telluric_linelist, \
                                                vel_telluric=vel_telluric, \
                                                check_derivatives=True, \
                                                discard_gaussian=False, \
                                                discard_voigt=True, \
                                                free_mu=True, \
                                                crossmatch_with_mu=False, \
                                                closest_match=False)
    
    # --- Discard bad masks ---
    
    #logging.info("Discarting bad masks...")                                                      
    
    #flux_peak = n_star_spectrum['flux'][chromosphere_line_regions['peak']]
    #flux_base = n_star_spectrum['flux'][chromosphere_line_regions['base']]
    #flux_top = n_star_spectrum['flux'][chromosphere_line_regions['top']]
    #bad_mask = np.logical_or(chromosphere_line_regions['wave_peak'] <= chromosphere_line_regions['wave_base'], chromosphere_line_regions['wave_peak'] >= chromosphere_line_regions['wave_top'])
    #bad_mask = np.logical_or(bad_mask, flux_peak >= flux_base)
    #bad_mask = np.logical_or(bad_mask, flux_peak >= flux_top)
    #chromosphere_line_regions = chromosphere_line_regions[~bad_mask]
    
    # Exclude lines that have not been successfully cross matched with the atomic data
    # because we cannot calculate the chemical abundance (it will crash the corresponding routines)
    
    rejected_by_atomic_line_not_found = (chromosphere_line_regions['wave_nm'] == 0)
    chromosphere_line_regions = chromosphere_line_regions[~rejected_by_atomic_line_not_found]
    
    # Exclude lines with EW equal to zero
    
    rejected_by_zero_ew = (chromosphere_line_regions['ew'] == 0)
    chromosphere_line_regions = chromosphere_line_regions[~rejected_by_zero_ew]
    
    # Exclude lines that may be affected by tellurics
    
    rejected_by_telluric_line = (chromosphere_line_regions['telluric_wave_peak'] != 0)
    chromosphere_line_regions = chromosphere_line_regions[~rejected_by_telluric_line]
    
    ew =  chromosphere_line_regions['ew']
    ew_err =  chromosphere_line_regions['ew_err']
    
    # Save linemasks (line masks + atomic cross-matched information + fit information)
    ispec.write_line_regions(chromosphere_line_regions, "Test_1.dat", extended=True)
    
# End fit_lines_determine_ew_and_crossmatch_with_atomic_data    
    
################################################################################

def calculate_stellar_parameters(n_star_spectrum, n_star_continuum_model, initial_R, line_regions, segments, \
                                 modeled_layers_pack, atomic_linelist, chemical_elements, isotopes, solar_abundances):
    
    # --- Calculate Stellar Parameters ---
    
    logging.info("Calculating stellar parameters for %s..." % star)    

    if calculate_logg == True:

        if estimate_vmic_and_vmac == True:
        
            # --- Calculate Teff & [M/H] ---
            
            logging.info("Calculating Teff & [M/H] for %s..." % star)
            
            free_params = ["teff", "MH"]
            
            obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
                ispec.model_spectrum(n_star_spectrum, n_star_continuum_model, \
                                     modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, \
                                     initial_teff = initial_teff, initial_logg = initial_logg, initial_MH = initial_MH, \
                                     initial_alpha = initial_alpha, initial_vmic = initial_vmic, initial_vmac = initial_vmac, \
                                     initial_vsini = initial_vsini, \
                                     initial_limb_darkening_coeff = initial_limb_darkening_coeff, \
                                     initial_R = initial_R, initial_vrad = initial_vrad, \
                                     free_params = free_params, segments = segments, linemasks = line_regions, \
                                     enhance_abundances = enhance_abundances, \
                                     use_errors = use_errors, \
                                     vmic_from_empirical_relation = True, \
                                     vmac_from_empirical_relation = True, \
                                     max_iterations = max_iterations, \
                                     tmp_dir = None, \
                                     code = code)

        elif estimate_vmic_and_vmac == False:
            
            # --- Calculate Teff & [M/H] ---
            
            logging.info("Calculating Teff & [M/H] for %s..." % star)
            
            free_params = ["teff", "MH"]
            
            obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
                ispec.model_spectrum(n_star_spectrum, n_star_continuum_model, \
                                     modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, \
                                     initial_teff = initial_teff, initial_logg = initial_logg, initial_MH = initial_MH, \
                                     initial_alpha = initial_alpha, initial_vmic = initial_vmic, initial_vmac = initial_vmac, \
                                     initial_vsini = initial_vsini, \
                                     initial_limb_darkening_coeff = initial_limb_darkening_coeff, \
                                     initial_R = initial_R, initial_vrad = initial_vrad, \
                                     free_params = free_params, segments = segments, linemasks = line_regions, \
                                     enhance_abundances = enhance_abundances, \
                                     use_errors = use_errors, \
                                     vmic_from_empirical_relation = False, \
                                     vmac_from_empirical_relation = False, \
                                     max_iterations = max_iterations, \
                                     tmp_dir = None, \
                                     code = code)

        final_teff = params['teff']
        e_teff = errors['teff']
        final_MH = params['MH']
        e_MH = errors['MH']
        final_alpha = params['alpha']
        e_alpha = errors['alpha']

        if estimate_vmic_and_vmac == True:
            
            # --- Calculate log g ---
            
            logging.info("Calculating log g for %s..." % star)
            
            free_params = ["logg"]
            
            obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
                ispec.model_spectrum(n_star_spectrum, n_star_continuum_model, \
                                     modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, \
                                     initial_teff = final_teff, initial_logg = initial_logg, initial_MH = final_MH, \
                                     initial_alpha = initial_alpha, initial_vmic = initial_vmic, initial_vmac = initial_vmac, \
                                     initial_vsini = initial_vsini, \
                                     initial_limb_darkening_coeff = initial_limb_darkening_coeff, \
                                     initial_R = initial_R, initial_vrad = initial_vrad, \
                                     free_params = free_params, segments = segments, linemasks = line_regions, \
                                     enhance_abundances = enhance_abundances, \
                                     use_errors = use_errors, \
                                     vmic_from_empirical_relation = True, \
                                     vmac_from_empirical_relation = True, \
                                     max_iterations = max_iterations, \
                                     tmp_dir = None, \
                                     code = code)

        elif estimate_vmic_and_vmac == False:
            
            # --- Calculate log g ---
            
            logging.info("Calculating log g for %s..." % star)
            
            free_params = ["logg"]
            
            obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
                ispec.model_spectrum(n_star_spectrum, n_star_continuum_model, \
                                     modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, \
                                     initial_teff = final_teff, initial_logg = initial_logg, initial_MH = final_MH, \
                                     initial_alpha = initial_alpha, initial_vmic = initial_vmic, initial_vmac = initial_vmac, \
                                     initial_vsini = initial_vsini, \
                                     initial_limb_darkening_coeff = initial_limb_darkening_coeff, \
                                     initial_R = initial_R, initial_vrad = initial_vrad, \
                                     free_params = free_params, segments = segments, linemasks = line_regions, \
                                     enhance_abundances = enhance_abundances, \
                                     use_errors = use_errors, \
                                     vmic_from_empirical_relation = False, \
                                     vmac_from_empirical_relation = False, \
                                     max_iterations = max_iterations, \
                                     tmp_dir = None, \
                                     code = code)

        final_logg = params['logg']
        e_logg = errors['logg']        
        
    elif calculate_logg == False:
        
        if estimate_vmic_and_vmac == True:

            # --- Calculate Teff & [M/H] ---
            
            logging.info("Calculating Teff & [M/H] for %s..." % star)
            
            free_params = ["teff", "MH"]
            
            obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
                ispec.model_spectrum(n_star_spectrum, n_star_continuum_model, \
                                     modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, \
                                     initial_teff = initial_teff, initial_logg = initial_logg, initial_MH = initial_MH, \
                                     initial_alpha = initial_alpha, initial_vmic = initial_vmic, initial_vmac = initial_vmac, \
                                     initial_vsini = initial_vsini, \
                                     initial_limb_darkening_coeff = initial_limb_darkening_coeff, \
                                     initial_R = initial_R, initial_vrad = initial_vrad, \
                                     free_params = free_params, segments = segments, linemasks = line_regions, \
                                     enhance_abundances = enhance_abundances, \
                                     use_errors = use_errors, \
                                     vmic_from_empirical_relation = True, \
                                     vmac_from_empirical_relation = True, \
                                     max_iterations = max_iterations, \
                                     tmp_dir = None, \
                                     code = code)

        elif estimate_vmic_and_vmac == False:
            
            # --- Calculate Teff & [M/H] ---
            
            logging.info("Calculating Teff & [M/H] for %s..." % star)
            
            free_params = ["teff", "MH"]
            
            obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
                ispec.model_spectrum(n_star_spectrum, n_star_continuum_model, \
                                     modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, \
                                     initial_teff = initial_teff, initial_logg = initial_logg, initial_MH = initial_MH, \
                                     initial_alpha = initial_alpha, initial_vmic = initial_vmic, initial_vmac = initial_vmac, \
                                     initial_vsini = initial_vsini, \
                                     initial_limb_darkening_coeff = initial_limb_darkening_coeff, \
                                     initial_R = initial_R, initial_vrad = initial_vrad, \
                                     free_params = free_params, segments = segments, linemasks = line_regions, \
                                     enhance_abundances = enhance_abundances, \
                                     use_errors = use_errors, \
                                     vmic_from_empirical_relation = False, \
                                     vmac_from_empirical_relation = False, \
                                     max_iterations = max_iterations, \
                                     tmp_dir = None, \
                                     code = code)
        
        final_teff = params['teff']
        e_teff = errors['teff']
        final_logg = params['logg']
        e_logg = errors['logg']
        final_MH = params['MH']
        e_MH = errors['MH']
        final_alpha = params['alpha']
        e_alpha = errors['alpha']
    
    final_vmic = params['vmic']
    e_vmic = errors['vmic']
    final_vmac = params['vmac']
    e_vmac = errors['vmac']
        
    # --- Calculate vsini ---
    
    logging.info("Calculating vsini for %s..." % star)
    
    free_params = ["vsini"]
    
    obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
        ispec.model_spectrum(n_star_spectrum, n_star_continuum_model, \
                             modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, \
                             initial_teff = final_teff, initial_logg = final_logg, initial_MH = final_MH, \
                             initial_alpha = final_alpha, initial_vmic = final_vmic, initial_vmac = final_vmac, \
                             initial_vsini = initial_vsini, \
                             initial_limb_darkening_coeff = initial_limb_darkening_coeff, \
                             initial_R = initial_R, initial_vrad = initial_vrad, \
                             free_params = free_params, segments = segments, linemasks = line_regions, \
                             enhance_abundances = enhance_abundances, \
                             use_errors = use_errors, \
                             vmic_from_empirical_relation = False, \
                             vmac_from_empirical_relation = False, \
                             max_iterations = max_iterations, \
                             tmp_dir = None, \
                             code = code)

    final_vsini = params['vsini']
    e_vsini = errors['vsini'] 

    # --- Calculate Fe abundance ---
    
    logging.info("Calculating Fe abundance for %s..." % star)    
    
    element_name = "Fe"
    new_free_abundances = ispec.create_free_abundances_structure([element_name], chemical_elements, solar_abundances)
    new_free_abundances['Abund'] += initial_MH # Scale to metallicity

    free_params = []
    
    obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
        ispec.model_spectrum(n_star_spectrum, n_star_continuum_model, \
                             modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, new_free_abundances, linelist_free_loggf, \
                             initial_teff = final_teff, initial_logg = final_logg, initial_MH = final_MH, \
                             initial_alpha = final_alpha, initial_vmic = final_vmic, initial_vmac = final_vmac, \
                             initial_vsini = final_vsini, \
                             initial_limb_darkening_coeff = initial_limb_darkening_coeff, \
                             initial_R = initial_R, initial_vrad = initial_vrad, \
                             free_params = free_params, segments = segments, linemasks = line_regions, \
                             enhance_abundances = enhance_abundances, \
                             use_errors = use_errors, \
                             vmic_from_empirical_relation = False, \
                             vmac_from_empirical_relation = False, \
                             max_iterations = max_iterations, \
                             tmp_dir = None, \
                             code = code)

    final_FeH = abundances_found['[X/H]']
    e_FeH = abundances_found['e[X/H]']
    
    return final_teff, e_teff, final_logg, e_logg, final_MH, e_MH, final_alpha, e_alpha, final_vmic, e_vmic, final_vmac, e_vmac, final_vsini, e_vsini, final_FeH, e_FeH, status['rms']
    
# End of calculate_stellar_parameters

################################################################################

if __name__ == '__main__':
        
    star = os.path.basename(os.getcwd())

    # Define filenames

    coadded_spectrum_filename = "A_%s_%s.dat" % (star, channel)
    synth_filename = "iSPar_Synth_Spectrum_%s.dat" % star
    synth_continuum_regions_filename = "iSPar_Synth_Continuum_Regions_%s.dat" % star
    star_line_regions_filename = "iSPar_Line_Regions_%s.dat" % star
    filtered_line_regions_filename = "iSPar_Line_Regions_Filtered_%s.dat" % star
    star_abundances_filename = "iSPar_Individual_Abundances_%s.dat" % star
    lines_by_element_filename = "iSPar_Lines_By_Element_%s.dat" % star
    filtered_individual_abundance_filename = "iSPar_Individual_Abundances_Filtered_%s.dat" % star
    params_results_filename = "iSPar_Results_%s.dat" % star
    residual_results_filename = "iSPar_Residuals_Results_%s.dat" % star
    normalizations_dirname = "iSPar Normalizations %s" % star
    results_dirname = "%s iSPar Results" % star
    
    # --- Create resulting directory ---
    
    logging.info("Creating resulting directory...")
    
    if not os.path.exists(results_dirname):
        
        os.mkdir(results_dirname)
        logging.info("Resulting directory for %s created" % star)
        
    else:
        
        logging.info("Resulting directory for %s already exists!" % star)

    # -----------------------------------------------------------------------------#
        
    # --- Load atmospheres model ---
    
    logging.info("Loading atmospheres model...")
    
    modeled_layers_pack = ispec.load_modeled_layers_pack(atm_model)        
        
    # --- Read atomic linelist ---
    
    logging.info("Reading atomic linelist...")
    
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base, wave_top)
    
    chemical_elements_file = '%s/input/abundances/chemical_elements_symbols.dat' % ispec_dir
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
    
    # --- Select lines that have some minimal contribution in the Sun ---
    
    logging.info("Selecting lines that have some minimal contribution in the Sun (>= 0.01)...")
    
    atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01]
    
    # --- Read isotopes list ---
    
    logging.info("Reading isotopes list...")
    
    isotope_file = '%s/input/isotopes/SPECTRUM.lst' % ispec_dir
    isotopes = ispec.read_isotope_data(isotope_file)
    
    # --- Load solar abundances ---
    
    logging.info("Loading solar abundances...")
    
    solar_abundances = ispec.read_solar_abundances(solar_abundances_file)
    
    # Read atomic data
    
    #mask_file = '%s/input/linelists/CCF/Narval.Sun.370_1048nm/mask.lst' % ispec_dir
    #mask_file = '%s/input/linelists/CCF/Atlas.Sun.372_926nm/mask.lst' % ispec_dir
    #mask_file = '%s/input/linelists/CCF/Synthetic.Sun.350_1100nm/mask.lst' % ispec_dir
    mask_file = '%s/input/linelists/CCF/VALD.Sun.300_1100nm/mask.lst' % ispec_dir

    # --- Load telluric linelist ---
    
    logging.info("Loading telluric linelist...")
                
    telluric_linelist_file = '%s/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst' % ispec_dir
    telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.0)

    # Telluric and strong line regions
        
    strong_lines = ispec.read_line_regions(ispec_dir + "/input/regions/strong_lines/telluric_and_absorption_lines.txt")

    # Create file for stellar parameters results
    
    params_results_columns = ['teff', 'e_teff', 'logg', 'e_logg', 'MH', 'e_MH', 'alpha', 'e_alpha', 'vmic', 'e_vmic', 'vmac', 'e_vmac', 'vsini', 'e_vsini', 'FeH', 'e_FeH', 'rms']
    params_results = pd.DataFrame(columns = params_results_columns)

    # -----------------------------------------------------------------------------#
    
    # --- Start the loop ---

    for i in xrange(0, n_spectra, 1):
        
        # --- Look the spectrum by SNR ---
        
        logging.info("Looking the spectrum by SNR...")
        
        spectra_info_filename = '%s_INFO.dat' % star
        
        if path.exists(spectra_info_filename) == False:
        
            logging.info("Please execute DataRedu.sh first!")
            sys.exit()
            
        else:
            
            spectra_info = pd.read_csv(spectra_info_filename, sep='\t')

            channel_spectra_info = spectra_info[spectra_info['CHANNEL'] == channel]
            
            snr_spectra_info = channel_spectra_info.sort_values(by='SNR', ascending=False)
            
            spectrum_filename = '%s.dat' % snr_spectra_info.iloc[i]['SPECTRUM']

            if use_spectrum_resolution == True:
                
                resolution = int(snr_spectra_info.iloc[i]['RESOLUTION'])

            else:

                resolution = defined_resolution
                
            SNR = int(snr_spectra_info.iloc[i]['SNR'])
            
            logging.info("The spectrum file is: %s" % spectrum_filename)
            logging.info("The resolution is: %s" % resolution)
            logging.info("The SNR from TIGRE is: %s" % SNR)

            # -----------------------------------------------------------------------------#
            
            m_spectrum_filename = 'M_%s' % spectrum_filename
        
            if path.exists("%s/%s" % (results_dirname, m_spectrum_filename)) == False:
                
                # --- Reduce the spectrum (cut and RV correction) ---
                
                wave_step, m_star_spectrum = spectrum_reduction(spectrum_filename, results_dirname, resolution, m_spectrum_filename)
                
            else:
                
                # --- Reading spectra ---
                
                logging.info("Reading modified (RV + cut) spectrum")
                
                m_star_spectrum = ispec.read_spectrum("%s/%s" % (results_dirname, m_spectrum_filename))

                wave_step = np.median(np.abs(m_star_spectrum['waveobs'][1:] - m_star_spectrum['waveobs'][:-1]))

                logging.info("The wave step is: %.5f" % wave_step)

            if resolution_degradation == True:
                
                resolution = final_R
                
                logging.info("The final resolution is: %.0f" % resolution)
                
            # -----------------------------------------------------------------------------#

            if coadd_spectra == True:
                
                if n_spectra <= 6 and n_spectra > 1:                    
                    
                    if i == 0:                    
                        
                        spectrum_1 = m_star_spectrum
                        resolution_1 = resolution
                        SNR_1 = SNR
                        wave_step_1 = wave_step
                        
                    elif i == 1:
                        
                        spectrum_2 = m_star_spectrum
                        resolution_2 = resolution
                        SNR_2 = SNR
                        wave_step_2 = wave_step
                        
                    elif i == 2:
                        
                        spectrum_3 = m_star_spectrum
                        resolution_3 = resolution
                        SNR_3 = SNR
                        wave_step_3 = wave_step
                        
                    elif i == 3:
                        
                        spectrum_4 = m_star_spectrum
                        resolution_4 = resolution
                        SNR_4 = SNR
                        wave_step_4 = wave_step
                        
                    elif i == 4:
                        
                        spectrum_5 = m_star_spectrum
                        resolution_5 = resolution
                        SNR_5 = SNR
                        wave_step_5 = wave_step
                        
                    elif i == 5:
                        
                        spectrum_6 = m_star_spectrum
                        resolution_6 = resolution
                        SNR_6 = SNR
                        wave_step_6 = wave_step
                        
                else:
                        
                    logging.info("A minimum of 2 and a maximum of 6 spectra can be coadded!")
                    sys.exit()

            else:

                initial_R = resolution

    # -----------------------------------------------------------------------------#

    if coadd_spectra == True:
        
        # --- Coadd Spectra ---
        
        logging.info("Coadding spectra...")
        
        if n_spectra == 2:

            wave_step = np.round((wave_step_1 + wave_step_2)/2, 5)
            wavelengths = np.arange(wave_base, wave_top, wave_step)
            
            # --- Resampling and combining ---
            
            logging.info("Resampling and combining...")
            
            # Resample
            
            resampled_spectrum_1 = ispec.resample_spectrum(spectrum_1, wavelengths, method="bessel", zero_edges=False)
            resampled_spectrum_2 = ispec.resample_spectrum(spectrum_2, wavelengths, method="bessel", zero_edges=False)
            
            # Coadd
            
            coadded_spectrum = ispec.create_spectrum_structure(resampled_spectrum_1['waveobs'])
            coadded_spectrum['flux'] = resampled_spectrum_1['flux'] + resampled_spectrum_2['flux']
            coadded_spectrum['err'] = np.sqrt(np.power(resampled_spectrum_1['err'],2) + \
                                              np.power(resampled_spectrum_2['err'],2))

            initial_R = (resolution_1 + resolution_2)/n_spectra

            N_1 = 1/float(SNR_1)
            N_2 = 1/float(SNR_2)

            SNR = 2/math.sqrt(math.pow(N_1,2) + math.pow(N_2,2))
            
        elif n_spectra == 3:

            wave_step = np.round((wave_step_1 + wave_step_2 + wave_step_3)/3, 5)
            wavelengths = np.arange(wave_base, wave_top, wave_step)
            
            # --- Resampling and combining ---
            
            logging.info("Resampling and combining...")        
        
            # Resample
            
            resampled_spectrum_1 = ispec.resample_spectrum(spectrum_1, wavelengths, method="bessel", zero_edges=False)
            resampled_spectrum_2 = ispec.resample_spectrum(spectrum_2, wavelengths, method="bessel", zero_edges=False)
            resampled_spectrum_3 = ispec.resample_spectrum(spectrum_3, wavelengths, method="bessel", zero_edges=False)
        
            # Coadd
            
            coadded_spectrum = ispec.create_spectrum_structure(resampled_spectrum_1['waveobs'])
            coadded_spectrum['flux'] = resampled_spectrum_1['flux'] + resampled_spectrum_2['flux'] + resampled_spectrum_3['flux']
            coadded_spectrum['err'] = np.sqrt(np.power(resampled_spectrum_1['err'],2) + \
                                              np.power(resampled_spectrum_2['err'],2) + \
                                              np.power(resampled_spectrum_3['err'],2))

            initial_R = (resolution_1 + resolution_2 + resolution_3)/n_spectra

            N_1 = 1/float(SNR_1)
            N_2 = 1/float(SNR_2)
            N_3 = 1/float(SNR_3)

            SNR = 3/math.sqrt(math.pow(N_1,2) + math.pow(N_2,2) + math.pow(N_3,2))
        
        elif n_spectra == 4:

            wave_step = np.round((wave_step_1 + wave_step_2 + wave_step_3 + wave_step_4)/4, 5)
            wavelengths = np.arange(wave_base, wave_top, wave_step)
            
            # --- Resampling and combining ---
            
            logging.info("Resampling and combining...")        
            
            # Resample
            
            resampled_spectrum_1 = ispec.resample_spectrum(spectrum_1, wavelengths, method="bessel", zero_edges=False)
            resampled_spectrum_2 = ispec.resample_spectrum(spectrum_2, wavelengths, method="bessel", zero_edges=False)
            resampled_spectrum_3 = ispec.resample_spectrum(spectrum_3, wavelengths, method="bessel", zero_edges=False)
            resampled_spectrum_4 = ispec.resample_spectrum(spectrum_4, wavelengths, method="bessel", zero_edges=False)
            
            # Coadd
            
            coadded_spectrum = ispec.create_spectrum_structure(resampled_spectrum_1['waveobs'])
            coadded_spectrum['flux'] = resampled_spectrum_1['flux'] + resampled_spectrum_2['flux'] + \
                resampled_spectrum_3['flux'] + resampled_spectrum_4['flux']
            coadded_spectrum['err'] = np.sqrt(np.power(resampled_spectrum_1['err'],2) + \
                                              np.power(resampled_spectrum_2['err'],2) + \
                                              np.power(resampled_spectrum_3['err'],2) + \
                                              np.power(resampled_spectrum_4['err'],2))

            initial_R = (resolution_1 + resolution_2 + resolution_3 + resolution_4)/n_spectra

            N_1 = 1/float(SNR_1)
            N_2 = 1/float(SNR_2)
            N_3 = 1/float(SNR_3)
            N_4 = 1/float(SNR_4)

            SNR = 4/math.sqrt(math.pow(N_1,2) + math.pow(N_2,2) + math.pow(N_3,2) + math.pow(N_4,2))
            
        elif n_spectra == 5:

            wave_step = np.round((wave_step_1 + wave_step_2 + wave_step_3 + wave_step_4 + wave_step_5)/5, 5)
            wavelengths = np.arange(wave_base, wave_top, wave_step)
            
            # --- Resampling and combining ---
            
            logging.info("Resampling and combining...")        
            
            # Resample
            
            resampled_spectrum_1 = ispec.resample_spectrum(spectrum_1, wavelengths, method="bessel", zero_edges=False)
            resampled_spectrum_2 = ispec.resample_spectrum(spectrum_2, wavelengths, method="bessel", zero_edges=False)
            resampled_spectrum_3 = ispec.resample_spectrum(spectrum_3, wavelengths, method="bessel", zero_edges=False)
            resampled_spectrum_4 = ispec.resample_spectrum(spectrum_4, wavelengths, method="bessel", zero_edges=False)
            resampled_spectrum_5 = ispec.resample_spectrum(spectrum_5, wavelengths, method="bessel", zero_edges=False)
            
            # Coadd
            
            coadded_spectrum = ispec.create_spectrum_structure(resampled_spectrum_1['waveobs'])
            coadded_spectrum['flux'] = resampled_spectrum_1['flux'] + resampled_spectrum_2['flux'] + \
                resampled_spectrum_3['flux'] + resampled_spectrum_4['flux'] + resampled_spectrum_5['flux']
            coadded_spectrum['err'] = np.sqrt(np.power(resampled_spectrum_1['err'],2) + \
                                              np.power(resampled_spectrum_2['err'],2) + \
                                              np.power(resampled_spectrum_3['err'],2) + \
                                              np.power(resampled_spectrum_4['err'],2) + \
                                              np.power(resampled_spectrum_5['err'],2))

            initial_R = (resolution_1 + resolution_2 + resolution_3 + resolution_4 + resolution_5)/n_spectra

            N_1 = 1/float(SNR_1)
            N_2 = 1/float(SNR_2)
            N_3 = 1/float(SNR_3)
            N_4 = 1/float(SNR_4)
            N_5 = 1/float(SNR_5)
            
            SNR = 5/math.sqrt(math.pow(N_1,2) + math.pow(N_2,2) + math.pow(N_3,2) + math.pow(N_4,2) + math.pow(N_5,2))

        elif n_spectra == 6:

            wave_step = np.round((wave_step_1 + wave_step_2 + wave_step_3 + wave_step_4 + wave_step_5 + wave_step_6)/6, 5)
            wavelengths = np.arange(wave_base, wave_top, wave_step)
            
            # --- Resampling and combining ---
            
            logging.info("Resampling and combining...")        
            
            # Resample
            
            resampled_spectrum_1 = ispec.resample_spectrum(spectrum_1, wavelengths, method="bessel", zero_edges=False)
            resampled_spectrum_2 = ispec.resample_spectrum(spectrum_2, wavelengths, method="bessel", zero_edges=False)
            resampled_spectrum_3 = ispec.resample_spectrum(spectrum_3, wavelengths, method="bessel", zero_edges=False)
            resampled_spectrum_4 = ispec.resample_spectrum(spectrum_4, wavelengths, method="bessel", zero_edges=False)
            resampled_spectrum_5 = ispec.resample_spectrum(spectrum_5, wavelengths, method="bessel", zero_edges=False)
            resampled_spectrum_6 = ispec.resample_spectrum(spectrum_6, wavelengths, method="bessel", zero_edges=False)
            
            # Coadd
            
            coadded_spectrum = ispec.create_spectrum_structure(resampled_spectrum_1['waveobs'])
            coadded_spectrum['flux'] = resampled_spectrum_1['flux'] + resampled_spectrum_2['flux'] + \
                resampled_spectrum_3['flux'] + resampled_spectrum_4['flux'] + resampled_spectrum_5['flux'] + resampled_spectrum_6['flux']
            coadded_spectrum['err'] = np.sqrt(np.power(resampled_spectrum_1['err'],2) + \
                                              np.power(resampled_spectrum_2['err'],2) + \
                                              np.power(resampled_spectrum_3['err'],2) + \
                                              np.power(resampled_spectrum_4['err'],2) + \
                                              np.power(resampled_spectrum_5['err'],2) + \
                                              np.power(resampled_spectrum_6['err'],2))

            initial_R = (resolution_1 + resolution_2 + resolution_3 + resolution_4 + resolution_5 + resolution_6)/n_spectra

            N_1 = 1/float(SNR_1)
            N_2 = 1/float(SNR_2)
            N_3 = 1/float(SNR_3)
            N_4 = 1/float(SNR_4)
            N_5 = 1/float(SNR_5)
            N_6 = 1/float(SNR_6)
            
            SNR = 6/math.sqrt(math.pow(N_1,2) + math.pow(N_2,2) + math.pow(N_3,2) + math.pow(N_4,2) + math.pow(N_5,2) + math.pow(N_6,2))
            
        # --- Save coadded spectra ---
            
        logging.info("Saving coadded spectrum...")
            
        ispec.write_spectrum(coadded_spectrum, "%s/%s" % (results_dirname, coadded_spectrum_filename))
        
        logging.info("The coadded spectrum file is: %s" % coadded_spectrum_filename)
        logging.info("The resolution of coadded spectrum is: %s" % initial_R)
        logging.info("The SNR of coadded spectrum is: %s" % SNR)
        logging.info("The final wave step is: %.5f" % wave_step)

        m_spectrum_filename = coadded_spectrum_filename
        m_star_spectrum = coadded_spectrum

    # -----------------------------------------------------------------------------#
            
    if precompute_synth_grid == True:
        
        # --- Precompute synthetic grid ---
        
        logging.info("Precomputing synthetic grid...")
        
        precompute_synthetic_grid(precomputed_grid_dirname, modeled_layers_pack, atomic_linelist, isotopes, solar_abundances)
        
    # -----------------------------------------------------------------------------#

    # --- Telluric velocity shift determination from spectrum ---
    
    logging.info("Telluric velocity shift determination...")
    
    models, ccf = ispec.cross_correlate_with_mask(m_star_spectrum, telluric_linelist, \
                                                  lower_velocity_limit=-100, upper_velocity_limit=100, \
                                                  velocity_step=0.5, mask_depth=0.01, \
                                                  fourier=False, \
                                                  only_one_peak=True)
    
    vel_telluric = np.round(models[0].mu(), 2) # km/s
    vel_telluric_err = np.round(models[0].emu(), 2) # km/s
    
    logging.info("Telluric velocity shift: %.2f +/- %.2f" % (vel_telluric, vel_telluric_err)) 
    
    # -----------------------------------------------------------------------------#
    
    if estimate_initial_parameters == True:
        
        # --- Estimate initial parameters ---
        
        logging.info("Estimating initial parameters for %s..." % star)
        
        estimated_teff, estimated_logg, estimated_MH, estimated_alpha = estimate_initial_parameters_with_precomputed_grid(precomputed_grid_dirname, m_star_spectrum, initial_R, \
                                                                                                                          vel_telluric, telluric_linelist, normalization=True)
        
        initial_teff = estimated_teff
        initial_logg = estimated_logg
        initial_MH = estimated_MH
        initial_alpha = estimated_alpha
        initial_vsini = 1.6 # Rotation velocity defined by default
        
        # Empirical relation for vmic & vmac
        
        initial_vmic = ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
        initial_vmac = ispec.estimate_vmac(initial_teff, initial_logg, initial_MH)
        
        logging.info("Initial parameter estimation:")
        logging.info("Teff = %.2f, log g = %.2f, [M/H] = %.2f, alpha = %.2f" % (initial_teff, initial_logg, initial_MH, initial_alpha))
        logging.info("vmic = %.2f, vmac = %.2f, vsini = %.2f" % (initial_vmic, initial_vmac, initial_vsini))
        
    elif estimate_initial_parameters == False:
        
        initial_teff = defined_teff
        initial_logg = defined_logg
        initial_MH = defined_MH
        initial_alpha = defined_alpha
        initial_vsini = defined_vsini
        
        if estimate_vmic_and_vmac == False:
            
            initial_vmic = defined_vmic
            initial_vmac = defined_vmac
            
        elif estimate_vmic_and_vmac == True:
            
            # Empirical relation for vmic & vmac
            
            initial_vmic = ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
            initial_vmac = ispec.estimate_vmac(initial_teff, initial_logg, initial_MH)        
            
        # Validate parameters
        
        if not ispec.valid_atmosphere_target(modeled_layers_pack, {'teff':initial_teff, 'logg':initial_logg, 'MH':initial_MH, 'alpha':initial_alpha}):
            logging.info("ERROR: The specified effective temperature (Teff), gravity (logg) and metallicity (M/H) fall out of the atmospheric models!")
            sys.exit()

        logging.info("Initial parameters defined by user:")
        logging.info("Teff = %.2f, log g = %.2f, [M/H] = %.2f, alpha = %.2f" % (initial_teff, initial_logg, initial_MH, initial_alpha))
        logging.info("vmic = %.2f, vmac = %.2f, vsini = %.2f" % (initial_vmic, initial_vmac, initial_vsini))

    # Other stellar parameters are defined by default
    
    initial_limb_darkening_coeff = 0.6
    initial_vrad = 0
            
    # -----------------------------------------------------------------------------#
    
    if path.exists("%s/%s" % (results_dirname, synth_continuum_regions_filename)) == False:
        
        if path.exists("%s/%s" % (results_dirname, synth_filename)) == False:
            
            # --- Synthesize spectrum ---
            
            logging.info("Synthesizing spectrum for %s..." % star)
            
            synth_spectrum = synthesize_spectrum(initial_R, modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, results_dirname, synth_filename)
            
        else:
            
            # --- Read synthesized spectrum ---
            
            logging.info("Reading synthesized spectrum for %s..." % star)
            
            synth_spectrum = ispec.read_spectrum("%s/%s" % (results_dirname, synth_filename))
            
        # --- Find continuum regions ---
            
        logging.info("Finding continuum regions for %s..." % star)
            
        synth_continuum_regions = find_synth_continuum_regions(synth_spectrum, initial_R, SNR, results_dirname, synth_continuum_regions_filename)
            
    else:
        
        # --- Read synthesized spectrum ---
        
        logging.info("Reading synthesized spectrum for %s..." % star)
        
        synth_spectrum = ispec.read_spectrum("%s/%s" % (results_dirname, synth_filename))
        
        # --- Read continuum regions ---
        
        logging.info("Reading continuum regions for %s..." % star)
        
        synth_continuum_regions = ispec.read_continuum_regions("%s/%s" % (results_dirname, synth_continuum_regions_filename))
        
    # -----------------------------------------------------------------------------#
    
    if path.exists("%s/%s" % (results_dirname, residual_results_filename)) == False:
        
        # --- Calculate continuum normalizations ---
        
        logging.info("Calculting continuum normalizations for %s..." % star)
        
        n_spectrum_filename = "%s/%s/%s" % (results_dirname, normalizations_dirname, \
                                            continuum_normalization(m_star_spectrum, initial_R, SNR, synth_spectrum, \
                                                                    synth_continuum_regions, vel_telluric, telluric_linelist, atomic_linelist, \
                                                                    results_dirname, normalizations_dirname, residual_results_filename))
        
    else:
        
        # Read the residuals results file
        
        residuals = pd.read_csv("%s/%s" % (results_dirname, residual_results_filename), sep='\t')
        
        # --- Look best normalized spectrum ---
        
        logging.info("Looking best normalized spectrum...")
        
        best_residual_results = residuals.sort_values(by='rms', ascending=True)
        
        n_spectrum_filename = "%s/%s/%s" % (results_dirname, normalizations_dirname, best_residual_results.iloc[0]['nfile'])
        
        logging.info("The best normalized spectrum file is: %s" % n_spectrum_filename)
        logging.info("The RMS of best normalized spectrum file is: %.10f" % best_residual_results.iloc[0]['rms'])

    # -----------------------------------------------------------------------------#    
    
    # --- Read normalized spectra ---
    
    logging.info("Reading best normalized spectrum")
    
    n_star_spectrum = ispec.read_spectrum(n_spectrum_filename)
    
    # Use a fixed value because the spectrum is already normalized
    
    n_star_continuum_model = ispec.fit_continuum(n_star_spectrum, fixed_value=1.0, model="Fixed value")

    # -----------------------------------------------------------------------------#
    
    # --- Telluric velocity shift verification from spectrum ---
    
    logging.info("Telluric velocity shift verification...")
    
    models, ccf = ispec.cross_correlate_with_mask(n_star_spectrum, telluric_linelist, \
                                                  lower_velocity_limit=-100, upper_velocity_limit=100, \
                                                  velocity_step=0.5, mask_depth=0.01, \
                                                  fourier=False, \
                                                  only_one_peak=True)
    
    vel_telluric = np.round(models[0].mu(), 2) # km/s
    vel_telluric_err = np.round(models[0].emu(), 2) # km/s
    
    logging.info("Telluric velocity shift: %.2f +/- %.2f" % (vel_telluric, vel_telluric_err))

    # -----------------------------------------------------------------------------#    

    if estimate_initial_parameters == True:
        
        # --- Re-estimate initial parameters with the final normalization ---
        
        logging.info("Re-estimating initial parameters for %s..." % star)
        
        estimated_teff, estimated_logg, estimated_MH, estimated_alpha = estimate_initial_parameters_with_precomputed_grid(precomputed_grid_dirname, n_star_spectrum, initial_R, \
                                                                                                                          vel_telluric, telluric_linelist, normalization=False)
                
        initial_teff = estimated_teff
        initial_logg = estimated_logg
        initial_MH = estimated_MH
        initial_alpha = estimated_alpha
        
        # Empirical relation for vmic & vmac
        
        initial_vmic = ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
        initial_vmac = ispec.estimate_vmac(initial_teff, initial_logg, initial_MH)
        
        logging.info("Initial parameter re-estimation:")
        logging.info("Teff = %.2f, log g = %.2f, [M/H] = %.2f, alpha = %.2f" % (initial_teff, initial_logg, initial_MH, initial_alpha))
        logging.info("vmic = %.2f, vmac = %.2f, vsini = %.2f" % (initial_vmic, initial_vmac, initial_vsini))

    # -----------------------------------------------------------------------------#
    
    if find_line_regions == True:
        
        if path.exists("%s/%s" % (results_dirname, filtered_line_regions_filename)) == False:
            
            # --- Find line regions --- #
            
            logging.info("Finding line regions...")
            
            line_regions = find_line_regions_function(n_star_spectrum, n_star_continuum_model, initial_R, SNR, \
                                                      vel_telluric, telluric_linelist, atomic_linelist, \
                                                      results_dirname, star_line_regions_filename, star_abundances_filename)
            
        else:
        
            # --- Read line regions ---
            
            logging.info("Reading line regions...")
            
            line_regions = ispec.read_line_regions("%s/%s" % (results_dirname, filtered_line_regions_filename))
            
        # --- Create line segments ---
        
        logging.info("Creating line segments...")
        
        segments = ispec.create_segments_around_lines(line_regions, margin=0.10)

    else:
        
        if path.exists(line_regions_path) == False:
            
            logging.info("ERROR: The line regions file doesn't exist!")            
            sys.exit()
            
        else:
            
            # --- Read line regions ---
            
            logging.info("Reading line regions...")
            
            line_regions = ispec.read_line_regions(line_regions_path)
            line_regions = ispec.adjust_linemasks(n_star_spectrum, line_regions, max_margin=0.05, check_derivatives=True)
            
            # --- Create line segments ---
            
            logging.info("Creating line segments...")
            
            segments = ispec.create_segments_around_lines(line_regions, margin=0.10)        

    # -----------------------------------------------------------------------------#

    fit_lines_determine_ew_and_crossmatch_with_atomic_data(telluric_linelist, modeled_layers_pack, \
                                                           atomic_linelist, chemical_elements, isotopes, solar_abundances)
    sys.exit()
    # -----------------------------------------------------------------------------#
    
    if simple_calculation == True:

        temp_vsini = initial_vsini

        for vsini in frange(temp_vsini-1.0, temp_vsini+1.0, 1.0):
        #for vsini in frange(temp_vsini, temp_vsini, 1.0): # Only central value
            
            initial_vsini = vsini
            
            # --- Calculate stellar parameters ---
            
            logging.info("Calculting stellar parameters for %s..." % star)
            
            logging.info("The initial parameters are: Teff = %.2f, log g = %.2f, [M/H] = %.2f, vmic = %.2f, vmac = %.2f & vsini = %.2f " \
                         % (initial_teff, initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini))
            
            final_teff, e_teff, final_logg, e_logg, final_MH, e_MH, final_alpha, e_alpha, final_vmic, e_vmic, final_vmac, e_vmac, final_vsini, e_vsini, final_FeH, e_FeH, rms = \
                calculate_stellar_parameters(n_star_spectrum, n_star_continuum_model, initial_R, line_regions, segments, \
                                             modeled_layers_pack, atomic_linelist, chemical_elements, isotopes, solar_abundances)
            
            logging.info("The stellar parameter are:")    
            logging.info("Teff = %.2f +/- %.2f" % (final_teff, e_teff))
            logging.info("log g = %.2f +/- %.2f" % (final_logg, e_logg))
            logging.info("[M/H] = %.2f +/- %.2f" % (final_MH, e_MH))
            logging.info("alpha = %.2f +/- %.2f" % (final_alpha, e_alpha))
            logging.info("vmic = %.2f +/- %.2f" % (final_vmic, e_vmic))
            logging.info("vmac = %.2f +/- %.2f" % (final_vmac, e_vmac))
            logging.info("vsini = %.2f +/- %.2f" % (final_vsini, e_vsini))
            logging.info("[Fe/H] = %.2f +/- %.2f" % (final_FeH, e_FeH))
            logging.info("rms = %.6f" % rms)
            
            # --- Save results ---                
            
            params_results = params_results.append({'teff':final_teff, 'e_teff':e_teff, \
                                                    'logg':final_logg, 'e_logg':e_logg, \
                                                    'MH':final_MH, 'e_MH':e_MH, \
                                                    'alpha':final_alpha, 'e_alpha':e_alpha, \
                                                    'vmic':final_vmic, 'e_vmic':e_vmic, \
                                                    'vmac':final_vmac, 'e_vmac':e_vmac, \
                                                    'vsini':final_vsini, 'e_vsini':e_vsini, \
                                                    'FeH':''.join(format(x, ".2f") for x in final_FeH), 'e_FeH':''.join(format(x, ".2f") for x in e_FeH), \
                                                    'rms':rms}, ignore_index=True)

    else:

        temp_teff = initial_teff
        temp_logg = initial_logg
        temp_MH = initial_MH
        
        for teff in frange(temp_teff-100, temp_teff+100, 100):
            
            initial_teff = teff
            
            for logg in frange(temp_logg-0.5, temp_logg+0.5, 0.5):
                
                initial_logg = logg
                
                for MH in frange(temp_MH-0.25, temp_MH+0.25, 0.25):
                    
                    initial_MH = MH
                    
                    # Empirical relation for vmic & vmac
        
                    initial_vmic = ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
                    initial_vmac = ispec.estimate_vmac(initial_teff, initial_logg, initial_MH)
                    
                    # --- Calculate stellar parameters ---
                    
                    logging.info("Calculting stellar parameters for %s..." % star)
                    
                    logging.info("The initial parameters are: Teff = %.2f, log g = %.2f, [M/H] = %.2f, vmic = %.2f, vmac = %.2f & vsini = %.2f " \
                                 % (initial_teff, initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini))                            
                    
                    final_teff, e_teff, final_logg, e_logg, final_MH, e_MH, final_alpha, e_alpha, final_vmic, e_vmic, final_vmac, e_vmac, final_vsini, e_vsini, final_FeH, e_FeH, rms = \
                        calculate_stellar_parameters(n_star_spectrum, n_star_continuum_model, initial_R, line_regions, segments, \
                                                     modeled_layers_pack, atomic_linelist, chemical_elements, isotopes, solar_abundances)
                    
                    logging.info("The stellar parameter are:")    
                    logging.info("Teff = %.2f +/- %.2f" % (final_teff, e_teff))
                    logging.info("log g = %.2f +/- %.2f" % (final_logg, e_logg))
                    logging.info("[M/H] = %.2f +/- %.2f" % (final_MH, e_MH))
                    logging.info("alpha = %.2f +/- %.2f" % (final_alpha, e_alpha))
                    logging.info("vmic = %.2f +/- %.2f" % (final_vmic, e_vmic))
                    logging.info("vmac = %.2f +/- %.2f" % (final_vmac, e_vmac))
                    logging.info("vsini = %.2f +/- %.2f" % (final_vsini, e_vsini))
                    logging.info("[Fe/H] = %.2f +/- %.2f" % (final_FeH, e_FeH))
                    logging.info("rms = %.6f" % rms)
                    
                    # --- Save results ---                
                    
                    params_results = params_results.append({'teff':final_teff, 'e_teff':e_teff, \
                                                            'logg':final_logg, 'e_logg':e_logg, \
                                                            'MH':final_MH, 'e_MH':e_MH, \
                                                            'alpha':final_alpha, 'e_alpha':e_alpha, \
                                                            'vmic':final_vmic, 'e_vmic':e_vmic, \
                                                            'vmac':final_vmac, 'e_vmac':e_vmac, \
                                                            'vsini':final_vsini, 'e_vsini':e_vsini, \
                                                            'FeH':''.join(format(x, ".2f") for x in final_FeH), 'e_FeH':''.join(format(x, ".2f") for x in e_FeH), \
                                                            'rms':rms}, ignore_index=True)
            
    # -----------------------------------------------------------------------------#
            
    logging.info("Saving results for %s..." % star)
    
    params_results = params_results.sort_values(by='rms', ascending=True)
    params_results['rms'] = params_results['rms'].map(lambda x: '%.6f' % x)

    if estimate_vmic_and_vmac == True:

        params_results_filename = "ET_%s" % params_results_filename

    if estimate_vmic_and_vmac == False:

        params_results_filename = "DT_%s" % params_results_filename
    
    if calculate_logg == True:
        
        params_results_filename = "CG_%s" % params_results_filename
        
    if calculate_logg == False:

        params_results_filename = "DG_%s" % params_results_filename
        
    if code == "spectrum":

        codename = "S"

    if code == "turbospectrum":

        codename = "T"

    if code == "synthe":

        codename = "ST"

    if code == "moog":
        
        codename = "M"            
        
    params_results_filename = "%s_%s" % (codename, params_results_filename)
        
    if estimate_initial_parameters == True:

        params_results_filename = "A_%s" % params_results_filename

    if estimate_initial_parameters == False:

        params_results_filename = "D_%s" % params_results_filename
    
    params_results.to_csv ("./%s/%s" % (results_dirname, params_results_filename), index = False, header=True, sep='\t', float_format='%.2f')
    
    sys.exit()
