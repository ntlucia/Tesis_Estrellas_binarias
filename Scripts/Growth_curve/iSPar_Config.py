#iSPar_Config.py:

# --- iSpec directory ---

ispec_dir = "virtual/home/iSpec/" # iSpec directory

# --- Spectrum properties ---

n_spectra = 1 # 1. The best spectrum, 2. The second best spectra, etc (6 maximum)
coadd_spectra = False # If coadd spectra is true, n_spectra must be between 2 to 6
channel = "B" # Channel of HEROS B or R
wave_base = 375.0 # Lower wavelength to cut spectrum. for R channel 580.0 and B channel 375.0
wave_top = 570.0 # Uper wavelength to cut spectrum: for R channel 870.0 and B channel 570.0
defined_resolution = 21000 # Spectrum resolution defined by user or
use_spectrum_resolution = True # Use the resolution inside of spectrum header

# The line selection was built based on a solar and stellar spectrum with R=20,000 and
# VALD atomic linelist. Use this function to downgrade the resolution with high resolution spectra

resolution_degradation = False
final_R = 20000 # Recommended for TIGRE spectra

# --- Model parameters ---

precompute_synth_grid = False
precomputed_grid_dirname = "TURBOSPECTRUM_MARCS_Grevesse.2007_VALD.480_680nm"

# Calculate the initial parameters or use those defined below

estimate_initial_parameters = False

# Define initial stellar parameters

defined_teff = 4052
defined_logg = 1.02
defined_MH = 0.03
defined_alpha = 0.0
defined_vsini = 1.6

# Use empirical relations to estimate the initial vmic and vmac...

estimate_vmic_and_vmac = True

# Or use vmic and vmac defined below

defined_vmic = 1.0
defined_vmac = 3.0

# --- Calculation properties ---

# Calculate log g or use the value defined above (it won't be calculated)

calculate_logg = True

# Use errors in the spectrum file to calculate stellar parameters

use_errors = True

# Number of iterations to calculate stellar parameters

max_iterations = 10

# Simple calculation or variation around initial parameters (for Teff: 100, log g: 0.5, [M/H]: 0.1)

simple_calculation = True

# Save all normalizations and residuals or only the best normalization

save_all_residuals = False

# Find line regions (calculating or not the individual abundances) or define the line regions path

find_line_regions = False
calculate_individual_abundances = False
#line_regions_path = "%s/input/regions/iSPar_Line_Regions_Filtered_Moon.dat" % ispec_dir
#line_regions_path = "%s/input/regions/iSPar_Line_Regions_Filtered_HD124897.dat" % ispec_dir
line_regions_path = "%s/input/regions/iSPar_Line_Regions_Filtered_Combined.dat" % ispec_dir
#line_regions_path = ispec_dir + "/input/regions/47000_VALD/turbospectrum_synth_good_for_params_all.txt"

# Radiative transfer code

#code = "spectrum" # 1D plane-parallel model atmosphere geometry, written in C
#code = "synthe" # 1D plane-parallel model atmosphere geometry,  written in FORTRAN but precompiled (Intel)
#code = "moog" # 1D plane-parallel model atmosphere geometry, written in FORTRAN
code = "turbospectrum" # 1D spherical model atmosphere geometry, written in FORTRAN

# Stellar atmosphere model

#atm_model = "%s/input/atmospheres/MARCS/" % ispec_dir
atm_model = "%s/input/atmospheres/MARCS.GES/" % ispec_dir
#atm_model = "%s/input/atmospheres/MARCS.APOGEE/" % ispec_dir
#atm_model = "%s/input/atmospheres/ATLAS9.Castelli/" % ispec_dir
#atm_model = "%s/input/atmospheres/ATLAS9.Kurucz/" % ispec_dir

# Solar abundances

#solar_abundances_file = "%s/input/abundances/Grevesse.1998/stdatom.dat" % ispec_dir
solar_abundances_file = "%s/input/abundances/Grevesse.2007/stdatom.dat" % ispec_dir
#solar_abundances_file = "%s/input/abundances/Asplund.2005/stdatom.dat" % ispec_dir
#solar_abundances_file = "%s/input/abundances/Asplund.2009/stdatom.dat" % ispec_dir

# Atomic linelist

atomic_linelist_file = "%s/input/linelists/transitions/VALD.300_1100nm/atomic_lines.tsv" % ispec_dir
#atomic_linelist_file = "%s/input/linelists/transitions/VALD.1100_2400nm/atomic_lines.tsv" % ispec_dir
#atomic_linelist_file = "%s/input/linelists/transitions/GESv5_atom_hfs_iso.420_920nm/atomic_lines.tsv" % ispec_dir
#atomic_linelist_file = "%s/input/linelists/transitions/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.tsv" % ispec_dir

# Free parameters

enhance_abundances = False

# Free individual element abundance
    
free_abundances = None
linelist_free_loggf = None
