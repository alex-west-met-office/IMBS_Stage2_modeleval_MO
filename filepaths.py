'''Defines files and directories used in stage 2 of the IMB calculation'''

def filepaths():
    output = {
        'input_dir' : # Directory where all IMB data series produced by Stage 1 of the code are stored in netCDF (equivalent to 'output_dir' in stage 1)
        'monthly_mean_output_file' : # Directory where the dataset of monthly mean energy fluxes is stored
        'false_bottoms' : # List of monthly mean energy fluxes likely to be corrupted by false bottom formation. Provided with the code
        'temps_too_high_posthoc' : # List of monthly mean energy fluxes likely to be corrupted, in the case of nonzero salinity being assumed, by the presence of temperatures above the melting point. Provided with the code
        'other_problems' : # List of monthly mean energy fluxes likely to be corrupted for other reasons, given in the file. Provided with the code
        }
    return output
