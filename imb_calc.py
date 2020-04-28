'''This module contains all scientific equations and operations necessary for calclating monthly mean energy fluxes from the IMB data.'''

import scientific_constants as sc
import parameters as param
import numpy as np
import difference_functions as df

def region_code(ncid_in,index,aux_vars_non_salin):
    '''Determines whether a buoy was contained entirely with the North Pole or Beaufort Sea regions within a particular month of operation'''
    import sys
    sys.path.insert(0,'/home/h01/hadax/Python/PhD/Buoys/hadax-IceMassBalanceBuoys/')
    import spatial
    region_codes = [spatial.region_key([ncid_in.variables['longitude_rt'][:][index][item],ncid_in.variables['latitude_rt'][:][index][item]]) for item in range(index[0].shape[0])]
    if np.unique(np.array(region_codes)).shape[0] == 1:
        return np.unique(np.array(region_codes))
    else:
        return 0

def latitude(ncid_in,index,aux_vars_non_salin):
    '''Returns average latitude for month of operation of a particular buoy'''
    try:
        return np.mean(ncid_in.variables['latitude_rt'][:][index])
    except KeyError:
        return None

def longitude(ncid_in,index,aux_vars_non_salin):
    '''Returns average longitude for month of operation of a particular buoy. Should be able to handle instances of crossing the date line.'''
    try:
        lon_array = ncid_in.variables['longitude_rt'][:][index]
    except KeyError:
        return None

    maxlon = np.max(lon_array)
    minlon = np.min(lon_array)

    if minlon < -90 and maxlon > 90:
        index = np.where(lon_array < 0.)
        lon_array[index] = lon_array[index] + 360.

    if minlon < 90 and maxlon > 270:
        index = np.where(lon_array > 180.)
        lon_array[index] = lon_array[index] - 360.

    return np.mean(lon_array)

'''Here follow a series of functions which calculate the monthly mean variable with the same name, given auxiliary variables'''

def ice_thickness(ncid_in,index,aux_vars_non_salin,aux_vars_salin):
    try:
        return df.mmean(ncid_in.variables['ice_thickness'][:][index])
    except KeyError:
        return None

def snow_thickness(ncid_in,index,aux_vars_non_salin,aux_vars_salin):
    try:
        return df.mmean(ncid_in.variables['snow_thickness'][:][index])
    except KeyError:
        return None

def thermal_insulance(ncid_in,index,aux_vars_non_salin,aux_vars_salin):
    try:
        return df.mmean(ncid_in.variables['ice_thickness'][:][index] / sc.k_fresh_ice_mu + ncid_in.variables['snow_thickness'][:][index] / sc.k_snow)
    except KeyError:
        return None

def sfctemp(ncid_in,index,aux_vars_non_salin,aux_vars_salin):
    try:
        return df.mmean(ncid_in.variables['surface_temp_-0.000'][:][index])
    except KeyError:
        return None

def fcondtop(ncid_in,index,aux_vars_non_salin,aux_vars_salin):
    try:
        return df.mmean(ncid_in.variables['zgrad_top_snow_adjusted_-0.500'][:][index] * aux_vars_salin['top_conductivity_estimate'])
    except KeyError:
        return None

def fcondbot(ncid_in,index,aux_vars_non_salin,aux_vars_salin):
    try:
        return df.mmean(ncid_in.variables['base_'+param.layer_fcondbot_label+'_zgrad'][:][index] * aux_vars_salin['basal_conductivity_estimate'])
    except KeyError:
        return None
    
def sensible_heat_uptake_lowest_layer(ncid_in,index,aux_vars_non_salin,aux_vars_salin):
    try:
        return df.mmean(aux_vars_non_salin['temp_rate_of_change'] * aux_vars_salin['heat_capacity'] * param.rhoi) * (param.layer_shu_calc[1] - param.layer_shu_calc[0])
    except KeyError:
        return None

def basal_melt(ncid_in,index,aux_vars_non_salin,aux_vars_salin):
    try:
        return -1. * aux_vars_non_salin['basal_melt_mass'] * aux_vars_salin['energy_of_melting_base'] / aux_vars_non_salin['seconds_in_month']
    except KeyError:
        return None

def basal_growth(ncid_in,index,aux_vars_non_salin,aux_vars_salin):
    try:
        return -1. * aux_vars_non_salin['basal_growth_mass'] * aux_vars_salin['energy_of_melting_base'] / aux_vars_non_salin['seconds_in_month']
    except KeyError:
        return None

def ocean_heat_flux(*args):
    return sensible_heat_uptake_lowest_layer(*args) + basal_melt(*args) + basal_growth(*args) - fcondbot(*args)

def top_melt(ncid_in,index,aux_vars_non_salin,aux_vars_salin):
    try:
        return aux_vars_non_salin['top_melt_mass'] * aux_vars_salin['energy_of_melting_top'] / aux_vars_non_salin['seconds_in_month']
    except KeyError:
        return None

def top_decrease_not_melt(ncid_in,index,aux_vars_non_salin,aux_vars_salin):
    try:
        return aux_vars_non_salin['top_decrease_not_melt_mass'] * aux_vars_salin['energy_of_melting_top'] / aux_vars_non_salin['seconds_in_month']
    except KeyError:
        return None

def top_increase(ncid_in,index,aux_vars_non_salin,aux_vars_salin):
    try:
        return aux_vars_non_salin['top_increase_mass'] * aux_vars_salin['energy_of_melting_top'] / aux_vars_non_salin['seconds_in_month']
    except KeyError:
        return None

def calc_aux_vars_non_salin(ncid_in,period,buoy_dates):
    '''Calculates auxiliary variables, not dependent upon salinity, necessary for calculating energy fluxes, e.g. elevation changes'''
    import difference_functions as df
    import numpy as np

    index = np.where(np.logical_and(buoy_dates > period[0],buoy_dates <= period[1]))
    aux_vars_non_salin = {}

    aux_vars_non_salin['ice_snow_column_thickness'] = ncid_in.variables['ice_thickness'][:] + ncid_in.variables['snow_thickness'][:]
    aux_vars_non_salin['ice_too_thin_logical_for_fcondbot'] = ncid_in.variables['ice_thickness'][:] < param.thin_ice_threshold_for_fcondbot
    aux_vars_non_salin['column_too_thin_logical_for_fcondtop'] = aux_vars_non_salin['ice_snow_column_thickness'] < param.thin_snow_ice_threshold
    aux_vars_non_salin['ice_too_thin_logical_for_shu'] = ncid_in.variables['ice_thickness'][:] < param.thin_ice_threshold_for_shu

    aux_vars_non_salin['seconds_in_month'] = sc.seconds_in_day * (period[1].toordinal()-period[0].toordinal())
    aux_vars_non_salin['basal_melt_mass'] = df.sum_increase(ncid_in.variables['bottom_rt'][:]*param.rhoi,period,\
        buoy_dates)
    aux_vars_non_salin['basal_growth_mass'] = df.sum_decrease(ncid_in.variables['bottom_rt'][:]*param.rhoi,period,\
        buoy_dates)

    aux_vars_non_salin['top_difference'] = df.data_difference(ncid_in.variables['surface_rt'][:],period,buoy_dates)
    aux_vars_non_salin['top_roc_positive'] = aux_vars_non_salin['top_difference'] * (aux_vars_non_salin['top_difference'] >= 0.)
    aux_vars_non_salin['top_roc_negative'] = aux_vars_non_salin['top_difference'] * (aux_vars_non_salin['top_difference'] < 0.)

    aux_vars_non_salin['snow_indicator'] = ncid_in.variables['surface_rt'][index] > ncid_in.variables['interface_rt'][index] + param.snow_threshold_for_topmelt
    aux_vars_non_salin['rho_array'] = aux_vars_non_salin['snow_indicator'].astype('int') * param.rhos + (1. - aux_vars_non_salin['snow_indicator'].astype('int')) * param.rhoi
    aux_vars_non_salin['rho_array_c'] = ((aux_vars_non_salin['rho_array'] + np.roll(aux_vars_non_salin['rho_array'],1)) / 2.)[1:]
    aux_vars_non_salin['sfctemp_c'] = ((ncid_in.variables['surface_temp_-0.000'][index] + np.roll(ncid_in.variables['surface_temp_-0.000'][index],1)) / 2.)[1:]
    aux_vars_non_salin['sfctemp_cold_c'] = aux_vars_non_salin['sfctemp_c'] < param.sfctemp_melt_threshold

    aux_vars_non_salin['top_melt_mass'] = df.msum(aux_vars_non_salin['top_roc_negative']*aux_vars_non_salin['rho_array_c']*(1.-aux_vars_non_salin['sfctemp_cold_c'].astype('int')))
    aux_vars_non_salin['top_decrease_not_melt_mass'] = df.msum(aux_vars_non_salin['top_roc_negative']*aux_vars_non_salin['rho_array_c']*aux_vars_non_salin['sfctemp_cold_c'].astype('int'))
    aux_vars_non_salin['top_increase_mass'] = df.msum(aux_vars_non_salin['top_roc_positive']*aux_vars_non_salin['rho_array_c'])

    aux_vars_non_salin['temp_change'] = df.data_difference(ncid_in.variables['base_'+param.layer_shu_label+'_zmean'][:],period,buoy_dates)
    aux_vars_non_salin['temp_rate_of_change'] = df.rate_of_change_polyfit(ncid_in.variables['base_'+param.layer_shu_label+'_zmean'][:],period,buoy_dates)
   
    return aux_vars_non_salin

def calc_aux_vars_salin(ncid_in,period,buoy_dates,salinity):
    '''Calculates auxiliary variables, dependent upon salinity, necessary for calculating vertical energy fluxes (e.g. conductivity, heat capacity)'''
    import difference_functions as df
    import numpy as np
    import copy

    aux_vars_salin = {}
    index = np.where(np.logical_and(buoy_dates > period[0],buoy_dates <= period[1]))

    aux_vars_salin['melting_temperature'] = 0. - salinity * sc.liquidus_constant

    max_top_temp_array = copy.deepcopy(ncid_in.variables['surface_temp_-0.000'][:][index])
    top_exceed_index = np.where(max_top_temp_array > aux_vars_salin['melting_temperature'])
    max_top_temp_array[top_exceed_index] = aux_vars_salin['melting_temperature']

    max_base_temp_array = copy.deepcopy(ncid_in.variables['basal_temp_0.000'][:][index])
    base_exceed_index = np.where(max_base_temp_array > aux_vars_salin['melting_temperature'])
    max_base_temp_array[base_exceed_index] = aux_vars_salin['melting_temperature']

    if salinity != 0.:
        aux_vars_salin['energy_of_melting_base'] = - sc.L_fresh * (1 + sc.liquidus_constant * salinity / df.mmean(max_base_temp_array))
    else:
        aux_vars_salin['energy_of_melting_base'] = -sc.L_fresh
    if salinity != 0.:
        aux_vars_salin['energy_of_melting_top'] = - sc.L_fresh * (1 + sc.liquidus_constant * salinity / df.mmean(max_top_temp_array))
    else:
        aux_vars_salin['energy_of_melting_top'] = -sc.L_fresh

    aux_vars_salin['top_conductivity_estimate_mu'] = sc.k_fresh_ice_mu  + salinity * sc.beta_k_coefficient / df.mmean(max_top_temp_array)
    aux_vars_salin['basal_conductivity_estimate_mu'] = sc.k_fresh_ice_mu  + salinity * sc.beta_k_coefficient * ncid_in.variables['base_'+param.layer_fcondbot_label+'_z_recip_mean'][:][index]

    aux_vars_salin['top_conductivity_estimate_pringle'] = sc.k_fresh_ice_pringle  - sc.pringle_temp_cond * df.mmean(max_top_temp_array) + sc.pringle_saltemp_cond * salinity / df.mmean(max_top_temp_array)
    aux_vars_salin['basal_conductivity_estimate_pringle'] = sc.k_fresh_ice_pringle  - sc.pringle_temp_cond * df.mmean(max_top_temp_array) + sc.pringle_saltemp_cond * salinity / df.mmean(max_top_temp_array)

    if param.conductivity_formulation == 'maykut_untersteiner':
        aux_vars_salin['top_conductivity_estimate'] = \
            aux_vars_salin['top_conductivity_estimate_mu']
        aux_vars_salin['basal_conductivity_estimate'] = \
            aux_vars_salin['basal_conductivity_estimate_mu']
    elif param.conductivity_formulation == 'pringle_bubbly_brine':
        aux_vars_salin['top_conductivity_estimate'] = \
            aux_vars_salin['top_conductivity_estimate_pringle']
        aux_vars_salin['basal_conductivity_estimate'] = \
            aux_vars_salin['basal_conductivity_estimate_pringle']
    else:
        print('Conductivity formulation must be \'maykut_untersteiner\' or '+\
             '\'pringle_bubbly_brine\'')
        raise ValueError

    aux_vars_salin['heat_capacity'] = sc.c_fresh_ice + sc.gamma_c_coefficient * salinity * df.mmean(ncid_in.variables['base_'+param.layer_shu_label+'_z_recip_sq_mean'][:][index])

    # Calculate required additional error quantities
    aux_vars_salin['temps_too_high_in_lowest_layer'] = -1. / np.sqrt(ncid_in.variables['base_0.000_0.400_z_recip_sq_mean'][:]) > aux_vars_salin['melting_temperature']
    #most_temps_too_high = np.sum(ncid_in.variables['surface_temp_-0.000'][:][index] > aux_vars_salin['melting_temperature'] + param.melting_temp_tol_most) > index[0].shape[0] * param.melting_temp_most_proportion
    #aux_vars_salin['temp_too_high'] = all_temps_too_high or most_temps_too_high

    return aux_vars_salin



def handle_errors(ncid_out,record_number,isalin,error_codes,aux_vars_salin,
                  aux_vars_non_salin,index,missing_data,false_bottoms,
                  temps_too_high_posthoc,other_problems):
    '''Control for various common errors, or problems, in the processed IMB data that mean it is impossible, or undesirable, to produce a monthly mean estimate of a particular variable. The error code for that month, salinity and variable is set to a particular nonzero value that will subsequently prevent the data from being read into the output array'''
        
    #if True in aux_vars_non_salin['column_too_thin_logical_for_fcondtop'][index]:
    #        ncid_out.variables['fcondtop_err'][record_number,isalin] = error_codes['column_too_thin']
        
    #if True in aux_vars_non_salin['ice_too_thin_logical_for_fcondbot'][index]:
    #        ncid_out.variables['fcondbot_err'][record_number,isalin] = error_codes['ice_too_thin_for_fcondbot']
    #        ncid_out.variables['ocean_heat_flux_err'][record_number,isalin] = error_codes['ice_too_thin_for_fcondbot']

    #if True in aux_vars_non_salin['ice_too_thin_logical_for_shu'][index]:
    #        ncid_out.variables['sensible_heat_uptake_lowest_layer_err'][record_number,isalin] = error_codes['ice_too_thin_for_sensible_heat_uptake']
    #        ncid_out.variables['ocean_heat_flux_err'][record_number,isalin] = error_codes['ice_too_thin_for_sensible_heat_uptake']
 
    #if True in aux_vars_salin['temps_too_high_in_lowest_layer'][index]:
    #        ncid_out.variables['sensible_heat_uptake_lowest_layer_err'][record_number,isalin] = error_codes['temps_too_high_for_this_salinity']
    #        ncid_out.variables['ocean_heat_flux_err'][record_number,isalin] = error_codes['temps_too_high_for_this_salinity']
            
    for varname in [*missing_data.keys()]:
        if missing_data[varname]:
            ncid_out.variables[varname+'_err'][record_number,isalin] = error_codes['missing_data']

    record_vars = ['buoy_number','month','year']
    if tuple([ncid_out.variables[rvar][record_number] for rvar in record_vars]) in false_bottoms['ocean_heat_flux']:
        ncid_out.variables['ocean_heat_flux_err'][record_number,isalin] = error_codes['false_bottom_occurrence']

    if tuple([ncid_out.variables[rvar][record_number] for rvar in record_vars]) in other_problems['ocean_heat_flux']:
        ncid_out.variables['ocean_heat_flux_err'][record_number,isalin] = error_codes['other_problems']

    for varname in ['fcondtop','fcondbot','sensible_heat_flux_lowest_layer','ocean_heat_flux']:
        try:
            if tuple([ncid_out.variables[rvar][record_number] for rvar in record_vars]) in temps_too_high_posthoc[varname]:
                if ncid_out.variables['salinity'][isalin] != 0.:
                    ncid_out.variables[varname+'_err'][record_number,isalin] = error_codes['temps_too_high_posthoc']
        except KeyError:
            pass
