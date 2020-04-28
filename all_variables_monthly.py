'''This is the driving routine for calculating monthly mean values of energy fluxes from the processed IMB data'''

# Imports
import aux_functions
import netCDF4 as nc4
import numpy as np
import numpy.ma as ma
import os
import difference_functions as df
import time_functions as tf
import imb_calc
import nc_functions
import parameters
import filepaths


# Control settings
year_ref = 1978
salinity_values = np.arange(0,10.1,.5)
force_new = False
version = 'v4_clean'
buoy_list = aux_functions.buoylist()

writefile = filepaths.filepaths()['monthly_mean_output_file']

error_codes = {
    'missing_data'                          : 99,
    'temps_too_high_for_this_salinity'      : 1,
    'ice_too_thin_for_fcondbot'             : 2,
    'ice_too_thin_for_sensible_heat_uptake' : 3,
    'column_too_thin'                       : 4,
    'false_bottom_occurrence'               : 5,
    'temps_too_high_posthoc'                : 6,
    'other_problems'                        : 7}

vars_out_salin = ['fcondtop','fcondbot','top_increase','top_decrease_not_melt','top_melt','basal_growth','basal_melt','sensible_heat_uptake_lowest_layer','ice_thickness','snow_thickness','thermal_insulance','sfctemp','ocean_heat_flux']

vars_out_non_salin = ['latitude','longitude','region_code']

false_bottoms = aux_functions.read_error_file('false_bottoms')
temps_too_high_posthoc = aux_functions.read_error_file('temps_too_high_posthoc')
other_problems = aux_functions.read_error_file('other_problems')

for buoy_name in buoy_list:
    print('Calculating monthly variables for buoy '+buoy_name)

    datafile = filepaths.filepaths()['input_dir'] + \
        buoy_name + '_processed.nc'

    ncid_in = nc4.Dataset(datafile)
    time = ncid_in.variables['time']
    buoy_dates = nc4.num2date(time[:],time.units)
    full_period = [buoy_dates[0],buoy_dates[-1]]
    months = np.array([ddt.month for ddt in nc4.num2date(time[:],time.units)])
    years = np.array([ddt.year for ddt in nc4.num2date(time[:],time.units)])
    month_nos = months + (years - year_ref) * 12
    unique_month_nos = np.unique(month_nos)


    if force_new or not os.path.exists(writefile):
        ncid_out = nc_functions.create_monthly_file(writefile,vars_out_salin,vars_out_non_salin,salinity_values,parameters)
        start_record = 0

    else:
        ncid_out = nc4.Dataset(writefile,'a')
        start_record = ncid_out.variables['record_number'].shape[0]


    for (irec,month_no) in enumerate(unique_month_nos):
        month = (month_no-1) % 12 + 1
        year = (month_no-1) // 12 + year_ref
        if tf.ym_in_period(year,month,full_period):
            period = [tf.start_datetime(year,month),tf.end_datetime(year,month)]
            index = np.where(np.logical_and(buoy_dates > period[0],buoy_dates <= period[1]))
            record_number = start_record + irec
            ncid_out.variables['record_number'][record_number] = record_number
            ncid_out.variables['buoy_number'][record_number] = aux_functions.buoynumber(buoy_name)
            ncid_out.variables['month'][record_number] = month
            ncid_out.variables['year'][record_number] = year

            for varname in vars_out_non_salin:
                ncid_out.variables[varname+'_err'][record_number] = 0.

            for varname in vars_out_salin:
                ncid_out.variables[varname+'_err'][record_number,:] = 0.

            # Calculate required compound variables not dependent on salinity
            aux_vars_non_salin = imb_calc.calc_aux_vars_non_salin(ncid_in,period,buoy_dates)

            for varname in vars_out_non_salin:
                funcname = getattr(imb_calc,varname)
                ncid_out.variables[varname][record_number] = funcname(ncid_in,index,aux_vars_non_salin)

            estimates = {}
            for (isalin,salinity) in enumerate(salinity_values):
                # Calculate salinity & time-dependent parameters
                aux_vars_salin = imb_calc.calc_aux_vars_salin(ncid_in,period,buoy_dates,salinity)

                # Calculate final quantities
                for varname in vars_out_salin:
                    funcname = getattr(imb_calc,varname)
                    estimates[varname] = funcname(ncid_in,index,aux_vars_non_salin,aux_vars_salin)

                missing_data = {varname: estimates[varname] is None for varname in vars_out_salin}

                imb_calc.handle_errors(ncid_out,record_number,isalin,error_codes,
                                       aux_vars_salin,aux_vars_non_salin,index,
                                       missing_data,false_bottoms,
                                       temps_too_high_posthoc,other_problems)

                for varname in vars_out_salin:
                    if ncid_out.variables[varname+'_err'][record_number,isalin] == 0 and np.isfinite(estimates[varname]):
                        ncid_out.variables[varname][record_number,isalin] = \
                          estimates[varname]
    

    ncid_in.close()
    
    ncid_out.close()
