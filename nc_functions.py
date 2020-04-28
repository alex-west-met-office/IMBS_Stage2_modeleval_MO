'''This module contains all 'administrative' netCDF functions used in the production of monthly mean energy flux data from the IMBs'''

import netCDF4 as nc4

def create_monthly_file(writefile,vars_out_salin,vars_out_non_salin,salinity_values,parameters):
    '''Creates netCDF dataset for monthly mean output data with all required variables and dimensions. Output variables have two dimensions: a record dimension (with associated month, year and buoy number variables) and a salinity dimension, as salinity affects most of the output fluxes strongly via conductivity and heat capacity'''
    ncid_out = nc4.Dataset(writefile,'w')
    record_dim = ncid_out.createDimension('record_number')
    record_var = ncid_out.createVariable('record_number','int',dimensions=('record_number',))
    buoy_var = ncid_out.createVariable('buoy_number','int',dimensions=('record_number',))
    month_var = ncid_out.createVariable('month','int',dimensions=('record_number',))
    month = ncid_out.createVariable('year','int',dimensions=('record_number',))
    start_record = 0

    salinity_dim = ncid_out.createDimension('salinity',salinity_values.shape[0])
    salinity_var = ncid_out.createVariable('salinity','float64',\
        dimensions=('salinity',))
    salinity_var[:] = salinity_values

    for var in vars_out_non_salin:
        output_var = ncid_out.createVariable(var,'float64',dimensions=('record_number',))
        output_var_err = ncid_out.createVariable(var+'_err','int',dimensions=('record_number',))

    for var in vars_out_salin:
        output_var = ncid_out.createVariable(var,'float64',dimensions=('record_number','salinity'))
        output_var_err = ncid_out.createVariable(var+'_err','int',dimensions=('record_number','salinity'))

    for key in dir(parameters):
        if key[:2] != '__':
            ncid_out.setncattr(key,getattr(parameters,key))

    return ncid_out

def create_daily_file(writefile,vars_out,salinity_values):
    ncid_out = nc4.Dataset(writefile,'w')
    record_dim = ncid_out.createDimension('record_number')
    record_var = ncid_out.createVariable('record_number','int',dimensions=('record_number',))
    buoy_var = ncid_out.createVariable('buoy_number','int',dimensions=('record_number',))
    day_var = ncid_out.createVariable('day','int',dimensions=('record_number',))
    month_var = ncid_out.createVariable('month','int',dimensions=('record_number',))
    year_var = ncid_out.createVariable('year','int',dimensions=('record_number',))
    start_record = 0

    salinity_dim = ncid_out.createDimension('salinity',salinity_values.shape[0])
    salinity_var = ncid_out.createVariable('salinity','float64',\
        dimensions=('salinity',))
    salinity_var[:] = salinity_values

    for var in vars_out:
        output_var = ncid_out.createVariable(var,'float64',dimensions=('record_number','salinity'))
        output_var_err = ncid_out.createVariable(var+'_err','int',dimensions=('record_number','salinity'))

    return ncid_out      
