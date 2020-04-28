'''This module contains all functions used for calculating changes in time of processed IMB data series. This is important for estimating elevation changes, and sensible heat uptake, when calculating monthly mean IMB fluxes'''


def data_difference(data,period,dates):
    '''Given a data series and a period in time, returns differences between successive data points in that time period (inclusive of end but not of start)'''
    import numpy as np
    index = np.where(np.logical_and(dates > period[0],dates <= period[1]))
    data_r = data[index]
    difference = data_r - np.roll(data_r,1)
    return difference[1:]   

def time_difference(period,dates): 
    '''Given a data series and a period in time, returns differences in time between the measurement points of successive data points in that time period (inclusive of end but not of start)'''
    import numpy as np
    times = [ddt for ddt in dates if ddt > period[0] and ddt <= period[1]]
    time_differences = np.array([(ddt2-ddt1).total_seconds() for (ddt1,ddt2) in \
         zip(times[:-1],times[1:])]) 
    return time_differences

def rate_of_change(data,period,dates):
    '''Given a data series and a period in time, returns estimate of rate of change between successive data points in that period (inclusive of end but not of start)'''

    data_differences = data_difference(data,period,dates)
    time_differences = time_difference(period,dates)
    return data_differences / time_differences 

def rate_of_change_polyfit(data,period,dates,fit_length=1):
    '''This function calculates a rate of change in the same way as the old code did, using successive lines of best fit'''
    import aux_functions
    import numpy as np
    import numpy.ma as ma

    seconds_in_day = 8.64e4
    tolerance = 0.95   
    numbers = np.array([aux_functions.datetime_to_float(ddt) for ddt in dates])
    dates_within = [ddt for ddt in dates if ddt >= period[0] and ddt <= period[1]]
    output_series = ma.masked_array(np.zeros(len(dates_within)),
            mask=np.zeros(len(dates_within),dtype='bool'))

    for (idate,ddt) in enumerate(dates_within):
        number = aux_functions.datetime_to_float(ddt)
        index_period = np.where(np.logical_and(
                                               numbers >= number-fit_length,
                                               numbers <= number+fit_length))
        numbers_period = numbers[index_period]
        data_period = data[index_period]
        
        if np.sum(data_period.mask) / data_period.shape[0] > tolerance:
            output_series.mask.itemset(idate,True)
        else:
            rate_of_change_instance = ma.polyfit(numbers_period,data_period,1)[0] / seconds_in_day
            output_series.data.itemset(idate,rate_of_change_instance)


    return output_series

def sum_increase(data,period,dates):
    '''Given a data series and a time period, returns the sum of the positive differences between successive data points in that time period'''
    import numpy as np
    difference = data_difference(data,period,dates)
    difference_positive = difference * (difference >= 0)
    return msum(difference_positive)

def sum_decrease(data,period,dates):
    '''Given a data series and a time period, returns the sum of the negative differences between successive data points in that time period'''
    import numpy as np
    difference = data_difference(data,period,dates)
    difference_negative = difference * (difference <= 0)
    return msum(difference_negative)

def sum_total(data,period,dates):
    '''Given a data series and a time period, returns the total change in that data series over the time period'''
    import numpy as np
    difference = data_difference(data,period,dates)
    return msum(difference)

def wmean(array,buoy_dates,period):
    '''Returns 'mean' value of a data series in a time period, where the average is weighted such that data points at the exact start and end point of the period are given half the weight of points strictly within the period'''
    import numpy as np
    import datetime as dt
    ieq = np.where(np.logical_or(buoy_dates==period[0],buoy_dates==period[1]))
    iwithin = np.where(np.logical_and(buoy_dates>period[0],buoy_dates<period[1]))

    neq = ieq[0].shape[0]
    nwithin = iwithin[0].shape[0]

    if neq+nwithin==0:
        print('No data from this period')
        return None

    weights = np.zeros(neq+nwithin)
    if neq>0:
        weights[:neq] = .5

    if nwithin>0:
        weights[-nwithin:] = 1.

    weights = weights / np.sum(weights)

    array_m = 0.
    if neq>0:
        array_m = array_m + np.sum(array[ieq] * weights[:neq])
    
    if nwithin>0.:
        array_m = array_m + np.sum(array[iwithin] * weights[-nwithin:])

    return array_m
        

def msum(array,tolerance=0.05):
    return mmean(array) * array.size     

def mmean(array,tolerance=0.05):
    import numpy as np
    import numpy.ma as ma

    if np.sum(array.mask) / array.shape[0] > tolerance:
        return ma.masked_array([0.],mask=[True])
    else:
        return ma.masked_array(np.mean(array),mask=[False])
