'''This module contains all functions used for subsetting IMB records by time'''

def start_datetime_daily(year,month,day):
    '''Returns datetime at the start of a particular day'''
    import datetime as dt
    return dt.datetime(year,month,day,0,0,0)

def end_datetime_daily(year,month,day):
    '''Returns datetime at the end of a particular day'''
    import datetime as dt
    ddate = dt.date(year,month,day)
    number = ddate.toordinal() + 1
    ddate_next = dt.date.fromordinal(number)
    return dt.datetime(ddate_next.year,ddate_next.month,ddate_next.day,0,0,0)

def start_datetime(year,month):
    '''Returns datatime at the start of a particular month'''
    import datetime as dt
    return dt.datetime(year,month,1,0,0)

def end_datetime(year,month):
    '''Returns datetime at the end of a particular month'''
    import datetime as dt
    return dt.datetime(year+(month==12),month%12+1,1,0,0)

def ym_in_period(year,month,period):
    import datetime as dt
    dec_log = month==12
    if dec_log:
        emonth = 1
        eyear = year + 1
    else:
        emonth = month + 1
        eyear = year
    
    after_start = dt.datetime(year,month,1,0,0) >= period[0]
    before_end = dt.datetime(eyear,emonth,1,0,0) <= period[1]    
    return (after_start and before_end)
    
