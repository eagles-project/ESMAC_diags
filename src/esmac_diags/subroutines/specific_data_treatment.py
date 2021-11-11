"""
function of some specific data treatment
"""

import numpy as np
from ..subroutines.quality_control import qc_remove_neg

#%%
def avg_time_1d(time0, data0, time):
    """
    average 1d data into coarser time resolution

    Parameters
    ----------
    time0 : numpy array
        time dimension for input data
    data0 : numpy array
        input data
    time : numpy array
        time dimension for output data

    Returns
    -------
    data : output data

    """
    data0 = qc_remove_neg(data0)
    if data0.shape[0] != len(time0):
        raise ValueError("Arrays must have the same size")
    data = np.full((len(time)), np.nan)
    dt = (time[1]-time[0])/2
    for tt in range(len(time)):
        idx = np.logical_and(time0 >= time[tt]-dt, time0 <= time[tt] + dt)
        data[tt] = np.nanmean(data0[idx], axis = 0)
    return(data)

# 
def avg_time_2d(time0, data0, time):
    """
    average 2d data into coarser time resolution

    Parameters
    ----------
    time0 : numpy array
        time dimension for input data
    data0 : numpy array
        input data
    time : numpy array
        time dimension for output data

    Returns
    -------
    data : output data

    """
    data0 = qc_remove_neg(data0)
    if data0.shape[0] != len(time0):
        raise ValueError("the first dimension of input data must have the same size with time")
    data = np.full((len(time), data0.shape[1]), np.nan)
    dt = (time[1]-time[0])/2
    for tt in range(len(time)):
        idx = np.logical_and(time0 >= time[tt]-dt, time0 <= time[tt] + dt)
        data[tt, :] = np.nanmean(data0[idx, :], axis = 0)
    return(data)

#%%
def lwc2cflag(lwc, lwcunit):
    """
    estimate cloud flag based on LWC

    Parameters
    ----------
    lwc : numpy array
        liquid water content data
    lwcunit : string
        unit of lwc

    Returns
    -------
    cldflag : estimated cloud flag

    """
    if lwcunit == 'kg/m3':
        lwc = lwc*0.001
        lwcunit = 'g/m3'
    elif lwcunit == 'g m-3':
        lwcunit = 'g/m3'
    if lwcunit not in ['g/m3', 'gram/m3']:
        print('unit of LWC: ' + lwcunit)
        raise ValueError("unit of LWC should be gram/m3 or kg/m3")

    cldflag = 0*np.array(lwc)
    
    # set threshold of LWC to identify cloud
    cldflag[lwc > 0.02] = 1
    return(cldflag)

#%%
def mask_model_ps(timem, psm, legnum, campaign, shipmetpath):
    """
    set model masks if the difference of Ps with observation is too large

    Parameters
    ----------
    timem : numpy array
        time in model
    psm : numpy array
        surface pressure in model
    legnum : string
        leg number,  or trip number
    campaign : string
        name of field campaign
    shipmetpath : string
        file path of shipmet data

    Returns
    -------
    datamask : mask flag of large Ps difference

    """
    import glob
    from ..subroutines.read_ship import read_marmet
    from ..subroutines.read_ARMdata import read_met
    from ..subroutines.time_format_change import yyyymmdd2cday,  cday2mmdd
    
    if campaign == 'MAGIC':
        filenameo = shipmetpath + 'marmet' + legnum + '.txt'
        (shipdata, shipvarlist) = read_marmet(filenameo)
        year = [a[1] for a in shipdata]
        month = [a[2] for a in shipdata]
        day = [a[3] for a in shipdata]
        hh = [int(a[4]) for a in shipdata]
        mm = [int(a[5]) for a in shipdata]
        ss = [int(a[6]) for a in shipdata]
        yyyymmdd = [year[i] + month[i] + day[i] for i in range(len(year))]   # yyyymmdd
        # get time in calendar day
        time = np.array(hh)/24. + np.array(mm)/1440. + np.array(ss)/86400. 
        time = np.array([time[i] + yyyymmdd2cday(yyyymmdd[i], 'noleap') for i in range(len(time))])
        if time[-1] < time[0]:
            time[time <= time[-1]] = time[time <= time[-1]] + 365
        # get variables
        ps = np.array([float(a[shipvarlist.index('bp')]) for a in shipdata])    
        ps[ps == -999] = np.nan

    elif campaign == 'MARCUS':
        if legnum == '1':
            startdate = '2017-10-30'
            enddate = '2017-12-02'
        elif legnum == '2':
            startdate = '2017-12-13'
            enddate = '2018-01-11'
        elif legnum == '3':
            startdate = '2018-01-16'
            enddate = '2018-03-04'
        elif legnum == '4':
            startdate = '2018-03-09'
            enddate = '2018-03-22'
        cday1 = yyyymmdd2cday(startdate, 'noleap')
        cday2 = yyyymmdd2cday(enddate,  'noleap')
        if startdate[0:4] != enddate[0:4]:
            cday2 = cday2 + 365  # cover two years
        time = np.empty(0)
        ps = np.empty(0)
        for cc in range(cday1, cday2 + 1):
            if cc <= 365:
                yyyymmdd = startdate[0:4] + cday2mmdd(cc)
            else:
                yyyymmdd = enddate[0:4] + cday2mmdd(cc-365)
            lst0 = glob.glob(shipmetpath + 'maraadmetX1.b1.' + yyyymmdd + '*')
            (time0, ps0, timeunit, psunit, ps_long_name) = read_met(lst0[0], 'atmospheric_pressure')
            time = np.hstack((time, time0/86400. + cc))
            ps = np.hstack((ps, ps0))
        ps[ps <= -999] = np.nan

    if len(timem) != len(time):
        raise ValueError("model and obs have inconsistent size")
        
    datamask = (ps-psm) > 10
    return(datamask)