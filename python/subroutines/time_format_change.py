# -*- coding: utf-8 -*-

# timestr='16:30:40'

# datestr='2016-05-11'
#%%
def yyyymmdd2cday(datestr,calendar='leap'):
    # must be yyyy-mm-dd format or yyyymmdd format
    if len(datestr)==10:
        dstr=datestr.split('-')
        yr=int(dstr[0])
        mon=int(dstr[1])
        day=int(dstr[2])
    elif len(datestr)==8:
        yr=int(datestr[0:4])
        mon=int(datestr[4:6])
        day=int(datestr[6:8])
    else:
        print('Error: input date format must be yyyymmdd or yyyy-mm-dd')
        return()
        
    
    if calendar=='leap' and yr%4==0:
        mmm=[31,29,31,30,31,30,31,31,30,31,30,31]
    else:
        mmm=[31,28,31,30,31,30,31,31,30,31,30,31]
    
    
    cday=sum(mmm[0:mon-1])+day
    # if mon==1:
    #     cday=day
    # else:
    #     cday=sum(mmm[0:mon-1])+day
        
    return(cday)

#%%
def hhmmss2sec(timestr):
    #must be hh:mm:ss format
    tstr=timestr.split(':')
    hh=int(tstr[0])
    mm=int(tstr[1])
    ss=int(tstr[2])
    
    tsec=hh*3600+mm*60+ss
    
    return(tsec)


#%%
def cday2hhmmss(cday):
    
    time=86400*(cday-int(cday))
    hh=int(time/3600)
    mm=int(time/60-hh*60)
    ss=int(time - hh*3600 - mm*60)
    
    hhmmss=str(hh).zfill(2) +':'+ str(mm).zfill(2) +':'+ str(ss).zfill(2)
    
    return(hhmmss)

#%% 
def cday2mmdd(cday,calendar='noleap'):
    
    if calendar=='leap':
        mmm=[31,29,31,30,31,30,31,31,30,31,30,31]
    else:
        mmm=[31,28,31,30,31,30,31,31,30,31,30,31]
    for ii in range(13):
        if sum(mmm[0:ii])>=int(cday):
            mm = ii
            break
    dd = int(cday)-sum(mmm[0:mm-1])
    
    mmdd = str(mm).zfill(2) +str(dd).zfill(2)
    return(mmdd)

#%%

def timeunit2cday(timeunit,calendar='leap'):
    tstr=timeunit.split(' ')
    cday0 = yyyymmdd2cday(tstr[2],calendar) + hhmmss2sec(tstr[3])/86400
    return(cday0)

# print(yyyymmdd2cday(datestr))
# print(hhmmss2sec(timestr))

#%% change the model time along the aircraft track to aircraft measurement time
def get_obs_time(time_model):
    
    import numpy as np
    time_obs=np.array([time_model[x]-int(time_model[x]) for x in range(len(time_model))])*86400
    t_unique = np.unique(time_model)
    
    for tt in range(len(t_unique)):
        idx = time_model==t_unique[tt]
        time_part = time_obs[idx]
        t_num = len(time_part)
        if tt==0:
            for ii in range(t_num):
                time_part[ii] = time_part[ii] + 60*(ii+30-t_num-15)
        else:
            if round((t_unique[tt]-t_unique[tt-1])*86400)==1800:
                for ii in range(t_num):
                    time_part[ii] = time_part[ii] + 60*(ii-15)
            else:
                for ii in range(t_num):
                    time_part[ii] = time_part[ii] + 60*(ii+30-t_num-15)
        time_obs[idx]=time_part
        
    return(time_obs)