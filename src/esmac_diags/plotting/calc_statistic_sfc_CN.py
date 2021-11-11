"""
# calculate statistics (mean, bias, correlation, RMSE) of Aerosol number concentration
# for surface measurements
# compare models and CPC measurements
"""

# import matplotlib.pyplot as plt
import os
import glob
import numpy as np
import scipy.stats
from ..subroutines.time_format_change import yyyymmdd2cday, cday2mmdd, timeunit2cday
from ..subroutines.read_ARMdata import read_cpc,read_uhsas
from ..subroutines.read_netcdf import read_E3SM
from ..subroutines.quality_control import  qc_remove_neg, qc_mask_qcflag_cpc,qc_mask_qcflag
from ..subroutines.specific_data_treatment import avg_time_1d


def run_plot(settings):
    #%% variables from settings

    campaign = settings['campaign']
    Model_List = settings['Model_List']
    cpcsfcpath = settings['cpcsfcpath']
    cpcusfcpath = settings['cpcusfcpath']
    uhsassfcpath = settings['uhsassfcpath']
    start_date = settings['start_date']
    end_date = settings['end_date']
    E3SM_sfc_path = settings['E3SM_sfc_path']
    figpath_sfc_statistics = settings['figpath_sfc_statistics']

    IOP = settings.get('IOP', None)
    
    #%% other settings
    # set time range you want to average
    # change start date into calendar day
    cday1 = yyyymmdd2cday(start_date,'noleap')
    cday2 = yyyymmdd2cday(end_date,'noleap')
    if start_date[0:4]!=end_date[0:4]:
        raise ValueError('currently not support multiple years. please set start_date and end_date in the same year')
    year0 = start_date[0:4]
    
            
    if not os.path.exists(figpath_sfc_statistics):
        os.makedirs(figpath_sfc_statistics)
    
    missing_value = np.nan
    
    
    #%% read in obs data
    if campaign=='ACEENA':
        # cpc
        if IOP=='IOP1':
            lst = glob.glob(cpcsfcpath+'enaaoscpcfC1.b1.2017062*')+glob.glob(cpcsfcpath+'enaaoscpcfC1.b1.201707*')
        elif IOP=='IOP2':
            lst = glob.glob(cpcsfcpath+'enaaoscpcfC1.b1.201801*')+glob.glob(cpcsfcpath+'enaaoscpcfC1.b1.201802*')
        lst.sort()
        t_cpc=np.empty(0)
        cpc=np.empty(0)
        for filename in lst:
            (time,data,qc,timeunit,cpcunit)=read_cpc(filename)
            data = qc_mask_qcflag(data,qc)
            timestr=timeunit.split(' ')
            date=timestr[2]
            cday=yyyymmdd2cday(date,'noleap')
            # average in time for consistent comparison with model
            time2=np.arange(0,86400,3600)
            data2 = avg_time_1d(np.array(time),np.array(data),time2)
            t_cpc=np.hstack((t_cpc, cday+time2/86400))
            cpc=np.hstack((cpc, data2))
        # fill missing days
        t_cpc2=np.arange(cday1*24,cday2*24+0.01,1)/24.
        cpc2=avg_time_1d(t_cpc,cpc,t_cpc2)
        cpc=cpc2
        t_cpc=t_cpc2
        cpc = qc_remove_neg(cpc)
        # no cpcu
        t_cpcu = np.array(np.nan)
        cpcu = np.array(np.nan)
        
        # uhsas
        if IOP=='IOP1':
            lst = glob.glob(uhsassfcpath+'enaaosuhsasC1.a1.2017062*')+glob.glob(uhsassfcpath+'enaaosuhsasC1.a1.201707*')
        elif IOP=='IOP2':
            lst = glob.glob(uhsassfcpath+'enaaosuhsasC1.a1.201801*')+glob.glob(uhsassfcpath+'enaaosuhsasC1.a1.201802*')
        lst.sort()
        t_uhsas=np.empty(0)
        uhsas=np.empty(0)
        for filename in lst:
            (time,dmin,dmax,data,timeunit,uhsasunit,long_name)=read_uhsas(filename)
            # sum up uhsas data for size >100nm
            data=np.ma.filled(data,np.nan)
            idx100 = dmin>=100
            data1=np.nansum(data[:,idx100],1)
            # average in time for consistent comparison with model
            time2=np.arange(0,86400,3600)
            data2 = avg_time_1d(np.array(time),np.array(data1),time2)
            t_uhsas=np.hstack((t_uhsas, timeunit2cday(timeunit)+time2/86400))
            uhsas=np.hstack((uhsas, data2))
        # fill missing days
        t_uhsas2=np.arange(cday1*24,cday2*24+0.01,1)/24.
        uhsas2=avg_time_1d(t_uhsas,uhsas,t_uhsas2)
        uhsas=uhsas2
        t_uhsas=t_uhsas2
        
        
    elif campaign=='HISCALE':  
        # cpc
        if IOP=='IOP1':
            lst = glob.glob(cpcsfcpath+'sgpaoscpcC1.b1.201604*')+glob.glob(cpcsfcpath+'sgpaoscpcC1.b1.201605*')+glob.glob(cpcsfcpath+'sgpaoscpcC1.b1.201606*')
        elif IOP=='IOP2':
            lst = glob.glob(cpcsfcpath+'sgpaoscpcC1.b1.201608*')+glob.glob(cpcsfcpath+'sgpaoscpcC1.b1.201609*')
        lst.sort()
        t_cpc=np.empty(0)
        cpc=np.empty(0)
        if len(lst)==0:
            t_cpc = np.array(np.nan)
            cpc = np.array(np.nan)
        else:
            for filename in lst:
                (time,data,qc,timeunit,cpcunit)=read_cpc(filename)
                data = qc_mask_qcflag_cpc(data,qc)
                timestr=timeunit.split(' ')
                date=timestr[2]
                cday=yyyymmdd2cday(date,'noleap')
                t_cpc= np.hstack((t_cpc,cday+time/86400))
                cpc=np.hstack((cpc,data))
            # average in time for consistent comparison with model
            t_cpc2=np.arange(cday1*24,cday2*24+0.01,1)/24.
            cpc2=avg_time_1d(t_cpc,cpc,t_cpc2)
            cpc=cpc2
            t_cpc=t_cpc2
            cpc = qc_remove_neg(cpc)
      
        # cpcu
        if IOP=='IOP1':
            lst = glob.glob(cpcusfcpath+'sgpaoscpcuS01.b1.201604*')+glob.glob(cpcusfcpath+'sgpaoscpcuS01.b1.201605*')+glob.glob(cpcusfcpath+'sgpaoscpcuS01.b1.201606*')
        elif IOP=='IOP2':
            lst = glob.glob(cpcusfcpath+'sgpaoscpcuS01.b1.201608*')+glob.glob(cpcusfcpath+'sgpaoscpcuS01.b1.201609*')
        lst.sort()
        t_cpcu=np.empty(0)
        cpcu=np.empty(0)
        if len(lst)==0:
            t_cpcu = np.array(np.nan)
            cpcu = np.array(np.nan)
        else:
            for filename in lst:
                (time,data,qc,timeunit,cpcuunit)=read_cpc(filename)
                data = qc_mask_qcflag_cpc(data,qc)
                timestr=timeunit.split(' ')
                date=timestr[2]
                cday=yyyymmdd2cday(date,'noleap')
                # t_cpcu= np.hstack((t_cpcu,cday+time/86400))
                # cpcu=np.hstack((cpcu,data))
                # average in time for consistent comparison with model
                time2=np.arange(0,86400,3600)
                data2 = avg_time_1d(np.array(time),np.array(data),time2)
                t_cpcu=np.hstack((t_cpcu, cday+time2/86400))
                cpcu=np.hstack((cpcu, data2))
            cpcu = qc_remove_neg(cpcu)
            # # average in time for consistent comparison with model
            # t_cpcu2=np.arange(t_cpcu[0]*24,t_cpcu[-1]*24,1)/24.
            # cpcu2=avg_time_1d(t_cpcu,cpcu,t_cpcu2)
            # cpcu=cpcu2
            # t_cpcu=t_cpcu2
            
        # uhsas
        if IOP=='IOP1':
            lst = glob.glob(uhsassfcpath+'sgpaosuhsasS01.a1.201604*')+glob.glob(uhsassfcpath+'sgpaosuhsasS01.a1.201605*')+glob.glob(uhsassfcpath+'sgpaosuhsasS01.a1.201606*')
        elif IOP=='IOP2':
            lst = glob.glob(uhsassfcpath+'sgpaosuhsasS01.a1.201608*')+glob.glob(uhsassfcpath+'sgpaosuhsasS01.a1.201609*')
        lst.sort()
        t_uhsas=np.empty(0)
        uhsas=np.empty(0)
        for filename in lst:
            (time,dmin,dmax,data,timeunit,uhsasunit,long_name)=read_uhsas(filename)
            # sum up uhsas data for size >100nm
            data=np.ma.filled(data,np.nan)
            idx100 = dmin>=100
            data1=np.nansum(data[:,idx100],1)
            # average in time for consistent comparison with model
            time2=np.arange(0,86400,3600)
            data2 = avg_time_1d(np.array(time),np.array(data1),time2)
            t_uhsas=np.hstack((t_uhsas, timeunit2cday(timeunit)+time2/86400))
            uhsas=np.hstack((uhsas, data2))
        # fill missing days
        t_uhsas2=np.arange(cday1*24,cday2*24+0.01,1)/24.
        uhsas2=avg_time_1d(t_uhsas,uhsas,t_uhsas2)
        uhsas=uhsas2
        t_uhsas=t_uhsas2
        
        
    #%% read in models
    ncn100_m = []
    ncn_m = []
    nucn_m = []
    nmodels = len(Model_List)
    for mm in range(nmodels):
        tmp_ncn100=np.empty(0)
        tmp_ncn=np.empty(0)
        tmp_nucn=np.empty(0)
        timem=np.empty(0)
        for cday in range(cday1,cday2+1):
            mmdd=cday2mmdd(cday)
            date=year0+'-'+mmdd[0:2]+'-'+mmdd[2:4]
            
            filename_input = E3SM_sfc_path+'SFC_CNsize_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
            (time,ncn,timemunit,dataunit,long_name)=read_E3SM(filename_input,'NCN')
            (time,nucn,timemunit,dataunit,long_name)=read_E3SM(filename_input,'NUCN')
            (time,ncnall,timemunit,dataunit,long_name)=read_E3SM(filename_input,'NCNall')
            
            timem = np.hstack((timem,time))
            tmp_ncn = np.hstack((tmp_ncn,ncn*1e-6))
            tmp_nucn = np.hstack((tmp_nucn,nucn*1e-6))
            tmp_ncn100 = np.hstack((tmp_ncn100, np.sum(ncnall[100:,:],0)*1e-6))  
        
        ncn100_m.append(tmp_ncn100)
        ncn_m.append(tmp_ncn)
        nucn_m.append(tmp_nucn)
        
    #%% calculate statistics
    
    # only choose the prescribed time range
    idx = np.logical_and(t_cpc>=cday1, t_cpc<=cday2)
    cpc=cpc[idx]
    t_cpc=t_cpc[idx]
    idx = np.logical_and(t_cpcu>=cday1, t_cpcu<=cday2)
    cpcu=cpcu[idx]
    t_cpcu=t_cpcu[idx]
    idx = np.logical_and(t_uhsas>=cday1, t_uhsas<=cday2)
    uhsas=uhsas[idx]
    t_uhsas=t_uhsas[idx]
    idx = np.logical_and(timem>=cday1, timem<=cday2)
    for mm in range(nmodels):
        ncn100_m[mm]=ncn100_m[mm][idx]
        ncn_m[mm]=ncn_m[mm][idx]
        nucn_m[mm]=nucn_m[mm][idx]
    timem=timem[idx]
    
    
    
    # select only valid data in obs and the corresponding data in models
    idx100 = ~np.isnan(uhsas)
    idx10 = ~np.isnan(cpc)
    idx3 = ~np.isnan(cpcu)
    
    mean100 = [None]*(nmodels+1)
    mean10 = [None]*(nmodels+1)
    mean3 = [None]*(nmodels+1)
    std100 = [None]*(nmodels+1)
    std10 = [None]*(nmodels+1)
    std3 = [None]*(nmodels+1)
    bias100 = [None]*(nmodels)
    bias10 = [None]*(nmodels)
    bias3 = [None]*(nmodels)
    corr100 = [None]*(nmodels)
    corr10 = [None]*(nmodels)
    corr3 = [None]*(nmodels)
    rmse100 = [None]*(nmodels)
    rmse10 = [None]*(nmodels)
    rmse3 = [None]*(nmodels)
    p10_100 = [None]*(nmodels+1)
    p10_10 = [None]*(nmodels+1)
    p10_3 = [None]*(nmodels+1)
    p25_100 = [None]*(nmodels+1)
    p25_10 = [None]*(nmodels+1)
    p25_3 = [None]*(nmodels+1)
    p75_100 = [None]*(nmodels+1)
    p75_10 = [None]*(nmodels+1)
    p75_3 = [None]*(nmodels+1)
    p90_100 = [None]*(nmodels+1)
    p90_10 = [None]*(nmodels+1)
    p90_3 = [None]*(nmodels+1)
        
    if len(idx100)==0 or sum(idx100)/len(idx100)<0.1:   # two few observation available
        # for obs
        mean100[nmodels] = missing_value
        std100[nmodels] = missing_value
        p10_100[nmodels] = missing_value
        p25_100[nmodels] = missing_value
        p75_100[nmodels] = missing_value
        p90_100[nmodels] = missing_value
        # for models
        for mm in range(nmodels):
            mean100[mm] = np.nanmean(ncn100_m[mm][idx100])
            std100[mm] = np.nanstd(ncn100_m[mm][idx100])
            p10_100[mm] = np.nanpercentile(ncn100_m[mm][idx100],10)
            p25_100[mm] = np.nanpercentile(ncn100_m[mm][idx100],25)
            p75_100[mm] = np.nanpercentile(ncn100_m[mm][idx100],75)
            p90_100[mm] = np.nanpercentile(ncn100_m[mm][idx100],90)
            bias100[mm] =  missing_value
            corr100[mm] = [missing_value, missing_value]
            rmse100[mm] = missing_value
    else:
        # for obs
        mean100[nmodels] = np.nanmean(uhsas[idx100])
        std100[nmodels] = np.nanstd(uhsas[idx100])
        p10_100[nmodels] = np.nanpercentile(uhsas[idx100],10)
        p25_100[nmodels] = np.nanpercentile(uhsas[idx100],25)
        p75_100[nmodels] = np.nanpercentile(uhsas[idx100],75)
        p90_100[nmodels] = np.nanpercentile(uhsas[idx100],90)
        # for models
        for mm in range(nmodels):
            mean100[mm] = np.nanmean(ncn100_m[mm][idx100])
            std100[mm] = np.nanstd(ncn100_m[mm][idx100])
            p10_100[mm] = np.nanpercentile(ncn100_m[mm][idx100],10)
            p25_100[mm] = np.nanpercentile(ncn100_m[mm][idx100],25)
            p75_100[mm] = np.nanpercentile(ncn100_m[mm][idx100],75)
            p90_100[mm] = np.nanpercentile(ncn100_m[mm][idx100],90)
            bias100[mm] = mean100[mm] - mean100[nmodels]
            c100 = scipy.stats.pearsonr(ncn100_m[mm][idx100],uhsas[idx100])
            corr100[mm] = [c100[0],c100[1]]
            rmse100[mm] = np.sqrt(((ncn100_m[mm][idx100]-uhsas[idx100])**2).mean())
    
    if len(idx10)==0 or sum(idx10)/len(idx10)<0.1:   # two few observation available
        # for obs
        mean10[nmodels] = missing_value
        std10[nmodels] = missing_value
        p10_10[nmodels] = missing_value
        p25_10[nmodels] = missing_value
        p75_10[nmodels] = missing_value
        p90_10[nmodels] = missing_value
        # for models
        for mm in range(nmodels):
            mean10[mm] = np.nanmean(ncn_m[mm][idx10])
            std10[mm] = np.nanstd(ncn_m[mm][idx10])
            p10_10[mm] = np.nanpercentile(ncn_m[mm][idx10],10)
            p25_10[mm] = np.nanpercentile(ncn_m[mm][idx10],25)
            p75_10[mm] = np.nanpercentile(ncn_m[mm][idx10],75)
            p90_10[mm] = np.nanpercentile(ncn_m[mm][idx10],90)
            bias10[mm] =  missing_value
            corr10[mm] = [missing_value, missing_value]
            rmse10[mm] = missing_value
    else:
        # for obs
        mean10[nmodels] = np.nanmean(cpc[idx10])
        std10[nmodels] = np.nanstd(cpc[idx10])
        p10_10[nmodels] = np.nanpercentile(cpc[idx10],10)
        p25_10[nmodels] = np.nanpercentile(cpc[idx10],25)
        p75_10[nmodels] = np.nanpercentile(cpc[idx10],75)
        p90_10[nmodels] = np.nanpercentile(cpc[idx10],90)
        # for models
        for mm in range(nmodels):
            mean10[mm] = np.nanmean(ncn_m[mm][idx10])
            std10[mm] = np.nanstd(ncn_m[mm][idx10])
            p10_10[mm] = np.nanpercentile(ncn_m[mm][idx10],10)
            p25_10[mm] = np.nanpercentile(ncn_m[mm][idx10],25)
            p75_10[mm] = np.nanpercentile(ncn_m[mm][idx10],75)
            p90_10[mm] = np.nanpercentile(ncn_m[mm][idx10],90)
            bias10[mm] = mean10[mm] - mean10[nmodels]
            c10 = scipy.stats.pearsonr(ncn_m[mm][idx10],cpc[idx10])
            corr10[mm] = [c10[0],c10[1]]
            rmse10[mm] = np.sqrt(((ncn_m[mm][idx10]-cpc[idx10])**2).mean())
    
    if len(idx3)==0 or sum(idx3)/len(idx3)<0.1:   # two few observation available
        # for obs
        mean3[nmodels] = missing_value
        std3[nmodels] = missing_value
        p10_3[nmodels] = missing_value
        p25_3[nmodels] = missing_value
        p75_3[nmodels] = missing_value
        p90_3[nmodels] = missing_value
        # for models
        for mm in range(nmodels):
            mean3[mm] = np.nanmean(nucn_m[mm][idx3])
            std3[mm] = np.nanstd(nucn_m[mm][idx3])
            p10_3[mm] = np.nanpercentile(nucn_m[mm][idx3],10)
            p25_3[mm] = np.nanpercentile(nucn_m[mm][idx3],25)
            p75_3[mm] = np.nanpercentile(nucn_m[mm][idx3],75)
            p90_3[mm] = np.nanpercentile(nucn_m[mm][idx3],90)
            bias3[mm] =  missing_value
            corr3[mm] = [missing_value, missing_value]
            rmse3[mm] = missing_value
    else:
        # for obs
        mean3[nmodels] = np.nanmean(cpcu[idx3])
        std3[nmodels] = np.nanstd(cpcu[idx3])
        p10_3[nmodels] = np.nanpercentile(cpcu[idx3],10)
        p25_3[nmodels] = np.nanpercentile(cpcu[idx3],25)
        p75_3[nmodels] = np.nanpercentile(cpcu[idx3],75)
        p90_3[nmodels] = np.nanpercentile(cpcu[idx3],90)
        # for models
        for mm in range(nmodels):
            mean3[mm] = np.nanmean(nucn_m[mm][idx3])
            std3[mm] = np.nanstd(nucn_m[mm][idx3])
            p10_3[mm] = np.nanpercentile(nucn_m[mm][idx3],10)
            p25_3[mm] = np.nanpercentile(nucn_m[mm][idx3],25)
            p75_3[mm] = np.nanpercentile(nucn_m[mm][idx3],75)
            p90_3[mm] = np.nanpercentile(nucn_m[mm][idx3],90)
            bias3[mm] = mean3[mm] - mean3[nmodels]
            c3 = scipy.stats.pearsonr(nucn_m[mm][idx3],cpcu[idx3])
            corr3[mm] = [c3[0],c3[1]]
            rmse3[mm] = np.sqrt(((nucn_m[mm][idx3]-cpcu[idx3])**2).mean())
    
    
    #%% write out files
    
    outfile = figpath_sfc_statistics+'statistics_CN10nm_'+campaign+'_'+IOP+'.txt'
    print('write statistics to file '+outfile)
    
    with open(outfile, 'w') as f:
        f.write('statistics of Aerosol Number Concentration comparing with CPC(>10nm). sample size '+format(sum(idx10))+'\n')
        line1 = list(Model_List)
        line1.insert(0,' --- ')
        line1.append('OBS')
        for ii in range(len(line1)):
            f.write(format(line1[ii],'10s')+', ')
        # write mean
        f.write('\n mean,\t')
        for ii in range(len(mean10)):
            f.write(format(mean10[ii],'10.2f')+', ')
        # write std
        f.write('\n std. dev.,')
        for ii in range(len(std10)):
            f.write(format(std10[ii],'10.2f')+', ')
        # write percentiles
        f.write('\n 10% percentile: ')
        for ii in range(len(p10_10)):
            f.write(format(p10_10[ii],'10.2f')+', ')
        f.write('\n 25% percentile: ')
        for ii in range(len(p25_10)):
            f.write(format(p25_10[ii],'10.2f')+', ')
        f.write('\n 75% percentile: ')
        for ii in range(len(p75_10)):
            f.write(format(p75_10[ii],'10.2f')+', ')
        f.write('\n 90% percentile: ')
        for ii in range(len(p90_10)):
            f.write(format(p90_10[ii],'10.2f')+', ')
        # write bias
        f.write('\n bias,\t')
        for ii in range(len(bias10)):
            f.write(format(bias10[ii],'10.2f')+', ')
        # write rmse
        f.write('\n RMSE,\t')
        for ii in range(len(rmse10)):
            f.write(format(rmse10[ii],'10.2f')+', ')
        # write correlation
        f.write('\n corrcoef,\t')
        for ii in range(len(rmse10)):
            f.write(format(corr10[ii][0],'10.4f')+', ')
        # write p value of correlation
        f.write('\n P_corr,\t')
        for ii in range(len(rmse10)):
            f.write(format(corr10[ii][1],'10.2f')+', ')
            
            
    outfile = figpath_sfc_statistics+'statistics_CN3nm_'+campaign+'_'+IOP+'.txt'
    print('write statistics to file '+outfile)
    
    with open(outfile, 'w') as f:
        f.write('statistics of Aerosol Number Concentration comparing with CPC(>3nm). sample size '+format(sum(idx3))+'\n')
        line1 = list(Model_List)
        line1.insert(0,' --- ')
        line1.append('OBS')
        for ii in range(len(line1)):
            f.write(format(line1[ii],'10s')+', ')
        # write mean
        f.write('\n mean,\t')
        for ii in range(len(mean3)):
            f.write(format(mean3[ii],'10.2f')+', ')
        # write std
        f.write('\n std. dev.,')
        for ii in range(len(std3)):
            f.write(format(std3[ii],'10.2f')+', ')
        # write percentiles
        f.write('\n 10% percentile: ')
        for ii in range(len(p10_3)):
            f.write(format(p10_3[ii],'10.2f')+', ')
        f.write('\n 25% percentile: ')
        for ii in range(len(p25_3)):
            f.write(format(p25_3[ii],'10.2f')+', ')
        f.write('\n 75% percentile: ')
        for ii in range(len(p75_3)):
            f.write(format(p75_3[ii],'10.2f')+', ')
        f.write('\n 90% percentile: ')
        for ii in range(len(p90_3)):
            f.write(format(p90_3[ii],'10.2f')+', ')
        # write bias
        f.write('\n bias,\t')
        for ii in range(len(bias3)):
            f.write(format(bias3[ii],'10.2f')+', ')
        # write rmse
        f.write('\n RMSE,\t')
        for ii in range(len(rmse3)):
            f.write(format(rmse3[ii],'10.2f')+', ')
        # write correlation
        f.write('\n corrcoef,\t')
        for ii in range(len(rmse3)):
            f.write(format(corr3[ii][0],'10.4f')+', ')
        # write p value of correlation
        f.write('\n P_corr,\t')
        for ii in range(len(rmse3)):
            f.write(format(corr3[ii][1],'10.2f')+', ')
        
    
    outfile = figpath_sfc_statistics+'statistics_CN100nm_'+campaign+'_'+IOP+'.txt'
    print('write statistics to file '+outfile)
    
    with open(outfile, 'w') as f:
        f.write('statistics of Aerosol Number Concentration comparing with UHSAS (>100nm). sample size '+format(sum(idx100))+'\n')
        line1 = list(Model_List)
        line1.insert(0,' --- ')
        line1.append('OBS')
        for ii in range(len(line1)):
            f.write(format(line1[ii],'10s')+', ')
        # write mean
        f.write('\n mean,\t')
        for ii in range(len(mean100)):
            f.write(format(mean100[ii],'10.2f')+', ')
        # write std
        f.write('\n std. dev.,')
        for ii in range(len(std100)):
            f.write(format(std100[ii],'10.2f')+', ')
        # write percentiles
        f.write('\n 10% percentile: ')
        for ii in range(len(p10_100)):
            f.write(format(p10_100[ii],'10.2f')+', ')
        f.write('\n 25% percentile: ')
        for ii in range(len(p25_100)):
            f.write(format(p25_100[ii],'10.2f')+', ')
        f.write('\n 75% percentile: ')
        for ii in range(len(p75_100)):
            f.write(format(p75_100[ii],'10.2f')+', ')
        f.write('\n 90% percentile: ')
        for ii in range(len(p90_100)):
            f.write(format(p90_100[ii],'10.2f')+', ')
        # write bias
        f.write('\n bias,\t')
        for ii in range(len(bias100)):
            f.write(format(bias100[ii],'10.2f')+', ')
        # write rmse
        f.write('\n RMSE,\t')
        for ii in range(len(rmse100)):
            f.write(format(rmse100[ii],'10.2f')+', ')
        # write correlation
        f.write('\n corrcoef,\t')
        for ii in range(len(rmse100)):
            f.write(format(corr100[ii][0],'10.4f')+', ')
        # write p value of correlation
        f.write('\n P_corr,\t')
        for ii in range(len(rmse100)):
            f.write(format(corr100[ii][1],'10.2f')+', ')
