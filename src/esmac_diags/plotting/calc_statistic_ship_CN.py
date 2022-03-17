"""
# calculate statistics (mean, bias, correlation, RMSE) of Aerosol number concentration
# for ship measurements
# compare models and CPC/UHSAS measurements
"""

import os
import glob
import numpy as np
import scipy.stats
from ..subroutines.read_ARMdata import read_cpc, read_uhsas
from ..subroutines.read_netcdf import read_E3SM, read_ship_exhaustfree
from ..subroutines.time_format_change import  cday2mmdd
from ..subroutines.specific_data_treatment import mask_model_ps, avg_time_1d
from ..subroutines.quality_control import qc_mask_qcflag, qc_remove_neg,qc_cn_max


def run_plot(settings):
    #%% variables from settings

    campaign = settings['campaign']
    Model_List = settings['Model_List']
    shipcpcpath = settings['shipcpcpath']
    shipmetpath = settings['shipmetpath']
    shipuhsaspath = settings['shipuhsaspath']
    shipcn_exhaustfree_path = settings['shipcn_exhaustfree_path']
    E3SM_ship_path = settings['E3SM_ship_path']
    figpath_ship_statistics = settings['figpath_ship_statistics']

    #%% other settings
    if not os.path.exists(figpath_ship_statistics):
        os.makedirs(figpath_ship_statistics)
    missing_value = np.nan
    
    #%% read model
    lst = glob.glob(E3SM_ship_path+'Ship_CNsize_'+campaign+'_'+Model_List[0]+'_shipleg*.nc')
    lst.sort()
    
    nmodels = len(Model_List)
    ncn10all = []
    ncn100all = []
    timemall = np.empty(0)
    for mm in range(nmodels):
        ncn10all.append(np.empty(0))
        ncn100all.append(np.empty(0))
        
    for ll in range(len(lst)):
        
        if campaign=='MAGIC':
            legnum=lst[ll][-5:-3]
        elif campaign=='MARCUS':
            legnum=lst[ll][-4]
        else: 
            raise ValueError('please check campaign name: '+campaign)
        print('legnum '+format(legnum))
        

        for mm in range(nmodels):
            filenamem = E3SM_ship_path+'Ship_CNsize_'+campaign+'_'+Model_List[mm]+'_shipleg'+legnum+'.nc'
        
            (timem,NCNall,timeunitm,datamunit,datamlongname)=read_E3SM(filenamem,'NCNall')
            (timem,data,timeunitm,datamunit,datamlongname)=read_E3SM(filenamem,'NCN')
            
            # mask data where model grid is not at ocean surface (Ps is too different than obs)
            filenamem = E3SM_ship_path+'Ship_vars_'+campaign+'_'+Model_List[mm]+'_shipleg'+legnum+'.nc'
            (timem,psm,timeunitx,psmunit,psmlongname)=read_E3SM(filenamem,'PS')
            datamask = mask_model_ps(timem,0.01*psm,legnum,campaign,shipmetpath)
            
            datam = data*1e-6    # change unit from 1/m3 to 1/cm3
            datam[datamask]=np.nan
        
            # Calculate model aerosol number concentration for UHSAS size range
            datam2=np.nansum(NCNall[100:1000,:],0) * 1e-6    # change unit from 1/m3 to 1/cm3
            datam2[datamask]=np.nan

            ncn10all[mm] = np.hstack((ncn10all[mm],datam))
            ncn100all[mm] = np.hstack((ncn100all[mm],datam2))
            
            if campaign=='MAGIC':
                if int(legnum)<=9:
                    timemall = np.hstack((timemall, timem))
                else:
                    timemall = np.hstack((timemall, timem+365))
            elif campaign=='MARCUS':
                if int(legnum)<=2:
                    timemall = np.hstack((timemall, timem))
                else:
                    timemall = np.hstack((timemall, timem+365))
                
                
    #%% read in observation
    if campaign=='MAGIC':
        t_cpc = np.empty(0)
        t_uhsas = np.empty(0)
        cpcall = np.empty(0)
        uhsasall = np.empty(0)
            
        for ll in range(len(lst)):       
            # use lat/lon from extracted model data
            filenamem = lst[ll]
            legnum=lst[ll][-5:-3]
            (timem,[lat1,lon1],timeunitm,varmunit,varmlongname)=read_E3SM(filenamem,['lat','lon'])
                
            print('legnum '+format(legnum))
            # find the days related to the ship leg
            day = [int(a) for a in timem]
            day = list(set(day))
            day.sort()
            
            # read in CPC    
            for dd in day:
                if int(legnum)<=9:
                    if dd<=365:  # year 2012
                        filenameo = glob.glob(shipcpcpath+'magaoscpcfM1.a1.2012'+cday2mmdd(dd,calendar='noleap')+'.*')
                    else:
                        filenameo = glob.glob(shipcpcpath+'magaoscpcfM1.a1.2013'+cday2mmdd(dd-365,calendar='noleap')+'.*')
                    cday = dd
                else:
                    filenameo = glob.glob(shipcpcpath+'magaoscpcfM1.a1.2013'+cday2mmdd(dd,calendar='noleap')+'.*')
                    cday = dd+365
                if len(filenameo)==0:
                    continue  # some days may be missing
                (time,obs,qc,timeunit,dataunit)=read_cpc(filenameo[0])
                obs=np.ma.filled(obs)
                obs=qc_mask_qcflag(obs,qc)
                obs=qc_cn_max(obs,10)
                obs = qc_remove_neg(obs)
                
                time2 = np.arange(150., 86400., 300.)
                obs2 = avg_time_1d(time, obs, time2)
                t_cpc=np.hstack((t_cpc, cday+time2/86400))
                cpcall=np.hstack((cpcall, obs2))

            
            # read in UHSAS
            for dd in day:
                if int(legnum)<=9:
                    if dd<=365:  # year 2012
                        filenameo = glob.glob(shipuhsaspath+'magaosuhsasM1.a1.2012'+cday2mmdd(dd,calendar='noleap')+'.*.cdf')
                    else:
                        filenameo = glob.glob(shipuhsaspath+'magaosuhsasM1.a1.2013'+cday2mmdd(dd-365,calendar='noleap')+'.*.cdf')
                    cday = dd
                else:
                    filenameo = glob.glob(shipuhsaspath+'magaosuhsasM1.a1.2013'+cday2mmdd(dd,calendar='noleap')+'.*.cdf')
                    cday = dd+365

                if len(filenameo)==0:
                    continue  # some days may be missing
                if len(filenameo)>1:
                    raise ValueError('find too many files')
                    
                (time,dmin,dmax,obs,timeunit,uhunit,uhlongname)=read_uhsas(filenameo[0])
                obs=np.ma.filled(obs)
                obs=qc_remove_neg(obs)
                obs = np.nansum(obs[:, dmin>=100],1)
                obs=qc_cn_max(obs,100)
                time2 = np.arange(150., 86400., 300.)
                obs2 = avg_time_1d(time,obs,time2)
                t_uhsas = np.hstack((t_uhsas,time2/86400+cday))
                uhsasall = np.hstack((uhsasall, obs2))
                
                
    elif campaign=='MARCUS':
        (t_cpc, cpc, timeunit, cpcunit, cpc_longname) = read_ship_exhaustfree(shipcn_exhaustfree_path + \
                                                     'CPC_UHSAS_exhaustfree_1hr.nc', 'CPC')
        (t_uhsas, uhsas, timeunit, uhsasunit, uhsas_longname) = read_ship_exhaustfree(shipcn_exhaustfree_path + \
                                                     'CPC_UHSAS_exhaustfree_1hr.nc', 'UHSAS100')
        cpcall = qc_remove_neg(cpc)
        uhsasall = qc_remove_neg(uhsas)
            
        idx = np.logical_and(t_cpc>=timemall[0]-0.04, t_cpc<=timemall[-1]+0.04)
        t_cpc = t_cpc[idx]
        t_uhsas = t_uhsas[idx]
        cpcall = cpcall[idx]
        uhsasall = uhsasall[idx]
             
    #%% calculate statistics
    
    time0 = np.arange(timemall[0],timemall[-1],1./24)
    cpcall = avg_time_1d(t_cpc,cpcall,time0)
    uhsasall = avg_time_1d(t_uhsas,uhsasall,time0)
    for mm in range(nmodels):
        ncn10all[mm] = avg_time_1d(timemall,ncn10all[mm],time0)
        ncn100all[mm] = avg_time_1d(timemall,ncn100all[mm],time0)
        
    
    # select only valid data in obs and the corresponding data in models (all data are not NAN)
    idx10 = sum(np.vstack((~np.isnan(ncn10all),~np.isnan(cpcall))))==nmodels+1
    idx100 = sum(np.vstack((~np.isnan(ncn100all),~np.isnan(uhsasall))))==nmodels+1
    
    mean10 = [None]*(nmodels+1)
    mean100 = [None]*(nmodels+1)
    std10 = [None]*(nmodels+1)
    std100 = [None]*(nmodels+1)
    bias10 = [None]*(nmodels)
    bias100 = [None]*(nmodels)
    corr10 = [None]*(nmodels)
    corr100 = [None]*(nmodels)
    rmse10 = [None]*(nmodels)
    rmse100 = [None]*(nmodels)
    p10_100 = [None]*(nmodels+1)
    p10_10 = [None]*(nmodels+1)
    p25_100 = [None]*(nmodels+1)
    p25_10 = [None]*(nmodels+1)
    p75_100 = [None]*(nmodels+1)
    p75_10 = [None]*(nmodels+1)
    p90_100 = [None]*(nmodels+1)
    p90_10 = [None]*(nmodels+1)
        
    
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
            mean10[mm] = np.nanmean(ncn10all[mm][idx10])
            std10[mm] = np.nanstd(ncn10all[mm][idx10])
            p10_10[mm] = np.nanpercentile(ncn10all[mm][idx10],10)
            p25_10[mm] = np.nanpercentile(ncn10all[mm][idx10],25)
            p75_10[mm] = np.nanpercentile(ncn10all[mm][idx10],75)
            p90_10[mm] = np.nanpercentile(ncn10all[mm][idx10],90)
            bias10[mm] =  missing_value
            corr10[mm] = [missing_value, missing_value]
            rmse10[mm] = missing_value
    else:
        # for obs
        mean10[nmodels] = np.nanmean(cpcall[idx10])
        std10[nmodels] = np.nanstd(cpcall[idx10])
        p10_10[nmodels] = np.nanpercentile(cpcall[idx10],10)
        p25_10[nmodels] = np.nanpercentile(cpcall[idx10],25)
        p75_10[nmodels] = np.nanpercentile(cpcall[idx10],75)
        p90_10[nmodels] = np.nanpercentile(cpcall[idx10],90)
        # for models
        for mm in range(nmodels):
            mean10[mm] = np.nanmean(ncn10all[mm][idx10])
            std10[mm] = np.nanstd(ncn10all[mm][idx10])
            p10_10[mm] = np.nanpercentile(ncn10all[mm][idx10],10)
            p25_10[mm] = np.nanpercentile(ncn10all[mm][idx10],25)
            p75_10[mm] = np.nanpercentile(ncn10all[mm][idx10],75)
            p90_10[mm] = np.nanpercentile(ncn10all[mm][idx10],90)
            bias10[mm] = mean10[mm] - mean10[nmodels]
            c10 = scipy.stats.pearsonr(ncn10all[mm][idx10],cpcall[idx10])
            corr10[mm] = [c10[0],c10[1]]
            rmse10[mm] = np.sqrt(((ncn10all[mm][idx10]-cpcall[idx10])**2).mean())
    
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
            mean100[mm] = np.nanmean(ncn100all[mm][idx100])
            std100[mm] = np.nanstd(ncn100all[mm][idx100])
            p10_100[mm] = np.nanpercentile(ncn100all[mm][idx100],10)
            p25_100[mm] = np.nanpercentile(ncn100all[mm][idx100],25)
            p75_100[mm] = np.nanpercentile(ncn100all[mm][idx100],75)
            p90_100[mm] = np.nanpercentile(ncn100all[mm][idx100],90)
            bias100[mm] =  missing_value
            corr100[mm] = [missing_value, missing_value]
            rmse100[mm] = missing_value
    else:
        # for obs
        mean100[nmodels] = np.nanmean(uhsasall[idx100])
        std100[nmodels] = np.nanstd(uhsasall[idx100])
        p10_100[nmodels] = np.nanpercentile(uhsasall[idx100],10)
        p25_100[nmodels] = np.nanpercentile(uhsasall[idx100],25)
        p75_100[nmodels] = np.nanpercentile(uhsasall[idx100],75)
        p90_100[nmodels] = np.nanpercentile(uhsasall[idx100],90)
        # for models
        for mm in range(nmodels):
            mean100[mm] = np.nanmean(ncn100all[mm][idx100])
            std100[mm] = np.nanstd(ncn100all[mm][idx100])
            p10_100[mm] = np.nanpercentile(ncn100all[mm][idx100],10)
            p25_100[mm] = np.nanpercentile(ncn100all[mm][idx100],25)
            p75_100[mm] = np.nanpercentile(ncn100all[mm][idx100],75)
            p90_100[mm] = np.nanpercentile(ncn100all[mm][idx100],90)
            bias100[mm] = mean100[mm] - mean100[nmodels]
            c100 = scipy.stats.pearsonr(ncn100all[mm][idx100],uhsasall[idx100])
            corr100[mm] = [c100[0],c100[1]]
            rmse100[mm] = np.sqrt(((ncn100all[mm][idx100]-uhsasall[idx100])**2).mean())
    
    
    #%% write out files
    
    outfile = figpath_ship_statistics+'statistics_CN10nm_'+campaign+'.txt'
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
            
            
    outfile = figpath_ship_statistics+'statistics_CN100nm_'+campaign+'.txt'
    print('write statistics to file '+outfile)
    
    with open(outfile, 'w') as f:
        f.write('statistics of Aerosol Number Concentration comparing with UHSAS100(>100nm). sample size '+format(sum(idx100))+'\n')
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

