"""
# calculate statistics (mean, bias, correlation, RMSE) of Aerosol number concentration
# for aircraft measurements
# compare models and CPC measurements
"""

import os
import glob
import numpy as np
import scipy.stats
from ..subroutines.read_aircraft import read_cpc, read_RF_NCAR
from ..subroutines.read_netcdf import read_merged_size,read_extractflight
from ..subroutines.quality_control import qc_cpc_air, qc_remove_neg, qc_mask_takeoff_landing

def run_plot(settings):
    #%% variables from settings

    campaign = settings['campaign']
    Model_List = settings['Model_List']
    E3SM_aircraft_path = settings['E3SM_aircraft_path']
    figpath_aircraft_statistics = settings['figpath_aircraft_statistics']

    IOP = settings.get('IOP', None)
    cpcpath = settings.get('cpcpath', None)
    merged_size_path = settings.get('merged_size_path', None)
    RFpath = settings.get('RFpath', None)

    #%% other settings
    if not os.path.exists(figpath_aircraft_statistics):
        os.makedirs(figpath_aircraft_statistics)
       
    missing_value = np.nan
    
    #%% find files for flight information
    
    lst = glob.glob(E3SM_aircraft_path+'Aircraft_vars_'+campaign+'_'+Model_List[0]+'_*.nc')
    lst.sort()
    if len(lst)==0:
        raise ValueError('cannot find any file')
    # choose files for specific IOP
    if campaign=='HISCALE':
        if IOP=='IOP1':
            lst=lst[0:17]
        elif IOP=='IOP2':
            lst=lst[17:]
        elif IOP[0:4]=='2016':
            a=lst[0].split('_'+Model_List[0]+'_')
            lst = glob.glob(a[0]+'_'+Model_List[0]+'_'+IOP+'*')
            lst.sort()
    elif campaign=='ACEENA':
        if IOP=='IOP1':
            lst=lst[0:20]
        elif IOP=='IOP2':
            lst=lst[20:]
        elif IOP[0:4]=='2017' or IOP[0:4]=='2018':
            a=lst[0].split('_'+Model_List[0]+'_')
            lst = glob.glob(a[0]+'_'+Model_List[0]+'_'+IOP+'*')
            lst.sort()
            
    alldates = [x.split('_')[-1].split('.')[0] for x in lst]
        
    #%% read all data
    
    uhsas100_o = np.empty(0)    # large particles. UHSAS for CSET and SOCRATES, PCASP for ACEENA and HISCALE
    cpc10_o = np.empty(0)
    cpc3_o = np.empty(0)
    ncn100_m = []
    ncn10_m = []
    ncn3_m = []
    nmodels=len(Model_List)
    for mm in range(nmodels):
        ncn100_m.append(np.empty(0))
        ncn10_m.append(np.empty(0))
        ncn3_m.append(np.empty(0))
        
    print('reading '+format(len(alldates))+' files to calculate the statistics: ')
    
    for date in alldates:
        print(date)
        
        #%% read in Models
        for mm in range(nmodels):
            filename_m = E3SM_aircraft_path+'Aircraft_CNsize_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
        
            (timem,heightm,cpc_m,timeunitm,ncn_unit,ncn_longname)=read_extractflight(filename_m,'NCN')
            (timem,heightm,cpcu_m,timeunitm,ncnu_unit,ncnu_longname)=read_extractflight(filename_m,'NUCN')
            (timem,heightm,ncnall,timeunitm,ncnall_unit,ncnall_longname)=read_extractflight(filename_m,'NCNall')
            
            if campaign=='HISCALE':
                if IOP=='IOP1':  # PCASP for HISCALE IOP1 size from 0.12 to 3 um
                    ncn100_m[mm] = np.hstack((ncn100_m[mm], np.sum(ncnall[120:,:],0)*1e-6))
                elif IOP=='IOP2': # PCASP for HISCALE IOP1 size from 0.09 to 3 um
                    ncn100_m[mm] = np.hstack((ncn100_m[mm], np.sum(ncnall[90:,:],0)*1e-6))
            else:
                ncn100_m[mm] = np.hstack((ncn100_m[mm], np.sum(ncnall[100:,:],0)*1e-6))  
            ncn10_m[mm] = np.hstack((ncn10_m[mm], cpc_m*1e-6))    # #/m3 to #/cm3
            ncn3_m[mm] = np.hstack((ncn3_m[mm], cpcu_m*1e-6))    # #/m3 to #/cm3
            
            
        #%% read in flight measurements (CPC and PCASP) for HISCALE and ACEENA
        if campaign in ['HISCALE', 'ACEENA']:
            if date[-1]=='a':
                flightidx=1
            else:
                flightidx=2
            if campaign=='HISCALE':
                filename_c=glob.glob(cpcpath+'CPC_G1_'+date[0:8]+'*R2_HiScale001s.ict.txt')
                filename_merge = merged_size_path+'merged_bin_fims_pcasp_HISCALE_'+date+'.nc'
            elif campaign=='ACEENA':
                filename_c=glob.glob(cpcpath+'CPC_G1_'+date[0:8]+'*R2_ACEENA001s.ict')    
                filename_merge = merged_size_path+'merged_bin_fims_pcasp_opc_ACEENA_'+date+'.nc'
            filename_c.sort()
            
            # read in CPC
            if len(filename_c)==1 or len(filename_c)==2: # some days have two flights
                (cpc,cpclist)=read_cpc(filename_c[flightidx-1])
                # fill missing timestep
                if np.logical_and(campaign=='ACEENA', date=='20180216a'):
                    cpc=np.insert(cpc,1404,(cpc[:,1403]+cpc[:,1404])/2,axis=1) 
                elif np.logical_and(campaign=='HISCALE', date=='20160425a'):
                    cpc=np.insert(cpc,0,cpc[:,0],axis=1)
                    cpc[0,0]=cpc[0,0]-1
                time_cpc = cpc[0,:]
                cpc10 = cpc[1,:]
                cpc3 = cpc[2,:]
            elif len(filename_c)==0:
                time_cpc=timem
                cpc10=np.nan*np.empty([len(timem)])
                cpc3=np.nan*np.empty([len(timem)])
            else:
                raise ValueError('find too many files: '+filename_c)
            
            # some quality checks
            (cpc3,cpc10) = qc_cpc_air(cpc3, cpc10)
            
            # read in PCASP
            (time_merge,size,pcasp,timeunit,pcaspunit,pcasplongname)=read_merged_size(filename_merge,'totalnum_pcasp')
            pcasp=qc_remove_neg(pcasp)
            if len(time_merge)!=len(time_cpc):
                raise ValueError('time dimension is inconsistent ')
            
            # exclude 30min after takeoff and before landing
            cpc3 = qc_mask_takeoff_landing(time_cpc,cpc3)
            cpc10 = qc_mask_takeoff_landing(time_cpc,cpc10)
            pcasp = qc_mask_takeoff_landing(time_cpc,pcasp)
            
            cpc10_o=np.hstack((cpc10_o, cpc10))
            cpc3_o=np.hstack((cpc3_o, cpc3))
            uhsas100_o=np.hstack((uhsas100_o, pcasp))
        
        #%% read in flight data (for CSET and SOCRATES)
        elif campaign in ['CSET', 'SOCRATES']:
            filename = glob.glob(RFpath+'RF*'+date+'*.PNI.nc')
            if len(filename)==1 or len(filename)==2:  # SOCRATES has two flights in 20180217, choose the later one
                (time_cpc,cpc10,timeunit,cpc10unit,cpc10longname,cellsize,cellunit)=read_RF_NCAR(filename[-1],'CONCN')
                if campaign=='CSET':
                    (time_cpc,uhsas100,timeunit,uhsas100unit,uhsas100longname,cellsize,cellunit)=read_RF_NCAR(filename[-1],'CONCU100_RWOOU')
                elif campaign=='SOCRATES':
                    # there are two variables: CONCU100_CVIU and CONCU100_LWII
                    (time_cpc,uhsas100,timeunit,uhsas100unit,uhsas100longname,cellsize,cellunit)=read_RF_NCAR(filename[-1],'CONCU100_LWII')
            else:
                raise ValueError('find too many files: '+filename)
            
            # some quality checks
            uhsas100=qc_remove_neg(uhsas100)
            
            # exclude 30min after takeoff and before landing
            cpc10 = qc_mask_takeoff_landing(time_cpc,cpc10)
            uhsas100 = qc_mask_takeoff_landing(time_cpc,uhsas100)
            
            cpc10_o=np.hstack((cpc10_o, cpc10))
            uhsas100_o=np.hstack((uhsas100_o, uhsas100))
            
    
    #%% calculate statistics
    
    # select only valid data in obs and the corresponding data in models
    idx100 = ~np.isnan(uhsas100_o)
    idx10 = ~np.isnan(cpc10_o)
    idx3 = ~np.isnan(cpc3_o)
    
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
        
    if sum(idx10)/len(idx10)<0.1:   # two few observation available
        # for obs
        mean10[nmodels] = missing_value
        std10[nmodels] = missing_value
        p10_10[nmodels] = missing_value
        p25_10[nmodels] = missing_value
        p75_10[nmodels] = missing_value
        p90_10[nmodels] = missing_value
        # for models
        for mm in range(nmodels):
            mean10[mm] = np.nanmean(ncn10_m[mm][idx10])
            std10[mm] = np.nanstd(ncn10_m[mm][idx10])
            p10_10[mm] = np.nanpercentile(ncn10_m[mm][idx10],10)
            p25_10[mm] = np.nanpercentile(ncn10_m[mm][idx10],25)
            p75_10[mm] = np.nanpercentile(ncn10_m[mm][idx10],75)
            p90_10[mm] = np.nanpercentile(ncn10_m[mm][idx10],90)
            bias10[mm] =  missing_value
            corr10[mm] = [missing_value, missing_value]
            rmse10[mm] = missing_value
    else:
        # for obs
        mean10[nmodels] = np.nanmean(cpc10_o[idx10])
        std10[nmodels] = np.nanstd(cpc10_o[idx10])
        p10_10[nmodels] = np.nanpercentile(cpc10_o[idx10],10)
        p25_10[nmodels] = np.nanpercentile(cpc10_o[idx10],25)
        p75_10[nmodels] = np.nanpercentile(cpc10_o[idx10],75)
        p90_10[nmodels] = np.nanpercentile(cpc10_o[idx10],90)
        # for models
        for mm in range(nmodels):
            mean10[mm] = np.nanmean(ncn10_m[mm][idx10])
            std10[mm] = np.nanstd(ncn10_m[mm][idx10])
            p10_10[mm] = np.nanpercentile(ncn10_m[mm][idx10],10)
            p25_10[mm] = np.nanpercentile(ncn10_m[mm][idx10],25)
            p75_10[mm] = np.nanpercentile(ncn10_m[mm][idx10],75)
            p90_10[mm] = np.nanpercentile(ncn10_m[mm][idx10],90)
            bias10[mm] = mean10[mm] - mean10[nmodels]
            c10 = scipy.stats.pearsonr(ncn10_m[mm][idx10],cpc10_o[idx10])
            corr10[mm] = [c10[0],c10[1]]
            rmse10[mm] = np.sqrt(((ncn10_m[mm][idx10]-cpc10_o[idx10])**2).mean())
            
    if sum(idx100)/len(idx100)<0.1:   # two few observation available
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
        mean100[nmodels] = np.nanmean(uhsas100_o[idx100])
        std100[nmodels] = np.nanstd(uhsas100_o[idx100])
        p10_100[nmodels] = np.nanpercentile(uhsas100_o[idx100],10)
        p25_100[nmodels] = np.nanpercentile(uhsas100_o[idx100],25)
        p75_100[nmodels] = np.nanpercentile(uhsas100_o[idx100],75)
        p90_100[nmodels] = np.nanpercentile(uhsas100_o[idx100],90)
        # for models
        for mm in range(nmodels):
            mean100[mm] = np.nanmean(ncn100_m[mm][idx100])
            std100[mm] = np.nanstd(ncn100_m[mm][idx100])
            p10_100[mm] = np.nanpercentile(ncn100_m[mm][idx100],10)
            p25_100[mm] = np.nanpercentile(ncn100_m[mm][idx100],25)
            p75_100[mm] = np.nanpercentile(ncn100_m[mm][idx100],75)
            p90_100[mm] = np.nanpercentile(ncn100_m[mm][idx100],90)
            bias100[mm] = mean100[mm] - mean100[nmodels]
            c100 = scipy.stats.pearsonr(ncn100_m[mm][idx100],uhsas100_o[idx100])
            corr100[mm] = [c100[0],c100[1]]
            rmse100[mm] = np.sqrt(((ncn100_m[mm][idx100]-uhsas100_o[idx100])**2).mean())
            
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
            mean3[mm] = np.nanmean(ncn3_m[mm][idx3])
            std3[mm] = np.nanstd(ncn3_m[mm][idx3])
            p10_3[mm] = np.nanpercentile(ncn3_m[mm][idx3],10)
            p25_3[mm] = np.nanpercentile(ncn3_m[mm][idx3],25)
            p75_3[mm] = np.nanpercentile(ncn3_m[mm][idx3],75)
            p90_3[mm] = np.nanpercentile(ncn3_m[mm][idx3],90)
            bias3[mm] =  missing_value
            corr3[mm] = [missing_value, missing_value]
            rmse3[mm] = missing_value
    else:
        # for obs
        mean3[nmodels] = np.nanmean(cpc3_o[idx3])
        std3[nmodels] = np.nanstd(cpc3_o[idx3])
        p10_3[nmodels] = np.nanpercentile(cpc3_o[idx3],10)
        p25_3[nmodels] = np.nanpercentile(cpc3_o[idx3],25)
        p75_3[nmodels] = np.nanpercentile(cpc3_o[idx3],75)
        p90_3[nmodels] = np.nanpercentile(cpc3_o[idx3],90)
        # for models
        for mm in range(nmodels):
            mean3[mm] = np.nanmean(ncn3_m[mm][idx3])
            std3[mm] = np.nanstd(ncn3_m[mm][idx3])
            p10_3[mm] = np.nanpercentile(ncn3_m[mm][idx3],10)
            p25_3[mm] = np.nanpercentile(ncn3_m[mm][idx3],25)
            p75_3[mm] = np.nanpercentile(ncn3_m[mm][idx3],75)
            p90_3[mm] = np.nanpercentile(ncn3_m[mm][idx3],90)
            bias3[mm] = mean3[mm] - mean3[nmodels]
            c3 = scipy.stats.pearsonr(ncn3_m[mm][idx3],cpc3_o[idx3])
            corr3[mm] = [c3[0],c3[1]]
            rmse3[mm] = np.sqrt(((ncn3_m[mm][idx3]-cpc3_o[idx3])**2).mean())
    
    
    #%% write out files
    
    if campaign in ['HISCALE', 'ACEENA']:   
        outfile = figpath_aircraft_statistics+'statistics_CN10nm_'+campaign+'_'+IOP+'.txt'
    elif campaign in ['CSET', 'SOCRATES']:
        outfile = figpath_aircraft_statistics+'statistics_CN10nm_'+campaign+'.txt'
    
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
            
    
    if campaign in ['HISCALE', 'ACEENA']:   
        outfile = figpath_aircraft_statistics+'statistics_CN3nm_'+campaign+'_'+IOP+'.txt'
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
            
            
    if campaign in ['HISCALE', 'ACEENA']:   
        outfile = figpath_aircraft_statistics+'statistics_CN100nm_'+campaign+'_'+IOP+'.txt'
    elif campaign in ['CSET', 'SOCRATES']:
        outfile = figpath_aircraft_statistics+'statistics_CN100nm_'+campaign+'.txt'
    print('write statistics to file '+outfile)
    
    with open(outfile, 'w') as f:
        if campaign in ['CSET', 'SOCRATES']:
            f.write('statistics of Aerosol Number Concentration comparing with UHSAS(>100nm). sample size '+format(sum(idx100))+'\n')
        elif campaign=='ACEENA':
            f.write('statistics of Aerosol Number Concentration comparing with PCASP(>100nm). sample size '+format(sum(idx100))+'\n')
        elif campaign=='HISCALE':
            f.write('statistics of Aerosol Number Concentration comparing with PCASP(>120nm for IOP1, >90nm for IOP2). sample size '+format(sum(idx100))+'\n')
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