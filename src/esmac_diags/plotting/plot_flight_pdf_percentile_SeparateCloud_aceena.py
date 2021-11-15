"""
# plot_flight_pdf_percentile_SeparateCloud_aceena.py
# plot pdf and percentiles in several aerosol size bins for aircraft data
# separated by observed PBLH 
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
# from time_format_change import  hhmmss2sec,yyyymmdd2cday
from ..subroutines.read_aircraft import read_cpc
from ..subroutines.read_netcdf import read_merged_size,read_extractflight
from ..subroutines.quality_control import qc_remove_neg

def run_plot(settings):

    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    cpcpath = settings['cpcpath']
    merged_size_path = settings['merged_size_path']
    E3SM_aircraft_path = settings['E3SM_aircraft_path']
    figpath_aircraft_statistics = settings['figpath_aircraft_statistics']
    
    IOP = settings.get('IOP', None)

    #%% other settings
    if not os.path.exists(figpath_aircraft_statistics):
        os.makedirs(figpath_aircraft_statistics)
       
    # set final bin sizes
    binl = np.array([3, 15, 70, 300, 1000])
    binh = np.array([10, 70, 300, 1000, 3000])
    binm = (binl+binh)/2
    
    d_mam=np.arange(1,3001)
    blen = len(binm)
    
    # numbers of bins in merged size data
    b2len=67  
    
    #%% find files for flight information
    
    lst = glob.glob(merged_size_path+'merged_bin_*'+campaign+'*.nc')
    lst.sort()
    
    if len(lst)==0:
        raise ValueError('cannot find any file')
    
    # choose files for specific IOP
    if campaign=='ACEENA':
        if IOP=='IOP1':
            lst=lst[0:20]
        elif IOP=='IOP2':
            lst=lst[20:]
        elif IOP[0:4]=='2017' or IOP[0:4]=='2018':
            a=lst[0].split('_'+campaign+'_')
            lst = glob.glob(a[0]+'*'+IOP+'*')
            lst.sort()
    else:
        raise ValueError('this code is only for ACEENA, check the campaign settings')
    
    if len(lst)==0:
        raise ValueError('cannot find any file')
    
    #%% read all data
    
    # pdf average for legs
    pdf_sfc_obs=np.zeros([b2len,0])
    pdf_near_obs=np.zeros([b2len,0])
    pdf_above_obs=np.zeros([b2len,0])
    
    cpcdiff_sfc=np.zeros([0])
    cpcdiff_near=np.zeros([0])
    cpcdiff_above=np.zeros([0])
    
    nmodels=len(Model_List)
    pdf_sfc_model=[]
    pdf_near_model=[]
    pdf_above_model=[]
    for mm in range(nmodels):
        pdf_sfc_model.append(np.zeros([3000,0]))
        pdf_near_model.append(np.zeros([3000,0]))
        pdf_above_model.append(np.zeros([3000,0]))
    
    # pdf for the final bin sizes
    p2_sfc_obs = []
    p2_near_obs = []
    p2_above_obs = []
    p2_sfc_model=[]
    p2_near_model=[]
    p2_above_model=[]
    for mm in range(nmodels):
        p2_sfc_model.append([])
        p2_near_model.append([])
        p2_above_model.append([])
        
    print('reading '+format(len(lst))+' files to calculate the statistics: ')
    
    for filename in lst:
        
        # get date info:        
        date=filename[-12:-3]
        if date[-1]=='a':
            flightidx=1
        else:
            flightidx=2
        print(date)
        
        #%% read aerosol size distribution
        (time,size,cvi,timeunit,cunit,long_name)=read_merged_size(filename,'CVI_inlet')
        (time,size,cflag,timeunit,cunit,long_name)=read_merged_size(filename,'cld_flag')
        (time,size,legnum,timeunit,cunit,long_name)=read_merged_size(filename,'leg_number')
        (time,size,height,timeunit,zunit,long_name)=read_merged_size(filename,'height')
        (time,size,sizeh,timeunit,dataunit,long_name)=read_merged_size(filename,'size_high')
        (time,size,sizel,timeunit,dataunit,long_name)=read_merged_size(filename,'size_low')
        (time,size,merge,timeunit,dataunit,long_name)=read_merged_size(filename,'size_distribution_merged')
        time=np.ma.compressed(time)
        time=time/3600.
        size=np.ma.compressed(size)*1000  # um to nm
        sizel=sizel*1000
        sizeh=sizeh*1000
        merge=qc_remove_neg(merge)
        
        
    
        #%% read in CPC measurements
        
        if campaign=='ACEENA':
            filename_c=glob.glob(cpcpath+'CPC_G1_'+date[0:8]+'*R2_ACEENA001s.ict')    
        else:
            raise ValueError('this code is only for ACEENA, check the campaign settings')
        filename_c.sort()
        # read in data
        if len(filename_c)==1 or len(filename_c)==2: # some days have two flights
            (cpc,cpclist)=read_cpc(filename_c[flightidx-1])
            if np.logical_and(campaign=='ACEENA', date=='20180216a'):
                cpc=np.insert(cpc,1404,(cpc[:,1403]+cpc[:,1404])/2,axis=1)
            elif np.logical_and(campaign=='HiScale', date=='20160425a'):
                cpc=np.insert(cpc,0,cpc[:,0],axis=1)
                cpc[0,0]=cpc[0,0]-1
            time_cpc = cpc[0,:]/3600
            cpc10 = cpc[1,:]
            cpc3 = cpc[2,:]
        elif len(filename_c)==0:
            time_cpc=time
            cpc10=np.nan*np.empty([len(time)])
            cpc3=np.nan*np.empty([len(time)])
        else:
            raise ValueError('find too many files in ' + filename_c)
        
        cpcdiff = cpc3-cpc10
        cpcdiff=qc_remove_neg(cpcdiff)
        
        #%% read in Models
        datam2 = []
        for mm in range(nmodels):
            filename_m = E3SM_aircraft_path+'Aircraft_CNsize_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
            (timem,heightm,datam,timeunitm,datamunit,datamlongname)=read_extractflight(filename_m,'NCNall')
            datam2.append(datam*1e-6)    # #/m3 to #/cm3
        
        timem = (timem - int(timem[0]))*24
        
        if len(timem)!=len(time) or len(time)!=len(time_cpc):
            raise ValueError('time dimension for obs and/or model are not consistent')
        
        #%% get leg information near surface, near cloud base and above cloud
        # leg_sfc = np.ma.compressed (np.unique(legnum[height<=200])[1:])
        leg_sfc = list()
        leg_near = list()
        leg_above = list()
        leg_toomuchcld = list()
        leg_nocld = list()
        leg_nodata = list()
        
        for ii in range(1,max(legnum)+1):
        # for ii in range(5,7):
            idx_l = legnum==ii
            # if any(cflag[idx_l]==1):
            # make sure cloud flag less than 10% of time and FIMS is not all missing (e.g., TD mode)
            if np.sum(cflag[idx_l])/len(cflag[idx_l]) >= 0.1:
            # if np.sum(cflag[idx_l]) > 1:
                leg_toomuchcld.append(ii)
                continue
            if all(np.isnan(merge[idx_l,10])):
                leg_nodata.append(ii)
                continue
            
            legheight = np.mean(height[idx_l])
            # if legheight<=250:   # leg number near surface
            #     leg_sfc.append(ii)
                
            # find the mean cloud height within 1hr of the leg
            i = np.argwhere(legnum==ii)
            i_start = max(i[0][0]-3600, 0)
            i_end = min(i[-1][0]+3600, len(cflag))
            if all(cflag[i_start:i_end]!=1):
                leg_nocld.append(ii)
                if legheight>2500:
                    leg_above.append(ii)
                elif legheight<=250:   # leg number near surface
                    leg_sfc.append(ii)
                continue
            idx_c = cflag[i_start:i_end]==1
            cldheight = np.mean(height[i_start:i_end][idx_c])
            cldmax = np.max(height[i_start:i_end][idx_c])
            cldmin = np.min(height[i_start:i_end][idx_c])
            # if (legheight-cldheight)<=200 and (legheight-cldheight)>=-400:
            if legheight>=max(cldmin,250) and legheight<=cldmax:
                leg_near.append(ii)
            elif legheight<min(cldmin,250):   # leg number near surface
                leg_sfc.append(ii)
            # elif (legheight-cldheight)>500:
            elif legheight>cldmax:
                leg_above.append(ii)
    
        #%% calculate all pdfs
        for ii in range(len(leg_sfc)):
            idx = legnum==leg_sfc[ii]
            tmp_obs = np.nanmean(merge[idx,:],0)
            tmp_obs[tmp_obs==0]=np.nan
            pdf_sfc_obs = np.hstack((pdf_sfc_obs, np.reshape(tmp_obs,(b2len,1))))
            cpcdiff_sfc = np.hstack((cpcdiff_sfc, np.nanmean(cpcdiff[idx])))
            for mm in range(nmodels):
                tmp_model = np.nanmean(datam2[mm][:,idx],1)
                tmp_model[tmp_model==0]=np.nan
                pdf_sfc_model[mm] = np.hstack((pdf_sfc_model[mm], np.reshape(tmp_model,(3000,1))))
            
        for ii in range(len(leg_near)):
            idx = legnum==leg_near[ii]
            tmp_obs = np.nanmean(merge[idx,:],0)
            tmp_obs[tmp_obs==0]=np.nan
            pdf_near_obs = np.hstack((pdf_near_obs, np.reshape(tmp_obs,(b2len,1))))
            cpcdiff_near = np.hstack((cpcdiff_near, np.nanmean(cpcdiff[idx])))
            for mm in range(nmodels):
                tmp_model = np.nanmean(datam2[mm][:,idx],1)
                tmp_model[tmp_model==0]=np.nan
                pdf_near_model[mm] = np.hstack((pdf_near_model[mm], np.reshape(tmp_model,(3000,1))))
            
        for ii in range(len(leg_above)):
            idx = legnum==leg_above[ii]
            tmp_obs = np.nanmean(merge[idx,:],0)
            tmp_obs[tmp_obs==0]=np.nan
            pdf_above_obs = np.hstack((pdf_above_obs, np.reshape(tmp_obs,(b2len,1))))
            cpcdiff_above = np.hstack((cpcdiff_above, np.nanmean(cpcdiff[idx])))
            for mm in range(nmodels):
                tmp_model = np.nanmean(datam2[mm][:,idx],1)
                tmp_model[tmp_model==0]=np.nan
                pdf_above_model[mm] = np.hstack((pdf_above_model[mm], np.reshape(tmp_model,(3000,1))))
        
    
    #%% change to the pre-defined size bins
        
    for bb in range(blen):
        idx_m = np.logical_and(d_mam>=binl[bb], d_mam<=binh[bb])
        for mm in range(nmodels):
            p2_sfc_model[mm].append(np.nansum(pdf_sfc_model[mm][idx_m,:],0))
            p2_near_model[mm].append(np.nansum(pdf_near_model[mm][idx_m,:],0))
            p2_above_model[mm].append(np.nansum(pdf_above_model[mm][idx_m,:],0))
    
        if bb==0:
            p2_sfc_obs.append(cpcdiff_sfc[~np.isnan(cpcdiff_sfc)])
            p2_near_obs.append(cpcdiff_near[~np.isnan(cpcdiff_near)])
            p2_above_obs.append(cpcdiff_above[~np.isnan(cpcdiff_above)])
        else:
            idx_o = np.logical_and(sizel>=binl[bb], sizeh<=binh[bb])
            if any(idx_o):
                tmp_sfc = np.nansum(pdf_sfc_obs[idx_o,:],0)
                tmp_near = np.nansum(pdf_near_obs[idx_o,:],0)
                tmp_above = np.nansum(pdf_above_obs[idx_o,:],0)
                p2_sfc_obs.append(tmp_sfc[tmp_sfc!=0])
                p2_near_obs.append(tmp_near[tmp_near!=0])
                p2_above_obs.append(tmp_above[tmp_above!=0])
            else:
               raise ValueError("no sample is found in the size bin")
                
    #%% calculate dlnDp for dN/dlnDp
    d_mam=np.arange(1,3001)
    dlnDp_m=np.full(3000,np.nan)
    for bb in range(3000):
        dlnDp_m[bb]=np.log((bb+2)/(bb+1))
    dlnDp_o=np.empty(len(size))
    for bb in range(len(size)):
        dlnDp_o[bb]=np.log(sizeh[bb]/sizel[bb])
    
    
    #%% plot entire pdf below and above PBL
    figname = figpath_aircraft_statistics+'SeparateCloud_pdf_AerosolSize_ACEENA_'+IOP+'.png'
    print('plotting PDF figures to '+figname)
    
    fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(6,8))
    
    ax1.plot(size,np.nanmedian(pdf_above_obs,1)/dlnDp_o,color='k',linewidth=1,label='Obs')
    for mm in range(nmodels):
        ax1.plot(d_mam,np.nanmedian(pdf_above_model[mm],1)/dlnDp_m,color=color_model[mm],linewidth=1, label=Model_List[mm])
    ax1.tick_params(color='k',labelsize=14)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    
    ax2.plot(size,np.nanmedian(pdf_near_obs,1)/dlnDp_o,color='k',linewidth=1,label='Obs')
    for mm in range(nmodels):
        ax2.plot(d_mam,np.nanmedian(pdf_near_model[mm],1)/dlnDp_m,color=color_model[mm],linewidth=1, label=Model_List[mm])
    ax2.tick_params(color='k',labelsize=14)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    
    ax3.plot(size,np.nanmedian(pdf_sfc_obs,1)/dlnDp_o,color='k',linewidth=1,label='Obs')
    for mm in range(nmodels):
        ax3.plot(d_mam,np.nanmedian(pdf_sfc_model[mm],1)/dlnDp_m,color=color_model[mm],linewidth=1, label=Model_List[mm])
    ax3.tick_params(color='k',labelsize=14)
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    
    # ax0.set_xlim(5,4000)
    # ax1.set_xlim(5,4000)
    ax1.set_ylim(1e-3,1e5)
    ax2.set_ylim(1e-3,1e5)
    ax3.set_ylim(1e-3,1e5)
    
    ax2.set_ylabel('aerosol #/dlnDp (cm$^{-3}$)',fontsize=14)
    ax3.set_xlabel('Diameter (nm)',fontsize=14)
    l=ax3.legend(loc='lower center', shadow=False, fontsize='medium')
    
    ax1.set_title('size distribution for ACEENA '+IOP,fontsize=15)
    
    ax3.text(200,3000,'Near Surface ('+str(pdf_sfc_obs.shape[1])+' legs)',fontsize=12)
    ax2.text(200,3000,'Near Clouds ('+str(pdf_near_obs.shape[1])+' legs)',fontsize=12)
    ax1.text(200,3000,'Above Clouds ('+str(pdf_above_obs.shape[1])+' legs)',fontsize=12)
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    # plt.close()
    
    #%% plot percentile on sizes
    
    figname = figpath_aircraft_statistics+'SeparateCloud_percentile_AerosolSize_ACEENA_'+IOP+'.png'
    print('plotting percentile figures to '+figname)
    
    # set position shift so that models and obs are not overlapped
    p_shift = np.arange(nmodels+1)
    p_shift = (p_shift - p_shift.mean())*0.2
    
    fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(6,8))
        
    ax1.boxplot(p2_above_obs,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(blen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax1.boxplot(p2_above_model[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(blen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    
    ax2.boxplot(p2_near_obs,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(blen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax2.boxplot(p2_near_model[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(blen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
        
    ax3.boxplot(p2_sfc_obs,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(blen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax3.boxplot(p2_sfc_model[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(blen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
        
    ax3.tick_params(color='k',labelsize=12)
    ax2.tick_params(color='k',labelsize=12)
    ax1.tick_params(color='k',labelsize=12)
    ax3.set_yscale('log')
    ax2.set_yscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(-.5,blen-.5)
    ax2.set_xlim(-.5,blen-.5)
    ax3.set_xlim(-.5,blen-.5)
    
    ax3.set_xlabel('Diameter (nm)',fontsize=14)
    ax2.set_ylabel('aerosol # (cm$^{-3}$)',fontsize=14)
    ax1.set_title('percentile for ACEENA '+IOP,fontsize=15)
    
    ax3.text(2.4,4000,'Near Surface ('+str(pdf_sfc_obs.shape[1])+' legs)',fontsize=12)
    ax2.text(2.4,4000,'Near Clouds ('+str(pdf_near_obs.shape[1])+' legs)',fontsize=12)
    ax1.text(2.4,4000,'Above Clouds ('+str(pdf_above_obs.shape[1])+' legs)',fontsize=12)
    
    xlabel=[str(binl[x])+'-'+str(binh[x]) for x in range(blen)]
    ax1.set_xticks(range(len(binm)))
    ax1.set_xticklabels(xlabel)
    ax2.set_xticks(range(len(binm)))
    ax2.set_xticklabels(xlabel)
    ax3.set_xticks(range(len(binm)))
    ax3.set_xticklabels(xlabel)
    ax1.set_ylim(1e-3,1e5)
    ax2.set_ylim(1e-3,1e5)
    ax3.set_ylim(1e-3,1e5)
    
    # plot temporal lines for label
    ax3.plot([],c='k',label='Obs')
    for mm in range(nmodels):
        ax3.plot([],c=color_model[mm],label=Model_List[mm])
    ax3.legend(loc='lower left', shadow=False, fontsize='medium')
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    # plt.close()
    
