"""
# plot percentile of CCN number concentration with height
# for flight data in IOPs
# compare models and CCN measurements
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_aircraft import read_ccn_hiscale, read_ccn_socrates
from ..subroutines.read_ARMdata import read_ccn
from ..subroutines.read_netcdf import read_extractflight,read_merged_size
from ..subroutines.quality_control import qc_mask_qcflag,qc_mask_cloudflag,qc_remove_neg

def run_plot(settings):

    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    height_bin = settings['height_bin']
    E3SM_aircraft_path = settings['E3SM_aircraft_path']
    figpath_aircraft_statistics = settings['figpath_aircraft_statistics']
    
    if campaign in ['HISCALE', 'ACEENA']:
        IOP = settings.get('IOP', None)
        merged_size_path = settings.get('merged_size_path', None)
        ccnpath = settings.get('ccnpath', None)
    elif campaign in ['CSET', 'SOCRATES']:
        ccnpath = settings.get('ccnpath', None)
    else:
        raise ValueError('campaign name is not recognized: '+campaign)

    #%% other settings
    
    if not os.path.exists(figpath_aircraft_statistics):
        os.makedirs(figpath_aircraft_statistics)
       
        
    #%%
    z=height_bin
    dz = z[1]-z[0]
    zmin=z-np.insert((z[1:]-z[0:-1])/2,0,dz)
    zmax=z+np.append((z[1:]-z[0:-1])/2,dz)
    
    zlen=len(z)   
    
    
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
    
    height_all = []
    ccna_all = []
    ccnb_all = []
    SSa_all = np.array([])
    SSb_all = np.array([])
    ccn3_all = []
    ccn5_all = []
    nmodels=len(Model_List)
    for mm in range(nmodels):
        ccn3_all.append([])
        ccn5_all.append([])
        
    print('reading '+format(len(alldates))+' files to calculate the statistics: ')
    
    for date in alldates:
        print(date)
        
        #%% read in Models
        
        ccn3=[]
        ccn5=[]
        for mm in range(nmodels):
            filename_m = E3SM_aircraft_path+'Aircraft_vars_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
        
            (timem,heightm,ccn3_tmp,timeunitm,ccn3_unit,ccn3_longname)=read_extractflight(filename_m,'CCN3')
            (timem,heightm,ccn5_tmp,timeunitm,ccn5_unit,ccn5_longname)=read_extractflight(filename_m,'CCN5')
            ccn3.append(ccn3_tmp)
            ccn5.append(ccn5_tmp)
            
            
        # get supersaturation
        SS3 = ccn3_longname.split('=')[-1]
        SS5 = ccn5_longname.split('=')[-1]
        
        #%% read in flight data (for HISCALE)
        if campaign=='HISCALE':
            filename_ccn=glob.glob(ccnpath+'CCN_G1_'+date[0:8]+'*R2_HiScale001s.*')
            filename_ccn.sort()
            if date[-1]=='a':
                flightidx=1
            else:
                flightidx=2
            # read in data
            if len(filename_ccn)==1 or len(filename_ccn)==2:
                (data0,ccnlist)=read_ccn_hiscale(filename_ccn[flightidx-1])
                # only choose data quality is good (flag=0)
                flag = data0[7,:]
                time_ccn = data0[0,:]
                ccna = data0[10,:]
                ccnb = data0[11,:]
                SSa = data0[2,:]
                SSb = data0[5,:]
                ccna = qc_mask_qcflag(ccna,flag)
                ccnb = qc_mask_qcflag(ccnb,flag)
                SSa=qc_remove_neg(SSa)
                SSb=qc_remove_neg(SSb)
            elif len(filename_ccn)==0:
                time_ccn=timem
                ccna=np.nan*np.empty([len(timem)])
                ccnb=np.nan*np.empty([len(timem)])
                SSa=0.24*np.full(len(timem),1)
                SSb=0.46*np.full(len(timem),1)
            else:
                raise ValueError('find too many files: '+filename_ccn)
            timea=time_ccn
            timeb=time_ccn
            # cloud flag
            filename = merged_size_path+'merged_bin_fims_pcasp_'+campaign+'_'+date+'.nc'
            (time,size,cflag,timeunit,cunit,long_name)=read_merged_size(filename,'cld_flag')
            ccna = qc_mask_cloudflag(ccna,cflag)
            ccnb = qc_mask_cloudflag(ccnb,cflag)
            
        elif campaign=='ACEENA':
            filename_ccna=glob.glob(ccnpath+'enaaafccn2colaF1.b1.'+date[0:8]+'*.nc')
            filename_ccnb=glob.glob(ccnpath+'enaaafccn2colbF1.b1.'+date[0:8]+'*.nc')
            # read in data
            if len(filename_ccna)==1:
                (timea,timeunita,ccna,qcflag,ccnunit,SSa)=read_ccn(filename_ccna[0])
                ccna=qc_mask_qcflag(ccna,qcflag)
                ccna=qc_remove_neg(ccna)
                SSa=qc_remove_neg(SSa)
            elif len(filename_ccna)==0:
                # print('no CCN data found. set as NaN')
                timea=timem
                SSa=np.nan*np.empty([len(timem)])
                ccna=np.nan*np.empty([len(timem)])
            else:
                raise ValueError('find too many files: '+filename_ccna)
            if len(filename_ccnb)==1:
                (timeb,timeunitb,ccnb,qcflag,ccnunit,SSb)=read_ccn(filename_ccnb[0])
                ccnb=qc_mask_qcflag(ccnb,qcflag)
                ccnb=qc_remove_neg(ccnb)
                SSb=qc_remove_neg(SSb)
            elif len(filename_ccnb)==0:
                # print('no CCN data found. set as NaN')
                timeb=timem
                SSb=np.nan*np.empty([len(timem)])
                ccnb=np.nan*np.empty([len(timem)])
            else:
                raise ValueError('find too many files: '+filename_ccnb)
            # cloud flag
            filename = merged_size_path+'merged_bin_fims_pcasp_opc_'+campaign+'_'+date+'.nc'
            (time,size,cflag,timeunit,cunit,long_name)=read_merged_size(filename,'cld_flag')
            if date=='20170707a':
                time=np.delete(time,5247)
                cflag=np.delete(cflag,5247)
            elif date=='20180201a':
                time=np.delete(time,3635)
                cflag=np.delete(cflag,3635)
            if time[0]<timea[0]:
                cflag = cflag[np.where(time==timea[0])[0][0]:]
                time = time[np.where(time==timea[0])[0][0]:]
            elif time[0]>timea[0]:
                cflag = np.insert(cflag,np.full(int(time[0]-timea[0]),0), -9999)
                time = np.insert(time,np.full(int(time[0]-timea[0]),0), -9999)
            if time[-1]<timea[-1]:
                cflag = np.append(cflag,np.full(int(timea[-1]-time[-1]), -9999))
                time = np.append(time,np.full(int(timea[-1]-time[-1]), -9999))
            elif time[-1]>timea[-1]:
                cflag = cflag[0:np.where(time==timea[-1])[0][0]+1]
                time = time[0:np.where(time==timea[-1])[0][0]+1]
            ccna = qc_mask_cloudflag(ccna,cflag)
            ccnb = qc_mask_cloudflag(ccnb,cflag)
            
        # CSET does not have observed CCN
        elif campaign=='CSET':
            timea=timem
            SSa=np.nan*np.empty([len(timem)])
            ccna=np.nan*np.empty([len(timem)])
            timeb=timem
            SSb=np.nan*np.empty([len(timem)])
            ccnb=np.nan*np.empty([len(timem)])
            
        # SOCRATES
        elif campaign=='SOCRATES':
            filename_ccn=glob.glob(ccnpath+'CCNscanning_SOCRATES_GV_RF*'+date[0:8]+'_R0.ict')
            if len(filename_ccn)==1:
                (data0,ccnlist)=read_ccn_socrates(filename_ccn[0])
                time_ccn = data0[0,:]
                ccn = data0[1,:]
                SS = data0[3,:]
                ccn=qc_remove_neg(ccn)
                timea=time_ccn
                timeb=time_ccn
                ccna=np.array(ccn)
                ccnb=np.array(ccn)
                idxa=np.logical_and(SS>0.05, SS<0.15)
                ccna[idxa==False]=np.nan
                SSa=np.full((len(timea)),0.1)
                idxb=np.logical_and(SS>0.45, SS<0.55)
                ccnb[idxb==False]=np.nan
                SSb=np.full((len(timeb)),0.5)
            elif len(filename_ccn)==0:
                timea=timem
                SSa=np.nan*np.empty([len(timem)])
                ccna=np.nan*np.empty([len(timem)])
                timeb=timem
                SSb=np.nan*np.empty([len(timem)])
                ccnb=np.nan*np.empty([len(timem)])
            else:
                raise ValueError('find too many files: '+filename_ccn)
                
        if any(timea!=timeb):
            raise ValueError('time dimension is inconsistent')
            
        # exclude NaNs
        idx = np.logical_or(~np.isnan(ccna), ~np.isnan(ccnb))
        ccna_all.append(ccna[idx])
        ccnb_all.append(ccnb[idx])
        SSa_all=np.append(SSa_all,SSa[idx])
        SSb_all=np.append(SSb_all,SSb[idx])
        
        height2=np.interp(timea,timem,heightm)
        height_all.append(height2[idx])
        
        # for interpolation of model results
        timea=timea[idx]
        timeb=timeb[idx]
        
        # interpolate model results into observational time
        for mm in range(nmodels):
            ccn3_all[mm].append(np.interp(timea,timem,ccn3[mm])) 
            ccn5_all[mm].append(np.interp(timeb,timem,ccn5[mm])) 
             
    #%% calculate percentiles for each height bin
    
    ccna_z = list()
    ccnb_z = list()
    ccn3_z = []
    ccn5_z = []
    nmodels=len(Model_List)
    for mm in range(nmodels):
        ccn3_z.append([])
        ccn5_z.append([])
    for zz in range(zlen):
        ccna_z.append(np.empty(0))
        ccnb_z.append(np.empty(0))
        for mm in range(nmodels):
            ccn3_z[mm].append(np.empty(0))
            ccn5_z[mm].append(np.empty(0))
        
    ndays=len(height_all)
    for dd in range(ndays):
        height = height_all[dd]
        ccna = ccna_all[dd]
        ccnb = ccnb_all[dd]
        for zz in range(zlen):
            idx = np.logical_and(height>=zmin[zz], height<zmax[zz])
            ccna_z[zz]=np.append(ccna_z[zz],ccna[np.logical_and(idx,~np.isnan(ccna))])
            ccnb_z[zz]=np.append(ccnb_z[zz],ccnb[np.logical_and(idx,~np.isnan(ccnb))])
            for mm in range(nmodels):
                ccn3 = ccn3_all[mm][dd]
                ccn5 = ccn5_all[mm][dd]
                ccn3_z[mm][zz]=np.append(ccn3_z[mm][zz],ccn3[idx])
                ccn5_z[mm][zz]=np.append(ccn5_z[mm][zz],ccn5[idx])
            
    
    #%% make plot
    # set position shift so that models and obs are not overlapped
    p_shift = np.arange(nmodels+1)
    p_shift = (p_shift - p_shift.mean())*0.2
    
    if campaign in ['HISCALE', 'ACEENA']:
        figname = figpath_aircraft_statistics+'percentile_height_CCN_'+campaign+'_'+IOP+'.png'
    else:
        figname = figpath_aircraft_statistics+'percentile_height_CCN_'+campaign+'.png'
    print('plotting figures to '+figname)
    
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(8,8))   # figsize in inches
    # plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
        
    ax1.boxplot(ccna_z,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(zlen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=False, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax1.boxplot(ccn3_z[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(zlen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=False, patch_artist=True)    # need patch_artist to fill color in box
    ax1.tick_params(color='k',labelsize=16)
    #ax1.set_xscale('log')
    ax1.set_ylim(-1,zlen)
    ax1.set_yticks(range(zlen))
    ax1.set_yticklabels(z)
    # ax1.set_yticks([1,3,5,7,9,11,12,13,14,15,16])
    # ax1.set_yticklabels(range(400,4100,400))
    # plot temporal lines for label
    ax1.plot([],c='k',label='Obs ('+format(np.nanmean(SSa_all),'.2f')+'%)')
    for mm in range(nmodels):
        ax1.plot([],c=color_model[mm],label=Model_List[mm])
    ax1.legend(loc='upper right', fontsize='x-large')
        
    ax2.boxplot(ccnb_z,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(zlen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=False, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax2.boxplot(ccn5_z[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(zlen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=False, patch_artist=True)    # need patch_artist to fill color in box
    ax2.tick_params(color='k',labelsize=16)
    #ax2.set_xscale('log')
    ax2.set_ylim(-1,zlen)
    ax2.set_yticks(range(zlen))
    ax2.set_yticklabels([])
    # ax1.set_yticks(np.arange(0,20,2))
    # ax1.set_yticklabels(range(400,4100,400))
    # plot temporal lines for label
    ax2.plot([],c='k',label='Obs ('+format(np.nanmean(SSb_all),'.2f')+'%)')
    for mm in range(nmodels):
        ax2.plot([],c=color_model[mm],label=Model_List[mm])
    ax2.legend(loc='upper right', fontsize='x-large')
        
    # set xlimit consistent in subplots
    # xlim1 = ax1.get_xlim()
    # xlim2 = ax2.get_xlim()
    # ax1.set_xlim([min(xlim1[0],xlim2[0]), max(xlim1[1],xlim2[1])])
    # ax2.set_xlim([min(xlim1[0],xlim2[0]), max(xlim1[1],xlim2[1])])
    
    ax1.set_ylabel('Height (m MSL)',fontsize=16)
    fig.text(0.4,0.06, 'CCN number (cm$^{-3}$)', fontsize=16)
    ax1.set_title('SS = '+SS3,fontsize=16)
    ax2.set_title('SS = '+SS5,fontsize=16)
    if campaign in ['HISCALE', 'ACEENA']:
        fig.text(0.48,0.92, IOP, fontsize=18)
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    # plt.close()
        
