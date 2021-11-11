"""
# plot pdf and percentiles in several aerosol size bins for aircraft data
# separated by observed PBLH 
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.time_format_change import  hhmmss2sec,yyyymmdd2cday
from ..subroutines.read_ARMdata import read_pblhtmpl1
from ..subroutines.read_surface import read_dl_pblh
from ..subroutines.read_aircraft import read_cpc
from ..subroutines.read_netcdf import read_merged_size,read_extractflight
from ..subroutines.quality_control import qc_remove_neg, qc_mask_qcflag

def run_plot(settings):

    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    cpcpath = settings['cpcpath']
    pblhpath = settings['pblhpath']
    dlpath = settings['dlpath']
    merged_size_path = settings['merged_size_path']
    E3SM_aircraft_path = settings['E3SM_aircraft_path']
    figpath_aircraft_statistics = settings['figpath_aircraft_statistics']
    
    IOP = settings.get('IOP', None)

    #%% other settings
    if not os.path.exists(figpath_aircraft_statistics):
        os.makedirs(figpath_aircraft_statistics)
       
    # set final bin sizes
    binl = np.array([3, 15, 70, 400, 1000])
    binh = np.array([10, 70, 400, 1000, 3000])
    binm = (binl+binh)/2
    
    # set a range around PBLH (PBLH +/- heightdiff) that only data outside of the range are counted
    heightdiff = 100
       
    #%% read in doppler lidar data. this is all days in one file
    dl=read_dl_pblh(dlpath+'sgpdlC1_mlh_0.08.txt')
    
    mlh_dl = dl[6,:]*1000
    day_dl = np.array(mlh_dl[:])
    time_dl = np.array(mlh_dl[:])
    for tt in range(len(time_dl)):
        yyyymmdd=format(int(dl[0,tt]),'04d')+format(int(dl[1,tt]),'02d')+format(int(dl[2,tt]),'02d')
        hhmmss=format(int(dl[3,tt]),'02d')+':'+format(int(dl[4,tt]),'02d')+':'+format(int(dl[5,tt]),'02d')
        day_dl[tt]=yyyymmdd2cday(yyyymmdd)
        time_dl[tt]=hhmmss2sec(hhmmss)
    mlh_dl=qc_remove_neg(mlh_dl)
    
    
    #%% find files for flight information
    
    lst = glob.glob(merged_size_path+'merged_bin_*'+campaign+'*.nc')
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
            a=lst[0].split('_'+campaign+'_')
            lst = glob.glob(a[0]+'*'+IOP+'*')
            lst.sort()
    else:
        raise ValueError('this code is only for HISCALE, check the campaign settings')
    
    if len(lst)==0:
        raise ValueError('cannot find any file')
        
    #%% read all data
    
    # pdf average for legs
    pdf_below_obs=np.full([44,len(lst)*10],np.nan)
    pdf_above_obs=np.full([44,len(lst)*10],np.nan)
    
    cpcdiff_above=np.full([len(lst)*10],np.nan)
    cpcdiff_below=np.full([len(lst)*10],np.nan)
    
    nmodels=len(Model_List)
    pdf_below_model=[]
    pdf_above_model=[]
    for mm in range(nmodels):
        pdf_below_model.append(np.full([3000,len(lst)*10],np.nan))
        pdf_above_model.append(np.full([3000,len(lst)*10],np.nan))
    
    n_below=0
    n_above=0
    n_total=0
        
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
        merge=merge.T
        
        #%% read in CPC measurements
        
        if campaign=='HISCALE':
            filename_c=glob.glob(cpcpath+'CPC_G1_'+date[0:8]+'*R2_HiScale001s.ict.txt')
        else:
            raise ValueError('this code is only for HISCALE, check the campaign settings')
        filename_c.sort()
        # read in data
        if len(filename_c)==1 or len(filename_c)==2: # some days have two flights
            (cpc,cpclist)=read_cpc(filename_c[flightidx-1])
            if np.logical_and(campaign=='ACEENA', date=='20180216a'):
                cpc=np.insert(cpc,1404,(cpc[:,1403]+cpc[:,1404])/2,axis=1)
            elif np.logical_and(campaign=='HISCALE', date=='20160425a'):
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
            raise ValueError('find too many files')
        
        cpcdiff = cpc3-cpc10
        cpcdiff=qc_remove_neg(cpcdiff)
        
        
        #%% read in PBLH data from MPL
        filename_mpl=glob.glob(pblhpath+'sgppblhtmpl1sawyerliC1*'+date[0:8]+'*.nc')
        # read in data
        if len(filename_mpl)==1:
            (time_pblh,timeunit,mpl,qc_mpl) = read_pblhtmpl1(filename_mpl[0])
            mpl = qc_mask_qcflag(mpl, qc_mpl)
        elif len(filename_mpl)==0:
            print('no pblh file in this day. skip...')
            continue
        else:
            raise ValueError('find too many files: ' + filename_mpl)
        time_pblh=time_pblh/3600
        
        #%% choose the same time of DL. get pblh
        cday0=yyyymmdd2cday(date[0:8])
        idx_dl = day_dl==cday0
        time_dl2 = time_dl[idx_dl]/3600
        mlh_dl2 = mlh_dl[idx_dl]
        
        #%% read in Models
        datam2 = []
        for mm in range(nmodels):
            filename_m = E3SM_aircraft_path+'Aircraft_CNsize_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
            (timem,heightm,datam,timeunitm,datamunit,datamlongname)=read_extractflight(filename_m,'NCNall')
            datam2.append(datam*1e-6)    # #/m3 to #/cm3
        
        timem = timem/3600
        
        if len(timem)!=len(time) or len(time)!=len(time_cpc):
            raise ValueError('time dimension for obs and/or model are not consistent')
        
        #%% get pdf for legs below and above PBLH
        
        for ii in range(max(legnum)):
            # get the mean pblh for this leg
            time_leg=time[legnum==ii+1]
            cflag_leg=cflag[legnum==ii+1]
            if np.sum(cflag_leg==1)>1: #0.01*len(cflag_leg):
                continue  # don't use legs with >10% cloud flag
            
            idx_dl2 = np.logical_and(time_dl2>=time_leg[0], time_dl2<=time_leg[-1])
            if idx_dl2.any()==False:
                idx_dl2 = np.logical_and(time_dl2>=time_leg[0]-2, time_dl2<=time_leg[-1]+2) # extend time range
            if idx_dl2.any():
                pblh = np.nanmean(mlh_dl2[idx_dl2])
            else:# use MPL pblh
                idx_mpl = np.logical_and(time_pblh>=time_leg[0], time_pblh<time_leg[-1])
                if any(idx_mpl)==False:
                    idx_mpl=np.logical_and(time_pblh>=time_leg[0]-2, time_pblh<time_leg[-1]+2)
                pblh=np.mean(mpl[idx_mpl])
                
            
            # average for each legs first
            hmean = np.mean(height[legnum==ii+1])
            if hmean<pblh-heightdiff:  # below PBLH
                pdf_below_obs[:,n_below] = np.nanmean(merge[:,legnum==ii+1],1)
                cpcdiff_below[n_below] = np.nanmean(cpcdiff[legnum==ii+1])
                for mm in range(nmodels):
                    pdf_below_model[mm][:,n_below] = np.nanmean(datam2[mm][:,legnum==ii+1],1)
                n_below=n_below+1
            elif hmean>pblh+heightdiff:
                pdf_above_obs[:,n_above] = np.nanmean(merge[:,legnum==ii+1],1)
                cpcdiff_above[n_above] = np.nanmean(cpcdiff[legnum==ii+1])
                for mm in range(nmodels):
                    pdf_above_model[mm][:,n_above] = np.nanmean(datam2[mm][:,legnum==ii+1],1)
                n_above=n_above+1
                
                
    #%% change to the pre-defined size bins
            
    d_model=np.arange(1,3001)
    blen = len(binm)
    p2_below_obs = list()
    p2_above_obs = list()
    p2_above_model = list()
    p2_below_model = list()
    for mm in range(nmodels):
        p2_above_model.append([])
        p2_below_model.append([])
    
    for bb in range(blen):
        idx_m = np.logical_and(d_model>=binl[bb], d_model<=binh[bb])
        for mm in range(nmodels):
            data_below = np.nansum(pdf_below_model[mm][idx_m,:],0)
            data_above = np.nansum(pdf_above_model[mm][idx_m,:],0)
            # exclude pre-assigned data space that are not used
            p2_below_model[mm].append(data_below[range(n_below)])
            p2_above_model[mm].append(data_above[range(n_above)])
        if bb==0:
            p2_below_obs.append(cpcdiff_below[~np.isnan(cpcdiff_below)])
            p2_above_obs.append(cpcdiff_above[~np.isnan(cpcdiff_above)])
        else:
            idx_o = np.logical_and(sizel>=binl[bb], sizeh<=binh[bb])
            if any(idx_o):
                tmp_below = np.nansum(pdf_below_obs[idx_o,:],0)
                tmp_above = np.nansum(pdf_above_obs[idx_o,:],0)
                # exclude not used or not detected (0 value) data
                p2_below_obs.append(tmp_below[tmp_below!=0])
                p2_above_obs.append(tmp_above[tmp_above!=0])
            else:
                p2_below_obs.append(np.full([n_below],np.nan))
                p2_above_obs.append(np.full([n_above],np.nan))
    
    #%% change to dN/dlnDp
    # model
    dlnDp=np.empty(3000)
    for bb in range(3000):
        dlnDp[bb]=np.log((bb+2)/(bb+1))
    for nn in range(n_below):
        for mm in range(nmodels):
            pdf_below_model[mm][:,nn]=pdf_below_model[mm][:,nn]/dlnDp
    for nn in range(n_above):
        for mm in range(nmodels):
            pdf_above_model[mm][:,nn]=pdf_above_model[mm][:,nn]/dlnDp
        
    # Obs
    dlnDp=np.empty(len(size))
    for bb in range(len(size)):
        dlnDp[bb]=np.log(sizeh[bb]/sizel[bb])
    for nn in range(n_below):
        pdf_below_obs[:,nn]=pdf_below_obs[:,nn]/dlnDp  
    for nn in range(n_above):  
        pdf_above_obs[:,nn]=pdf_above_obs[:,nn]/dlnDp
        
         
    #%% plot entire pdf below and above PBL
    figname = figpath_aircraft_statistics+'SeparatePBLH_pdf_AerosolSize_HISCALE_'+IOP+'.png'
    print('plotting PDF figures to '+figname)
    
    fig,(ax0,ax1) = plt.subplots(2,1,figsize=(8,6))
    idx_v=range(n_above)
    h3=ax0.plot(size,np.nanmean(pdf_above_obs[:,idx_v],1),color='k',linewidth=1,label='Obs')
    for mm in range(nmodels):
        ax0.plot(np.arange(1,3001),np.nanmean(pdf_above_model[mm][:,idx_v],1),color=color_model[mm],linewidth=1, label=Model_List[mm])
    # ax0.legend(loc='lower center', shadow=False, fontsize='large')
    ax0.tick_params(color='k',labelsize=14)
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    
    idx_v=range(n_below)
    h3=ax1.plot(size,np.nanmean(pdf_below_obs[:,idx_v],1),color='k',linewidth=1,label='Obs')
    for mm in range(nmodels):
        ax1.plot(np.arange(1,3001),np.nanmean(pdf_below_model[mm][:,idx_v],1),color=color_model[mm],linewidth=1, label=Model_List[mm])
    ax1.legend(loc='lower left', shadow=False, fontsize='large')
    ax1.tick_params(color='k',labelsize=14)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    
    # ax0.set_xlim(5,4000)
    # ax1.set_xlim(5,4000)
    ax0.set_ylim(1e-3,1e5)
    ax1.set_ylim(1e-3,1e5)
    ax1.set_xlabel('Diameter (nm)',fontsize=14)
    ax0.set_ylabel('aerosol #/dlnDp (cm$^{-3}$)',fontsize=13)
    ax1.set_ylabel('aerosol #/dlnDp (cm$^{-3}$)',fontsize=13)
    ax0.set_title('size distribution for Hi-Scale '+IOP,fontsize=15)
    fig.text(.65,.83,'Above PBL ('+str(n_above)+' legs)',fontsize=12)
    fig.text(.65,.43,'Below PBL ('+str(n_below)+' legs)',fontsize=12)
    # fig.text(.68,.83,'Above PBL ('+format(n_above/n_total*100,'.1f')+'%)',fontsize=12)
    # fig.text(.68,.43,'Below PBL ('+format(n_below/n_total*100,'.1f')+'%)',fontsize=12)
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    # plt.close()
    
    #%% plot percentile on sizes
    
    figname = figpath_aircraft_statistics+'SeparatePBLH_percentile_AerosolSize_HISCALE_'+IOP+'.png'
    print('plotting percentile figures to '+figname)
    
    # set position shift so that models and obs are not overlapped
    p_shift = np.arange(nmodels+1)
    p_shift = (p_shift - p_shift.mean())*0.2
    
    fig,(ax0,ax1) = plt.subplots(2,1,figsize=(8,6))
        
    ax0.boxplot(p2_above_obs,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(blen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax0.boxplot(p2_above_model[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(blen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    
    ax1.boxplot(p2_below_obs,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(blen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax1.boxplot(p2_below_model[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(blen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
        
    ax0.tick_params(color='k',labelsize=12)
    ax1.tick_params(color='k',labelsize=14)
    # ax0.set_xscale('log')
    # ax1.set_xscale('log')
    ax0.set_yscale('log')
    ax1.set_yscale('log')
    ax0.set_xlim(-1,blen)
    ax1.set_xlim(-1,blen)
    ax1.set_xlabel('Diameter (nm)',fontsize=14)
    ax0.set_ylabel('aerosol # (cm$^{-3}$)',fontsize=14)
    ax1.set_ylabel('aerosol # (cm$^{-3}$)',fontsize=14)
    ax0.set_title('percentile for Hi-Scale '+IOP,fontsize=15)
    fig.text(.66,.83,'Above PBL ('+str(n_above)+' legs)',fontsize=12)
    fig.text(.66,.43,'Below PBL ('+str(n_below)+' legs)',fontsize=12)
    
    xlabel=[str(binl[x])+'-'+str(binh[x]) for x in range(blen)]
    ax0.set_xticks(range(len(binm)))
    ax0.set_xticklabels(xlabel)
    ax1.set_xticks(range(len(binm)))
    ax1.set_xticklabels(xlabel)
    # ax0.set_yticks([1,3,5,7,9,11,12,13,14,15,16])
    # ax0.set_yticklabels(range(400,4100,400))
    ax0.set_ylim(1e-3,1e5)
    ax1.set_ylim(1e-3,1e5)
    
    # plot temporal lines for label
    ax1.plot([],c='k',label='Obs')
    for mm in range(nmodels):
        ax1.plot([],c=color_model[mm],label=Model_List[mm])
    ax1.legend(loc='lower left', shadow=False, fontsize='large')
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    # plt.close()
