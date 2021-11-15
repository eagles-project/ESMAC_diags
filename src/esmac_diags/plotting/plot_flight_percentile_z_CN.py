"""
# plot percentile of aerosol number concentration (CN) with height
# for flight data in IOPs
# compare models and CPC measurements
"""

import glob
import os
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_aircraft import read_cpc, read_RF_NCAR
from ..subroutines.read_netcdf import read_merged_size,read_extractflight
from ..subroutines.quality_control import qc_mask_cloudflag,qc_remove_neg,qc_cpc_air

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
        cpcpath = settings.get('cpcpath', None)
    elif campaign in ['CSET', 'SOCRATES']:
        RFpath = settings.get('RFpath', None)
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
    cpc10_o = []
    cpc3_o = []
    uhsas100_o = []
    cpc100_m = []
    cpc10_m = []
    cpc3_m = []
    nmodels=len(Model_List)
    for mm in range(nmodels):
        cpc100_m.append([])
        cpc10_m.append([])
        cpc3_m.append([])
        
    print('reading '+format(len(alldates))+' files to calculate the statistics: ')
    
    for date in alldates:
        print(date)
        
        #%% read in Models
        for mm in range(nmodels):
            filename_m = E3SM_aircraft_path+'Aircraft_CNsize_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
        
            (timem,heightm,cpc_m,timeunitm,ncn_unit,ncn_longname)=read_extractflight(filename_m,'NCN')
            (timem,heightm,cpcu_m,timeunitm,ncnu_unit,ncnu_longname)=read_extractflight(filename_m,'NUCN')
            (timem,heightm,ncnall,timeunitm,ncnall_unit,ncnall_longname)=read_extractflight(filename_m,'NCNall')
            
            cpc100_m[mm].append(np.sum(ncnall[100:,:],0)*1e-6) # #/m3 to #/cm3
            cpc10_m[mm].append(cpc_m*1e-6)    # #/m3 to #/cm3
            cpc3_m[mm].append(cpcu_m*1e-6)    # #/m3 to #/cm3
        
        height_all.append(heightm)
        
        #%% read in CPC measurements
        if campaign in ['HISCALE', 'ACEENA']:
            if date[-1]=='a':
                flightidx=1
            else:
                flightidx=2
            if campaign=='HISCALE':
                filename_c=glob.glob(cpcpath+'CPC_G1_'+date[0:8]+'*R2_HiScale001s.ict.txt')
            elif campaign=='ACEENA':
                filename_c=glob.glob(cpcpath+'CPC_G1_'+date[0:8]+'*R2_ACEENA001s.ict')    
            filename_c.sort()
            # read in data
            if len(filename_c)==1 or len(filename_c)==2: # some days have two flights
                (cpc,cpclist)=read_cpc(filename_c[flightidx-1])
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
            
            # cloud flag
            if campaign=='HISCALE':
                filename = merged_size_path+'merged_bin_fims_pcasp_'+campaign+'_'+date+'.nc'
            elif campaign=='ACEENA':
                filename = merged_size_path+'merged_bin_fims_pcasp_opc_'+campaign+'_'+date+'.nc'
            (time,size,cflag,timeunit,cunit,long_name)=read_merged_size(filename,'cld_flag')
            cpc3=qc_mask_cloudflag(cpc3,cflag)
            cpc10=qc_mask_cloudflag(cpc10,cflag)
            
            # some quality checks
            (cpc3,cpc10) = qc_cpc_air(cpc3,cpc10)
            
            cpc10_o.append(cpc10)
            cpc3_o.append(cpc3)
        
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
                raise ValueError('cannot find any file or find too many files: '+filename)
            
            # some quality checks
            uhsas100=qc_remove_neg(uhsas100)
            cpc10=qc_remove_neg(cpc10)
            
            cpc10_o.append(cpc10)
            uhsas100_o.append(uhsas100)
            
    #%% calculate percentiles for each height bin
    
    uhsas100_o_z = list()
    cpc10_o_z = list()
    cpc3_o_z = list()
    cpc100_m_z = []
    cpc10_m_z = []
    cpc3_m_z = []
    nmodels=len(Model_List)
    for mm in range(nmodels):
        cpc100_m_z.append([])
        cpc10_m_z.append([])
        cpc3_m_z.append([])
    for zz in range(zlen):
        uhsas100_o_z.append(np.empty(0))
        cpc10_o_z.append(np.empty(0))
        cpc3_o_z.append(np.empty(0))
        for mm in range(nmodels):
            cpc100_m_z[mm].append(np.empty(0))
            cpc10_m_z[mm].append(np.empty(0))
            cpc3_m_z[mm].append(np.empty(0))
        
    ndays=len(height_all)
    for dd in range(ndays):
        height = height_all[dd]
        cpc10 = cpc10_o[dd]
        if campaign in ['HISCALE', 'ACEENA']:
            cpc3 = cpc3_o[dd]
        elif campaign in ['CSET', 'SOCRATES']:
            uhsas100 = uhsas100_o[dd]
        for zz in range(zlen):
            idx = np.logical_and(height>=zmin[zz], height<zmax[zz])
            cpc10_o_z[zz]=np.append(cpc10_o_z[zz],cpc10[np.logical_and(idx,~np.isnan(cpc10))])
            for mm in range(nmodels):
                model10 = cpc10_m[mm][dd]
                cpc10_m_z[mm][zz]=np.append(cpc10_m_z[mm][zz],model10[idx])
            if campaign in ['HISCALE', 'ACEENA']:
                cpc3_o_z[zz]=np.append(cpc3_o_z[zz],cpc3[np.logical_and(idx,~np.isnan(cpc3))])
                for mm in range(nmodels):
                    model3 = cpc3_m[mm][dd]
                    cpc3_m_z[mm][zz]=np.append(cpc3_m_z[mm][zz],model3[idx])
            elif campaign in ['CSET', 'SOCRATES']:
                uhsas100_o_z[zz]=np.append(uhsas100_o_z[zz],uhsas100[np.logical_and(idx,~np.isnan(uhsas100))])
                for mm in range(nmodels):
                    model100 = cpc100_m[mm][dd]
                    cpc100_m_z[mm][zz]=np.append(cpc100_m_z[mm][zz],model100[idx])
    
    #%% make plot
    # set position shift so that models and obs are not overlapped
    p_shift = np.arange(nmodels+1)
    p_shift = (p_shift - p_shift.mean())*0.2
    
    if campaign in ['HISCALE', 'ACEENA']:
        figname = figpath_aircraft_statistics+'percentile_height_CN_'+campaign+'_'+IOP+'.png'
    else:
        figname = figpath_aircraft_statistics+'percentile_height_CN_'+campaign+'.png'
    print('plotting figures to '+figname)
    
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(8,8))   # figsize in inches
    # plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
        
    ax1.boxplot(cpc10_o_z,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(zlen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=False, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax1.boxplot(cpc10_m_z[mm],whis=(5,95),showmeans=False,showfliers=False,
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
    ax1.plot([],c='k',label='CPC(>10nm)')
    for mm in range(nmodels):
        ax1.plot([],c=color_model[mm],label=Model_List[mm])
    ax1.legend(loc='upper right', fontsize='x-large')
        
    if campaign in ['HISCALE', 'ACEENA']:
        ax2.boxplot(cpc3_o_z,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(zlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=False, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax2.boxplot(cpc3_m_z[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(zlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=False, patch_artist=True)    # need patch_artist to fill color in box
        ax2.plot([],c='k',label='CPC(>3nm)')
    elif campaign in ['CSET', 'SOCRATES']:
        ax2.boxplot(uhsas100_o_z,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(zlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=False, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax2.boxplot(cpc100_m_z[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(zlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=False, patch_artist=True)    # need patch_artist to fill color in box
        ax2.plot([],c='k',label='UHSAS(>100nm)')
        
    ax2.tick_params(color='k',labelsize=16)
    ax2.set_ylim(-1,zlen)
    ax2.set_yticks(range(zlen))
    ax2.set_yticklabels([])
    # ax1.set_yticks(np.arange(0,20,2))
    # ax1.set_yticklabels(range(400,4100,400))
    # plot temporal lines for label
    for mm in range(nmodels):
        ax2.plot([],c=color_model[mm],label=Model_List[mm])
    ax2.legend(loc='upper right', fontsize='x-large')
        
    # set xlimit consistent in subplots
    # xlim1 = ax1.get_xlim()
    # xlim2 = ax2.get_xlim()
    # ax1.set_xlim([min(xlim1[0],xlim2[0]), max(xlim1[1],xlim2[1])])
    # ax2.set_xlim([min(xlim1[0],xlim2[0]), max(xlim1[1],xlim2[1])])
    if campaign=='HISCALE':
        ax1.set_xlim(-200,15000)
        ax2.set_xlim(-200,15000)
    elif campaign=='ACEENA':
        ax1.set_xlim(0,3000)
        ax2.set_xlim(0,3000)
    elif campaign=='CSET':
        ax1.set_xscale('log')
        ax2.set_xscale('log')
        ax1.set_xlim(1,1e4)
        ax2.set_xlim(1,1e4)
    elif campaign=='SOCRATES':
        ax1.set_xscale('log')
        ax2.set_xscale('log')
        ax1.set_xlim(1,1e4)
        ax2.set_xlim(1,1e4)
    
    ax1.set_ylabel('Height (m MSL)',fontsize=16)
    fig.text(0.4,0.06, 'Aerosol number (cm$^{-3}$)', fontsize=16)
    if campaign in ['HISCALE', 'ACEENA']:
        fig.text(0.48,0.9, IOP, fontsize=18)
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    # plt.close()
    
    #%% plot sample numbers
    # num_sample = [len(a) for a in cpc10_o_z]
    # fig,ax = plt.subplots(figsize=(2,8))
    # ax.plot(num_sample,z,color='k',linewidth=1,linestyle='-')
    # ax.set_xscale('log')
    # ax.set_xlabel('Sample Number (#)',fontsize=16)
    # ax.set_ylabel('Height (m MSL)',fontsize=16)
    # ax.set_xticks([10,1e3,1e5])
    # ax.set_yticks(z)
    # ax.tick_params(color='k',labelsize=16)
    
    # if campaign in ['HISCALE', 'ACEENA']:
    #     figname = figpath_aircraft_statistics+'samplenumber_'+campaign+'_'+IOP+'.png'
    # else:
    #     figname = figpath_aircraft_statistics+'samplenumber_'+campaign+'.png'
    # fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)