"""
# plot percentile of Aerosol Sulfate and Organics with height
# for flight data in IOPs
# compare models and aircraft measurements
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.time_format_change import hhmmss2sec
from ..subroutines.read_aircraft import read_ams,read_iwg1
from ..subroutines.read_netcdf import read_merged_size,read_extractflight
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
        amspath = settings.get('amspath', None)
        iwgpath = settings.get('iwgpath', None)
    elif campaign in ['CSET', 'SOCRATES']:
        raise ValueError('CSET or SOCRATES do not have aerosol composition data')
    else:
        raise ValueError('check campaign name setting: '+campaign)

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
    elif campaign=='ACEENA':
        if IOP=='IOP1':
            lst=lst[0:20]
        elif IOP=='IOP2':
            lst=lst[20:]
        elif IOP[0:4]=='2017' or IOP[0:4]=='2018':
            a=lst[0].split('_'+campaign+'_')
            lst = glob.glob(a[0]+'*'+IOP+'*')
            lst.sort()
    else:
        raise ValueError('check campaign name setting: '+campaign)
    
    if len(lst)==0:
        raise ValueError('cannot find any file')
    
    #%% read all data
    
    height_all = []
    so4_o_all = []
    org_o_all = []
    so4_m_all = []
    org_m_all = []
    nmodels=len(Model_List)
    for mm in range(nmodels):
        so4_m_all.append([])
        org_m_all.append([])
        
    print('reading '+format(len(lst))+' files to calculate the statistics: ')
    
    for filename in lst:
        
        # get date info:        
        date=filename[-12:-3]
        if date[-1]=='a':
            flightidx=1
        else:
            flightidx=2
        print(date)
        
        #% read in flight information
        (time,size,cvi,timeunit,cunit,long_name)=read_merged_size(filename,'CVI_inlet')
        (time,size,cflag,timeunit,cunit,long_name)=read_merged_size(filename,'cld_flag')
        (time,size,height,timeunit,zunit,long_name)=read_merged_size(filename,'height')
        time=np.ma.compressed(time)
        
        #%% read T and P from iwg
        filename_i=glob.glob(iwgpath+'aaf.iwg*.'+date+'*txt')
        filename_i.sort()
        # read in data
        if len(filename_i)==1: 
            (iwg,iwgvars)=read_iwg1(filename_i[0])
            timelen = len(iwg)
            if np.logical_and(campaign=='ACEENA', date=='20180216a'):
                iwg.insert(1403,list(iwg[1403]))
                tstr=iwg[1403][1]
                tstr=tstr[0:-1]+str(int(tstr[-1])-1)
                iwg[1403][1]=tstr
                del iwg[-1]
            # get variables
            time_iwg=np.empty(timelen)
            T_iwg=np.empty(timelen)
            P_iwg=np.empty(timelen)
            for t in range(timelen):
                T_iwg[t]=float(iwg[t][20])+273.15
                P_iwg[t]=float(iwg[t][23])*100
                timestr=iwg[t][1].split(' ')
                time_iwg[t]=hhmmss2sec(timestr[1])
        else:
            raise ValueError('cannot find any file or find too many files: ' + filename_i)
        # remove cloud flag
        T_iwg=qc_mask_cloudflag(T_iwg,cflag)
        P_iwg=qc_mask_cloudflag(P_iwg,cflag)
        
        #%% read aerosol composition in AMS
        
        filename_ams=glob.glob(amspath+'*'+date[0:8]+'*')
        filename_ams.sort()
        
        if len(filename_ams)==1 or len(filename_ams)==2:
            (ams,amslist)=read_ams(filename_ams[flightidx-1])
            time_ams=ams[0,:]
            flag=ams[-1,:]
            orgaaf=ams[1,:]
            so4aaf=ams[5,:]
            orgaaf=qc_mask_qcflag(orgaaf,flag)
            so4aaf=qc_mask_qcflag(so4aaf,flag)
        elif len(filename_ams)==0:
            time_ams = time_iwg
            orgaaf = np.full(len(time_ams),np.nan)
            so4aaf = np.full(len(time_ams),np.nan)
        else:
            raise ValueError('find too many files: ' + filename_ams)
        
        # change values from standardize condition to ambient condition
        T_ams = np.interp(time_ams,time,T_iwg)
        P_ams = np.interp(time_ams,time,P_iwg)
        so4aaf = so4aaf * (296.15/T_ams) * (P_ams/101325.)
        orgaaf = orgaaf * (296.15/T_ams) * (P_ams/101325.)
        
        so4aaf = qc_remove_neg(so4aaf)
        orgaaf = qc_remove_neg(orgaaf)
    
        # exclude NaNs
        idx = np.logical_and(~np.isnan(so4aaf), ~np.isnan(orgaaf))
        so4_o_all.append(so4aaf[idx])
        org_o_all.append(orgaaf[idx])
        
        height2=np.interp(time_ams,time,height)
        height_all.append(height2[idx])
        
        # for interpolation of model results
        time_ams=time_ams[idx]
    
        #%% read in Models
        for mm in range(nmodels):
            filename_m = E3SM_aircraft_path+'Aircraft_vars_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
        
            (timem,heightm,soa_a1,timeunitm,soaunit,soaname)=read_extractflight(filename_m,'soa_a1')
            (timem,heightm,soa_a2,timeunitm,soaunit,soaname)=read_extractflight(filename_m,'soa_a2')
            (timem,heightm,soa_a3,timeunitm,soaunit,soaname)=read_extractflight(filename_m,'soa_a3')
            (timem,heightm,so4_a1,timeunitm,so4unit,so4name)=read_extractflight(filename_m,'so4_a1')
            (timem,heightm,so4_a2,timeunitm,so4unit,so4name)=read_extractflight(filename_m,'so4_a2')
            (timem,heightm,so4_a3,timeunitm,so4unit,so4name)=read_extractflight(filename_m,'so4_a3')
            (timem,heightm,pom_a1,timeunitm,pomunit,pomname)=read_extractflight(filename_m,'pom_a1')
            (timem,heightm,pom_a3,timeunitm,pomunit,pomname)=read_extractflight(filename_m,'pom_a3')
            (timem,heightm,pom_a4,timeunitm,pomunit,pomname)=read_extractflight(filename_m,'pom_a4')
            (timem,heightm,mom_a1,timeunitm,momunit,momname)=read_extractflight(filename_m,'mom_a1')
            (timem,heightm,mom_a2,timeunitm,momunit,momname)=read_extractflight(filename_m,'mom_a2')
            (timem,heightm,mom_a3,timeunitm,momunit,momname)=read_extractflight(filename_m,'mom_a3')
            (timem,heightm,mom_a4,timeunitm,momunit,momname)=read_extractflight(filename_m,'mom_a4')
            
            # add nucleation mode if available
            try:
                (timem,heightm,soa_a5,timeunitm,soaunit,soaname)=read_extractflight(filename_m,'soa_a5')
                model_org = soa_a1+soa_a2+soa_a3+soa_a5 + pom_a1+pom_a3+pom_a4 + mom_a1+mom_a2+mom_a3+mom_a4
            except:
                model_org = soa_a1+soa_a2+soa_a3 + pom_a1+pom_a3+pom_a4 + mom_a1+mom_a2+mom_a3+mom_a4
            try:
                (timem,heightm,so4_a5,timeunitm,so4unit,so4name)=read_extractflight(filename_m,'so4_a5')
                model_so4 = so4_a1+so4_a2+so4_a3+so4_a5
            except:
                model_so4 = so4_a1+so4_a2+so4_a3
            
            # change E3SM unit from kg/kg to ug/m3 
            rho = P_iwg/T_iwg/287.06
            model_so4=model_so4*1e9*rho
            model_org=model_org*1e9*rho
            
            # interpolate into observational time
            so4_m_all[mm].append(np.interp(time_ams,timem,model_so4)) 
            org_m_all[mm].append(np.interp(time_ams,timem,model_org)) 
        
    #%% calculate percentiles for each height bin
    
    so4_o_z = list()
    org_o_z = list()
    so4_m_z = []
    org_m_z = []
    for mm in range(nmodels):
        so4_m_z.append([])
        org_m_z.append([])
    for zz in range(zlen):
        so4_o_z.append(np.empty(0))
        org_o_z.append(np.empty(0))
        for mm in range(nmodels):
            so4_m_z[mm].append(np.empty(0))
            org_m_z[mm].append(np.empty(0))
        
    ndays=len(height_all)
    for dd in range(ndays):
        height = height_all[dd]
        so4_o = so4_o_all[dd]
        org_o = org_o_all[dd]
        for zz in range(zlen):
            idx = np.logical_and(height>=zmin[zz], height<zmax[zz])
            so4_o_z[zz]=np.append(so4_o_z[zz],so4_o[idx])
            org_o_z[zz]=np.append(org_o_z[zz],org_o[idx])
            for mm in range(nmodels):
                so4_m = so4_m_all[mm][dd]
                org_m = org_m_all[mm][dd]
                so4_m_z[mm][zz]=np.append(so4_m_z[mm][zz],so4_m[idx])
                org_m_z[mm][zz]=np.append(org_m_z[mm][zz],org_m[idx])
            
    
    #%% make plot
    # set position shift so that models and obs are not overlapped
    p_shift = np.arange(nmodels+1)
    p_shift = (p_shift - p_shift.mean())*0.2
    
        
    figname = figpath_aircraft_statistics+'percentile_height_AerosolComposition_'+campaign+'_'+IOP+'.png'
    print('plotting figures to '+figname)
    
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(8,8))   # figsize in inches
    # plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
        
    ax1.boxplot(so4_o_z,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(zlen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=False, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax1.boxplot(so4_m_z[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(zlen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=False, patch_artist=True)    # need patch_artist to fill color in box
    ax1.tick_params(color='k',labelsize=16)
    # ax1.set_xscale('log')
    ax1.set_ylim(-1,zlen)
    ax1.set_yticks(range(zlen))
    ax1.set_yticklabels(z)
    # ax1.set_yticks([1,3,5,7,9,11,12,13,14,15,16])
    # ax1.set_yticklabels(range(400,4100,400))
    # plot temporal lines for label
    ax1.plot([],c='k',label='Obs')
    for mm in range(nmodels):
        ax1.plot([],c=color_model[mm],label=Model_List[mm])
    ax1.legend(loc='upper right', fontsize='x-large')
        
    ax2.boxplot(org_o_z,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(zlen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=False, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax2.boxplot(org_m_z[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(zlen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=False, patch_artist=True)    # need patch_artist to fill color in box
    ax2.tick_params(color='k',labelsize=16)
    # ax2.set_xscale('log')
    ax2.set_ylim(-1,zlen)
    ax2.set_yticks(range(zlen))
    ax2.set_yticklabels([])
    # ax1.set_yticks(np.arange(0,20,2))
    # ax1.set_yticklabels(range(400,4100,400))
    # plot temporal lines for label
    ax2.plot([],c='k',label='Obs')
    for mm in range(nmodels):
        ax2.plot([],c=color_model[mm],label=Model_List[mm])
    ax2.legend(loc='upper right', fontsize='x-large')
        
    # set xlimit consistent in subplots
    # xlim1 = ax1.get_xlim()
    # xlim2 = ax2.get_xlim()
    # ax1.set_xlim([min(xlim1[0],xlim2[0]), max(xlim1[1],xlim2[1])])
    # ax2.set_xlim([min(xlim1[0],xlim2[0]), max(xlim1[1],xlim2[1])])
    
    ax1.set_ylabel('Height (m MSL)',fontsize=16)
    fig.text(0.46,0.06, '$\mu$g/m$^3$', fontsize=16)
    ax1.set_title('Sulfate',fontsize=16)
    ax2.set_title('Organic',fontsize=16)
    fig.text(0.48,0.92, IOP, fontsize=18)
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    # plt.close()
    