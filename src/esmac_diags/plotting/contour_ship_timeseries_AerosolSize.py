# plot timeseries of surface aerosol size distribution along each ship leg

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_ARMdata import read_uhsas
from ..subroutines.read_netcdf import read_E3SM
from ..subroutines.time_format_change import cday2mmdd
from ..subroutines.specific_data_treatment import  avg_time_2d
from ..subroutines.quality_control import qc_remove_neg

def run_plot(settings):
    #%% variables from settings

    campaign = settings['campaign']
    Model_List = settings['Model_List']
    shipuhsaspath = settings['shipuhsaspath']
    E3SM_ship_path = settings['E3SM_ship_path']
    figpath_ship_timeseries = settings['figpath_ship_timeseries']

    #%% other settings
    
    if not os.path.exists(figpath_ship_timeseries):
        os.makedirs(figpath_ship_timeseries)
    
    lst = glob.glob(E3SM_ship_path+'Ship_CNsize_'+campaign+'_'+Model_List[0]+'_shipleg*.nc')
    lst.sort()
    
    for ll in range(len(lst)):
        if campaign=='MAGIC':
            legnum=lst[ll][-5:-3]
        elif campaign=='MARCUS':
            legnum=lst[ll][-4]
        
    
        #%% read in model
        nmodels=len(Model_List)
        datam = list()
        for mm in range(nmodels):
            filenamem = E3SM_ship_path+'Ship_CNsize_'+campaign+'_'+Model_List[mm]+'_shipleg'+legnum+'.nc'
        
            (timem0,data,timeunitm,datamunit,datamlongname)=read_E3SM(filenamem,'NCNall')
        
            # if time expands two years, add 365 days to the second year
            if timem0[0]>timem0[-1]:
                timem0[timem0<=timem0[-1]]=timem0[timem0<=timem0[-1]]+365
            
            # average in time for quicker plot
            timem=np.arange(timem0[0]-0.1,timem0[-1]+0.1,1/24.)
            data2 = avg_time_2d(timem0,data.T,timem)
            data2 = data2.T
            
            # change to dN/dlnDp
            for bb in range(3000):
                dlnDp=np.log((bb+2)/(bb+1))
                data2[bb,:]=data2[bb,:]/dlnDp
            datam.append(data2*1e-6)    # change unit from 1/m3 to 1/cm3
            
        year0 = str(int(timeunitm.split()[2][0:4])+1)
        
        #%% read in observations
        # find the days related to the ship leg
        day = [int(a) for a in timem]
        day = list(set(day))
        day.sort()
        
        nbins = 99 # for UHSAS at MAGIC
        t_uh=np.empty(0)
        uhsasall=np.empty((0,nbins))
        for dd in day:
            if campaign=='MAGIC':
                if int(legnum)<=9:
                    if dd<=365:  # year 2012
                        filenameo = glob.glob(shipuhsaspath+'magaosuhsasM1.a1.2012'+cday2mmdd(dd,calendar='noleap')+'.*.cdf')
                    else:
                        filenameo = glob.glob(shipuhsaspath+'magaosuhsasM1.a1.2013'+cday2mmdd(dd-365,calendar='noleap')+'.*.cdf')
                else:
                    filenameo = glob.glob(shipuhsaspath+'magaosuhsasM1.a1.2013'+cday2mmdd(dd,calendar='noleap')+'.*.cdf')
            elif campaign=='MARCUS':
                if int(legnum)<=2:
                    if dd<=365:  # year 2012
                        filenameo = glob.glob(shipuhsaspath+'maraosuhsasM1.a1.2017'+cday2mmdd(dd,calendar='noleap')+'.*')
                    else:
                        filenameo = glob.glob(shipuhsaspath+'maraosuhsasM1.a1.2018'+cday2mmdd(dd-365,calendar='noleap')+'.*')
                else:
                    filenameo = glob.glob(shipuhsaspath+'maraosuhsasM1.a1.2018'+cday2mmdd(dd,calendar='noleap')+'.*')
            
            if len(filenameo)==0:
                continue  # some days may be missing
            if len(filenameo)>1:
                raise ValueError('find too many files: ' + filenameo)
            
            
            (time,dmin,dmax,uhsas,timeunit,uhunit,uhlongname)=read_uhsas(filenameo[0])
            
            uhsas=np.ma.filled(uhsas)
            uhsas=qc_remove_neg(uhsas)
            
            # average in time for quicker plot
            time2=np.arange(1800,86400,3600)
            data2 = avg_time_2d(time,uhsas,time2)
            uhsasall=np.vstack((uhsasall, data2))
            t_uh = np.hstack((t_uh,time2/86400+dd))
            
        # if no obs available, fill one data with NaN
        if len(t_uh)==0:
            t_uh=[timem[0],timem[1]]
            uhsasall=np.full((2,nbins),np.nan)
            
        # if time expands two years, add 365 days to the second year
        if t_uh[0]>t_uh[-1]:
            t_uh[t_uh<=t_uh[-1]]=t_uh[t_uh<=t_uh[-1]]+365
            
        size_u = (dmin+dmax)/2
        dsize_u = dmax-dmin
        
        uhsasall=qc_remove_neg(uhsasall)
        
        # change to dN/dlnDp
        dlnDp_u=np.empty(nbins)
        for bb in range(len(size_u)):
            dlnDp_u[bb]=np.log(dmax[bb]/dmin[bb])
            uhsasall[:,bb]=uhsasall[:,bb]/dlnDp_u[bb]
        
        #%% make plot
            
        figname = figpath_ship_timeseries+'timeseries_AerosolSize_'+campaign+'_ship'+legnum+'.png'
        print('plotting figures to '+figname)
        
        #fig = plt.figure()
        fig,ax = plt.subplots(nmodels+1,1,figsize=(8,2*(nmodels+1)))   # figsize in inches
        plt.tight_layout(h_pad=1.1)   #pad=0.4, w_pad=0.5, h_pad=1.0
        plt.subplots_adjust(right=0.9,bottom=0.1)
        
        leveltick=[0.1,1,10,100,1000,10000,100000]
        levellist=np.arange(np.log(leveltick[0]),12,.5)
        
        uhsasall[uhsasall<0.01]=0.01
        h1 = ax[0].contourf(t_uh,size_u,np.log(uhsasall.T),levellist,cmap=plt.get_cmap('jet'))
        
        size_m=np.arange(1,3001)
        h2=[]
        for mm in range(nmodels):
            data = datam[mm]
            data[data<0.01]=0.01
            h_m = ax[mm+1].contourf(timem,size_m,np.log(data),levellist,cmap=plt.get_cmap('jet'))
            h2.append(h_m)
    
        # colorbar
        cax = plt.axes([0.92, 0.2, 0.02, 0.6])
        cbar=fig.colorbar(h2[0], cax=cax, ticks=np.log(leveltick))
        cbar.ax.set_yticklabels(leveltick, fontsize=14)
        
        # set axis
        for ii in range(nmodels+1):
            ax[ii].set_xlim(timem[0],timem[-1])
            ax[ii].set_yscale('log')
            ax[ii].set_ylim(5, 3000)
            ax[ii].set_yticks([10,100,1000])
            ax[ii].tick_params(color='k',labelsize=14)
            if ii==0:
                ax[ii].text(0.01, 0.94, 'OBS', fontsize=14,transform=ax[ii].transAxes, verticalalignment='top')
            else:
                ax[ii].text(0.01, 0.94, Model_List[ii-1], fontsize=14,transform=ax[ii].transAxes, verticalalignment='top')
            
        ax[1].set_ylabel('Diameter (nm)',fontsize=14)
        ax[0].set_title('Size Distribution (#/dlnDp, cm-3)',fontsize=15)
        ax[nmodels].set_xlabel('Calendar Day in '+year0,fontsize=14)
        
        fig.text(.08, .97,'ship leg '+legnum, fontsize=12)
        
        fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        