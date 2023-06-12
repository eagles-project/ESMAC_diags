"""
script to generate all plots for HISCALE surface data

"""
import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import esmac_diags.plotting.plot_esmac_diags as plot
import esmac_diags.plotting.calc_statistics as calc
import matplotlib.dates as mdates

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings
# set site name and datapath

# set site name.
site = 'HISCALE'

prep_model_path = '../prep_data/'+site+'/model/'
prep_sfc_path = '../prep_data/'+site+'/surface/'
prep_sat_path = '../prep_data/'+site+'/satellite/'

# time_hiscale = pd.date_range(start='2016-04-25', end='2016-05-21', freq="3600s")
# IOP = 'IOP2'
for IOP in ['IOP1','IOP2']:

        
        if IOP=='IOP1':
            time_hiscale = pd.date_range(start='2016-04-25', end='2016-05-21', freq="3600s")
        elif IOP=='IOP2':
            time_hiscale = pd.date_range(start='2016-08-28', end='2016-09-23', freq="3600s") 
        
                    
        # path of output figures
        figpath= '../figures/RRM/'+site+'/sfc_toa/'
        if not os.path.exists(figpath):
            os.makedirs(figpath)
        
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # read in data
        obsdata = xr.open_dataset(prep_sfc_path + 'sfc_ACSM_'+site+'.nc')
        org = obsdata['org'].load()
        so4 = obsdata['so4'].load()
        nh4 = obsdata['nh4'].load()
        no3 = obsdata['no3'].load()
        chl = obsdata['chl'].load()
        obsdata.close()
        org_hiscale = org.sel(time=time_hiscale)
        so4_hiscale = so4.sel(time=time_hiscale)
        nh4_hiscale = nh4.sel(time=time_hiscale)
        no3_hiscale = no3.sel(time=time_hiscale)
        chl_hiscale = chl.sel(time=time_hiscale)
        
        obsdata = xr.open_dataset(prep_sfc_path + 'totcld_'+site+'.nc')
        cld_arscl = obsdata['tot_cld_arscl'].load()
        cld_tsi = obsdata['tot_cld_tsi'].load()
        cld_visst = obsdata['tot_cld_visst'].load()
        obsdata.close()
        cld_arscl_hiscale = cld_arscl.sel(time=time_hiscale)
        cld_tsi_hiscale = cld_tsi.sel(time=time_hiscale)
        cld_visst_hiscale = cld_visst.sel(time=time_hiscale)
        obsdata = xr.open_dataset(prep_sfc_path + 'cloud_2d_'+site+'.nc')
        height_o = obsdata['height'].load()
        cloud_2d = obsdata['cloud'].load()
        obsdata.close()
        cloud_2d_hiscale = cloud_2d.sel(time=time_hiscale)
        obsdata = xr.open_mfdataset(sorted(glob.glob(prep_sfc_path + 'sfc_CCN_'+site+'_*.nc')))
        ccn2 = obsdata['CCN2'].load()
        obsdata.close()
        ccn2_hiscale = ccn2.sel(time=time_hiscale)
        obsdata = xr.open_dataset(prep_sfc_path + 'sfc_CPC_'+site+'.nc')
        cpc10 = obsdata['cpc10'].load()
        cpc3 = obsdata['cpc3'].load()
        obsdata.close()
        cpc3_hiscale = cpc3.sel(time=time_hiscale)
        cpc10_hiscale = cpc10.sel(time=time_hiscale)
        obsdata = xr.open_dataset(prep_sfc_path + 'cod_'+site+'.nc')
        cod = obsdata['cod'].load()
        obsdata.close()
        cod_hiscale = cod.sel(time=time_hiscale)
        obsdata = xr.open_dataset(prep_sfc_path + 'LWP_'+site+'.nc')
        lwp_armbe = obsdata['lwp_armbe'].load()
        lwp_mfrsr = obsdata['lwp_mfrsr'].load()
        obsdata.close()
        lwp_armbe_hiscale = lwp_armbe.sel(time=time_hiscale)
        lwp_mfrsr_hiscale = lwp_mfrsr.sel(time=time_hiscale)
        obsdata = xr.open_dataset(prep_sfc_path + 'Ndrop_'+site+'.nc')
        ndrop = obsdata['cdnc'].load()
        obsdata.close()
        ndrop_hiscale = ndrop.sel(time=time_hiscale)
        obsdata = xr.open_dataset(prep_sfc_path + 'reff_'+site+'.nc')
        reff = obsdata['reff'].load()
        obsdata.close()
        reff_hiscale = reff.sel(time=time_hiscale)
        obsdata = xr.open_dataset(prep_sfc_path + 'precip_'+site+'.nc')
        precip = obsdata['precip'].load()
        obsdata.close()
        precip_hiscale = precip.sel(time=time_hiscale)
        
        obsdata = xr.open_dataset(prep_sfc_path + 'sfc_radiation_'+site+'.nc')
        lwdnsfc = obsdata['lwdn'].load()
        swdnsfc = obsdata['swdn'].load()
        lwupsfc = obsdata['lwup'].load()
        swupsfc = obsdata['swup'].load()
        obsdata.close()
        lwnetsfc = lwupsfc - lwdnsfc
        swnetsfc = swdnsfc - swupsfc
        lwnetsfc_hiscale = lwnetsfc.sel(time=time_hiscale)
        swnetsfc_hiscale = swnetsfc.sel(time=time_hiscale)
        
        obsdata = xr.open_dataset(prep_sfc_path + 'sfc_UHSAS_'+site+'.nc')
        size_uhsas = obsdata['size'].load()
        dmin_hiscale = obsdata['size_low'].load()
        dmax_hiscale = obsdata['size_high'].load()
        uhsas100 = obsdata['uhsas100'].load()
        uhsas_all = obsdata['uhsas_all'].load()
        obsdata.close()
        uhsas100_hiscale = uhsas100.sel(time=time_hiscale)
        uhsasall_hiscale = uhsas_all.sel(time=time_hiscale)
        dlogDp_uhsas = np.mean(np.log10(dmax_hiscale/dmin_hiscale))
        
        obsdata = xr.open_dataset(prep_sfc_path + 'sfc_SMPS_'+site+'_IOP1.nc')
        time1 = obsdata['time'].load()
        smps100_1 = obsdata['smps100_dlogDp'].load()
        smpsall_1 = obsdata['dN_dlogDp'].load()
        size1 = obsdata['size'].load()
        obsdata.close()
        obsdata = xr.open_dataset(prep_sfc_path + 'sfc_SMPS_'+site+'_IOP2.nc')
        time2 = obsdata['time'].load()
        smps100_2 = obsdata['smps100_dlogDp'].load()
        smpsall_2 = obsdata['dN_dlogDp'].load()
        size2 = obsdata['size'].load()
        obsdata.close()
        time_smps = xr.concat((time1,time2),dim='time')
        size_smps = size1
        smps100 = xr.concat((smps100_1,smps100_2),dim='time')
        smps_all = xr.concat((smpsall_1,smpsall_2.interp(size=smpsall_1['size'])),dim='time')
        # SMPS data is already dN/dlogDp, total number concentration must multiply by dlogDp
        dlogDp_smps = np.mean(np.log10(size_smps[1:].data/size_smps[0:-1].data))
        smps100 = smps100 * dlogDp_smps
        smps_all = smps_all * dlogDp_smps
        smps100_hiscale = smps100.sel(time=time_hiscale)
        smpsall_hiscale = smps_all.sel(time=time_hiscale)
        
        obsdata = xr.open_mfdataset(prep_sfc_path + 'cloudheight_ARSCL_'+site+'.nc')
        cth = obsdata['cth'].load()
        cbh = obsdata['cbh'].load()
        cths = obsdata['cths'].load()
        obsdata.close()
        cth_hiscale = cth.sel(time=time_hiscale)
        cbh_hiscale = cbh.sel(time=time_hiscale)
        
        # satellite data
        obsdata = xr.open_dataset(prep_sat_path + 'albedo_VISSTgrid_'+site+'.nc')
        albedo_sat = obsdata['albedo'].load()
        obsdata.close()
        albedo_hiscale = albedo_sat.sel(time=time_hiscale)
        obsdata = xr.open_dataset(prep_sat_path + 'cloudfraction_VISSTgrid_'+site+'.nc')
        cfall_sat = obsdata['cldtot'].load()
        cflow_sat = obsdata['cldlow'].load()
        obsdata.close()
        cfall_sat_hiscale = cfall_sat.sel(time=time_hiscale)
        obsdata = xr.open_dataset(prep_sat_path + 'cloudtop_VISSTgrid_'+site+'.nc')
        ctt_sat = obsdata['ctt'].load()
        cth_sat = obsdata['cth'].load()
        obsdata.close()
        ctt_sat_hiscale = ctt_sat.sel(time=time_hiscale)
        cth_sat_hiscale = cth_sat.sel(time=time_hiscale)
        obsdata = xr.open_dataset(prep_sat_path + 'cod_VISSTgrid_'+site+'.nc')
        cod_sat = obsdata['cod'].load()
        obsdata.close()
        cod_sat_hiscale = cod_sat.sel(time=time_hiscale)
        obsdata = xr.open_dataset(prep_sat_path + 'Reff_VISSTgrid_'+site+'.nc')
        reff_sat = obsdata['reff'].load()
        obsdata.close()
        reff_sat_hiscale = reff_sat.sel(time=time_hiscale)
        obsdata = xr.open_dataset(prep_sat_path + 'LWP_VISSTgrid_'+site+'.nc')
        lwp_sat = obsdata['lwp'].load()
        obsdata.close()
        lwp_sat_hiscale = lwp_sat.sel(time=time_hiscale)
        obsdata = xr.open_dataset(prep_sat_path + 'Nd_VISSTgrid_'+site+'.nc')
        nd_sat = obsdata['Nd'].load()
        obsdata.close()
        nd_sat_hiscale = nd_sat.sel(time=time_hiscale)
        obsdata = xr.open_dataset(prep_sat_path + 'lwflx_VISSTgrid_'+site+'.nc')
        lwnettoa = obsdata['lwnettoa'].load()
        obsdata.close()
        lwnettoa_hiscale = lwnettoa.sel(time=time_hiscale)
        obsdata = xr.open_dataset(prep_sat_path + 'swflx_VISSTgrid_'+site+'.nc')
        swnettoa = obsdata['swnettoa'].load()
        obsdata.close()
        swnettoa_hiscale = swnettoa.sel(time=time_hiscale)
        
        # E3SM data
        filename = prep_model_path + 'hiscale_RRMne30_'+IOP+'_sfc.nc'
        modeldata = xr.open_dataset(filename)
        time_m = modeldata['time'].load()
        bc_m = modeldata['bc'].load()
        dst_m = modeldata['dst'].load()
        mom_m = modeldata['mom'].load()
        pom_m = modeldata['pom'].load()
        ncl_m = modeldata['ncl'].load()
        so4_m = modeldata['so4'].load()
        soa_m = modeldata['soa'].load()
        ccn2_m = modeldata['CCN4'].load()
        ncn3_m = modeldata['NCN3'].load()
        ncn10_m = modeldata['NCN10'].load()
        ncn100_m = modeldata['NCN100'].load()
        CNsize_m = modeldata['NCNall'].load()
        cod_m = modeldata['cod'].load()
        reff_m = modeldata['reff'].load()
        lwp_m = modeldata['TGCLDLWP'].load()
        nd_m = modeldata['Nd_mean'].load()
        precip_m = modeldata['PRECT'].load()
        cld_m = modeldata['CLDTOT'].load()
        cbh_m = modeldata['cbh'].load()
        cth_m = modeldata['cth'].load()
        Hcld_m = modeldata['clddepth'].load()
        cldlow_m = modeldata['CLDLOW'].load()
        cldmid_m = modeldata['CLDMED'].load()
        cldhgh_m = modeldata['CLDHGH'].load()
        lwdnsfc_m = modeldata['FLDS'].load()
        lwnetsfc_m = modeldata['FLNS'].load()
        lwnettoa_m = modeldata['FLNT'].load()
        lwuptoa_m = modeldata['FLUT'].load()
        swdnsfc_m = modeldata['FSDS'].load()
        swnetsfc_m = modeldata['FSNS'].load()
        swdntoa_m = modeldata['SOLIN'].load()
        swnettoa_m = modeldata['FSNT'].load()
        swuptoa_m = modeldata['FSUTOA'].load()
        modeldata.close()
        lwupsfc_m = lwnetsfc_m + lwdnsfc_m
        swupsfc_m = swdnsfc_m - swnetsfc_m
        albedo_m = swuptoa_m/swdntoa_m*100
        org_m = pom_m + mom_m + soa_m
        bc_m_hiscale = bc_m.sel(time=time_hiscale)
        dst_m_hiscale = dst_m.sel(time=time_hiscale)
        org_m_hiscale = org_m.sel(time=time_hiscale)
        so4_m_hiscale = so4_m.sel(time=time_hiscale)
        ncl_m_hiscale = ncl_m.sel(time=time_hiscale)
        ccn2_m_hiscale = ccn2_m.sel(time=time_hiscale)
        ncn3_m_hiscale = ncn3_m.sel(time=time_hiscale)
        ncn10_m_hiscale = ncn10_m.sel(time=time_hiscale)
        ncn100_m_hiscale = ncn100_m.sel(time=time_hiscale)
        CNsize_m_hiscale = CNsize_m.sel(time=time_hiscale)
        cod_m_hiscale = cod_m.sel(time=time_hiscale)
        reff_m_hiscale = reff_m.sel(time=time_hiscale)
        lwp_m_hiscale = lwp_m.sel(time=time_hiscale)
        nd_m_hiscale = nd_m.sel(time=time_hiscale)
        precip_m_hiscale = precip_m.sel(time=time_hiscale)
        cld_m_hiscale = cld_m.sel(time=time_hiscale)
        cbh_m_hiscale = cbh_m.sel(time=time_hiscale)
        cth_m_hiscale = cth_m.sel(time=time_hiscale)
        Hcld_m_hiscale = Hcld_m.sel(time=time_hiscale)
        cldlow_m_hiscale = cldlow_m.sel(time=time_hiscale)
        cldmid_m_hiscale = cldmid_m.sel(time=time_hiscale)
        cldhgh_m_hiscale = cldhgh_m.sel(time=time_hiscale)
        lwnetsfc_m_hiscale = lwnetsfc_m.sel(time=time_hiscale)
        lwnettoa_m_hiscale = lwnettoa_m.sel(time=time_hiscale)
        swnetsfc_m_hiscale = swnetsfc_m.sel(time=time_hiscale)
        swnettoa_m_hiscale = swnettoa_m.sel(time=time_hiscale)
        albedo_m_hiscale = albedo_m.sel(time=time_hiscale)
        
        filename = prep_model_path + 'hiscale_RRM3km_'+IOP+'_sfc.nc'
        modeldata = xr.open_dataset(filename)
        time_m2 = modeldata['time'].load()
        bc_m2 = modeldata['bc'].load()
        dst_m2 = modeldata['dst'].load()
        mom_m2 = modeldata['mom'].load()
        pom_m2 = modeldata['pom'].load()
        ncl_m2 = modeldata['ncl'].load()
        so4_m2 = modeldata['so4'].load()
        soa_m2 = modeldata['soa'].load()
        ccn2_m2 = modeldata['CCN4'].load()
        ncn3_m2 = modeldata['NCN3'].load()
        ncn10_m2 = modeldata['NCN10'].load()
        ncn100_m2 = modeldata['NCN100'].load()
        CNsize_m2 = modeldata['NCNall'].load()
        cod_m2 = modeldata['cod'].load()
        reff_m2 = modeldata['reff'].load()
        lwp_m2 = modeldata['TGCLDLWP'].load()
        nd_m2 = modeldata['Nd_mean'].load()
        precip_m2 = modeldata['PRECT'].load()
        cld_m2 = modeldata['CLDTOT'].load()
        cbh_m2 = modeldata['cbh'].load()
        cth_m2 = modeldata['cth'].load()
        Hcld_m2 = modeldata['clddepth'].load()
        cldlow_m2 = modeldata['CLDLOW'].load()
        cldmid_m2 = modeldata['CLDMED'].load()
        cldhgh_m2 = modeldata['CLDHGH'].load()
        lwdnsfc_m2 = modeldata['FLDS'].load()
        lwnetsfc_m2 = modeldata['FLNS'].load()
        lwnettoa_m2 = modeldata['FLNT'].load()
        lwuptoa_m2 = modeldata['FLUT'].load()
        swdnsfc_m2 = modeldata['FSDS'].load()
        swnetsfc_m2 = modeldata['FSNS'].load()
        swdntoa_m2 = modeldata['SOLIN'].load()
        swnettoa_m2 = modeldata['FSNT'].load()
        swuptoa_m2 = modeldata['FSUTOA'].load()
        modeldata.close()
        lwupsfc_m2 = lwnetsfc_m2 + lwdnsfc_m2
        swupsfc_m2 = swdnsfc_m2 - swnetsfc_m2
        albedo_m2 = swuptoa_m2/swdntoa_m2*100
        org_m2 = pom_m2 + mom_m2 + soa_m2
        bc_m2_hiscale = bc_m2.sel(time=time_hiscale)
        dst_m2_hiscale = dst_m2.sel(time=time_hiscale)
        org_m2_hiscale = org_m2.sel(time=time_hiscale)
        so4_m2_hiscale = so4_m2.sel(time=time_hiscale)
        ncl_m2_hiscale = ncl_m2.sel(time=time_hiscale)
        ccn2_m2_hiscale = ccn2_m2.sel(time=time_hiscale)
        ncn3_m2_hiscale = ncn3_m2.sel(time=time_hiscale)
        ncn10_m2_hiscale = ncn10_m2.sel(time=time_hiscale)
        ncn100_m2_hiscale = ncn100_m2.sel(time=time_hiscale)
        CNsize_m2_hiscale = CNsize_m2.sel(time=time_hiscale)
        cod_m2_hiscale = cod_m2.sel(time=time_hiscale)
        reff_m2_hiscale = reff_m2.sel(time=time_hiscale)
        lwp_m2_hiscale = lwp_m2.sel(time=time_hiscale)
        nd_m2_hiscale = nd_m2.sel(time=time_hiscale)
        precip_m2_hiscale = precip_m2.sel(time=time_hiscale)
        cld_m2_hiscale = cld_m2.sel(time=time_hiscale)
        cbh_m2_hiscale = cbh_m2.sel(time=time_hiscale)
        cth_m2_hiscale = cth_m2.sel(time=time_hiscale)
        Hcld_m2_hiscale = Hcld_m2.sel(time=time_hiscale)
        cldlow_m2_hiscale = cldlow_m2.sel(time=time_hiscale)
        cldmid_m2_hiscale = cldmid_m2.sel(time=time_hiscale)
        cldhgh_m2_hiscale = cldhgh_m2.sel(time=time_hiscale)
        lwnetsfc_m2_hiscale = lwnetsfc_m2.sel(time=time_hiscale)
        lwnettoa_m2_hiscale = lwnettoa_m2.sel(time=time_hiscale)
        swnetsfc_m2_hiscale = swnetsfc_m2.sel(time=time_hiscale)
        swnettoa_m2_hiscale = swnettoa_m2.sel(time=time_hiscale)
        albedo_m2_hiscale = albedo_m2.sel(time=time_hiscale)
        
        filename = prep_model_path + 'hiscale_RRMne30_'+IOP+'_profiles.nc'
        modeldata = xr.open_dataset(filename)
        height_m = modeldata['height'].load()
        cf_e3sm = modeldata['cloud_z'].load()
        modeldata.close()
        cloud_m_hiscale = cf_e3sm.sel(time=time_hiscale)
        
        filename = prep_model_path + 'hiscale_RRM3km_'+IOP+'_profiles.nc'
        modeldata = xr.open_dataset(filename)
        height_m2 = modeldata['height'].load()
        cf_e3sm2 = modeldata['cloud_z'].load()
        modeldata.close()
        cloud_m2_hiscale = cf_e3sm2.sel(time=time_hiscale)
        
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # specific data treatments
        
        # divided by dlogDp in size distribution
        dlogDp_e3sm = np.log10(np.arange(2,3002)/np.arange(1,3001))
        CNsize_m_hiscale = CNsize_m_hiscale.T/dlogDp_e3sm
        CNsize_m2_hiscale = CNsize_m2_hiscale.T/dlogDp_e3sm
        smpsall_hiscale = smpsall_hiscale / dlogDp_smps
        uhsasall_hiscale = uhsasall_hiscale / dlogDp_uhsas
        
        pdf_uhsas_hiscale = np.nanmean(uhsasall_hiscale,axis=0)
        pdf_smps_hiscale = np.nanmean(smpsall_hiscale,axis=0)
        pdf_m_hiscale = np.nanmean(CNsize_m_hiscale,axis=0)
        pdf_m2_hiscale = np.nanmean(CNsize_m2_hiscale,axis=0)
        
        ndrop_hiscale[ndrop_hiscale<10] = np.nan
        nd_sat_hiscale[nd_sat_hiscale<10] = np.nan
        nd_m_hiscale[nd_m_hiscale<10] = np.nan
        nd_m2_hiscale[nd_m2_hiscale<10] = np.nan
        ndrop_hiscale[ndrop_hiscale>500] = np.nan
        nd_sat_hiscale[nd_sat_hiscale>500] = np.nan
        nd_m_hiscale[nd_m_hiscale>500] = np.nan
        nd_m2_hiscale[nd_m2_hiscale>500] = np.nan
        
        lwp_armbe_hiscale[lwp_armbe_hiscale<20] = np.nan
        lwp_mfrsr_hiscale[lwp_mfrsr_hiscale<20] = np.nan
        lwp_sat_hiscale[lwp_sat_hiscale<20] = np.nan
        lwp_m_hiscale[lwp_m_hiscale<20] = np.nan
        lwp_m2_hiscale[lwp_m2_hiscale<20] = np.nan
        
        cod_hiscale[cod_hiscale<2] = np.nan
        cod_sat_hiscale[cod_sat_hiscale<2] = np.nan
        cod_m_hiscale[cod_m_hiscale<2] = np.nan
        cod_m2_hiscale[cod_m2_hiscale<2] = np.nan
        cod_hiscale[cod_hiscale>100] = np.nan
        cod_sat_hiscale[cod_sat_hiscale>100] = np.nan
        cod_m_hiscale[cod_m_hiscale>100] = np.nan
        cod_m2_hiscale[cod_m2_hiscale>100] = np.nan
        
        # unit change:
        precip_m_hiscale = precip_m_hiscale*3600*1000   # m/s to mm/hr
        precip_m2_hiscale = precip_m2_hiscale*3600*1000   # m/s to mm/hr
        cloud_m_hiscale = cloud_m_hiscale*100  # fraction to %
        cloud_m2_hiscale = cloud_m2_hiscale*100  # fraction to %
        height_o = height_o.data*0.001   # m to km
        height_m = height_m.data*0.001   # m to km
        height_m2 = height_m2.data*0.001   # m to km
        
        # set a small threshold of E3SM precipitation
        precip_m_hiscale[precip_m_hiscale<0.02] = 0
        precip_m2_hiscale[precip_m2_hiscale<0.02] = 0
        
        
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        #%% bar plot
        datagroup0 = [org_hiscale,so4_hiscale,nh4_hiscale,no3_hiscale,chl_hiscale, [], []]
        datagroup1 = [org_m_hiscale, so4_m_hiscale, [], [], [], bc_m_hiscale, dst_m_hiscale]
        datagroup2 = [org_m2_hiscale, so4_m2_hiscale, [], [], [], bc_m2_hiscale, dst_m2_hiscale]
        dataall=[datagroup0,datagroup1, datagroup2,]
        labelall = ['Organic', 'SO$_4$', 'NH$_4$', 'NO$_3$', 'Chl', 'BC', 'Dust']
        colorall = ['limegreen', 'red', 'lightblue', 'orange', 'cyan', 'k', 'silver']
        fig,ax = plot.bar(dataall, datalabel=['Obs','ne30','3km_regrid'], xlabel=None, ylabel='unit: $\mu$g/m$^3$', 
                          title='Aerosol Composition  '+site+' '+IOP, varlabel= labelall, colorall=colorall)
        fig.savefig(figpath+'bar_composition_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        #%% mean size distribution
        # fig,ax = plot.mean_size([size_uhsas,size_smps,np.arange(1,3001),np.arange(1,3001)], 
        #             [pdf_uhsas_hiscale, pdf_smps_hiscale, pdf_m_hiscale, pdf_m2_hiscale], 
        #             legend = ['UHSAS','SMPS','ne30','3km_regrid'],color=['k','gray','r','b'], 
        #             marker=['o','+',None,None], linestyles=['none','none','-','-'],
        #             xlimit=(2, 2e3), ylimit=(1e-2,2e4), xlabel='Diameter (nm)', ylabel='dN/dlogDp (cm$^{-3}$)', 
        #             title = 'Mean Aerosol Size Distribution '+site+' '+IOP)
        # fig.savefig(figpath+'mean_aerosol_size_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.mean_size_witherror([size_uhsas,size_smps,np.arange(1,3001),np.arange(1,3001)], 
                                [uhsasall_hiscale,smpsall_hiscale,CNsize_m_hiscale, CNsize_m2_hiscale], 
                                legend = ['UHSAS','SMPS','ne30','3km_regrid'],color=['k','gray','r','b'], 
                                marker=['o','+',None,None], linestyles=['none','none','-','-'],
                          xlimit=(2, 2e3), ylimit=(1e-2,3e4), 
                          xlabel='Diameter (nm)', ylabel='dN/dlogDp (cm$^{-3}$)', 
                            title = 'Mean Aerosol Size Distribution '+site+' '+IOP)
        fig.savefig(figpath+'mean_aerosol_size_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        
        #%% timeseries
        fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [org_hiscale,org_m_hiscale,org_m2_hiscale], 
                                  legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                  title='Total Organic '+site+' '+IOP, xlabel=None, ylabel='${\mu}$g/m$^{3}$')
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'timeseries_org_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [so4_hiscale,so4_m_hiscale,so4_m2_hiscale], 
                                  legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'],  
                                  title='Sulfate '+site+' '+IOP, xlabel=None, ylabel='${\mu}$g/m$^{3}$')
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'timeseries_so4_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [ccn2_hiscale,ccn2_m_hiscale,ccn2_m2_hiscale], 
                                  legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                title='0.2%CCN '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'timeseries_CCN2_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [cpc3_hiscale,ncn3_m_hiscale,ncn3_m2_hiscale], 
                                  legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                  title='CN(>3nm) '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'timeseries_CPC3_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [cpc10_hiscale,ncn10_m_hiscale,ncn10_m2_hiscale], 
                                  legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                  title='CN(>10nm) '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'timeseries_CPC10_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale,time_hiscale], [smps100_hiscale,uhsas100_hiscale,ncn100_m_hiscale,ncn100_m2_hiscale], 
                                legend = ['SMPS','UHSAS','ne30','3km_regrid'], color=['k','gray','r','b'],
                                title='CN(>100nm) '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'timeseries_CN100_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale,time_hiscale], [cod_hiscale,cod_sat_hiscale,cod_m_hiscale,cod_m2_hiscale], 
                                  legend = ['MFRSR','Satellite','ne30','3km_regrid'], color=['k','gray','r','b'], #marker='.',
                                title='cloud optical depth '+site+' '+IOP, xlabel=None, ylabel=None)
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'timeseries_cod_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale,time_hiscale], [lwp_armbe_hiscale, lwp_sat_hiscale, lwp_m_hiscale, lwp_m2_hiscale], 
                                legend = ['ARMBE','Satellite','ne30','3km_regrid'], color=['k','gray','r','b'],
                                title='LWP '+site+' '+IOP,xlabel=None, ylabel="g/m$^2$")
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'timeseries_LWP_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale,time_hiscale], [ndrop_hiscale,nd_sat_hiscale, nd_m_hiscale, nd_m2_hiscale], 
                                  legend = ['Ndrop','Satellite','ne30','3km_regrid'], color=['k','gray','r','b'], #marker='.',
                                  title='Nd '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'timeseries_Nd_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale,time_hiscale], [reff_hiscale,reff_sat_hiscale,reff_m_hiscale,reff_m2_hiscale],  
                                legend = ['MFRSR','Satellite','ne30','3km_regrid'], color=['k','gray','r','b'],marker='.',
                                title='Reff '+site+' '+IOP,xlabel=None, ylabel='$\mu$m')
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'timeseries_reff_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [precip_hiscale,precip_m_hiscale,precip_m2_hiscale],  
                                  legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                title='Precip '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='mm/hr')
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'timeseries_precip_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [lwnetsfc_hiscale,lwnetsfc_m_hiscale,lwnetsfc_m2_hiscale], 
                                  legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                title='Sfc. net LW Flux '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='W/m$^2$')
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'timeseries_LWsfc_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [swnetsfc_hiscale,swnetsfc_m_hiscale,swnetsfc_m2_hiscale], 
                                  legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                title='Sfc. net SW Flux '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='W/m$^2$')
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'timeseries_SWsfc_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [lwnettoa_hiscale,lwnettoa_m_hiscale,lwnettoa_m2_hiscale], 
                                  legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                title='TOA. net LW Flux '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='W/m$^2$')
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'timeseries_LWtoa_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [swnettoa_hiscale,swnettoa_m_hiscale,swnettoa_m2_hiscale], 
                                  legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                title='TOA. net SW Flux '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='W/m$^2$')
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'timeseries_SWtoa_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale,time_hiscale], [cld_arscl_hiscale,cld_visst_hiscale,cld_m_hiscale,cld_m2_hiscale], 
                                legend = ['ARSCL','Satellite','ne30','3km_regrid'], color=['k','gray','r','b'],
                                title='Cloud fraction '+site+' '+IOP,xlabel=None, ylabel="%")
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'timeseries_totcld_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        #%%
        fig,ax = plot.timeseries_size([time_hiscale,time_hiscale,time_hiscale,time_hiscale], 
                                      [size_uhsas,size_smps, np.arange(1,3001), np.arange(1,3001)], 
                                      [uhsasall_hiscale.T.data, smpsall_hiscale.T.data, CNsize_m_hiscale.T.data, CNsize_m2_hiscale.T.data], 
                                      legend = ['UHSAS','SMPS','ne30','3km_regrid'],
                                  ylabel='Diameter (nm)', ylimit=(3,1000),
                                  title = 'Aerosol Size Distribution (dN/dlogDp, cm$^{-3}$)')
        # for ax_i in ax:
        #     ax_i.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        #     ax_i.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'aerosol_size_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.timeseries_2d([time_hiscale,time_hiscale,time_hiscale], 
                                    [height_o, height_m, height_m2], 
                                    [cloud_2d_hiscale.T.data, cloud_m_hiscale.T.data, cloud_m2_hiscale.T.data], 
                                      yticks=[0,3,6,9,12], ylimit=(0,12), ylabel='Height (m)',cmap='jet', #ylimit=(3,1000),
                                      legend = ['Obs','ne30','3km_regrid'], title = 'Cloud Fraction (%)')
        # for ax_i in ax:
        #     ax_i.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        #     ax_i.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        fig.savefig(figpath+'cloud_2d_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        
        #%% diurnal cycle
        fig,ax = plot.diurnalcycle([org_hiscale,org_m_hiscale,org_m2_hiscale], legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                  title='Organic '+site+' '+IOP, xlabel='Time (UTC)', ylabel='${\mu}$g/m$^{3}$')
        fig.savefig(figpath+'diurnalcycle_org_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        fig,ax = plot.diurnalcycle([so4_hiscale,so4_m_hiscale,so4_m2_hiscale], legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                  title='Sulfate '+site+' '+IOP, xlabel='Time (UTC)', ylabel='${\mu}$g/m$^{3}$')
        fig.savefig(figpath+'diurnalcycle_so4_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.diurnalcycle([ccn2_hiscale,ccn2_m_hiscale,ccn2_m2_hiscale], legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                title='0.2%CCN '+site+' '+IOP, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
        fig.savefig(figpath+'diurnalcycle_CCN2_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.diurnalcycle([cpc3_hiscale,ncn3_m_hiscale,ncn3_m2_hiscale], legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                title='CN(>3nm) '+site+' '+IOP, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
        fig.savefig(figpath+'diurnalcycle_CPC3_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        fig,ax = plot.diurnalcycle([cpc10_hiscale,ncn10_m_hiscale,ncn10_m2_hiscale], legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                title='CN(>10nm) '+site+' '+IOP, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
        fig.savefig(figpath+'diurnalcycle_CPC10_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        fig,ax = plot.diurnalcycle([smps100_hiscale,uhsas100_hiscale,ncn100_m_hiscale,ncn100_m2_hiscale], legend = ['SMPS100','UHSAS100','ne30','3km_regrid'], 
                                title='CN(>100nm) '+site+' '+IOP, color=['k','gray','r','b'], xlabel='Time (UTC)',ylabel='cm$^{-3}$')
        fig.savefig(figpath+'diurnalcycle_CN100_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.diurnalcycle( [cod_hiscale, cod_sat_hiscale, cod_m_hiscale, cod_m2_hiscale], 
                                    legend = ['MFRSR','Satellite','ne30','3km_regrid'], color=['k','gray','r','b'], 
                                title='Cloud optical depth '+site+' '+IOP, xlabel='Time (UTC)', ylabel=None)
        fig.savefig(figpath+'diurnalcycle_cod_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.diurnalcycle([lwp_armbe_hiscale,lwp_sat_hiscale, lwp_m_hiscale, lwp_m2_hiscale], 
                                    legend = ['ARMBE','Satellite','ne30','3km_regrid'], color=['k','gray','r','b'],
                                title='LWP '+site+' '+IOP,  xlabel='Time (UTC)',ylabel="g/m$^2$")
        fig.savefig(figpath+'diurnalcycle_LWP_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.diurnalcycle([ndrop_hiscale, nd_sat_hiscale, nd_m_hiscale,nd_m2_hiscale], 
                                    legend = ['Ndrop', 'Satellite','ne30','3km_regrid'], color=['k','gray','r','b'], 
                                  title='Nd '+site+' '+IOP, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
        fig.savefig(figpath+'diurnalcycle_Nd_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.diurnalcycle([reff_hiscale, reff_sat_hiscale, reff_m_hiscale, reff_m2_hiscale], 
                                    legend = ['MFRSR','Satellite','ne30','3km_regrid'], color=['k','gray','r','b'],
                                title='droplet effective radius '+site+' '+IOP, xlabel='Time (UTC)', ylabel='$\mu$m')
        fig.savefig(figpath+'diurnalcycle_reff_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.diurnalcycle( [precip_hiscale,precip_m_hiscale,precip_m2_hiscale], 
                                    legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                nozero_percentile=True, title='Precipitation '+site+' '+IOP, xlabel='Time (UTC)',ylabel='mm/hr')
        fig.savefig(figpath+'diurnalcycle_precip_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.diurnalcycle([lwnetsfc_hiscale,lwnetsfc_m_hiscale,lwnetsfc_m2_hiscale], legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                title='Sfc. net LW Flux '+site+' '+IOP, xlabel='Time (UTC)',ylabel='W/m$^2$')
        fig.savefig(figpath+'diurnalcycle_LWsfc_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        fig,ax = plot.diurnalcycle([swnetsfc_hiscale,swnetsfc_m_hiscale,swnetsfc_m2_hiscale], legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                title='Sfc. net SW Flux '+site+' '+IOP, xlabel='Time (UTC)', ylabel='W/m$^2$')
        fig.savefig(figpath+'diurnalcycle_SWsfc_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        fig,ax = plot.diurnalcycle([lwnettoa_hiscale,lwnettoa_m_hiscale,lwnettoa_m2_hiscale], legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                title='TOA. net LW Flux '+site+' '+IOP, xlabel='Time (UTC)', ylabel='W/m$^2$')
        fig.savefig(figpath+'diurnalcycle_LWtoa_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        fig,ax = plot.diurnalcycle([swnettoa_hiscale,swnettoa_m_hiscale,swnettoa_m2_hiscale], legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'], 
                                title='TOA. net SW Flux '+site+' '+IOP, xlabel='Time (UTC)', ylabel='W/m$^2$')
        fig.savefig(figpath+'diurnalcycle_SWtoa_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.diurnalcycle([cld_arscl_hiscale,cld_visst_hiscale,cld_m_hiscale,cld_m2_hiscale],
                                    legend = ['ARSCL','Satellite','ne30','3km_regrid'], color=['k','gray','r','b'],
                                    title='Total cloud fraction '+site+' '+IOP, xlabel='Time (UTC)', ylabel="%")
        fig.savefig(figpath+'diurnalcycle_totcld_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.diurnalcycle_2d([uhsasall_hiscale.T, smpsall_hiscale.T, CNsize_m_hiscale.T, CNsize_m2_hiscale.T], 
                                      y=[size_uhsas,size_smps, np.arange(1,3001), np.arange(1,3001)], 
                                      title= ['UHSAS','SMPS','ne30','3km_regrid'],
                                      levellist=np.arange(0,11500,200), xlabel='Time (UTC)', ylabel='Diameter (nm)', 
                                      ylimit=(3,1000),cmap='jet')
        for ax_i in ax:
            ax_i.set_yscale('log')
        fig.savefig(figpath+'diurnalcycle_aerosol_size_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.diurnalcycle_2d([cloud_2d_hiscale, cloud_m_hiscale, cloud_m2_hiscale], 
                                      y = [height_o, height_m, height_m2],
                                yticks=[0,3,6,9,12], ylimit=(0,12), ylabel='Height (km)',  cmap='jet',
                                levellist=np.arange(0,45,1),
                                  title= ['Obs','ne30','3km_regrid'])
        fig.savefig(figpath+'diurnalcycle_cloud2d_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        #%% 1d histogram
        
        w1 = np.ones_like(org_hiscale)/sum(~np.isnan(org_hiscale.data))
        w2 = np.ones_like(org_m_hiscale)/sum(~np.isnan(org_m_hiscale.data))
        w3 = np.ones_like(org_m2_hiscale)/sum(~np.isnan(org_m2_hiscale.data))
        fig,ax = plot.hist([org_hiscale,org_m_hiscale,org_m2_hiscale], weights=[w1,w2,w3], bins=np.arange(0,10,0.3),
                            legend =['Obs','ne30','3km_regrid'], color=['k','r','b'],
                            title = 'Total Organic '+site+' '+IOP, ylabel='Fraction', xlabel='${\mu}$g/m$^{3}$')
        fig.savefig(figpath+'hist_org_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        w1 = np.ones_like(so4_hiscale)/sum(~np.isnan(so4_hiscale.data))
        w2 = np.ones_like(so4_m_hiscale)/sum(~np.isnan(so4_m_hiscale.data))
        w3 = np.ones_like(so4_m2_hiscale)/sum(~np.isnan(so4_m2_hiscale.data))
        fig,ax = plot.hist([so4_hiscale,so4_m_hiscale,so4_m2_hiscale], weights=[w1,w2,w3], bins=np.arange(0,6,0.2),
                            legend =['Obs','ne30','3km_regrid'], color=['k','r','b'],
                            title = 'Sulfate '+site+' '+IOP, ylabel='Fraction', xlabel='${\mu}$g/m$^{3}$')
        fig.savefig(figpath+'hist_SO4_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        w1 = np.ones_like(ccn2_hiscale)/sum(~np.isnan(ccn2_hiscale.data))
        w2 = np.ones_like(ccn2_m_hiscale)/sum(~np.isnan(ccn2_m_hiscale.data))
        w3 = np.ones_like(ccn2_m2_hiscale)/sum(~np.isnan(ccn2_m2_hiscale.data))
        fig,ax = plot.hist([ccn2_hiscale,ccn2_m_hiscale,ccn2_m2_hiscale], weights=[w1,w2,w3], bins=np.arange(0,1800,50),
                            legend =['Obs','ne30','3km_regrid'], color=['k','r','b'],
                            title = 'CCN (SS=0.2%) '+site+' '+IOP, ylabel='Fraction', xlabel='cm$^{-3}$')
        fig.savefig(figpath+'hist_CCN2_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        
        w0 = np.ones_like(cpc10_hiscale)/sum(~np.isnan(cpc10_hiscale.data))
        w1 = np.ones_like(ncn10_m_hiscale)/sum(~np.isnan(ncn10_m_hiscale.data))
        w2 = np.ones_like(ncn10_m2_hiscale)/sum(~np.isnan(ncn10_m2_hiscale.data))
        fig,ax = plot.hist([cpc10_hiscale,ncn10_m_hiscale,ncn10_m2_hiscale], weights=[w0,w1,w2], bins=np.arange(0,22000,1000),
                            legend = ['Obs','ne30','3km_regrid'], color=['k','r','b'],
                            title='Aerosol number (>10nm) '+site+' '+IOP,ylabel='Fraction', xlabel='cm$^{-3}$')
        fig.savefig(figpath+'hist_CPC10_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        w0 = np.ones_like(smps100_hiscale)/sum(~np.isnan(smps100_hiscale.data))
        w1 = np.ones_like(ncn100_m_hiscale)/sum(~np.isnan(ncn100_m_hiscale.data))
        w2 = np.ones_like(ncn100_m2_hiscale)/sum(~np.isnan(ncn100_m2_hiscale.data))
        fig,ax = plot.hist([smps100_hiscale,ncn100_m_hiscale,ncn100_m2_hiscale], weights=[w0,w1,w2], bins=np.arange(0,2100,100),
                            legend = ['SMPS100','ne30','3km_regrid'], color=['k','r','b'],
                            title='Aerosol number (>100nm) '+site+' '+IOP,  ylabel='Fraction', xlabel='cm$^{-3}$')
        fig.savefig(figpath+'hist_CN100_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        w0 = np.ones_like(cod_hiscale)/sum(~np.isnan(cod_hiscale.data))
        w00 = np.ones_like(cod_sat_hiscale)/sum(~np.isnan(cod_sat_hiscale.data))
        w1 = np.ones_like(cod_m_hiscale)/sum(~np.isnan(cod_m_hiscale.data))
        w2 = np.ones_like(cod_m2_hiscale)/sum(~np.isnan(cod_m2_hiscale.data))
        fig,ax = plot.hist( [cod_hiscale, cod_sat_hiscale, cod_m_hiscale, cod_m2_hiscale], weights=[w0,w00,w1,w2], 
                            legend = ['MFRSR','Satellite','ne30','3km_regrid'], color=['k','gray','r','b'],
                            title='Cloud Optical Depth '+site+' '+IOP, bins=np.arange(0,61,3), ylabel='Fraction', xlabel='N/A')
        fig.savefig(figpath+'hist_cod_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        w0 = np.ones_like(lwp_armbe_hiscale)/sum(~np.isnan(lwp_armbe_hiscale.data))
        # w0 = np.ones_like(lwp_mfrsr)/sum(~np.isnan(lwp_mfrsr.data))
        w00 = np.ones_like(lwp_sat_hiscale)/sum(~np.isnan(lwp_sat_hiscale.data))
        w1 = np.ones_like(lwp_m_hiscale)/sum(~np.isnan(lwp_m_hiscale.data))
        w2 = np.ones_like(lwp_m2_hiscale)/sum(~np.isnan(lwp_m2_hiscale.data))
        fig,ax = plot.hist([lwp_mfrsr_hiscale, lwp_sat_hiscale, lwp_m_hiscale, lwp_m2_hiscale], weights=[w0,w00,w1,w2], 
                            legend = ['ARMBE','Satellite','ne30','3km_regrid'], color=['k','gray','r','b'],
                            title='LWP '+site+' '+IOP, bins=np.arange(10,410,20), ylabel='Fraction', xlabel="g/m$^2$")
        fig.savefig(figpath+'hist_LWP_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        w0 = np.ones_like(ndrop_hiscale)/sum(~np.isnan(ndrop_hiscale.data))
        w00 = np.ones_like(nd_sat_hiscale)/sum(~np.isnan(nd_sat_hiscale.data))
        w1 = np.ones_like(nd_m_hiscale)/sum(~np.isnan(nd_m_hiscale.data))
        w2 = np.ones_like(nd_m2_hiscale)/sum(~np.isnan(nd_m2_hiscale.data))
        fig,ax = plot.hist([ndrop_hiscale,nd_sat_hiscale,nd_m_hiscale,nd_m2_hiscale],  weights=[w0,w00,w1,w2], 
                            legend = ['Ndrop','Satellite','ne30','3km_regrid'], color=['k','gray','r','b'],
                            title = 'Nd '+site+' '+IOP, bins=np.arange(0,410,20), ylabel='Fraction', xlabel='cm$^{-3}$')
        fig.savefig(figpath+'hist_Nd_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        w0 = np.ones_like(reff_hiscale)/sum(~np.isnan(reff_hiscale.data))
        w00 = np.ones_like(reff_sat_hiscale)/sum(~np.isnan(reff_sat_hiscale.data))
        w1 = np.ones_like(reff_m_hiscale)/sum(~np.isnan(reff_m_hiscale.data))
        w2 = np.ones_like(reff_m2_hiscale)/sum(~np.isnan(reff_m2_hiscale.data))
        fig,ax = plot.hist([reff_hiscale,reff_sat_hiscale,reff_m_hiscale,reff_m2_hiscale], weights=[w0,w00,w1,w2], 
                            legend = ['MFRSR','Satellite','ne30','3km_regrid'], color=['k','gray','r','b'],
                            title = 'Cloud Effective Radius '+site+' '+IOP, bins=np.arange(4,28,1), ylabel='Fraction', xlabel='$\mu$m')
        fig.savefig(figpath+'hist_reff_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        pr0 = precip_hiscale[precip_hiscale!=0]
        prm = precip_m_hiscale[precip_m_hiscale!=0]
        prm2 = precip_m_hiscale[precip_m2_hiscale!=0]
        w0 = np.ones_like(pr0)/sum(~np.isnan(pr0.data))
        w1 = np.ones_like(prm)/sum(~np.isnan(prm.data))
        w2 = np.ones_like(prm2)/sum(~np.isnan(prm2.data))
        fig,ax = plot.hist( [pr0,prm,prm2], weights=[w0,w1,w2], legend = ['Obs','ne30','3km_regrid'], 
                            color=['k','r','b'],  bins=np.arange(0,5,.1), 
                            title = 'Precipitation '+site+' '+IOP, ylabel='Fraction', xlabel='mm/hr')
        fig.savefig(figpath+'hist_precip_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        w0 = np.ones_like(cld_arscl_hiscale)/sum(~np.isnan(cld_arscl_hiscale.data))
        w00 = np.ones_like(cld_visst_hiscale)/sum(~np.isnan(cld_visst_hiscale.data))
        w1 = np.ones_like(cld_m_hiscale)/sum(~np.isnan(cld_m_hiscale.data))
        w2 = np.ones_like(cld_m2_hiscale)/sum(~np.isnan(cld_m2_hiscale.data))
        fig,ax = plot.hist([cld_arscl_hiscale,cld_visst_hiscale,cld_m_hiscale,cld_m2_hiscale], 
                            weights=[w0,w00,w1,w2],  bins=np.arange(0,101,5), 
                            legend = ['ARMBE','Satellite','ne30','3km_regrid'], color=['k','gray','r','b'],
                              title = 'Cloud Fraction '+site+' '+IOP, ylabel='Fraction', xlabel="%")
        fig.savefig(figpath+'hist_totcld_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        
        #%% calculate statistics
        calc.mean_std_percentiles([org_hiscale,org_m_hiscale,org_m2_hiscale],legend=['Obs','ne30','3km_regrid'], 
                                  outfile=figpath+'statistics_1var_ORG_'+site+'_'+IOP+'.txt')
        calc.mean_std_percentiles([so4_hiscale, so4_m_hiscale, so4_m2_hiscale],legend=['Obs','ne30','3km_regrid'], 
                                  outfile=figpath+'statistics_1var_SO4_'+site+'_'+IOP+'.txt')
        calc.mean_std_percentiles([ccn2_hiscale,ccn2_m_hiscale,ccn2_m2_hiscale],legend=['Obs','ne30','3km_regrid'], 
                                  outfile=figpath+'statistics_1var_CCN2_'+site+'_'+IOP+'.txt')
        calc.mean_std_percentiles([cpc3_hiscale,ncn3_m_hiscale,ncn3_m2_hiscale],legend=['Obs','ne30','3km_regrid'], 
                                  outfile=figpath+'statistics_1var_CPC3_'+site+'_'+IOP+'.txt')
        calc.mean_std_percentiles([cpc10_hiscale,ncn10_m_hiscale,ncn10_m2_hiscale],legend=['Obs','ne30','3km_regrid'], 
                                  outfile=figpath+'statistics_1var_CPC10_'+site+'_'+IOP+'.txt')
        calc.mean_std_percentiles([uhsas100_hiscale, smps100_hiscale, ncn100_m_hiscale, ncn100_m2_hiscale],legend=['UHSAS','SMPS','ne30','3km_regrid'], 
                                  outfile=figpath+'statistics_1var_CN100_'+site+'_'+IOP+'.txt')
        calc.mean_std_percentiles([cod_hiscale,cod_sat_hiscale, cod_m_hiscale, cod_m2_hiscale],legend=['MFRSR','Satellite','ne30','3km_regrid'],
                                  outfile=figpath+'statistics_1var_COD_'+site+'_'+IOP+'.txt')
        calc.mean_std_percentiles([reff_hiscale,reff_sat_hiscale,reff_m_hiscale,reff_m2_hiscale],legend=['MFRSR','Satellite','ne30','3km_regrid'],
                                  outfile=figpath+'statistics_1var_Reff_'+site+'_'+IOP+'.txt')
        calc.mean_std_percentiles([lwp_mfrsr_hiscale,lwp_armbe_hiscale,lwp_sat_hiscale,lwp_m_hiscale,lwp_m2_hiscale],
                                  legend=['MFRSR','ARMBE','Satellite','ne30','3km_regrid'],
                                  outfile=figpath+'statistics_1var_LWP_'+site+'_'+IOP+'.txt')
        calc.mean_std_percentiles([ndrop_hiscale,nd_sat_hiscale,nd_m_hiscale,nd_m2_hiscale],legend=['Ndrop','Nd_satellite','ne30','3km_regrid'],
                                  outfile=figpath+'statistics_1var_Nd_'+site+'_'+IOP+'.txt')
        calc.mean_std_percentiles([precip_hiscale,precip_m_hiscale,precip_m2_hiscale],legend=['Obs','ne30','3km_regrid'],
                                  outfile=figpath+'statistics_1var_Precip_'+site+'_'+IOP+'.txt')
        calc.mean_std_percentiles([cld_arscl_hiscale,cld_visst_hiscale,cld_tsi_hiscale,cld_m_hiscale,cld_m2_hiscale],
                                  legend=['ARSCL','Satellite','TSI','ne30','3km_regrid'],
                                  outfile=figpath+'statistics_1var_totcld_'+site+'_'+IOP+'.txt')
        
        
        calc.bias_corrcoef_RMSE(org_hiscale,org_m_hiscale,label1='Obs',label2='RRMne30', 
                                outfile=figpath+'statistics_ORG_RRMne30vsOBS_'+site+'_'+IOP+'.txt')
        calc.bias_corrcoef_RMSE(org_hiscale,org_m2_hiscale,label1='Obs',label2='RRM3km', 
                                outfile=figpath+'statistics_ORG_RRM3kmvsOBS_'+site+'_'+IOP+'.txt')
        
        calc.bias_corrcoef_RMSE(so4_hiscale, so4_m_hiscale,label1='Obs',label2='RRMne30', 
                                outfile=figpath+'statistics_SO4_RRMne30vsOBS_'+site+'_'+IOP+'.txt')
        calc.bias_corrcoef_RMSE(so4_hiscale, so4_m2_hiscale,label1='Obs',label2='RRM3km', 
                                outfile=figpath+'statistics_SO4_RRM3kmvsOBS_'+site+'_'+IOP+'.txt')
        
        calc.bias_corrcoef_RMSE(ccn2_hiscale,ccn2_m_hiscale,label1='Obs',label2='RRMne30', 
                                outfile=figpath+'statistics_CCN2_RRMne30vsOBS_'+site+'_'+IOP+'.txt')
        calc.bias_corrcoef_RMSE(ccn2_hiscale,ccn2_m2_hiscale,label1='Obs',label2='RRM3km', 
                                outfile=figpath+'statistics_CCN2_RRM3kmvsOBS_'+site+'_'+IOP+'.txt')
        
        if IOP=='IOP1':
            calc.bias_corrcoef_RMSE(cpc3_hiscale,ncn3_m_hiscale,label1='Obs',label2='RRMne30', 
                                outfile=figpath+'statistics_CN3nm_RRMne30vsOBS_'+site+'_'+IOP+'.txt')
            calc.bias_corrcoef_RMSE(cpc3_hiscale,ncn3_m2_hiscale,label1='Obs',label2='RRM3km', 
                                outfile=figpath+'statistics_CN3nm_RRM3kmvsOBS_'+site+'_'+IOP+'.txt')
            calc.bias_corrcoef_RMSE(cpc10_hiscale,ncn10_m_hiscale,label1='Obs',label2='RRMne30', 
                                outfile=figpath+'statistics_CN10nm_RRMne30vsOBS_'+site+'_'+IOP+'.txt')
            calc.bias_corrcoef_RMSE(cpc10_hiscale,ncn10_m2_hiscale,label1='Obs',label2='RRM3km', 
                                outfile=figpath+'statistics_CN10nm_RRM3kmvsOBS_'+site+'_'+IOP+'.txt')
        
        calc.bias_corrcoef_RMSE(smps100_hiscale, ncn100_m_hiscale,label1='Obs',label2='RRMne30', 
                                outfile=figpath+'statistics_CN100_RRMne30vsOBS_'+site+'_'+IOP+'.txt')
        calc.bias_corrcoef_RMSE(smps100_hiscale, ncn100_m2_hiscale,label1='Obs',label2='RRM3km', 
                                outfile=figpath+'statistics_CN100_RRM3kmvsOBS_'+site+'_'+IOP+'.txt')
        
        calc.bias_corrcoef_RMSE(lwp_armbe_hiscale, lwp_m_hiscale,label1='ARMBE',label2='RRMne30', 
                                outfile=figpath+'statistics_lwp_RRMne30vsOBS_'+site+'_'+IOP+'.txt')
        calc.bias_corrcoef_RMSE(lwp_armbe_hiscale, lwp_m2_hiscale,label1='ARMBE',label2='RRM3km', 
                                outfile=figpath+'statistics_lwp_RRM3kmvsOBS_'+site+'_'+IOP+'.txt')
        
        calc.bias_corrcoef_RMSE(ndrop_hiscale, nd_m_hiscale,label1='Ndrop',label2='RRMne30', 
                                outfile=figpath+'statistics_Nd_RRMne30vsOBS_'+site+'_'+IOP+'.txt')
        calc.bias_corrcoef_RMSE(ndrop_hiscale, nd_m2_hiscale,label1='Ndrop',label2='RRM3km', 
                                outfile=figpath+'statistics_Nd_RRM3kmvsOBS_'+site+'_'+IOP+'.txt')
        
        calc.bias_corrcoef_RMSE(nd_sat_hiscale, nd_m_hiscale,label1='Satellite',label2='RRMne30', 
                                outfile=figpath+'statistics_Nd_E3SMv2vsSat_'+site+'_'+IOP+'.txt')
        calc.bias_corrcoef_RMSE(nd_sat_hiscale, nd_m2_hiscale,label1='Satellite',label2='RRM3km', 
                                outfile=figpath+'statistics_Nd_RRM3kmvsSat_'+site+'_'+IOP+'.txt')
        
        # calc.bias_corrcoef_RMSE(reff_hiscale, reff_m_hiscale,label1='Obs',label2='RRMne30', 
        #                         outfile=figpath+'statistics_Reff_RRMne30vsOBS_'+site+'_'+IOP+'.txt')
        # calc.bias_corrcoef_RMSE(reff_hiscale, reff_m2_hiscale,label1='Obs',label2='RRM3km', 
        #                         outfile=figpath+'statistics_Reff_RRM3kmvsOBS_'+site+'_'+IOP+'.txt')
        
        # calc.bias_corrcoef_RMSE(reff_sat_hiscale, reff_m_hiscale,label1='Satellite',label2='RRMne30', 
        #                         outfile=figpath+'statistics_Reff_E3SMv2vsSat_'+site+'_'+IOP+'.txt')
        # calc.bias_corrcoef_RMSE(reff_sat_hiscale, reff_m2_hiscale,label1='Satellite',label2='RRM3km', 
        #                         outfile=figpath+'statistics_Reff_E3SMv2'+modeltype+'vsSat_'+site+'_'+IOP+'.txt')
        
        #%% joint histogram
        fig,ax = plot.jointhist([uhsas100_hiscale,ncn100_m_hiscale,ncn100_m2_hiscale], [ccn2_hiscale,ccn2_m_hiscale,ccn2_m2_hiscale], 
                            xedges=np.arange(0,800,40),yedges=np.arange(0,800,40), normalize_x=True,
                            xlabel='CN (>100nm) (cm$^{-3}$)', ylabel='CCN (SS=0.2%) (cm$^{-3}$)', vmax=0.5,
                            title=['Ground','ne30','3km_regrid'])
        fig.savefig(figpath+'jointhist_CN100_CCN2_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.jointhist([ccn2_hiscale,ccn2_hiscale,ccn2_m_hiscale,ccn2_m2_hiscale], 
                                [ndrop_hiscale,nd_sat_hiscale,nd_m_hiscale,nd_m2_hiscale],
                            xedges=np.arange(0,500,30),yedges=np.arange(0,300,20), normalize_x=True,
                            xlabel='CCN (SS=0.2%) (cm$^{-3}$)', ylabel='Nd (cm$^{-3}$)', vmax=0.4,
                            title=['Ground','Satellite','ne30','3km_regrid'])
        fig.savefig(figpath+'jointhist_CCN2_Nd_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.jointhist([ndrop_hiscale,nd_sat_hiscale,nd_m_hiscale,nd_m2_hiscale],
                                [lwp_sat_hiscale,lwp_sat_hiscale,lwp_m_hiscale,lwp_m2_hiscale], 
                            xedges=np.arange(0,300,20),yedges=np.arange(0,300,20), normalize_x=True,
                            xlabel='Nd (cm$^{-3}$)', ylabel='LWP (g/m$^2$)', vmax=0.4,
                            title=['Ground','Satellite','ne30','3km_regrid'])
        fig.savefig(figpath+'jointhist_LWP_Nd_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.jointhist([ndrop_hiscale,nd_sat_hiscale,nd_m_hiscale,nd_m2_hiscale],
                                [reff_hiscale,reff_sat_hiscale,reff_m_hiscale,reff_m2_hiscale],
                            xedges=np.arange(0,300,20),yedges=np.arange(4,25,1), normalize_x=True,
                            xlabel='Nd (cm$^{-3}$)', ylabel='Reff ($\mu$m)', vmax=0.25,
                            title=['Ground','Satellite','ne30','3km_regrid'])
        fig.savefig(figpath+'jointhist_Reff_Nd_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        fig,ax = plot.jointhist([cod_sat_hiscale,cod_sat_hiscale,cod_m_hiscale,cod_m2_hiscale],[lwp_armbe_hiscale,lwp_sat_hiscale,lwp_m_hiscale,lwp_m2_hiscale], 
                            xedges=np.arange(0,40,3),yedges=np.arange(0,300,20), normalize_x=True,
                            xlabel='Cloud Optical Depth (N/A)', ylabel='LWP (g/m$^2$)', vmax=0.25,
                            title=['Ground','Satellite','ne30','3km_regrid'])
        fig.savefig(figpath+'jointhist_COD_Nd_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        #%% scatter plot
        
        fig,ax = plot.scatter([ndrop_hiscale.data, nd_sat_hiscale.data,nd_m_hiscale.data,nd_m2_hiscale.data], 
                              [ccn2_hiscale.data,ccn2_hiscale.data,ccn2_m_hiscale.data,ccn2_m2_hiscale.data],
                              xlimit=(0,300), ylimit=(0,600),
                            xlabel='Nd (cm$^{-3}$)', ylabel='Surface CCN (SS=0.2%) (cm$^{-3}$)', 
                            title=['Ground','Satellite','ne30','3km_regrid'],
                        linear_fit=True, intercept=False)
        fig.savefig(figpath+'scatter_Nd_CCN2_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        fig,ax = plot.scatter([smps100_hiscale.data,ncn100_m_hiscale.data,ncn100_m2_hiscale.data], 
                              [ccn2_hiscale.data,ccn2_m_hiscale.data,ccn2_m2_hiscale.data],
                              xlimit=(0,800), ylimit=(0,800),
                            xlabel='Surface CN (>100nm) (cm$^{-3}$)', ylabel='Surface CCN (SS=0.2%) (cm$^{-3}$)', 
                            title=['Ground','ne30','3km_regrid'],
                        linear_fit=True, intercept=True)
        fig.savefig(figpath+'scatter_CN100_CCN2_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        #%% heatmaps
        
        # xedges=np.exp(np.arange(np.log(10),6.5,0.5))
        # yedges=np.exp(np.arange(np.log(10),6.5,0.5))
        fig,ax = plot.heatmap([ndrop_hiscale.data, nd_sat_hiscale.data,nd_m_hiscale.data,nd_m2_hiscale.data],
                              [lwp_sat_hiscale,lwp_sat_hiscale,lwp_m_hiscale,lwp_m2_hiscale],
                              [albedo_hiscale,albedo_hiscale,albedo_m_hiscale,albedo_m2_hiscale],vmax=60,
                            xedges=np.arange(0,300,20), yedges=np.arange(10,300,20),
                            # xedges=xedges, yedges=yedges, 
                            xlabel='Nd (cm$^{-3}$)', ylabel='LWP (g/m$^2$)', zlabel='TOA Albedo (%)',
                            title=['Ground','Satellite','ne30','3km_regrid'])
        fig.savefig(figpath+'heatmap_Albedo_vs_Nd_LWP_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
