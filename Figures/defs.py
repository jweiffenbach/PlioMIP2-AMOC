import numpy as np
import xarray as xr
import xskillscore

def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    #Code from Oldeman (2021)
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]

def regularAtlmask(ds):
    
    #Create Atlantic mask on regular 1x1 grid and exclude Mediterranean
    import regionmask

    Atlantic = np.array([[-110,30],[-75,50],[-75, 65],[-90,75],[-80,80],[-40,80], [25,80], [25,70],[10,60],[-10,35], [-10,10], [20,10], [20,-34], [-65,-34], [-40, -10], [-70,5], [-100,20]])
    region = regionmask.Regions([Atlantic])

    # define lat/ lon grid
    dsgrid = ds.copy(deep=True)
    mask = region.mask(ds.lon.values, ds.lat.values)
    return(mask)

def lon180(ds):
    "Rewrite 0-360 longitude to -180-180"
    ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
    ds = ds.sortby(ds.lon)
    return ds

def amocstrength(model, experiment):
    "Computes AMOC strength by taking the maximum streamfunction below 500 m depth and north of 0N"
    folder = '/Volumes/External/DataPlioMIP2/Data/Processed/'
    dsmoc = xr.open_dataset(folder+model+'/'+experiment+'/AMOC_annual_100yr.nc')
    amoc = dsmoc.AMOC.where(dsmoc.AMOC<1e10).where(dsmoc.AMOC>-1e10).where(dsmoc.z>500).where(dsmoc.lat>0).max(dim=['z','lat'])
    return amoc 

def amocstrength_100yr(model, experiment):
    folder = '/Volumes/External/DataPlioMIP2/Data/Processed/'
    dsmoc = xr.open_dataset(folder+model+'/'+experiment+'/AMOC_100yr.nc')
    amoc = dsmoc.AMOC.where(dsmoc.AMOC<1e10).where(dsmoc.AMOC>-1e10).where(dsmoc.z>500).where(dsmoc.lat>0).max(dim=['z','lat'])
    return amoc 
    
def detrend_dim(da, dim, deg=1):
    # detrend along a single dimension
    p = da.polyfit(dim=dim, deg=deg)
    fit = xr.polyval(da[dim], p.polyfit_coefficients)
    return da - fit

def corsst(model, experiment):
    "Computes Pearson correlation coefficient and associated p-value between SST field and AMOC strength"
    folder = '/Volumes/External/DataPlioMIP2/Data/Processed/'
    dstos = xr.open_dataset(folder+model+'/'+experiment+'/SST_annual_100yr.nc')
    
    sst = dstos.sst
    amoc = amocstrength(model, experiment)
    
    #Shorten timeseries for CESM1.2 and EC-Earth3-LR to correct for mismatch in years in time series
    if model == 'EC-Earth3-LR' and experiment == 'Eoi400':
        amoc = detrend_dim(amoc[3:], "time", deg=1)
        sst = detrend_dim(sst[:-3], "time", deg=1)
        cor = xr.corr(amoc, sst, dim='time')
        pvalue = xskillscore.pearson_r_eff_p_value(amoc, sst, dim='time', skipna=True, keep_attrs=False)
    else:
        amoc = detrend_dim(amoc, "time", deg=1)
        sst = detrend_dim(sst, "time", deg=1)
        cor = xr.corr(amoc, sst, dim='time')
        pvalue = xskillscore.pearson_r_eff_p_value(amoc, sst, dim='time', skipna=True, keep_attrs=False)
        
    return cor, pvalue

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def corsst_10rm(model, experiment):
    "Computes Pearson correlation coefficient and associated p-value between SST field and AMOC strength"
    folder = '/Volumes/External/DataPlioMIP2/Data/Processed/'
    dstos = xr.open_dataset(folder+model+'/'+experiment+'/SST_annual_100yr.nc')
    
    sst = dstos.sst
    amoc = amocstrength(model, experiment)
    
    #Shorten timeseries for CESM1.2 and EC-Earth3-LR to correct for mismatch in years in time series
    if model == 'EC-Earth3-LR' and experiment == 'Eoi400':
        amoc = detrend_dim(amoc[3:], "time", deg=1)
        amoc = amoc.rolling(time=10,center=True).mean()
        sst = detrend_dim(sst[:-3], "time", deg=1)
        sst = sst.rolling(time=10,center=True).mean()
        cor = xr.corr(amoc, sst, dim='time')
        pvalue = xskillscore.pearson_r_eff_p_value(amoc, sst, dim='time', skipna=True, keep_attrs=False)
    else:
        amoc = detrend_dim(amoc, "time", deg=1)
        amoc = amoc.rolling(time=10,center=True).mean()
        sst = detrend_dim(sst, "time", deg=1)
        sst = sst.rolling(time=10,center=True).mean()
        cor = xr.corr(amoc, sst, dim='time')
        pvalue = xskillscore.pearson_r_eff_p_value(amoc, sst, dim='time', skipna=True, keep_attrs=False)
        
    return cor, pvalue

def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]

def pvalmask(ds, CI):
    "Creates mask using p-values with value 1 where significant, 0 where not significant"
    mask_p = ds.pval.where(ds.pval<CI)
    mask_p = (mask_p/mask_p).fillna(0)
    return mask_p

def makedz(ds):
    dz = np.zeros(len(ds.z))
    dz[0] = ds.z[0].values*2
    depth = dz[0]
    for i in range(1,len(ds.z)):
        dz[i] = (ds.z[i]-depth)*2
        depth = depth+dz[i] #Depth spacing. Units: cm
    return xr.DataArray(data=dz, dims=("z"))

def difmodel(dsE280var, dsEoi400var, latmin, latmax):
    lats = np.arange(-89, 90, 0.5)
    var1 = dsE280var.interp(lat = lats, method='linear')
    var2 = dsEoi400var.interp(lat = lats, method = 'linear')
    dif = var2.where(var2.lat>latmin).where(var2.lat<latmax).mean(dim='lat')-var1.where(var1.lat>latmin).where(var1.lat<latmax).mean(dim='lat')
    return dif

def oht_amoc_anomaly(modellist, latmin, latmax):
    folder = '/Volumes/External/DataPlioMIP2/Data/Processed/'
    lats = np.arange(-89, 90, 0.5)
    difamoc = np.zeros(len(modellist))
    difov = np.zeros(len(modellist))
    difoht = np.zeros(len(modellist))

    for i in range(len(modellist)):
        model = modellist[i]
        ds1 = xr.open_dataset(folder+model+'/E280/decomOHT_100yr.nc')
        ds2 = xr.open_dataset(folder+model+'/Eoi400/decomOHT_100yr.nc')
        difov[i] = difmodel(ds1.OHTov, ds2.OHTov, latmin, latmax)

        ds1 = xr.open_dataset(folder+model+'/E280/OHT_100yr.nc')
        ds2 = xr.open_dataset(folder+model+'/Eoi400/OHT_100yr.nc')
        difoht[i] = difmodel(ds1.OHT, ds2.OHT, latmin, latmax)

        ds1 = xr.open_dataset(folder+model+'/E280/AMOC_100yr.nc')
        ds2 = xr.open_dataset(folder+model+'/Eoi400/AMOC_100yr.nc')
        amoc1 = ds1.AMOC.where(ds1.AMOC<1e10).where(ds1.AMOC>-1e10).where(ds1.lat>0).where(ds1.lat<60).where(ds1.z>500).max(dim=['lat','z'])
        amoc2 = ds2.AMOC.where(ds2.AMOC<1e10).where(ds2.AMOC>-1e10).where(ds2.lat>0).where(ds2.lat<60).where(ds2.z>500).max(dim=['lat','z'])
        difamoc[i] = amoc2 - amoc1
        
    return difamoc, difov, difoht

def merid_gradient(variable, lat_bnds, z_bot):
    weights = makedz(variable)
    N_avg = variable.where(variable.lat>lat_bnds[0]).where(variable.lat<lat_bnds[1]).where(variable.z<z_bot).weighted(weights).mean()
    S_avg = variable.where(variable.lat>lat_bnds[2]).where(variable.lat<lat_bnds[3]).where(variable.z<z_bot).weighted(weights).mean()
    gradient = N_avg-S_avg
    return(gradient)

def N_avg(variable, lat_bnds, z_bot):
    weights = makedz(variable)
    N_avg = variable.where(variable.lat>lat_bnds[0]).where(variable.lat<lat_bnds[1]).where(variable.z<z_bot).weighted(weights).mean()
    return(N_avg)

def S_avg(variable, lat_bnds, z_bot):
    weights = makedz(variable)
    S_avg = variable.where(variable.lat>lat_bnds[2]).where(variable.lat<lat_bnds[3]).where(variable.z<z_bot).weighted(weights).mean()
    return(S_avg)