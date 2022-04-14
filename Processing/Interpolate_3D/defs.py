import numpy as np
import xarray as xr
import xskillscore

def regularAtlmask(ds):
    
    #Create Atlantic mask on regular 1x1 grid and exclude Mediterranean
    import regionmask

    Atlantic = np.array([[-110,30],[-75,50],[-75, 65],[-90,75],[-80,80],[-40,80], [25,80], [25,70],[10,60],[-10,35], [-10,10], [20,10], [20,-34], [-65,-34], [-40, -10], [-70,5], [-100,20]])
    region = regionmask.Regions([Atlantic])

    # define lat/ lon grid
    dsgrid = ds.copy(deep=True)
    mask = region.mask(ds.longitude.values, ds.latitude.values)
    mask = mask.rename({'lon': 'longitude', 'lat': 'latitude'})
    return(mask)

def lon180(ds):
    "Rewrite 0-360 longitude to -180-180"
    ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
    ds = ds.sortby(ds.lon)
    return ds

def amocstrength(model, experiment):
    "Computes AMOC strength by taking the maximum streamfunction below 500 m depth and north of 0N"
    folder = '/Users/6497241/surfdrive/Documents/PlioMIP2-OHT/Data/Processed/'
    dsmoc = xr.open_dataset(folder+model+'/'+experiment+'/AMOC_annual_100yr.nc')
    amoc = dsmoc.AMOC.where(dsmoc.z>500).where(dsmoc.lat>30).max(dim=['z','lat'])
    return amoc 

def amocstrength_100yr(model, experiment):
    folder = '/Users/6497241/surfdrive/Documents/PlioMIP2-OHT/Data/Processed/'
    dsmoc = xr.open_dataset(folder+model+'/'+experiment+'/AMOC_100yr.nc')
    amoc = dsmoc.AMOC.where(dsmoc.z>500).where(dsmoc.lat>30).max(dim=['z','lat'])
    return amoc 
    

def corsst(model, experiment):
    "Computes Pearson correlation coefficient and associated p-value between SST field and AMOC strength"
    folder = '/Users/6497241/surfdrive/Documents/PlioMIP2-OHT/Data/Processed/'
    dstos = xr.open_dataset(folder+model+'/'+experiment+'/SST_annual_100yr.nc')
    
    sst = dstos.sst
    amoc = amocstrength(model, experiment)
    
    #Shorten timeseries for CESM1.2 and EC-Earth3-LR to correct for mismatch in years in time series
    if model == 'EC-Earth3-LR' and experiment == 'Eoi400':
        cor = xr.corr(amoc[3:], sst[:-3], dim='time')
        pvalue = xskillscore.pearson_r_p_value(amoc[3:], sst[:-3], dim='time', weights=None, skipna=True, keep_attrs=False)
    else:
        cor = xr.corr(amoc, sst, dim='time')
        pvalue = xskillscore.pearson_r_p_value(amoc, sst, dim='time', weights=None, skipna=True, keep_attrs=False)
        
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