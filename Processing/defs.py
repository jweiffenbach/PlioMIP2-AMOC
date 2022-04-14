import numpy as np
import xarray as xr
import xskillscore

def regularAtlmask(ds):
    
    #Create Atlantic mask on regular 1x1 grid and exclude Mediterranean
    #This mask is used for computing zonal mean Atlantic salinity/temperature
    import regionmask

    Atlantic = np.array([[-110,30],[-75,50],[-75, 65],[-90,75],[-80,80],[-40,80], [25,80], [25,70],[10,60],[-10,35], [-10,10], [20,10], [20,-34], [-65,-34], [-40, -10], [-70,5], [-100,20]])
    region = regionmask.Regions([Atlantic])

    # define lat/ lon grid
    dsgrid = ds.copy(deep=True)
    mask = region.mask(ds.longitude.values, ds.latitude.values)
    mask = mask.rename({'lon': 'longitude', 'lat': 'latitude'})
    return(mask)
