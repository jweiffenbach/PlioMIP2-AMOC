import numpy as np
import xarray as xr
import xskillscore


def lon180(ds):
    "Rewrite 0-360 longitude to -180-180 for regular 1x1 grid"
    ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
    ds = ds.sortby(ds.lon)
    return ds

def lon180cesm(ds):
    "Rewrite 0-360 longitude to -180-180 for CESM models"
    ds.coords['TLONG'] = (ds.coords['TLONG'] + 180) % 360 - 180
    #ds = ds.sortby(ds.TLONG)
    return ds

def lon180ipsl(ds):
    "Rewrite 0-360 longitude to -180-180 for IPSL models"
    ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
    #ds = ds.sortby(ds.TLONG)
    return ds

def salmask(a=1):
    #Create Atlantic mask to exclude points in Mediterranean
    import regionmask
    model = 'CCSM4'
    Atlantic = np.array([[-110,30],[-75,50],[-75, 65],[-90,75],[-80,80],[-40,80], [25,80], [25,70],[10,60],[-10,35], [-10,10], [20,10], [20,-34], [-70,-34], [-40, -10], [-70,5],[-100,20]])
    region = regionmask.Regions([Atlantic])

    # define lat/ lon grid using standard 1x1 grid
    ds = xr.open_dataset('/Users/6497241/surfdrive/Documents/PlioMIP2-OHT/Data/Processed/CCSM4/E280/SST_annual_100yr.nc')
    ds = lon180(ds)

    mask_Atl_reg = region.mask(ds.lon.values, ds.lat.values)
    return mask_Atl_reg

def latweights(a=1):
    #Latitudinal weights for SST averaging using standard grid
    ds = xr.open_dataset('/Users/6497241/surfdrive/Documents/PlioMIP2-OHT/Data/Processed/CCSM4/E280/SST_annual_100yr.nc')
    weights = np.cos(np.deg2rad(ds.lat))
    weights.name = "weights"
    return weights


def makedz(ds):
    #Spacing of the vertical grid
    dz = np.zeros(len(ds.z))
    dz[0] = ds.z[0].values*2
    depth = dz[0]
    for i in range(1,len(ds.z)):
        dz[i] = (ds.z[i]-depth)*2
        depth = depth+dz[i] #Depth spacing. Units: cm
    return xr.DataArray(data=dz, dims=("z"))