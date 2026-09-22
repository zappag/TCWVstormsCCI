import argparse
import os
import xarray as xr
import numpy as np
from cdo import Cdo, CDOException
cdo = Cdo()
import matplotlib.pyplot as plt
import pandas as pd
import glob
import gc
from datetime import datetime

def solar_noon_utc_hour(longitude):
    """
    Returns the UTC hour at which solar noon occurs for a given longitude.
    Parameters:
        longitude (float): Longitude in degrees. East is positive, west is negative.
    Returns:
        float: UTC hour (in 24-hour format) when the local solar noon occurs.
    """
    longitude = ((longitude + 180) % 360) - 180
    return 12 - (longitude / 15.0)

def create_time_24(xarray_dataset):
    """
    Create a new time coordinate in the xarray dataset that represents 24-hour time.
    Parameters:
        xarray_dataset (xarray.Dataset): The input dataset with a 'time' coordinate.
    Returns:
        xarray.Dataset: The dataset with a new 'time_24' coordinate.
    """
    xarray_newtime = xarray_dataset.copy()
    times_0 = xarray_dataset['time'].values[xarray_dataset['time'].dt.hour == 0]
    print(times_0)
    new_times = []
    for t in times_0[:]:
        prev_day = np.datetime64(t) - np.timedelta64(1, 'D')
        print(prev_day)
        new_time = np.datetime64(str(prev_day)[:10] + 'T23:30')
        print(new_time)
        new_times.append(new_time)
    data_at_0 = xarray_dataset.sel(time=xarray_dataset['time'].dt.hour == 0)
    data_at_0 = data_at_0.assign_coords(time=('time', new_times))
    xarray_newtime = xr.concat([xarray_dataset, data_at_0], dim='time')
    xarray_newtime = xarray_newtime.sortby('time')
    return xarray_newtime

def save_seasons_to_disk(ds, outdir, prefix="season",varname='TCWV'):
    """
    Split an xarray Dataset/DataArray into meteorological seasons and save each to disk.
    Only complete seasons are saved.
    Parameters:
        ds (xr.Dataset or xr.DataArray): Input data with a 'time' coordinate.
        outdir (str): Output directory.
        prefix (str): Prefix for output files.
    """
    season_months = {
        "DJF": [12, 1, 2],
        "MAM": [3, 4, 5],
        "JJA": [6, 7, 8],
        "SON": [9, 10, 11]
    }

    time = ds['time'].to_index()
    years = time.year
    months = time.month
    season_labels = []
    season_years = []
    for t in time:
        y, m = t.year, t.month
        if m == 12:
            season = "DJF"
            season_year = y
        elif m in [1, 2]:
            season = "DJF"
            season_year = y - 1
        elif m in [3, 4, 5]:
            season = "MAM"
            season_year = y
        elif m in [6, 7, 8]:
            season = "JJA"
            season_year = y
        elif m in [9, 10, 11]:
            season = "SON"
            season_year = y
        season_labels.append(season)
        season_years.append(season_year)
    ds = ds.assign_coords(season=("time", season_labels), season_year=("time", season_years))
    import pandas as pd
    df = pd.DataFrame({"season": season_labels, "season_year": season_years}, index=time)
    grouped = df.groupby(["season_year", "season"])
    complete_seasons = []
    for (y, s), group in grouped:
        months_in_season = season_months[s]
        if s == "DJF":
            months_needed = [12, 1, 2]
            years_needed = [y, y+1, y+1]
            months_present = [(t.year, t.month) for t in group.index]
            if all((yy, mm) in months_present for yy, mm in zip(years_needed, months_needed)):
                complete_seasons.append((y, s))
        else:
            if set(months_in_season).issubset(set(group.index.month)):
                complete_seasons.append((y, s))
    for y, s in complete_seasons:
        season_ds = ds.sel(time=(ds.season == s) & (ds.season_year == y))
        if season_ds.time.size >= 3*28:
            outname = f"{prefix}_{y}_{s}.nc"
            season_ds=season_ds.drop_vars(['season', 'season_year'])
            nlat= season_ds.sizes['lat']
            nlon= season_ds.sizes['lon']
            print(nlat,nlon)
            encoding = {varname: {'chunksizes': (1, nlat, nlon)}}
            season_ds.to_netcdf(os.path.join(outdir, outname), mode='w', encoding=encoding, format='NETCDF4', engine='netcdf4')
            print(f"Saved {outname}")

def file_has_date_in_range(ncfile, date_start, date_end):
    """
    Returns True if any date in the NetCDF file is within [date_start, date_end].
    """
    import subprocess
    from datetime import datetime
    try:
        result = subprocess.run(
            ["cdo", "showdate", ncfile],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True
        )
        dates = result.stdout.strip().split()
        for d in dates:
            if date_start <= d <= date_end:
                return True
        return False
    except Exception as e:
        print(f"Error running cdo showdate on {ncfile}: {e}")
        return False

def compute_esa_masks(ds_esa,version):
    if version == '3':
        # compute masks based on land type and quality flag
        land_coast_ice_mask = ds_esa['surface_type_flag'].isin([0, 2, 4, 5, 6, 7])
        ocean_mask = ds_esa['surface_type_flag'].isin([1, 3])
        land_mask = ds_esa['surface_type_flag'].isin([0, 2, 6])
        land_clear_mask = ds_esa['surface_type_flag'].isin([0])
        valid_data_mask= ds_esa['tcwv_quality_flag'] == 0

    elif version == '4':
        land_coast_ice_mask = ds_esa['surface_type_flag'].isin([0, 2, 3, 4, 5, 6])
        ocean_mask = ds_esa['surface_type_flag'].isin([1])
        land_mask = ds_esa['surface_type_flag'].isin([0])
        land_clear_mask = ds_esa['atmospheric_conditions_flag'].isin([1])
        valid_data_mask= ds_esa['tcwv_quality_flag'] == 0

    return {
        "land_coast_ice_mask": land_coast_ice_mask,
        "ocean_mask": ocean_mask,
        "land_mask": land_mask,
        "land_clear_mask": land_clear_mask,
        "valid_data_mask": valid_data_mask
    }

def load_esa_and_masks(date_start, date_end, version):
    """
    Load ESA CCI WV data for the given date range and return (ds_esa, masks_dict).
    """
    if version == '3':
        #esadir='/work/users/clima/zappa/ESACCI-WV/cdr2/raw/'
        esadir='/mnt/naszappa/ESACCI-WV-draftpaper/cdr2/raw/'
    elif version == '4':
        esadir='/mnt/naszappa/ESACCI-WV_v4/cdr2/raw2/'

    ds_esa = xr.open_mfdataset(f"{esadir}/HTWdm*.nc", combine='by_coords')
    ds_esa = ds_esa.sel(time=slice(date_start, date_end))
    ds_esa['time'] = ds_esa['time'].dt.floor('D') + np.timedelta64(12, 'h')
    ds_esa = ds_esa.assign_coords(lon=(ds_esa['lon'] + 360) % 360)
    ds_esa = ds_esa.sortby('lon')
    ds_esa = ds_esa.assign_coords(
        time=ds_esa['time'],
        lon=ds_esa['lon'],
        lat=ds_esa['lat'],
    )
    # keep data as is where mask true, and set rest (north pole) to 7 (ice)
    if version == '3':
        ds_esa['surface_type_flag'] = ds_esa['surface_type_flag'].where(np.abs(ds_esa['lat']) < 88.5, 7)
    elif version == '4':
        ds_esa['surface_type_flag'] = ds_esa['surface_type_flag'].where(np.abs(ds_esa['lat']) < 88.5, 6)
    
    masks = compute_esa_masks(ds_esa,version)

    return ds_esa, masks

def load_esa_and_masks3(date_start, date_end, version):

    if version == '3':
        #esadir='/work/users/clima/zappa/ESACCI-WV/cdr2/raw/'
        esadir='/mnt/naszappa/ESACCI-WV/cdr2/raw/'
    elif version == '4':
        esadir='/mnt/naszappa/ESACCI-WV_v4/cdr2/raw2/'

    # 1️⃣ List all files
    all_files = sorted(glob.glob(f"{esadir}/HTWdm*.nc"))

    def extract_time_from_filename(path):
        basename = os.path.basename(path)
        date_str = basename[5:13]  # '20171214'
        return pd.to_datetime(date_str, format='%Y%m%d')


    # 3️⃣ Filter files by desired time range
    times = [extract_time_from_filename(f) for f in all_files]
    file_df = pd.DataFrame({"file": all_files, "time": times})
    sel_files = file_df[
        (file_df["time"] >= pd.to_datetime(date_start)) &
        (file_df["time"] <= pd.to_datetime(date_end))
    ]["file"].tolist()

    # 4️⃣ Open only relevant files
    ds_esa = (
        xr.open_mfdataset(sel_files, combine="by_coords", parallel=False)
        .assign_coords(
            lon=lambda ds: (ds.lon + 360) % 360,
            time=lambda ds: ds.time.dt.floor("D") + np.timedelta64(12, "h"),
        )
        .sortby("lon")     
    )

    ds_esa = ds_esa.assign_coords(
        time=ds_esa['time'],
        lon=ds_esa['lon'],
        lat=ds_esa['lat'],
    )

    # keep data as is where mask true, and set rest (north pole) to 7 (ice)
    if version == '3':
        ds_esa['surface_type_flag'] = ds_esa['surface_type_flag'].where(np.abs(ds_esa['lat']) < 88.5, 7)
    elif version == '4':
        ds_esa['surface_type_flag'] = ds_esa['surface_type_flag'].where(np.abs(ds_esa['lat']) < 88.5, 6)

    masks = compute_esa_masks(ds_esa,version)

    return ds_esa, masks


def compute_low_data_points(ds_esa, threshd=10):
    """
    Compute low data points where TCWV is below a certain threshold.
    Parameters:
        ds_esa (xarray.Dataset): ESA dataset with 'tcwv' variable.
        threshold (float): Threshold value for low data points.
    Returns:
        xarray.DataArray: Boolean mask of low data points.
    """
    # select valid data points based on quality flag
    valid_data_mask= ds_esa['tcwv_quality_flag'] == 0

    # Load it into memory
    valid_data_mask = valid_data_mask.load()

    # sum grouped by season (DJF, MAM, JJA, SON)
    time_sum = valid_data_mask.groupby('time.season').sum(dim='time')

    # number of days in each season
    season_counts = ds_esa['time'].groupby('time.season').count()

    # select points with less than 10% of data in each season
    threshd= 10
    low_data_points = time_sum < threshd/100 * season_counts

    return low_data_points


if __name__ == "__main__":
    # script code here

    # script requires as input reanalysis data - ideally hourly - interpolated to ESA CCI grid

    parser = argparse.ArgumentParser(description="Process TCWV data.")
    parser.add_argument('--fdata', type=str, default='ERA5',
                    help="Data type to be added (e.g., ERA5, MERRA2, ESA)")
    parser.add_argument('--year', type=str, default=None,
                    help="")
    parser.add_argument('--version', type=str, default=None,
                    help="")
    parser.add_argument('--fullrun', action='store_true',
                    help="If set, process the full dataset.")
    parser.add_argument('--crossmask', action='store_true', default=False,
                    help="If true, runs v4 data with v3 masks (for sensitivity test).")

    # Reanalysis Data
    dataset= parser.parse_args().fdata
    tyear= parser.parse_args().year
    version= parser.parse_args().version
    fullrun= parser.parse_args().fullrun
    crossmask= parser.parse_args().crossmask

    if version == '3':
        basedir='/mnt/naszappa/ESACCI-WV/'
    elif version == '4':
        basedir='/mnt/naszappa/ESACCI-WV_v4/'

    if crossmask == True:
        crossext='cross'
        if version == '4':
            version_mask = '3'
        elif version == '3':
            version_mask = '4'
    else:
        crossext=''    
        version_mask = version

    if dataset == 'ERA5':
        varname='TCWV'
        rean_dir='/mnt/naszappa/ERA5/store/1hr/'
        rean_file='ERA5_total_column_water_vapour_1hr_full_sfc_*.nc'
        outdir_hourly= basedir + '/ERA5/tcwv/hourly_50km/'
        outdir_6h_full= basedir + '/ERA5/tcwv/6h_v2/'
        outdir_12UTC_full= basedir + '/ERA5/tcwv/12UTC_v2/'
        outdir_6h= basedir + '/ERA5/tcwv/6h_50km_v2/'
        outdir_day= basedir + '/ERA5/tcwv/day_50km_v2/'
        outdir_day_masked= basedir + '/ERA5/tcwv/day_50km_masked_v2/'
    elif dataset == 'ESA':
        if version == '3':
            outdir_day= '/work/users/clima/zappa/ESACCI-WV/cdr2/day_v2/'
            outdir_day_masked= f'/work/users/clima/zappa/ESACCI-WV/cdr2/day_masked_v2{crossext}/'
            outdir_masks= basedir + '/cdr2/masks/'
            esadir=basedir + '/cdr2/raw/'
        elif version == '4':
            outdir_day= basedir + '/cdr2/day_v2/'
            outdir_day_masked= basedir + f'/cdr2/day_masked_v2{crossext}/'
            outdir_masks= basedir + '/cdr2/masks/'
            esadir=basedir + '/cdr2/raw2/'
    # elif dataset == 'MERRA2-nas':
    #     varname='TQV'
    #     rean_dir='/home/zappa/naszappa/MERRA2/TQV/hourly/'
    #     rean_file='MERRA2_TQV_1hr*'
    #     outdir_hourly= '/home/zappa/naszappa/ESACCI-WV/MERRA2/tcwv/hourly_50km/'
    #     outdir_6h_full= '/home/zappa/naszappa/ESACCI-WV/MERRA2/tcwv/6h_v2/'
    #     outdir_12UTC_full= '/home/zappa/naszappa/ESACCI-WV/MERRA2/tcwv/12UTC_v2/'
    #     outdir_6h= '/home/zappa/naszappa/ESACCI-WV/MERRA2/tcwv/6h_50km_v2/'
    #     outdir_day= '/home/zappa/naszappa/ESACCI-WV/MERRA2/tcwv/day_50km_v2/'
    #     outdir_day_masked= '/home/zappa/naszappa/ESACCI-WV/MERRA2/tcwv/day_50km_masked_v2/'
    elif dataset == 'MERRA2':
        varname='TQV'
        rean_dir='/mnt/naszappa/MERRA2/TQV/hourly/'
        rean_file='MERRA2_*.nc'
        outdir_hourly= basedir + '/MERRA2/tcwv/hourly_50km/'
        outdir_6h_full= basedir + '/MERRA2/tcwv/6h_v2/'
        outdir_12UTC_full= basedir + '/MERRA2/tcwv/12UTC_v2/'
        outdir_6h= basedir + '/MERRA2/tcwv/6h_50km_v2/'
        outdir_day= basedir + '/MERRA2/tcwv/day_50km_v2/'
        outdir_day_masked= basedir + '/MERRA2/tcwv/day_50km_masked_v2/'
    elif dataset == 'JRA3Q':
        varname='tciwv-col-an-gauss'
        rean_dir='/mnt/naszappa/JRA3Q/tciwv/6hourly/'
        rean_file='jra3q*tciwv*.nc'
        outdir_hourly= basedir + '/JRA3Q/tcwv/hourly_50km/'
        outdir_6h_full= basedir + '/JRA3Q/tcwv/6h_v2/'
        outdir_12UTC_full= basedir + '/JRA3Q/tcwv/12UTC_v2/'
        outdir_6h= basedir + '/JRA3Q/tcwv/6h_50km_v2/'
        outdir_day= basedir + '/JRA3Q/tcwv/day_50km_v2/'
        outdir_day_masked= basedir + '/JRA3Q/tcwv/day_50km_masked_v2/'
    else:
        raise Exception(f"Dataset {dataset} not recognized. Choose from ERA5, MERRA2, ESA, JRA3Q.")

    # if era5 data is used, create output directories
    if dataset == 'ERA5' or dataset == 'MERRA2' or dataset == 'MERRA2-work' or dataset == 'JRA3Q':
        os.system('mkdir -p ' + outdir_hourly)
        os.system('mkdir -p ' + outdir_6h_full)
        os.system('mkdir -p ' + outdir_6h)
        os.system('mkdir -p ' + outdir_12UTC_full)
        os.system('mkdir -p ' + outdir_day)
        os.system('mkdir -p ' + outdir_day_masked)
    elif dataset == 'ESA':
        os.system('mkdir -p ' + outdir_day)
        os.system('mkdir -p ' + outdir_day_masked)


    # ESA CCI WV Data
    # time period to process

    if tyear == None:
        date_start= '2002-12-01'
        if version == '3':
            date_end= '2017-11-30'
        elif version == '4' and crossmask == False:
            date_end= '2022-11-30'
        elif version == '4' and crossmask == True:
            date_end= '2017-11-30'
    else:
        tyearp= str(int(tyear)+1)
        print(tyear,tyearp)
        date_start= f"{tyear}-12-01" 
        date_end= f"{tyearp}-11-30"
        #date_start= '2007-12-01'
    #date_end= '2008-02-28'

    # Read ESA CCI WV data and prepare masks
    print(datetime.now().strftime("%H:%M:%S"))
    print('Reading ESA data...')
    #ds_esa, masks = load_esa_and_masks3(date_start, date_end, version)  # ok with raw/ filenames
    ds_esa, masks = load_esa_and_masks(date_start, date_end, version)   # test with yearly batches (raw2/)

    if crossmask == True:
        print(f"Using version {version_mask} masks for version {version} data (crossmask=True).")
        _ , masks = load_esa_and_masks(date_start, date_end, version_mask)

    # # check for duplicates
    # time = ds_esa['time'].values
    # import pandas as pd
    # t = pd.DatetimeIndex(time)
    # print(t[t.duplicated(keep=False)])  # shows the actual duplicate timestamps

    # tm = pd.DatetimeIndex(masks['valid_data_mask']['time'].values)
    # print(tm[tm.duplicated(keep=False)])
    # raise Exception("Finished loading ESA data and masks, stopping here for inspection.")

    if dataset == 'ESA':
        # save low data points
        threshd= 10
        low_data_points = compute_low_data_points(ds_esa, threshd=threshd)
        os.system('mkdir -p ' + outdir_masks)
        low_data_points.to_netcdf(os.path.join(outdir_masks, f'low_data_points_{threshd}p_{date_start}_{date_end}.nc'))
        print('Saved low data points mask to disk.')

        # split into seasons
        ds_esa = ds_esa.rename({'tcwv': 'TCWV'})
        save_seasons_to_disk(ds_esa['TCWV'], outdir_day, 'tcwv_ESA')
        
        # apply valid ma
        ds_esa_valid= ds_esa.where(masks['valid_data_mask'].compute())
        save_seasons_to_disk(ds_esa_valid['TCWV'], outdir_day_masked, 'tcwv_ESA')

        ds_esa_land_coast_ice_masked = ds_esa_valid.where(masks['land_coast_ice_mask'].compute())
        save_seasons_to_disk(ds_esa_land_coast_ice_masked['TCWV'], outdir_day_masked, 'tcwv_ESA_LandCoastIce')

        ds_esa_land_clear_masked = ds_esa_valid.where(masks['land_clear_mask'].compute())
        save_seasons_to_disk(ds_esa_land_clear_masked['TCWV'], outdir_day_masked, 'tcwv_ESA_LandClear')

        ds_esa_ocean_masked = ds_esa_valid.where(masks['ocean_mask'].compute())
        save_seasons_to_disk(ds_esa_ocean_masked['TCWV'], outdir_day_masked, 'tcwv_ESA_Ocean')
        raise Exception("ESA processing complete, stopping here.") 


    # clean memory 
    del ds_esa
    gc.collect()

    if fullrun:

        # Convert reanalyses files (as they are) to ESA spatial resolution for later processing
        print(datetime.now().strftime("%H:%M:%S"))
        print('Remap reanalysis to ESA spatial resolution...')
        grid_esa=basedir + "/grid_esa.txt"
        for ff in glob.glob(os.path.join(rean_dir,rean_file)):
            if not file_has_date_in_range(ff, date_start, date_end):
                print(f"Skipping {ff} as it does not contain data in the specified date range.")
                continue
            ffname= os.path.basename(ff)
            outfile= os.path.join(outdir_hourly, ffname)
            if os.path.exists(outfile):
                print("File already exists: " + outfile)
                continue
            #cdo.remapycon(grid_esa, input=ff, output=outfile, options='-f nc4 -z zip')
            cdo.remapycon(grid_esa, input=f"-selname,{varname} " + ff, output=outfile)


        selected_hours_6h = [0, 6, 12, 18]
        # # 6 hourly full resolution split by season (6h)
        print(datetime.now().strftime("%H:%M:%S"))
        print('Extract 6 hours samples split be season...')
        ds_rean_full = xr.open_mfdataset(os.path.join(rean_dir, rean_file), combine='by_coords').sel(time=slice(date_start, date_end))
        if 'TCWV' not in ds_rean_full:
            ds_rean_full = ds_rean_full.rename({varname: 'TCWV'})
        time_mask = (ds_rean_full['time'].dt.hour.isin(selected_hours_6h))
        ds_rean_full = ds_rean_full['TCWV']
        if dataset == 'MERRA2':
            # drop latitude values in magnitude smaller than 0.0001 (drop wierd equator)
            ds_rean_full = ds_rean_full.assign_coords(lon=(ds_rean_full['lon'] + 360) % 360)
            ds_rean_full = ds_rean_full.where(np.abs(ds_rean_full['lat']) > 0.0001, drop=True)
            ds_rean_full = ds_rean_full.sortby('lon')
        ds_rean_full_6h_samples = ds_rean_full.sel(time=time_mask).compute()
        ds_rean_full_6h_samples = ds_rean_full_6h_samples.assign_coords(
            time=ds_rean_full_6h_samples['time'],
            lon=ds_rean_full_6h_samples['lon'],
            lat=ds_rean_full_6h_samples['lat'],
        )
        save_seasons_to_disk(ds_rean_full_6h_samples, outdir_6h_full, f'tcwv_{dataset}')


        # # for radial composite and comparison at fixed storm maximum (12UTC)
        print(datetime.now().strftime("%H:%M:%S"))
        print('Extract 12 UTC samples split be season...')                                                                             
        selected_hours_12utc = [12]
        time_mask = (ds_rean_full['time'].dt.hour.isin(selected_hours_12utc))
        ds_rean_full_12UTC_samples = ds_rean_full.sel(time=time_mask).compute()
        ds_rean_full_12UTC_samples = ds_rean_full_12UTC_samples.assign_coords(
            time=ds_rean_full_12UTC_samples['time'],
            lon=ds_rean_full_12UTC_samples['lon'],
            lat=ds_rean_full_12UTC_samples['lat'],
        )
        save_seasons_to_disk(ds_rean_full_12UTC_samples, outdir_12UTC_full, f'tcwv_{dataset}')


        # 6hourly ESA RESOUTION split by season (6h_50km)
        print(datetime.now().strftime("%H:%M:%S"))
        print('Extract 6 hourly samples split be season at ESA resolution...') 
        ds_rean_50km = xr.open_mfdataset(os.path.join(outdir_hourly, rean_file), combine='by_coords').sel(time=slice(date_start, date_end))
        if 'TCWV' not in ds_rean_50km:
            ds_rean_50km = ds_rean_50km.rename({varname: 'TCWV'})

        time_mask = (ds_rean_50km['time'].dt.hour.isin(selected_hours_6h))
        ds_rean_50km = ds_rean_50km['TCWV'].compute()
        ds_rean_50km_6h_samples = ds_rean_50km.sel(time=time_mask).compute()
        ds_rean_50km_6h_samples = ds_rean_50km_6h_samples.assign_coords(
            time=ds_rean_50km_6h_samples['time'],
            lon=ds_rean_50km_6h_samples['lon'],
            lat=ds_rean_50km_6h_samples['lat'],
        )
        save_seasons_to_disk(ds_rean_50km_6h_samples, outdir_6h, f'tcwv_{dataset}')

    ## clean memory
    #del ds_rean_50km_6h_samples, ds_rean_full_6h_samples, ds_rean_full_12UTC_samples
    #gc.collect()

    # rename: hourly 50km dataset to ds_rean
    print(datetime.now().strftime("%H:%M:%S"))
    # print('Interpolate JRA3Q to hourly...')

    if fullrun:
        ds_rean=ds_rean_50km
    else:
        ds_rean_50km = xr.open_mfdataset(os.path.join(outdir_hourly, rean_file), combine='by_coords').sel(time=slice(date_start, date_end))
        if 'TCWV' not in ds_rean_50km:
            ds_rean_50km = ds_rean_50km.rename({varname: 'TCWV'})
        ds_rean=ds_rean_50km
    
    # if dataset == 'JRA3Q':
    #     # Assume ds is your 6-hourly xarray Dataset/DataArray with a 'time' coordinate
    #     # Create a new hourly time index covering the same range
    #     new_time = pd.date_range(ds_rean['time'].min().item(), ds_rean['time'].max().item(), freq='1h')
        
    #     # Interpolate to hourly using nearest neighbor
    #     ds_rean = ds_rean.interp(time=new_time, method="nearest").load()

    print(datetime.now().strftime("%H:%M:%S"))
    print('Prepare for ESA samplingy...')    
    time_12utc = ds_rean['time'].sel(time=ds_rean['time'].dt.hour == 12)
    ds_rean_land2=ds_rean.copy(deep=True)
    ds_rean_ocean=ds_rean.copy(deep=True)

    # Compute the solar noon hour for each longitude (vectorized)
    print(datetime.now().strftime("%H:%M:%S"))
    print('Compute solar noon timestep...')
    solar_noon_hours = 12 - ((ds_rean['lon'] + 180) % 360 - 180) / 15.0 

    if dataset == 'JRA3Q': # 6 hourly
        solar_noon_hours_rounded = (6 * np.round(solar_noon_hours / 6)) % 24
        solar_noon_hours_rounded = solar_noon_hours_rounded.astype(int)
    else:
        solar_noon_hours_rounded = np.round(solar_noon_hours).astype(int) % 24

    #solar_noon_hours_rounded = np.round(solar_noon_hours).astype(int) % 24
    lon = ds_rean['lon']
    time = ds_rean['time']
    hour_2d = xr.DataArray(
        np.broadcast_to(solar_noon_hours_rounded.values[:, np.newaxis], (len(lon), len(time))),
        dims=('lon', 'time'),
        coords={'lon': lon, 'time': time}
    )
    time_hours = ds_rean['time'].dt.hour
    time_hours_2d = xr.DataArray(
        np.broadcast_to(time_hours.values, (len(lon), len(time))),
        dims=('lon', 'time'),
        coords={'lon': lon, 'time': time}
    )
    mask = (time_hours_2d == hour_2d)
    ds_rean_land2 = ds_rean.where(mask)
    ds_rean_land2['time'] = ds_rean_land2['time'].dt.floor('D') + np.timedelta64(12, 'h')
    ds_rean_land2 = ds_rean_land2.groupby('time').sum('time')

    ds_rean_land2 = ds_rean_land2.assign_coords(
        time=time_12utc,
        lon=ds_rean_land2['lon'],
        lat=ds_rean_land2['lat'],
    )

    # compute daily-mean from hourly data
    print(datetime.now().strftime("%H:%M:%S"))
    print('Compute daily means...')
    ds_rean_ocean= ds_rean.resample(time='1D').mean()
    ds_rean_ocean['time'] = ds_rean_ocean['time'].dt.floor('D') + np.timedelta64(12, 'h')
    ds_rean_ocean=ds_rean_ocean
    ds_rean_ocean = ds_rean_ocean.transpose('time', 'lat', 'lon')
    ds_rean_ocean = ds_rean_ocean.assign_coords(
        time=time_12utc,
        lon=ds_rean_ocean['lon'],
        lat=ds_rean_ocean['lat'],
    )
    #ds_rean_ocean.attrs = {}

    # # compute masks based on land type and quality flag
    # land_coast_ice_mask = ds_esa['surface_type_flag'].isin([0, 2, 4, 5, 6, 7])
    # ocean_mask = ds_esa['surface_type_flag'].isin([1, 3])
    # land_mask = ds_esa['surface_type_flag'].isin([0, 2, 6])
    # valid_data_mask= ds_esa['tcwv_quality_flag'] == 0

    # # Plot both land and ocean masks for the same time on the same figure
    # fig, axes = plt.subplots(1, 3, figsize=(12, 5))
    # land_mask.sel(time='2002-12-01').plot(ax=axes[0])
    # axes[0].set_title('Land Mask (2002-12-01)')
    # ocean_mask.sel(time='2002-12-01').plot(ax=axes[1])
    # axes[1].set_title('Ocean Mask (2002-12-01)')
    # land_coast_ice_mask.sel(time='2002-12-01').plot(ax=axes[2])
    # axes[1].set_title('Land Coast Ice Mask (2002-12-01)')
    # plt.tight_layout()
    # plt.show()


    ds_rean_land_only= ds_rean_land2.where(masks['land_mask'].compute())
    ds_rean_land_clear_only= ds_rean_land2.where(masks['land_clear_mask'].compute())
    ds_rean_land_coast_ice_only= ds_rean_land2.where(masks['land_coast_ice_mask'].compute())
    ds_rean_ocean_only= ds_rean_ocean.where(masks['ocean_mask'].compute())
    ds_rean_combined = ds_rean_land_coast_ice_only.fillna(0) + ds_rean_ocean_only.fillna(0)

    # clean memory
    del ds_rean_land2, ds_rean_ocean
    gc.collect()

    # fig, axes = plt.subplots(2, 2, figsize=(12, 5))
    # ds_rean_land_coast_ice_only.isel(time=0).plot(ax=axes[0,0])
    # ds_rean_land_only.isel(time=0).plot(ax=axes[0,1])
    # ds_rean_ocean_only.isel(time=0).plot(ax=axes[1,0])
    # ds_rean_combined.isel(time=0).plot(ax=axes[1,1])


    print(datetime.now().strftime("%H:%M:%S"))
    print('Compute and save land/ocean split data...')
    print(ds_rean_land_only)
    
    save_seasons_to_disk(ds_rean_land_only, outdir_day, prefix=f"tcwv_{dataset}_Land")
    save_seasons_to_disk(ds_rean_land_clear_only, outdir_day, prefix=f"tcwv_{dataset}_LandClear")
    save_seasons_to_disk(ds_rean_land_coast_ice_only, outdir_day, prefix=f"tcwv_{dataset}_LandCoastIce")
    save_seasons_to_disk(ds_rean_ocean_only, outdir_day, prefix=f"tcwv_{dataset}_Ocean")
    save_seasons_to_disk(ds_rean_combined, outdir_day, prefix=f"tcwv_{dataset}")

    ds_rean_land_masked= ds_rean_land_only.where(masks['valid_data_mask'].compute())
    ds_rean_land_clear_masked= ds_rean_land_clear_only.where(masks['valid_data_mask'].compute())
    ds_rean_land_coast_ice_masked= ds_rean_land_coast_ice_only.where(masks['valid_data_mask'].compute())
    ds_rean_ocean_masked= ds_rean_ocean_only.where(masks['valid_data_mask'].compute())
    ds_rean_combined_masked= ds_rean_combined.where(masks['valid_data_mask'].compute())

    print(ds_rean_land_masked)

    # clean memory
    del ds_rean_land_only, ds_rean_land_coast_ice_only, ds_rean_ocean_only, ds_rean_combined
    gc.collect()

    print(datetime.now().strftime("%H:%M:%S"))
    print('Compute and save masked data...')
    save_seasons_to_disk(ds_rean_land_masked, outdir_day_masked, prefix=f"tcwv_{dataset}_Land")
    save_seasons_to_disk(ds_rean_land_clear_masked, outdir_day_masked, prefix=f"tcwv_{dataset}_LandClear")
    save_seasons_to_disk(ds_rean_land_coast_ice_masked, outdir_day_masked, prefix=f"tcwv_{dataset}_LandCoastIce")
    save_seasons_to_disk(ds_rean_ocean_masked, outdir_day_masked, prefix=f"tcwv_{dataset}_Ocean")
    save_seasons_to_disk(ds_rean_combined_masked, outdir_day_masked, prefix=f"tcwv_{dataset}")
