import track_wrapper
import os
import argparse


# Parse command-line arguments
parser = argparse.ArgumentParser(description="Process TCWV data.")
parser.add_argument('--fdata', type=str, default='ESA',
                    help="Data type to be added (e.g., ERA5, MERRA2, ESA, ESAcross)")
parser.add_argument('--fres', type=str, default='',
                    help="Data type to be added (e.g.,'', 50km, 50km-day, 50km-day-masked")
parser.add_argument('--fsubmask', type=str, default='',
                    help="Submask for data (e.g., '', 'LandCoastIce', 'Ocean'). Use '' for no submask.")
parser.add_argument('--radius', type=str, default='5',
                    help="Radius for mean/max/min field (default: 5).")
parser.add_argument('--dosubsample', action='store_true', help="Run time subsampling step")
parser.add_argument('--doaddtcwv', action='store_true', help="Run addtcwv step")
parser.add_argument('--dostats', action='store_true', help="Run stats step")
parser.add_argument('--doradial', action='store_true', help="Run radial step")
parser.add_argument('--testing', action='store_true', default=False,
                    help="Run in testing mode (default: False). If True, only processes one season and year.")
parser.add_argument('--version', type=str, default='3')
parser.add_argument('--feature', type=str, default='cyc', help="cyc/anticyc")
#parser.add_argument('--crossmask', action='store_true', default=False, help="True/False")


args = parser.parse_args()

# CONFIGURATION 
# seasons and years to be processed
if args.testing == False:
    aseas=('MAM','JJA','SON','DJF')
    if args.version == '3':
        seasdict={'MAM':(2003,2017),'JJA':(2003,2017),'SON':(2003,2017),'DJF':(2002,2016)}
    elif args.version == '4' and args.fdata == 'ESAcross':
        seasdict={'MAM':(2003,2017),'JJA':(2003,2017),'SON':(2003,2017),'DJF':(2002,2016)}
    elif args.version == '4' :
        seasdict={'MAM':(2003,2022),'JJA':(2003,2022),'SON':(2003,2022),'DJF':(2002,2021)}
    else:
        raise Exception('No acceptable version option')
elif args.testing == True:
    aseas=('DJF',)
    seasdict={'MAM':(2006,2006),'JJA':(2003,2017),'SON':(2003,2003),'DJF':(2002,2002)}

rd=args.radius # radius for the mean field
rdd=float(rd)

mask1=args.fsubmask
if mask1 != '':
    mask1='_'+mask1
res=args.fres
if res != '':
    res1='-'+res
else:
    res1=''

# operations to be performed
run_time_subsampling=args.dosubsample
run_addtcwv=args.doaddtcwv
run_stats=args.dostats
run_radial=args.doradial

run_ffi=False
run_ff2d=False
if args.feature == 'cyc':
    fftype0='ff_trs_neg'
    run_ffi=True
    print('Analysing cyclonic features')
elif args.feature == 'anticyc':
    run_ffi=True
    fftype0='ff_trs_pos'
    print('Analysing anticyclonic features')
elif args.feature == 'anticyc2d':
    run_ff2d=True
    fftype0='tr_trs_pos'
    print('Analysing all anticyclonic features')
elif args.feature == 'cyc2d':
    run_ff2d=True
    fftype0='tr_trs_neg'
    print('Analysing all cyclonic features')
else:
    raise Exception("wrong feature type")

if args.version == '3':
    basedir='/mnt/naszappa/ESACCI-WV'
    print('Using ESACCI version 3')
elif args.version == '4':
    basedir='/mnt/naszappa/ESACCI-WV_v4'
    print('Using ESACCI version 4')

# data to be added
fdataset = args.fdata  # Use the passed argument
fdata = args.fdata + res1
#fdata='ERA5-50km-day-masked' # ERA5-50km, ERA5-50km-day, ERA5-50km-day-masked, ESA

if fdataset not in ['ERA5','MERRA2','ESA','JRA3Q','ESAcross']:
    raise Exception(f"{fdataset} not implemented yet")

if res not in ['','50km','50km-day','50km-day-masked']:
    raise Exception(f"{res1} not implemented yet")

if fdata not in [f"{fdataset}-50km-day",f"{fdataset}-50km-day-masked","ESA",'ESAcross'] and mask1 != '':
    raise Exception(f"Submask {mask1} only applicable for ERA5-50km-day, ERA5-50km-day-masked and ESA data")

# path to directories containing data to be added to the tracks
fmiss=True # flag for missing values (NOTE: saying that missing values are present, cause general search is otherwise not working)

# # path to directory containing tracks
# if args.crossmask == False:
#     crossext=''
#     trackdir=f"{basedir}/era5_tracks_v3/"
# elif args.crossmask == True:
#     crossext='cross'
#     trackdir=f"{basedir}/era5_tracks_v3cross/"
#

trackdir=f"{basedir}/era5_tracks_v3/"

if fdata=='ESA':
    adddatadir=f"{basedir}/cdr2/day_masked_v2/"
elif fdata=='ESAcross':
    adddatadir=f"{basedir}/cdr2/day_masked_v2cross/"
elif fdata=='ESA-orig':
    adddatadir=f"{basedir}/cdr2/day_v2/"
elif res1=='':
    adddatadir=f"{basedir}/{fdataset}/tcwv/6h_v2/"
elif res1=='-50km':
    adddatadir=f"{basedir}/{fdataset}/tcwv/6h_50km_v2/"
elif res1=='-50km-day':
    adddatadir=f"{basedir}/{fdataset}/tcwv/day_50km_v2/"
elif res1=='-50km-day-masked':
    adddatadir=f"{basedir}/{fdataset}/tcwv/day_50km_masked_v2/"


## PROCESSING STARTS HERE: 
# loop over years and seasons
for season in aseas:

    # season=aseas[s]
    sy=seasdict[season][0]
    ly=seasdict[season][1]
    for y in range(sy,ly+1):
        print('Processing year '+str(y)+' season '+season)

        trackdiry=trackdir+'/'+season+'/msl/NH_ERA5_msl_6hr_'+str(y)+'_'+season
        trackfile=trackdiry+'/'+fftype0
        print('Track file to be analysed: '+trackfile)

        # unzip
        if os.path.exists(f"{trackfile}.gz") and not os.path.exists(f"{trackfile}"):
            os.system(f"gzip -d {trackfile}.gz")
        
        if fdataset=='ERA5':
            adddatafile=adddatadir+'tcwv_ERA5_'+str(y)+'_'+season+'.nc'
        elif fdataset=='ESA' or fdataset=='ESAcross':
            adddatafile=adddatadir+'tcwv_ESA_'+str(y)+'_'+season+'.nc'
        elif fdataset=='MERRA2':
            adddatafile=adddatadir+'tcwv_MERRA2_'+str(y)+'_'+season+'.nc'
        elif fdataset=='JRA3Q':
            adddatafile=adddatadir+'tcwv_JRA3Q_'+str(y)+'_'+season+'.nc'

        print('Added field to be added: '+adddatafile)

        if not os.path.exists(adddatafile):
            print(f"WARNING: file {adddatafile} does not exist")
            continue

        if run_ffi:
            print('Filtering ff_trs based on intensity')
            track_wrapper.filt_trtrs(trackfile,0,5,"ff5i")
            trackfile=trackfile.replace("ff_trs","ff5i_trs")
            fftype=fftype0.replace("ff_trs","ff5i_trs")
        elif run_ff2d:
            print('Filtering tr_trs tracks to two days (no displacement)')
            track_wrapper.filt_trtrs(trackfile,8,10,"ff2d10i")
            trackfile=trackfile.replace("tr_trs","ff2d10i_trs")
            fftype=fftype0.replace("tr_trs","ff2d10i_trs")
        else:
            fftype=fftype0
            
        # temporal subsampling
        if run_time_subsampling:
            if res=='' or res=='50km-day' or res=='50km-day-masked' or fdataset=='ESA' or fdataset=='ESAcross':
                print('<-- Added field has daily time resolution. Subsampling tracks to 12UTC daily values')
                # if clean
                # os.system('rm '+trackdiry+'/subsampled_day/ff_trs_neg')

                if not os.path.exists(trackdiry+'/subsampled_day12h/'+fftype):
                    track_wrapper.subsample_timesteps(trackfile,4,2)
            
        if  res=='50km-day' or res=='50km-day-masked' or fdataset=='ESA' or fdataset=='ESAcross':
            trackfile=trackdiry+'/subsampled_day12h/'+fftype

        if run_addtcwv:
            print('Adding file to tracks: started')               
            newtrack=trackdiry+'/'+fftype+'_'+fdata
            os.system('cp '+trackfile+' '+newtrack)
            track_wrapper.add_mean_field(adddatafile,newtrack,rdd,"TCWV",scaling=1,missing=fmiss,interpolation='nearest')
            track_wrapper.add_mean_field(adddatafile,f"{newtrack}.TCWV{rd}mean",rdd,"TCWV",scaling=1.,operation='min',missing=fmiss,interpolation='nearest')
            track_wrapper.add_mean_field(adddatafile,f"{newtrack}.TCWV{rd}mean.TCWV{rd}min",rdd,"TCWV",scaling=1.,operation='max',missing=fmiss,interpolation='nearest')
            
            # clean
            #os.system('rm '+newtrack + ' ' + newtrack+'.TCWV5mean' + ' ' + newtrack+'.TCWV5mean.TCWV5min')
            os.system(f"rm {newtrack} {newtrack}.TCWV{rd}mean {newtrack}.TCWV{rd}mean.TCWV{rd}min")

            print('Adding file to tracks: finished')

if run_stats:
    for s in range(len(aseas)):
        season=aseas[s]
        sy=seasdict[season][0]
        ly=seasdict[season][1] 

        print('Processing STATS for season '+ season)
        trackdirin=trackdir+'/'+season+'/msl/'

        fname=f"{fftype}_{fdata}.TCWV{rd}mean.TCWV{rd}min.TCWV{rd}max"
        track_wrapper.stats(trackdirin,fftype,'std',sy,ly) # not working, right?
        track_wrapper.stats(trackdirin,fname,'add1',sy,ly)
        track_wrapper.stats(trackdirin,fname,'add2',sy,ly)
        track_wrapper.stats(trackdirin,fname,'add3',sy,ly)

        # rename and combine for simplicity
        fstats0=f'{trackdirin}/stats/{fftype}_scl.std_{sy}-{ly}_1.nc'
        fstats1=f'{trackdirin}/stats/{fname}_scl.add1_{sy}-{ly}_1.nc'
        os.system(f'ncrename -v mstr,meantcwv {fstats1}')
        fstats2=f'{trackdirin}/stats/{fname}_scl.add2_{sy}-{ly}_1.nc'
        os.system(f'ncrename -v mstr,mintcwv {fstats2}')
        fstats3=f'{trackdirin}/stats/{fname}_scl.add3_{sy}-{ly}_1.nc'
        os.system(f'ncrename -v mstr,maxtcwv {fstats3}')
        fstatscomb=f'{trackdirin}/stats/{fname}_scl.stats_{sy}-{ly}_1.nc'
        if os.path.exists(fstatscomb):
            os.system(f'rm {fstatscomb}')
        os.system(f'cdo merge {fstats0} -selvar,meantcwv {fstats1} -selvar,mintcwv {fstats2} -selvar,maxtcwv {fstats3} {fstatscomb}')

        os.system('mkdir -p '+trackdirin+'/stats/unscaled')
        os.system(f'mv {trackdirin}/stats/{fname}.add*{sy}-{ly}_1.nc {trackdirin}/stats/unscaled/.')
        os.system(f'mv {trackdirin}/stats/{fftype}.std*{sy}-{ly}_1.nc {trackdirin}/stats/unscaled/.')

        os.system('mkdir -p '+trackdirin+'/stats/combined')
        os.system(f'mv {trackdirin}/stats/combined*{fname}*{sy}-{ly}* {trackdirin}/stats/combined/.')
        os.system(f'mv {trackdirin}/stats/combined*{fftype}*{sy}-{ly}* {trackdirin}/stats/combined/.')

        os.system('mkdir -p '+trackdirin+'/stats/indat')
        os.system(f'mv {trackdirin}/stats/STATS_*{fname}*{sy}-{ly}*.in {trackdirin}/stats/indat/.')
        os.system(f'mv {trackdirin}/stats/STATS_standard_{fftype}*{sy}-{ly}*.in {trackdirin}/stats/indat/.')
        #os.system(f'mv {trackdirin}/stats/combine*{fname}*{sy}-{ly}*.in {trackdirin}/stats/indat/.')        
        os.system(f'mv {trackdirin}/stats/combine_{sy}-{ly}*.in {trackdirin}/stats/indat/.')        

        # clean
        os.system(f'rm {fstats0} {fstats1} {fstats2} {fstats3}')

f_timestorm='12utc' # 12utc/any                                                                                                       
if run_radial:
    for s in range(len(aseas)):
        season=aseas[s]
        sy=seasdict[season][0]
        ly=seasdict[season][1] 

        print('Processing RADIAL for season '+ season)
        trackdirin=trackdir+'/'+season+'/msl/'
        if fdataset in ['ERA5','MERRA2','JRA3Q'] and res in ['','50km']:
            if f_timestorm == '12utc':
                trackfilein='subsampled_day12h/'+fftype
                if res == '':
                    adddatadir1=adddatadir.replace('6h_v2/','12UTC_v2/')
                elif res == '50km': # still using 6hourly track information here
                    trackfilein=fftype
                    adddatadir1=f"{adddatadir}"
                    #adddatadir1=adddatadir.replace('6h_50km_v2/','12UTC_50km_v2/')
                    #print('datainput 12UTC_50km_v2 NOT READY: CONTINUE')
            elif f_timestorm == 'any':
                trackfilein=fftype
                adddatadir1=f"{adddatadir}"
        #elif fdata in ['ERA5-50km-day','ERA5-50km-day-masked','ESA','ESA-orig']:
        else:
            trackfilein='subsampled_day12h/'+fftype
            adddatadir1=f"{adddatadir}"

        if fdataset=='ESA' or fdataset=='ESAcross':
            filenc=adddatadir1+'tcwv_ESA'+mask1+'_YYYY_'+season+'_filled.nc'
        elif res in ['','50km']:
            filenc=adddatadir1+'tcwv_'+fdataset+'_YYYY_'+season+'.nc'
        elif res in ['50km-day','50km-day-masked']:
            filenc=adddatadir1+'tcwv_'+fdataset+''+mask1+'_YYYY_'+season+'_filled.nc'

        if mask1 != '': 
            for i in range(sy,ly+1):
                filenc1=filenc.replace('YYYY',str(i))
                if not os.path.exists(filenc1):
                    filenc1_unfilled=filenc1.replace('_filled','')
                    os.system(f"cdo setmisstoc,1000000000000000000 {filenc1_unfilled} {filenc1}")

        fdatamask=f"{fdata}{mask1}"

        if args.feature == 'cyc':
            extin=fdatamask
        elif args.feature == 'anticyc':
            extin=f"anticyclones_{fdatamask}"
        elif args.feature == 'cyc2d':
            extin=f"2d_{fdatamask}"
        elif args.feature == 'anticyc2d':
            extin=f"2d_anticyclones_{fdatamask}"

        track_wrapper.radial_maps_2(trackdirin,trackfilein,filenc,'TCWV',intensity=(5,100),rotate=0,sy=sy,ly=ly,interpolation='nearest',ext=extin)
        track_wrapper.radial_maps_2(trackdirin,trackfilein,filenc,'TCWV',intensity=(10,100),rotate=0,sy=sy,ly=ly,interpolation='nearest',ext=extin)
        #track_wrapper.radial_maps_2(trackdirin,trackfilein,filenc,'TCWV',intensity=(10,100),rotate=0,sy=sy,ly=ly,interpolation='nearest',ext=extin,selection='singlestorms')
        

# clean and move anticyclone data for simplicity
if args.feature == 'anticyc':
    for s in range(len(aseas)):
        season=aseas[s]

        trackdirin=trackdir+'/'+season+'/msl/'
        os.system('mkdir -p '+trackdirin+'/radial_maps_anticyclones')
        os.system('mkdir -p '+trackdirin+'/radial_maps_anticyclones/ff_trs')
        os.system('mkdir -p '+trackdirin+'/radial_maps_anticyclones/ff_trs_reg')

        os.system(f'mv {trackdirin}/radial_maps/*anticyclones* {trackdirin}/radial_maps_anticyclones/.')
        os.system(f'mv {trackdirin}/radial_maps/ff_trs/*anticyclones* {trackdirin}/radial_maps_anticyclones/ff_trs/.')
        os.system(f'mv {trackdirin}/radial_maps/ff_trs_reg/*anticyclones* {trackdirin}/radial_maps_anticyclones/ff_trs_reg/.')

        os.system('mkdir -p '+trackdirin+'/stats_anticyclones')
        os.system('mkdir -p '+trackdirin+'/stats_anticyclones/combined')
        os.system('mkdir -p '+trackdirin+'/stats_anticyclones/indat')
        os.system('mkdir -p '+trackdirin+'/stats_anticyclones/scaled')

        os.system(f'mv {trackdirin}/stats/*_pos_* {trackdirin}/stats_anticyclones/.')
        os.system(f'mv {trackdirin}/stats/combined/*_pos_* {trackdirin}/stats_anticyclones/combined/.')
        os.system(f'mv {trackdirin}/stats/scaled/*_pos_* {trackdirin}/stats_anticyclones/scaled/.')
        os.system(f'mv {trackdirin}/stats/indat/*_pos_* {trackdirin}/stats_anticyclones/indat/.')


