#################################################################################################################################
###  CREATING OF 12-DAY STATCUBES of S1 DATA FOR S14AMAZONAS PROJECT
### Script written by Dr. Neha Hunka, dated 10.02.2021
# The statistics are generated using the available mosaic files.
### Script execute with command -
# python3 SARbackscatter_statcubes.py --in_tile "21LYG" --mosaic_folder ".../output/Mosaic"
# --statcubes_folder ".../output/StatCubes"
# --in_multipolygon_string "MultiPolygon(((-55.1650695389999 -11.753601353,-54.158203488 -11.7453511209999,-54.147570924 -12.737064393,-55.15819755 -12.746032956,-55.1650695389999 -11.753601353)))"
# #################################################################################################################################
import logging
import os
import glob
from csv import reader
import json
from time import gmtime
import time
import subprocess
import argparse

from pathlib import Path
from tempfile import mkdtemp

try:
   import gdal
except:
   from osgeo import gdal
import numpy as np
import scipy
from scipy import stats
import pandas as pd
from osgeo import ogr, osr

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

HLINE = '-----------------------------------------------------------------------------------------------'
log = logging.getLogger(__name__)

def init_logging():
    logging.Formatter.converter = gmtime
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter(fmt="{levelname} {message}", style="{"))
    root_log = logging.getLogger()
    root_log.addHandler(console_handler)
    root_log.setLevel(logging.NOTSET)
    log.info("Logging has been started.")

def raster2array(rasterfn):
    raster = gdal.Open(rasterfn)
    band = raster.GetRasterBand(1).ReadAsArray().astype('float')
    return band

def getNoDataValue(rasterfn):
    raster = gdal.Open(rasterfn)
    band = raster.GetRasterBand(1)
    return band.GetNoDataValue()

def getDate(rasterfn):
    acqdate = os.path.basename(rasterfn).split('_')[3].split('.')[0]
    return acqdate

def getDateList(files):
    Date_list = list()
    for i in range(0, len(files)):
        Date = getDate(files[i])
        Date_list.append(Date)
    return Date_list

def getRasterDimensions(rasterfn):
    raster = gdal.Open(str(rasterfn))
    cols = raster.RasterXSize
    rows = raster.RasterYSize
    return [rows,cols]

def array2raster(rasterfn,ORX,ORY,newRasterfn,array):
    raster = gdal.Open(str(rasterfn))
    geotransform = raster.GetGeoTransform() 
    originX = ORX
    originY = ORY
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    cols = raster.RasterXSize
    rows = raster.RasterYSize

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(str(newRasterfn), cols, rows, 1, gdal.GDT_Float32,['COMPRESS=DEFLATE', 'TILED=YES','BIGTIFF=IF_NEEDED'])
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(raster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

Stack_size = 5
VH_pattern = '*BAC*VH*.tif'
POLS = ['VH','VV']

#################################### THE START ##################################################
def run(tile, mosaic_folder, statcubes_folder, multipolygon_string):

    init_logging()
    log.debug("Tile :{}".format(tile))

    # folder inside the backscatter mosaic where the tile mosaic are stored
    mosaic_folder_tile = Path(mosaic_folder).joinpath(tile)
    # the orbit directions for the tile
    ORBIT_DIRECTIONS = [x.stem for x in mosaic_folder_tile.iterdir() if x.is_dir()]
    log.debug("orbit directions {}".format(ORBIT_DIRECTIONS))

    # folder where the backscatter statistics will be stored
    statcubes_folder = Path(statcubes_folder)
    log.debug(HLINE)

    # find the origin of the output tif files
    aoi_geom = ogr.CreateGeometryFromWkt(multipolygon_string)
    (ORX, _, _, ORY) = aoi_geom.GetEnvelope()

    # create the tile folder within the statcubes folder
    Out_folder_tile_stat = statcubes_folder.joinpath(tile)
    if not Out_folder_tile_stat.exists():
        os.makedirs(Out_folder_tile_stat)

    SRC_folder_start = mosaic_folder
    log.debug("output folder in :{}".format(Out_folder_tile_stat))

    for ORBIT_DIRECTION in ORBIT_DIRECTIONS:
        # find statcubes for each orbit direction
        log.debug("Orbit Direction : {}".format(ORBIT_DIRECTION))
        if os.path.exists(Path(SRC_folder_start).joinpath(tile, ORBIT_DIRECTION)):

            # create the orbit folder within the tile-statcube folder
            Out_folder = Path(Out_folder_tile_stat).joinpath(ORBIT_DIRECTION)
            Container_Out_folder = Out_folder
            if not os.path.exists(Out_folder):
                os.makedirs(Out_folder)

            # find mosaic files with VV and VH polarisation
            MOS_folder_tile_orbit_VH = Path(SRC_folder_start).joinpath(tile, ORBIT_DIRECTION, VH_pattern)
            MOS_folder_tile_orbit_VV = Path(SRC_folder_start).joinpath(tile, ORBIT_DIRECTION, VH_pattern.replace('VH', 'VV'))
            ####### GENERATE A LIST OF FILES TO BE ANALYZED
            VH_Files = glob.glob(str(MOS_folder_tile_orbit_VH))
            VH_Files_datelist = [int(x) for x in getDateList(VH_Files)]
            globals()['VH_Files'] = [x for _,x in sorted(zip(VH_Files_datelist,VH_Files))]
            VV_Files = glob.glob(str(MOS_folder_tile_orbit_VV))
            VV_Files_datelist = [int(x) for x in getDateList(VV_Files)]
            globals()['VV_Files'] = [x for _,x in sorted(zip(VV_Files_datelist,VV_Files))]
            log.debug("vh files {}".format(VH_Files_datelist))
            log.debug("vv files {}".format(VV_Files_datelist))

            if len(VH_Files_datelist) < Stack_size or len(VV_Files_datelist) < Stack_size:
               log.debug("vv and vh files less than stack size, skipping the orbit")
               continue

            # find dates which has both VV and VH polarisation mosaics
            Common_dates = list(set(VH_Files_datelist) & set(VV_Files_datelist))
            # create list of mosaic files with these dates
            for POL in POLS:
                pol_r_files_list = [str(Path(SRC_folder_start).joinpath(tile, ORBIT_DIRECTION,
                                                                        tile + '_BAC_' + POL + '_' + str(
                                                                        comm_value) + '.tif')) for comm_value in
                                    Common_dates]
                globals()[(POL + '_r_files')] = pol_r_files_list
            VH_r_files = [x for _,x in sorted(zip(Common_dates,globals()['VH_r_files']))]
            VV_r_files = [x for _,x in sorted(zip(Common_dates,globals()['VV_r_files']))]
            S1_images_list = np.hstack((VH_Files, VV_Files, VH_r_files, VV_r_files))

            MASTER_GRID_TIFF = S1_images_list[0]
            log.debug("MASTER_GRID_TIFF ... {}".format(MASTER_GRID_TIFF))

            ######## GENERATE NEAR-REAL TIME LAST DATES BACKSCATTER STATS ##############
            for POL in POLS:
                for i in range(len((globals()[(POL + '_Files')])) - Stack_size + 1, len((globals()[(POL + '_Files')]))):
                    Date = getDate((globals()[(POL + '_Files')])[i])
                    log.debug( "Currently processing backscatter stats for date: {} and for POL: {}".format(Date, POL))
                    if not os.path.exists(Path(Out_folder).joinpath(POL + '_pmin_' + Date + '.tif')):
                      ### GENERATE PAST STACK
                      rasterarray = []
                      Stack_p_list = list()
                      for p in range(i-Stack_size,i):
                          Stack_p_list.append('raster2array(' + POL + '_Files[' + str(p) + '])')
                          rasterarray.append(raster2array(globals()[(POL + '_Files')][p]))
                      Stack_p    = np.stack(rasterarray, axis=0)

                      Stack_p[Stack_p <= -99] = 'nan'
                      Stack_p_MIN = np.nanmean(Stack_p, axis=0)
                      array2raster(MASTER_GRID_TIFF, ORX, ORY, Path(Out_folder).joinpath(POL + '_pmin_' + Date + '.tif'),Stack_p_MIN)  # (globals()[(POL + '_Files')])[i]
                      log.debug("writing file: {}".format(str(POL + '_pmin_' + Date + '.tif')))

                      ####### GENERATE SINGLE-DATE FUTURE STACK
                      BAC_state = raster2array((globals()[(POL + '_Files')])[i])  # UPDATE IF MULTIPLE ORBITS ARE PROCESSED SEPARATELY
                      BAC_mean_change = np.subtract(BAC_state, Stack_p_MIN)
                      BAC_mean_change[BAC_mean_change <= -99] = 'nan'
                      array2raster(MASTER_GRID_TIFF, ORX, ORY, Path(Out_folder).joinpath(POL + '_mean_change_' + Date + '.tif'),BAC_mean_change)
                      log.debug("writing file: {}".format(str(POL + '_mean_change_' + Date + '.tif')))
                      del Stack_p_MIN, BAC_state, BAC_mean_change

            ######## GENERATE BACKSCATTER STATS ##############
            for POL in POLS:
                for i in range(Stack_size, len((globals()[(POL + '_Files')])) - Stack_size + 1):
                    Date = getDate((globals()[(POL + '_Files')])[i])
                    log.debug("CURRENTLY PROCESSING BACKSCATTER STATS for DATE: {} and for POL: {} -  ".format(Date, POL))
                    if not os.path.exists(Path(Out_folder).joinpath(POL + '_pmin_' + Date + '.tif')):
                      ### GENERATE PAST STACK
                      rasterarray = []
                      Stack_p_list = list()
                      for p in range(i-Stack_size,i):
                          Stack_p_list.append('raster2array(' + POL + '_Files[' + str(p) + '])')
                          rasterarray.append(raster2array(globals()[(POL + '_Files')][p]))
                      Stack_p    = np.stack(rasterarray, axis=0)

                      Stack_p[Stack_p <= -99] = 'nan'
                      Stack_p_MIN = np.nanmean(Stack_p, axis=0)
                      array2raster(MASTER_GRID_TIFF, ORX, ORY, Path(Out_folder).joinpath(POL + '_pmin_' + Date + '.tif'), Stack_p_MIN)  # (globals()[(POL + '_Files')])[i]

                      log.debug("writing file: {}".format(str(POL + '_pmin_' + Date + '.tif')))

                      ### GENERATE (INCLUSIVE) FUTURE STACK
                      rasterarray = []
                      Stack_f_list = list()
                      for f in range(i, i+Stack_size):
                          Stack_f_list.append('raster2array(' + POL + '_Files[' + str(f) + '])')
                          rasterarray.append(raster2array(globals()[(POL + '_Files')][f]))
                      Stack_f    = np.stack(rasterarray, axis=0)

                      Stack_f[Stack_f <= -99] = 'nan'
                      Stack_f_MIN = np.nanmean(Stack_f, axis=0)
                      array2raster(MASTER_GRID_TIFF, ORX, ORY, Path(Out_folder).joinpath(POL + '_fmin_' + Date + '.tif'),Stack_f_MIN)  # (globals()[(POL + '_Files')])[i]

                      log.debug("writing file: {}".format(str(POL + '_fmin_' + Date + '.tif')))

                      ### GENERATE STATISTICS OF PAST AND FUTURE STACK
                      POL_std = np.nanstd(np.vstack((Stack_p,Stack_f)),axis=0)
                      array2raster(MASTER_GRID_TIFF, ORX, ORY,Path(Out_folder).joinpath(POL + '_std_' + Date + '.tif'),POL_std)      #(globals()[(POL + '_Files')])[i]

                      log.debug("writing file: {}".format(str(POL + '_std_' + Date + '.tif')))
                      del POL_std

                      ### MOVING WINDOW TTEST ON STACK
                      POL_mean_change = np.subtract(np.nanmean(Stack_f, axis=0),np.nanmean(Stack_p, axis=0))
                      ttest_pvalue = (scipy.stats.ttest_ind(Stack_p, Stack_f, axis=0,nan_policy='omit')).pvalue
                      ttest_tstatistic = (scipy.stats.ttest_ind(Stack_p, Stack_f, axis=0,nan_policy='omit')).statistic
                      array2raster(MASTER_GRID_TIFF, ORX, ORY, Path(Out_folder).joinpath(POL + '_mean_change_' + Date + '.tif'), POL_mean_change)
                      array2raster(MASTER_GRID_TIFF, ORX, ORY, Path(Out_folder).joinpath(POL + '_ttest_pvalue_' + Date + '.tif'), ttest_pvalue)
                      array2raster(MASTER_GRID_TIFF, ORX, ORY, Path(Out_folder).joinpath(POL + '_ttest_tstatistic_' + Date + '.tif'), ttest_tstatistic)      #(globals()[(POL + '_Files')])[i]

                      log.debug("writing file: {}".format(str(POL + '_mean_change_' + Date + '.tif')))
                      log.debug("writing file: {}".format(str(POL + '_ttest_pvalue_' + Date + '.tif')))
                      log.debug("writing file: {}".format(str(POL + '_ttest_tstatistic_' + Date + '.tif')))
                      del Stack_p,Stack_f,ttest_pvalue,POL_mean_change,ttest_tstatistic,Stack_p_MIN,Stack_f_MIN

                log.debug(HLINE)

            ######## GENERATE BACKSCATTER RATIO STATS ################]
            for i in range(Stack_size, len(VH_r_files) - Stack_size + 1):
                Date = getDate(VH_r_files[i])
                log.debug("CURRENTLY PROCESSING BACKSCATTER RATIO STATS for DATE ... - {} ".format(Date))
                ### GENERATE VV-VH RATIO PAST PLUS FUTURE STACK - i.e. 2 x Stack Size - FOR LINEAR MODEL FIT
                if not os.path.exists(Path(Container_Out_folder).joinpath('VVVHRatio_ttest_tstatistic_' + Date + '.tif')):
                    Stack_r = np.empty([Stack_size*2, getRasterDimensions(VH_r_files[i])[0], getRasterDimensions(VH_r_files[i])[1]])
                    count = -1
                    Date_list = np.empty(Stack_size*2)
                    for r in range(i-Stack_size,i+Stack_size):
                        count +=1
                        Date_list[count] = getDate(VH_r_files[r])
                        VV_for_ratio = raster2array(VV_r_files[r])
                        VV_for_ratio[VV_for_ratio <= -99] = 'nan'
                        VH_for_ratio = raster2array(VH_r_files[r])
                        VH_for_ratio[VH_for_ratio <= -99] = 'nan'
                        Stack_r[count,:,:] = np.subtract(VV_for_ratio,VH_for_ratio)
                        del VV_for_ratio,VH_for_ratio
                    Date_list = [pd.to_datetime(str(int(x)))for x in Date_list]
                    Date_deltas = [dates - Date_list[0] for dates in Date_list]

                    ######## APPLY LINEAR MODEL FIT TO STACK X 2 SIZE VV-VH STACK ########
                    if not os.path.exists(Path(Container_Out_folder).joinpath('VVVHRatio_r_squared_' + Date + '.tif')):
                        x = np.array([int(str(n).split(' ')[0]) for n in Date_deltas])
                        A = np.c_[x, np.ones_like(x)]
                        y = np.reshape(Stack_r,(Stack_size*2, getRasterDimensions(VH_r_files[i])[0]*getRasterDimensions(VH_r_files[i])[1]))
                        col_mean = np.nanmean(y, axis=0)
                        inds = np.where(np.isnan(y))
                        y[inds] = np.take(col_mean,inds[1])
                        m,resid,rank,s = np.linalg.lstsq(A, y,rcond=None)
                        slope_r = np.reshape(m[0],(getRasterDimensions(VH_r_files[i])[0],getRasterDimensions(VH_r_files[i])[1]))
                        array2raster(MASTER_GRID_TIFF, ORX, ORY, Path(Out_folder).joinpath('VVVHRatio_slope_' + Date + '.tif'), slope_r)
                        r_square = 1-resid/(x.size*np.var(y,axis=0))
                        r_square = np.reshape(r_square, (getRasterDimensions(VH_r_files[i])[0], getRasterDimensions(VH_r_files[i])[1]))
                        array2raster(MASTER_GRID_TIFF, ORX, ORY, Path(Out_folder).joinpath('VVVHRatio_r_squared_' + Date + '.tif'), r_square)

                        log.debug("writing file: {}".format(str('VVVHRatio_r_squared_' + Date + '.tif')))
                        del A,Stack_r,y,col_mean,slope_r,resid,m,inds

                    ######## GENERATE VV-VH RATIO PAST STACK
                    if not os.path.exists(Path(Container_Out_folder).joinpath('VVVHRatio_ttest_tstatistic_' + Date + '.tif')):
                        Stack_r_p = np.empty([Stack_size, getRasterDimensions(VH_r_files[i])[0], getRasterDimensions(VH_r_files[i])[1]])
                        count = -1
                        Date_list = np.empty(Stack_size)
                        for r in range(i-Stack_size,i):
                            count +=1
                            Date_list[count] = getDate(VH_r_files[r])
                            VV_for_ratio = raster2array(VV_r_files[r])
                            VV_for_ratio[VV_for_ratio <= -99] = 'nan'
                            VH_for_ratio = raster2array(VH_r_files[r])
                            VH_for_ratio[VH_for_ratio <= -99] = 'nan'
                            Stack_r_p[count,:,:] = np.subtract(VV_for_ratio,VH_for_ratio)
                            del VV_for_ratio
                            del VH_for_ratio

                        ######## GENERATE VV-VH RATIO FUTURE STACK
                        Stack_r_f = np.empty([Stack_size, getRasterDimensions(VH_r_files[i])[0], getRasterDimensions(VH_r_files[i])[1]])
                        count = -1
                        Date_list = np.empty(Stack_size)
                        for r in range(i,i+Stack_size):
                            count +=1
                            Date_list[count] = getDate(VH_r_files[r])
                            VV_for_ratio = raster2array(VV_r_files[r])
                            VV_for_ratio[VV_for_ratio <= -99] = 'nan'
                            VH_for_ratio = raster2array(VH_r_files[r])
                            VH_for_ratio[VH_for_ratio <= -99] = 'nan'
                            Stack_r_f[count,:,:] = np.subtract(VV_for_ratio,VH_for_ratio)
                            del VV_for_ratio
                            del VH_for_ratio

                        ######## MOVING WINDOW TTEST ON VV-VH RATIO STACK
                        Ratio_mean_change = np.subtract(np.nanmean(Stack_r_f, axis=0), np.nanmean(Stack_r_p, axis=0))
                        ttest_pvalue = (scipy.stats.ttest_ind(Stack_r_p, Stack_r_f, axis=0,nan_policy='omit')).pvalue
                        ttest_tstatistic = (scipy.stats.ttest_ind(Stack_r_p, Stack_r_f, axis=0,nan_policy='omit')).statistic
                        array2raster(MASTER_GRID_TIFF, ORX, ORY, Path(Out_folder).joinpath('VVVHRatio_mean_change_' + Date + '.tif'), Ratio_mean_change)
                        array2raster(MASTER_GRID_TIFF, ORX, ORY, Path(Out_folder).joinpath('VVVHRatio_ttest_pvalue_' + Date + '.tif'),ttest_pvalue)
                        array2raster(MASTER_GRID_TIFF, ORX, ORY, Path(Out_folder).joinpath('VVVHRatio_ttest_tstatistic_' + Date + '.tif'),ttest_tstatistic)

                        log.debug("writing file: {}".format(str('VVVHRatio_mean_change_' + Date + '.tif')))
                        del Stack_r_p, Stack_r_f, Ratio_mean_change, ttest_pvalue,ttest_tstatistic
        log.debug(HLINE)
        log.debug(HLINE)

    log.debug("SAR statistics script has finished")

if __name__ == "__main__":

        parser = argparse.ArgumentParser(description="Deforestation detection tool.")
        parser.add_argument("--in_tile", type=str, help="Tile name.", required=True, default=None)
        parser.add_argument("--mosaic_folder", type=str, help="Mosaic folder inside output folder.", required=True, default=None)
        parser.add_argument("--statcubes_folder", type=str, help="StatCubes folder inside output folder.", required=True, default=None)
        parser.add_argument("--in_multipolygon_string", type=str, help="Area of interest in multipolygon/polygon string.",
                            required=True, default=None)
        args = parser.parse_args()

        run(tile=args.in_tile,
            mosaic_folder=args.mosaic_folder,
            statcubes_folder=args.statcubes_folder,
            multipolygon_string=args.in_multipolygon_string
            )
