#################################################################################################################################
### Deforestation detection S1 DATA FOR S14AMAZONAS PROJECT
### Script written by Dr. Neha Hunka, dated 10.02.2021
# The available statcubes are used to detect doforestation events. When the statistics are above a certain
# predefined threshold, it is recorded as a deforestation event.
### Script execute with command -
# python3 SARbackscatter_detection.py --in_tile "21LYG" --in_timeperiod "20150601-20171231"
# --output_dir".../output"
# --aux_dir ".../aux_data"
# #################################################################################################################################

import os
import sys
import logging
import time
import datetime
from datetime import datetime
from datetime import  timedelta
from time import gmtime
from shutil import copyfile
import shutil
from pathlib import Path
import subprocess
from tempfile import mkdtemp
from csv import reader
import json
import argparse
import glob

from scipy import ndimage as nd

try:
    import gdal
except:
    from osgeo import gdal

import numpy as np
from osgeo import ogr, osr

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
logger = logging.getLogger("deforestation")

HLINE = '-----------------------------------------------------------------------------------------------'
DETECTION_SET = 'Detections'

GENERATE_MASK = False

def setup_logging():
    logger.setLevel(logging.DEBUG)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    console_formatter = logging.Formatter("%(levelname)s, %(message)s", "%Y-%m-%d %H:%M:%S")
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

def raster2array(rasterfn):
    raster = gdal.Open(rasterfn)
    band = raster.GetRasterBand(1).ReadAsArray().astype('float')
    return band

def getNoDataValue(rasterfn):
    raster = gdal.Open(rasterfn)
    band = raster.GetRasterBand(1)
    return band.GetNoDataValue()

def getDate(rasterfn):
    acqdate = os.path.basename(rasterfn).split('_')[-1].split('.')[0]
    return acqdate

def getDateList(files):
    Date_list = list()
    for i in range(0, len(files)):
        Date = getDate(files[i])
        Date_list.append(Date)
    return Date_list

def getRasterDimensions(rasterfn):
    raster = gdal.Open(rasterfn)
    band = raster.GetRasterBand(1)
    cols = raster.RasterXSize
    rows = raster.RasterYSize
    return [rows,cols]

def getRasterExtent(rasterfn):
    raster = gdal.Open(rasterfn)
    ulx, xres, xskew, uly, yskew, yres = raster.GetGeoTransform()
    lrx = ulx + (raster.RasterXSize * xres)
    lry = uly + (raster.RasterYSize * yres)
    return [ulx,uly,lrx,lry]

def array2raster(rasterfn,newRasterfn,array):
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    cols = raster.RasterXSize
    rows = raster.RasterYSize

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float32,['COMPRESS=DEFLATE', 'TILED=YES','BIGTIFF=IF_NEEDED'])
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(raster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

def dist(data, invalid=None):
    if invalid is None: invalid = np.isnan(data)
    dis = nd.distance_transform_edt(invalid,return_distances=True,return_indices=False)
    return dis

def fill(data, invalid=None):
    if invalid is None: invalid = np.isnan(data)
    ind = nd.distance_transform_edt(invalid,return_distances=False,return_indices=True)
    return data[tuple(ind)]

def open_json(path):
    with open(path, 'r') as file:
        return json.load(file)

def run(tile, timeperiod, statcubes_dir, output_dir, aux_dir):

    config = open_json('threshold.json')

    # the folder where the statcube of tile is stored
    output_dir = Path(output_dir)
    # the folder where the SRTM, forest cover, paths etc are stored. Refer the next section for all the auxiliary data
    aux_dir = Path(aux_dir)
    # path to the statcubes dir
    statcubes_dir = Path(statcubes_dir)

    start_time = time.time()
    ####################################################################################################
    setup_logging()

    logger.debug("Deforestation main run has started.")

    GFW_treecover = aux_dir.joinpath('CCI', 'ESACCI-LC-L4-LCCS-Map-300m-2015_Amazon_only30to90_20m.tif')
    if GFW_treecover.exists():
       logger.debug("GFW treecover exists")
    else:
       logger.debug("GFW treecover doesnt exists")

    SRTM_path = aux_dir.joinpath('SRTM')
    if SRTM_path.exists():
       logger.debug("SRTM path exists")
    else:
       logger.debug("SRTM path doesnt exists")

    HYDRO_layer = aux_dir.joinpath('hydrography_biome', 'hydrography_biome_holes1_dissolved_003deg_buffer.shp')
    if HYDRO_layer.exists():
       logger.debug("Hydro layer exists")
    else:
       logger.debug("Hydro layer doesnt exists")

    THRESHOLD_LIST = aux_dir.joinpath('geom_data', 'S2_grid_AmazonBasin_detections_thresholds.csv')
    if THRESHOLD_LIST.exists():
       logger.debug("Threshold list exists")
    else:
       logger.debug("Threshold list doesnt exists")

    PATHS_LIST = aux_dir.joinpath('geom_data', 'S2_grid_AmazonBasin_detections_paths.csv')
    if PATHS_LIST.exists():
       logger.debug("Path list exists")
    else:
       logger.debug("Path list doesnt exists")

    # find whether the tile has one or two orbit direction
    mosaic_folder_tile = Path(statcubes_dir).joinpath(tile)
    ORBIT_DIRECTIONS = [x.stem for x in mosaic_folder_tile.iterdir() if x.is_dir()]
    # find the corresponding SRTM for the tile
    with open(PATHS_LIST, 'r') as read_paths:
        txt_reader = reader(read_paths)
        for each_Ama_tile in txt_reader:
            Ama_tile = str(each_Ama_tile[0])
            if Ama_tile == str(tile):
                Ama_path = str(each_Ama_tile[1])
                SRTM_layer = Path(SRTM_path).joinpath(str(each_Ama_tile[3]))
                logger.debug("USE SRTM layer: {}".format(SRTM_layer))

    # find deforestation events for each orbit and then merge them
    Tile_count = -1
    for ORBIT_DIRECTION in ORBIT_DIRECTIONS:
        logger.debug(HLINE)
        logger.debug(HLINE)

        logger.debug("CURRENT ON ORBIT {}".format(ORBIT_DIRECTION))

        # Thresholds
        THR_pvalue = config["THR_pvalue"]["value"]
        THR_tstat = config["THR_tstat"]["value"]

        VV_THR_mean_change = config["VV_THR_mean_change"]["value"]
        VH_THR_mean_change = config["VH_THR_mean_change"]["value"]

        VH_THR_std = config["VH_THR_std"]["value"]
        VV_THR_std = config["VV_THR_std"]["value"]

        sieve_size = config["sieve_size"]["value"]

        VH_THR_fmin = config["VH_THR_fmin"]["value"]
        VV_THR_fmin = config["VV_THR_fmin"]["value"]

        VH_THR_pmin = config["VH_THR_pmin"]["value"]
        VV_THR_pmin = config["VV_THR_pmin"]["value"]

        THR_BASEMaxSTD = config["THR_BASEMaxSTD"]["value"]

        VVVHR_THR_slope = config["VVVHR_THR_slope"]["value"]
        VVVH_THR_R2 = config["VVVH_THR_R2"]["value"]
        VVVHR_THR_mean_change = config["VVVHR_THR_mean_change"]["value"]

        COUNT_THR = config["COUNT_THR"]["value"] # Number of StatCubes that should be valid in 1 day for a detection to be TRUE
        NUM_Changes_THR = config["NUM_Changes_THR"]["value"]  # Min. number of times the change can be seen in whole time series
        DISTANCE_to_WATER_THR = config["DISTANCE_to_WATER_THR"]["value"]  # This is a placeholder, in case polygons very close to water have to be removed

        THR_SRTM_elevation_mountains = config["THR_SRTM_elevation_mountains"]["value"]  # Areas above this height are mountains and are masked out.
        THR_SRTM_elevation = config["THR_SRTM_elevation"]["value"]  # Areas below elevation of 40 m are masked out. These are usually flood plains + swamps
        THR_SRTM_elevation_cutoff = config["THR_SRTM_elevation_cutoff"]["value"]  # Areas below elevation of 25 m are masked out with no exceptions.

        VH_THR_fmin_water = config["VH_THR_fmin_water"]["value"]
        VV_THR_fmin_water = config["VV_THR_fmin_water"]["value"]

        ### CREATE A GRID-LIST OF DATES EVERY 6/12 DAYS (for 12-day acquisition frequency) ###
        acq_frequency = 12

        stack_size = -5 # Be careful to specify the stack size used in generating the statcubes!! The negative number means the last few days are marked.

        OUT_folder   = output_dir.joinpath('StatCubes',tile, ORBIT_DIRECTION)
        OUT_d_folder = output_dir.joinpath(DETECTION_SET, tile, ORBIT_DIRECTION)
        work_dir_detection_folder = output_dir.joinpath(DETECTION_SET)

        statcubes_files = os.listdir(OUT_folder)
        statcubes_datelist = [int(x) for x in getDateList(statcubes_files)]
        startdate_statcubes_str = str(sorted(statcubes_datelist)[0])

        startdatestr = timeperiod.split('-')[0]
        enddatestr = timeperiod.split('-')[1]

        # Adjust the start date so that it matches the available statcube files
        start_date = datetime.strptime(startdate_statcubes_str, '%Y%m%d')
        end_date = datetime.strptime(enddatestr, '%Y%m%d')
        days_interval = np.ndarray.tolist(
            np.arange(start_date, end_date, timedelta(days=acq_frequency), dtype='datetime64[D]'))
        logger.debug("performing detection for date above {}".format(start_date))
        year = start_date.year
        last_days = days_interval[stack_size:]
        change_start = np.ndarray.tolist(np.arange(datetime(2017, 1, 1), datetime(2017, 1, 2), timedelta(days=acq_frequency), dtype='datetime64[D]'))

        logger.debug("creating Statcubes folder and Detections folder")
        if not os.path.exists(OUT_folder):
            os.makedirs(OUT_folder)
        if not os.path.exists(OUT_d_folder):
            os.makedirs(OUT_d_folder)

        POLS = ['VH','VV']

        logger.debug(HLINE)
        logger.debug("-- BASEMin --")
        logger.debug(HLINE)

        VV_BASEMIN_array_tif_path = Path(OUT_d_folder).joinpath('VV_BASEMIN_array.tif')
        VH_BASEMIN_array_tif_path = Path(OUT_d_folder).joinpath('VH_BASEMIN_array.tif')

        ################### GENERATE BASE YEAR FOREST/NON-FOREST MASKS ##################
        basemin_array_pol_orbit_created = 0
        basemin_array_pol_orbit_copied = 0
        rclone_baseminmask_registry = []
        if not VV_BASEMIN_array_tif_path.exists() or not VH_BASEMIN_array_tif_path.exists() or GENERATE_MASK: # or not os.path.exists(OUT_d_folder + 'VH_Greyzone.tif')
           logger.debug("VV VH basemin array doesnt exist. Generating them")
           for root, dirs, files in os.walk(OUT_folder):
               for POL in POLS:
                   count = -1
                   for file in files:
                       if file.__contains__(POL + '_pmin_2015') or file.__contains__(POL + '_pmin_2016') or file.__contains__(POL + '_fmin_2015') or file.__contains__(POL + '_fmin_2016'):
                           count += 1
                           if count == 0:
                               MASTER_GRID_TIFF = os.path.join(OUT_folder,file)
                               logger.debug("Reading ... {}".format(file))
                               MINS = raster2array(os.path.join(root, file))
                           else:
                               logger.debug("Reading ... {}".format(file))
                               MINS_NEXT = raster2array(os.path.join(root,file))
                               MINS = np.nanmin(np.stack((MINS,MINS_NEXT),axis=0),axis=0)

                   if 'MINS' in globals() or 'MINS' in locals():
                       array2raster(MASTER_GRID_TIFF, str(OUT_d_folder.joinpath(POL + '_BASEMIN_array.tif')), MINS)
                       del MINS

        else:
           logger.debug("VV or VH basemin array exists.  Moving on to basemaxstd")

        # at this stage BASEMIN array tif files will be generated. Check if both VV and VH are present and make a copy if missing
        VV_basemin_array_tif_path_local = Path(OUT_d_folder).joinpath('VV_BASEMIN_array.tif')
        VH_basemin_array_tif_path_local = Path(OUT_d_folder).joinpath('VH_BASEMIN_array.tif')
        if not VH_basemin_array_tif_path_local.exists() and VV_basemin_array_tif_path_local.exists():
           logger.debug("VH basemin array is not generated and VV basemin is gen. Copying {} to {}".format(str(VV_basemin_array_tif_path_local), str(VH_basemin_array_tif_path_local)))
           copyfile(str(VV_basemin_array_tif_path_local), str(VH_basemin_array_tif_path_local))

        for POL in POLS:
            if set(ORBIT_DIRECTIONS) == set(['descending','ascending']) and ORBIT_DIRECTION == ORBIT_DIRECTIONS[1]:
               logger.debug("{} checks.".format(ORBIT_DIRECTIONS[1]))
               basemin_array_orbit0 = output_dir.joinpath(DETECTION_SET, tile, ORBIT_DIRECTIONS[0], POL + '_BASEMIN_array.tif')
               basemin_array_orbit1 = output_dir.joinpath(DETECTION_SET, tile, ORBIT_DIRECTIONS[1], POL + '_BASEMIN_array.tif')
               if basemin_array_orbit0.exists() and not basemin_array_orbit1.exists():
                  logger.debug("{} exists and {} doesnt exists. Copy first to second".format(basemin_array_orbit0, basemin_array_orbit1))
                  copyfile(basemin_array_orbit0, basemin_array_orbit1)

        ####################################################### BASEMaxSTD arrays ######################################################
        logger.debug(HLINE)
        logger.debug("-- BASEMaxSTD --")
        logger.debug(HLINE)

        rclone_mask_registry = []
        VV_BASEMaxSTD_array_tif_path = Path(OUT_d_folder).joinpath('VV_BASEMaxSTD_array.tif')
        VH_BASEMaxSTD_array_tif_path = Path(OUT_d_folder).joinpath('VH_BASEMaxSTD_array.tif')
        if not VV_BASEMaxSTD_array_tif_path.exists() or not VH_BASEMaxSTD_array_tif_path.exists() or GENERATE_MASK:
           logger.debug("CREATING BASEMaxSTD_array ....")
           for root, dirs, files in os.walk(OUT_folder):
              for POL in POLS:
                  count = -1
                  for file in files:
                      if file.__contains__(POL + '_std_2015') or file.__contains__(POL + '_std_2016'):
                         count += 1
                         if count == 0:
                            MASTER_GRID_TIFF = os.path.join(root,file)
                            logger.debug("Reading :{} ".format(file))
                            MINS = raster2array(os.path.join(root, file))
                         else:
                            logger.debug("Reading :{} ".format(file))
                            MINS_NEXT = raster2array(os.path.join(root,file))
                            MINS = np.nanmax(np.stack((MINS,MINS_NEXT),axis=0),axis=0)
                  if 'MINS' in globals() or 'MINS' in locals():
                     array2raster(MASTER_GRID_TIFF, str(OUT_d_folder.joinpath(POL + '_BASEMaxSTD_array.tif')), MINS)
                     del MINS
                     logger.debug("WRITING FILE: {}".format(str(OUT_d_folder.joinpath(POL + '_BASEMaxSTD_array.tif'))))

        # Allow exception if there is not enough VH data to make a basemap
        VV_basemaxstd_array_tif_path_local = Path(OUT_d_folder).joinpath('VV_BASEMaxSTD_array.tif')
        VH_basemaxstd_array_tif_path_local = Path(OUT_d_folder).joinpath('VH_BASEMaxSTD_array.tif')
        if not VH_basemaxstd_array_tif_path_local.exists() and VV_basemaxstd_array_tif_path_local.exists():
           logger.debug("VV exist  while VH mask was not created. Copying {} to {}".format(str(VV_basemaxstd_array_tif_path_local), str(VH_basemaxstd_array_tif_path_local)))
           copyfile(str(VV_basemaxstd_array_tif_path_local),str(VH_basemaxstd_array_tif_path_local))

        for POL in POLS:
            if set(ORBIT_DIRECTIONS) == set(['descending','ascending']) and ORBIT_DIRECTION == ORBIT_DIRECTIONS[1]:
               basemaxstd_array_orbit0 = output_dir.joinpath(DETECTION_SET, tile, ORBIT_DIRECTIONS[0], POL + '_BASEMaxSTD_array.tif')
               basemaxstd_array_orbit1 = output_dir.joinpath(DETECTION_SET, tile, ORBIT_DIRECTIONS[1], POL + '_BASEMaxSTD_array.tif')
               if basemaxstd_array_orbit0.exists() and not basemaxstd_array_orbit1.exists():
                  copyfile(basemaxstd_array_orbit0, basemaxstd_array_orbit1)

        ####################################################### GFWTC mask ######################################################
        logger.debug(HLINE)
        logger.debug("-- GFWTC mask --")
        logger.debug(HLINE)

        GFWTC_mask_orbit_tif_path = Path(OUT_d_folder).joinpath(tile + '_GFWTC_mask.tif')
        if not GFWTC_mask_orbit_tif_path.exists() or GENERATE_MASK:
            logger.debug("Generating GFWTC_file_mask for orbit:{}".format(ORBIT_DIRECTION))
            if GFW_treecover.exists():
                MASTER_GRID_TIFF = OUT_d_folder.joinpath('VV_BASEMIN_array.tif')
                rows = getRasterDimensions(str(MASTER_GRID_TIFF))[0]
                cols = getRasterDimensions(str(MASTER_GRID_TIFF))[1]
                extent = getRasterExtent(str(MASTER_GRID_TIFF))

                cmd_GFWTC_tile = ["gdal_translate",
                                  "-ot", "Int32",
                                  "-r", "nearest", "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES",
                                  "-co", "BIGTIFF=IF_NEEDED",
                                  "-projwin",
                                  "{}".format(extent[0]), "{}".format(extent[1]), "{}".format(extent[2]), "{}".format(extent[3]),
                                  "-outsize",
                                  "{}".format(str(cols)), "{}".format(str(rows)),
                                  "{}".format(str(GFW_treecover)), "{}".format(str(OUT_d_folder.joinpath(tile + '_GFWTC_mask.tif')))]
                cmd_output = subprocess.run(cmd_GFWTC_tile, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                logger.debug("exit code {} --> {}".format(cmd_output.returncode, cmd_GFWTC_tile))

                GFWTC_file_mask = raster2array(str(OUT_d_folder.joinpath(tile + '_GFWTC_mask.tif')))
                np.nan_to_num(GFWTC_file_mask)
                array2raster(str(MASTER_GRID_TIFF), str(OUT_d_folder.joinpath(tile + '_GFWTC_mask.tif')), GFWTC_file_mask)

        ################################################### SRTM elevation mask ####################################################
        logger.debug(HLINE)
        logger.debug("-- SRTM elevation mask --")
        logger.debug(HLINE)

        SRTM_Elevation_mask_orbit_tif_path = Path(OUT_d_folder).joinpath(tile + '_SRTM_Elevation_mask.tif')
        if not SRTM_Elevation_mask_orbit_tif_path.exists() or GENERATE_MASK: # and os.path.exists(OUT_d_folder + 'VH_Greyzone.tif')
            logger.debug("Creating ... SRTM_Elevation_mask")
            MASTER_GRID_TIFF = OUT_d_folder.joinpath('VV_BASEMIN_array.tif')
            rows = getRasterDimensions(str(MASTER_GRID_TIFF))[0]
            cols = getRasterDimensions(str(MASTER_GRID_TIFF))[1]
            extent = getRasterExtent(str(MASTER_GRID_TIFF))
            if not OUT_d_folder.joinpath(tile + '_SRTM.tif').exists():
                cmd_SRTM_tile = ["gdal_translate",
                                  "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES",
                                  "-co", "BIGTIFF=IF_NEEDED",
                                  "-r", "average",
                                  "-projwin",
                                  "{}".format(extent[0]), "{}".format(extent[1]), "{}".format(extent[2]), "{}".format(extent[3]),
                                  "-outsize",
                                  "{}".format(str(cols)), "{}".format(str(rows)),
                                  "{}".format(str(SRTM_layer)), "{}".format(str(OUT_d_folder.joinpath(tile + '_SRTM.tif')))]
                cmd_output = subprocess.run(cmd_SRTM_tile, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                logger.debug("exit code {} --> {}".format(cmd_output.returncode, cmd_SRTM_tile))

            SRTM_ELEV = raster2array(str(OUT_d_folder.joinpath(tile + '_SRTM.tif')))
            SRTM_GZ   = raster2array(str(OUT_d_folder.joinpath(tile + '_SRTM.tif')))
            SRTM_GZ[SRTM_GZ > THR_SRTM_elevation_mountains] = 1
            SRTM_GZ[SRTM_GZ <= THR_SRTM_elevation] = 1
            SRTM_GZ[SRTM_GZ > THR_SRTM_elevation] = 0

            if GFWTC_mask_orbit_tif_path.exists():
                GFWTC_file_mask = raster2array(str(GFWTC_mask_orbit_tif_path))
                SRTM_GZ[(GFWTC_file_mask == 1)] = 1  # THIS IS AN EXTERNAL MASK - remove areas that were not forest
            array2raster(str(MASTER_GRID_TIFF), str(OUT_d_folder.joinpath(tile + '_SRTM_Elevation_mask.tif')), SRTM_GZ)

        ################# GENERATE ITERATIVE WATER MASK FOR ALL YEARS INCLUDING BASE YEARS ##################
        logger.debug(HLINE)
        logger.debug("-- Water mask --")
        logger.debug(HLINE)

        watermask_orbit_tif_path = Path(OUT_d_folder).joinpath(tile + '_WATER_mask' + '_s500-8LT.tif')
        watermask_tif = Path(OUT_d_folder).joinpath(tile + '_WATER_mask.tif')
        if not watermask_orbit_tif_path.exists() or GENERATE_MASK:
            logger.debug("Sieving ITERATIVE WATER_mask and adding STRM and PRODES TERRABRAZILIS HYDRO LAYER ....")
            if not watermask_tif.exists() or GENERATE_MASK:
                logger.debug("CREATING ITERATIVE WATER_mask ....")
                for root, dirs, files in os.walk(OUT_folder):
                    for POL in POLS:
                        for file in files:
                            Tile_count += 1
                            if Tile_count == 0:
                                logger.debug("Creating empty WATER MASK with ...{}".format(file))
                                rows = getRasterDimensions(str(OUT_folder.joinpath(file)))[0]
                                cols = getRasterDimensions(str(OUT_folder.joinpath(file)))[1]
                                WATER_mask = np.zeros(shape=(rows, cols))
                                MASTER_GRID_TIFF = str(OUT_folder.joinpath(file))
                            else:
                                WATER_mask = WATER_mask
                            if file.__contains__(POL + '_fmin_'):
                                logger.debug("Reading file...{}".format(file))
                                FMIN = raster2array(str(OUT_folder.joinpath(file)))
                                WATER_mask[FMIN <= locals()[POL + '_THR_fmin_water']] = 1
                array2raster(str(MASTER_GRID_TIFF), str(OUT_d_folder.joinpath(tile + '_WATER_mask.tif')), WATER_mask)
                logger.debug("Writing file: {}".format(str(OUT_d_folder.joinpath(tile + '_WATER_mask.tif'))))

            ###### MERGE OUR LARGE_WATER_BODIES and SRTM MASK LAYER #########
            MASTER_GRID_TIFF = os.path.join(OUT_d_folder.joinpath('VV_BASEMIN_array.tif'))
            logger.debug("MERGE OUR LARGE_WATER_BODIES and SRTM MASK LAYER")
            SRTM_GZ = raster2array(str(OUT_d_folder.joinpath(tile + '_SRTM_Elevation_mask.tif')))
            WATER_mask[SRTM_GZ == 1] = 1
            array2raster(MASTER_GRID_TIFF, str(OUT_d_folder.joinpath(tile + '_WATER_mask.tif')), WATER_mask)

            ######## SIEVE FILTER WATER MASK TO KEEP ONLY LARGE WATER BODIES ############
            logger.debug("APPLYING SIEVE FILTER TO MAKE LARGE WATER BODIES")
            cmd_sieve = ["gdal_sieve.py",
                         "-st", "500", "-8",
                         "-nomask", "-of", "GTiff",
                         "{}".format(str(OUT_d_folder.joinpath(tile + '_WATER_mask.tif'))),
                         "{}".format(str(OUT_d_folder.joinpath(tile + '_WATER_mask' + '_s500-8LT.tif')))]
            cmd_output = subprocess.run(cmd_sieve, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            logger.debug("exit code {} --> {}".format(cmd_output.returncode, cmd_sieve))

            ###### MERGE OUR LARGE_WATER_BODIES and PRODES TERRABRAZILIS HYDRO LAYER #########
            if HYDRO_layer.exists():
                logger.debug("MERGE OUR LARGE_WATER_BODIES and PRODES TERRABRAZILIS HYDRO LAYER")
                OUT_water_file = OUT_d_folder.joinpath(tile + '_WATER_mask' + '_s500-8LT.tif')
                cmd_water = ["gdal_rasterize",
                             "-burn", "1", "-at",
                             "{}".format(str(HYDRO_layer)),
                             "{}".format(str(OUT_water_file))]
                cmd_output = subprocess.run(cmd_water, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                logger.debug("exit code {} --> {}".format(cmd_output.returncode, cmd_water))

        logger.debug(HLINE)
        logger.debug("--- Going through the date for detection ---")
        logger.debug(HLINE)

        datepoint = 0
        for date in days_interval:
            if date >= change_start[0]:
                logger.debug(HLINE)
                count = -1
                datepoint += 1
                date_pattern = [str(date).replace('-', '')]
                logger.debug("--- date pattern : {} ---".format(date_pattern))
                logger.debug("threshold checks")
                ##################### APPLY THRESHOLDS TO STATCUBES RASTERS ##################
                for root, dirs, files in os.walk(OUT_folder):
                    for file in files:
                        count += 1
                        if count == 0:
                            rows = getRasterDimensions(os.path.join(OUT_folder,file))[0]
                            cols = getRasterDimensions(os.path.join(OUT_folder,file))[1]
                            extent = getRasterExtent(os.path.join(OUT_folder, file))
                            DEC_array = np.zeros(shape=(rows, cols))
                            MASTER_GRID_TIFF = os.path.join(OUT_folder,file)
                        if file.__contains__(str(date).replace('-', '')):
                            logger.debug("processing file :{}".format(file))
                            for POL in POLS:
                                if file.__contains__(POL + '_std'):
                                    logger.debug("- in std")
                                    STD = raster2array(os.path.join(OUT_folder,file))
                                    STD[STD < locals()[(POL + '_THR_std')]] = 0
                                    STD[np.isnan(STD)] = 0
                                    STD[STD >= locals()[(POL + '_THR_std')]] = 1
                                    DEC_array = DEC_array + STD
                                if file.__contains__(POL + '_mean_change'):
                                    logger.debug("- in mean_change")
                                    Mean_change = raster2array(os.path.join(OUT_folder, file))
                                    Mean_change[Mean_change > locals()[(POL + '_THR_mean_change')]] = 0
                                    Mean_change[Mean_change <= locals()[(POL + '_THR_mean_change')]] = 1
                                    Mean_change[np.isnan(Mean_change)] = 0
                                    DEC_array = DEC_array + Mean_change
                                    if os.path.exists(os.path.join(OUT_folder, file.replace('mean_change', 'ttest_tstatistic'))):
                                        logger.debug("- in ttest_tstatistic")
                                        TSTAT = raster2array(os.path.join(OUT_folder, file.replace('mean_change', 'ttest_tstatistic')))
                                        TSTAT[TSTAT < THR_tstat] = 0
                                        TSTAT[np.isnan(TSTAT)] = 0
                                        TSTAT[TSTAT >= THR_tstat] = 1
                                        TSTAT[Mean_change == 0] = 0
                                        DEC_array = DEC_array + TSTAT
                                        if os.path.exists(os.path.join(OUT_folder,file.replace('mean_change','ttest_pvalue'))):
                                            logger.debug("- in ttest_pvalue")
                                            PVAL = raster2array(os.path.join(OUT_folder,file.replace('mean_change','ttest_pvalue')))
                                            PVAL[PVAL > THR_pvalue] = 2
                                            PVAL[PVAL <= THR_pvalue] = 1
                                            PVAL[np.isnan(PVAL)] = 2
                                            PVAL[PVAL == 2] = 0
                                            PVAL[TSTAT == 0] = 0
                                            DEC_array = DEC_array + PVAL
                            if file.__contains__('VVVHRatio_slope'):
                                logger.debug("- in vvvh ratio")
                                RSLOPE = raster2array(os.path.join(OUT_folder, file))
                                RSLOPE[RSLOPE >= VVVHR_THR_slope] = 1
                                RSLOPE[RSLOPE < VVVHR_THR_slope] = 0
                                RSLOPE[np.isnan(RSLOPE)] = 0
                                DEC_array = DEC_array + RSLOPE
                            if file.__contains__('VVVHRatio_r_squared'):
                                logger.debug("- in vvvh ratio sq")
                                RR2 = raster2array(os.path.join(OUT_folder, file))
                                RR2[RR2 >= VVVH_THR_R2] = 1
                                RR2[RR2 < VVVH_THR_R2] = 0
                                RR2[np.isnan(RR2)] = 0
                                DEC_array = DEC_array + RR2
                            if file.__contains__('VVVHRatio_mean_change'):
                                logger.debug("- in vvvh mean change")
                                R_Mean_change = raster2array(os.path.join(OUT_folder, file))
                                R_Mean_change[R_Mean_change < VVVHR_THR_mean_change] = 0
                                R_Mean_change[R_Mean_change >= VVVHR_THR_mean_change] = 1
                                R_Mean_change[np.isnan(R_Mean_change)] = 0
                                DEC_array = DEC_array + R_Mean_change
                                if os.path.exists(os.path.join(OUT_folder, file.replace('mean_change', 'ttest_tstatistic'))):
                                    logger.debug("- in vvvh ttest tstatistic")
                                    RTSTAT = raster2array(os.path.join(OUT_folder, file.replace('mean_change', 'ttest_tstatistic')))
                                    RTSTAT[RTSTAT <= THR_tstat] = 1
                                    RTSTAT[RTSTAT > THR_tstat] = 0
                                    RTSTAT[np.isnan(RTSTAT)] = 0
                                    RTSTAT[R_Mean_change == 0] = 0
                                    DEC_array = DEC_array + RTSTAT
                                    if os.path.exists(os.path.join(OUT_folder, file.replace('mean_change','ttest_pvalue'))):
                                        logger.debug("- in vvvh ttest pvalue")
                                        RPVAL = raster2array(os.path.join(OUT_folder, file.replace('mean_change','ttest_pvalue')))
                                        RPVAL[RPVAL > THR_pvalue] = 2
                                        RPVAL[RPVAL <= THR_pvalue] = 1
                                        RPVAL[np.isnan(RPVAL)] = 2
                                        RPVAL[RPVAL == 2] = 0
                                        RPVAL[RTSTAT == 0] = 0
                                        DEC_array = DEC_array + RPVAL

                ######### APPLY PAST AND FUTURE FOREST/NON-FOREST/WATER MASK ##############
                for root, dirs, files in os.walk(OUT_folder):
                    for file in files:
                        if file.__contains__(str(date).replace('-', '')):
                            for POL in POLS:
                                if file.__contains__(POL + '_pmin_'):
                                    PMIN = raster2array(os.path.join(OUT_folder, file))
                                    PMIN[PMIN >= locals()[(POL + '_THR_pmin')]] = 1
                                    PMIN[PMIN < locals()[(POL + '_THR_pmin')]] = 0
                                    PMIN[np.isnan(PMIN)] = 1
                                    DEC_array = DEC_array*PMIN
                                if file.__contains__(POL + '_fmin_'):
                                    FMIN = raster2array(os.path.join(OUT_folder, file))
                                    FMIN[FMIN > locals()[(POL + '_THR_fmin')]] = 0
                                    FMIN[FMIN <= locals()[(POL + '_THR_fmin')]] = 1
                                    FMIN[np.isnan(FMIN)] = 1
                                    DEC_array = DEC_array*FMIN

                ################## APPLY BASE YEAR FOREST MASK #########################
                for POL in POLS:
                    BASE_Forestmask = raster2array(str(OUT_d_folder.joinpath(POL + '_BASEMIN_array.tif')))
                    BASE_Forestmask[BASE_Forestmask >= locals()[(POL + '_THR_pmin')]] = 1
                    BASE_Forestmask[BASE_Forestmask < locals()[(POL + '_THR_pmin')]] = 0
                    BASE_Forestmask[np.isnan(BASE_Forestmask)] = 1
                    DEC_array = DEC_array * BASE_Forestmask

                    BASE_Forestmask_STD = raster2array(str(OUT_d_folder.joinpath(POL + '_BASEMaxSTD_array.tif')))
                    BASE_Forestmask_STD[BASE_Forestmask_STD < THR_BASEMaxSTD] = 1
                    BASE_Forestmask_STD[BASE_Forestmask_STD >= THR_BASEMaxSTD] = 0
                    BASE_Forestmask_STD[np.isnan(BASE_Forestmask_STD)] = 1
                    DEC_array = DEC_array * BASE_Forestmask_STD

                ################## APPLY MIN STATCUBES CUTOFF FOR CHANGE ###############
                if date in last_days:
                    COUNT_THR = 0
                DEC_array_before_COUNT_THR = DEC_array
                DEC_array[DEC_array <= COUNT_THR] = 0
                DEC_array[DEC_array > COUNT_THR] = 1
                array2raster(MASTER_GRID_TIFF, str(OUT_d_folder.joinpath('DEC_array_' + str(date).replace('-', '') + '.tif')), DEC_array)
                ################## APPLY SIEVE FILTER ####################################
                OUT_d_file = str(OUT_d_folder.joinpath('DEC_array_' + str(date).replace('-', '') + '_s6-8LT.tif'))
                cmd_sieve = ["gdal_sieve.py", "-st",
                             "{}".format(str(sieve_size)),
                             "-8", "-nomask", "-of", "GTiff",
                             "{}".format(str(OUT_d_folder.joinpath('DEC_array_' + str(date).replace('-', '') + '.tif'))),
                             "{}".format(OUT_d_file)]
                cmd_output = subprocess.run(cmd_sieve, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                logger.debug("exit code {} --> {}".format(cmd_output.returncode, cmd_sieve))

                ################## APPLY GAP FILLING ####################################
                DEC_array_sieve_ORIGINAL = raster2array(OUT_d_file)
                DEC_array_sieve = raster2array(OUT_d_file)
                DEC_array_bool = DEC_array_sieve == 0
                DEC_DIS = dist(DEC_array_sieve,DEC_array_bool)

                array2raster(MASTER_GRID_TIFF, str(OUT_d_folder.joinpath('DEC_array_' + str(date).replace('-', '') + '_s6-8LT_INT.tif')), DEC_array_sieve)
                if os.path.exists(str(OUT_d_folder.joinpath('DEC_array_' + str(date).replace('-', '') + '.tif'))):
                    os.remove(str(OUT_d_folder.joinpath('DEC_array_' + str(date).replace('-', '') + '.tif')))
                if os.path.exists(str(OUT_d_folder.joinpath('DEC_array_' + str(date).replace('-', '') + '_s6-8LT.tif'))):
                    os.remove(str(OUT_d_folder.joinpath('DEC_array_' + str(date).replace('-', '') + '_s6-8LT.tif')))

                logger.debug("--- {} STAT CUBE class done".format(str(date).replace('-', '')))

        ################ DO REALITY CHECK ON EACH DATE ##############################
        logger.debug(HLINE)
        logger.debug("--- REMOVE DATES WITH POOR VALUES AND SUSPICIUS DETECTIONS ---")
        logger.debug(HLINE)

        FILES_VALUES = []
        if not OUT_d_folder.joinpath(tile + '_DEC_array_NOMASK' + '.tif').exists():
            for root, dirs, files in os.walk(OUT_d_folder):
                for file in files:
                    if file.__contains__('INT.tif') and (file.__contains__('_2017') or file.__contains__('_2018') or file.__contains__('_2019') or file.__contains__('_2020') or file.__contains__('_2021')):
                        logger.debug("processing file :{}".format(file))
                        MASTER_GRID_TIFF = os.path.join(OUT_d_folder, file)
                        FILE_VAL = raster2array(os.path.join(OUT_d_folder, file))
                        FILE_SUM = np.sum(FILE_VAL)
                        FILES_VALUES.append(FILE_SUM)
            logger.debug("Files values {}".format(FILES_VALUES))
            FILES_CUTOFF = np.mean(FILES_VALUES) + np.std(FILES_VALUES)
            logger.debug("Files cutoff {}".format(FILES_CUTOFF))

            for root, dirs, files in os.walk(OUT_d_folder):
                for file in files:
                    if file.__contains__('INT.tif') and (file.__contains__('_2017') or file.__contains__('_2018') or file.__contains__('_2019') or file.__contains__('_2020') or file.__contains__('_2021')):
                        logger.debug("checking file {}".format(file))
                        MASTER_GRID_TIFF = os.path.join(OUT_d_folder, file)
                        FILE_val = raster2array(os.path.join(OUT_d_folder, file))
                        FILE_SUM = np.sum(FILE_val)
                        logger.debug("file sum {}".format(FILE_SUM))
                        FILE_Z = (FILE_SUM - np.mean(FILES_VALUES)) /np.std(FILES_VALUES)
                        logger.debug("file Z {}".format(FILE_Z))

                        if (FILE_Z > 2.0):
                            array2raster(MASTER_GRID_TIFF, str(OUT_d_folder.joinpath(file.split(".tif")[0] + "_ERROR.tif")), FILE_val)
                            logger.debug("-marked error:{} --".format(file))
                            os.remove(os.path.join(OUT_d_folder, file))

        logger.debug("REMOVED SUSPICIUS DETECTIONS")
        logger.debug(HLINE)

        ################ CREATE ALL YEARS CHANGE MAP #################################
        logger.debug(HLINE)
        logger.debug("------- CREATE ALL YEARS CHANGE MAP -----------------")
        logger.debug(HLINE)
        YEARS = ['_2017','_2018','_2019', '_2020', '_2021']
        count = -1
        if not OUT_d_folder.joinpath(tile + '_DEC_array_NOMASK' + '.tif').exists():
            for root, dirs, files in os.walk(OUT_d_folder):
                for file in files:
                    if file.__contains__('INT.tif') and (file.__contains__('_2017') or file.__contains__('_2018') or file.__contains__('_2019') or file.__contains__('_2020') or file.__contains__('_2021')):
                        logger.debug("processing file :{}".format(file))
                        count += 1
                        if count == 0:
                            MASTER_GRID_TIFF = str(OUT_d_folder.joinpath(file))
                            logger.debug("Reading ... {}".format(file))
                            DEC_YEAR = raster2array(str(OUT_d_folder.joinpath(file)))
                        else:
                            logger.debug("Reading ... {}".format(file))
                            DEC = raster2array(str(OUT_d_folder.joinpath(file)))
                            DEC_YEAR = DEC_YEAR + DEC
            array2raster(MASTER_GRID_TIFF, str(OUT_d_folder.joinpath(tile + '_DEC_array_NOMASK' + '.tif')), DEC_YEAR)

        ###### APPLY OUR LARGE_WATER_BODIES, PRODES TERRABRAZILIS HYDRO LAYER, SRTM_GZ (OUT_water_file),and GFWTC TO ALL_YEAR MASK #########
        DEC_YEAR = raster2array(str(OUT_d_folder.joinpath(tile + '_DEC_array_NOMASK' + '.tif')))
        MASTER_GRID_TIFF = str(OUT_d_folder.joinpath(tile + '_DEC_array_NOMASK' + '.tif'))
        logger.debug("---- APPLY OUR LARGE_WATER_BODIES, SRTM_GZ, and PRODES TERRABRAZILIS HYDRO LAYER TO MASK -----")
        OUT_water_file = str(OUT_d_folder.joinpath(tile + '_WATER_mask' + '_s500-8LT.tif'))
        LARGE_WATER_BODIES = raster2array(OUT_water_file)
        DEC_YEAR[LARGE_WATER_BODIES == 1] = 0
        array2raster(MASTER_GRID_TIFF, str(OUT_d_folder.joinpath(tile + '_DEC_array' + '.tif')), DEC_YEAR)

        ################## CLEAN ALL YEARS CHANGE MAP - ITERATIVE CHANGE + WATER MASK #################################
        ########## THE FOR LOOP HERE IS JUST A PLACEHOLDER, IT IS NOT USED IN THIS VERSION OF THE SCRIPT ##############
        MASTER_GRID_TIFF = str(OUT_d_folder.joinpath(tile + '_DEC_array' + '.tif'))
        CHANGE_ORG = raster2array(str(OUT_d_folder.joinpath(tile + '_DEC_array' + '.tif')))
        OUT_water_file = str(OUT_d_folder.joinpath(tile + '_WATER_mask' + '_s500-8LT.tif'))
        DEC_YEAR = raster2array(str(OUT_d_folder.joinpath(tile + '_DEC_array' + '.tif')))
        for n in range(DISTANCE_to_WATER_THR):
            logger.debug('Running cleaning mask, iteration {}'.format(str(n)))
            if n==0:
                globals()[('CHANGE_3_' + str(n))] = raster2array(str(OUT_d_folder.joinpath(tile + '_DEC_array' + '.tif')))
                globals()[('LARGE_WATER_BODIES_' + str(n))] = raster2array(OUT_water_file)
            else:
                globals()[('CHANGE_3_' + str(n))] = globals()[('CHANGE_3_' + str(n-1))]
                globals()[('LARGE_WATER_BODIES_' + str(n))] = globals()[('LARGE_WATER_BODIES_' + str(n-1))]
            Bool_water = globals()[('LARGE_WATER_BODIES_' + str(n))] == 0
            DIS_water = dist(globals()[('LARGE_WATER_BODIES_' + str(n))],Bool_water)
            globals()[('LARGE_WATER_BODIES_' + str(n))][(DIS_water <= 1.0) & (globals()[('CHANGE_3_' + str(n))] > 0)] = 1
            globals()[('CHANGE_3_' + str(n))][DIS_water <= 1.0] = 0
            globals()[('CHANGE_3_' + str(n))][globals()[('CHANGE_3_' + str(n))] < NUM_Changes_THR] = 0
            globals()[('CHANGE_3_' + str(n))][globals()[('CHANGE_3_' + str(n))] >= NUM_Changes_THR] = 1
            CHANGE_3_bool = globals()[('CHANGE_3_' + str(n))] == 0
            CHANGE_3_DIS = dist(globals()[('CHANGE_3_' + str(n))],CHANGE_3_bool)
            globals()[('CHANGE_3_' + str(n))][(CHANGE_3_DIS <= 1.0) & ((CHANGE_ORG == 1.0) | (CHANGE_ORG == 2.0))] = 1
            globals()[('CHANGE_3_' + str(n))][globals()[('CHANGE_3_' + str(n))] == 1] = NUM_Changes_THR
            if n > 1:
                del globals()[('CHANGE_3_' + str(n-2))]
                del globals()[('LARGE_WATER_BODIES_' + str(n-2))]
        DEC_YEAR[globals()[('CHANGE_3_' + str(n))] == 0] = 0
        DEC_YEAR[globals()[('LARGE_WATER_BODIES_' + str(n))] == 1] = 0
        DEC_YEAR[DEC_YEAR >= 1] = 1
        array2raster(MASTER_GRID_TIFF, str(OUT_d_folder.joinpath(tile + '_DEC_array_MASK.tif')),DEC_YEAR)
        array2raster(MASTER_GRID_TIFF, str(OUT_d_folder.joinpath(tile + '_WATER_mask' + '_s500-8LT_C.tif')), globals()[('LARGE_WATER_BODIES_' + str(n))])

        ################## GENERATE AND APPLY LARGE-ERRORS SIEVE FILTER ####################################
        dec_mask_tifpath = str(Path(OUT_d_folder).joinpath(tile + '_DEC_array_MASK.tif'))
        dec_mask_largesieve_tifpath = str(Path(OUT_d_folder).joinpath(tile + '_DEC_array_MASK_largesieve.tif'))
        cmd_sieve = ["gdal_sieve.py", "-st", "{}".format(50000),
                     "-8", "-nomask", "-of", "GTiff",
                     "{}".format(dec_mask_tifpath),
                     "{}".format(dec_mask_largesieve_tifpath)]
        cmd_output = subprocess.run(cmd_sieve, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logger.debug("exit code {} --> {}".format(cmd_output.returncode, cmd_sieve))

        CHANGE_VALIDITY = raster2array(dec_mask_tifpath)
        CHANGE_VALIDITY_largesieve = raster2array(dec_mask_largesieve_tifpath)
        CHANGE_VALIDITY[CHANGE_VALIDITY_largesieve == 1] = 0
        array2raster(MASTER_GRID_TIFF, dec_mask_tifpath, CHANGE_VALIDITY)
        if os.path.exists(dec_mask_largesieve_tifpath):
            os.remove(dec_mask_largesieve_tifpath)

        ################# APPLY ALL YEARS CHANGE MAP TO EACH DATE ################
        logger.debug(HLINE)
        logger.debug("--- Apply all years change map to each data ---")
        logger.debug(HLINE)

        CHANGE_VALIDITY = raster2array(str(OUT_d_folder.joinpath(tile + '_DEC_array_MASK.tif')))
        for root, dirs, files in os.walk(OUT_d_folder):
            for file in files:
                if file.__contains__('INT.tif') and not file.__contains__('_2015') and not file.__contains__('_2016'):
                    MASTER_GRID_TIFF = os.path.join(OUT_d_folder, file)
                    logger.debug("Validating pixels in: {}".format(file))
                    DEC = raster2array(os.path.join(OUT_d_folder, file))
                    DEC[CHANGE_VALIDITY == 0] = 0
                    array2raster(MASTER_GRID_TIFF, str(OUT_d_folder.joinpath(file.split('.tif')[0] + '_C.tif')), DEC)
                    os.remove(os.path.join(OUT_d_folder, file))

        ################# APPLY ALL YEARS CHANGE MAP TO NO_MASK LAYER ################
        MASTER_GRID_TIFF = str(Path(OUT_d_folder).joinpath(tile + '_DEC_array_MASK.tif'))
        NO_M = raster2array(str(Path(OUT_d_folder).joinpath(tile + '_DEC_array.tif')))
        DEC_LAYER = raster2array(str(Path(OUT_d_folder).joinpath(tile + '_DEC_array_MASK.tif')))
        NO_M[DEC_LAYER == 0] = 0
        NO_M[(NO_M == 0) & (DEC_LAYER == 1)] = 1
        array2raster(MASTER_GRID_TIFF, str(Path(OUT_d_folder).joinpath(tile + '_DEC_array_C.tif')), NO_M)

    logger.debug(HLINE)
    logger.debug("--- Merging files in one/two orbit(s) ---")
    logger.debug(HLINE)
    OUT_d_folder = output_dir.joinpath(DETECTION_SET, tile)

    if set(ORBIT_DIRECTIONS) == set(['descending','ascending']):
        files_orbit0_list = []
        files_orbit1_list = []
        files_orbits_list = []
        for root, dirs, files in os.walk(Path(OUT_d_folder).joinpath(ORBIT_DIRECTIONS[0])):
            for file in files:
                if file.__contains__('INT_C.tif') or file.__contains__('DEC_array_C.tif') or file.__contains__('DEC_array_MASK.tif') or file.__contains__('DEC_array_NOMASK.tif') or file.__contains__('DEC_array.tif'):
                   files_orbit0_list.append(file)
        logger.debug("orbit 0 files: {}".format(sorted(files_orbit0_list)))
        for root, dirs, files in os.walk(Path(OUT_d_folder).joinpath(ORBIT_DIRECTIONS[1])):
            for file in files:
                if file.__contains__('INT_C.tif') or file.__contains__('DEC_array_C.tif') or file.__contains__('DEC_array_MASK.tif') or file.__contains__('DEC_array_NOMASK.tif') or file.__contains__('DEC_array.tif'):
                   files_orbit1_list.append(file)
        logger.debug("orbit 1 files: {}".format(sorted(files_orbit1_list)))

        files_orbit0_set = set(files_orbit0_list)
        files_orbit1_set = set(files_orbit1_list)
        files_orbits_set = files_orbit0_set.union(files_orbit1_set)
        files_orbits_list = list(files_orbits_set)
        logger.debug("Files to merge or transfer to detection")
        logger.debug("{}".format(files_orbits_list))

        for files_orbits_list_item in files_orbits_list:
            files_orbits_list_item_path_orbit0 = str(Path(OUT_d_folder).joinpath(ORBIT_DIRECTIONS[0], files_orbits_list_item))
            files_orbits_list_item_path_orbit1 = str(Path(OUT_d_folder).joinpath(ORBIT_DIRECTIONS[1], files_orbits_list_item))
            files_orbitscombined_tile = str(Path(OUT_d_folder).joinpath(files_orbits_list_item))

            if os.path.exists(files_orbits_list_item_path_orbit0) and os.path.exists(files_orbits_list_item_path_orbit1):
               logger.debug("Asc and Des has {}".format(files_orbits_list_item))
               FILE_1 = raster2array(files_orbits_list_item_path_orbit0)
               FILE_2 = raster2array(files_orbits_list_item_path_orbit1)
               if not files_orbits_list_item_path_orbit0.__contains__('DEC_array_C.tif') and not files_orbits_list_item_path_orbit1.__contains__('DEC_array_C.tif'):
                  FILE_C = np.nansum(np.dstack((FILE_1, FILE_2)), 2)
                  FILE_C[FILE_C > 1] = 1
               else:
                  FILE_C = np.nanmean(np.dstack((FILE_1, FILE_2)), 2)
               array2raster(files_orbits_list_item_path_orbit0, files_orbitscombined_tile, FILE_C)
            elif os.path.exists(files_orbits_list_item_path_orbit0) and not os.path.exists(files_orbits_list_item_path_orbit1):
               logger.debug("Des has {}".format(files_orbits_list_item))
               copyfile(files_orbits_list_item_path_orbit0, files_orbitscombined_tile)
            elif not os.path.exists(files_orbits_list_item_path_orbit0) and os.path.exists(files_orbits_list_item_path_orbit1):
               logger.debug("Asc has {}".format(files_orbits_list_item))
               copyfile(files_orbits_list_item_path_orbit1, files_orbitscombined_tile)
            else:
               logger.error("This shouldnt be happening: {} doesnt meet the conditions".format(files_orbits_list_item))
               raise Exception("Check merge part of code")

    else:
        for root, dirs, files in os.walk(Path(OUT_d_folder).joinpath(ORBIT_DIRECTIONS[0])):
            for file in files:
                if file.__contains__('INT_C.tif') or file.__contains__('DEC_array_C.tif') or file.__contains__('DEC_array_MASK.tif') or file.__contains__('DEC_array_NOMASK.tif') or file.__contains__('DEC_array.tif'):
                   files_orbits_list_item_path_orbit0 = str(Path(OUT_d_folder).joinpath(ORBIT_DIRECTIONS[0], file))
                   files_orbitscombined_tile = str(Path(OUT_d_folder).joinpath(file))
                   copyfile(files_orbits_list_item_path_orbit0, files_orbitscombined_tile)


    for root, dirs, files in os.walk(Path(OUT_d_folder)):
        for file in files:
            if file.__contains__('DEC_array_C.tif'):
                MASTER_GRID_TIFF = str(os.path.join(OUT_d_folder,file))
                CONF_LAYER = raster2array(str(os.path.join(OUT_d_folder,file)))
                CONF_LAYER[CONF_LAYER>10] = 10
                CONF_LAYER = CONF_LAYER*10
                conf_filepath = Path(OUT_d_folder).joinpath(tile + '_DEC_array_CONF.tif')
                array2raster(MASTER_GRID_TIFF, str(conf_filepath), CONF_LAYER)

    stop_time = time.time()
    seconds = stop_time - start_time
    minutes = seconds/60
    hours = minutes/60
    processing_time = 'Processing time: %02d:%02d:%02d - %s [s] \n' % (hours, minutes % 60, seconds % 60, stop_time - start_time)
    logger.debug("Processing_time is {}".format(processing_time))
    logger.debug("Deforestation detection finished successfully")

if __name__ == "__main__":

        parser = argparse.ArgumentParser(description="Deforestation detection tool.")
        parser.add_argument("--in_tile", type=str, help="Tile name.", required=True, default=None)
        parser.add_argument("--in_timeperiod", type=str, help="time period.", required=True, default=None)
        parser.add_argument("--statcubes_folder", type=str, help="StatCubes folder inside output folder.",
                            required=True, default=None)
        parser.add_argument("--output_dir", type=str, help="Output directory.", required=True, default=None)
        parser.add_argument("--aux_dir", type=str, help="Auxiliary directory.", required=True, default=None)
        args = parser.parse_args()

        run(tile=args.in_tile,
            timeperiod=args.in_timeperiod,
            statcubes_dir=args.statcubes_folder,
            output_dir=args.output_dir,
            aux_dir=args.aux_dir
            )
