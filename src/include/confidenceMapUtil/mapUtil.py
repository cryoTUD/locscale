from . import FDRutil
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import argparse, os, sys
import subprocess
import math
import gc
import os.path
from time import sleep


#Author: Maximilian Beckers, EMBL Heidelberg, Sachse Group (2017)


#------------------------------------------------------------------------------------------------------
def localFiltration(map, locResMap, apix, localVariance, windowSize, boxCoord, ECDF):

	#**************************************************
	#**** function to perform a local filtration ******
	#****** according to the local resolution *********
	#**************************************************

	#some initialization
	mapSize = map.shape;
	numX = mapSize[0];
	numY = mapSize[1];
	numZ = mapSize[2];
	
	mean = np.zeros((numX, numY, numZ));
	var = np.zeros((numX, numY, numZ));
	ECDFmap = np.ones((numX, numY, numZ));
	filteredMapData = np.zeros((numX, numY, numZ));

	#transform to numpy array
	locResMapData = np.copy(locResMap);

	#set all resoltuon lower than 2.1 to 2.1
	#locResMapData[locResMapData > 2.5] = 2.5;
	
	locResMapData[locResMapData == 0.0] = 100.0;
	locResMapData[locResMapData >= 100.0] = 100.0;

	#transform to abosulte frequency units(see http://sparx-em.org/sparxwiki/absolute_frequency_units)
	locResMapData = np.divide(apix, locResMapData);
	
	#round to 3 decimals
	locResMapData = np.around(locResMapData, 3);	

	#set resolution search range, 3 decimals exact
	locResArray = np.arange(0, 0.5+0.001 , 0.001);
	
	#set maximum resolution, important as ResMap is masking
	limRes = np.min(locResMapData);
	counter = 0;
	numRes = len(locResArray);	

	#get initial noise statistics
	initMapData = np.copy(map);
	initMean, initVar, _ = FDRutil.estimateNoiseFromMap(initMapData, windowSize, boxCoord);
	noiseMapData = np.random.normal(initMean, math.sqrt(initVar), (100, 100, 100));

	#do FFT of the respective map
	mapFFT = np.fft.rfftn(map);

	#get frequency map
	frequencyMap = FDRutil.calculate_frequency_map(map);

	# Initial call to print 0% progress
	#printProgressBar(counter, numRes, prefix = 'Progress:', suffix = 'Complete', bar_length = 50)
	print("Start local filtering. This might take a few minutes ...");

	counterRes = 0;
	for tmpRes in locResArray:   
		counterRes = counterRes + 1;
		progress = counterRes/float(numRes);
		if counterRes%(int(numRes/20.0)) == 0:
			output = "%.1f" %(progress*100) + "% finished ..." ;
			print(output);
		
		#get indices of voxels with the current resolution	
		indices = np.where(locResMapData == tmpRes);
	
		if (indices[0].size == 0):
			#this resolution is obviously not in the map, so skip
			counter = counter + 1;
			continue;
		elif math.fabs(tmpRes - limRes) < 0.0000001:
			xInd, yInd, zInd = indices[0], indices[1], indices[2];
			
			#do local filtration
			tmpFilteredMapData = FDRutil.lowPassFilter(mapFFT, frequencyMap, tmpRes, map.shape);

			#set the filtered voxels
			filteredMapData[xInd, yInd, zInd] = tmpFilteredMapData[xInd, yInd, zInd];

		else:
			xInd, yInd, zInd = indices[0], indices[1], indices[2];
			#do local filtration
			tmpFilteredMapData = FDRutil.lowPassFilter(mapFFT, frequencyMap, tmpRes, map.shape);
			#set the filtered voxels
			filteredMapData[xInd, yInd, zInd] = tmpFilteredMapData[xInd, yInd, zInd];
			if localVariance == True:
				#estimate and set noise statistic

				if ECDF == 1:
					#if ecdf shall be used, use if to p-vals
					tmpECDF, sampleSort = FDRutil.estimateECDFFromMap(tmpFilteredMapData, windowSize, boxCoord);
					vecECDF = np.interp(tmpFilteredMapData[xInd, yInd, zInd], sampleSort, tmpECDF, left=0.0, right=1.0);
					ECDFmap[xInd, yInd, zInd] = vecECDF; 
				else:
					ECDFmap = 0;

				tmpMean, tmpVar, _ = FDRutil.estimateNoiseFromMap(tmpFilteredMapData, windowSize, boxCoord);
				mean[xInd, yInd, zInd] = tmpMean;
				var[xInd, yInd, zInd] = tmpVar;

	print("Local filtering finished ...");

	return filteredMapData, mean, var, ECDFmap;

#------------------------------------------------------------------------------------------------
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, bar_length = 100):
     
    #******************************************
    #** progress bar for local visualization **
    #******************************************	
    """
    Call in a loop to create terminal progress bar
    params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    str_format = "{0:." + str(decimals) + "f}";
    percents = str_format.format(100 * (iteration / float(total)));
    filled_length = int(round(bar_length * iteration / float(total)));
    bar = '#' * filled_length + '-' * (bar_length - filled_length);

    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix));

    if iteration == total:
        sys.stdout.write('\n');
    sys.stdout.flush();

    return;
#---------------------------------------------------------------------------------
def makeCircularMask(map, sphereRadius):

	#some initialization
	mapSize = map.shape;

	x = np.linspace(-math.floor(mapSize[0]/2.0), -math.floor(mapSize[0]/2.0) + mapSize[0], mapSize[0]);
	y = np.linspace(-math.floor(mapSize[1]/2.0), -math.floor(mapSize[1]/2.0) + mapSize[1], mapSize[1]);
	z = np.linspace(-math.floor(mapSize[2]/2.0), -math.floor(mapSize[2]/2.0) + mapSize[2], mapSize[2]);

	xx, yy, zz = np.meshgrid(x, y, z, indexing='ij');

	mask = np.sqrt(xx**2 + yy**2 + zz**2);

	mask[mask>sphereRadius] = 0.0;

	mask[mask>0.0] = 1.0;

	return mask;
	


