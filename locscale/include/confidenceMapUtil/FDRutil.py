import numpy as np
import subprocess
import math
import gc
import os
import sys


#Author: Maximilian Beckers, EMBL Heidelberg, Sachse Group

#-------------------------------------------------------------------------------------
def estimateNoiseFromMap(map, windowSize, boxCoord):

	#**************************************************
	#****** function to estimate var an mean from *****
	#**** nonoverlapping boxes outside the particle ***
	#**************************************************

	if boxCoord == 0:
		#extract a sample of pure noise from the map
		sizeMap = map.shape;
		sizePatch = np.array([windowSize, windowSize, windowSize]);
		center = np.array([0.5*sizeMap[0], 0.5*sizeMap[1], 0.5*sizeMap[2]]);
		
		sampleMap1 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		int(0.02*sizeMap[1]):(int(0.02*sizeMap[1]) + sizePatch[1]),
		(int(center[2]-0.5*sizePatch[2])):(int((center[2]-0.5*sizePatch[2]) + sizePatch[2]))];

		sampleMap2 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		int(0.98*sizeMap[1] - sizePatch[1]):(int(0.98*sizeMap[1])),
		(int(center[2]-0.5*sizePatch[2])):(int((center[2]-0.5*sizePatch[2]) + sizePatch[2]))];
	
		sampleMap3 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		(int(center[1]-0.5*sizePatch[1])):(int((center[1]-0.5*sizePatch[1]) + sizePatch[1])), 
		int(0.02*sizeMap[2]):(int(0.02*sizeMap[2]) + sizePatch[2])];

		sampleMap4 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		(int(center[1]-0.5*sizePatch[1])):(int((center[1]-0.5*sizePatch[1]) + sizePatch[1])), 
		int(0.98*sizeMap[2]) - sizePatch[2]:(int(0.98*sizeMap[2]))];

		#concatenate the two samples
		sampleMap = np.concatenate((sampleMap1, sampleMap2, sampleMap3, sampleMap4), axis=0);

	else:
		sizePatch = np.array([windowSize, windowSize, windowSize]);
		center = np.array(boxCoord);
		sampleMap = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		int(center[1]-0.5*sizePatch[1]):(int(center[1]-0.5*sizePatch[1]) + sizePatch[1]),
		(int(center[2]-0.5*sizePatch[2])):(int((center[2]-0.5*sizePatch[2]) + sizePatch[2]))];	
	
	#estimate variance and mean from the sample
	mean = np.mean(sampleMap);
	var = np.var(sampleMap);

	if var == 0.0:
		print("Variance is estimated to be 0. You are probably estimating noise in a masked region. Exit ...")
		sys.exit();

	return mean, var, sampleMap;
#-------------------------------------------------------------------------------------
def estimateNoiseFromMapInsideMask(map, mask):

	#**************************************************
	#****** function to estimate var an mean from *****
	#******* map outside the user provided mask *******
	#**************************************************

	mask[mask<=0.5] = 0.0;
	mask[mask>0.0] = 1000.0;
	mask[mask<1000.0] = 1.0;
	mask[mask==1000.0] = 0.0;

	sampleMap = np.copy(map)*mask;
	sampleMap = sampleMap[sampleMap != 0.0];
	
	#estimate variance and mean from the sample
	mean = np.mean(sampleMap);
	var = np.var(sampleMap);
	
	
	return mean, var, sampleMap;

#-------------------------------------------------------------------------------------
def estimateECDFFromMap(map, windowSize, boxCoord):

	#**************************************************
	#****** function to estimate empirical cumul. *****
	#**** distribution function from solvent area *****
	#**************************************************

	if boxCoord == 0:
		#extract a sample of pure noise from the map
		sizeMap = map.shape;
		sizePatch = np.array([windowSize, windowSize, windowSize]);
		center = np.array([0.5*sizeMap[0], 0.5*sizeMap[1], 0.5*sizeMap[2]]);
		sampleMap1 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		int(0.02*sizeMap[1]):(int(0.02*sizeMap[1]) + sizePatch[1]),
		(int(center[2]-0.5*sizePatch[2])):(int((center[2]-0.5*sizePatch[2]) + sizePatch[2]))];

		sampleMap2 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		int(0.98*sizeMap[1] - sizePatch[1]):(int(0.98*sizeMap[1])),
		(int(center[2]-0.5*sizePatch[2])):(int((center[2]-0.5*sizePatch[2]) + sizePatch[2]))];
	
		sampleMap3 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		(int(center[1]-0.5*sizePatch[1])):(int((center[1]-0.5*sizePatch[1]) + sizePatch[1])), 
		int(0.02*sizeMap[2]):(int(0.02*sizeMap[2]) + sizePatch[2])];

		sampleMap4 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		(int(center[1]-0.5*sizePatch[1])):(int((center[1]-0.5*sizePatch[1]) + sizePatch[1])), 
		int(0.98*sizeMap[2]) - sizePatch[2]:(int(0.98*sizeMap[2]))];

		#conatenate the two samples
		sampleMap = np.concatenate((sampleMap1, sampleMap2, sampleMap3, sampleMap4), axis=0);
	
	elif boxCoord == -1:
		sampleMap = map;
	
	else:
		sizePatch = np.array([windowSize, windowSize, windowSize]);
		center = np.array(boxCoord);
		sampleMap = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		int(center[1]-0.5*sizePatch[1]):(int(center[1]-0.5*sizePatch[1]) + sizePatch[1]),
		(int(center[2]-0.5*sizePatch[2])):(int((center[2]-0.5*sizePatch[2]) + sizePatch[2]))];	
	
	#estimate ECDF from map
	sampleMap = sampleMap.flatten();	
	
	#downsize the sample
	finalSampleSize = min(100000, sampleMap.size);
	sampleMap = np.random.choice(sampleMap, finalSampleSize, replace = False);
	numSamples = sampleMap.size;
	sampleMapSort = np.sort(sampleMap);

	minX = sampleMapSort[0];
	maxX = sampleMapSort[numSamples-1];
	numInterval = numSamples;
	spacingX = (maxX - minX)/(float(numInterval));

	ECDF = np.zeros(numInterval);	
	for index in range(numInterval):
		val = sampleMapSort[index];
		ECDF[index] = ((sampleMapSort[sampleMapSort<= val]).size)/float(numSamples);
	
	return ECDF, sampleMapSort;

#------------------------------------------------------------------------------------
def getCDF(x, ECDF, sampleMapSort):

	#****************************************************
	#********* get the value of the CDF at point x ******
	#******* CDF : Cumulative distribution function *****
	#****************************************************
	
	numSamples = sampleMapSort.size;
	minX = sampleMapSort[0];
	maxX = sampleMapSort[numSamples-1];


	if x >= maxX:
		CDFval = 1.0;
	elif x <= minX:
		CDFval = 0.0;
	else:
		#get the index in the ECDF array		
		index = np.searchsorted(sampleMapSort, x) - 1;	
		CDFval = ECDF[index];

	return CDFval;
#------------------------------------------------------------------------------------
def AndersonDarling(sample):

	#********************************************
	#*** Anderson-Darling test for normality ****
	#********************************************

	sample = np.random.choice(sample, min(10000,sample.size), replace=False);
	sampleMapSort = np.sort(sample);
	numSamples = sampleMapSort.size;

	Ad = -1.0*numSamples;
	for i in range(numSamples):
		CDF_Yi = 0.5 * (1.0 + math.erf(sampleMapSort[i]/math.sqrt(2.0)));
		CDF_Yn = 0.5 * (1.0 + math.erf(sampleMapSort[numSamples-i-1]/math.sqrt(2.0)));
		if CDF_Yi == 0:
			CDF_Yi = 0.000001;
		if CDF_Yi == 1:
			CDF_Yi = 0.999999;
		if CDF_Yn == 0:
			CDF_Yn = 0.000001;
		if CDF_Yn == 1:
			CDF_Yn = 0.999999;

		#calculate the Anderson-Darling test statistic
		Ad = Ad - (1.0/float(numSamples)) * (2*(i+1)-1)*(math.log(CDF_Yi) + (math.log(1.0-CDF_Yn)));

	#do adjustment for estimation of mean and variance, as unknown before
	Ad = Ad*(1 + 0.75/float(numSamples) + 2.25/float(numSamples*numSamples));

	#calculate p-values
	# R.B. D'Augostino and M.A. Stephens, Eds., 1986, Goodness-of-Fit Techniques, Marcel Dekker
	try:
		if Ad >= 0.6:
			pVal = math.exp(1.2937 - 5.709*(Ad) + 0.0186*Ad*Ad);
		elif 0.34<Ad<0.6:
			pVal = math.exp(0.9177 - 4.279*Ad - 1.38 * Ad*Ad);
		elif 0.2 < Ad <= 0.34:
			pVal = 1 - math.exp(-8.318 + 42.796*Ad - 59.938*Ad*Ad);
		else:
			pVal = 1 - math.exp(-13.436 + 101.14 * Ad - 223.73 * Ad*Ad);
	except:
		pVal = -1.0;

	return Ad, pVal, numSamples;


#------------------------------------------------------------------------------------
def KolmogorowSmirnow(ECDF, sampleMapSort):

	#***********************************************
	#***** KS test by supremum of distance *********
	#*********** between CDF and ECDF **************
	#***********************************************

	#some initialization
	numSamples = sampleMapSort.size;
	X = np.linspace(-5, 5, 200000);
	vectorizedErf = np.vectorize(math.erf);

	#maximum distances between CDF and ECDF over the whole defintion region
	Y_stdNorm = 0.5 * (1.0 + vectorizedErf(X/math.sqrt(2.0)));
	Y_ecdf = np.interp(X, sampleMapSort, ECDF, left=0.0, right=1.0);
	Dn = np.amax(np.absolute(np.subtract(Y_stdNorm,Y_ecdf)));

	#get Kolmogorow-Smirnow test statistic
	KS_testStat = math.sqrt(numSamples)*Dn;

	#maximum distances between CDF and ECDF for tail regions
	X_tail_right = X[X>2.0];
	X_tail_left = X[X<-2.0];
	X_tail = np.concatenate((X_tail_right, X_tail_left));
	Y_stdNorm = 0.5 * (1.0 + vectorizedErf(X_tail/math.sqrt(2.0)));
	Y_ecdf = np.interp(X_tail, sampleMapSort, ECDF, left=0.0, right=1.0);
	Dn_tail = np.amax(np.absolute(np.subtract(Y_stdNorm,Y_ecdf)));

	return KS_testStat, Dn, Dn_tail, numSamples;

#-----------------------------------------------------------------------------------
def checkNormality(map, windowSize, boxCoord):

	#***************************************
	#** check normal distribution ass. *****
	#***************************************

	print('Checking the normal distribution assumption ...');

	mean, var, _ = estimateNoiseFromMap(map, windowSize, boxCoord);

	map = np.subtract(map, mean);
	tMap = np.multiply(map, (1.0/(math.sqrt(var))));
	map = np.copy(tMap);
		
	#get maximum distances between ECDF and CDF
	ECDFvals, sampleSort = estimateECDFFromMap(map, windowSize, boxCoord);
	KSstat, Dn, Dn_tail, n = KolmogorowSmirnow(ECDFvals, sampleSort);
	output = "Maximum Distance Dn between ECDF and CDF: Dn=" + " %.4f" %Dn + ", in Tail:" + " %.4f" %Dn_tail + ". Sample size used: " + repr(n);
	print(output);

	#do Anderson-Darling test for normality
	AnDarl, pVal, n = AndersonDarling(sampleSort);
	output = "Anderson-Darling test summary: " + repr(AnDarl) + ". p-Value: " + "%.4f" %pVal + ". Sample size used: " + repr(n);
	if pVal != -1.0:
		print(output);
	else:
		pVal = -1.0;		

	if (Dn_tail > 0.01):
		output = "WARNING: Deviation in the tail areas between the normal distribution and the empircal CDF is higher than 1%. If boxes for background noise estimation are set properly, please consider using the flag -ecdf to use the empirical CDF instead of the normal distribution."
		print(output);

#------------------------------------------------------------------------------------
def studentizeMap(map, mean, var):
	
	#****************************************
	#********* normalize map ****************
	#****************************************
	
	if np.isscalar(var):
		studMap = np.subtract(map, mean);
		studMap = np.multiply(studMap, (1.0/(math.sqrt(var))));
	else: #if local variances are known, use them
		var[var == 0] = 1000;
		studMap = np.subtract(map, mean);
		studMap = np.divide(studMap, np.sqrt(var));
		var[var == 1000] = 0.0;
		studMap[var == 0.0] = 0.0;
		
	return studMap;


#-----------------------------------------------------------------------------------
def calcQMap(map, mean, var, ECDF, windowSize, boxCoord, mask, method, test):

	#*****************************************
	#***** generate qMap of a 3D density *****
	#*****************************************

	#get some map data
	sizeMap = map.shape;

	#calculate the test statistic
	if np.isscalar(var):
		#map[map == 0.0] = -100000000;
		tmap = np.subtract(map, mean);   
		tMap = np.multiply(tmap, (1.0/(math.sqrt(var))));
		map = np.copy(tMap);
	else:
		var[var==0.0] = 1000.0; #just to avoid division by zero
		tmap = np.subtract(map, mean); 	
		tMap = np.divide(tmap, np.sqrt(var));
		var[var==1000.0] = 0.0;
		#upadte the mask, necessary as resmap is masking as well
		mask = np.multiply(mask,var);

	if test == 'rightSided':
		tMap[tMap<0] = -10000000.0;
	elif test == 'leftSided':
		tMap[tMap>0] = 10000000.0;

	#calculate the p-Values
	print('Calculating p-Values ...');
	pMap = np.zeros(sizeMap);
	
	if np.isscalar(ECDF): 
		if ECDF == 0: 
			vectorizedErf = np.vectorize(math.erf);
			erfMap = vectorizedErf(tMap/math.sqrt(2.0));
			#erf2Map = special.erf(tMap/math.sqrt(2.0));

			pMapRight = 1.0 - (0.5*(1.0 + erfMap)); 
			pMapLeft = (0.5*(1.0 + erfMap));
	
		else:
			#if ecdf shall be used, use if to p-vals
			ECDF, sampleSort = estimateECDFFromMap(map, windowSize, boxCoord);
			print('start ECDF calculation ...');	
			vecECDF = np.interp(map, sampleSort, ECDF, left=0.0, right=1.0);
			pMapRight = 1.0 - vecECDF;
			pMapLeft = vecECDF;

			#************************************
			#****** make diagnostic plots *******
			#************************************
			X = np.linspace(-5, 5, 10000);
			vectorizedErf = np.vectorize(math.erf);
			Y1 = np.interp(X, sampleSort, ECDF, left=0.0, right=1.0);
			Y2 = (0.5*(1.0 + vectorizedErf(X/math.sqrt(2.0))));
			np.savetxt('ECDF.txt', Y1, delimiter=';');
			import matplotlib.pyplot as plt
			from matplotlib.backends.backend_pdf import PdfPages
			import matplotlib.gridspec as gridspec
			plt.rc('xtick', labelsize=8);    # fontsize of the tick labels
			plt.rc('ytick', labelsize=8);    # fontsize of the tick labels
			pp = PdfPages('ECDF_vs_CDF.pdf');
			gs = gridspec.GridSpec(1, 2);
			
			ax1 = plt.subplot(gs[0]);
			ax1.set_title('Full CDF');
			ax1.plot(X, Y1);
			ax1.plot(X, Y2);
	    	#do test plots
			X = np.linspace(2, 5, 10000);
			Y1 = np.interp(X, sampleSort, ECDF, left=0.0, right=1.0);
			Y2 = (0.5*(1.0 + vectorizedErf(X/math.sqrt(2.0))));
			ax2 = plt.subplot(gs[1]);
			ax2.set_title('Tail CDF');
			ax2.plot(X, Y1);
			ax2.plot(X, Y2);
			pp.savefig();
			pp.close();
			plt.close();
	else:
		pMapRight = 1 - ECDF;
		pMapLeft = ECDF;

	if test == 'twoSided':
		pMap = np.minimum(pMapLeft, pMapRight);
	elif test == 'rightSided':
		pMap = pMapRight;
		pMap[tMap==-10000000.0] = 1.0
	elif test == 'leftSided':
		pMap = pMapLeft;
		pMap[tMap==10000000.0] = 1.0;
	
	#take the p-values in the mask	
	binaryMask = np.copy(mask);
	binaryMask[ binaryMask != 0.0 ] = 1.0;
	binaryMask[ binaryMask ==  0.0] = np.nan;

	pMap = pMap * binaryMask;
	pFlat = pMap.flatten();
	pInBall = pFlat[~np.isnan(pFlat)];

	#do FDR control, i.e. calculate the qMap
	print('Start FDR control ...');
	pAdjValues = pAdjust(pInBall, method);	

	pAdjFlat = np.copy(pFlat);
	pAdjFlat[~np.isnan(pAdjFlat)] = pAdjValues;
	pAdjMap = np.reshape(pAdjFlat, sizeMap);
	pAdjMap[np.isnan(pAdjMap)] = 1.0;

	return pAdjMap;

#---------------------------------------------------------------------------------
def pAdjust(pValues, method):

	#***********************************
	#***** FDR correction methods ******
	#***********************************

	numPVal = len(pValues);

	print("Sorting p-values ...")
	pSortInd = np.argsort(pValues);
	pSort = pValues[pSortInd];
	
	pAdjust = np.zeros(numPVal);
	prevPVal = 1.0;

	#use expansion for harmonic series
	Hn = math.log(numPVal) + 0.5772 + 0.5/numPVal - 1.0/(12*numPVal**2) + 1.0/(120*numPVal**4);    

	print("Adjusting p-values ...");
	if method =='BH': #do benjamini-hochberg procedure
		for i in range(numPVal-1, -1, -1):
			pAdjust[i] = min(prevPVal, pSort[i]*numPVal/(i+1.0));
			prevPVal = pAdjust[i];	
	
	elif method == 'BY': #do benjamini yekutieli procedure
		for i in range(numPVal-1, -1, -1):
			pAdjust[i] =  min(prevPVal, pSort[i]*(numPVal/(i+1.0))*Hn);
			prevPVal = pAdjust[i];
	elif method == "Holm":
		prevPVal = 0.0;
		for i in range(numPVal):
			tmpPVal = (numPVal - i)*pSort[i];
			pAdjust[i] = max(prevPVal, tmpPVal);
			prevPVal = pAdjust[i];
		pAdjust[pAdjust>1.0] = 1.0;
	elif method == "Hochberg":
		for i in range(numPVal-1, -1, -1):
			pAdjust[i] = min(prevPVal, pSort[i]*(numPVal-i));
			prevPVal = pAdjust[i];		
	else:
		print('Please specify a method. Execution is stopped ...');
		quit();	

	#sort back to the original order
	pSortIndOrig = np.argsort(pSortInd);

	return pAdjust[pSortIndOrig];

#---------------------------------------------------------------------------------
def binarizeMap(map, threshold):

	#***********************************
	#*binarize map at given threshold **
	#***********************************

	binMap = np.array(map);
	binMap[binMap <= threshold] = 0;
	binMap[binMap > threshold] = 1;

	#finally invert the map
	binMap = np.subtract(np.ones(binMap.shape), binMap);
  	
	return binMap;

#---------------------------------------------------------------------------------
def tanh_filter(x, cutoff):

	#**********************************
	#******* tanh lowpass filter ******
	#**********************************

	filter_fall_off = 0.1;

	fx = 1.0 - (1.0 - 0.5*(np.tanh((np.pi*(x+cutoff)/(2*filter_fall_off*cutoff))) - np.tanh((np.pi*(x-cutoff)/(2*filter_fall_off*cutoff)))));

	return fx;
#---------------------------------------------------------------------------------
def calculate_frequency_map(map):

    sizeMap = map.shape;

    #calc frequency for each voxel
    freqi = np.fft.fftfreq(sizeMap[0], 1.0);
    freqj = np.fft.fftfreq(sizeMap[1], 1.0);
    freqk = np.fft.rfftfreq(sizeMap[2], 1.0);

    sizeFFT = np.array([freqi.size, freqj.size, freqk.size]);
    FFT = np.zeros(sizeFFT);

    freqMapi = np.copy(FFT);
    for j in range(sizeFFT[1]):
        for k in range(sizeFFT[2]):
            freqMapi[:,j,k] = freqi*freqi;

    freqMapj = np.copy(FFT);
    for i in range(sizeFFT[0]):
        for k in range(sizeFFT[2]):
            freqMapj[i,:,k] = freqj*freqj;

    freqMapk = np.copy(FFT);
    for i in range(sizeFFT[0]):
        for j in range(sizeFFT[1]):
            freqMapk[i,j,:] = freqk*freqk;
    
    frequencyMap = np.sqrt(freqMapi + freqMapj + freqMapk);

    return frequencyMap;

#---------------------------------------------------------------------------------
def lowPassFilter(mapFFT, frequencyMap, cutoff, shape):

	#**********************************
	#*** filter in fourier domain *****
	#**********************************

	sizeMap = mapFFT.shape;

	#do filtering of the map
	filterMap = tanh_filter(frequencyMap, cutoff);
	filteredftMap = filterMap*mapFFT;
	
	#do iverse FFT
	try: #optional with pyfftw
		import pyfftw
		import multiprocessing

		# get number of cpus
		numCores = multiprocessing.cpu_count();

		fftObject = pyfftw.builders.irfftn(filteredftMap, shape, threads = numCores);
		filteredMap = fftObject();

	except:
		filteredMap = np.fft.irfftn(filteredftMap, shape);


	filteredMap = np.real(filteredMap);

	return filteredMap;

#---------------------------------------------------------------------------------
def printSummary(args, time):

	#***********************************
	#** print a Summary of the job *****
	#***********************************

	print("*** Done ***");
	print(" ");
	print("******** Summary ********");

	output = "Elapsed Time: " + repr(time);
	print(output);

	#print input map filename
	splitFilename = os.path.splitext(os.path.basename(args.em_map));
	output = "Input EM-map: " + splitFilename[0] + ".mrc"; 
	print(output);

	#print model map filename
	if args.model_map is not None:
		splitFilename = os.path.splitext(os.path.basename(args.model_map));
		output = "LocScale was done with the input model-map: " + splitFilename[0] + ".mrc";
		print(output);

		#print window size
		if args.window_size is not None:
			w = int(math.ceil(args.window_size / 2.) * 2);
		else:
			w = int(round(7 * 3 * args.apix));
		output = "Window size: " + repr(w);

	#print local resolution map filename
	if args.locResMap is not None:
		splitFilename = os.path.splitext(os.path.basename(args.locResMap));
		output = "Input local resolution map: " + splitFilename[0] + ".mrc";
		print(output);
	
	#print output filenames
	if args.outputFilename is not None:
		splitFilename = os.path.splitext(os.path.basename(args.outputFilename));
	else:
		splitFilename = os.path.splitext(os.path.basename(args.em_map));
	
	#print specifications for locScale and local filtering
	if args.model_map is not None:
		output = "Output LocScale map: " + splitFilename[0] + "_scaled" + ".mrc";
		print(output);
		output = "Output confidence Map: " + splitFilename[0] + "_confidenceMap" + ".mrc";
		print(output);
		
		if args.stepSize is not None:
			output = "Step size used for the sliding window moves: " + repr(args.stepSize);
		else:
			output = "Step size used for the sliding window moves: " + repr(5);
		print(output);

	else:
		if args.locResMap is not None:
			output = "Output locally filtered map: " + splitFilename[0] + "_locFilt" + ".mrc";
			print(output);
			output = "Output confidence Map: " + splitFilename[0] + "_confidenceMap" + ".mrc";
			print(output);
		else:
			output = "Output confidence Map: " + splitFilename[0] + "_confidenceMap" + ".mrc";
			print(output);

	#print method used for FDR-control
	if args.method is None:
		output = "Multiple Testing was controlled with: BY";			
	else:
		output = "Multiple Testing  was controlled with: " + args.method;
	print(output);


