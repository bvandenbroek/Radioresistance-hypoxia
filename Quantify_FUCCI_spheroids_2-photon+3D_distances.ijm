/*
 * Macro to quantify the green and red signals from fluorescent FUCCI cells in 3D spheroids
 * 
 * Brief workflow:
 * - Rolling ball background subtraction is performed on both green and red channels
 * - The green channel is normalized to the mean intensity of the (thresholded) red channel stack, in order to equalize the two signal strengths
 * - The two channels are added to combine all the cells into a single 3D stack.
 * - StarDist nuclei segmentation is performed in separate 2D slices of this stack.
 * - The outline of each slice is determined by thresholding the 2D image, after sequentially applying a median filter and Gaussian blurred 2D.
 * - The 2D distances from the centroids of all nuclei to the spheroid section rim are computed, 
 *   and then recalculated to the shortest 3D distances to the spheroid edge, assuming a spherical shape.
 * - For every segmented nucleus, the green and red intensity g and r are measured.
 * - Positive green cells are defined as normalized green signal g/(g+r) being larger than a user-defined relative threshold.
 * - A histogram of the fractions of green and red cells vs distance is constructed by grouping them in distance bins.
 * 
 * Required Fiji update sites:
 * - StarDist
 * - CSBDeep
 * - CLIJ
 * - CLIJ2
 * - SCF MPI CBG
 * 
 * Author: Bram van den Broek, Netherlands Cancer Institute, December 2020 - March 2021
 * b.vd.broek@nki.nl
 *
 * N.B. Filtering causes mistakes with the 3D calculations. It's not really necessary anyway. :-) 
 * 
 */

#@ File (label = "Input file", style = "file") input
#@ File (label = "Output directory", style = "directory") output
#@ Integer (label = "Green channel", value = 1) chGreen
#@ Integer (label = "Red channel", value = 2) chRed
#@ Integer (label = "Rolling ball radius (background subtraction)", value = 100) rollingBallRadius
#@ Float (label = "StarDist nucleus probability threshold", maxvalue = 0.3) probabilityThreshold
#@ Integer (label = "Remove nulei with diameter smaller than (um)", value = 4) minNucleusSize_setting
#@ Integer (label = "Remove nulei with diameter larger than (um)", value = 40) maxNucleusSize_setting
#@ Integer (label = "Start Z slice", min = 1, value = 1) startZ
#@ Integer (label = "Z step size", value = 1) stepZ
#@ Integer (label = "End Z slice", value = 99) endZ
#@ Integer (label = "Crop edge to circumvent stitching artifacts (pixels)", value = 20) cropBorder
#@ Integer (label = "Blur size for spheroid detection", value = 5) sigma
#@ Integer (label = "Distance zone step size (um)", value = 10) binSize
#@ Integer (label = "Maximum distance (um)", value = 200) maxDistance
#@ Integer (label = "Estimated spheroid radius (um)", value = 200) spheroidRadius
#@ Integer (label = "Slice at the bottom of the spheroid", value = 0) offset
#@ Float (label = "Relative threshold intensity for green cells", min=0, max=1, value = 0.5) thresholdIntensityGreen
#@ Boolean (label = "Only re-analyze Results with new threshold for green cells", value=false) reanalyze


nBins = maxDistance/binSize;
minDistance = 0;

filename = File.getNameWithoutExtension(input);
//CLose previous tables
close(filename + "_greenFraction.tsv");

saveSettings();

if(reanalyze==false) {

run("Conversions...", " ");
run("Colors...", "foreground=white background=black selection=gray");
run("Set Measurements...", "area mean standard min median stack redirect=None decimal=3");
roiManager("reset");

run("Close All");
open(input);


setBatchMode(true);

original = getTitle();
getDimensions(width, height, channels, slices, frames);
if(cropBorder>0) {
	makeRectangle(cropBorder, cropBorder, width-2*cropBorder, height-2*cropBorder);
	wait(50);
	run("Crop");
}
getDimensions(width, height, channels, slices, frames);
getVoxelSize(pixelWidth, pixelHeight, voxelDepth, unit);
minNucleusSize = PI*pow((minNucleusSize_setting / pixelWidth / 2),2);	//Calculate the nucleus area as if it were a circle
maxNucleusSize = PI*pow((maxNucleusSize_setting / pixelWidth / 2),2);

green = "green";
red = "red";
selectWindow(original);
setBatchMode("show");
run("Duplicate...", "title="+green+" duplicate channels="+chGreen);
run("Green");
run("32-bit");
run("Subtract Background...", "rolling="+rollingBallRadius+" stack");
selectWindow(original);
run("Duplicate...", "title="+red+" duplicate channels="+chRed);
run("Red");
run("32-bit");
run("Subtract Background...", "rolling="+rollingBallRadius+" stack");

greenInt = getIntensity(green);
redInt = getIntensity(red);
ratioGR = greenInt/redInt;
print("-------------------------------");
print("Mean green intensity: "+greenInt);
print("Mean red intensity: "+redInt);
print("Adjusting green intensity to match red intensity, with factor "+1/ratioGR);
selectWindow(green);
run("Divide...", "value="+ratioGR+" stack");	//normalize green to red intensity
imageCalculator("Add create stack", green,red);
rename("combined_cells");
run("Grays");
Stack.setSlice(slices/2);
run("Enhance Contrast", "saturated=0.35");
Stack.setSlice(1);
//setBatchMode("show");
//run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");	//Make it look like a time-lapse for StarDist

//Initialize GPU
run("CLIJ2 Macro Extensions", "cl_device=");
Ext.CLIJ2_clear();
// In case another GPU needs to be selected:
//Ext.CLIJ2_listAvailableGPUs();
//availableGPUs = Table.getColumn("GPUName");
//run("CLIJ2 Macro Extensions", "cl_device=" + availableGPUs[1]);

//Analyze slice by slice

MeasurementTable = "Measurements";
Table.create(MeasurementTable);
Table.setLocationAndSize(100, 100, 300, 600);
Table.reset(MeasurementTable);

for(z=startZ; z<=minOf(slices,endZ); z+=stepZ) {
	showProgress((z-startZ)/(endZ-startZ));
	selectWindow("combined_cells");
	Stack.setSlice(z);
	run("Duplicate...", "title=slice_"+z);
	slice = getTitle();
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'slice_"+z+"', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.60000000000001', 'probThresh':'"+probabilityThreshold+"', 'nmsThresh':'0.3', 'outputType':'ROI Manager', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Stack', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");

	//rename ROIs
	for(i=0; i<roiManager("count"); i++) {
		roiManager("select",i);
		roiManager("rename",i+1);
	}
	run("Select None");
	
	run("ROI Manager to LabelMap(2D)");
	labelmap = getTitle();
	Ext.CLIJ2_push(labelmap);
	close(labelmap);

	//Display ROIs on combined cells image during processing
	selectWindow(original);
	Stack.setSlice(z);
	roiManager("Show All without Labels");
	
	//Filter on area
	run("Clear Results");
	Ext.CLIJ2_statisticsOfBackgroundAndLabelledPixels(labelmap, labelmap); //Somehow if you put (image, labelmap) as arguments the pixel count is wrong
	Ext.CLIJ2_pushResultsTableColumn(area, "PIXEL_COUNT");

	Ext.CLIJ2_excludeLabelsWithValuesOutOfRange(area, labelmap, labelmap_filtered, minNucleusSize, maxNucleusSize);
	Ext.CLIJ2_release(labelmap);


	Ext.CLIJ2_centroidsOfLabels(labelmap_filtered, pointlist);
	//Ext.CLIJ2_pull(pointlist);
	//setBatchMode("show");

	// make another similar image which contains the coordinates + another line with intensities (1, 2, ... n) 
	Ext.CLIJ2_getDimensions(pointlist, number_of_labels, dimensionality, _);
	Ext.CLIJ2_create2D(coordinates_and_index, number_of_labels + 1, dimensionality + 1, 32);
	Ext.CLIJ2_setRampX(coordinates_and_index);
	Ext.CLIJ2_paste2D(pointlist, coordinates_and_index, 1, 0);
	
	// generate an output image, set it to 0 everywhwere
	Ext.CLIJ2_getDimensions(labelmap_filtered, width, height, depth);
	Ext.CLIJ2_create2D(pointmap, width, height, 16);
	Ext.CLIJ2_set(pointmap, 0);
	
	// at every pixel position defined in the coordinate list above, write a number
	Ext.CLIJ2_writeValuesToPositions(coordinates_and_index, pointmap);
	
	/* visualize the pointmap output
	Ext.CLIJ2_pull(pointmap);
	setBatchMode("show");
	setMinAndMax(0, number_of_labels);
	run("glasbey_on_dark");
	*/

	Ext.CLIJ2_release(pointlist);
	Ext.CLIJ2_release(coordinates_and_index);

	// Generate spheroid distance map
	selectWindow(slice);
	run("Duplicate...", "title=median_"+z);
	median = "median_"+z;
	run("Median...", "radius="+sigma);	//Perform median on the CPU, because the GPU version often crashes.
	Ext.CLIJ2_push(median);
	close(median);
//	Ext.CLIJ2_median2DBox(slice, median, sigma, sigma);
//	Ext.CLIJ2_release(slice);	
	Ext.CLIJ2_gaussianBlur2D(median, blurred, sigma, sigma);
	Ext.CLIJ2_release(median);
	Ext.CLIJ2_thresholdHuang(blurred, mask);
	Ext.CLIJ2_binaryFillHoles(mask, mask_filled);
	Ext.CLIJ2_pullBinary(mask_filled);
//	run("EDM Binary Operations", "iterations="+sigma+" operation=close stack");	//BioVoxxel Toolbox close (better than the built-in one)
//To do: maybe first erode a bit
	rename("spheroid_mask_"+z);

	//Estimate 2D spheroid radius
	setAutoThreshold("Default dark");
	run("Create Selection");
	run("Fit Circle");
	getSelectionBounds(x, y, sliceDiameter, sliceDiameter);
	resetThreshold();
	run("Fill");
	run("Select None");
	
	Ext.CLIJ2_distanceMap(mask_filled, distancemap);
	Ext.CLIJ2_release(mask_filled);

	run("Clear Results");
	Ext.CLIJ2_statisticsOfBackgroundAndLabelledPixels(distancemap, pointmap);
	distance_ = Table.getColumn("MEAN_INTENSITY","Results");
	distance_ = Array.deleteIndex(distance_, 0);	//remove the background value (0) in the array
	distance_ = multiplyArraywithScalar(distance_, pixelWidth);		//multiply with pixel size
	Ext.CLIJ2_release(pointmap);
	Ext.CLIJ2_release(distancemap);

	//Create a visual distance map of all nuclei
	setResult("MEAN_INTENSITY", 0, 0);	//Set background distance to 0
	Ext.CLIJ2_pushResultsTableColumn(distance_vector, "MEAN_INTENSITY");
//Ext.CLIJ2_pull(distance_vector);
//	Ext.CLIJ2_generateParametricImage(labelmap_filtered, distance_vector, distance_map_nuclei);
//	Ext.CLIJ2_replaceIntensities(labelmap_filtered, distance_vector, distance_map_nuclei);
//	Ext.CLIJ2_pull(distance_map_nuclei);
//	setBatchMode("show");

	//Measure green
	Ext.CLIJ2_pull(labelmap_filtered);
	selectWindow(green);
	Stack.setSlice(z);
	run("Duplicate...", "title=green_"+z);
	run("Intensity Measurements 2D/3D", "input=green_"+z+" labels="+labelmap_filtered+" mean median");	//Measure with MorphoLibJ
	green_ = Table.getColumn("Mean", "green_"+z+"-intensity-measurements");
	close("green_"+z);
	close("green_"+z+"-intensity-measurements");

	//Measure red
	selectWindow(red);
	Stack.setSlice(z);
	run("Duplicate...", "title=red_"+z);
	run("Intensity Measurements 2D/3D", "input=red_"+z+" labels="+labelmap_filtered+" mean median");	//Measure with MorphoLibJ
	red_ = Table.getColumn("Mean", "red_"+z+"-intensity-measurements");
	close("red_"+z);
	close("red_"+z+"-intensity-measurements");

	close(labelmap_filtered);

	cell_nr_ = Array.getSequence(distance_.length+1);
	cell_nr_ = Array.deleteIndex(cell_nr_, 0);	//remove element 0

	selectWindow(MeasurementTable);
	n=0;
	tableSize = Table.size;
	for(i=tableSize; i<(tableSize+distance_.length); i++) {
		Table.set("slice",i ,z, MeasurementTable);
		Table.set("nucleus", i, cell_nr_[n], MeasurementTable);
		Table.set("distance_2D", i, distance_[n], MeasurementTable);
		//Estimate 3D distance to the edge, assuming a spherical structure, sitting on the bottom (z=1)
		distance_[n] = maxOf(0, spheroidRadius - sqrt(pow((sliceDiameter/2*pixelWidth - distance_[n]),2) + pow((spheroidRadius - (z-offset)*voxelDepth),2)) );
		Table.set("distance", i, distance_[n], MeasurementTable);
		Table.set("green", i, green_[n], MeasurementTable);
		Table.set("red", i, red_[n], MeasurementTable);
		Table.set("g/(g+r)", i, green_[n]/(red_[n]+green_[n]), MeasurementTable);
		n++;
	}
	Table.update;
	
	Array.rotate(distance_, 1);
	distance_[0] = 0;
	Ext.CLIJ2_pushArray(distance_vector, distance_, distance_.length, 1, 1);
	Ext.CLIJ2_replaceIntensities(labelmap_filtered, distance_vector, distance_map_nuclei);
	Ext.CLIJ2_pull(distance_map_nuclei);
	rename("distance_map_nuclei_"+z);	//These images will later be merged into the analyzed image
	run("Fire");
	selectWindow(slice);

	Ext.CLIJ2_clear();
}
run("Images to Stack", "method=[Copy (center)] name=Combined_nuclei_analyzed_slices title=slice_ use");
run("16-bit");
run("Images to Stack", "method=[Copy (center)] name=Distance_map_nuclei title=distance_map_nuclei use");
run("Images to Stack", "method=[Copy (center)] name=Spheroid_mask title=mask use");
run("16-bit");
setBatchMode("exit and display");

Table.rename(MeasurementTable,"Results");

run("Merge Channels...", "c3=Combined_nuclei_analyzed_slices c4=Spheroid_mask c5=Distance_map_nuclei create");
getDimensions(width, height, channels, slices, frames);
Stack.setSlice(slices/2);
Stack.setChannel(3);
setMinAndMax(0, maxDistance);
Stack.setChannel(2);
setMinAndMax(0, 1523);
run("Green");
Stack.setChannel(1);
saveAs("Tiff", output + File.separator + filename + "_analyzed");

} //End of image analysis - hereafter only data analysis


allDistances = Table.getColumn("distance", "Results");
allGreenIntensityFractions = Table.getColumn("g/(g+r)", "Results");

distanceHist = Array.getSequence(nBins);
distanceHist = multiplyArraywithScalar(distanceHist, binSize);
distanceHist = addScalarToArray(distanceHist, binSize/2);
totalCells = newArray(nBins);
greenCells = newArray(nBins);

for (i=0 ; i<nResults ; i++) {
	distance = getResult("distance", i);
	currentBin = floor(((distance - minDistance)/(maxDistance - minDistance))*nBins);
	if(distance>=maxDistance) currentBin = nBins-1;
	else if(distance<=minDistance) currentBin = 0;
	if(allGreenIntensityFractions[i] > thresholdIntensityGreen) greenCells[currentBin] += 1;
	totalCells[currentBin] += 1;
}

fractionGreen = divideArrays(greenCells, totalCells);

Table.create("Green fraction");
Table.setLocationAndSize(100, 100, 400, 600);
Table.setColumn("distance", distanceHist);
Table.setColumn("green cells", greenCells);
Table.setColumn("total cells", totalCells);
Table.setColumn("fraction green cells", fractionGreen);

//Save results
if(reanalyze == true) filename = File.getNameWithoutExtension(input);
selectWindow("Results");
saveAs("Results", output + File.separator + filename + "_allResults.tsv");
selectWindow("Green fraction");
saveAs("Results", output + File.separator + filename + "_greenFraction.tsv");

// Create plot
Plot.create("Green cell fraction", "distance (um)", "fraction");
Plot.add("circle", allDistances, allGreenIntensityFractions, "Green intensity fraction");
thresholdLine = newArray(distanceHist.length);
Array.fill(thresholdLine, thresholdIntensityGreen);
Plot.add("line",distanceHist,thresholdLine, "Green cell threshold");
Plot.add("circle",distanceHist, fractionGreen, "Fraction cells with green signal");
Plot.setStyle(0, "#77a0ff,none,2.0,Dot");
Plot.setStyle(1, "#a00000,none,2.0,Line");
//Plot.setStyle(2, "#00a000,green,2.0,Connected");
Plot.setStyle(2, "#00a000,none,3.0,Separated Bars");
Plot.setLimits(0, maxDistance, 0, 1);
Plot.setFrameSize(600, 300);
Plot.addLegend("Green fractional intensity [g/(g+r)]\nRelative threshold for green signal\nFraction cells with green signal", "Top-Right Bottom-To-Top");
//saveAs("Tiff", output + File.separator + filename + "_Plot");


restoreSettings;


function getIntensity(image) {
	selectWindow(image);
	getDimensions(width, height, channels, slices, frames);
	Stack.setSlice(slices/2);
	setAutoThreshold("Otsu dark stack");
	mean = getValue("Mean limit");
	resetThreshold;
	return mean;
}

function show(string) {
	setBatchMode("show");
	waitForUser(string);
}

//Multiplies all elements of an array with a scalar
function multiplyArraywithScalar(array, scalar) {
	multiplied_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		multiplied_array[a]=array[a]*scalar;
	}
	return multiplied_array;
}

//Divides the elements of two arrays and returns the new array
function divideArrays(array1, array2) {
	divArray=newArray(lengthOf(array1));
	for (a=0; a<lengthOf(array1); a++) {
		divArray[a]=array1[a]/array2[a];
	}
	return divArray;
}

//Adds a scalar to all elements of an array
function addScalarToArray(array, scalar) {
	added_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		added_array[a]=array[a] + scalar;
	}
	return added_array;
}
