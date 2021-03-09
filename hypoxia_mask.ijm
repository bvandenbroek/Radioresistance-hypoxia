
saveSettings();

print("\\Clear");
run("Clear Results");
run("Input/Output...", "jpeg=85 gif=-1 file=.csv save_column");

var nr_series;
var file_name;
var format;
var default_pixel_size = 0.4900;

var median_radius = 5;
var min_size = 100;				//Minimum size of the spheroid (in units^2)
var use_hematoxylin = true;		//Use hematoxylin staining for spheroid detection? Otherwise pimo is used.
var distance_map_scale = 4;		//scaling factor for determining the distance map. A value of 4 allows a maximum distance of 4x255 = 1020 um from the periphery
var cluster_nr = 1;				//Number of the cluster that contains the hypoxic region (k-means method can sometimes switch cluster numbers).
var all_files = false;

var zone_width = 10;	//Width of the zones in microns
var fill_holes = true;

if(nImages>0) run("Close All");
path = File.openDialog("Select any file in the folder to be processed");

setBatchMode(true);

run("Bio-Formats Macro Extensions");
//run("Bio-Formats Importer", "open=["+path+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_1");
Ext.setId(path);
Ext.getFormat(path, format);
//print(format);

//getDimensions(width, height, channels, slices, frames);

file_name = File.getName(path);
dir = File.getParent(path)+"\\";
savedir= dir+"\\results\\";
if(!File.exists(savedir)) File.makeDirectory(savedir);

extension = substring(file_name, lastIndexOf(file_name,".")+1);

file_list = getFileList(dir); //get filenames of directory
//make a list of images with 'extension' as extension.
j=0;
image_list=newArray(file_list.length);	//Dynamic array size doesn't work on some computers, so first make image_list the maximal size and then trim.
image_list_no_extension=newArray(file_list.length);
for(i=0; i<file_list.length; i++){
	if (endsWith(file_list[i],extension)) {
		image_list[j] = file_list[i];
		image_list_no_extension[j] = substring(file_list[i], 0, lastIndexOf(file_list[i],"."));
		j++;
	}
}
image_list = Array.trim(image_list, j);	//Trimming the array of images
image_list_no_extension = Array.trim(image_list_no_extension, j);	//Trimming the array of images

//Dialog
Dialog.create("Settings");
Dialog.addNumber("Minimum spheroid size",min_size,0,4,"um");	//there is no check for units yet, so we just assume it is um.
Dialog.addCheckbox("Use hematoxylin staining for spheroid detection?",use_hematoxylin);
Dialog.addNumber("Width of the zones",zone_width,0,4,"um");
Dialog.addCheckbox("Fill holes in hypoxic area?",fill_holes);
Dialog.addSlider("Hypoxic zone cluster nr",1,2,cluster_nr);
Dialog.addCheckbox("Process all "+image_list.length+" ."+extension+" files in this folder?",all_files);
Dialog.show();
min_size = Dialog.getNumber();
use_hematoxylin = Dialog.getCheckbox();
zone_width = Dialog.getNumber();
fill_holes = Dialog.getCheckbox();
cluster_nr = Dialog.getNumber();
all_files = Dialog.getCheckbox();

print("\\Clear");
print("Directory contains "+file_list.length+" files, of which "+image_list.length+" ."+extension+" ("+format+") files.");
run("Set Measurements...", "area mean standard integrated median area_fraction stack limit redirect=None decimal=3");
run("Clear Results");
roiManager("Reset");

current_image_nr=0;
do {
	if(all_files==true) {
		run("Close All");
		file_name = image_list[current_image_nr];		//retrieve file name from image list
	}
	else file_name = File.getName(path);
	Ext.setId(dir+file_name);
	Ext.getSeriesCount(nr_series);

	for(i=0;i<nr_series;i++) {
		print("Processing file "+current_image_nr+1+"/"+image_list.length+", series "+i+1+"/"+nr_series+": "+dir+file_name+"...");
		
		run("Bio-Formats Importer", "open=["+dir+file_name+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_"+i);
		//open(dir+file_name);
		getPixelSize(unit, pixelWidth, pixelHeight);
		if(unit == "microns" && pixelWidth>5) {
			print("WARNING: Pixel calibration is probably wrong. Reverting to the default value of "+default_pixel_size+" um/pixel.");
			run("Properties...", "unit=um pixel_width="+default_pixel_size+" pixel_height="+default_pixel_size+" voxel_depth=1.0000000");
		}
		else if(unit == "inches" || unit == "cm" && pixelWidth>0.0005) {
			print("WARNING: Pixel calibration is probably wrong. Reverting to the default value of "+default_pixel_size+" um/pixel.");
			run("Properties...", "unit=um pixel_width="+default_pixel_size+" pixel_height="+default_pixel_size+" voxel_depth=1.0000000");
		}
		getPixelSize(unit, pixelWidth, pixelHeight);
		
		image = getTitle();
		image = replace(image,"/","-");	//replace slashes by dashes in the name
		image = substring(image, 0,lengthOf(image)-lengthOf(extension)-1);	//remove extension
		rename(image);

		roiManager("Reset");
		run("Clear Results");

		//waitForUser("Select Pimo image and click OK");
		if(bitDepth!=24) {
			run("RGB Color");
			close(image);
			rename(image);
		}
		setBatchMode("show");
		image = getTitle();
		run("Colour Deconvolution", "vectors=[H AEC] hide");
		close(image+"-(Colour_3)");
		selectWindow(image+"-(Colour_2)");
		run("Median...", "radius="+median_radius);
		//run("Gaussian Blur...", "sigma=10");
		run("8-bit");
		//setAutoThreshold("Default");
		//run("Threshold...");
		run("k-means Clustering ...", "number_of_clusters=3 cluster_center_tolerance=0.00010000 enable_randomization_seed randomization_seed=48");
		setThreshold(1,2);
		if(use_hematoxylin == true) {
			selectWindow(image+"-(Colour_1)");
			run("Median...", "radius="+median_radius);
			setAutoThreshold("Otsu");
		}
		run("Analyze Particles...", "size="+min_size+"-Infinity include add");

		//Create combined ROI in case of multiple selections
		if(roiManager("count")>1) {
			roiManager("Select All");
			roiManager("Combine");
			roiManager("Reset");
			roiManager("Add");
		}
		
		//Ask for manual ROI and merge with automatic ROI
		selectWindow(image);
		roiManager("Show All without labels");
		setTool("freehand");
		waitForUser("Select region of interest around the spheroid or just click OK to select the whole image.");
		if(selectionType == -1) run("Select All");
		else if(selectionType<0 || selectionType>3) {
			do {
				wait(50);
				waitForUser("Please select a single area.");
				wait(50);
			} while (selectionType<0 || selectionType>3);
		}
		roiManager("Add");

		roiManager("Select All");
		roiManager("AND");
		roiManager("Reset");
		roiManager("Add");
		roiManager("Select",0);
		roiManager("Rename", "Spheroid");

		run("Create Mask");
		rename("mask_spheroid");
		
		run("Select None");

		//detect the necrotic part
		selectWindow(image+"-(Colour_1)");
		run("Select None");
		setAutoThreshold("Mean"); //Yen also works nicely
		run("Create Selection");
		roiManager("add");
		roiManager("Select All");
		roiManager("XOR");
		roiManager("Add");
		roiManager("Select",1);
		roiManager("Delete");
		roiManager("Select",1);
		roiManager("Rename", "Necrotic+outside");
		roiManager("Select All");
		roiManager("AND");
		roiManager("Add");
		roiManager("Select",1);
		roiManager("Delete");
		roiManager("Select",1);
		roiManager("Rename", "Necrotic");

		//Create distance map by subsequent downscaling and upscaling
		selectWindow("mask_spheroid");
		run("Scale...", "x="+1/distance_map_scale+" y="+1/distance_map_scale+" interpolation=None average create");
		setThreshold(1, 255);
		run("Convert to Mask");		//make binary again
		run("Distance Map");
		run("16-bit");				//convert to 16-bit to allow larger values than 255
		run("Multiply...", "value="+distance_map_scale*pixelWidth);	//distance in units
		run("Scale...", "x="+distance_map_scale+" y="+distance_map_scale+" interpolation=Bilinear average create");
		resetMinAndMax();
		rename("distance_map");

		//Define hypoxic area
		selectWindow("Clusters");
		run("Select None");
		setThreshold(cluster_nr,cluster_nr);
		run("Convert to Mask");
		if(fill_holes==true) {
			run("Fill Holes");
		}
		run("Create Selection");
		roiManager("add");
		roiManager("Select",roiManager("count")-1);
		roiManager("Rename", "Hypoxia");
		run("Create Mask");
		rename("mask_hypoxia");
		//Make necrotic area NaN
		run("32-bit");
		roiManager("Select",1);
		changeValues(0, 0, NaN);

		//Create and empty results table with no additional columns
		run("Measure");
		selectWindow("Results");
		run("Close");

		run("Set Measurements...", "area mean standard area_fraction redirect=None decimal=3");
		selectWindow("distance_map");
		getMinAndMax(min, max);
		if((max/zone_width)%1 !=0) count = floor(max/zone_width)+1;
		else count = floor(max/zone_width);					//if maximum distance is just at the end of the zone.
		for(z=0;z<count;z++) {
				selectWindow("distance_map");
				setThreshold(z*zone_width+1 , (z+1)*(zone_width));
				//print("threshold: "+z*zone_width+1 + " - "+(z+1)*(zone_width));
				run("Create Selection");
				changeValues(0, 99999, (z+1)*zone_width);
				roiManager("add");
				roiManager("Select",roiManager("count")-1);
				roiManager("Rename", "zone_"+z+1);
				getRawStatistics(total_area);
				roiManager("Deselect");
				selectWindow("mask_hypoxia");
				roiManager("Select",roiManager("count")-1);
				List.setMeasurements();
				area_perc = List.getValue("%Area");
				//viable_area = List.getValue("Area");
				getRawStatistics(viable_area);
				setResult("zone ("+zone_width+" um)", z, z+1);
				setResult("% area hypoxia", z, area_perc);
				setResult("total area ("+unit+"^2)", z, total_area*pixelWidth);
				setResult("viable area ("+unit+"^2)", z, viable_area*pixelWidth);
				setResult("% area viable", z, viable_area/total_area*100);
		}
		selectWindow("distance_map");
		resetThreshold();
		run("Select None");
		run("Fire");
		resetMinAndMax();
		run("Select None");
		setBatchMode("show");
		roiManager("Show None");
	}
	close("mask_spheroid");
	close("mask_spheroid-1");
	close("clusters");
	close("mask_hypoxia");
	close(image+"-(Colour_1)");
	close(image+"-(Colour_2)");


	//save images
	roiManager("save", savedir+image+"_ROIs.zip");
	selectWindow("distance_map");
	saveAs("Tiff", savedir+image+"_distance_map");

	selectWindow(image);
	setLineWidth(3);
	roiManager("Select", 0);	//spheroid
	setForegroundColor(255,0,0);
	run("Draw");
	roiManager("Select", 1);	//hypoxia
	setForegroundColor(0,0,255);
	run("Draw");
	roiManager("Deselect");
	saveAs("Tiff", savedir+image+"_overlay");

	updateResults();
	//Save results
	saveAs("Results", savedir+image+"_results.txt");

	current_image_nr++;
	if(current_image_nr<image_list.length && all_files == true) waitForUser("Press OK to process the next image.");

} while (current_image_nr<image_list.length && all_files==true);


restoreSettings();

