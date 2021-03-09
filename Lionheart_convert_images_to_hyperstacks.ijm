/*
 * Macro to bundle separate .tif files from the Lionheart microscope into hyperstacks
 * output: A single multi-page tif file per well, with as prefix the folder name.
 * 
 * The number of slices, channels, frames and tiles is automatically obtained from the file names.
 * (This is easier than reading from the metadata.)
 * Additionally, it pads the well names with a leading zero (e.g. B2 -> B02)
 * 
 * Author: Bram van den Broek, The Netherlands Cancer Institute (b.vd.broek@nki.nl)
 * 
 * version 1.0, March 2019
 * version 2.0, August 2020
 * 	- Optimized code for metadata (z,c,f) reading
 * 	- Added objective magnification
 * 	- Does not rename raw files on the disk any more
 * 
 */

#@ File (label = "Input directory", style = "directory") input
#@ String (label = "File suffix", value = ".tif") suffix
#@ String (label = "Objective magnification", choices={10,20,40,60,100}) objectiveMagnification
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "output file name prefix (leave empty to inherit the folder name)", value = "experimentName") outputPrefix

pixelSize = 12.9/objectiveMagnification;	//For our Lionheart system (?)

if(output == input) {
	showMessage("Warning: Output directory cannot be the same as Input directory. Creating a subdirectory 'Output' in "+input);
	File.makeDirectory(input + File.separator + "Output");
	output = input + File.separator + "Output"; 
}

//Initialize
z=0;
f=0;
c=0;

list = getFileList(input);
list = Array.sort(list);
wellNames = newArray(list.length);
w=0;

setBatchMode(true);
print("\\Clear");

//Determine max values for wells, z, c, anf f.
for (i = 0; i < list.length; i++) {
	if(endsWith(list[i], suffix)) {
		if(!File.isDirectory(input + File.separator + list[i])) {
			name = list[i];
			metadata = split(name, "_.");	//Split the name on both '_' and '.'
			if(i==0) wellNames[w] = metadata[0];	//first well name
			if(i>0) {
				if (metadata[0] != wellNames[w]) {	//compare well names
					w++;
					wellNames[w] = metadata[0];
				}
			}

			c = maxOf(parseInt(metadata[2]),c);
			z = maxOf(metadata[3],z);
			f = maxOf(metadata[5],f);
		}
	}
}
w = w+1;

wellNames = Array.trim(wellNames,w);
Array.print(wellNames);
print("\\Clear");
print("Detected properties from file names:");
print("nr of wells = "+w);
print("z-slices = "+z);
print("frames = "+f);
print("channels = "+c);

//The actual processing
for(i=0;i<w;i++) {
	convertWell(input, output, wellNames[i]);
	showProgress(i/wellNames.length);
}
print("\nProcessing finished.");



function convertWell(input, output, well) {
	//Importing images, creating a hyperstack and saving as a multi-page tiff file

	run("Image Sequence...", "open=["+input+"] file="+well+"_ sort");

	//Fix pixel calibration
	run("Properties...", "unit=um pixel_width="+pixelSize+" pixel_height="+pixelSize+" voxel_depth=1");

	tiles = nSlices/(z*f*c);	//determine the number of tiles
	//print("tiles = "+tiles);
	
	run("Stack to Hyperstack...", "order=xytzc channels="+c+" slices="+z+" frames="+f*tiles+" display=Color");
	outputName = outputPrefix;
	name = getTitle();
	print("processing well "+well+"...");
	if (outputPrefix == "") outputName = name + "_" + substring(well,0,1) + IJ.pad(substring(well,1), 2);
	else outputName = outputPrefix + "_" + substring(well,0,1) + IJ.pad(substring(well,1), 2);
	stack = getTitle;
	if(tiles>1) {
		for(j=1;j<=tiles;j++) {
			selectWindow(stack);
			if(c==1 && z==1) run("Duplicate...", "duplicate range="+f*(j-1)+1+"-"+f*j);
			else  run("Duplicate...", "duplicate frames="+f*(j-1)+1+"-"+f*j);
			
//			if(c==1 && z==1) run("Make Substack...", "channels=1-"+c+" slices=1-"+z+" frames="+f*(j-1)+1+"-"+f*j);
//			else if (z==1) run("Make Substack...", "channels=1-"+c+" frames="+f*(j-1)+1+"-"+f*j);
//			else if (c==1) run("Make Substack...", "slices=1-"+z+" frames="+f*(j-1)+1+"-"+f*j);
//			else if (f==1) run("Make Substack...", "frames="+f*(j-1)+1+"-"+f*j);
			print(name + "_tile"+j+"...");
			saveAs("Tiff", output + File.separator + outputName + "_tile"+j);
			close();
		}
	}
	else saveAs("Tiff", output + File.separator + outputName);
	run("Close All");
}

