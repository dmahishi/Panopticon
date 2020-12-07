// Pano.ijm
// 
// This ImageJ macro extracts objects from stacks of images. For further processing of the results use PanDataHandleList.m in Matlab.
//
// Author: Tilman Triphan, tilman.triphan@uni-leipzig.de
// License: GNU General Public License v3.0
// Copyright (C) 2020  Tilman Triphan
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

//Default settings
batchsize = 1800;
maxbatch = 200; // frames taken at beginning and end
blur_sigma = 3; // settings for Gaussian blur
// settings for particle analysis
min_diff_threshold = 60;
max_diff_threshold = 255;
min_size_threshold = 150;
max_size_threshold = 1200;
min_circ_threshold = 0.3;
path = ""; 	//defines folder to process and save data
fn_roi = ""; //rois to define innerBorder, outerBorder, Arena

//Process a batch of images
function processBatch(fn,start,batchnumber) {
	run("Image Sequence...", "open="+fn+" number="+d2s(batchsize,0)+" starting="+d2s(i*batchsize+1,0)+" sort");
	title = getTitle();
	sliceCount = nSlices;
	id = getImageID();

	// Create MIN-projection for each batch stack (shows fly movement as the fly is the darkest object)
	fn_min = title+"_"+batchnumber+"_MIN.tif";
	run("Z Project...", "projection=[Min Intensity]");
	saveAs("Tiff", path+fn_min);
	close();
	
	if (sliceCount < maxbatch) {
		selectWindow(title);
		close();
		return;
	}	
	
	selectWindow(title);
	getDimensions(m_width, m_height, m_channels, m_slices, m_frames);

	fn_max = title+"_MAX.tif";
	fn_max_save = title+"_"+batchnumber+"_MAX.tif";
	fn_max_stack = title+"_stack_MAX.tif";
	
	//Make MAX INTENSITY stack locally, using frames at start and end
	newImage(fn_max_stack, "8-bit white", m_width, m_height, maxbatch*2);
	idm = getImageID();
	for (i = 1; i <= maxbatch; i++) {
		selectImage(id);
		setSlice(i);
		run("Copy");
		selectImage(idm);
		setSlice(i);
		run("Paste");
	}	
	for (i = maxbatch+1; i <= maxbatch*2; i++) {
		selectImage(id);
		setSlice(m_slices-maxbatch*2+i);
		run("Copy");
		selectImage(idm);
		setSlice(i);
		run("Paste");
	}
	
	selectImage(idm);
	run("Z Project...", "projection=[Max Intensity]");	
	saveAs("Tiff", path+fn_max_save);
	rename(fn_max);
	selectImage(idm);
	close();
	
	selectWindow(fn_max);
	run("Select All");
	run("Copy");
	newImage(fn_max_stack, "8-bit white", m_width, m_height, sliceCount);
	selectWindow(fn_max_stack);
	for (i = 1; i <= sliceCount; i++) {
		setSlice(i);
		run("Paste");
	}

	imageCalculator("Subtract stack", fn_max_stack, title);
	selectWindow(fn_max_stack);
	//Enhance contrast to find flies on the dark side
	run("Enhance Contrast", "saturated=0.35");
	run("Apply LUT", "stack");
	run("Gaussian Blur...", "sigma="+blur_sigma+" stack");	
	setThreshold(min_diff_threshold, max_diff_threshold);

	//Use a roi to define the area in which to track objects
	roiManager("reset");
	roiManager("open", fn_roi);
	roiManager("select", 2);
	
	run("Analyze Particles...", "size="+min_size_threshold+"-"+max_size_threshold+" circularity="+min_circ_threshold+"-1.00 display clear add stack");

	selectWindow("Results");
	fn_results = title+"_"+batchnumber+"_results.xls";
	saveAs("Results", path+fn_results);	
	fn_roi = title+"_"+batchnumber+"_roiset.zip";
	roiManager("Save", path+fn_roi);

	selectWindow(title);
	close();
	selectWindow(fn_max_stack);
	close();
}

//Process a folder containing several batches of images
function processFolder(input) {
	starttime = getTime();
	print("Processing folder: " + input);
	list = getFileList(input);
	print("frames: " + list.length);
	filecount = list.length;
	batchcount = floor(filecount/batchsize)+1;
	print ("batches: " + batchcount);

	//create and process batches
	for (i = 0; i < batchcount; i++) {
		print("Current batch: " + i);
		processBatch(input,i*batchsize+1,IJ.pad(d2s(i,0),2));
		print("batch " + i + " is done!");
	}
	endtime = getTime();
	print((endtime - starttime)/60000 + " minutes");
}

//Read textfile containing paths to folders and start analysis
run("Set Measurements...", "area centroid shape stack redirect=None decimal=3");
setBatchMode(true);
path = getDirectory("Please select folder to analyze");
fn_roi = path+"RoiSet.zip";
processFolder(path);
print("all done!");
setBatchMode(false);
