% Block 1: 
% Tracked data and arena info (center, areas, borders) in .csv format is 
% loaded into MATLAB, followed by assigning and splitting into individual 
% sectors\animal(s). After adding metadata info (.xlsx), processed data is 
% saved, and XY positions are plotted.
% 
% Block 2:
% Quality check includes:
% Temporal distribution and total percentage of NaNs
% Walking traces (with time on z axis) for visual inspection
% Sector-wise traces of longest movements, to identify tracing artifacts
% 
% Block 3:
% Loading processed data either independently or combined
%
% Each block 1-3 can be run independently as well.

%% Loading and processing of tracked data
% here, FIJI-generated .csv files are loaded and pre-processed.
% 
% path to folder containing FIJI files of single experiments
% !! metaData file name must match expt folder name !!
datapath = 'D:\Panopticon\DatasetA';
% batch size
batchsize = 1800;
% data structure / list for all files is created
pano = PanDataHandleList.CreateFromFolder(datapath);
% sectors (1-8) are assigned
pano.assignSector
% original data is split into different sectors 
pFixed = pano.splitAndRepair('Sector',batchsize);
% NaNs are handled more efficiently
pFixed.fixGaps
% metadata is added
pFixed.addDataPath(datapath)
pFixed.addImagePath;
[~,b,~] = fileparts(datapath);
xlspath = [datapath,'\',b,'.xlsx'];
pFixed.addMetaData(xlspath);
pFixed.mergeInfoFields({'Genotype','Sex','Substrate'},'Name');
% data is saved
pFixed.saveToFile([datapath, '.mat']);
% XY positions are plotted
pFixed.plotXY

%% Data quality check
% percentage of NaNs for every individual fly data
pFixed.plotNans
% individual fly tracks in 3D (time on z axis)
pFixed.plotWalkingTrace
% images with traces of longest movements (default: 10 frames)
pFixed.showLongDistanceMoves(nFrames)

%% Load processed data
% processed data can be loaded in multiple ways
% single .mat file
inputfile = 'D:\Panopticon\DatasetA.mat';
DatasetA = PanDataHandleList.LoadFromFile(inputfile);
% text file with several .mat files combined
% text file format: 
% D:\Panopticon\DatasetB1.mat
% D:\Panopticon\DatasetB2.mat 
inputtextfile = 'D:\Panopticon\DatasetB.txt'; 
DatasetB = PanDataHandleList.LoadFromTextFile(inputtextfile);

%% Analysis and plotting
% example for plotting total activity (%) with default parameters.
% !! for customization, refer to PanoPlot.m !!
PanoPlot({DatasetA, DatasetB},{'defaultPlot','Activity'}) 
%
%% Additional plotting functions
% optional, uncomment for use

% fly location probability - use RESSCALE variable to change resolution
% pFixed.plotResidency

% images with movement traces (Frame subset range)
% pFixed.plotTracksOnImage(1:20)

% temporal distribution of patch area visits (areaID equal or greater 2)
% pFixed.plotOnArea

% fly positions in 3D (time on z axis)
% pFixed.plotPositions
