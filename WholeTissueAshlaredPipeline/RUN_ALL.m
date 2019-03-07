clear all
% This code is supposed to run the THIRD AND FOURTH step of the t-CycIF Ashlar
% Analysis Pipeline
%
% t-CycIF ASHLAR ANALYSIS PIPELINE 
% Step 0: Stitch and register the fields into an "ome.tif" by using Ashlar 
% Step 1:(RUN_Step1_new_fields_preilastik_crops.m) Cut the "ome.tif" into 
% fields of a size specified by the user and cut crops out of each field 
% for Ilastik training. Create full stacks of all cycles and channels. 
% Omits fields that have no cells.   
% Step 2: Create segmentation probabilities using Ilastik
% (RUN_ALL_segon.m)Runs segmentation and measurements scripts (Step 3 & 4) 
% Step 3:(RUN_Step3_segmentfromilastik.m) Segment based on segmentation
% probabilities produce by Ilastik
% Step 4: (RUN_Step4_CycIF_measurments_ilastik.m) Makes measurements of
% signal and foci segmentation 
% Step 5: (RUN_Step5_ROI.m) Finds pixel location of ROIs and indexes of
% centroids within the ROI 
% 
% In order to begin the analysis a set of parameters are needed to be
% defined by the user 
% 1) where the files are located and how the names are formatted
% 2) parameters to choose the size of the field and number of crops per
% field 

%%% OUTPUTS FILES AND LOCATIONS 
% Step 1:
% FullStacks: AllTheRawData\OmeTifName\FullStacks\OmeTifName_Field_row_column.tif
% CroppedStacks: AllTheRawData\OmeTifName\CroppedData\OmeTifName_Field_row_column_#ofCrop.tif
% Coordinates: AllTheRawData\Coordinates_FullStacks\OmeTifName.mat 
% y= # of pixels in 1 column; x = # of pixels in 1 row ; t = size of field desired 
% Coordinates Matrix: Coordinates.Field(row,column) = [keep field (0= no, 1=yes), x1, x2, y1, y2]
% Montage: AllTheRawData\MontageforROI\OmeTifName_montage.tif 
% Step 3:
% Segmented Images: AllTheRawData\ANALYSIS\OmeTifName\Ilastik_Segmentation\OmeTifName_Field_row_column_Seg.tif
% Check Segmented Images: AllTheRawData\ANALYSIS\OmeTifName\Ilastik_Segmentation\OmeTifName_Field_row_column_checkseg.tif
% Step 4: 
% Nucleus & Cytoplasm SegmentedImages: AllTheRawData\ANALYSIS\OmeTifName\Ilastik_Segmentation\OmeTifName_Field_row_column_NucCytSeg.tif
% Foci Segmented Images: AllTheRawData\ANALYSIS\OmeTifName\Foci_Segmentation\OmeTifName_Field_row_column__HSF1_FociSeg.tif
% Foci Check Segmented Images: AllTheRawData\ANALYSIS\OmeTifName\Foci_Segmentation\OmeTifName_Field_row_column__HSF1_FociSeg_check.tif
% Measurments file: AllTheRawData\ANALYSIS\Analysis_Results\Results_data.mat 
% Step 5:
% ROI_pixels: AllTheRawData\ROIpixels\OmeTifName\ROI_pixels.mat: matrix of pixel locations (x,y) of ROIs, vector of row indexes of Centroids that are within the ROI  


%%% INPUT FILE LOCATION AND FORMATTING
%
% The code will run on the Ilastik Probabilities file and FullStacks file 
% The expected file structure is that a master folder will contain all of
% the "ome.tif" data and within that folder, a folder called "ANALYSIS\"
% will contain all the Ilastik Probabilities and FullStacks data 
%
%%% IMPORTANT: Ilastik Output Filename Format for Export Options
% {dataset_dir}/Ilastik_Probabilities/{nickname}_Probabilities.tif

% eg. AllTheRawData\ANALYSIS\OmeTifName\FullStacks\Ilastik_Probabilities\SlideName_Field_row_column_Probabilities.tif 
% eg. AllTheRawData\ANALYSIS\OmeTifName\FullStacks\SlideName_Field_row_column.tif 

% ROI:
% Must be drawn on ImageJ on the the montage created in RUN_Step1_new_fields_preilastik_crops
% Location: Montage: AllTheRawData\MontageforROI\OmeTifName_montage.tif 
% ROIs must be save in "AllTheRawData\MontageforROI\" and named with "OmeTifName" 

% 1) THE MASTER OME.TIF DATA FOLDER

filename.folders.main = 'W:\Analysis\Ashlared\';

% 2) SLIDE SPECIFIC PARAMETERS:
filename.folders.fols = {'AJ0160_P2', 'AJ0160_P3'}; %Name of ASHLARED image without '.ome.tif' ending 
filename.cycles = 8; % total # of cycles

% 3) USER DESIRED PARAMETERS 
filename.sizefield = 6000; %Size of field desired for a square 
filename.crops = 2; %# of cropped fields desired per field 


% 4) OPTIONS  
% Step 3: SEGMENTATION OPTIONS 
options.nuc = 1;    %Channel nucleus Ilastik probability is in
options.cyt = 3; %Channel cytoplasm Ilastik probability is in
options.backgr = 2; %Channel background Ilastik probability is in 
options.cellsize = 20;
options.max_prob = 65535; %Maximum Ilastik probability, usually does not change

% Step 4: MEASUREMENTS AND FOCI OPTIONS 
options.date = '20190223';
options.focithr = 2; %Foci threshold 
options.hsf1round = 0; %HSF1 cycle 

% 5) OUTPUT PARAMETERS: DO NOT EDIT 
filename.folders.output = 'ANALYSIS\'; 
filename.folders.fullstacks = 'FullStacks\';
filename.folders.cropfol = 'CroppedData\';
filename.folders.coordinates = 'Coordinates_FullStacks\';
filename.folders.ilastikprob= 'Ilastik_Probabilities\';
filename.folders.ilastikseg = 'Ilastik_Segmentation\';
filename.ilastiksuffix = '_Probabilities.tif';
filename.folders.results = 'Analysis_Results\';
filename.folders.fociseg = 'Foci_Segmentation\'; 
filename.folders.ROI = 'ROIpixels\'; 
filename.suffix = '.tif';
filename.folders.montage = 'MontageforROI\'; 
filename.dim = ['%02d']; %Delimeter 

%% MAKE SURE TO COMMENT OUT THE OTHER STEPS THAT YOU AREN'T AT YET!!
% Runs Step 1: (Creating FullStacks, CroppedStacks, Montages) 
RUN_Step1_new_fields_preilastik_crops(filename,options) 

%% Runs Step 3 and 4: (Segmentation and Measurments) 
RUN_Step3_segmentfromilastik(filename, options) 
RUN_Step4_CycIF_measurements_ilastik(filename, options)  

%% Run Step 5: (ROI Selection)
RUN_Step5_ROI(filename, options)