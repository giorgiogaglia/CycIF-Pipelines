clear all
% This Code is supposed to run all the steps of the t-CycIF Analysis
% Pipeline in MatLab

% In order to begin the analysis a set of parameters are needed to be
% defined by the user 
% 1) where the files are located and how the names are formatted
% 2) parameters to choose the type of segmentation and details like size of
% cells

%%% FILE LOCATION AND FORMATTING
%
% The code will run on the raw .tiff files as most platforms will export or
% allow to out the .tiff files
% The expected file structure is that a master folder will contain all of
% the raw image files from the cycles with each slide having a separate 
% subfolder, followed by cycle having a separate subfolder with all the
% .tiff files in it
%
% eg. AllTheRawData\SlideXX\CycleXX\Field_1_Channel_1.tif

% 1) THE MASTER RAW DATA FOLDER

master_folderIN = 'D:\Myeloma_HSF1_2018-07-27\RawData';

% 2) THE RAW DATA SUBFOLDERS FOR EACH SLIDES AND HOW MANY TILES THEY HAVE

slides_folders = {'Slide_1','Slide_2','Slide_3'};   % names of the folders for each slides
slides_rowtils = [4, 5, 11];                        % number of rows in each slide
slides_coltils = [4, 4,  2];                        % number of cols in each slide
slides_maxflds = slides_rowtils.*slides_coltils;    % number of fields imaged

% 3) THE RAW DATA SUBFOLDER FOR EACH CYCLE WITHIN EACH SLIDE

cycles_folders = {'Cycle_1','Cycle_2','Cycle_3','Cycle_4'};

% 4) THE DATA FOLDER FOR BASIC CORRECTION IMAGES
%Use imagej_basic_ashlar.py to create Basic Correction Images
basic_folderloc = 'D:\Myeloma_HSF1_2018-07-27\Analysis\AnalysisWBaSiCcorrection\BasicCorrectionImages\Myeloma_14824\Cycle_';

% 5) MONTAGES FOR ROI
montfolder = 'X:\Carmen\Myeloma_Montages\ROIs_Montages\';   % base folder for montages
montfols = {'ROI_1', 'ROI_2', 'ROI_3'};                     % names of the ROI for each montage, include extention (.zip or .roi) 


% 6) USER INPUT PARAMTERS
cycles = [2 3 4 5]; %Cycle number is based on number in cycles_folders 
dates = '\2018-07-27_'; %Date to differentiate between repeat analysis 
marker_name = {'HSF1','Hsp70','Hsp90','CD138','CytoC'}; %List of markers
antibody_type = {'N', 'C', 'C','C'}; %Type of marker for each antibody; N= Nucleus; C= Cytoplasm 
antibody_rounds = [3, 6, 14, 7, 15]; %Numbering is 1=DAPI cycle 1; 2= 488 cycle 1; 3=555 cycle 1; 4=647 cycle 1 
marker_cycle = [1 1 2 3 4]; %Cycle the antibody is in 1 being the first cycle in cycles variable 
canceround = 10; %Used to filter the data based on the cancer marker

% THRESHOLD VALUES
cellsize = 23.5;     %Threshold for cell segmentation 
sol_thresh = 0.9;    %Threshold to filter data based on cancer marker 

%% Running analysis
[base_folder,filename, ext] = fileparts(master_folderIN); 
basefolder = [base_folder '\Analysis\']; 
basicfolderloc = [basic_folderloc '\Cycle_'];
maxcycle = length(cycles); 
marker_name = ['not used', marker_name]; 
tilesize = [2048 2048];
prefix1 = {'A - 1(fld '};
wavelengths = {'wv Blue - FITC','wv Green - dsRed','wv Red - Cy5'};
DAPIslice = []; 
x=0; 
for i2 = 1:maxcycle
    if i2 == 1
        x = 1;
    else
        x = y + 4; 
    end
    DAPIslice = [DAPIslice, x]; 
    y = x; %Previous DAPI slice 
end 
    
disp('Run 1: Making Stacks')
RUN_1_makecorestacks_Validation(infolderloc, outfolderloc,basicfolderloc, folders, maxfield, cycles, tilesize, prefix1, wavelengths) 
disp('Run 2: Segmenting Cells')
RUN_2_segmentstacks_Validation(basefolder, slides_folders, slides_maxflds, prefix1,DAPIslice, cellsize)
disp('Run 3: Making Measurements')
RUN_3_CycIF_measurements_Validation(basefolder,folders, maxfield, maxcycle, dates, DAPIslice)    
disp('Removing ROIs')
remove_ROI(slides_folders, slides_rowtils, slides_coltils, slides_maxflds,montfolder, montfols) 
disp('Run 4: Creating Graphs')
RUN_4_CellStateCompare_Path_Validation(basefolder, folders, dates,marker_name, antibody_rounds, marker_cycle, HSF1round, canceround, antibody_type, sol_thresh) 

