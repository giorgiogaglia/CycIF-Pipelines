% clear all
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

master_folderIN = 'Y:\IN Cell Analyzer 6000\Carmen\Test';

% 2) SLIDE AND PLATE SPECIFIC PARAMETERS: Comment out the other section

% 2A) WHOLE TISSUE SLIDE: THE RAW DATA SUBFOLDERS FOR EACH SLIDES AND HOW MANY TILES THEY HAVE

slides_folders = {'A_5B_dHSF2','A-1C_WT','B_5B_dHSF2', 'B-1C_WT'};   % names of the folders for each slides
slides_rowtils = [5,5,5,5];                         % number of rows in each slide
slides_coltils = [2,2,2,2];                        % number of cols in each slide
maxfields = slides_rowtils.*slides_coltils;       % number of fields imaged
%DO NOT CHANGE
    prefix1 = {'A - 1(fld '};    % prefix for the name of each file
    max_rows = 1;                       % total number of plate rows; ex. 2 for rows B & C
    well_nums = 1;                   % vector of number for wells; ex. 2 to 10 for columns 2 through 10 
    well_lets = {'A'};
    DAPIwv = 'wv UV - DAPI';    
    wavelengths = {'wv Blue - FITC','wv Green - dsRed','wv Red - Cy5'};

% 2B) IF YOU HAVE A TMA 
%CARMEN ADD THIS IN! 

% % 2B) IF YOU HAVE A PLATE
% 
% slides_folders = {'Slide_1','Slide_2','Slide_3'};   % names of the folders for each plate 
% maxfields = 9;                        % number of fields imaged 
% prefix1 = {'B - 02(fld ','B - 03(fld ','B - 04(fld ','B - 05(fld ','B - 06(fld ','B - 07(fld ','B - 08(fld ','B - 09(fld ','B - 10(fld '; ...     %Prefix for the name of each file, MUST START A NEW ROW FOR A NEW ROW 
%     'C - 02(fld ','C - 03(fld ','C - 04(fld ','C - 05(fld ','C - 06(fld ','C - 07(fld ','C - 08(fld ','C - 09(fld ','C - 10(fld '};
% well_nums = 2:10;                   % Vector of number for wells; ex. 2 to 10 for columns 2 through 10 
% well_lets = {'B', 'C'};             % Letters that correspond to the name of the row 
% DAPIwv = 'wv UV - DAPI z 07';    % Change if suffix is different  
% wavelengths = {'wv Blue - FITC z 07','wv Green - dsRed z 07','wv Red - Cy5 z 07'};
% % DO NOT CHANGE
%     slides_rowtils = 0;                        
%     slides_coltils = 0;     
%     max_rows = 2=length(well_lets);                       %Total number of plate rows; ex. 2 for rows B & C

% 3) THE RAW DATA SUBFOLDER FOR EACH CYCLE WITHIN EACH SLIDE

cycles_folders = {'Cycle_1','Cycle_2','Cycle_3'};

% 4) THE DATA FOLDER FOR BASIC CORRECTION IMAGES
%Use imagej_basic_ashlar.py to create Basic Correction Images
basic_folderloc = 'Y:\IN Cell Analyzer 6000\Carmen\Test\BasicCorrection\Cycle_';

% 5) MONTAGES FOR ROI: If you do not have montages then leave both
% variables as empty
montfolder = '';   % base folder for montages
montfols = {''};                     % names of the ROI for each montage 


% 6) USER INPUT PARAMTERS
cycles = [1 2 3]; %Cycle number is based on number in cycles_folders 
dates = '\2018-10-19'; 
% marker_name = {'HSF1','Hsp70','Hsp90','CD138','CytoC'}; %List of markers
% antibody_type = {'N', 'C', 'C','C'}; %Type of marker for each antibody; N= Nucleus; C= Cytoplasm 
% antibody_rounds = [3, 6, 14, 7, 15]; %Numbering is 1=DAPI cycle 1; 2= 488 cycle 1; 3=555 cycle 1; 4=647 cycle 1 
% marker_cycle = [1 1 2 3 4]; %Cycle the antibody is in 1 being the first cycle in cycles variable 
% canceround = 10; %Used to filter the data based on the cancer marker
% sol_thresh = 0.9;    %Threshold to filter data based on cancer marker 

%7) USER INPUT PARAMETERS FOR SEGMENTATION
options.cellsize = 31;     %Threshold for cell segmentation (must be odd) 
options.sigma = 1; 
options.findminvec_method = 3;
options.findminvec_quantile = 0.3;
options.findminvec_folddiff = 3;
options.writeflag = 0;
options.fixbrokencellsflag = 0;
sigma = 1;



% 7) IF THERE ARE Z-STACKS
z_cycle = 0;                    % Number of cycle the z-stack is in, default is 0 for no z-stacks
z_wave = 'wv Green - dsRed';  % Wavelength of the z-stack, possible wavelengths: {'wv Blue - FITC','wv Green - dsRed','wv Red - Cy5'}
z_num = 0;                      % Number of z-stacks 

% 8) IF YOU NEED FOCI SEGMENTATION 
HSF1cycle = 1; %Cycle is from 1... 
HSF1round = 3; %Numbering is 1=DAPI cycle 1; 2= 488 cycle 1; 3=555 cycle 1; 4=647 cycle 1 
thr = 2; %Foci threshold 

%% Running analysis 
basefolder = [master_folderIN '\Analysis\']; 
basicfolderloc = [basic_folderloc '\Cycle_'];
maxcycle = length(cycles); 
tilesize = [2048 2048];

DAPIslice = []; 
CYCLEslice = {};
x=0; 
for i2 = 1:maxcycle
    if i2 == 1
        x = 1;
    elseif cycles(i2-1) == z_cycle
        x = y + 3 + z_num; %3 is for the other 2 channels
    else
        x = y + 4; 
    end
    DAPIslice = [DAPIslice, x]; 
    y = x; %Previous DAPI slice 
    channels = {}; 
    if cycles(i2) == z_cycle
        channels{1} = 'wv UV - DAPI'; 
        for i3 = 1:3
            if isequal(wavelengths{i3}, z_wave) 
                for i4 = 1:z_num 
                    channels = [channels, z_wave];
                end
            else 
                channels = [channels, wavelengths{i3}];
            end 
        end     
    else 
        channels = {'wv UV - DAPI', 'wv Blue - FITC','wv Green - dsRed','wv Red - Cy5'};
    end 
    CYCLEslice{i2} = channels; 
end 


%%
disp('Run 1: Making Stacks')
RUN_1_makecorestacks(master_folderIN, basefolder,basicfolderloc, slides_folders, maxfields, cycles, tilesize, prefix1, wavelengths, z_wave, z_cycle, max_rows, well_nums , well_lets, DAPIwv)
disp('Run 2: Segmenting Cells')
RUN_2_segmentstacks(basefolder, slides_folders, maxfields,DAPIslice, options, max_rows, well_nums, well_lets)
disp('Run 3: Foci Segmentation & making Measurements')
RUN_3_CycIF_measurements(basefolder,slides_folders, maxfields, maxcycle, dates, DAPIslice, z_cycle, z_num, thr, CYCLEslice, prefix1,max_rows, well_nums, well_lets, cycles)
if ~isempty(montfols)
    disp('Removing ROIs')
    remove_ROI(slides_folders, slides_rowtils, slides_coltils, maxfields,montfolder, montfols) 
end

% % disp('Run 4: Creating Graphs')
% RUN_4_CellStateCompare_Path(basefolder, slides_folders, dates,marker_name, antibody_rounds, marker_cycle, HSF1round, canceround, antibody_type, max_rows, well_nums)

%% OVERVIEW AND EXPLANATION

% CODE REQUIRED TO RUN ANALYSIS:
% 1) RUN_1_makecorestacks
% 2) RUN_2_segmentstacks
% 3) RUN_3_CycIF_measurements
% 4) CycIF_Segmentation_normfit_opencellsizeover5_medianpercell_v2
% 5) findmininvect
% 6) saveastiff

% ADDITIONAL CODES
% - imagej_basic_ashlar.py : Required if planning to use BaSic Correction
% - ReadImageJROI : Required if planning to remove portions of tissue from analysis 
% - remove_ROI : Required if planning to remove portions of tissue from analysis 
% - RUN_4_CellStateCompare_Path : Required if planning to create graphs
% from foci analysi 

% EXPLANATION:
% 1) RUN_1_makecorestacks: Creates two stacks; skipping fields that have no
% cells in first cycle 
        % A) Only DAPI channel
        % B) All Channels for all cycles  
% 2) RUN_2_segmentstacks: Creates one stack- all segmented DAPI from all
% cycles
% 3) RUN_3_CycIF_measurements: Creates one mat file per slide of single
% cell measurments 

