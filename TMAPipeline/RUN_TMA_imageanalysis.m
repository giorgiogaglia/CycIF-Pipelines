%This program runs the image analysis for TMAs

% Step 1: Making the core montage using 1 10x simage- finds each core and outputs a RealCore tif file for each core

% Step 2: Makes the CoreStacks for each core: ouputs DAPI_Core (DAPIStack) & Core (FullStack)

% Step 3: Segments each Core: outputs TrackedField stack of segmented DAPI images 

%%% INPUT PARAMETERS:

%% A. Parameters that depend on the microscopy setup 
microscope.tilesize = [2048 2048];

%% B. Parameters that define where the files to be analysed are located 
% It is assumed that the files are organized as follows
%
% "folderloc" = each slide/TMA group has a master folder where the DATA and the ANALYSIS will be contained
%
% -- the master folder will have two subfolders called DATA and ANALYSIS
%
% --- the DATA subfolder will have subfolder called [cycprefix][cycle number]
% ie either "Cycle_1" or "Round1" or "iamterribleatnamingfolders_HAHA_1"
% the REAL CORE images should be called with the same prefix and cycle number == 0 
%
% ---- the TMA cores will have a "prefix1" (usally a letter representing 
% the row) and "prefix2" (usually a number representing the column) and a
% "suffix" which has to be common to all cores.

filename.folderloc  = 'D:\Giorgio\TMA_BRC15010_3_6Jun18\';
filename.datafolder = [folderloc 'DATA' filesep];      % no need to change this
filename.analfolder = [folderloc 'ANALYSIS' filesep];  % no need to change this
filename.cycprefix  = 'Cycle_';



filename.prefix1  = {'A','B','C','D','E','F','G','H','I','J'};  % cols of cores
filename.midfix1  = ' - ';
filename.prefix2  = linspace(1,15,15); % rows of cores
filename.midfix2  = '(fld ';
filename.suffix   = '.tif';
filename.channels = {' wv UV - DAPI)',' wv Blue - FITC)',' wv Green - dsRed)', ' wv Red - Cy5)'};
filename.cycles   = [1 2 3 4 5 6 7 8];  % cycles to analyse
filename.dim      = '%02d'; %deliminator - if the cores are called A - 1, set this to '%01d'
filename.rowspercore = 2; % images per each core
filename.colspercore = 2; % images per each core
filename.maxtile = filename.rowspercore*filename.colspercore;

filename.DAPIslices = (filename.cycles-1)*4+1;


%% C. Analysis parameters - need to be organized in an "options" structure

%For 10X Montage
% rows = 1;
% cols = 1;
% maxtile = rows*cols;

options.cut = 100;
options.background = 1000;
options.radius10X = 700;
options.MagDiff10x = 2;
options.offsetthr = microscope.tilesize/4;


% if this is only an update modify this parameters
options.num_old_DAPI = 0;

options.writeflag = 0;
options.fixbrokencellsflag = 0;
options.sigma = 2;
options.cellsize = 40;
options.smallmorphs = 2;
options.findminvect_method = 3;
options.findminvect_quantile = 0.5;
options.findminvect_fold = 8;
options.cellfracoverlap = 0.4;


%% Running image Analysis

% first we make the "RealCores" from the 10X imaging we did on the first
% round

RUN_Step1_makecoremontage_single10X(microscope,filename,options)

RUN_Step2_makecorestacks_update(microscope,filename,options)

RUN_Step3_segmentstacks(microscope,filename,options) 

RUN_Step4_CycIF_measurements(folderloc, prefix2, prefix1, maxcycle, dates, DAPIslice)











