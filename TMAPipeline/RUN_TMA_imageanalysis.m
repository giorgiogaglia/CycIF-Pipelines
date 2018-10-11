%This program runs the image analysis for TMAs
%Step1: Making the core montage using 1 10x simage- finds each core and
%outputs a RealCore tif file for each core
%Step 2: Makes the CoreStacks for each core: ouputs DAPI_Core (DAPIStack) &
%Core (FullStack)
%Step 3: Segments each Core: outputs TrackedField stack of segmented DAPI
%images 

% parameters that should not change
tilesize = [2048 2048];
wavelengths = {'wv Blue - FITC','wv Green - dsRed','wv Red - Cy5'};
thr = max(tilesize);

%Change All parameters below 
folderloc = 'D:\Giorgio\TMA_BRC15010_3_6Jun18\';
prefix1 = {'A','B','C','D','E','F','G','H','I','J'};  % cols of cores
prefix2 = linspace(1,15,15); % rows of cores
cycles = [1 2 3 4 5 6 7 8];  % cycles to analyse
dim = ['%02d']; %deliminator 
dates = '\2018-07-27_'; 

%For 10X Montage
% rows = 1;
% cols = 1;
% maxtile = rows*cols;
cut = 100;
background = 1000;
radius10X = 700;

% images per each core
rows2 = 2;
cols2 = 2;
maxtile2 = rows2*cols2;

% if this is only an update modify this parameters
num_old_DAPI = 0;


%% Running image Analysis

RUN_Step1_makecoremontage_single10X(folderloc, prefix1, prefix2, tilesize, dim, background, radius10X)
RUN_Step2_makecorestacks_update(folderloc,prefix1, prefix2, wavelengths, thr, cycles, dim, num_old_DAPI, rows2, cols2)
RUN_Step3_segmentstacks(folderloc, prefix1,prefix2, DAPIslice, cellsize) 
RUN_Step4_CycIF_measurements(folderloc, prefix2, prefix1, maxcycle, dates, DAPIslice)
