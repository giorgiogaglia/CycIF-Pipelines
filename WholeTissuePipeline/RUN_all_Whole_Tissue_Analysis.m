%Runs the whole tissue analysis 
clear all
%Please change all these inputs 
infolderloc = 'D:\Myeloma_HSF1_2018-07-27\RawData';
outfolderloc = 'D:\Myeloma_HSF1_2018-07-27\Test\';
basicfolderloc = 'D:\Myeloma_HSF1_2018-07-27\Analysis\AnalysisWBaSiCcorrection\BasicCorrectionImages\Myeloma_BS16_14824\Cycle_';
folders = {'Myeloma_BS16_49921_P1'};
maxfield = [16];
rowstiles = 4;
colstiles = 4;
cycles = [2 3 4 5];
maxcycle = 4; 
HSF1cycle = 1; %Cycle is from 1... 
znum = 4; %number of zstacks 
DAPIslice = [1, 6, 10, 14];
cellsize = 23.5; 
dates = '\2018-07-27_'; 
thr = 2; %Foci threshold 
montfols = {'n/a'};
montfolder ='X:\Carmen\Myeloma_Montages\ROIs_Montages\';
marker_name = {'not used', 'HSF1','Hsp70','Hsp70_2','Hsp90','Hsp90_2','CD138','CytoC'}; %HSF1 should always be 1st marker  
antibody_rounds = [3, 6, 14, 7, 15,10, 8]; %Numbering is 1=DAPI cycle 1; 2= 488 cycle 1; 3=555 cycle 1; 4=647 cycle 1 
marker_cycle = [1 1 2 3 2 3 3 2];
HSF1round = 3; %Numbering is 1=DAPI cycle 1; 2= 488 cycle 1; 3=555 cycle 1; 4=647 cycle 1 
canceround = 10; %Cancer marker 

%Don't need to change 
tilesize = [2048 2048];
prefix1 = {'A - 1(fld '};
wavelengths = {'wv Blue - FITC','wv Green - dsRed','wv Red - Cy5'};
basefolder = outfolderloc;

%% Running analysis 
x= input('Are there Z-stacks? ''Y''/''N'' ');
antibody_type = [];
for i1=1:length(antibody_rounds)
            str= sprintf('Is %s a cytoplasm or nucleus marker? (''C''/''N'')', marker_name{i1+1}); 
            str = input(str); 
            antibody_type = [antibody_type, str];
end 
if isequal(x, 'Y')
    [indx,tf] = listdlg('ListString',wavelengths, 'PromptString', 'Which wavelength is the Z-stacks in?: ', 'SelectionMode', 'single');
    indx1 = input('Which cycle is it in? (1 being the first cycle of cycles): ');
    disp('Run 1: Making Stacks')
    RUN_1_makecorestacks_z(infolderloc, outfolderloc,basicfolderloc, folders, maxfield, cycles, tilesize, prefix1, wavelengths, indx, indx1) 
    disp('Run 2: Segmenting Cells')
    RUN_2_segmentstacks(basefolder, folders, maxfield, prefix1,DAPIslice, cellsize) 
    disp('Run 3: Foci Segmentation & making Measurements') 
    RUN_3_CycIF_measurements_Z(basefolder,folders, maxfield, maxcycle, dates, DAPIslice, HSF1cycle, znum, thr)   
    disp('Removing ROIs')
    remove_ROI(folders, rowstiles, colstiles, maxfield,montfolder, montfols) %Montages should have same name as folders
    disp('Run 4: Creating Graphs') 
    RUN_4_CellStateCompare_Path(basefolder, folders, dates,marker_name, antibody_rounds, marker_cycle, HSF1round, canceround, antibody_type) 
else 
    disp('Run 1: Making Stacks')
    RUN_1_makecorestacks(infolderloc, outfolderloc,basicfolderloc, folders, maxfield, cycles, tilesize, prefix1, wavelengths) 
    disp('Run 2: Segmenting Cells')
    RUN_2_segmentstacks(basefolder, folders, maxfield, prefix1,DAPIslice, cellsize) 
    disp('Run 3: Foci Segmentation & making Measurements') 
    RUN_3_CycIF_measurements(basefolder,folders, maxfield, maxcycle, dates, DAPIslice)
    disp('Removing ROIs')
    remove_ROI(folders, rowstiles, colstiles, maxfield,montfolder, montfols) %Montages should have same name as folders
    disp('Run 4: Creating Graphs') 
    RUN_4_CellStateCompare_Path(basefolder, folders, dates,marker_name, antibody_rounds, marker_cycle, HSF1round, canceround, antibody_type) 
end 




