% This program saves the smallest montage made by ASHLAR for ROI drawing on
% ImageJ 
filename.folders.main = 'Z:\sorger\data\IN_Cell_Analyzer_6000\Giorgio\2019-01-03_Mouse_Lung_Tumors_Round2\Analysis\Ashlared\';
filename.folders.fols = {'AJ0160_P3'};
filename.folders.output = 'ANALYSIS\'; 
filename.cycles = 8;

for i = 1:length(filename.folders.fols)
    tic
    fold = filename.folders.fols{i};
    disp(fold) 
    
    %Making folders
    outputfolder = [filename.folders.main filename.folders.output 'MontageforROI\'];
    addpath(filename.folders.main)
    mkdir(outputfolder)
  
    filenamesv = [outputfolder fold '_montage.tif'];
    omifile = [ filename.folders.main filename.folders.fols{i} '.ome.tif'];
    ominfo = imfinfo(omifile); 
    slices = length(ominfo); %Total number of slices in the omi.tif file 
    
    channels = filename.cycles*4; 
    
    smallslice = slices - (channels-1); 
    
    DAPImontage = imread(omifile, smallslice); 
    
    imwrite(DAPImontage, filenamesv) 
    toc 
end 

