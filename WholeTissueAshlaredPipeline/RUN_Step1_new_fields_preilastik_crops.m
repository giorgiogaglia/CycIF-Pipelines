clear all
% This code is supposed to run the FIRST step of the t-CycIF Ashlar
% Analysis Pipeline
%
%t-CycIF ASHLAR ANALYSIS PIPELINE 
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
%
% In order to begin the analysis a set of parameters are needed to be
% defined by the user 
% 1) where the files are located and how the names are formatted
% 2) parameters to choose the size of the field and number of crops per
% field 

%%% OUTPUTS 
% FullStacks: AllTheRawData\OmeTifName\FullStacks\OmeTifName_Field_row_column.tif
% CroppedStacks: AllTheRawData\OmeTifName\CroppedData\OmeTifName_Field_row_column_#ofCrop.tif
% Coordinates: AllTheRawData\Coordinates_FullStacks\OmeTifName.mat 
% Coordinates Matrix: Coordinates.Field(row,column) = [keep field (0= no, 1=yes), x1, x2, y1, y2]
%
%%% FILE LOCATION AND FORMATTING
%
% The code will run on the "ome.tif" file produced by ASHLAR 
% The expected file structure is that a master folder will contain all of
% the "ome.tif" files from each slide
%
% eg. AllTheRawData\SlideName.ome.tif 

% 1) THE MASTER OME.TIF DATA FOLDER

filename.folders.main = 'W:\Analysis\Ashlared\';

% 2) SLIDE SPECIFIC PARAMETERS:
filename.folders.fols = {'AJ0160_P2', 'AJ0160_P3'}; %Name of ASHLARED image without '.ome.tif' ending 
filename.cycles = 8; % total # of cycles

% 3) USER DESIRED PARAMETERS 
filename.sizefield = 6000; %Size of field desired for a square 
filename.crops = 2; %# of cropped fields desired per field 
filename.dim = ['%02d']; %Delimeter 

% 4) OPTIONS  
options.max_prob = 65535;
options.cellsize = 20;
options.date = '20190223';

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
filename.suffix = '.tif';

outputfold = [filename.folders.main filename.folders.output]; 
addpath(filename.folders.main)
mkdir(outputfold) 

%% Splitting into size of field(t) desired 
for i = 1:length(filename.folders.fols)
    t= filename.sizefield ; %minus 1 to get the field to desired size
    tic
    Coordinates = {};
    fold = filename.folders.fols{i};
    disp(fold) 
    
    %Making folders
    outputfolder = [outputfold fold '\'];
    addpath(outputfold)
    mkdir([outputfolder filename.folders.fullstacks])
    mkdir([outputfold filename.folders.coordinates]) 
    mkdir([outputfolder filename.folders.cropfol])
    
    filenamesv = [outputfold filename.folders.coordinates fold];
    omifile = [ filename.folders.main filename.folders.fols{i} '.ome.tif'];
    
    DAPICycle0 = imread(omifile, 1);
    [y x] = size(DAPICycle0); %y= # of pixels in 1 column; x = # of pixels in 1 row
    column = 1;
    ix = 1;
    iy = 1; % Intital x and y coordinates
    toc
    %Loops while first coordinate is less than x-t for first row
    while ix < x-t
        fieldname = [fold '_Field_' num2str(1 , filename.dim) '_' num2str(column , filename.dim)];
        disp(fieldname)
        outfilename = [outputfolder filename.folders.fullstacks fieldname '.tif'];
        if ix == 1
            xcorr = t; %t to make each field (t+1)*(t+1)
        else
            xcorr = ix + t; %t to make each field (t+1)*(t+1)
        end 
        
        field = DAPICycle0(1:t, ix:xcorr);
        if prctile(field(:),97.5) < 100
            kf = 0; %kf= keep field?
        else
            kf= 1;
            imwrite(field, outfilename)
        end
        Coordinates.Field{1,column} = [kf, ix, xcorr, 1, t]; %Saving coordinates
        
        ix = xcorr + 1; %Save previous corrdinate to intialize next field +1 cause start at next matrix
        column = column + 1; %column number
    end
    
    %Create last field in first row
    fieldname = [fold '_Field_' num2str(1 , filename.dim) '_' num2str(column , filename.dim)];
    outfilename = [outputfolder filename.folders.fullstacks fieldname '.tif'];
    xcorr = x; %Last pixel in row
    field = DAPICycle0(1:t, ix:xcorr);
    if prctile(field(:),97.5) < 100
        kf = 0; %kf= keep field?
    else
        kf= 1;
        imwrite(field, outfilename)
    end
    Coordinates.Field{1,column} = [kf, ix, xcorr, 1, t]; %Saving coordinates
    toc
    
    %Loops while first coordinate is less than y-t for columns
    row = 1;
    iy = 1;
    while iy < y-t
        fieldname = [fold '_Field_' num2str(row , filename.dim) '_' num2str(1 , filename.dim)];
        disp(fieldname)
        outfilename = [outputfolder filename.folders.fullstacks fieldname '.tif'];
        if iy == 1
            ycorr = t;
        else 
            ycorr = iy + t; %t to make (t+1)*(t+1)
        end 
        field = DAPICycle0(iy:ycorr, 1:t);
        if prctile(field(:),97.5) < 100
             kf = 0; %kf= keep field?
        else
            kf= 1; 
            imwrite(field, outfilename)
        end
        Coordinates.Field{row,1} = [kf, 1,t, iy, ycorr]; %Saving coordinates
        iy = ycorr + 1; %Save previous corrdinate to intialize next field +1 cause start at next matrix
        row = row + 1; %column number
    end
    
    %Create last field in first Column 
    fieldname = [fold '_Field_' num2str(row , filename.dim) '_' num2str(1 , filename.dim)];
    outfilename = [outputfolder filename.folders.fullstacks fieldname '.tif'];
    ycorr = y; %Last pixel in row
    field = DAPICycle0(iy:ycorr, 1:t);
    if prctile(field(:),97.5) < 100
        kf = 0; %kf= keep field?
        else
            kf= 1;
            imwrite(field, outfilename)
    end 
    Coordinates.Field{row,1} = [kf, 1, t, iy, ycorr]; %Saving coordinates
    toc
    % Loop through column and row to create all fields
    for j1 = 2:column
        coordvect = Coordinates.Field{1,j1};
        xcoord = coordvect(2:3); %Intial, final
        for j2 = 2:row
            coordvect = Coordinates.Field{j2,1};
            ycoord = coordvect(4:5); %Intial, final
            field = DAPICycle0( ycoord(1):ycoord(2),xcoord(1):xcoord(2));
            fieldname = [fold '_Field_' num2str(j2 , filename.dim) '_' num2str(j1 , filename.dim)];
            disp(fieldname)
            outfilename = [outputfolder filename.folders.fullstacks fieldname '.tif'];
            if prctile(field(:),97.5) < 100
                kf = 0; %kf= keep field?
            else
                kf= 1;
                imwrite(field, outfilename)
            end
            Coordinates.Field{j2,j1} = [kf, xcoord, ycoord]; %Saving coordinates
        end
    end
    toc
    save(filenamesv, 'Coordinates', 'x', 'y', 't')
    %% Loop through channels to create stacks
    for ch = 2:(filename.cycles*4) 
        fprintf('Channel = %d', ch)
        Cycle = imread(omifile, ch);
        fprintf('Channel %d', ch)
        if size(Cycle)==size(DAPICycle0)
            for j1 = 1:column
                for j2 = 1:row
                    coordvect = Coordinates.Field{j2,j1};
                    ycoord = coordvect(4:5); %Intial, final
                    xcoord = coordvect(2:3); %Intial, final
                    kf = coordvect(1); %keep field?
                    field = Cycle( ycoord(1):ycoord(2), xcoord(1):xcoord(2));
                    fieldname = [fold '_Field_' num2str(j2 , filename.dim) '_' num2str(j1 , filename.dim)];
                    disp(fieldname)
                    outfilename = [outputfolder filename.folders.fullstacks fieldname '.tif'];
                    if kf == 1
                        imwrite(field, outfilename, 'WriteMode','append')
                    end
                end
            end
        else
            disp('Channels done')
            disp(ch)
        end
    end
    toc
    
    %% Loop through channels and fields to create CROPPED images
    
    for j1 = 1:column
        for j2 = 1:row
            for k1 = 1:filename.crops
                rowstart = randi(t-250);    %Cropped image coordiantes
                colstart = randi(t-250);
                for ch = 1:(filename.cycles*4) 
                    fprintf('Channel = %d\n', ch)
                    coordvect = Coordinates.Field{j2,j1};
                    kf = coordvect(1); %keep field?
                    fieldname = [fold '_Field_' num2str(j2 , filename.dim) '_' num2str(j1 , filename.dim)];
                    disp(fieldname)
                    outfilename = [outputfolder filename.folders.fullstacks fieldname '.tif'];
                    crfilename = [outputfolder filename.folders.cropfol fieldname '_' num2str(k1,'%02d') '.tif'];
                    if kf == 1
                        % Making Cropped image
                        orifield = imread(outfilename, ch);
                        try
                            crfield = orifield(rowstart:(rowstart+250), colstart:(colstart+250));
                        catch
                            continue
                        end 
                        if ch ==1
                            if prctile(crfield(:),97.5) < 100 %Skipping sections that don't have cells 
                                break 
                            else 
                                imwrite(crfield, crfilename);
                            end 
                        else
                            imwrite(crfield, crfilename, 'WriteMode','append');
                        end
                    end
                end
            end
        end
    end

end