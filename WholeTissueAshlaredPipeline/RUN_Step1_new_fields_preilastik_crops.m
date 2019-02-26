%Program splits ASHLAR image into different sized fields and into cropped
%images for Ilastik 
%CordMatrix = [ field# row column x y coordinates];
%Outputs "Field_row_column.tif" files 
filename.folders.ashlared = 'Ashlared'; 
filename.folders.output = 'ANALYSIS\'; 
filename.folders.fullstacks = 'FullStacks\';
filename.folders.ilastikprob= 'Ilastik_Probabilities\';
filename.folders.ilastikseg = 'Ilastik_Segmentation\';
filename.folders.results = 'Analysis_Results\';
filename.folders.coordinates = 'Coordinates_FullStacks\';
filename.folders.cropfol = 'CroppedData\';

%% Edit below 
filename.folders.main = 'W:\Analysis\Ashlared\';
filename.folders.fols = {'AJ0160_P2',...
    'AJ0160_P3',...
    'AJ0160_P5',...
    'AJ0160_P6',...
    'AJ0176_P1',...
    'AJ0176_P2',...
    'AJ0176_P3',...
    'AJ0176_P4',...
    'AJ0179_P1',...
    'AJ0179_P2',...
    'AJ0179_P4',...
    'AJ0179_P5',...
    'AJ0179_P6',...
    'AJ0188_P1',...
    'AJ0188_P2',...
    'AJ0188_P4',...
    'AJ0188_P5',...
    'AJ0188_P6',...
    'AJ0189_P2',...
    'AJ0189_P4',...
    'AJ0189_P5',...
    'AJ0189_P6',...
    'AJ0205_P1',...
    'AJ0205_P3',...
    'AJ0205_P4',...
    'AJ0205_P5',...
    'AJ0205_P6'}; %Name of Ashlared image without '.ome.tif' ending 

filename.dim = ['%02d']; 
filename.sizefield = 6000; %Size of field desired 
filename.cycles = 8; % total # of cycles
filename.crops = 2; %# of cropped fields desired per field 

options.max_prob = 65535;
options.cellsize = 20;
options.date = '20190223';

%% Do not edit 
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