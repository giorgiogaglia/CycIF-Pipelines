%function RUN_Step6_PlottingResults(filename, options) 

for i3= 1:length(filename.folders.fols) %Looping through each slide
    fold = filename.folders.fols{i3};
    coordmat = [filename.folders.main filename.folders.output filename.folders.coordinates fold '.mat'];
    load(coordmat, 'Coordinates')
    [r, c] = size(Coordinates.Field);
    filename.folders.resultfile = [filename.folders.results fold '_Results_' options.date '.mat'];
    filename_results = [filename.folders.main filename.folders.output filename.folders.resultfile];
    filename_roi = [filename.folders.main filename.folders.output filename.folders.ROI fold '_ROIpixels.mat'];
    load(filename_results, 'Field')
    load(filename_roi, 'ROI_pixels')
    for k = 1:length(Field)
        disp(Field(k).Name)
        Centroids =[];
        ROI_loc = [];
        CytMatrix =[];
        NucMatrix =[];
        Solidity = [];
        FociMatrix = [];
        cyt = [];
        nuc = [];
        foci = [];
        if Field(k).Name == ROI_pixels(k).Name
            temp_sol = Field(k).Solidity;
            Centroids(:,1) = round(Field(k).CentroidRow);%Centroid locations
            Centroids(:,2) = round(Field(k).CentroidCol);
            try
                ROI_loc = ROI_pixels(k).PixelList.PixelList(:,:); %ROI locations
            catch
                [r1, c1]= size(Centroids);
                ROI_loc = zeros(r1, c1);
            end
            index=ismember(Centroids, ROI_loc, 'rows'); %index of centroids within the ROI- specifies based on rows
            index1 = ismember(Centroids(:,1), ROI_loc(:,1));
            index2 = ismember(Centroids(:,2), ROI_loc(:,2));
            %foci(:,1) = Field(k).HSF1Foci_Sign(index,1)./Field(k).HSF1Foci_Sign(index,2); % HSF1 foci fraction
            foci(:,2) = Field(k).MedianNucSign(index,options.hsf1round);    % HSF1 Foci signal
            for j1 = 1:(filename.cycles*4)
                cyt(:,j1) = Field(k).MedianCytSign(index,j1); %Median Cytoplasm Signal
                nuc(:,j1) = Field(k).MedianNucSign(index,j1); %Median Nucleus Signal
            end
            FociMatrix = [FociMatrix; foci];
            CytMatrix = [CytMatrix; cyt];
            NucMatrix = [NucMatrix, nuc];
            Solidity = [Solidity; temp_sol];
            
            %IF YOU WANT TO DOUBLE CHECK THE ROIS 
            figure
            img = ['Z:\sorger\data\IN_Cell_Analyzer_6000\Giorgio\2019-01-03_Mouse_Lung_Tumors_Round2\Analysis\ANALYSIS\FullStacks\' fold '_' Field(k).Name '.tif'];
            img2 = imread(img,1); %DAPI first cycle 
            imshow(img2 ,[])
            hold on 
            %plot(ROI_loc(:,1), ROI_loc(:,2), 'b.')
            plot(Centroids(:,1), Centroids(:,2), 'b.')
            hold on
            %plot(ROI_loc(:,1), ROI_loc(:,2), 'y.')
            plot(Centroids(index,1), Centroids(index,2), 'y.')
            hold off
        end
    end
end


                    
            
            