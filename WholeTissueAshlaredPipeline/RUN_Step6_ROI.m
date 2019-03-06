%function RUN_Step5_ROI(slides_folders, slides_rowtils, slides_coltils, maxfields,montfolder, montslides_folders) 
filename.folders.main = 'Z:\sorger\data\IN_Cell_Analyzer_6000\Giorgio\2019-01-03_Mouse_Lung_Tumors_Round2\Analysis\Ashlared\';
filename.folders.fols = {'AJ0160_P2'};
filename.folders.output = 'ANALYSIS\'; 
filename.cycles = 8;
filename.montfolder = 'Z:\sorger\data\IN_Cell_Analyzer_6000\Giorgio\2019-01-03_Mouse_Lung_Tumors_Round2\Analysis\ANALYSIS\MontageforROI\';
% ROI sets must be saved as the name of the ome.tif it corresponds to 
    
%% Creating structure of ROIs
for i2= 1:length(filename.folders.fols) %Looping through each slide
    if isequal(filename.montfolder, 'n/a')
        ROI_Fields(i2).Name = filename.folders.fols{i2};
        ROI_Fields(i2).NotCancer = [];
    else
        fold = filename.folders.fols{i2};
        filename_res = [filename.folders.main filename.folders.output 'ROIpixels\' fold '_ROIpixels.mat']; 
        coordmat = [filename.folders.main filename.folders.output filename.folders.coordinates fold];
        load(coordmat) 
        orowsize = x; %Size of largest image of omi tif 
        ocolsize = y;
        
        disp(fold) 
        %montfile = [filename.folders.main filename.folders.output filename.folders.montage fold '_montage.tif'];
        montfile = [filename.montfolder fold '_montage.tif'];
        montinfo = imfinfo(montfile); 
        mrowsize = montinfo.Width;
        mcolsize = montinfo.Height; 
        
         

        mscale1 = orowsize/mrowsize; % Scale of image based on row
        mscale2 = ocolsize/mcolsize; % Scale of image based on column 
        
        if round(mscale1) ~= round(mscale2)
            disp("Scale is not matching!!") 
            continue
        else 
            mscale = round(mscale1); %Scale of image 
        end 
        
        %% Getting ROI coordinates in decreased slide Montage and filling in ROI 
        ROI_File =[filename.montfolder filename.folders.fols{i2} '.roi']; %one ROI 
        if exist(ROI_File, 'file') == 2
            ROI_mont = ReadImageJROI(ROI_File);
        else 
            ROI_File =[filename.montfolder filename.folders.fols{i2} '.zip']; %Multiple ROI's are saved in a zip file 
            ROI_mont = ReadImageJROI(ROI_File);
        end 
        
        ROI_img = zeros(mcolsize(i2), mrowsize(i2)); %Creating zero matrix of montage
        rows = [];
        cols = []; 
        for i3 = 1:length(ROI_mont)
            if length(ROI_mont) ==1
                rows = ROI_mont.mnCoordinates(:,2);
                rows(find(rows==0))=1;
                cols = ROI_mont.mnCoordinates(:,1);
                cols(find(cols==0))=1;
            else
                rows = ROI_mont{1,i3}.mnCoordinates(:,2);
                rows(find(rows==0))=1;
                cols = ROI_mont{1,i3}.mnCoordinates(:,1);
                cols(find(cols==0))=1;
            end
            rows = [rows(end); rows];   %Adding last point to beginning of vector in order for all the points to be connected to each other
            cols = [cols(end); cols];
            ROI_img(sub2ind([mcolsize(i2) mrowsize(i2)], rows(rows.*cols>0), cols(rows.*cols>0)))= 1;
            %Connecting two points in the matrix and making connection = 1
            for j1 = 1:length(rows)-1
                xline = round(linspace(rows(j1), rows(j1+1),50));
                yline = round(linspace(cols(j1), cols(j1+1),50));
                idx=sub2ind([mcolsize(i2) mrowsize(i2)],xline,yline);
                ROI_img(idx)=1;
            end
        end
        ROI_img = imdilate(ROI_img,strel('disk',6,8));
        ROI_img = imfill(ROI_img,'holes');
%         figure(1)
%         imshow(ROI_img, [])
        ROI_img = imresize(ROI_img, mscale); 
%         figure(2)
%         imshow(ROI_img, [])
        ROI_pixels = regionprops(ROI_img,'PixelList');
        pixs = cat(1, ROI_pixels.PixelList);
        figure(3)
        imshow(ROI_img,[])
        hold on 
        plot(pixs(:,1),pixs(:,2), 'b*')
    end 
%% Dividing ROI_pixels into matrices corresponding to each field in the slide 
ROI_pixels = {}; 
field = 0; 
for j1 = 1:filename.prefix1
    for j2 = 1:filename.prefix2
        coordinfo = Coordinates.Field{j1,j2}; %Coordinates matrix info for one field: [keep field (0= no, 1=yes), x1, x2, y1, y2]
        if coordinfo(1) == 1 %Field was kept 
            filename_res = [fold '_Field_' num2str(j1, filename.dim) '_' num2str(j2 , filename.dim) '_ROIs.mat']; 
            x1 = coordinfo(2);
            x2 = coordinfo(3);
            y1 = coordinfo(4); 
            y2 = coordinfo(5);
            roifield = ROI_img(y1:y2, x1:x2); %y= column size, x= row size
            field = field + 0; 
            ROI_pixels(field).Name = ['Field_' num2str(j1, filename.dim) '_' num2str(j2 , filename.dim)]; 
            ROI_pixels(field).PixelList= regionprops(ROI_img,'PixelList'); 
        end 
    end
    save(filename_res, 'ROI_pixels')
end
end 
