function remove_ROI(slides_folders, slides_rowtils, slides_coltils, maxfields,montfolder, montslides_folders) 
%This program find the fields that correspond to ROI's that are incorrectly
%indentified as cancer by CD138 
filename_res = 'Incorrect_ROIs.mat';
save(filename_res)
    
%% Creating structure of ROIs
for i2= 1:length(slides_folders) %Looping through each slide
    if isequal(montslides_folders{i2}, 'n/a')
        ROI_Fields(i2).Name = slides_folders{i2};
        ROI_Fields(i2).NotCancer = [];
    else
        montfile = imread(montslides_folders{i2}); 
        [colsize, rowsize] = size(montfile); 
        npix = rowsize(i2)/slides_coltils(i2); %# of pixels the reduced field is
        colf = 0:npix:rowsize(i2);   %Making vectors of the indexes of each field
        rowf = 0:npix:colsize(i2);
        colf(1) = 1;
        rowf(1) = 1;
        ROI_File =[montfolder montslides_folders{i2}];
        ROI_mont = ReadImageJROI(ROI_File);
        ROI_img = zeros(colsize(i2), rowsize(i2)); %Creating zero matrix of montage
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
            ROI_img(sub2ind([colsize(i2) rowsize(i2)], rows(rows.*cols>0), cols(rows.*cols>0)))= 1;
            %Connecting two points in the matrix and making connection = 1
            for j1 = 1:length(rows)-1
                xline = round(linspace(rows(j1), rows(j1+1),50));
                yline = round(linspace(cols(j1), cols(j1+1),50));
                idx=sub2ind([colsize(i2) rowsize(i2)],xline,yline);
                ROI_img(idx)=1;
            end
        end
        ROI_img = imdilate(ROI_img,strel('disk',6,8));
        ROI_img = imfill(ROI_img,'holes');
        naf = [];
        x=[];
        y=[];
        for i5 = 1:maxfields(i2)
            r1 = ceil(i5/slides_coltils(i2)); %Row number
            c1 =slides_coltils(i2)+(i5-(slides_coltils(i2))*r1); %Col number
            rind = rowf(r1:r1+1);
            cind = colf(c1:c1+1);
            if any(any(ROI_img(rind(1):rind(2), cind(1):cind(2))))~=0
                naf = [naf, i5]; %Vector of fields that are not cancer
                x=[x rind];
                y=[y cind];
            end
            
        end
        ROI_Fields(i2).Name = slides_folders{i2};
        ROI_Fields(i2).NotCancer = naf;
        figure
        imshow(ROI_img)
        %Checking if correct fields have been removed
        mont = ones(slides_rowtils(i2), slides_coltils(i2))';
        mont(naf) = 20;
        figure
        image(mont') %transpose cause matlab indexes by consecutive columns not rows
    end
    save(filename_res, 'ROI_Fields', '-append')
end
