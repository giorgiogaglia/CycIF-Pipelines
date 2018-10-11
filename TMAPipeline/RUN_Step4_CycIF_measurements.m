function RUN_Step4_CycIF_measurements(folderloc, prefix2, prefix1, maxcycle, dates, DAPIslice)  
%% Comment out to prevent rewriting over saved mat files 

% x= input('STOP! Are you sure you want to write over the MAT file????? Y/N?');
% if x == 'Y'

filename_res = [folderloc dates  '_Results_Measurements.mat'];
mkdir([folderloc'\FociSeg'])
save(filename_res)

% end
trackfile_seg = '\TrackedImages\TrackedField';
core_rawimage = '\FullStacks\Core';
%% Measuring CycIF 

for i1 = 1:length(prefix1)

filename_res = [folderloc dates  '_Results_Measurements.mat'];

load(filename_res)

tic
Field = [];
% open the tiff file of the DAPI images and segment them separately
field = 0;
   

for i2 = 1:length(prefix2)
    
tracking_stack = [folderloc filesep 'TrackedImages\TrackedField' prefix1{i1} num2str(prefix2(i2),dim) '.tif']; %Tracked field 
rawimage_stack = [folderloc filesep 'FullStacks\Core' prefix1{i1}  num2str(prefix2(i2),dim) '.tif'];%Full Stack 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  if %want to restart program to continue to appending to Field of a slide
%  that is in complete comment out this portion or else it will write over
%  Field 
field = field + 1;
Field(field).Name = prefix2(i2);
Field(field).Area = [];
Field(field).Solidity = [];
Field(field).CentroidRow = [];
Field(field).CentroidCol = [];
Field(field).MedianNucSign = [];
Field(field).MedianCytSign = [];
Field(field).MeanNucSign = [];
Field(field).MeanCytSign = [];
Field(field).HSF1Foci_Area = [];
Field(field).HSF1Foci_Sign = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    disp(tracking_stack)
    A = imread(tracking_stack,'Index',1);
    clear A
catch
    continue
end


for cycle = 1:maxcycle
    toc

    % load images 
    try
        lb_Nuc_Image = uint16(imread(tracking_stack,'Index',cycle));
    catch
        disp(['Cycle ' num2str(cycle) 'was not found'])
        tracking_stack
        continue
    end
    
    DAPI_Image = uint16(imread(rawimage_stack,'Index',DAPIslice(cycle)));
    Image{1} = uint16(imread(rawimage_stack,'Index',DAPIslice(cycle)+1));
    Image{2} = uint16(imread(rawimage_stack,'Index',DAPIslice(cycle)+2));
    Image{3} = uint16(imread(rawimage_stack,'Index',DAPIslice(cycle)+3));
    
    % correct shift between colors 
    x_shift = 2;
    y_shift = 2;
    for j2 = 1:3
        try 
            Image_temp{j2} = padarray(Image{j2},[y_shift x_shift],0,'pre');
        catch
            continue
        end 
    end

    for j3 = 1:3
            try
                Image{j3} = Image_temp{j3}(1:length(Image{1}(:,1)),1:length(Image{1}(1,:)));
            catch
                continue
            end
    end

    clear Image_temp 
    
    % subtract background locally - ie correct for flatfield
    DAPI_BK = DAPI_Image; %-imopen(DAPI_Image, strel('disk',10));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for j4 = 1:3
        try
            Image_BK{j4} = Image{j4};
        catch
            continue
        end
    end

    
    SegImage = lb_Nuc_Image > 0;


    % create cytoplasmic mask
    lb_Cyt_Image = imdilate(lb_Nuc_Image,offsetstrel('ball',5,0));
    lb_Cyt_Image(lb_Nuc_Image>0)=0;
    CytImage = lb_Cyt_Image > 0;

  
    stats_NucImage = regionprops(lb_Nuc_Image,'Area','Solidity','Centroid');
    
    if cycle == 1 %initializing size of struct 
    totcells = length(stats_NucImage);
    Field(field).Name = prefix2(i2);
    Field(field).Area = zeros(totcells,maxcycle)+NaN;
    Field(field).Solidity = zeros(totcells,maxcycle)+NaN;
    Field(field).CentroidRow = zeros(totcells,maxcycle)+NaN;
    Field(field).CentroidCol = zeros(totcells,maxcycle)+NaN;
    Field(field).MedianNucSign = zeros(totcells,maxcycle*4)+NaN;
    Field(field).MedianCytSign = zeros(totcells,maxcycle*4)+NaN;
    Field(field).MeanNucSign = zeros(totcells,maxcycle*4)+NaN;
    Field(field).MeanCytSign = zeros(totcells,maxcycle*4)+NaN;
    end
        
    Area = {stats_NucImage.Area};
    Solidity = {stats_NucImage.Solidity};
    Centroid = {stats_NucImage.Centroid};
    CentroidVect = cell2mat(Centroid);
    CentroidMat = reshape(CentroidVect,2,length(CentroidVect)/2);
    
    totcyclecells = length(Area);
   
    Field(field).Area(1:totcyclecells,cycle) = cell2mat(Area);
    Field(field).Solidity(1:totcyclecells,cycle) = cell2mat(Solidity);
    Field(field).CentroidRow(1:totcyclecells,cycle) = CentroidMat(1,:);
    Field(field).CentroidCol(1:totcyclecells,cycle) = CentroidMat(2,:); 
    
    Field(field).Area(isnan(Field(field).Area))=0;
    Field(field).Solidity(isnan(Field(field).Solidity))=0;
    Field(field).CentroidRow(isnan(Field(field).CentroidRow))=0;
    Field(field).CentroidCol(isnan(Field(field).CentroidCol))=0;
    
    Nuclei_DAPI_stats = regionprops(lb_Nuc_Image,DAPI_BK,'PixelValues');
    Cytopl_DAPI_stats = regionprops(lb_Cyt_Image,DAPI_BK,'PixelValues');
  
    Nuclei_A488_stats = regionprops(lb_Nuc_Image,Image_BK{1},'PixelValues');
    Nuclei_A555_stats = regionprops(lb_Nuc_Image,Image_BK{2},'PixelValues');
    Nuclei_A647_stats = regionprops(lb_Nuc_Image,Image_BK{3},'PixelValues');
    
    Cytopl_A488_stats = regionprops(lb_Cyt_Image,Image_BK{1},'PixelValues');
    Cytopl_A555_stats = regionprops(lb_Cyt_Image,Image_BK{2},'PixelValues');
    Cytopl_A647_stats = regionprops(lb_Cyt_Image,Image_BK{3},'PixelValues');
    
    Field(field).MedianNucSign(1:totcyclecells,4*(cycle-1)+1) = cellfun(@median,{Nuclei_DAPI_stats.PixelValues});
    Field(field).MedianNucSign(1:totcyclecells,4*(cycle-1)+2) = cellfun(@median,{Nuclei_A488_stats.PixelValues});
    Field(field).MedianNucSign(1:totcyclecells,4*(cycle-1)+3) = cellfun(@median,{Nuclei_A555_stats.PixelValues});
    Field(field).MedianNucSign(1:totcyclecells,4*(cycle-1)+4) = cellfun(@median,{Nuclei_A647_stats.PixelValues});
    
    Field(field).MedianCytSign(1:totcyclecells,4*(cycle-1)+1) = cellfun(@median,{Cytopl_DAPI_stats.PixelValues});
    Field(field).MedianCytSign(1:totcyclecells,4*(cycle-1)+2) = cellfun(@median,{Cytopl_A488_stats.PixelValues});
    Field(field).MedianCytSign(1:totcyclecells,4*(cycle-1)+3) = cellfun(@median,{Cytopl_A555_stats.PixelValues});
    Field(field).MedianCytSign(1:totcyclecells,4*(cycle-1)+4) = cellfun(@median,{Cytopl_A647_stats.PixelValues});
    
    Field(field).MeanNucSign(1:totcyclecells,4*(cycle-1)+1) = cellfun(@mean,{Nuclei_DAPI_stats.PixelValues});
    Field(field).MeanNucSign(1:totcyclecells,4*(cycle-1)+2) = cellfun(@mean,{Nuclei_A488_stats.PixelValues});
    Field(field).MeanNucSign(1:totcyclecells,4*(cycle-1)+3) = cellfun(@mean,{Nuclei_A555_stats.PixelValues});
    Field(field).MeanNucSign(1:totcyclecells,4*(cycle-1)+4) = cellfun(@mean,{Nuclei_A647_stats.PixelValues});
    
    Field(field).MeanCytSign(1:totcyclecells,4*(cycle-1)+1) = cellfun(@mean,{Cytopl_DAPI_stats.PixelValues});
    Field(field).MeanCytSign(1:totcyclecells,4*(cycle-1)+2) = cellfun(@mean,{Cytopl_A488_stats.PixelValues});
    Field(field).MeanCytSign(1:totcyclecells,4*(cycle-1)+3) = cellfun(@mean,{Cytopl_A555_stats.PixelValues});
    Field(field).MeanCytSign(1:totcyclecells,4*(cycle-1)+4) = cellfun(@mean,{Cytopl_A647_stats.PixelValues});
    
    Field(field).MedianNucSign(isnan(Field(field).MedianNucSign))=0;
    Field(field).MedianCytSign(isnan(Field(field).MedianCytSign))=0;
    
    Field(field).MeanNucSign(isnan(Field(field).MeanNucSign))=0;
    Field(field).MeanCytSign(isnan(Field(field).MeanCytSign))=0;        

end
save(filename_res,'Field','-append') %Saving matrix 
end
end
end