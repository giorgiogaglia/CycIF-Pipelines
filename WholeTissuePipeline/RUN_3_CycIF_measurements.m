function RUN_3_CycIF_measurements(basefolder,slides_folders, maxfields, maxcycle, dates, DAPIslice, z_cycle, znum, thr, CYCLEslice, prefix1,max_rows, well_nums, well_lets)
%% Comment out to prevent rewriting over saved mat files

for folder = 1:length(slides_folders)
    filename_res = [basefolder slides_folders{folder} dates slides_folders{folder} '_Results.mat'];
    mkdir([basefolder slides_folders{folder} '\FociSeg'])
    try
        load(filename_res)
    catch
        save(filename_res)
        continue
    end
end


%% Measuring CycIF
trackfile_seg = '\TrackedImages\TrackedField';


for folder = 1:length(slides_folders)
    
    filename_res = [basefolder slides_folders{folder} dates slides_folders{folder} '_Results.mat'];
    
    load(filename_res)
    
    tic
    Field = [];
    % open the tiff file of the DAPI images and segment them separately
    
    prefix2 = linspace(1,maxfields(folder),maxfields(folder));
    field = 0;
    
    
    
    for i1 = 1:max_rows
        for i3 = 1:length(well_nums)
            
            for i2 = 1:length(prefix2)
                tracking_stack = [basefolder slides_folders{folder} trackfile_seg well_lets{i1} num2str(well_nums(i3)) '_TrackedField' num2str(prefix2(i2),'%04d') '.tif']; %Tracked field
                core_rawimage = ['\FullStacks\' well_lets{i1} num2str(well_nums(i3)) '_Field'];
                rawimage_stack = [basefolder slides_folders{folder} core_rawimage num2str(prefix2(i2),'%04d') '.tif']; %Field
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %  if %want to restart program to continue to appending to Field of a slide
                %  that is in complete comment out this portion or else it will write over
                %  Field
                field = field + 1;
                Plate{i1, well_nums(i3)}.Fields(field).Name = prefix2(i2);
                Plate{i1, well_nums(i3)}.Fields(field).Area = [];
                Plate{i1, well_nums(i3)}.Field(field).Solidity = [];
                Plate{i1, well_nums(i3)}.Field(field).CentroidRow = [];
                Plate{i1, well_nums(i3)}.Field(field).CentroidCol = [];
                Plate{i1, well_nums(i3)}.Field(field).MedianNucSign = [];
                Plate{i1, well_nums(i3)}.Field(field).MedianCytSign = [];
                Plate{i1, well_nums(i3)}.Field(field).MeanNucSign = [];
                Plate{i1, well_nums(i3)}.Field(field).MeanCytSign = [];
                Plate{i1, well_nums(i3)}.Field(field).HSF1Foci_Area = [];
                Plate{i1, well_nums(i3)}.Field(field).HSF1Foci_Sign = [];
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
                    channels = CYCLEslice{cycle};
                    % load images
                    try
                        lb_Nuc_Image = uint16(imread(tracking_stack,'Index',cycle));
                    catch
                        disp(['Cycle ' num2str(cycle) 'was not found'])
                        tracking_stack
                        continue
                    end
                    % Needs to change depending on z-stack
                    DAPI_Image = uint16(imread(rawimage_stack,'Index',DAPIslice(cycle)));
                    for i2 = 1:length(channels)-1
                        Image{i2} = uint16(imread(rawimage_stack,'Index',DAPIslice(cycle)+i2));
                    end
                    
                    % correct shift between colors
                    x_shift = 2;
                    y_shift = 2;
                    for j2 = 1:length(channels)-1
                        try
                            Image_temp{j2} = padarray(Image{j2},[y_shift x_shift],0,'pre');
                        catch
                            continue
                        end
                    end
                    
                    for j3 = 1:length(channels)-1
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
                    %Making Maximum Intensity Projection of Z-stacks
                    Image_Zstack = []; %Stack of each Z-stack
                    
                    if cycles(cycle) == z_cycle
                        for j4 = 1:znum
                            try
                                Image_BK{j4} = Image{j4};
                                Image_Zstack = cat(3,Image_Zstack,Image{j4});
                            catch
                                continue
                            end
                        end
                        MP_Zstack = max(Image_Zstack, [], 3); %Creates Maximum intensity projection of all Z-stacks
                    else
                        for j4 = 1:3
                            try
                                Image_BK{j4} = Image{j4};
                            catch
                                continue
                            end
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
                        Plate{i1, well_nums(i3)}.Field(field).Name = prefix2(i2);
                        Plate{i1, well_nums(i3)}.Field(field).Area = zeros(totcells,maxcycle)+NaN;
                        Plate{i1, well_nums(i3)}.Field(field).Solidity = zeros(totcells,maxcycle)+NaN;
                        Plate{i1, well_nums(i3)}.Field(field).CentroidRow = zeros(totcells,maxcycle)+NaN;
                        Plate{i1, well_nums(i3)}.Field(field).CentroidCol = zeros(totcells,maxcycle)+NaN;
                        Plate{i1, well_nums(i3)}.Field(field).MedianNucSign = zeros(totcells,maxcycle*4)+NaN;
                        Plate{i1, well_nums(i3)}.Field(field).MedianCytSign = zeros(totcells,maxcycle*4)+NaN;
                        Plate{i1, well_nums(i3)}.Field(field).MeanNucSign = zeros(totcells,maxcycle*4)+NaN;
                        Plate{i1, well_nums(i3)}.Field(field).MeanCytSign = zeros(totcells,maxcycle*4)+NaN;
                    end
                    
                    Area = {stats_NucImage.Area};
                    Solidity = {stats_NucImage.Solidity};
                    Centroid = {stats_NucImage.Centroid};
                    CentroidVect = cell2mat(Centroid);
                    CentroidMat = reshape(CentroidVect,2,length(CentroidVect)/2);
                    
                    totcyclecells = length(Area);
                    
                    Plate{i1, well_nums(i3)}.Field(field).Area(1:totcyclecells,cycle) = cell2mat(Area);
                    Plate{i1, well_nums(i3)}.Field(field).Solidity(1:totcyclecells,cycle) = cell2mat(Solidity);
                    Plate{i1, well_nums(i3)}.Field(field).CentroidRow(1:totcyclecells,cycle) = CentroidMat(1,:);
                    Plate{i1, well_nums(i3)}.Field(field).CentroidCol(1:totcyclecells,cycle) = CentroidMat(2,:);
                    
                    Plate{i1, well_nums(i3)}.Field(field).Area(isnan(Field(field).Area))=0;
                    Plate{i1, well_nums(i3)}.Field(field).Solidity(isnan(Field(field).Solidity))=0;
                    Plate{i1, well_nums(i3)}.Field(field).CentroidRow(isnan(Field(field).CentroidRow))=0;
                    Plate{i1, well_nums(i3)}.Field(field).CentroidCol(isnan(Field(field).CentroidCol))=0;
                    
                    Nuclei_DAPI_stats = regionprops(lb_Nuc_Image,DAPI_BK,'PixelValues');
                    
                    Cytopl_DAPI_stats = regionprops(lb_Cyt_Image,DAPI_BK,'PixelValues');
                    
                    if cycles(cycle) == z_cycle
                        %stats is for Maximum Intensity Projection of all Z-stacks
                        Nuclei_A488_stats = regionprops(lb_Nuc_Image,Image_BK{1},'PixelValues');
                        Nuclei_A555_stats = regionprops(lb_Nuc_Image,MP_Zstack,'PixelValues');
                        Nuclei_A647_stats = regionprops(lb_Nuc_Image,Image_BK{3},'PixelValues');
                        
                        Cytopl_A488_stats = regionprops(lb_Cyt_Image,Image_BK{1},'PixelValues');
                        Cytopl_A555_stats = regionprops(lb_Cyt_Image,MP_Zstack,'PixelValues');
                        Cytopl_A647_stats = regionprops(lb_Cyt_Image,Image_BK{3},'PixelValues');
                        
                        Plate{i1, well_nums(i3)}.Field(field).MedianNucSign(1:totcyclecells,4*(cycle-1)+1) = cellfun(@median,{Nuclei_DAPI_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MedianNucSign(1:totcyclecells,4*(cycle-1)+2) = cellfun(@median,{Nuclei_A488_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MedianNucSign(1:totcyclecells,4*(cycle-1)+3) = cellfun(@median,{Nuclei_A555_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MedianNucSign(1:totcyclecells,4*(cycle-1)+4) = cellfun(@median,{Nuclei_A647_stats.PixelValues});
                        
                        Plate{i1, well_nums(i3)}.Field(field).MedianCytSign(1:totcyclecells,4*(cycle-1)+1) = cellfun(@median,{Cytopl_DAPI_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MedianCytSign(1:totcyclecells,4*(cycle-1)+2) = cellfun(@median,{Cytopl_A488_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MedianCytSign(1:totcyclecells,4*(cycle-1)+3) = cellfun(@median,{Cytopl_A555_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MedianCytSign(1:totcyclecells,4*(cycle-1)+4) = cellfun(@median,{Cytopl_A647_stats.PixelValues});
                        
                        Plate{i1, well_nums(i3)}.Field(field).MeanNucSign(1:totcyclecells,4*(cycle-1)+1) = cellfun(@mean,{Nuclei_DAPI_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MeanNucSign(1:totcyclecells,4*(cycle-1)+2) = cellfun(@mean,{Nuclei_A488_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MeanNucSign(1:totcyclecells,4*(cycle-1)+3) = cellfun(@mean,{Nuclei_A555_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MeanNucSign(1:totcyclecells,4*(cycle-1)+4) = cellfun(@mean,{Nuclei_A647_stats.PixelValues});
                        
                        Plate{i1, well_nums(i3)}.Field(field).MeanCytSign(1:totcyclecells,4*(cycle-1)+1) = cellfun(@mean,{Cytopl_DAPI_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MeanCytSign(1:totcyclecells,4*(cycle-1)+2) = cellfun(@mean,{Cytopl_A488_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MeanCytSign(1:totcyclecells,4*(cycle-1)+3) = cellfun(@mean,{Cytopl_A555_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MeanCytSign(1:totcyclecells,4*(cycle-1)+4) = cellfun(@mean,{Cytopl_A647_stats.PixelValues});
                    else
                        
                        Nuclei_A488_stats = regionprops(lb_Nuc_Image,Image_BK{1},'PixelValues');
                        Nuclei_A555_stats = regionprops(lb_Nuc_Image,Image_BK{2},'PixelValues');
                        Nuclei_A647_stats = regionprops(lb_Nuc_Image,Image_BK{3},'PixelValues');
                        
                        Cytopl_A488_stats = regionprops(lb_Cyt_Image,Image_BK{1},'PixelValues');
                        Cytopl_A555_stats = regionprops(lb_Cyt_Image,Image_BK{2},'PixelValues');
                        Cytopl_A647_stats = regionprops(lb_Cyt_Image,Image_BK{3},'PixelValues');
                        
                        Plate{i1, well_nums(i3)}.Field(field).MedianNucSign(1:totcyclecells,4*(cycle-1)+1) = cellfun(@median,{Nuclei_DAPI_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MedianNucSign(1:totcyclecells,4*(cycle-1)+2) = cellfun(@median,{Nuclei_A488_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MedianNucSign(1:totcyclecells,4*(cycle-1)+3) = cellfun(@median,{Nuclei_A555_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MedianNucSign(1:totcyclecells,4*(cycle-1)+4) = cellfun(@median,{Nuclei_A647_stats.PixelValues});
                        
                        Plate{i1, well_nums(i3)}.Field(field).MedianCytSign(1:totcyclecells,4*(cycle-1)+1) = cellfun(@median,{Cytopl_DAPI_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MedianCytSign(1:totcyclecells,4*(cycle-1)+2) = cellfun(@median,{Cytopl_A488_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MedianCytSign(1:totcyclecells,4*(cycle-1)+3) = cellfun(@median,{Cytopl_A555_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MedianCytSign(1:totcyclecells,4*(cycle-1)+4) = cellfun(@median,{Cytopl_A647_stats.PixelValues});
                        
                        Plate{i1, well_nums(i3)}.Field(field).MeanNucSign(1:totcyclecells,4*(cycle-1)+1) = cellfun(@mean,{Nuclei_DAPI_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MeanNucSign(1:totcyclecells,4*(cycle-1)+2) = cellfun(@mean,{Nuclei_A488_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MeanNucSign(1:totcyclecells,4*(cycle-1)+3) = cellfun(@mean,{Nuclei_A555_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MeanNucSign(1:totcyclecells,4*(cycle-1)+4) = cellfun(@mean,{Nuclei_A647_stats.PixelValues});
                        
                        Plate{i1, well_nums(i3)}.Field(field).MeanCytSign(1:totcyclecells,4*(cycle-1)+1) = cellfun(@mean,{Cytopl_DAPI_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MeanCytSign(1:totcyclecells,4*(cycle-1)+2) = cellfun(@mean,{Cytopl_A488_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MeanCytSign(1:totcyclecells,4*(cycle-1)+3) = cellfun(@mean,{Cytopl_A555_stats.PixelValues});
                        Plate{i1, well_nums(i3)}.Field(field).MeanCytSign(1:totcyclecells,4*(cycle-1)+4) = cellfun(@mean,{Cytopl_A647_stats.PixelValues});
                    end
                    
                    Plate{i1, well_nums(i3)}.Field(field).MedianNucSign(isnan(Field(field).MedianNucSign))=0;
                    Plate{i1, well_nums(i3)}.Field(field).MedianCytSign(isnan(Field(field).MedianCytSign))=0;
                    
                    Plate{i1, well_nums(i3)}.Field(field).MeanNucSign(isnan(Field(field).MeanNucSign))=0;
                    Plate{i1, well_nums(i3)}.Field(field).MeanCytSign(isnan(Field(field).MeanCytSign))=0;
                    
                    if cycles(cycle) == z_cycle
                        % only use HSF1cycle to segment the foci
                        if totcyclecells > 0
                            % initialize results variables
                            Field(field).HSF1Foci_Area = zeros(totcyclecells,2)+NaN;
                            Field(field).HSF1Foci_Sign = zeros(totcyclecells,2)+NaN;
                            Dilate_Nuc_Image = lb_Nuc_Image;
                            FociSeg_Cumulative= 0; % = [];
                            FociZstack = []; %Stack of all z-stack images
                            for foci_im = 1:znum
                                % segement HSF1 foci
                                try
                                    HSF1Foci_BK = Image{foci_im};%-imopen(Image{foci_im}, strel('disk',15));
                                catch
                                    continue
                                end
                                
                                SpatVar_Raw = stdfilt(HSF1Foci_BK)./mean(Image{foci_im}( Image{foci_im}(:)>10  ));
                                
                                % FUNDAMENTAL PARAMETER HERE!!!!
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                FociSeg1 = SpatVar_Raw>thr;       % we increase this parameter to be more conservative (from 0.8 to 2)-make very conservative; FociSeg1 is segmentation for single z-stack
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %Continue analysis on 1 HSF1 Z-stack
                                lb_FociImage1 = uint16(FociSeg1).*Dilate_Nuc_Image;
                                lb_FociImage1 = imfill(lb_FociImage1,'holes');
                                FociSegImage1 = lb_FociImage1 > 0;
                                FociSegImage1 = FociSegImage1 - bwareaopen(FociSegImage1,150);
                                L = bwlabel(FociSegImage1);
                                stats_foci = regionprops(L,'Area','Solidity');
                                idx = find([stats_foci.Solidity] > 0.7);
                                L = ismember(L,idx);
                                FociSegImage1 = L > 0;
                                FociCheckImage1 = uint8(SegImage);
                                FociCheckImage1(FociSegImage1>0)=2;
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                FociSeg_Cumulative = FociSeg_Cumulative + uint16(FociSegImage1); %adding all the foci segmentation for all z-stacks together
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %Save stacked image of all z-stacks
                                FociZstack = cat(3,uint16(FociZstack),uint16(FociCheckImage1));
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % THIS VERSION OF THE ANALYSIS HAS NO DILATION!!!!
                            %             FociSeg = FociSeg.*imdilate(uint16(SegImage),offsetstrel('ball',3,0));
                            %FociSeg = FociSeg > 0;
                            FociSeg = FociSeg_Cumulative>0;
                            lb_FociImage = uint16(FociSeg).*Dilate_Nuc_Image;
                            
                            FociCheckImage = uint8(SegImage);
                            FociCheckImage(FociSeg>0)=2;
                            
                            DirectCheckIm = cat(3,uint16(FociCheckImage),uint16(MP_Zstack));
                            DirectCheckIm = cat(3,DirectCheckIm,uint16(DAPI_BK));
                            
                            imwrite(uint16(DirectCheckIm),[basefolder slides_folders{folder} '\FociSeg\Field' num2str(prefix2(i2),'%04d') '_FociSeg_SV1_normfit.tif'])  %Foci segmentation to check
                            
                            FociZstack = cat(3,uint16(FociCheckImage),uint16(FociZstack)); %Stack of 1st=Cumulative Foci segmentation, 2nd... Each z-stack
                            filez = [basefolder slides_folders{folder} '\FociSeg\' prefix1(i1) 'Field' num2str(prefix2(i2),'%04d') '_FociSeg_Z-Stack.tif']; %File name of z-stack file
                            
                            saveastiff(uint16(FociZstack),filez);   %Saving tiff file of Foci segmentation for all z-stacks
                            
                            Foci_stats = regionprops(lb_FociImage,MP_Zstack,'PixelValues'); %uses maximum intensity projection of all images
                            Foci_tot_stats = regionprops(Dilate_Nuc_Image,MP_Zstack,'PixelValues');
                            
                            Plate{i1, well_nums(i3)}.Field(field).HSF1Foci_Area(1:max(lb_FociImage(:)),1) = cellfun(@length,{Foci_stats.PixelValues});
                            Plate{i1, well_nums(i3)}.Field(field).HSF1Foci_Area(1:max(Dilate_Nuc_Image(:)),2) = cellfun(@length,{Foci_tot_stats.PixelValues});
                            Plate{i1, well_nums(i3)}.Field(field).HSF1Foci_Sign(1:max(lb_FociImage(:)),1) = cellfun(@sum,{Foci_stats.PixelValues});
                            Plate{i1, well_nums(i3)}.Field(field).HSF1Foci_Sign(1:max(Dilate_Nuc_Image(:)),2) = cellfun(@sum,{Foci_tot_stats.PixelValues});
                            
                            Plate{i1, well_nums(i3)}.Field(field).HSF1Foci_Area(isnan(Field(field).HSF1Foci_Area))=0;
                            Plate{i1, well_nums(i3)}.Field(field).HSF1Foci_Sign(isnan(Field(field).HSF1Foci_Sign))=0;
                            
                        end
                    end
                end
                save(filename_res,'Plate','-append') %Saving matrix
            end
        end
    end
end
