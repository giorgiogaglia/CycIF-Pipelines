function RUN_Step4_CycIF_measurements_ilastik(filename, options)  



addpath(filename.folders.main)
mkdir([filename.folders.main  filename.folders.output 'Analysis_Results\']);
HSF1cycle = options.hsf1round;     
tic
core = 0;

for k = 1:length(filename.folders.fols)
    fold = filename.folders.fols{k};
    filename.folders.resultfile = [filename.folders.results fold '_Results_' options.date '.mat'];
    addpath(filename.folders.main)
    mkdir([filename.folders.main filename.folders.output fold filesep filename.folders.fociseg]); 
    for i1 = 1:filename.prefix1
        for i2 = 1:filename.prefix2
            
            CoreFile = [fold '_Field_' num2str(i1 , filename.dim) '_' num2str(i2 , filename.dim)];
            FileTif = [filename.folders.main filename.folders.output fold filesep filename.folders.fullstacks CoreFile  filename.suffix];
            try
                DAPI = imread(FileTif,'Index',1);
            catch
                continue
            end
            disp(CoreFile)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  if %want to restart program to continue to appending to Field of a slide
            %  that is in complete comment out this portion or else it will write over
            %  Field
            core = core + 1; 
            Field(core).Name = ['Field_' num2str(i1 , filename.dim) '_' num2str(i2 , filename.dim)];
            Field(core).Area = [];
            Field(core).Solidity = [];
            Field(core).CentroidRow = [];
            Field(core).CentroidCol = [];
            Field(core).MedianNucSign = [];
            Field(core).MedianCytSign = [];
            Field(core).MeanNucSign = [];
            Field(core).MeanCytSign = [];
            Field(core).HSF1Foci_Area = [];
            Field(core).HSF1Foci_Sign = [];
            
            % load nuclear mask, create cytoplasmic mask and save them
            segfile = [filename.folders.main filename.folders.output fold filesep filename.folders.ilastikseg CoreFile '_Seg.tif'];
            
            if exist(segfile,'file') == 2
                disp('yes') 
                NucMask = (imread(segfile));
                lb_NucMask = uint16(bwlabel(NucMask));
                
                % create cytoplasmic mask
                lb_CytMask = imdilate(lb_NucMask,offsetstrel('ball',5,0));
                lb_CytMask(lb_NucMask>0)=0;
                CytImage = lb_CytMask > 0;
                
                cyt_segfile = [filename.folders.main filename.folders.output fold filesep filename.folders.ilastikseg CoreFile '_NucCytSeg.tif'];
                DAPI_img = uint16(imread(FileTif,'Index',1));
                NucCytSeg = cat(3,lb_NucMask,DAPI_img,lb_CytMask);
                imwrite(NucCytSeg,cyt_segfile);
                
                % measure properties and assign them to results
                stats_nuc = regionprops(lb_NucMask,'Area','Centroid','Solidity');
                
                Area = {stats_nuc.Area};
                Solidity = {stats_nuc.Solidity};
                Centroid = {stats_nuc.Centroid};
                CentroidVect = cell2mat(Centroid);
                CentroidMat = reshape(CentroidVect,2,length(CentroidVect)/2);
                
                Field(core).Area = cell2mat(Area);
                Field(core).Solidity = cell2mat(Solidity);
                Field(core).CentroidRow = CentroidMat(1,:);
                Field(core).CentroidCol = CentroidMat(2,:);
                
                
                for i3 = 1:filename.cycles*4
                    % load image
                    FluorImage = imread(FileTif,'Index',i3);
                    FluorImage_BK = FluorImage - imopen(imclose(FluorImage,strel('disk',round(options.cellsize/10))),strel('disk',round(options.cellsize*3)));
                    
                    Nuclei_stats = regionprops(lb_NucMask,FluorImage_BK,'PixelValues');
                    Cytopl_stats = regionprops(lb_CytMask,FluorImage_BK,'PixelValues');
                    
                    if length(Nuclei_stats)>0
                        Field(core).MedianNucSign(:,i3) = cellfun(@median,{Nuclei_stats.PixelValues});
                        Field(core).MedianCytSign(:,i3) = cellfun(@median,{Cytopl_stats.PixelValues});
                        Field(core).MeanNucSign(:,i3) = cellfun(@mean,{Nuclei_stats.PixelValues});
                        Field(core).MeanCytSign(:,i3) = cellfun(@mean,{Cytopl_stats.PixelValues});
                        
                         if i3 == HSF1cycle % only use HSF1cycle to segment the foci
                             % segement HSF1 foci
                             try
                                 HSF1Foci_BK = FluorImage{HSF1round};%-imopen(Image{foci_im}, strel('disk',15));
                             catch
                                 continue
                             end
                             
                             SpatVar_Raw = stdfilt(HSF1Foci_BK)./mean(FluorImage{HSF1round}(FluorImage{HSF1round}(:)>10  ));
                             
                             % FUNDAMENTAL PARAMETER HERE!!!!
                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                             FociSeg1 = SpatVar_Raw>options.focithr;       % we increase this parameter to be more conservative (from 0.8 to 2)-make very conservative; FociSeg1 is segmentation for single z-stack
                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                             lb_FociImage1 = uint16(FociSeg1).*lb_NucMask;
                             lb_FociImage1 = imfill(lb_FociImage1,'holes');
                             FociSegImage1 = lb_FociImage1 > 0;
                             FociSegImage1 = FociSegImage1 - bwareaopen(FociSegImage1,150);
                             L = bwlabel(FociSegImage1);
                             stats_foci = regionprops(L,'Area','Solidity');
                             idx = find([stats_foci.Solidity] > 0.7);
                             L = ismember(L,idx);
                             FociSegImage1 = L > 0;
                             FociCheckImage1 = uint8(NucMask);
                             FociCheckImage1(FociSegImage1>0)=2;
                             FociSeg_Cumulative = uint16(FociSegImage1); 
                             FociSeg = FociSeg_Cumulative>0;
                             lb_FociImage = uint16(FociSeg).*lb_NucMask;
                             
                             FociCheckImage = uint8(NucMask);
                             FociCheckImage(FociSeg>0)=2;
                             
                             DirectCheckIm = cat(3,uint16(FociCheckImage),uint16(FluorImage{HSF1round}));
                             DirectCheckIm = cat(3,DirectCheckIm,uint16(DAPI_img));
                             
                             imwrite(uint16(DirectCheckIm),[filename.folders.main filename.folders.output fold filesep filename.folders.fociseg CoreFile '_HSF1_FociSeg_check.tif'])  %Foci segmentation to check
                             
                             FociZstack = cat(3,uint16(FociCheckImage),uint16(FociCheckImage1));
                             filez = [filename.folders.main filename.folders.output fold filesep filename.folders.fociseg CoreFile '_HSF1_FociSeg.tif']; 
                             imwrite(uint16(FociSeg), filez); %Saving TIFF file of foci segmentation 

                             Foci_stats = regionprops(lb_FociImage,FluorImage{HSF1round},'PixelValues'); %uses maximum intensity projection of all images
                             Foci_tot_stats = regionprops(lb_NucMask,FluorImage{HSF1round},'PixelValues');
                             
                             Field(core).HSF1Foci_Area(1:max(lb_FociImage(:)),1) = cellfun(@length,{Foci_stats.PixelValues});
                             Field(core).HSF1Foci_Area(1:max(lb_NucMask(:)),2) = cellfun(@length,{Foci_tot_stats.PixelValues});
                             Field(core).HSF1Foci_Sign(1:max(lb_FociImage(:)),1) = cellfun(@sum,{Foci_stats.PixelValues});
                             Field(core).HSF1Foci_Sign(1:max(lb_NucMask(:)),2) = cellfun(@sum,{Foci_tot_stats.PixelValues});
                         end 
                    end
                    
                end
            end
            save([filename.folders.main filename.folders.output filename.folders.resultfile],'Field') %Saving matrix
        end
    end
end