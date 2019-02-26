function RUN_Step4_CycIF_measurements_ilastik(filename, options)  

% clear all  
% filename.folders.main = 'W:\Analysis\AJ0160_P2\';
% filename.folders.fullstacks = 'FullStacks\';
% filename.folders.ilastikprob= 'Ilastik_Probabilities\';
% filename.folders.ilastikseg = 'Ilastik_Segmentation\';
% filename.folders.fols = {'AJ0160_P2'};
% filename.folders.results = 'Analysis_Results\';
% 
% 
% filename.dim = ['%02d']; 
% filename.prefix1 = 10; % rows of Fields
% filename.prefix2 = 10; % columns of Fields 
% filename.cycles = 8; % max cycle
% filename.suffix = '.tif';
% filename.ilastiksuffix = '_Probabilities.tif';
% 
% options.max_prob = 65535;
% options.cellsize = 20;
% options.date = '20190223';

filename.folders.resultfile = [filename.folders.results 'Results_' options.date '.mat'];

addpath(filename.folders.main)
mkdir([filename.folders.main 'Analysis_Results\']);
    
tic
core = 0;

for k = 1:length(filename.folders.fols)
    fold = filename.folders.fols{k};
    for i1 = 1:filename.prefix1
        for i2 = 1:filename.prefix2
            
            CoreFile = [fold '_Field_' num2str(i1 , filename.dim) '_' num2str(i2 , filename.dim)];
            FileTif = [filename.folders.main filename.folders.output fold filesep filename.folders.fullstacks CoreFile  filename.suffix];
            disp(FileTif)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  if %want to restart program to continue to appending to Field of a slide
            %  that is in complete comment out this portion or else it will write over
            %  Field
            core = (i1-1)*length(filename.prefix2)+i2;
            Field(core).Name = ['Field_' num2str(i1 , filename.dim) '_' num2str(i2 , filename.dim)];
            Field(core).Area = [];
            Field(core).Solidity = [];
            Field(core).CentroidRow = [];
            Field(core).CentroidCol = [];
            Field(core).MedianNucSign = [];
            Field(core).MedianCytSign = [];
            Field(core).MeanNucSign = [];
            Field(core).MeanCytSign = [];
            
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
                    end
                end
            end
            save([filename.folders.main filename.folders.output filename.folders.resultfile],'Field') %Saving matrix
        end
    end
end