function RUN_Step3_segmentfromilastik(filename, options) 
% open the tiff file of the full stack images and segment them separately
% Inputs "Field_row_column.tif" files 
% clear all  
% filename.folders.main = 'X:\sorger\data\IN_Cell_Analyzer_6000\Giorgio\2019-01-03_Mouse_Lung_Tumors_Round2\Analysis\AJ0160_P2\';
% filename.folders.output = 'Analysis\';
% filename.folders.fullstacks = 'FullStacks\';
% filename.folders.ilastikprob= 'Ilastik_Probabilities\';
% filename.folders.ilastikseg = 'Ilastik_Segmentation\';
% filename.folders.fols = {'AJ0160_P2'};
% 
% filename.dim = ['%02d']; 
% filename.prefix1 = 20; % rows of Fields
% filename.prefix2 = 20; % columns of Fields 
% filename.cycles = 1:8; % cycles 
% filename.suffix = '.tif';
% filename.ilastiksuffix = '_Probabilities.tif';
% 
% options.max_prob = 65535;
% options.cellsize = 20;

addpath(filename.folders.main)
mkdir([filename.folders.main 'Ilastik_Segmentation\'])
    
tic
count2 = 0;

for k = 1:length(filename.folders.fols) 
    fold = filename.folders.fols{k}; 
    for i1 = 1:filename.prefix1
        for i2 = 1:filename.prefix2
            core = [fold '_Field_' num2str(i1 , filename.dim) '_' num2str(i2 , filename.dim)];
            FileTif = [filename.folders.main filename.folders.output fold filesep filename.folders.fullstacks core  filename.suffix];
            disp(FileTif)
            
            
            IlastikTif = [filename.folders.main filename.folders.output fold filesep filename.folders.ilastikprob core filename.ilastiksuffix];
            if exist(IlastikTif,'file') == 2
                try
                    IlastikRGB = imread(IlastikTif);
                catch
                    continue
                end
            else
                try
                    IlastikTif = [filename.folders.main fold filesep core filesep filename.folders.ilastikprob filename.ilastiksuffix];
                    IlastikRGB = imread(IlastikTif);
                catch
                    disp('Field does not exist')
                    continue
                end
            end
            
           NucProb = IlastikRGB(:,:,1);
           CytProb = IlastikRGB(:,:,2);
           BkgProb = IlastikRGB(:,:,3);
           DAPI = imread(FileTif,'Index',1);
            
           % define the nuclear area less stringently
           
           ThreshImg = NucProb > options.max_prob*0.2 & CytProb < options.max_prob*0.95 & BkgProb < options.max_prob*0.50;
           NucVect = DAPI(ThreshImg);
            if length(NucVect)>200
               
                ThreshImg = imfill(ThreshImg,'holes');
           
                % watershed image
                coImage=imopen(imclose(DAPI,strel('disk',round(options.cellsize/10))),strel('disk',round(options.cellsize*2)));
                sImage = DAPI-coImage;
            
                % now define seeds very stringently
                seeds_img = NucProb > options.max_prob*0.5 & CytProb < options.max_prob*0.5;
                seeds_img = imfill(seeds_img,'holes');
                seeds=seeds_img&ThreshImg;
                
                  
                L=watershed(imimposemin(-sImage,seeds));
                Nuclei=ThreshImg&L;
                Nuclei=bwareaopen(Nuclei,options.cellsize);
                
%                 figure
%                 imshow(Nuclei+seeds_img,[])
%            dddd
            else
                Nuclei = zeros(size(DAPI));
                seeds = Nuclei;
            end
           
           check_img = cat(3,uint16(Nuclei),uint16(DAPI),uint16(seeds));
           check_file = [filename.folders.main filename.folders.output fold filesep  filename.folders.ilastikseg core 'checkseg.tif'];
           imwrite(check_img,check_file)
           
           segfile = [filename.folders.main filename.folders.output fold filesep filename.folders.ilastikseg core '_Seg.tif'];
           imwrite(uint16(Nuclei),segfile)
           
           
           
        end
    end
end 
