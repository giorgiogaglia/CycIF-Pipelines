function RUN_Step3_segmentfromilastik(filename, options) 
% open the tiff file of the full stack images and segments them based on
% Ilastik Probabilities file 
% Inputs FullStacks and Ilastik Probabilities files 

tic

for k = 1:length(filename.folders.fols) 
    fold = filename.folders.fols{k}; 
    addpath(filename.folders.main)
    mkdir([filename.folders.main filename.folders.output fold filesep  filename.folders.ilastikseg])
    for i1 = 1:filename.prefix1
        for i2 = 1:filename.prefix2
            core = [fold '_Field_' num2str(i1 , filename.dim) '_' num2str(i2 , filename.dim)];
            FileTif = [filename.folders.main filename.folders.output fold filesep filename.folders.fullstacks core  filename.suffix];
            try
                DAPI = imread(FileTif,'Index',1);
            catch
                continue
            end
            disp(FileTif)
            
            IlastikTif = [filename.folders.main filename.folders.output fold filesep filename.folders.fullstacks filename.folders.ilastikprob core filename.ilastiksuffix];
            try
                IlastikRGB = imread(IlastikTif);
            catch
                continue
            end

           NucProb = IlastikRGB(:,:,options.nuc);
           CytProb = IlastikRGB(:,:,options.cyt);
           BkgProb = IlastikRGB(:,:,options.backgr);
           
            
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
                
            else
                Nuclei = zeros(size(DAPI));
                seeds = Nuclei;
            end
           
           check_img = cat(3,uint16(Nuclei),uint16(DAPI),uint16(seeds));
           check_file = [filename.folders.main filename.folders.output fold filesep filename.folders.ilastikseg core '_checkseg.tif'];
           imwrite(check_img,check_file)
           
           segfile = [filename.folders.main filename.folders.output fold filesep filename.folders.ilastikseg core '_Seg.tif'];
           imwrite(uint16(Nuclei),segfile)
           
           
           
        end
    end
end 
