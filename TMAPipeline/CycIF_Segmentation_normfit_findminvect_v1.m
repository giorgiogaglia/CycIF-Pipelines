function OutImage = CycIF_Segmentation_normfit_findminvect_v1(image_in,filename_out,options)

writeflag = options.writeflag;
fixbrokencellsflag = options.fixbrokencellsflag;
sigma = options.sigma;
cellsize = options.cellsize;
background = options.backgorund/10;
method = options.findminvect_method;
quantile = options.findminvect_quantile;
fold = options.findminvect_fold;
smallmorphs = options.smallmorphs;


    % load image to be segmented if not loaded already
    if ischar(image_in)
        RawImage = double(imread(image_in));
    else
        RawImage = double(image_in);
    end

    
    %       *IMAGE PROCESSING GOES HERE.*
    
    % --------------------------- FINDING NUCLEI ------------------------------
    % takes the DNA image and gives back a "Nuclei" image
    
    
    % first do local background subtraction
    coImage=imopen(imclose(RawImage,strel('disk',round(cellsize/20))),strel('disk',round(cellsize*1)));
    sImage = RawImage-coImage;
   
    
    % then determine the threshold by fitting two gaussians
    TwoNormModel = fitgmdist(log2(sImage(sImage>background)),2); 
    [~,m]=max(TwoNormModel.mu);
    [~,m_min]=min(TwoNormModel.mu);
    thr = 2^(TwoNormModel.mu(m)-sigma*TwoNormModel.Sigma(:,:,m));

    if thr < 2^(TwoNormModel.mu(m_min))
        thr = 2^(TwoNormModel.mu(m_min)+sigma*TwoNormModel.Sigma(:,:,m_min));
    end
    
    
    % apply a gaussian filter to smooth image and pick up cells
    seeds = imclose(sImage,strel('disk',ceil(smallmorphs))); 
    seeds=imgaussfilt(seeds,smallmorphs*2.375,'FilterSize',round(cellsize*1.4)); %4.75 is best
    seeds = imregionalmax(seeds);
    seeds=seeds&sImage>thr;
    seeds=imdilate(seeds,strel('disk',ceil(smallmorphs*2)));
    
    
    % threshold and watershed image
    opencells = imopen(imclose(sImage,strel('disk',round(cellsize/20))),strel('disk',round(cellsize/10)));
    ThresImage=opencells>thr;
    L=watershed(imimposemin(-sImage,seeds));
    Nuclei=ThresImage&L;

    Nuclei=bwareaopen(Nuclei,ceil(cellsize));
    Nuclei = imfill(Nuclei,'holes');
    
    
    %%%%%%%%%%%%%% FIXING BROKEN CELLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if fixbrokencellsflag == 1
        LineSegments=ThresImage&~L;
        LineSegments=bwareaopen(LineSegments,ceil(smallmorphs*2.5));
        LineSegments=bwlabel(LineSegments,8);
        
        for Seg_i=1:max(LineSegments(:))
            
            FullCell=imreconstruct(LineSegments==Seg_i,ThresImage);
            BrokenCells=bwlabel(FullCell&~LineSegments>0);
            BrokenStats=regionprops(BrokenCells,'solidity','area');
            FullStats=regionprops(FullCell,'solidity','area');
            
            for Seg_k=1:max(BrokenCells(:))
                BrokenArea(Seg_k)=BrokenStats(Seg_k).Area;
                BrokenSolidity(Seg_k)=BrokenStats(Seg_k).Solidity;
            end
            FullArea=FullStats(1).Area;
            FullSolidity=FullStats(1).Solidity;
            if sum(FullSolidity>BrokenSolidity)
                
                Nuclei=Nuclei|LineSegments==Seg_i;
            end
            
            Nuclei=imopen(Nuclei,strel('disk',ceil(smallmorphs)));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % erode out a little bit of cells
    NucleiComplement = ~Nuclei;
    NucleiComplement = imdilate(NucleiComplement,strel('disk',ceil(smallmorphs/2)));
    Nuclei = Nuclei - NucleiComplement > 0;
    Nuclei = bwlabel(Nuclei,8);
    
    if max(Nuclei(:)>0)
        % shave off the excess for each cell
        pixel_values = regionprops(Nuclei, sImage, 'PixelValues');
        quant = cellfun(@(x) findmininvect(log2(x),method,quantile,fold), {pixel_values.PixelValues});
        
        if ~isnan(quant)
            
            quant = round(2.^quant);
            quantileMask = Nuclei;
            quantileMask(Nuclei > 0) = quant(Nuclei(Nuclei > 0));
            Nuclei_th = sImage > quantileMask;
            Nuclei = Nuclei_th & Nuclei > 0;
            Nuclei = imfill(Nuclei,'holes');
        end
    end

    SegmentedImage = bwlabel(Nuclei,8);
    if writeflag == 1
        imwrite(uint16(SegmentedImage),filename_out)
    end
    
       OutImage = double(SegmentedImage);