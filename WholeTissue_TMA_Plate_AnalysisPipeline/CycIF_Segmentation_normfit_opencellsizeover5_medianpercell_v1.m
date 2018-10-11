function OutImage = CycIF_Segmentation_normfit_opencellsizeover5_medianpercell_v1(image_in,filename_out,options)
%For Myeloma
% clear all
% basefolder = 'D:\Myeloma_HSF1_2018-07-27\Myeloma_Segmentation\';
% imagefolder = 'Test_files\';
% %core = ['Field' num2str(prefix2(i2),dim) 'test.tif'];
% core = ['Field0002test.tif','Field0018test.tif','Field0032test.tif','Field0033test.tif','Field0039test.tif','Field0052test.tif',];
% prefix2= {'0002', '0018' '0032' '0033' '0039' '0052'};
% %prefix2= {'051'};
% maxcycle = 4;
% %DAPIslice = (linspace(1,maxcycle,maxcycle)-1)*4+1; %Normal
% 
% %DAPIslice= [1, 3, 7, 11];
% 
% sigma = 1;
% cellsize = 23.5;  %Change first!!!
% writeflag = 1;
% fixbrokencellsflag = 0;
% output=[basefolder 'Seg_results\']; 
% testvar= 'T7_16_r'; 
% filename_out = [basefolder imagefolder 'test2.tif'];

%for i2=1:length(prefix2)
%     core = ['Field' prefix2{i2} 'test.tif'];
%     if i2 == 2
%         DAPIslice = [1, 3, 7, 11]; 
%     else
%         DAPIslice = [1, 6, 10, 14];
%     end
    
%     FileTif= [basefolder imagefolder core];
%     DAPIImage=imread(FileTif,'Index',DAPIslice(1));
%     image_in = DAPIImage;
%     
    % load image to be segmented if not loaded already
writeflag = options(1);
fixbrokencellsflag = options(2);
sigma = options(3);
cellsize = options(4); 
    if ischar(image_in)
        RawImage = double(imread(filename_in));
    else
        RawImage = double(image_in);
    end
    
    % figure(1)
    % imshow(RawImage,[0 200])
    
    %       *IMAGE PROCESSING GOES HERE.*
    
    % --------------------------- FINDING NUCLEI ------------------------------
    % takes the DNA image and gives back a "Nuclei" image
    
    
    % first do local background subtraction
    coImage=imopen(imclose(RawImage,strel('disk',round(cellsize/20))),strel('disk',round(cellsize*1)));
    sImage = RawImage-coImage;
    
    % figure(2)
    % imshow(RawImage,[])
    % figure(12)
    % imshow(coImage,[])
    % figure(13)
    % imshow(sImage,[])
    
    % then determine the threshold
    TwoNormModel = fitgmdist(log2(sImage(sImage>150)),2); 
    [y,m]=max(TwoNormModel.mu);
    thr = 2^(TwoNormModel.mu(m)-sigma*TwoNormModel.Sigma(:,:,m));
    % [n,h]=ksdensity(log2(sImage(sImage>8)));
    % [pks,locs] = findpeaks(n);
    % if length(locs)~=3
    %     thr = prctile(sImage(~imclose(sImage<0,strel('disk',5))),80)
    % else
    % %     plot(h(locs(length(locs)-1):locs(length(locs))),n(locs(length(locs)-1):locs(length(locs))));
    %     [x,thr_loc] = min(n(locs(length(locs)-1):locs(length(locs))));
    %     thr = h(locs(length(locs)-1)+thr_loc)+1;
    %     thr = 2^(thr);
    % end
    
    % figure(fignum)
    % plot(h,n,log2([thr thr]),[0 0.3])
    % hold on
    % plot(h,normpdf(h,TwoNormModel.mu(1),TwoNormModel.Sigma(:,:,1)),h,normpdf(h,TwoNormModel.mu(2),TwoNormModel.Sigma(:,:,2)))
    % hold on
    % plot(h(locs(length(locs)-1):locs(length(locs))),n(locs(length(locs)-1):locs(length(locs))));
    
    % apply a gaussian filter to smooth image and pick up cells
%   Test7 code:
    seeds = imclose(sImage,strel('disk',2)); 
%     figure(3)
%     imshow(seeds,[])
    seeds=imgaussfilt(seeds,4.75,'FilterSize',round(cellsize*1.4)); %4.75 is best
%     figure(4)
%     imshow(seeds,[])
    seeds = imregionalmax(seeds);
%     figure(5)
%     imshow(max(sImage,seeds*60000),[])
    seeds=seeds&sImage>thr;
    seeds=imdilate(seeds,strel('disk',4));
    
    % % optional parts to respace the seeds, probably not useful
    % stats = regionprops(seeds,'centroid');
    % centroids = round(cat(1, stats.Centroid));
    % if ~isempty(centroids)
    %     seeds = zeros(size(seeds));
    %     seeds(sub2ind(size(seeds), centroids(:,2), centroids(:,1))) = 1;
    %     seeds=imdilate(seeds,strel('disk',4));
    % end
    
    % threshold and watershed image
    % sImage2 = imclose(sImage,strel('disk',round(1)));
    opencells = imopen(imclose(sImage,strel('disk',round(cellsize/20))),strel('disk',round(cellsize/6)));
    % ThresImage=sImage>thr;
    ThresImage=opencells>thr;
    %ThresImage = imfill(ThresImage,'holes');%uncommented
    L=watershed(imimposemin(-sImage,seeds));
    Nuclei=ThresImage&L;
    
%     figure(6)
%     imshow(Nuclei,[])
    
    Nuclei=bwareaopen(Nuclei,cellsize/20);
    Nuclei = imfill(Nuclei,'holes');
    
    
    %%%%%%%%%%%%%% FIXING BROKEN CELLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if fixbrokencellsflag == 1
        LineSegments=ThresImage&~L;
        LineSegments=bwareaopen(LineSegments,5);
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
            
            Nuclei=imopen(Nuclei,strel('disk',2));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     figure(7)
%     imshow(Nuclei,[])
    
    % erode out a little bit of cells
    NucleiComplement = ~Nuclei;
    NucleiComplement = imdilate(NucleiComplement,strel('disk',1));
    Nuclei = Nuclei - NucleiComplement > 0;
    Nuclei = bwlabel(Nuclei,8);
    
    if max(Nuclei(:)>0)
        % shave off the excess for each cell
        pixel_values = regionprops(Nuclei, sImage, 'PixelValues');
        quant = cellfun(@(x) findmininvect(log2(x),3,0.3,3), {pixel_values.PixelValues});
        
        if ~isnan(quant)
            
            quant = round(2.^quant);
            quantileMask = Nuclei;
            quantileMask(Nuclei > 0) = quant(Nuclei(Nuclei > 0));
%             figure(8); imshow(quantileMask,[])
            Nuclei_th = sImage > quantileMask;
%             figure(9); imshow(Nuclei_th,[])
            Nuclei = Nuclei_th & Nuclei > 0;
%             figure(10); imshow(Nuclei,[])
            Nuclei = imfill(Nuclei,'holes');
        end
    end
%     figure(11)
%     imshow(Nuclei,[])

    SegmentedImage = bwlabel(Nuclei,8);
    if writeflag == 1
        imwrite(uint16(SegmentedImage),filename_out)
    end
    
       OutImage = double(SegmentedImage);
%     imwrite(uint16(SegmentedImage),[output testvar 'Seg_' prefix2{i2} '.tif'])
%     
%     overlay = cat(3,uint16(sImage),uint16(Nuclei));
%     %overlay = cat(3,overlay,(zeros(size(RawImage))));
%     overlay = cat(3,overlay, max(sImage,seeds*60000));
% % %     figure(12)
% %     imshow(SegmentedImage,[])
%     
%    imwrite(overlay,[output testvar 'Seg_overlay_' prefix2{i2} '.tif'])
%end
%disp('Done!')