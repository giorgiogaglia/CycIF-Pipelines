function RUN_Step3_segmentstacks(microscope,filename,options) folderloc, prefix1,prefix2, DAPIslice, cellsize) 
% open the tiff file of the full stack images and segment them separately

fullstackfolder = [filename.analfolder filesep 'FullStacks' filesep];
trackedstackfolder = [filename.analfolder filesep 'TrackedStacks' filesep];
mkdir(trackedstackfolder)

    for i1 = length(filename.prefix1) 
        for i2 = 1:length(filename.prefix2)
            core = ['Core' filename.prefix1{i1} filename.midfix1 num2str(filename.prefix2(i2),filename.dim) filename.suffix];
            FileTif = [fullstackfolder core];
           
            disp(FileTif) %Display 
            
            clear Matrix_1
            
            track_stack = [trackedstackfolder 'TrackedCore' filename.prefix1{i1} num2str(filename.prefix2(i2),filename.dim) filename.suffix];
            
            if exist(track_stack, 'file') == 2
                disp(['Segmentation Error: ' core ' has already been tracked - skipping it!'])
                continue
            end
            
            try
                DAPIImage=imread(FileTif,'Index',filename.DAPIslices(1)); %Run Segmentation for 1st cycle 
            catch
                disp(['Segmentation Error: ' core ' not found'])
                continue
            end
          
            % run segmentation for cycle 1
            SegImage1 = CycIF_Segmentation_normfit_findminvect_v1(DAPIImage,0,options); %Test segmentation first 
            Cycle1Stats = regionprops(SegImage1,'centroid','area','solidity','PixelList','PixelIdxList');
            TrackImage_base = zeros(size(SegImage1));
            TrackImage = [];
            SegImage1 = SegImage1 > 0;
            
            if max(SegImage1(:)) > 0    %Checking that the field has cells 
                for i = 1:length(Cycle1Stats)
                    Matrix_1(i,1) = i;
                    Matrix_1(i,2) = Cycle1Stats(i).Centroid(1);
                    Matrix_1(i,3) = Cycle1Stats(i).Centroid(2);
                    Matrix_1(i,4) = Cycle1Stats(i).Area;
                    Matrix_1(i,5) = Cycle1Stats(i).Solidity;
                end
      
                for cycle=1:length(filename.DAPIslices)
                    clear Matrix_N CycleNStats centroid_N
                    DAPIImage=imread(FileTif,'Index',filename.DAPIslices(cycle));
                    if max(DAPIImage(:))>0
                        SegImageN = CycIF_Segmentation_normfit_findminvect_v1(DAPIImage,0,options);    %Segementing each cycle 
                    else
                        SegImageN = zeros(size(DAPIImage));
                        disp(['Segmentation Error: ' core ' DAPI cycle ' num2str((filename.DAPIslices+3)/4) ' image is empty'])
                    end
                    
                    CycleNStats = regionprops(SegImageN,'centroid','area','solidity','PixelList','PixelIdxList');
                    if ~isempty(CycleNStats)
                        for i = 1:length(CycleNStats)
                            Matrix_N(i,1) = i;
                            Matrix_N(i,2) = CycleNStats(i).Centroid(1);
                            Matrix_N(i,3) = CycleNStats(i).Centroid(2);
                            Matrix_N(i,4) = CycleNStats(i).Area;
                            Matrix_N(i,5) = CycleNStats(i).Solidity;
                        end
                        [IDX,D] = knnsearch(Matrix_1(:,2:3),Matrix_N(:,2:3));
                        
                        centroid_N = zeros(size(SegImage1));
                        centroid_N(sub2ind(size(SegImage1),round(Matrix_N(:,3)),round(Matrix_N(:,2)))) = 1;
                        centroid_N = imdilate(centroid_N,strel('disk',5));
                       
                        temp_id = [IDX zeros(size(IDX))+NaN];
                        
                        TrackN = zeros(size(SegImageN));
                        
                        test = zeros(size(TrackN));
                        
                        check = [];

                        for i = 1:length(Matrix_N(:,1))
                            pix_overlap = intersect(CycleNStats(i).PixelList,Cycle1Stats(IDX(i)).PixelList,'rows');
                            overlap = length(pix_overlap(:,1))/length(CycleNStats(i).PixelList(:,1));
                            temp_id(i,2)=overlap;
                            test(sub2ind(size(SegImage1),CycleNStats(i).PixelList(:,2),CycleNStats(i).PixelList(:,1)))= D(i);
                            if overlap > options.cellfracoverlap
                                TrackN(sub2ind(size(SegImage1),CycleNStats(i).PixelList(:,2),CycleNStats(i).PixelList(:,1)))=IDX(i);
                            else
                                check = [check i];
                                TrackN(CycleNStats(i).PixelIdxList)=0;
                            end
                        end
                        
                        TrackImage = cat(3,TrackImage,TrackN);
                    else
                        TrackN = zeros(size(SegImageN));
                    end
                end
                saveastiff(uint16(TrackImage),track_stack); % track cells based off cycle 1
                toc 
            end
        end
    end
    
