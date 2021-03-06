function RUN_2_segmentstacks_Validation(basefolder, folders, maxfield, prefix1,DAPIslice, cellsize) 
% open the tiff file of the DAPI images and segment them separately
filenameout=0; 
sigma = 1; 
writeflag = 0;
fixbrokencellsflag = 0;
options = [writeflag,fixbrokencellsflag,sigma,cellsize]; 
imagefolder = 'FullStacks\';
count2 = 0;
dim = ['%04d']; 

for folder = 1:length(folders)
    prefix2 = linspace(1,maxfield(folder),maxfield(folder));
    
    for i1 = 1:length(prefix1)
        for i2 = 1:length(prefix2)
            core = ['Field' num2str(prefix2(i2),dim) 'test.tif'];
            FileTif = [basefolder folders{folder} filesep imagefolder core];
           
            disp(FileTif) %Display 
            
            clear Matrix_1
            
            track_stack = [basefolder folders{folder} filesep 'TrackedImages\TrackedField' num2str(prefix2(i2),dim) '.tif'];
            
            if exist(track_stack, 'file') == 2
                continue
            end
            
            try
                DAPIImage=imread(FileTif,'Index',DAPIslice(1)); %Run Segmentation for HSF1 cycle- Change DAPIslice to correct one for HSF1 cycle or 1st cycle 
            catch
                continue
            end
            
            count2 = count2 + 1;
            
            % run segmentation for cycle 1
            
            filename_out = [basefolder folders{folder} filesep 'SegImages\Seg' num2str(prefix2(i2),dim) '.tif'];
            SegImage1 = CycIF_Segmentation_normfit_opencellsizeover5_medianpercell_v1(DAPIImage,filenameout,options); %Test segmentation first 
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
                
                for cycle=1:length(DAPIslice)
                    clear Matrix_N CycleNStats centroid_N
                    DAPIImage=imread(FileTif,'Index',DAPIslice(cycle));
                    if max(DAPIImage(:))>0
                        SegImageN = CycIF_Segmentation_normfit_opencellsizeover5_medianpercell_v1(DAPIImage,filenameout,options);    %Segementing each cycle 
                    else
                        SegImageN = zeros(size(DAPIImage));
                        disp('This image is empty')
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
                        count = 0;
                        for i = 1:length(Matrix_N(:,1))
                            pix_overlap = intersect(CycleNStats(i).PixelList,Cycle1Stats(IDX(i)).PixelList,'rows');
                            overlap = length(pix_overlap(:,1))/length(CycleNStats(i).PixelList(:,1));
                            temp_id(i,2)=overlap;
                            test(sub2ind(size(SegImage1),CycleNStats(i).PixelList(:,2),CycleNStats(i).PixelList(:,1)))= D(i);
                            if overlap > 0.5
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
            end
        end
    end
    
end
end 
