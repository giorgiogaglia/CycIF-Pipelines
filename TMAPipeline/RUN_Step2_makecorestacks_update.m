function RUN_Step2_makecorestacks_update(microscope,filename,options)

folderloc,prefix1, prefix2, wavelengths, thr, cycles, dim, num_old_DAPI, rows, cols)  

realcorefolder =  [filename.analfolder filesep 'RealCores' filesep];
fullstackfolder = [filename.analfolder filesep 'FullStacks' filesep];
dapistackfolder = [filename.analfolder filesep 'DapiStacks' filesep];
mkdir(fullstackfolder)
mkdir(dapistackfolder)

% cycle through the cores
for i1 = 1:length(filename.prefix1)
    for i2 = 1:length(filename.prefix2)
        
        final_stack = [filename.analfolder filesep 'FullStacks' filesep 'Core'      filename.prefix1{i1}  num2str(filename.prefix2(i2),filename.dim) filename.suffix ];
        DAPI_stack  = [filename.analfolder filesep 'DAPIStacks' filesep 'DAPI_Core' filename.prefix1{i1}  num2str(filename.prefix2(i2),filename.dim) filename.suffix ];
        DAPIStack = [];
        RunningStack = [];
        
        img_name = [realcorefolder filesep 'RealCore' filesep filename.prefix1{i1} filename.midfix1  num2str(filename.prefix2(i2),filename.dim) filename.suffix];
        try
            DAPICycle0 = imread(img_name);
        catch
            disp(['MakeCoreStacks Error: RealCore ' filename.prefix1{i1} filename.midfix1  num2str(filename.prefix2(i2),filename.dim) ' not found'])
            continue
        end

        try
            for i = 1:options.num_old_DAPI
                DAPIStack(:,:,i) = imread(DAPI_stack,'Index',i);
                RunningStack(:,:,i) = imread(final_stack,'Index',i);
            end
        catch
            RunningStack = [];
            DAPIStack = [];
        end
        
        
        
        
        
        for i3 = 1:length(filename.cycles) 
            clear mont_row mont_col tile_row tile_col flag
            % initialize new DAPI images with all zeros
            DAPI{i3} = zeros(size(DAPICycle0,1),size(DAPICycle0,2),'uint16');

            % load the DAPI images from the next cycle
            for j1 = 1:filename.rowspercore
                for j2 = 1:filename.colspercore
                    tiles_img{j1,j2} = imread([filename.datafolder filename.cycprefix num2str(filename.cycles(i3)) filesep ...
                                               filename.prefix1{i1} filename.midfix1  num2str(filename.prefix2(i2),filename.dim) ...
                                               filename.midfix2 num2str((j1-1)*2+j2) filename.channels{1} filename.suffix]);
                
                    c = normxcorr2(tiles_img{j1,j2},DAPICycle0);
                    [ypeak, xpeak] = find(c==max(c(:)));
                    tileshift_y(j1,j2) = ypeak-size(tiles_img{j1,j2},1);
                    tileshift_x(j1,j2) = xpeak-size(tiles_img{j1,j2},2);
                
                    if abs(tileshift_y(j1,j2)-(j1-1)*size(tiles_img{j1,j2},1)) < options.offsetthr && abs(tileshift_x(j1,j2)-(j2-1)*size(tiles_img{j1,j2},2)) < options.offsetthr
                        flag = 1;
                        mont_row{j1,j2} = 1+max(0,tileshift_y(j1,j2)):min(tileshift_y(j1,j2)+size(tiles_img{j1,j2},1),size(DAPI{i3},1));
                        mont_col{j1,j2} = 1+max(0,tileshift_x(j1,j2)):min(tileshift_x(j1,j2)+size(tiles_img{j1,j2},2),size(DAPI{i3},2));
                        tile_row{j1,j2} = (1-min(0,tileshift_y(j1,j2))):(size(tiles_img{j1,j2},1)-max(0,tileshift_y(j1,j2)+size(tiles_img{j1,j2},1)-size(DAPI{i3},1)));
                        tile_col{j1,j2} = (1-min(0,tileshift_x(j1,j2))):(size(tiles_img{j1,j2},2)-max(0,tileshift_x(j1,j2)+size(tiles_img{j1,j2},2)-size(DAPI{i3},2)));
                
                    
                        DAPI{i3}(mont_row{j1,j2},mont_col{j1,j2}) = tiles_img{j1,j2}(tile_row{j1,j2},tile_col{j1,j2});                                 
                
                    else
                        abs(tileshift_y(j1,j2)-(j1-1)*size(tiles_img{j1,j2},1))
                        abs(tileshift_x(j1,j2)-(j2-1)*size(tiles_img{j1,j2},2))
                    end
                end
            end

        DAPIStack = cat(3,DAPIStack,DAPI{i3});
        RunningStack = cat(3,RunningStack,DAPI{i3});
        
        
        for wv1 = 1:length(filename.channels)
            temp_mont = zeros(size(DAPICycle0,1),size(DAPICycle0,2),'uint16');
            if flag == 1
            for j1 = 1:size(tile_col,1)
                for j2 = 1:size(tile_col,2)
                    temp_img{j1,j2} = imread([filename.datafolder filename.cycprefix num2str(filename.cycles(i3)) filesep ...
                                              filename.prefix1{i1} filename.midfix1  num2str(filename.prefix2(i2),filename.dim) ...
                                              filename.midfix2 num2str((j1-1)*2+j2) filename.channels{wv1} filename.suffix]);
                    
                    temp_mont(mont_row{j1,j2},mont_col{j1,j2}) = temp_img{j1,j2}(tile_row{j1,j2},tile_col{j1,j2}); 
                end
            end
            end
                    
            RunningStack = cat(3,RunningStack,temp_mont);
   
            clear  temp_mont_cropped temp_mont temp_img
        end

            
        end
                
    final_stack = [fullstackfolder 'Core'      filename.prefix1{i1} filename.midfix1 num2str(filename.prefix2(i2),filename.dim) filename.suffix];
    DAPI_stack = [dapistackfolder  'DAPI_Core' filename.prefix1{i1} filename.midfix1 num2str(filename.prefix2(i2),filename.dim) filename.suffix];
    
    saveastiff(RunningStack,final_stack)
    saveastiff(DAPIStack,DAPI_stack)
    clear RunningStack DAPIStack DAPI_img DAPICycle0
        
    end
end


