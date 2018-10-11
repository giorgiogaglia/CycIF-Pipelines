function RUN_Step2_makecorestacks_update(folderloc,prefix1, prefix2, wavelengths, thr, cycles, dim, num_old_DAPI, rows, cols)  

tic

% % parameters that should not change
% tilesize = [2048 2048];
% wavelengths = {'wv Blue - FITC','wv Green - dsRed','wv Red - Cy5'};
% thr = max(tilesize);
% 
% folderloc = 'D:\Giorgio\TMA_BRC15010_3_6Jun18\';
% prefix1 = {'A','B','C','D','E','F','G','H','I','J'};  % cols of cores
% prefix2 = linspace(1,15,15); % rows of cores
% cycles = [1 2 3 4 5 6 7 8];  % cycles to analyse
% 
% % images per each core
% rows = 2;
% cols = 2;
% maxtile = rows*cols;
% dim = ['%02d'];
% 
% % if this is only an update modify this parameters
% num_old_DAPI = 0;
% num_old_Running = num_old_DAPI*4;





% cycle through the cores
for i1 = 1:length(prefix1)
    for i2 = 1:length(prefix2)
        
        final_stack = [folderloc filesep 'FullStacks\Core' prefix1{i1}  num2str(prefix2(i2),dim) '.tif'];
        DAPI_stack = [folderloc filesep 'DAPIStacks\DAPI_Core' prefix1{i1}  num2str(prefix2(i2),dim) '.tif'];
        DAPIStack = [];
        RunningStack = [];
        
        try
            DAPICycle0 = imread([folderloc filesep 'RealCore_' prefix1{i1} ' - '  num2str(prefix2(i2),dim) '.tif']);
        catch
            continue
        end

        try
            for i = 1:num_old_DAPI
                DAPIStack(:,:,i) = imread(DAPI_stack,'Index',i);
            end
            for i = 1:num_old_DAPI
                RunningStack(:,:,i) = imread(final_stack,'Index',i);
            end
        catch
            RunningStack = [];
            DAPIStack = [];
        end
        
        
        
        
        
        for i3 = 1:length(cycles) 
            clear mont_row mont_col tile_row tile_col flag
            % initialize new DAPI images with all zeros
            DAPI{i3} = zeros(size(DAPICycle0,1),size(DAPICycle0,2),'uint16');

            % load the DAPI images from the next cycle
            for j1 = 1:rows
                for j2 = 1:cols
                    tiles_img{j1,j2} = imread([folderloc filesep 'Cycle' num2str(cycles(i3)) '_20X' ...
                                        filesep prefix1{i1} ' - '  num2str(prefix2(i2),dim) ...
                                        '(fld ' num2str((j1-1)*2+j2) ' wv UV - DAPI).tif']);
                
                    c = normxcorr2(tiles_img{j1,j2},DAPICycle0);
                    [ypeak, xpeak] = find(c==max(c(:)));
                    tileshift_y(j1,j2) = ypeak-size(tiles_img{j1,j2},1);
                    tileshift_x(j1,j2) = xpeak-size(tiles_img{j1,j2},2);
                
                    if abs(tileshift_y(j1,j2)-(j1-1)*size(tiles_img{j1,j2},1)) < thr && abs(tileshift_x(j1,j2)-(j2-1)*size(tiles_img{j1,j2},2)) < thr
                        flag = 1;
                        mont_row{j1,j2} = 1+max(0,tileshift_y(j1,j2)):min(tileshift_y(j1,j2)+size(tiles_img{j1,j2},1),size(DAPI{i3},1));
                        mont_col{j1,j2} = 1+max(0,tileshift_x(j1,j2)):min(tileshift_x(j1,j2)+size(tiles_img{j1,j2},2),size(DAPI{i3},2));
                        tile_row{j1,j2} = (1-min(0,tileshift_y(j1,j2))):(size(tiles_img{j1,j2},1)-max(0,tileshift_y(j1,j2)+size(tiles_img{j1,j2},1)-size(DAPI{i3},1)));
                        tile_col{j1,j2} = (1-min(0,tileshift_x(j1,j2))):(size(tiles_img{j1,j2},2)-max(0,tileshift_x(j1,j2)+size(tiles_img{j1,j2},2)-size(DAPI{i3},2)));
                
                    
                        DAPI{i3}(mont_row{j1,j2},mont_col{j1,j2}) = tiles_img{j1,j2}(tile_row{j1,j2},tile_col{j1,j2});                                 
                
%                     figure(10+(j1-1)*3+j2)
%                     imshow(DAPI{i3},[])
                    else
                        abs(tileshift_y(j1,j2)-(j1-1)*size(tiles_img{j1,j2},1))
                        abs(tileshift_x(j1,j2)-(j2-1)*size(tiles_img{j1,j2},2))
                    end
                end
            end
        i3
            toc

        DAPIStack = cat(3,DAPIStack,DAPI{i3});
        RunningStack = cat(3,RunningStack,DAPI{i3});
        
        
        for wv1 = 1:length(wavelengths)
            temp_mont = zeros(size(DAPICycle0,1),size(DAPICycle0,2),'uint16');
            if flag == 1
            for j1 = 1:size(tile_col,1)
                for j2 = 1:size(tile_col,2)
                    temp_img{j1,j2} = imread([folderloc filesep 'Cycle' num2str(cycles(i3)) '_20X' ...
                                        filesep prefix1{i1} ' - '  num2str(prefix2(i2),dim) ...
                                        '(fld ' num2str((j1-1)*2+j2) ' '  wavelengths{wv1} ').tif']);
                    
                    temp_mont(mont_row{j1,j2},mont_col{j1,j2}) = temp_img{j1,j2}(tile_row{j1,j2},tile_col{j1,j2}); 
                end
            end
            end
                    
%             temp_mont_cropped = imcrop(temp_mont,rect);
%             if sum(abs(size(temp_mont_cropped)-size(DAPICycle0_cropped))) > 0 
%                 temp_mont_cropped = padarray(temp_mont_cropped,size(DAPICycle0_cropped)-size(temp_mont_cropped),0,'post');
%             end
            RunningStack = cat(3,RunningStack,temp_mont);
   
            clear  temp_mont_cropped temp_mont temp_img
        end

            
        end
        
% % %         

        
        
    final_stack = [folderloc filesep 'Core' prefix1{i1}  num2str(prefix2(i2),dim) '.tif'];
    DAPI_stack = [folderloc filesep 'DAPI_Core' prefix1{i1}  num2str(prefix2(i2),dim) '.tif'];
% % %     
    saveastiff(RunningStack,final_stack)
    saveastiff(DAPIStack,DAPI_stack)
    clear RunningStack DAPIStack DAPI_img DAPICycle0
        
    end
end


