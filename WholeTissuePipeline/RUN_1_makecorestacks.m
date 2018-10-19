function RUN_1_makecorestacks(master_folderIN, basefolder,basicfolderloc, slides_folders, maxfields, cycles, tilesize, prefix1, wavelengths, z_wave, z_cycle, max_rows, well_nums , well_lets, DAPIwv)

tic
rows = 1;
cols = 1;
thr = max(tilesize);
dim1 = '%04d';
for folder = 1:length(slides_folders)
    if maxfields(folder)<=100 && maxfields(folder)>10 %Changing dim to get correct file name
        dim = '%02d';
    elseif maxfields(folder)>1000
        dim = '%04d';
    elseif maxfields(folder)<=10
        dim = '%1d';
    else
        dim = '%03d';
    end
    prefix2 = linspace(1,maxfields(folder),maxfields(folder));
    
    for x2 = 1:length(cycles)
        for x1 = 1:4
            try
                DG{x2,x1} = double(imread([basicfolderloc num2str(cycles(x2)) '-dfp-basic.tif'],'Index',x1));  %
                FF{x2,x1} = double(imread([basicfolderloc num2str(cycles(x2)) '-ffp-basic.tif'],'Index',x1));  %
            catch
                disp([x2,x1])
                DG{x2,x1} = zeros(tilesize);    %No basic correction if it didn't work
                FF{x2,x1} = zeros(tilesize)+1;
            end
        end
    end
    
    
    % cycle through the cores
    for i1 = 1:max_rows
        for i3 = 1:length(well_nums)
            for i2 = 1:length(prefix2)
                %disp([slides_folders{folder} ' ' num2str(i2)])
                try
                    DAPICycle0 = imread([master_folderIN filesep slides_folders{folder} filesep 'Cycle_' num2str(cycles(1)) filesep prefix1{i1, i3} num2str(prefix2(i2),dim) ' ' DAPIwv ').tif']);    %
                    DAPICycle0 = uint16((double(DAPICycle0)-DG{1,1})./FF{1,1});
                catch
                    disp([slides_folders{folder} ' ' num2str(i2)])
                    continue
                end
                RunningStack = [];
                DAPIStack = [];
                
                if prctile(DAPICycle0(:),97.5) < 100    %To skip fields that have no signal
                    continue
                elseif exist([basefolder filesep slides_folders{folder} filesep 'DAPIStacks' filesep  well_lets{i1} well_nums(i3) 'DAPI_Field' num2str(prefix2(i2),dim1) 'test.tif'], 'file') == 2
                    continue
                end
                
                DAPIStack = cat(3,DAPIStack,DAPICycle0);    %Putting DAPICycle0 as 1st in stack- 2nd image is same as 1st in stack
                
                for c1 = 1:length(cycles)
                    clear mont_row mont_col tile_row tile_col flag
                    % initialize new DAPI images with all zeros
                    DAPI{c1} = zeros(size(DAPICycle0,1),size(DAPICycle0,2),'uint16');
                    
                    % load the DAPI images from the next cycle
                   for g1 = 1:rows
                        for g2 = 1:cols
                            tiles_img{g1,g2} = imread([master_folderIN filesep slides_folders{folder} filesep 'Cycle_' num2str(cycles(c1)) filesep prefix1{i1, i3} num2str(prefix2(i2),dim) ' ' DAPIwv ').tif']);   %
                            tiles_img{g1,g2} = uint16((double(tiles_img{g1,g2})-DG{c1,1})./FF{c1,1});
                            
                            c = normxcorr2(tiles_img{g1,g2},DAPICycle0);
                            [ypeak, xpeak] = find(c==max(c(:)));
                            tileshift_y(g1,g2) = ypeak-size(tiles_img{g1,g2},1);
                            tileshift_x(g1,g2) = xpeak-size(tiles_img{g1,g2},2);
                            
                            if abs(tileshift_y(g1,g2)-(g1-1)*size(tiles_img{g1,g2},1)) < thr && abs(tileshift_x(g1,g2)-(g2-1)*size(tiles_img{g1,g2},2)) < thr
                                flag = 1;
                                mont_row{g1,g2} = 1+max(0,tileshift_y(g1,g2)):min(tileshift_y(g1,g2)+size(tiles_img{g1,g2},1),size(DAPI{c1},1));
                                mont_col{g1,g2} = 1+max(0,tileshift_x(g1,g2)):min(tileshift_x(g1,g2)+size(tiles_img{g1,g2},2),size(DAPI{c1},2));
                                tile_row{g1,g2} = (1-min(0,tileshift_y(g1,g2))):(size(tiles_img{g1,g2},1)-max(0,tileshift_y(g1,g2)+size(tiles_img{g1,g2},1)-size(DAPI{c1},1)));
                                tile_col{g1,g2} = (1-min(0,tileshift_x(g1,g2))):(size(tiles_img{g1,g2},2)-max(0,tileshift_x(g1,g2)+size(tiles_img{g1,g2},2)-size(DAPI{c1},2)));
                                
                                
                                DAPI{c1}(mont_row{g1,g2},mont_col{g1,g2}) = tiles_img{g1,g2}(tile_row{g1,g2},tile_col{g1,g2});
                                
                            else
                                abs(tileshift_y(g1,g2)-(g1-1)*size(tiles_img{g1,g2},1))
                                abs(tileshift_x(g1,g2)-(g2-1)*size(tiles_img{g1,g2},2))
                            end
                        end
                    end
                    
                    toc
                    disp(['Field ' num2str(i2) ' Cycle ' num2str(cycles(c1))])
                    
                    DAPIStack = cat(3,DAPIStack,DAPI{c1});
                    RunningStack = cat(3,RunningStack,DAPI{c1});
                    
                    for wv1 = 1:length(wavelengths)
                        temp_mont = zeros(size(DAPICycle0,1),size(DAPICycle0,2),'uint16');

                            if  cycles(c1) == z_cycle && isequal(z_wave, wavelengths{wv1}) %z-stacks
                                %disp('z-stacks')
                                for z1 = 1:5
                                    for j1 = 1:size(tile_col,1)
                                        for j2 = 1:size(tile_col,2)
                                            try
                                                temp_img{j1,j2} = imread([master_folderIN filesep slides_folders{folder} filesep 'Cycle_' num2str(cycles(c1)) filesep prefix1{i1, i3} num2str(prefix2(i2),dim) ' ' wavelengths{wv1} ' z ' num2str(z1, '%1d') ').tif']);
                                            catch
                                                temp_mont=[];
                                                continue
                                            end
                                            try
                                                temp_img{j1,j2} = uint16((double(temp_img{j1,j2})-DG{c1,wv1+1})./FF{c1,wv1+1});
                                                temp_mont(mont_row{j1,j2},mont_col{j1,j2}) = temp_img{j1,j2}(tile_row{j1,j2},tile_col{j1,j2});
                                            catch
                                                continue
                                            end
                                        end
                                    end
                                    RunningStack = cat(3,RunningStack,temp_mont);
                                    clear  temp_mont_cropped temp_mont temp_img
                                end
                            else
                               for y1 = 1:size(tile_col,1)
                                    for y2 = 1:size(tile_col,2)
                                        try
                                            temp_img{y1,y2} = imread([master_folderIN filesep slides_folders{folder} filesep 'Cycle_' num2str(cycles(c1)) filesep prefix1{i1, i3} num2str(prefix2(i2),dim) ' '  wavelengths{wv1} ').tif']);  
                                        catch
                                            continue
                                        end
                                        try
                                            temp_img{y1,y2} = uint16((double(temp_img{y1,y2})-DG{c1,wv1+1})./FF{c1,wv1+1});
                                            temp_mont(mont_row{y1,y2},mont_col{y1,y2}) = temp_img{y1,y2}(tile_row{y1,y2},tile_col{y1,y2});
                                        catch
                                            continue
                                        end
                                    end
                                end
                                RunningStack = cat(3,RunningStack,temp_mont);
                                clear  temp_mont_cropped temp_mont temp_img
                            end
                        
                    end
                    
                end
                
                
                final_stack = [basefolder filesep slides_folders{folder} filesep 'FullStacks' filesep well_lets{i1} num2str(well_nums(i3)) '_Field' num2str(prefix2(i2),dim1) '.tif']; 
                DAPI_stack  = [basefolder filesep slides_folders{folder} filesep 'DAPIStacks' filesep  well_lets{i1} num2str(well_nums(i3)) 'DAPI_Field' num2str(prefix2(i2),dim1) '.tif'];    
                
                saveastiff(RunningStack,final_stack);
                saveastiff(DAPIStack,DAPI_stack);
                clear RunningStack DAPIStack DAPI_img DAPICycle0
                
            end
        end
    end
end
end

