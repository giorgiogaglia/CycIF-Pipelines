function RUN_1_makecorestacks(infolderloc, outfolderloc,basicfolderloc, folders, maxfield, cycles, tilesize, prefix1, wavelengths, indx, indx1, max_rows, well_nums , well_lets, DAPIwv)

tic
rows = 1;
cols = 1;
thr = max(tilesize);

for folder = 1:length(folders)
    if maxfield(folder)<=100 && maxfield(folder)>10 %Changing dim to get correct file name
        dim = '%02d';
    elseif maxfield(folder)>1000
        dim = '%04d';
    elseif maxfield(folder)<=10
        dim = '%1d';
    else
        dim = '%03d';
    end
    prefix2 = linspace(1,maxfield(folder),maxfield(folder));
    
    for i3 = 1:length(cycles)
        for j1 = 1:4
            try
                DG{i3,j1} = double(imread([basicfolderloc num2str(cycles(i3)) '-dfp-basic.tif'],'Index',j1));  %
                FF{i3,j1} = double(imread([basicfolderloc num2str(cycles(i3)) '-ffp-basic.tif'],'Index',j1));  %
            catch
                [i3,j1]
                DG{i3,j1} = zeros(tilesize);    %No basic correction if it didn't work
                FF{i3,j1} = zeros(tilesize)+1;
            end
        end
    end
    
    
    % cycle through the cores
    for i1 = 1:max_rows
        for i3 = 1:length(well_nums)
            for i2 = 1:length(prefix2)
                %disp([folders{folder} ' ' num2str(i2)])
                try
                    DAPICycle0 = imread([infolderloc filesep folders{folder} filesep 'Cycle' num2str(cycles(1)) filesep prefix1{i1, i3} num2str(prefix2(i2),dim) ' ' DAPIwv ').tif']);    %
                    DAPICycle0 = uint16((double(DAPICycle0)-DG{1,1})./FF{1,1});
                catch
                    disp([folders{folder} ' ' num2str(i2)])
                    continue
                end
                RunningStack = [];
                DAPIStack = [];
                
                if prctile(DAPICycle0(:),97.5) < 100    %To skip fields that have no signal
                    continue
                elseif exist([outfolderloc filesep folders{folder} filesep 'DAPIStacks' filesep  well_lets{i1} well_nums(i3) 'DAPI_Field' num2str(prefix2(i2),dim1) 'test.tif'], 'file') == 2
                    continue
                end
                
                DAPIStack = cat(3,DAPIStack,DAPICycle0);    %Putting DAPICycle0 as 1st in stack- 2nd image is same as 1st in stack
                
                for i3 = 1:length(cycles)
                    clear mont_row mont_col tile_row tile_col flag
                    % initialize new DAPI images with all zeros
                    DAPI{i3} = zeros(size(DAPICycle0,1),size(DAPICycle0,2),'uint16');
                    
                    % load the DAPI images from the next cycle
                    for j1 = 1:rows
                        for j2 = 1:cols
                            tiles_img{j1,j2} = imread([infolderloc filesep folders{folder} filesep 'Cycle' num2str(cycles(1)) filesep prefix1{i1, i3} num2str(prefix2(i2),dim) ' ' DAPIwv ').tif']);   %
                            tiles_img{j1,j2} = uint16((double(tiles_img{j1,j2})-DG{i3,1})./FF{i3,1});
                            
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
                                
                            else
                                abs(tileshift_y(j1,j2)-(j1-1)*size(tiles_img{j1,j2},1))
                                abs(tileshift_x(j1,j2)-(j2-1)*size(tiles_img{j1,j2},2))
                            end
                        end
                    end
                    
                    toc
                    disp(['Field ' num2str(i2) ' Cycle ' num2str(i3)])
                    
                    DAPIStack = cat(3,DAPIStack,DAPI{i3});
                    RunningStack = cat(3,RunningStack,DAPI{i3});
                    
                    for wv1 = 1:length(wavelengths)
                        temp_mont = zeros(size(DAPICycle0,1),size(DAPICycle0,2),'uint16');
                        if flag == 1
                            if  i3 == indx1 && wv1 == indx %z-stacks
                                %disp('z-stacks')
                                for z1 = 1:5
                                    for j1 = 1:size(tile_col,1)
                                        for j2 = 1:size(tile_col,2)
                                            try
                                                temp_img{j1,j2} = imread([infolderloc filesep folders{folder} filesep 'Cycle' num2str(cycles(1)) filesep prefix1{i1, i3} num2str(prefix2(i2),dim) ' ' wavelengths{wv1} ' z ' num2str(z1, '%1d') ').tif'];
                                            catch
                                                temp_mont=[];
                                                continue
                                            end
                                            try
                                                temp_img{j1,j2} = uint16((double(temp_img{j1,j2})-DG{i3,wv1+1})./FF{i3,wv1+1});
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
                                for j1 = 1:size(tile_col,1)
                                    for j2 = 1:size(tile_col,2)
                                        try
                                            temp_img{j1,j2} = imread([infolderloc filesep folders{folder} filesep 'Cycle' num2str(cycles(1)) filesep prefix1{i1, i3} num2str(prefix2(i2),dim) ' '  wavelengths{wv1} ').tif']);   
                                        catch
                                            continue
                                        end
                                        try
                                            temp_img{j1,j2} = uint16((double(temp_img{j1,j2})-DG{i3,wv1+1})./FF{i3,wv1+1});
                                            temp_mont(mont_row{j1,j2},mont_col{j1,j2}) = temp_img{j1,j2}(tile_row{j1,j2},tile_col{j1,j2});
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
                    
                end
                dim1 = ['%04d'];
                
                final_stack = [outfolderloc filesep folders{folder} filesep 'FullStacks' filesep well_lets{i1} well_nums(i3) '_Field' num2str(prefix2(i2),dim1) 'test.tif']; 
                DAPI_stack  = [outfolderloc filesep folders{folder} filesep 'DAPIStacks' filesep  well_lets{i1} well_nums(i3) 'DAPI_Field' num2str(prefix2(i2),dim1) 'test.tif'];    
                
                saveastiff(RunningStack,final_stack);
                saveastiff(DAPIStack,DAPI_stack);
                clear RunningStack DAPIStack DAPI_img DAPICycle0
                
            end
        end
    end
end
end
