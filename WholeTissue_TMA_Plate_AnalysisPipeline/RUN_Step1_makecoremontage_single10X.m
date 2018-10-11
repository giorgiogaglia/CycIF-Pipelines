function RUN_Step1_makecoremontage_single10X(folderloc, prefix1, prefix2, tilesize, dim, background, radius10X)
tic
% 
% folderloc = 'D:\Giorgio\TMA_BRC15010_3_6Jun18\Cycle_1_10X\';
% prefix1 = {'A','B','C','D','E','F','G','H','I','J'};  % cols of cores
% prefix2 = linspace(1,15,15); % rows of cores
% tilesize = [2048 2048];
% 
% rows = 1;
% cols = 1;
%maxtile = rows*cols;
% dim = ['%02d'];
% cut = 100;
% background = 1000;
% radius10X = 700;


% cycle through the cores
for i1 = 1:length(prefix1)
    for i2 = 1:length(prefix2)

        img_name = [folderloc '\Cycle_1_10X' filesep 'Cycle1_10X' filesep prefix1{i1} ' - '  num2str(prefix2(i2),dim) '(wv DAPI - DAPI).tif'];
        try
            Cycle0_DAPI = imread(img_name);
            clear A
        catch
            disp(num2str(i1*100+i2))
            continue
        end

        opened_img = imclose(Cycle0_DAPI,strel('disk',50));
        core_outline = opened_img > background;
        core_outline = bwareaopen(core_outline,5000);
        stats = regionprops(core_outline,'Centroid','Area');
        area = cat(1,stats.Area);
        [val,ind] = max(area);
        if length(area)>1
            disp(['Core' prefix1{i1} ' - '  num2str(prefix2(i2),dim) 'might not have worked'])
        end
        row_core = max(1,round(stats(ind).Centroid(1,2)-radius10X)):min(round(stats(ind).Centroid(1,2)+radius10X),tilesize(1));
        col_core = max(1,round(stats(ind).Centroid(1,1)-radius10X)):min(round(stats(ind).Centroid(1,1)+radius10X),tilesize(2));
        realcore = Cycle0_DAPI(row_core,col_core);
        options.message = 'false';
        realcore = imresize(realcore, 2);
        saveastiff(uint16(realcore),[folderloc filesep 'RealCore_' prefix1{i1} ' - '  num2str(prefix2(i2),dim) '.tif'],options)

    end
end
        
  
    
    
    
    
    