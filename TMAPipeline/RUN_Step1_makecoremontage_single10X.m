function RUN_Step1_makecoremontage_single10X(microscope,filename,options)

realcorefolder = [filename.analfolder filesep 'RealCores' filesep];
mkdir(realcorefolder);

% cycle through the cores
for i1 = 1:length(filename.prefix1)
    for i2 = 1:length(filename.prefix2)

        img_name = [filename.datafolder filename.cycprefix filesep filename.prefix1{i1} ' - '  num2str(filename.prefix2(i2),filename.dim) filename.channels{1} filename.suffix ];
        try
            Cycle0_DAPI = imread(img_name);
        catch
            disp(['RealCore Error: Core ' filename.prefix1{i1} ' - '  num2str(filename.prefix2(i2),filename.dim) ' not found'])
            continue
        end

        opened_img = imclose(Cycle0_DAPI,strel('disk',round(options.cellsize*options.MagDiff10x)));
        core_outline = opened_img > options.background;
        core_outline = bwareaopen(core_outline,round((microscope.tilesize/3)^2));
        stats = regionprops(core_outline,'Centroid','Area');
        area = cat(1,stats.Area);
        [~,ind] = max(area);
        if length(area)>1
            disp(['RealCore Error: Core ' filename.prefix1{i1} ' - '  num2str(filename.prefix2(i2),filename.dim) 'might not have worked'])
        end
        row_core = max(1,round(stats(ind).Centroid(1,2)-options.radius10X)):min(round(stats(ind).Centroid(1,2)+options.radius10X),microscope.tilesize(1));
        col_core = max(1,round(stats(ind).Centroid(1,1)-options.radius10X)):min(round(stats(ind).Centroid(1,1)+options.radius10X),microscope.tilesize(2));
        realcore = Cycle0_DAPI(row_core,col_core);
        optionstiff.message = 'false';
        realcore = imresize(realcore, options.MagDiff10x);
        saveastiff(uint16(realcore),[realcorefolder filesep 'RealCore' filesep filename.prefix1{i1} ' - '  num2str(filename.prefix2(i2),filename.dim) filename.suffix],optionstiff)
        disp(['RealCore ' filename.prefix1{i1} ' - '  num2str(filename.prefix2(i2),filename.dim) ' Done'])

    end
end
        
  
    
    
    
    
    