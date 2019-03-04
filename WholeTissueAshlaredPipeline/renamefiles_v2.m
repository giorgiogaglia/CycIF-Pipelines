
basefolder = 'W:\Analysis\';
slideNames = {'AJ0188_P1', 'AJ0188_P2', 'AJ0188_P3', 'AJ0188_P4', 'AJ0188_P5','AJ0188_P6', 'AJ0189_P1', 'AJ0189_P2', 'AJ0189_P3', 'AJ0189_P4', 'AJ0189_P5', 'AJ0189_P6',...
    'AJ0205_P1', 'AJ0205_P2', 'AJ0205_P3', 'AJ0205_P4', 'AJ0205_P5','AJ0205_P6'}; 
for i = 1:length(slideNames)
    for i1 = 1:10
        for i2 = 1:10
            outputfilename = [basefolder slideNames{i} '\FullStacks\Field_' num2str(i1 , filename.dim) '_' num2str(i2, filename.dim) '.tif']; 
            
            rename = [basefolder slideNames{i} '\FullStacks\' slideNames{i} '_Field_' num2str(i1 , filename.dim) '_' num2str(i2, filename.dim) '.tif']; 
            disp(rename)
            try
            movefile(outputfilename, rename) 
            catch
                disp('Not renamed') 
                continue 
            end 
        end 
    end 
end 