function RUN_4_CellStateCompare_Path(basefolder, slides_folders, dates,marker_name, antibody_rounds, marker_cycle, HSF1round, canceround, antibody_type, max_rows, well_nums, sol_thresh) 
%This program chooses the fields with good staining then plots the
%frequency of HSF1 foci per core using the good fields 
% PLOT THE FREQUENCY OF HSF1 FOCI PER CORE

x_plot= [];
MovingMean = []; 
for folder = 1:length(slides_folders) 
openfile = [basefolder slides_folders{folder} dates slides_folders{folder} '_Results.mat'];
load(openfile)

openfile1 = 'Incorrect_ROIs.mat'; 
try
    load(openfile1)
    bad_fields = ROI_Fields(folder).NotCancer;    %Vector of fields that are not Cancer but have CD138 signal 
catch
    bad_fields = [];    %Vector of fields that are not Cancer but have CD138 signal 
end

Matrix = [];
Solidity = [];
temp = [];
temp_sol = [];
a = [];

for i1 = 1:max_rows
    for i3 = 1:length(well_nums)
for k = 1:length(Field)
    temp = [];
       try
        index = Field(k).Solidity(:,1)> sol_thresh;
       catch
           continue
       end 
       if ~ismember(k,bad_fields) %If field is not a bad field
        if length(index)==length(Plate{i1, well_nums(i3)}.Field(k).HSF1Foci_Sign) && length(index)>0 && length(Plate{i1, well_nums(i3)}.Field(k).MedianNucSign)==length(Plate{i1, well_nums(i3)}.Field(k).HSF1Foci_Sign) 
        temp_sol = Field(k).Solidity;
        % foci cells count
        temp(:,1) = Plate{i1, well_nums(i3)}.Field(k).HSF1Foci_Sign(:,1)./Plate{i1, well_nums(i3)}.Field(k).HSF1Foci_Sign(:,2);
        temp(:,2) = Plate{i1, well_nums(i3)}.Field(k).MedianNucSign(:,HSF1round);
        for i1=1:length(antibody_rounds)
            if isequal(antibody_type(i1), 'C')
                temp(:,i1+2) = Plate{i1, well_nums(i3)}.Field(k).MedianCytSign(:,antibody_rounds(i1));
            elseif isequal(antibody_type(i1), 'N')
                temp(:,i1+2) = Plate{i1, well_nums(i3)}.Field(k).MedianNucSign(:,antibody_rounds(i1));
            end 
        end 
        Matrix = [Matrix; temp];
        Solidity = [Solidity; temp_sol];
        else
            a = [a length(index)];
        end
       end 
end

Matrix = sortrows(Matrix,1);

figure('Name', [fols{folder} '_10'], 'NumberTitle', 'off')
for j1 = 2:length(marker_name)
% filter on Cancer value
cancol = find(canceround==antibody_rounds)+1; %Index of cancer antibody 
index1 = Matrix(:,cancol) > 2^7; %Filtering 
index2 = Solidity(:,marker_cycle(j1))>sol_thresh;

hold on 
submatrix = [Matrix(index1.*index2>0,1) Matrix(index1.*index2>0,j1-1)];
submatrix = submatrix(submatrix(:,2)>0,:);
try
    [x,idx] = datasample(submatrix(:,1),50000);
catch
    continue
end 
y = submatrix(idx,2);
subplot(1,7,j1-1)
plot(x,(y),'.')
if j1 == 2
    ylabel('Fluor AU')
elseif j1 == 4
    xlabel('% of HSF1 in Foci')
end
title(marker_name{j1})
hold on
bins = 201;
x = linspace(0,1,bins);
for ii = 2:bins
    if sum(submatrix(:,1)>=x(ii-1) & submatrix(:,1)<x(ii)) < 10
        x_plot(ii-1) = NaN;
        MovingMean(ii-1) = NaN;
    else
        x_plot(ii-1) = x(ii)-(x(2)/2);
        MovingMean(ii-1) = mean(submatrix(submatrix(:,1)>=x(ii-1) & submatrix(:,1)<x(ii),2));
    end
end
   
plot(x_plot,MovingMean,'r','LineWidth',3)
end

end
end 
end 

%% Plotting All the Data Together 

x_plot= [];
MovingMean = []; 

for folder = 1:length(slides_folders) 
openfile = [basefolder slides_folders{folder} dates slides_folders{folder} '_Results.mat'];
load(openfile)

try
    load(openfile1)
    bad_fields = ROI_Fields(folder).NotCancer;    %Vector of fields that are not Cancer but have CD138 signal 
catch
    bad_fields = [];    %Vector of fields that are not Cancer but have CD138 signal 
end

Matrix = [];
Solidity = [];
temp = [];
temp_sol = [];
a = [];

for i1 = 1:max_rows
    for i3 = 1:length(well_nums)
for k = 1:length(Field)
    temp = [];
       try
        index = Field(k).Solidity(:,1)> sol_thresh;
       catch
           continue
       end 
        if ~ismember(k,bad_fields) %If field is not a bad field
        if length(index)==length(Plate{i1, well_nums(i3)}.Field(k).HSF1Foci_Sign) && length(index)>0 && length(Plate{i1, well_nums(i3)}.Field(k).MedianNucSign)==length(Plate{i1, well_nums(i3)}.Field(k).HSF1Foci_Sign) 
        temp_sol = Plate{i1, well_nums(i3)}.Field(k).Solidity;
        % foci cells count
        temp(:,1) = Plate{i1, well_nums(i3)}.Field(k).HSF1Foci_Sign(:,1)./Plate{i1, well_nums(i3)}.Field(k).HSF1Foci_Sign(:,2);
        temp(:,2) = Plate{i1, well_nums(i3)}.Field(k).MedianNucSign(:,HSF1round);
        for i1=1:length(antibody_rounds)
            if isequal(antibody_type(i1), 'C')
                temp(:,i1+2) = Plate{i1, well_nums(i3)}.Field(k).MedianCytSign(:,antibody_rounds(i1));
            elseif isequal(antibody_type(i1), 'N')
                temp(:,i1+2) = Plate{i1, well_nums(i3)}.Field(k).MedianNucSign(:,antibody_rounds(i1));
            end 
        end 
        Matrix = [Matrix; temp];
        Solidity = [Solidity; temp_sol];
        else
            a = [a length(index)];
        end
       end 
end
    end 
end 
end 
Matrix = sortrows(Matrix,1);

%%

h=figure('Name', 'All Data_All_Antibodies', 'NumberTitle', 'off');
for j1 = 2:length(marker_name)
% filter on Cancer value
cancol = find(canceround==antibody_rounds)+1; %Index of cancer antibody 
index1 = Matrix(:,cancol) > 2^7; %Filtering 
index2 = Solidity(:,marker_cycle(j1))>sol_thresh;
    
submatrix = [Matrix(index1.*index2>0,1) Matrix(index1.*index2>0,j1)];
submatrix = submatrix(submatrix(:,2)>0,:);
try
    [x,idx] = datasample(submatrix(:,1),50000);
catch
    continue
end 
y = submatrix(idx,2);
subplot(1,length(marker_name)-1,j1-1)
plot(x,(y),'.')
if j1 == 2
    ylabel('Fluor AU')
elseif j1 == 4
    xlabel('% of HSF1 in Foci')
end
title(marker_name{j1})
hold on

bins = 201;
x1 = linspace(0,1,bins);
for ii = 2:bins
    if sum(submatrix(:,1)>=x1(ii-1) & submatrix(:,1)<x1(ii)) < 100
        x_plot(ii-1) = NaN;
        MovingMean(ii-1) = NaN;
    else
        x_plot(ii-1) = x1(ii)-(x1(2)/2);
        MovingMean(ii-1) = mean(submatrix(submatrix(:,1)>=x1(ii-1) & submatrix(:,1)<x1(ii),2));
    end
end
   
plot(x_plot,MovingMean,'r','LineWidth',3)

end
savefig(h, 'HSP_All_Plot_path_filt.fig')

clear MovingMean
end

 
