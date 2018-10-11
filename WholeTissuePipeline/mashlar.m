CycleNStats = regionprops(SegImageN,'centroid','area','solidity','PixelList','PixelIdxList');
for i = 1:length(CycleNStats)
    Matrix_N(i,1) = i;
    Matrix_N(i,2) = CycleNStats(i).Centroid(1);
    Matrix_N(i,3) = CycleNStats(i).Centroid(2);
    Matrix_N(i,4) = CycleNStats(i).Area;
    Matrix_N(i,5) = CycleNStats(i).Solidity;
end
TrackN = zeros(size(SegImageN));
[IDX,D] = knnsearch(Matrix_1(:,2:3),Matrix_N(:,2:3));
for i = 1:length(Matrix_N(:,1))
TrackN(sub2ind(size(SegImage1),CycleNStats(i).PixelList(:,2),CycleNStats(i).PixelList(:,1)))=IDX(i);
end 