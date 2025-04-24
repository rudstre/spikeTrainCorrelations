function dataPairs = convertPairsFromSqFm3d(sqform)
for i = 1:size(sqform,3)
    data = squeeze(sqform(:,:,i));
    data(isnan(data)) = 0;
    pairs = squareform(data);
    pairs(~pairs) = nan;
    dataPairs(:,i) = pairs;
end