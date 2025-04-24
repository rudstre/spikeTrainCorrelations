function convertPairsFromSqFm(sqform)
sqform(isnan(sqform)) = 0;
pairs = squareform(sqform);
pairs(~pairs) = nan;