function res = min_empt(vec)
if isempty(vec)
    res = inf;
else
    res = min(vec);
end