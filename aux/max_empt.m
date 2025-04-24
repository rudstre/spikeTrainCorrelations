function res = max_empt(vec)
if isempty(vec)
    res = -inf;
else
    res = max(vec);
end

end