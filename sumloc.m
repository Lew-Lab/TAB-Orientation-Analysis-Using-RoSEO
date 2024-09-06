sumlocs = 0;
for idx = 1:5000
    tmp = loc_data{idx};
    sumlocs = sumlocs+size(tmp,1);
end