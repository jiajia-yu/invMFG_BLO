function gradR = comp_gradR_p(g,g_true,mask)

gradR = cell(3,1);
total = sum(mask(:));
for i = 1:3
    gi = g{i};
    gtruei = g_true{i};
    gradR{i} = (gi-gtruei).*(mask/total);
end

end