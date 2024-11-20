function valR = comp_valR_p(g,g_true,mask)

valR = 0;
total = sum(mask(:));
for i = 1:3
    gi = g{i};
    gtruei = g_true{i};
    valR = valR + 0.5/total* sum((gi(mask)-gtruei(mask)).^2);
end
