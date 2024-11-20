function valR = comp_valR_p(g,g_true,mask)

total = sum(mask(:));
valR = 0.5/total* sum((g(mask)-g_true(mask)).^2);

end