function gradR = comp_gradR_p(g,g_true,mask)

total = sum(mask(:));
gradR = (g-g_true).*(mask/total);

end