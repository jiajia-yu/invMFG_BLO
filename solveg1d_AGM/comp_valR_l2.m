function valR = comp_valR_l2(g)

nx = size(g,1);
grad1 = (g(2:end)-g(1:end-1));
valR = (0.5*nx)*sum(grad1.^2,'all') ;

end