function gradR = comp_gradR_l2(g)

nx = size(g,1);

grad1 = (g(2:end,:)-g(1:end-1,:));
gradR = conv2(grad1,[-1;1],'full') * nx;
    


end