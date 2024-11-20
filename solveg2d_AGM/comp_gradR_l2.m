function gradR = comp_gradR_l2(g)


gradR = cell(3,1);
nx = size(g{1},1);
ny = size(g{1},2);
for i=1:3
    gi = g{i};

    grad1 = (gi(2:end,:)-gi(1:end-1,:));
    grad2 = (gi(:,2:end)-gi(:,1:end-1));

    gradR1 = conv2(grad1,[-1;1],'full') * (nx/ny);
    gradR2 = conv2(grad2,[-1,1],'full') * (ny/nx);
    
    gradR{i} = gradR1 + gradR2;
end


end