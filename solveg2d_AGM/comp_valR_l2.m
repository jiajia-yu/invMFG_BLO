function valR = comp_valR_l2(g)

valR = 0;
nx = size(g{1},1);
ny = size(g{1},2);
for i=1:3
    gi = g{i};

    grad1 = (gi(2:end,:)-gi(1:end-1,:));
    grad2 = (gi(:,2:end)-gi(:,1:end-1));
    valR = valR + 0.5*((nx/ny)*sum(grad1.^2,'all') ...
                     + (ny/nx)*sum(grad2.^2,'all'));
end

end