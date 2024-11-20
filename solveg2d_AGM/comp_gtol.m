function gtol = comp_gtol(g)
gxx = g{1}; 
gxy = g{2};
gyy = g{3};
gtol = 1;

[nx,ny] = size(gxx);
for i = 1:nx
    for j = 1:ny
        D = eig([gxx(i,j),gxy(i,j); gxy(i,j),gyy(i,j)]);
        gtol = min(min(D),gtol);
    end
end

end