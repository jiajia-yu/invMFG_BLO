function g = proj_g(g,tol)

gxx = g{1}; 
gxy = g{2};
gyy = g{3};

[nx,ny] = size(gxx);
for i = 1:nx
    for j = 1:ny
        [V,D] = eig([gxx(i,j),gxy(i,j); gxy(i,j),gyy(i,j)]);
        D = diag(max(diag(D),tol));
        V = V*D*V';
        gxx(i,j) = V(1,1);
        gxy(i,j) = V(1,2);
        gyy(i,j) = V(2,2);
    end
end
g = {gxx;gxy;gyy};

end