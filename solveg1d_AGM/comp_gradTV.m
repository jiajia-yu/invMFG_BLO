function gradg = comp_gradTV(g)

gradg = zeros(size(g));
diffg = diff(g);
gradg(1:end-1) = -sign(diffg);
gradg(2:end) = gradg(2:end) + sign(diffg);

end