function b = proj_b(b)
% b = proj_b(b)
% solves the projection such that sum(b(:))=0

b = b-mean(b(:));

end