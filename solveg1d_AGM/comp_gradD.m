function [gradrhoD_cell,gradmD_cell] ...
    = comp_gradD(rho_cell,m_cell,rhotilde_cell,mtilde_cell)
% [gradrhoD_cell,gradmD_cell] ...
%     = comp_gradD(rho_cell,m_cell,rhotilde_cell,mtilde_cell)
% compute the gradient of upper level objective D
% with respect to density rho and flux m
% the value is not normalized by meshgrid in this function

fun = @(a,atilde) a-atilde;
gradrhoD_cell = cellfun(fun,rho_cell,rhotilde_cell,'UniformOutput',false);
gradmD_cell = cellfun(fun,m_cell,mtilde_cell,'UniformOutput',false);

end

