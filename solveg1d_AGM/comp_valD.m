function valD = comp_valD(dd,rho_cell,m_cell,rhotilde_cell,mtilde_cell)
% valD = comp_valD(dd,rho_cell,m_cell,rhotilde_cell,mtilde_cell)
% compute the value of upper level objective D
% the value is normalized by meshgrid in this function 
%       (dd is the volume of one infinitestimal element)

fun = @(a,atilde) sum((a-atilde).^2,'all');
valrho = cellfun(fun,rho_cell,rhotilde_cell);
valm = cellfun(fun,m_cell,mtilde_cell);

valD = 0.5*dd*(sum(valrho(:)) + sum(valm(:)));

end