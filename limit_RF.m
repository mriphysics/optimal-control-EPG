function [c,ceq,gradc,gradceq] = limit_RF(params,max_val,Nch)
% nonlin constraint to limit RF 
Nt = length(params)/Nch/2;
c = zeros(1,Nt*Nch);
for j=1:Nt*Nch
    c(j) = params(j)^2 + params(j+Nt*Nch)^2- max_val.^2;
end
ceq = [];
if nargout > 1
gradc = zeros(Nt*Nch*2,Nt*Nch);
for j = 1:Nt*Nch
    gradc(j,j) = 2*params(j);
    gradc(j+Nt*Nch,j) = 2*params(j+Nt*Nch);
end
    
gradceq = [];
end