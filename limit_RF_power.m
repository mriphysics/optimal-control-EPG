function [c,ceq,gradc,gradceq] = limit_RF_power(params,max_power,Nch)
% nonlin constraint to limit RF power 
Nt = length(params)/Nch/2;
C = zeros(1,Nt*Nch);
for j=1:Nt*Nch
    C(j) = params(j)^2 + params(j+Nt*Nch)^2;
end
c = sum(C(:))-max_power;
ceq = [];
if nargout > 1
gradc = zeros(Nt*Nch*2,1);
for j = 1:Nt*Nch
    gradc(j,1) = 2*params(j);
    gradc(j+Nt*Nch,1) = 2*params(j+Nt*Nch);
end
    
gradceq = [];
end