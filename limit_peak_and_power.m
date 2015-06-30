function [c,ceq,gradc,gradceq] = limit_peak_and_power(params,max_power,max_peak,Nch,frequencies)
% nonlin constraint to limit peak RF and total power 
Nt = length(params)/Nch/2;
% power contraints
Cpow = zeros(1,Nt*Nch);
for j=1:Nt*Nch
    freq = rem(j,Nt);
    if freq == 0
        freq = Nt;
    end
    Cpow(j) = frequencies(freq)*(params(j)^2 + params(j+Nt*Nch)^2);
end
% peak contraints
Cpeak = zeros(1,Nt*Nch);
for j=1:Nt*Nch
    Cpeak(j) = params(j)^2 + params(j+Nt*Nch)^2- max_peak.^2;
end

c = [sum(Cpow(:))-max_power, Cpeak];

ceq = [];
% gradient

% power
if nargout > 1
gradcPow = zeros(Nt*Nch*2,1);
for j = 1:Nt*Nch
    freq = rem(j,Nt);
    if freq == 0
        freq = Nt;
    end
    gradcPow(j,1) = 2*params(j)*frequencies(freq);
    gradcPow(j+Nt*Nch,1) = 2*params(j+Nt*Nch)*frequencies(freq);
end

% peak
gradcPeak = zeros(Nt*Nch*2,Nt*Nch);
for j = 1:Nt*Nch
    gradcPeak(j,j) = 2*params(j);
    gradcPeak(j+Nt*Nch,j) = 2*params(j+Nt*Nch);
end

gradc = [gradcPow, gradcPeak];
gradceq = [];
end