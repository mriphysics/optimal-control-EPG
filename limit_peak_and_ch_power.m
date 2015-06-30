function [c,ceq,gradc,gradceq] = limit_peak_and_ch_power(params,max_chpower,max_peak,Nch,frequencies)
% nonlin constraint to limit 
% peak RF and 
% total power CHANNEL-WISE! 
%
% max_chpower is hard constraints vector 1 x Nch
Nt = length(params)/Nch/2;
% power contraints
Cpow = zeros(1,Nt*Nch);
x = reshape(params(1:Nt*Nch),Nt,Nch); %real parts
y = reshape(params(Nt*Nch+1:end),Nt,Nch); % imaginary parts
Cpow = diag(frequencies)*(x.^2+y.^2);
sumPow = sum(Cpow,1);
% peak contraints
Cpeak = zeros(1,Nt*Nch);
for j=1:Nt*Nch
    Cpeak(j) = params(j)^2 + params(j+Nt*Nch)^2- max_peak.^2;
end

c = [sumPow-max_chpower, Cpeak];

ceq = [];
% gradient

% power
if nargout > 1
gradcPow = zeros(Nt*Nch*2,Nch);
j = 0;
for n = 1:Nt
    for ch = 1:Nch
        j = j+1;
        freq = rem(j,Nt);
        if freq == 0
            freq = Nt;
        end
        gradcPow(j,ch) = 2*params(j)*frequencies(freq);
        gradcPow(j+Nt*Nch,ch) = 2*params(j+Nt*Nch)*frequencies(freq);
    end
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