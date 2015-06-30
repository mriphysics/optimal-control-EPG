function full = red2full(reduced,frequencies,Nch);
% transform vector of tip and phases from reduced form to full form
np = length(reduced)/(2*Nch);
xopt = reshape(reduced(1:np*Nch),np,Nch).'; %real parts
yopt = reshape(reduced(np*Nch+1:end),np,Nch).'; % imaginary parts
params0 = xopt +1i*yopt;
params1 = zeros(Nch,sum(frequencies));
sumfreq = cumsum(frequencies);
for j = (length(frequencies)):-1:2
    params1(:,sumfreq(j-1)+1:sumfreq(j)) = repmat(params0(:,j),[1 frequencies(j)]);
end
params1(:,1:frequencies(1)) = repmat(params0(:,1),[1 frequencies(1)]);
params1 = params1.';
full = [real(params1(:).') , imag(params1(:).')];