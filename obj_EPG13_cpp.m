function [obj,gr] = obj_EPG13_cpp(params,ESP,T1,T2,c,B1,target,frequencies,klim)
% call armadillo-c++ equivalent of obj_EPG13
%% construct angles -> params2
[nt] = size(target,1);
[Ns, Nch] = size(B1);
np = length(params)/(2*Nch);
x = reshape(params(1:np*Nch),np,Nch).'; %real parts
y = reshape(params(np*Nch+1:end),np,Nch).'; % imaginary parts
params0 = x +1i*y;
params1 = zeros(Nch,sum(frequencies));
sumfreq = cumsum(frequencies);
for j = (length(frequencies)):-1:2
    params1(:,sumfreq(j-1)+1:sumfreq(j)) = repmat(params0(:,j),[1 frequencies(j)]);
end
params1(:,1:frequencies(1)) = repmat(params0(:,1),[1 frequencies(1)]);
params1 = params1.';
params2 = [real(params1(:).') , imag(params1(:).')];

%% write data to file

Nt = nt-1;

header = [nargout;Nt;Ns;Nch;klim];
vector = [header;params2(:);ESP;T1;T2;c(:);real(B1(:));imag(B1(:));real(target(:));imag(target(:))];
save input_data.out vector -ASCII
%% call C++ function
system('./obj_EPG13');
%% read output of C++ function
data = load('output_data.out');
obj = data(1);
GR = data(2:end);
%% construct gradient: gr = A*GR

A = zeros(length(frequencies),sum(frequencies));
for j = 2:length(frequencies)
    A(j,sumfreq(j-1)+1:sumfreq(j-1)+frequencies(j)) = 1;
end
A(1,1:frequencies(1)) = 1;
A = sparse(A);
A = kron(eye(2*Nch),A);
gr = A*GR;
