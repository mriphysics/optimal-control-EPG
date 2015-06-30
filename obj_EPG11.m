%% EPG forward for FSE
%% objective function. maximizes the signal 
function [obj, grad,FF] = obj_EPG11(params,ESP,T1,T2,c,B1,target)
% efficient implementation of blockdiagonal matrix vec multiplication
% w.r.t. obj_EPG8
%
% c is vector of {0,1} samplings, for instance: c(t) = 1 counts, c(t) = 0
% does not count.
% includes optimization w.r.t. phi (phase)
% includes 2D spatially resolved , Nch channels
% B1 :  Nvoxels x Nch array
% params: [real(theta1) real(theta2) ... imag(theta1) imag(theta2) ...];
% where index stands for channel
% target: Nt x Ns target states
%
% NOTE: objective w.r.t. real and imaginary parts of theta!
%
% Alessandro Sbrizzi March 2015 (based on Shaihan Malik EPG implementation)
[Ns, Nch] = size(B1);
F0 = 0;

np = length(params)/(2*Nch);

x = reshape(params(1:np*Nch),np,Nch).'; %real parts
y = reshape(params(np*Nch+1:end),np,Nch).'; % imaginary parts
th = zeros(Ns,np,Nch);
for ch = 1:Nch
    th(:,:,ch) = B1(:,ch)*(x(ch,:)+1i*y(ch,:)); % effective theta
end
th = sum(th,3);

%th(:,2:end) = th(:,2:end)*exp(1i*pi/2); % not working!!
alph = abs(th);% effective alpha
ph = angle(th); % effective phi

% add CPMG phase
ph(:,2:end) = ph(:,2:end) + pi/2; % working!
xeff = real(th); % effective real
yeff = imag(th); % effective imag


kmax = 2*np - 1; % up to just before the next RF pulse
N = 3*(kmax-1)/2; % number of states in total
% split magnitude and phase




Npathway = inf;
klimit=false;
flipback=false;

% enforce pathway limit
if (N>Npathway)&&~klimit
    N=Npathway;
end

if klimit
    nr = size(th,2)-1; % number of refocus pulses
    kmax = nr - 1 + mod(nr,2);
    N = 1.5*(nr+mod(nr,2)); % number of states in total
    if mod(nr,2) 
        % odd
        KMAX = [1:2:kmax (kmax-2):-2:1];
    else
        %even
        KMAX = [1:2:kmax kmax:-2:1];
    end
    NMAX = 1.5*(KMAX+1);
else
    % previous implementation
    NMAX = 3:3:3*(np-1);
    NMAX(NMAX>N)=(N-mod(N,3));
end
    
    
%% ==== build Shift matrix, S with indices 
S = zeros([N N]);
%%% F(k>1) look @ states just BEFORE np+1 pulse
kidx = 4:3:N; % miss out F1+
sidx = kidx-3;
%idx = sub2ind([N N],kidx,sidx);
idx = kidx + N*(sidx-1);
S(idx)=1;

%%% F(k<1) 
kidx = 2:3:N;
kidx(end)=[];% most negative state relates to nothing; related row is empty
sidx = kidx+3;
% ix = sub2ind([N N],kidx,sidx);
ix = kidx + N*(sidx-1);
S(ix)=1;

%%% Z states
kidx = 3:3:N;
%ix = sub2ind([N N],kidx,kidx);
ix = kidx + N*(kidx-1);
S(ix)=1;

%%% finally F1+ - relates to F-1-
S(1,2)=1;
S(end-1,end-2) = 1;
S = sparse(S);
%% Relaxation =====: note, Z0 regrowth not handled (No Z0 state here)
E1=exp(-ESP/T1);
E2=exp(-ESP/T2);
%R = diag(repmat([E2 E2 E1],[1 kmax+1]));
R=eye(N);
ii = 1:(N/3);
R(3*N*(ii-1)+3*(ii-1)+1)=E2;
R(3*N*ii-2*N+3*(ii-1)+2)=E2;
R(3*N*ii-N+3*(ii-1)+3)=E1;

%%%% composites
RS=sparse(R*S);

%% F matrix (many elements zero, not efficient)
FF = zeros(N*2,np+1,Ns);

% Now the other states
for ns = 1:Ns % loop over space
    
    F = zeros([N*2 np+1]); %% records state 
    % initial state:
    F(3,1) = 1; 
    for jj=2:np+1 %loop over time

        if jj == 2 % excitation
            % Excitation pulse
            A = Trot_fun(alph(ns,jj-1),ph(ns,jj-1));
            F([1:3 N+1:N+3],jj) = [real(A); imag(A)]*F(1:3,jj-1); %<---- state straight after excitation [F0 F0* Z0]
        elseif jj == 3 % first refocusing
            A = Trot_fun(alph(ns,jj-1),ph(ns,jj-1));
            T = build_T_matrix_sub(A,3);
            F(1,jj) = exp(-0.5*ESP/T2)*F(1,jj-1);
            F(N+1,jj) = exp(-0.5*ESP/T2)*F(N+1,jj-1);
            F([1:3 N+1:N+3],jj) = [real(T), -imag(T);imag(T),real(T)]*[F(1:3,jj);F(N+1:N+3,jj)];    
        else
            A = Trot_fun(alph(ns,jj-1),ph(ns,jj-1));
            temp1 = RS*F(1:N,jj-1);
            temp2 = RS*F(N+1:2*N,jj-1);
            temp2(1) = -temp2(1);% D*temp;
            F(:,jj) = [blockdiag_mult(real(A),temp1)-blockdiag_mult(imag(A),temp2); ...
                       blockdiag_mult(imag(A),temp1)+blockdiag_mult(real(A),temp2)];
        end
    end
    FF(:,:,ns) = F; 
end
state = FF(2,:,:)+1i*FF(N+2,:,:);
state = diag(c)*(squeeze(state)-target);
obj = 0.5*norm(state(:))^2;

%% now adjoint states
if nargout > 1 % need gradient 

%C = zeros(2*N,2*N);
%C(2,2) = 1; % sampling matrix
%C(N+2,N+2) = 1; % sampling matrix

LL = zeros(np+1,N*2,Ns);
RRSS = [RS, RS*0;RS*0, RS];
RRSS(N+1,:) = -RRSS(N+1,:);
% Now the other states
for ns = 1:Ns % loop over space
    L = zeros(np+1,N*2);
    for jj=np:-1:1

        if jj == 1 % 
            A = Trot_fun(alph(ns,jj+1),ph(ns,jj+1));
            T = build_T_matrix_sub(A,3);
            d = zeros(6,6);
            d(1) = exp(-0.5*ESP/T2);d(3+1,3+1) = exp(-0.5*ESP/T2);
            L(jj,[1:3 N+1:N+3]) = L(jj+1,[1:3 N+1:N+3])*[real(T) -imag(T);imag(T) real(T)]*d+[0,c(jj+1)*(FF(2,jj+1,ns)-real(target(jj+1,ns))),0,0,c(jj+1)*(FF(N+2,jj+1,ns)-imag(target(jj+1,ns))),0];
        
        elseif jj == np
        
            temp = zeros(1,2*N);
            temp(2) = c(jj+1)*(FF(2,jj+1,ns)-real(target(jj+1,ns)));
            temp(N+2) = c(jj+1)*(FF(N+2,jj+1,ns)-imag(target(jj+1,ns)));
            L(jj,:) = temp; 
        else

            A = Trot_fun(alph(ns,jj+1),ph(ns,jj+1));
            % NMAX is maximum state index
            kidx = 1:NMAX(end);
            temp = [blockdiag_mult(real(A).',L(jj+1,1:N).')+blockdiag_mult(imag(A).',L(jj+1,N+1:end).');...
                   -blockdiag_mult(imag(A).',L(jj+1,1:N).')+blockdiag_mult(real(A).',L(jj+1,N+1:end).')].';
            temp = temp*RRSS;
            temp1 = zeros(1,2*N);
            temp1(2) = c(jj+1)*(FF(2,jj+1,ns)-real(target(jj+1,ns)));
            temp1(N+2) = c(jj+1)*(FF(N+2,jj+1,ns)-imag(target(jj+1,ns)));
            L(jj,:) = temp+temp1;%-(F(:,jj+1))'*C;

        end
    end
    LL(:,:,ns) = L; 
end

%% now gradient

GRAD = zeros(2*np*Nch,Ns);
for ns = 1:Ns % loop over space
    for jj = 2:np+1
        dT_daeff = dTrot_fun_da(alph(ns,jj-1),ph(ns,jj-1));
        dT_dpeff = dTrot_fun_dp(alph(ns,jj-1),ph(ns,jj-1));
        d = zeros(2*N,2*N);
        d(1,1) = exp(-0.5*ESP/T2);d(N+1,N+1) = exp(-0.5*ESP/T2);
        temp1 = d*FF(:,jj-1,ns);
        temp2 = RRSS*FF(:,jj-1,ns);
        TEMP1a = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_daeff),temp1(1:N))-blockdiag_mult(imag(dT_daeff),temp1(N+1:end));...
                                blockdiag_mult(imag(dT_daeff),temp1(1:N))+blockdiag_mult(real(dT_daeff),temp1(N+1:end))];
        TEMP1p = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_dpeff),temp1(1:N))-blockdiag_mult(imag(dT_dpeff),temp1(N+1:end));...
                                blockdiag_mult(imag(dT_dpeff),temp1(1:N))+blockdiag_mult(real(dT_dpeff),temp1(N+1:end))];
        TEMP2a = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_daeff),temp2(1:N))-blockdiag_mult(imag(dT_daeff),temp2(N+1:end));...
                                blockdiag_mult(imag(dT_daeff),temp2(1:N))+blockdiag_mult(real(dT_daeff),temp2(N+1:end))];
        TEMP2p = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_dpeff),temp2(1:N))-blockdiag_mult(imag(dT_dpeff),temp2(N+1:end));...
                                blockdiag_mult(imag(dT_dpeff),temp2(1:N))+blockdiag_mult(real(dT_dpeff),temp2(N+1:end))];
        
        for ch = 1:Nch
            dadxj = xeff(ns,jj-1)/alph(ns,jj-1)*real(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)*imag(B1(ns,ch));
            dpdxj = xeff(ns,jj-1)/alph(ns,jj-1)^2*imag(B1(ns,ch))-yeff(ns,jj-1)/alph(ns,jj-1)^2*real(B1(ns,ch));
            dadyj = -xeff(ns,jj-1)/alph(ns,jj-1)*imag(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)*real(B1(ns,ch));
            dpdyj = xeff(ns,jj-1)/alph(ns,jj-1)^2*real(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)^2*imag(B1(ns,ch));
            if jj == 3
                GRAD(jj-1+(ch-1)*np,ns) = dadxj*TEMP1a+dpdxj*TEMP1p;
                GRAD(np*Nch+jj-1+(ch-1)*np,ns) = dadyj*TEMP1a+dpdyj*TEMP1p;
            else
                GRAD(jj-1+(ch-1)*np,ns) = dadxj*TEMP2a+dpdxj*TEMP2p; 
                GRAD(np*Nch+jj-1+(ch-1)*np,ns) = dadyj*TEMP2a+dpdyj*TEMP2p; 
            end
        end
    end
end

grad = real(sum(GRAD,2));
grad = grad(:);
end
%%
%sig = exp(-0.5*ESP/T2)*conj(F(2,:))*exp(1i*pi/2); % conjugate because it's F-1*
 
%% Flipback/DRIVE
if flipback
    f0 = conj(F(2,end))*exp(-0.5*ESP/T2);% 8-2-12: Add conj here
        
    % rotation matrix for flip back pulse
    Af = Trot_fun(abs(alpha_fb),angle(alpha_fb));

    % state directly after flip back
    F0 = Af*[f0;conj(f0);0];% IGNORE Z0 here
end

    % reduced version
    function T =  build_T_matrix_sub(AA,nn)
         T=zeros([nn nn]);
%          for ii = 1:(nn/3);
%             T(3*nn*(ii-1)+3*(ii-1)+1)=AA(1);
%             T(3*nn*(ii-1)+3*(ii-1)+2)=AA(2);
%             T(3*nn*(ii-1)+3*(ii-1)+3)=AA(3);
%             T(3*nn*ii-2*nn+3*(ii-1)+1)=AA(4);
%             T(3*nn*ii-2*nn+3*(ii-1)+2)=AA(5);
%             T(3*nn*ii-2*nn+3*(ii-1)+3)=AA(6);
%             T(3*nn*ii-nn+3*(ii-1)+1)=AA(7);
%             T(3*nn*ii-nn+3*(ii-1)+2)=AA(8);
%             T(3*nn*ii-nn+3*(ii-1)+3)=AA(9);
%          end
            ind = 1:(nn/3);
            T(3*nn*(ind-1)+3*(ind-1)+1)=AA(1);
            T(3*nn*(ind-1)+3*(ind-1)+2)=AA(2);
            T(3*nn*(ind-1)+3*(ind-1)+3)=AA(3);
            T(3*nn*ind-2*nn+3*(ind-1)+1)=AA(4);
            T(3*nn*ind-2*nn+3*(ind-1)+2)=AA(5);
            T(3*nn*ind-2*nn+3*(ind-1)+3)=AA(6);
            T(3*nn*ind-nn+3*(ind-1)+1)=AA(7);
            T(3*nn*ind-nn+3*(ind-1)+2)=AA(8);
            T(3*nn*ind-nn+3*(ind-1)+3)=AA(9);

        T = sparse(T);
%        T = sparse(kron(eye(nn/3),AA));% slow
    end

    function w = blockdiag_mult(A,v)
        %given A square matrix, it returns w =  kron(A,eye(N))*v efficiently
        nA = size(A,1); %square
        V = reshape(v,nA,length(v)/nA);
        W = A*V;
        w = W(:);
    end


    % Rotation matrix direct definition
    function T = Trot_fun(a,p)
        
        T = zeros([3 3]);
        T(1) = cos(a/2).^2;
        T(2) = exp(-2*1i*p)*(sin(a/2)).^2;
        T(3) = -0.5*1i*exp(-1i*p)*sin(a);
        T(4) = conj(T(2));
        T(5) = T(1);
        T(6) = 0.5*1i*exp(1i*p)*sin(a);
        T(7) = -1i*exp(1i*p)*sin(a);
        T(8) = 1i*exp(-1i*p)*sin(a);
        T(9) = cos(a);
    end
    
    function T = dTrot_fun_da(a,p)
        % derivative Rotation matrix w.r.t. a
        T = zeros([3 3]);
        T(1) = -sin(a)/2;%cos(a/2).^2;
        T(2) = 0.5*exp(-2*1i*p)*sin(a);%exp(-2*1i*p)*(sin(a/2)).^2;
        T(3) = -0.5*1i*exp(-1i*p)*cos(a);%-0.5*1i*exp(-1i*p)*sin(a);
        T(4) = conj(T(2));%conj(T(2));
        T(5) = T(1);%T(1);
        T(6) = 0.5*1i*exp(1i*p)*cos(a);%0.5*1i*exp(1i*p)*sin(a);
        T(7) = -1i*exp(1i*p)*cos(a);%-1i*exp(1i*p)*sin(a);
        T(8) = 1i*exp(-1i*p)*cos(a);%1i*exp(-1i*p)*sin(a);
        T(9) = -sin(a);%cos(a);
    end
    function T = dTrot_fun_dp(a,p)
        % derivative Rotation matrix w.r.t. p
        T = zeros([3 3]);
        T(1) = 0;%cos(a/2).^2;
        T(2) = -2*1i*exp(-2*1i*p)*(sin(a/2)).^2;
        T(3) = -0.5*exp(-1i*p)*sin(a);
        T(4) = conj(T(2));
        T(5) = T(1);%T(1);
        T(6) = -0.5*exp(1i*p)*sin(a);
        T(7) = exp(1i*p)*sin(a);
        T(8) = exp(-1i*p)*sin(a);
        T(9) = 0;%cos(a);
    end

end

