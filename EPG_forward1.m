%% EPG forward for FSE

function [F, grad] = EPG_forward1(theta,varargin)

% F0 is the FID/Z0 created by the flipback pulse. Initialize in case not set
F0 = 0;

np = length(theta);
kmax = 2*np - 1; % up to just before the next RF pulse
N = 3*(kmax-1)/2; % number of states in total



% split magnitude and phase
alpha = abs(theta);phi=angle(theta);

% add CPMG phase
%alpha(2:end) = alpha(2:end)*exp(1i*pi/2);
phi(2:end) = phi(2:end) + pi/2;

%% get variables
T1 = inf;
T2 = inf;
ESP=10;
Npathway = inf;
klimit=false;
flipback=false;

for ii=1:length(varargin)
    
    if strcmp(varargin{ii},'T1')
        T1 = varargin{ii+1};
    end
    if strcmp(varargin{ii},'T2')
        T2 = varargin{ii+1};
    end
    if strcmp(varargin{ii},'target')
        target = varargin{ii+1};
    end
        if strcmp(varargin{ii},'Fend')
        Fend = varargin{ii+1};
    end
    if strcmp(varargin{ii},'ESP')||strcmp(varargin{ii},'TE')
        ESP = varargin{ii+1};
    end
    % # of coherence pathways to consider (default is inf)
    if strcmp(varargin{ii},'Npathway')||strcmp(varargin{ii},'npath')
        Npathway = varargin{ii+1};
    end
    % more drastic version of above (see lab book 9-7-11)
    if strcmp(varargin{ii},'klimit')
        klimit=true;
    end
    % Allow user to define starting Equilibrium magnetization
    if strcmp(varargin{ii},'E')||strcmp(varargin{ii},'M0')
        E = varargin{ii+1};
    end
    % Simulate Flip-back pulse at the time of the last echo
    if strcmpi(varargin{ii},'flipback')||strcmpi(varargin{ii},'drive')
        flipback = true;
        alpha_fb = varargin{ii+1}; % COMPLEX
    end
end

% enforce pathway limit
if (N>Npathway)&&~klimit
    N=Npathway;
end

if klimit
    nr = length(theta)-1; % number of refocus pulses
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
RS=R*S;

%% F matrix (many elements zero, not efficient)
F = sparse(zeros([N np+1])); %% records state 
% initial state:
F(3,1) = 1; 
% Now the other states
for jj=2:np+1 
    if jj == 2 % excitation
        % Excitation pulse
        A = Trot_fun(alpha(jj-1),phi(jj-1));
        F(1:3,jj) = A*F(1:3,jj-1); %<---- state straight after excitation [F0 F0* Z0]
    elseif jj == 3 % first refocusing
        A = Trot_fun(alpha(jj-1),phi(jj-1));
        T = build_T_matrix_sub(A,3);
        F(1,jj) = exp(-0.5*ESP/T2)*F(1,jj-1);
        F(1:3,jj) = T*F(1:3,jj);        
    else
        A = Trot_fun(alpha(jj-1),phi(jj-1));
        % NMAX is maximum state index
        kidx = 1:NMAX(end);
        T = build_T_matrix_sub(A,NMAX(end));
        %FF(1)=conj(FF(1)); %<---- complex conjugate links the +/- states
        temp = RS*F(kidx,jj-1);
        temp(1)=conj(temp(1));
        F(kidx,jj) = T*temp;
    end
end
%% now adjoint states
if nargout > 1 % need gradient
C = zeros(N,N);
C(2,2) = 1; % sampling matrix
C = sparse(C);
L = sparse(zeros(np+1,N));
% Now the other states
for jj=np:-1:1
    if jj == 1 % 
        A = Trot_fun(alpha(jj+1),phi(jj+1));
        T = build_T_matrix_sub(A,3);
        d = zeros(3,3);
        d(1) = exp(-0.5*ESP/T2);d = sparse(d);
        L(jj,1:3) = L(jj+1,1:3)*T*d+[0,conj(F(2,jj+1)-target(2,jj+1)),0];
    elseif jj == np
        L(jj,:) = (F(:,jj+1)-target(:,jj+1))'*C;
    else

        A = Trot_fun(alpha(jj+1),phi(jj+1));
        % NMAX is maximum state index
        kidx = 1:NMAX(end);
        T = build_T_matrix_sub(A,NMAX(end));
        temp = L(jj+1,:)*T*RS;
        %temp(1) = conj(temp(1));
        % First evolve, then conj, then flip
        L(jj,:) = temp+(F(:,jj+1)-target(:,jj+1))'*C;

    end
end

%% now gradient
grad = zeros(np,1);
for jj = 2:np+1
    dT = build_T_matrix_sub(dTrot_fun_da(alpha(jj-1),phi(jj-1)),NMAX(end));
    if jj == 3
        d = zeros(N,N);
        d(1) = exp(-0.5*ESP/T2);d = sparse(d);
        grad(jj-1) = real(L(jj-1,:)*dT*d*F(:,jj-1));
    else
        grad(jj-1) = real(L(jj-1,:)*dT*RS*F(:,jj-1));
    end
end
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
        ii = 1:(nn/3);
        T(3*nn*(ii-1)+3*(ii-1)+1)=AA(1);
        T(3*nn*(ii-1)+3*(ii-1)+2)=AA(2);
        T(3*nn*(ii-1)+3*(ii-1)+3)=AA(3);
        T(3*nn*ii-2*nn+3*(ii-1)+1)=AA(4);
        T(3*nn*ii-2*nn+3*(ii-1)+2)=AA(5);
        T(3*nn*ii-2*nn+3*(ii-1)+3)=AA(6);
        T(3*nn*ii-nn+3*(ii-1)+1)=AA(7);
        T(3*nn*ii-nn+3*(ii-1)+2)=AA(8);
        T(3*nn*ii-nn+3*(ii-1)+3)=AA(9);
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
    % derivative Rotation matrix w.r.t. alpha
    function T = dTrot_fun_da(a,p)
        
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

end

