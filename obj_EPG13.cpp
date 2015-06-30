// Example EPG simulator

#include <iostream>
//#include <math.h>
#include "armadillo"
#include <algorithm>
#include <ctime>
using namespace arma;
using namespace std;

arma::uvec index_vec(int start, int interval, int finish) // makes vector of INTEGER indeces as in matlab: vector = [start:interval:finish]
	{
	uvec vector = zeros<uvec>((finish-start)/interval+1);
	
	for (u32 ind=0;ind<vector.n_elem;++ind){
		vector(ind) = start+ind*interval;
	}
	return vector;
}

arma::sp_mat makeS(int N) // return shift matrix
	{
 	using namespace arma;
	//mat S = zeros<mat>(N, N);
 	sp_mat S(N, N);
 	for(u32 ind=0; ind<N/3.0-1.0; ++ind){
		S(ind*3+2.0,ind*3+2.0)=1.0; // diag elems
		S(ind*3+1.0,ind*3+4.0)=1.0;
		S(ind*3+3.0,ind*3+0.0)=1.0;
	}
	//S(0,0) = 1.0;
	S(0,1) = 1.0;
	S(N-1.0,N-1.0) = 1.0;
	//S(N-2.0,N-3.0) = 1.0;
 	//uvec kidx = ones<uvec>(N/3.0);
 	
 	//S.cols(kidx) += 1.0;
 	return S;
}
arma::sp_mat makeR(int N, double E1, double E2) // return decay matrix
	{
 	using namespace arma;
 	sp_mat R(N, N);
 	for(u32 ind=0; ind<N/3.0-1.0; ++ind){
		R(ind*3,ind*3)=E2;
		R(ind*3+1.0,ind*3+1.0)=E2;
		R(ind*3+2.0,ind*3+2.0)=E1;
	}
	R(N-1.0,N-1.0) = E1;
	R(N-2.0,N-2.0) = E2;
	R(N-3.0,N-3.0) = E2;
 	return R;
}

arma::cx_mat Trot_fun(double a, double p) // % Rotation matrix direct definition
	{
 	using namespace arma;
 	cx_mat T = zeros<cx_mat>(3, 3);
 	std::complex<double> ii (0.0,1.0);
	T(0,0) = pow(cos(a/2.0),2);
	T(1,0) = exp(-2.0*ii*p)*pow(sin(a/2.0),2);
	T(2,0) = -0.5*ii*exp(-ii*p)*sin(a);
	T(0,1) = conj(T(1,0));
	T(1,1) = T(0,0);
	T(2,1) = 0.5*ii*exp(ii*p)*sin(a);
	T(0,2) = -ii*exp(ii*p)*sin(a);
	T(1,2) = ii*exp(-ii*p)*sin(a);
	T(2,2) = cos(a);
 	return T;
}

arma::cx_mat dTrot_fun_da(double a, double p) // % derivative Rotation matrix w.r.t. a
	{
 	using namespace arma;
    // derivative Rotation matrix w.r.t. a
    cx_mat T = zeros<cx_mat>(3, 3);
    std::complex<double> ii (0.0,1.0);
    T(0,0) = -sin(a)/2.0;
    T(1,0) = 0.5*exp(-2.0*ii*p)*sin(a);
    T(2,0) = -0.5*ii*exp(-ii*p)*cos(a);
    T(0,1) = conj(T(1,0));
    T(1,1) = T(0,0);
    T(2,1) = 0.5*ii*exp(ii*p)*cos(a);
    T(0,2) = -ii*exp(ii*p)*cos(a);
    T(1,2) = ii*exp(-ii*p)*cos(a);
    T(2,2) = -sin(a);
    
 	return T;
}

arma::cx_mat dTrot_fun_dp(double a, double p) // % derivative Rotation matrix w.r.t. p
	{
 	using namespace arma;
    // derivative Rotation matrix w.r.t. p
    cx_mat T = zeros<cx_mat>(3, 3);
    std::complex<double> ii (0.0,1.0);
    T(0,0) = 0.0;
    T(1,0) = -2.0*ii*exp(-2.0*ii*p)*pow(sin(a/2.0),2);
    T(2,0) = -0.5*exp(-ii*p)*sin(a);
    T(0,1) = conj(T(1,0));
    T(1,1) = T(0,0);
    T(2,1) = -0.5*exp(ii*p)*sin(a);
    T(0,2) = exp(ii*p)*sin(a);
    T(1,2) = exp(-ii*p)*sin(a);
    T(2,2) = 0.0;
    
 	return T;
}


arma::vec blockmat_vec(arma::mat A, arma::vec v) 
	// % compute kron(A,eye(n))*v to replace large block diagonal sparse operator
	// A supposed to be squared m x m
	{	
 	using namespace arma;
 	int n = v.n_elem;
 	int m = A.n_rows;
 	colvec w;
 	mat W;
 	mat V = reshape(v,m,n/m);
 	W = A*V;
 	w = vectorise(W);
 	return w;
 }

int main(int argc, char** argv)
  	{
  	// read data ========================================================================================================
  	colvec data;
  	clock_t begin = clock();
	data.load("input_data.out", raw_ascii);  // default save format is arma_binary
	clock_t end = clock();
	
	int nargout = data(0);
	int Nt = data(1);
	int Ns = data(2);
	int Nch = data(3);
    int klim = data(4); // maximum order of k
	int start_params = 5;
	int end_params   = start_params-1+Nt*Nch*2.0;
	double ESP = data(end_params+1);
	double T1 = data(end_params+2);
	double T2 = data(end_params+3);
	int start_c = end_params+4;
	int end_c   = start_c+Nt;
	int start_B1r = end_c+1;
	int end_B1r   = start_B1r-1+Ns*Nch;
	int start_B1i = end_B1r+1;
	int end_B1i   = start_B1i-1+Ns*Nch;
	int start_Tr  = end_B1i+1;
	int end_Tr    = start_Tr-1+(Nt+1)*Ns;
	int start_Ti  = end_Tr+1;
	int end_Ti    = start_Ti-1+(Nt+1)*Ns;
	
	rowvec params = trans(data(span(start_params,end_params))); // input tip angles (real and imag parts)
	// finished loading ======================================================================================================
	vec c = data(span(start_c,end_c));
	// make B1+
	mat B1r = reshape(data(span(start_B1r,end_B1r)),Ns,Nch);
	mat B1i = reshape(data(span(start_B1i,end_B1i)),Ns,Nch);
	cx_mat B1 = cx_mat(B1r,B1i);
	// make TARGET
	mat Tr = reshape(data(span(start_Tr,end_Tr)),Nt+1,Ns);
	mat Ti = reshape(data(span(start_Ti,end_Ti)),Nt+1,Ns);
	cx_mat target = cx_mat(Tr,Ti);
	
  	int np;
  	const double pi = 3.141592653589793;
  	std::complex<double> ii (0.0,1.0);
  	//ii = (0.0,1.0); // imaginary unit
	
	np = params.n_cols/(2*Nch);
	mat x = trans(reshape(params.cols(0,np*Nch-1.0),np,Nch));
	mat y = trans(reshape(params.cols(np*Nch,2.0*np*Nch-1.0),np,Nch));
	cx_mat xy = cx_mat(x,y); // assemble complex
	cx_mat th = B1*xy;
	
	mat alph = abs(th);
	mat ph = imag(log(th)); // angle function (!)
	
	// add CPMG phase
	ph.cols(1,ph.n_cols-1) = ph.cols(1,ph.n_cols-1) + pi/2;
	
	mat xeff = real(th); // effective real
    mat yeff = imag(th); // effective imag
	int kmax = min(2*np - 1,2*klim - 1); // up to just before the next RF pulse
	int N = 3*(kmax-1)/2; 
	// ==== build Shift matrix, S with indices 
	sp_mat S = makeS(N);
	//Relaxation =====: note, Z0 regrowth not handled (No Z0 state here)
	double E1 = exp(-ESP/T1);
	double E2 = exp(-ESP/T2);
	sp_mat R = makeR(N,E1,E2);
	// composites
	sp_mat RS = R*S;
	vec temp;
	vec temp1;
	// F matrix (many elements zero, not efficient)
	cube FF = zeros<cube>(N*2,np+1,Ns);
	// Now the other states
	uvec colind;
    uvec rowind;
    colvec TEMP1;
    colvec TEMP2;
    rowind << 0<<1<<2<<N<<N+1<<N+2;
    // make locations for sparse matrix TT
    begin = clock();
    cx_mat A(3,3);
	for (u32 ns = 0; ns<Ns; ++ns){// loop over space
		mat F = zeros<mat>(N*2,np+1); // records state 
		// initial state:
    	F(2,0) = 1.0; 
    	for (u32 jj=1; jj<np+1; ++jj){ //loop over time
    		colind <<jj;
    		if (jj == 1){ // excitation
    			A = Trot_fun(alph(ns,jj-1),ph(ns,jj-1));
    			F(rowind,colind) = join_cols(real(A),imag(A))*F(span(0,2),span(jj-1,jj-1)); // <---- state straight after excitation [F0 F0* Z0]
    			
    		}
    		
    		else if (jj == 2){ // first refocusing
    			A = Trot_fun(alph(ns,jj-1),ph(ns,jj-1));
    			F(0,jj) = exp(-0.5*ESP/T2)*F(0,jj-1);
    			F(N,jj) = exp(-0.5*ESP/T2)*F(N,jj-1);
    			TEMP1 = blockmat_vec(real(A),F(span(0,2),span(jj,jj)))-blockmat_vec(imag(A),F(span(N,N+2),span(jj,jj)));
    			TEMP2 = blockmat_vec(imag(A),F(span(0,2),span(jj,jj)))+blockmat_vec(real(A),F(span(N,N+2),span(jj,jj)));
    			F(rowind,colind) = join_cols(TEMP1,TEMP2);
    			
    		}
    		else{
    			A = Trot_fun(alph(ns,jj-1),ph(ns,jj-1));
    			temp = join_cols(RS*F(span(0,N-1),span(jj-1,jj-1)),RS*F(span(N,2*N-1),span(jj-1,jj-1)));
    			temp(N) = -temp(N);// D*temp;
    			TEMP1 = blockmat_vec(real(A),temp(span(0,N-1)))-blockmat_vec(imag(A),temp(span(N,2*N-1)));
    			TEMP2 = blockmat_vec(imag(A),temp(span(0,N-1)))+blockmat_vec(real(A),temp(span(N,2*N-1)));
    			F.col(jj) = join_cols(TEMP1,TEMP2);    			
    		}    		
    	}
    	FF.slice(ns) = F;
	}
	mat stater = FF(span(1),span(),span());
	mat statei = FF(span(N+1),span(),span());
	cx_mat state(stater,statei);
	cx_mat residual = diagmat(c)*(state-target);
	mat obj = 0.5*abs((trans(vectorise(residual))*vectorise(residual)));
	
	colvec grad;
	
	if (nargout == 2){// nargout = 2 means that gradient is needed
	//===============================================================================================================
	// Adjoint states calculation
	cube LL = zeros<cube>(np+1,N*2,Ns);
	sp_mat RRSS = join_cols(join_rows(RS, RS*0.0),join_rows(RS*0.0, RS));
	RRSS(span(N),span()) = -RRSS(span(N),span());
	
	// Now the other states
	colind << 0<<1<<2<<N<<N+1<<N+2;
	rowvec temp_vec;
	rowvec tempv = zeros<rowvec>(2*N);
    mat tempm;
    mat d = zeros<mat>(6,6);
    rowvec Temp1 = zeros<rowvec>(2*N);
	for (u32 ns = 0; ns<Ns; ++ns){// loop over space
		mat L = zeros<mat>(np+1,N*2);
		for (int jj=np-1; jj>-1; --jj){ //loop over time
		
			if (jj == 0){ 
				A = Trot_fun(alph(ns,jj+1),ph(ns,jj+1));
				d.zeros();
				d(0,0) = exp(-0.5*ESP/T2);d(3,3) = exp(-0.5*ESP/T2); 
				rowind<<jj;
				temp_vec<<0.0<<c(jj+1)*(FF(2,jj+1,ns)-real(target(jj+1,ns)))<<0.0<<0.0<<c(jj+1)*(FF(N+2,jj+1,ns)-imag(target(jj+1,ns)))<<0.0;
				L(rowind,colind) = L(rowind+1,colind)*join_cols(join_rows(real(A),-imag(A)),join_rows(imag(A),real(A)))*d+temp_vec;
			}
			else if (jj == np-1){
				tempv.zeros();
				tempv(1) = c(jj+1)*(FF(1,jj+1,ns)-real(target(jj+1,ns)));
				tempv(N+1) = c(jj+1)*(FF(N+1,jj+1,ns)-imag(target(jj+1,ns)));
				L.row(jj) = tempv; 
			}
    		else{
                A = Trot_fun(alph(ns,jj+1),ph(ns,jj+1));
                tempm = trans(join_cols(
                        blockmat_vec(trans(real(A)),trans(L(jj+1,span(0,N-1))))+blockmat_vec(trans(imag(A)),trans(L(jj+1,span(N,2*N-1)))),
                       -blockmat_vec(trans(imag(A)),trans(L(jj+1,span(0,N-1))))+blockmat_vec(trans(real(A)),trans(L(jj+1,span(N,2*N-1))))
                       ));	
                tempm = tempm*RRSS;
                Temp1.zeros();
                Temp1(1) = c(jj+1)*(FF(1,jj+1,ns)-real(target(jj+1,ns)));
                Temp1(N+1) = c(jj+1)*(FF(N+1,jj+1,ns)-imag(target(jj+1,ns)));
                L.row(jj) = tempm+Temp1;
			}
		
		}
        LL.slice(ns) = L; 
	
	}
    // now gradient ========================================================================
    cx_mat dT_daeff;
    cx_mat dT_dpeff;
    mat dd = zeros<mat>(2*N,2*N);
    mat GRAD = zeros<mat>(2*np*Nch,Ns);
    vec vec1;
    vec temp2;
    rowvec LL1;
    mat TEMP1a,TEMP1p,TEMP2a,TEMP2p;
    double dadxj,dpdxj,dadyj,dpdyj;
    for (u32 ns = 0; ns<Ns; ++ns){// loop over space
        for (u32 jj=1; jj<np+1; ++jj){ //loop over time
            dT_daeff = dTrot_fun_da(alph(ns,jj-1),ph(ns,jj-1));
            dT_dpeff = dTrot_fun_dp(alph(ns,jj-1),ph(ns,jj-1));
            dd.zeros();
            dd(0,0) = exp(-0.5*ESP/T2);dd(N,N) = exp(-0.5*ESP/T2);
            vec1 = FF(span(),span(jj-1),span(ns));
            temp1 = dd*vec1;
            temp2 = RRSS*vec1;
            LL1 = LL(span(jj-1),span(),span(ns));
            TEMP1a = LL1*join_cols(blockmat_vec(real(dT_daeff),temp1(span(0,N-1)))-blockmat_vec(imag(dT_daeff),temp1(span(N,2*N-1))),
                                   blockmat_vec(imag(dT_daeff),temp1(span(0,N-1)))+blockmat_vec(real(dT_daeff),temp1(span(N,2*N-1))));
            TEMP1p = LL1*join_cols(blockmat_vec(real(dT_dpeff),temp1(span(0,N-1)))-blockmat_vec(imag(dT_dpeff),temp1(span(N,2*N-1))),
                                   blockmat_vec(imag(dT_dpeff),temp1(span(0,N-1)))+blockmat_vec(real(dT_dpeff),temp1(span(N,2*N-1))));
            TEMP2a = LL1*join_cols(blockmat_vec(real(dT_daeff),temp2(span(0,N-1)))-blockmat_vec(imag(dT_daeff),temp2(span(N,2*N-1))),
                                   blockmat_vec(imag(dT_daeff),temp2(span(0,N-1)))+blockmat_vec(real(dT_daeff),temp2(span(N,2*N-1))));
            TEMP2p = LL1*join_cols(blockmat_vec(real(dT_dpeff),temp2(span(0,N-1)))-blockmat_vec(imag(dT_dpeff),temp2(span(N,2*N-1))),
                                   blockmat_vec(imag(dT_dpeff),temp2(span(0,N-1)))+blockmat_vec(real(dT_dpeff),temp2(span(N,2*N-1))));
            for (u32 ch = 0; ch<Nch; ++ch){
                dadxj = xeff(ns,jj-1)/alph(ns,jj-1)*real(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)*imag(B1(ns,ch));
                dpdxj = xeff(ns,jj-1)/pow(alph(ns,jj-1),2)*imag(B1(ns,ch))-yeff(ns,jj-1)/pow(alph(ns,jj-1),2)*real(B1(ns,ch));
                dadyj = -xeff(ns,jj-1)/alph(ns,jj-1)*imag(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)*real(B1(ns,ch));
                dpdyj = xeff(ns,jj-1)/pow(alph(ns,jj-1),2)*real(B1(ns,ch))+yeff(ns,jj-1)/pow(alph(ns,jj-1),2)*imag(B1(ns,ch));
                if (jj == 2){ 
                    GRAD(jj-1+(ch)*np,ns)= as_scalar(dadxj*TEMP1a+dpdxj*TEMP1p);
                    GRAD(np*Nch+jj-1+(ch)*np,ns) = as_scalar(dadyj*TEMP1a+dpdyj*TEMP1p);
                }
                else{
                    GRAD(jj-1+(ch)*np,ns) = as_scalar(dadxj*TEMP2a+dpdxj*TEMP2p); 
                    GRAD(np*Nch+jj-1+(ch)*np,ns) = as_scalar(dadyj*TEMP2a+dpdyj*TEMP2p); 
                }
            }
        }
    }
    grad = real(sum(GRAD,1));
    } // end gradient calculation if nargout = 2;
    if (nargout == 2){
    	colvec output = join_cols(obj,grad);
    	output.save("output_data.out", raw_ascii);
    }
    else{
    	obj.save("output_data.out", raw_ascii);
    }
    	
	end = clock();
	
	
  	return 0;
  }
