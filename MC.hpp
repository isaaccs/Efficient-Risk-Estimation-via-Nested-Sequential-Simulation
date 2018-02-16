#include "outils.hpp"

template <int m,int n>
class scenario{
    private: double* s;
    
    public: 
    scenario();
    //scenario (scenario<n> const & p);//constructeur par copie
    ~scenario (){delete [] s;}
    double& operator()(int j,int i);
    double operator()(int j,int i) const;
    
};



template <int m,int n>
scenario<m, n>::scenario(){
	    s=new double[n*m];
}

template <int m,int n>
double& scenario<m, n>::operator()(int j,int i){
    return s[j*n+i];
}

template <int m,int n>
double scenario<m, n>::operator()(int j,int i) const{
    return s[j*n+i];
}

////////////////////////////////////////////////
////////////Uniforme avec template/////////////
//////////////////////////////////////////////

//////////////////////////////////////////////
////////////Version Gaussienne///////////////
/////////////////////////////////////////////

template <int m,int n>
void Uniforme_gauss (scenario<m,n> & s,double c){
    
    double a;
    double *L=new double [n];
    double *M=new double [n];
    
    scenario<m-1,n> Z;
    for (int i=0;i<n;i++){
	//On simule un premier scénario jusqu'au temps tau.
	s(0,i)=scenario1_gauss();
	L[i]=0;
	for (int j=1;j<m;j++){
	    //Puis un deuxième jusqu'au temps T.
	    s(j,i)=scenario2_gauss();
	    //on calcule les pertes de ce dexuieme scénario.
	    Z(j-1,i)=Z_bar_gauss(s(0,i),s(j,i));
	    L[i]+=Z(j-1,i);
	}
	//Puis on fait la moyenne.
	L[i]/=m;
    }
    //calcul de alpha chapeau
    a=moyenne_indicatrice_sup(L,c,n);
    //calcul de la variance
    double var=Var(a,n);
    //calcul du biais.
    double biais;
        if(c==C_G_10){
	    biais=pow((a-0.1),2);
    }else{
	if(c==C_G_1){
	    biais=pow((a-0.01),2);
	}else{
	    biais=pow((a-0.001),2);
	}
    }
    //calcul de la MSE
    double MSE=var+biais;
    cout<<"variance Gauss "<<var<<endl;
    cout<<"biais Gauss "<<biais<<endl;
    cout<<"MSE "<<MSE<<endl;
    cout<<"alpha Gauss "<<a*100<<endl;
    cout<<"\n"<<endl;
    delete [] L;
}

template <int m,int n>
void Uniforme_Put(scenario<m,n> & s,double c){
    double a;
    double *L=new double [n];
    scenario<m-1,n> Z;
    for (int i=0;i<n;i++){
	s(0,i)=scenario1_Put();;
	L[i]=0;
	for (int j=1;j<m;j++){ 
	    s(j,i)=scenario2_Put(s(0,i));
	    Z(j-1,i)=Z_bar_Put(s(j,i));
	    L[i]+=Z(j-1,i);
	}
    }
    for(int i=0;i<n;i++){
	L[i]=L[i]/m;
    }
    a=moyenne_indicatrice_sup(L,c,n);
    double var=Var(a,n);
    double biais;
    if(c==C_P_10){
	biais=pow((a-0.1),2);
    }else{
	if(c==C_P_1){
	    biais=pow((a-0.01),2);
	}else{
	    biais=pow((a-0.001),2);
	}
    }
    double MSE=var+biais;
    cout<<"variance PUT "<<var<<endl;
    cout<<"biais PUT "<<biais<<endl;
    cout<<"MSE PUT "<<MSE<<endl;
    cout<<"alpha Put "<<a*100<<endl;
    cout<<"\n"<<endl;
    delete []L;
}


///////////////////////////////////////////////////////
////////////// Fin Uniforme Template///////////////////
///////////////////////////////////////////////////////






////////////////////////////////
/////Sequential////////////////
///////////////////////////////


///////Version Gauss//////

//Version 10%
 double Sequential_gauss_10 (int m,int n,int m0){
    double c=C_G_10;
    double a;
    double* S=new double [n];
    double *L=new double [n];
    double *Z=new double [n];
    double *M=new double [n];
    double S2;
    for (int i=0;i<n;i++){
	S[i]=scenario1_gauss();
	L[i]=0;
	Z[i]=0;
	for (int j=0;j<m0;j++){
	    S2=scenario2_gauss();
	    Z[i]+=Z_bar_gauss(S[i],S2);
	}
	M[i]=m0;
	L[i]=Z[i]/M[i];
    }
    int idx=0;
    double *y=new double [n];
    for(int i=0;i<n;i++){
	y[i]=f(L[i],c,sig_inner,M[i]);
    }
    while (somme(M,n)<n*m){
	if(fmod(somme(M,n)*100/(n*m),10.)==0.){
	    cout<<somme(M,n)*100/(n*m)<<"%"<<endl;
	}
	y[idx]=f(L[idx],c,sig_inner,M[idx]);
	idx=argmin(y,n);
	S2=scenario2_gauss();
	Z[idx]+=Z_bar_gauss(S[idx],S2);
	M[idx]+=1;
	L[idx]=Z[idx]/M[idx];
    }
    a=moyenne_indicatrice_sup(L,c,n);
    double var=Var(a,n);
    double biais=pow((a-0.10),2);
    double MSE=var+biais;
    cout<<"variance Seq Gauss "<<var<<endl;
    cout<<"biais Seq Gauss "<<biais<<endl;
    cout<<"MSE Seq Gauss "<<MSE<<endl;
    cout<<"alpha Seq Gauss "<<a*100<<endl;
    cout<<moyenne(M,n)<<endl;
    cout<<"\n"<<endl;
    delete [] L;
    return a;
}

//Version 1%

double Sequential_gauss_1 (int m,int n,int m0){
    ofstream sortie;
    sortie.open("gauss");
    double c=C_G_1;
    double a;
    double* S=new double [n];
    double *L=new double [n];
    double *Z=new double [n];
    double *M=new double [n];
    double S2;
    for (int i=0;i<n;i++){
	S[i]=scenario1_gauss();
	L[i]=0;
	Z[i]=0;
	for (int j=0;j<m0;j++){
	    S2=scenario2_gauss();
	    Z[i]+=Z_bar_gauss(S[i],S2);
	}
	M[i]=m0;
	L[i]=Z[i]/M[i];
    }
    int idx=0;
    double *y=new double [n];
    for(int i=0;i<n;i++){
	y[i]=f(L[i],c,sig_inner,M[i]);
    }
    while (somme(M,n)<n*m){
	if(fmod(somme(M,n)*100/(n*m),10.)==0.){
	    cout<<somme(M,n)*100/(n*m)<<"%"<<endl;
	}
	y[idx]=f(L[idx],c,sig_inner,M[idx]);
	idx=argmin(y,n);
	S2=scenario2_gauss();
	Z[idx]+=Z_bar_gauss(S[idx],S2);
	M[idx]+=1;
	L[idx]=Z[idx]/M[idx];
    }
    for (int i=0;i<n;i++){
	sortie<<L[i]<<"\t"<<M[i]<<endl;
    }
    a=moyenne_indicatrice_sup(L,c,n);
    double var=Var(a,n);
    double biais=pow((a-0.010),2);
    double MSE=var+biais;
    cout<<"variance Seq Gauss"<<var<<endl;
    cout<<"biais Seq Gauss"<<biais<<endl;
    cout<<"MSE Seq Gauss "<<MSE<<endl;
    cout<<"alpha Seq Gauss "<<a*100<<endl;
    cout<<moyenne(M,n)<<endl;
    cout<<"\n"<<endl;
    delete [] L;
    return a;
    sortie.close();
}

//Version 0.1%

double Sequential_gauss_01 (int m,int n,int m0){
    double c=C_G_01;
    double a;
    double* S=new double [n];
    double *L=new double [n];
    double *Z=new double [n];
    double *M=new double [n];
    double S2;
    for (int i=0;i<n;i++){
	S[i]=scenario1_gauss();
	L[i]=0;
	Z[i]=0;
	for (int j=0;j<m0;j++){
	    S2=scenario2_gauss();
	    Z[i]+=Z_bar_gauss(S[i],S2);
	}
	M[i]=m0;
	L[i]=Z[i]/M[i];
    }
    int idx=0;
    double *y=new double [n];
    for(int i=0;i<n;i++){
	y[i]=f(L[i],c,sig_inner,M[i]);
    }
    while (somme(M,n)<n*m){
	if(fmod(somme(M,n)*100/(n*m),10.)==0.){
	    cout<<somme(M,n)*100/(n*m)<<"%"<<endl;
	}
	y[idx]=f(L[idx],c,sig_inner,M[idx]);
	idx=argmin(y,n);
	S2=scenario2_gauss();
	Z[idx]+=Z_bar_gauss(S[idx],S2);
	M[idx]+=1;
	L[idx]=Z[idx]/M[idx];
    }
    a=moyenne_indicatrice_sup(L,c,n);
    double var=Var(a,n);
    double biais=pow((a-0.0010),2);
    double MSE=var+biais;
    cout<<"variance Seq Gauss "<<var<<endl;
    cout<<"biais Seq Gauss "<<biais<<endl;
    cout<<"MSE Seq Gauss "<<MSE<<endl;
    cout<<"alpha Seq Gauss"<<a*100<<endl;
    cout<<moyenne(M,n)<<endl;
    cout<<"\n"<<endl;
    delete [] L;
    return a;
}










//////Version Put///////

//Version 10%
double Sequential_PUT_10 (int m,int n,int m0){
    double c=C_P_10;
    double a;
    double* S=new double [n];
    double *L=new double [n];
    double *Z=new double [n];
    double *M=new double [n];
    double S2;
    for (int i=0;i<n;i++){
	S[i]=scenario1_Put();
	L[i]=0;
	Z[i]=0;
	for (int j=0;j<m0;j++){
	    S2=scenario2_Put(S[i]);
	    Z[i]+=Z_bar_Put(S2);
	}
	M[i]=m0;
	L[i]=Z[i]/M[i];
    }
    int idx=0;
    double *y=new double [n];
    for(int i=0;i<n;i++){
	y[i]=f(L[i],c,sig,M[i]);
    }
    while (somme(M,n)<n*m){
	if(fmod(somme(M,n)*100/(n*m),10.)==0.){
	    cout<<somme(M,n)*100/(n*m)<<"%"<<endl;
	}
	y[idx]=f(L[idx],c,sig,M[idx]);
	idx=argmin(y,n);
	double S2=scenario2_Put(S[idx]);
	Z[idx]+=Z_bar_Put(S2);
	M[idx]+=1;
	L[idx]=Z[idx]/M[idx];
    }
    a=moyenne_indicatrice_sup(L,c,n);
    double var=Var(a,n);
    double biais=pow((a-0.10),2);
    double MSE=var+biais;
    cout<<"variance PUT Seq "<<var<<endl;
    cout<<"biais PUT Seq "<<biais<<endl;
    cout<<"MSE PUT Seq "<<MSE<<endl;
    cout<<"alpha PUT Seq "<<a*100<<endl;
    cout<<moyenne(M,n)<<endl;
    cout<<"\n"<<endl;
    delete [] L;
    return a;
}

//Version 1%
double Sequential_PUT_1 (int m,int n,int m0){
    ofstream sortie;
    sortie.open("seq");
    double c=C_P_1;
    double a;
    double* S=new double [n];
    double *L=new double [n];
    double *Z=new double [n];
    double *M=new double [n];
    double S2;
    for (int i=0;i<n;i++){
	S[i]=scenario1_Put();
	L[i]=0;
	Z[i]=0;
	for (int j=0;j<m0;j++){
	    S2=scenario2_Put(S[i]);
	    Z[i]+=Z_bar_Put(S2);
	}
	M[i]=m0;
	L[i]=Z[i]/M[i];
    }
    int idx=0;
    double *y=new double [n];
    for(int i=0;i<n;i++){
	y[i]=f(L[i],c,sig,M[i]);
    }
    while (somme(M,n)<n*m){
	if(fmod(somme(M,n)*100/(n*m),10.)==0.){
	    cout<<somme(M,n)*100/(n*m)<<"%"<<endl;
	}
	y[idx]=f(L[idx],c,sig,M[idx]);
	idx=argmin(y,n);
	double S2=scenario2_Put(S[idx]);
	Z[idx]+=Z_bar_Put(S2);
	M[idx]+=1;
	L[idx]=Z[idx]/M[idx];
    }
    for (int i=0;i<n;i++){
	sortie<<M[i]<<"\t"<<L[i]<<endl;
    }
    a=moyenne_indicatrice_sup(L,c,n);
    double var=Var(a,n);
    double biais=pow((a-0.010),2);
    double MSE=var+biais;
    cout<<"variance PUT Seq "<<var<<endl;
    cout<<"biais PUT Seq "<<biais<<endl;
    cout<<"MSE PUT Seq "<<MSE<<endl;
    cout<<"alpha PUT Seq "<<a*100<<endl;
    cout<<moyenne(M,n)<<endl;
    cout<<"\n"<<endl;
    sortie.close();
    delete [] L;
    return a;
}

//Version 0.1%
double Sequential_PUT_01 (int m,int n,int m0){
    double c=C_P_01;
    double a;
    double* S=new double [n];
    double *L=new double [n];
    double *Z=new double [n];
    double *M=new double [n];
    double S2;
    for (int i=0;i<n;i++){
	S[i]=scenario1_Put();
	L[i]=0;
	Z[i]=0;
	for (int j=0;j<m0;j++){
	    S2=scenario2_Put(S[i]);
	    Z[i]+=Z_bar_Put(S2);
	}
	M[i]=m0;
	L[i]=Z[i]/M[i];
    }
    int idx=0;
    double *y=new double [n];
    for(int i=0;i<n;i++){
	y[i]=f(L[i],c,sig,M[i]);
    }
    while (somme(M,n)<n*m){
	if(fmod(somme(M,n)*100/(n*m),10.)==0.){
	    cout<<somme(M,n)*100/(n*m)<<"%"<<endl;
	}
	y[idx]=f(L[idx],c,sig,M[idx]);
	idx=argmin(y,n);
	double S2=scenario2_Put(S[idx]);
	Z[idx]+=Z_bar_Put(S2);
	M[idx]+=1;
	L[idx]=Z[idx]/M[idx];
    }
    a=moyenne_indicatrice_sup(L,c,n);
    double var=Var(a,n);
    double biais=pow((a-0.0010),2);
    double MSE=var+biais;
    cout<<"variance PUT Seq "<<var<<endl;
    cout<<"biais PUT Seq "<<biais<<endl;
    cout<<"MSE PUT Seq "<<MSE<<endl;
    cout<<"alpha PUT Seq "<<a*100<<endl;
    cout<<moyenne(M,n)<<endl;
    cout<<"\n"<<endl;
    delete [] L;
    return a;
}

////////////////////////////////////////
//////////// Fin Sequential/////////////
////////////////////////////////////////







void Threshold_gauss (int n,double gamma,double c){
    double a;
    double* S=new double [n];
    double *L=new double [n];
    double *Z=new double [n];
    double *M=new double [n];
    double S2;
    for (int i=0;i<n;i++){
	cout<<i/n<<"%"<<endl;
	S[i]=scenario1_gauss();
	L[i]=0;
	Z[i]=0;
	M[i]=0;
	do{
	    S2=scenario2_gauss();
	    Z[i]+=Z_bar_gauss(S[i],S2);
	    M[i]+=1;
	    L[i]=Z[i]/M[i];
	    cout<<"while "<<f(L[i],c,sig_inner,M[i])<<" %"<<endl;
	}
	while(!(f(L[i],c,sig_inner,M[i])>=gamma));
    }
    a=moyenne_indicatrice_sup(L,c,n);
    double var=Var(a,n);
    double biais=pow((a-0.1),2);
    double MSE=var+biais;
    cout<<"variance Threshold Gauss "<<var<<endl;
    cout<<"biais Threshold Gauss "<<biais<<endl;
    cout<<"MSE Threshold Gauss "<<MSE<<endl;
    cout<<"alpha Threshold Gauss "<<a*100<<endl;
    cout<<"\n"<<endl;
    delete []L;
    
}







void Adaptative_Gauss(int m0,int n0,double tau_e,double k,double c){
    vector<double> L(n0);
    vector<double> Z(n0);
    vector<double> S(n0);
    double S2=0;
    double a=0;
    
    for (int i=0;i<n0;i++){
	S[i]=scenario1_gauss();
	}
	
    int n=n0;
    
    vector<double> y(n);
    vector<double> M(n);

    for (int i=0;i<n0;i++){
	Z[i]=0;
	M[i]=0;
	L[i]=0;
	for (int j=0;j<m0;j++){
	    S2=scenario2_gauss();
	    Z[i]+=Z_bar_gauss(S[i],S2);
	}
	M[i]=m0;
	L[i]=Z[i]/M[i];
	y[i]=f(L[i],c,sig_inner,M[i]);
    }
    
    double m_bar=0;
    double m_bar_4=0;
    double V=0;
    double B=0;
    double B_2=0;
    int n_prime=0;
    double frac1=0;
    double frac=0;
    int idx=0;
    double a_bar=0;
    a=moyenne_indicatrice_sup_vector(L,c,n);
    cout<<"aaaa "<<a<<endl;
    for(int l=1;l<partie_entiere_exces(k/tau_e);l++){
	a_bar=alpha_bar(n,sig_inner,M,c,L);
	cout<<" a_bar "<<a_bar<<endl;
	a=moyenne_indicatrice_sup_vector(L,c,n);
	cout<<" a "<<a<<endl;
	m_bar=moyenne_vector(M,n);
	cout<<" M bar "<<m_bar<<endl;
	m_bar_4=pow(m_bar,4);
	V=Var(a_bar,n);
	cout<<" V "<<V<<endl;
	B=Biais(a,a_bar);
	B_2=pow(B,2);
	cout<<" B "<<B<<endl;
	frac1=(V*n)/(4*B_2*m_bar_4);
	frac=pow((pow(m_bar*n+tau_e,4))*frac1,0.2);
	cout<<"P4 "<<pow(m_bar*n+tau_e,4)*frac1<<endl;
	cout<<"frac "<<frac<<endl;
	cout<<" n "<<n<<endl;
	n_prime=partie_entiere(min(max(frac,(double)n),(double)(n+tau_e)));
	cout<<n_prime<<endl;
	/*M=change_dim(M,n,n_prime);
	L=change_dim(L,n,n_prime);
	S=change_dim(S,n,n_prime);
	Z=change_dim(Z,n,n_prime);
	y=change_dim(y,n,n_prime);*/
	if(n<n_prime){
	L.resize(n_prime);
	M.resize(n_prime);
	S.resize(n_prime);
	Z.resize(n_prime);
	y.resize(n_prime);
	for(int i=n;i<n_prime;i++){
	    S[i]=scenario1_gauss();
	    M[i]=0;
	    L[i]=0;
	    y[i]=f(L[i],c,sig_inner,M[i]);
	}
	}
	n=n_prime;
	cout<<" n "<<n<<" n prime "<<n_prime<<endl;
	
	
	while(somme_vector(M,n)<l*tau_e){
	    if(min_tab(M,n)<m0){
		idx=argmin_vector(M,n);
	    }else{
		idx=argmin_vector(y,n);
	    }
	    double G=gauss();
	    Z[idx]+=Z_bar_gauss(S[idx],G);
	    M[idx]+=1;
	    L[idx]=Z[idx]/M[idx];
	    y[idx]=f(L[idx],c,sig_inner,M[idx]);
	}
    }

    a=moyenne_indicatrice_sup_vector(L,c,n);
    double var=(a*(1-a))/n;
    double biais=pow((a-0.1),2);
    double MSE=var+biais;
    cout<<"variance Gauss Seq "<<var<<endl;
    cout<<"biais Gauss Seq "<<biais<<endl;
    cout<<"MSE Seq "<<MSE<<endl;
    cout<<"alpha Gauss Seq "<<a*100<<endl;
    cout<<moyenne_vector(M,n)<<endl;
    cout<<"\n"<<endl;
}
