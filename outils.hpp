#include <iostream>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <vector>

using namespace std;
const double sig = 0.2; //Valeur de la vol
const double r=0.03;
const double S0=100;
const double K_put=95;
const double tau=1./52;
const double X0_put=1.669;
const double X0_gauss=0;
const double T=0.25;
const double mu = 0.08;
const double sig_inner=5;
const double C_G_10=1.282;
const double C_P_10=0.859;
const double C_G_1=2.326;
const double C_P_1=1.221;
const double C_G_01=3.09;
const double C_P_01=1.39;

double uniforme(void)
{
	return (1.+rand())/(1.+RAND_MAX);
	// On renvoie un nombre aleatoire dans ]0,1]. Au numerateur, on ecrit 1.+rand() pour ne jamais renvoyer la valeur 0 (singularite du logarithme)
}

double gauss(void)
{
	// Methode de Box-Muller pour les Gaussiennes
	return sqrt( -2*log(uniforme()) )  *  cos( 2*M_PI*uniforme() );
}

double* indicatrice_sup(double *x, double y,int z){
    for (int i=0;i<z;i++){
	
	if(x[i]>=y){
	    x[i]=1;
	}else{
	    x[i]=0;
	}
    }
    return x;
}

vector<double> indicatrice_sup_vector(vector <double>x, double y,int z){
    for (int i=0;i<z;i++){
	
	if(x[i]>=y){
	    x[i]=1;
	}else{
	    x[i]=0;
	}
    }
    return x;
}
//Fonction renvoyant la moyenne des valeurs d'un tableau.
double moyenne (double *a,int n){
    double s=0;
    for (int i=0;i<n;i++){
	s+=a[i];
    }
    s=s/n;
    return s;
}

double moyenne_vector (vector<double> a,int n){
    double s=a[0];
    for (int i=1;i<n;i++){
	s+=a[i];
    }
    s=s/n;
    return s;
}

double moyenne_indicatrice_sup(double *x, double y,int z){
    double cpt=0;
    for (int i=0;i<z;i++){
	if(x[i]>=y){
	    cpt+=1;
	}
    }
    return cpt/z;
}

double moyenne_indicatrice_sup_vector(vector<double> x, double y,int z){
    double cpt=0;
    for (int i=0;i<z;i++){
	if(x[i]>=y){
	    cpt+=1;
	}
    }
    return cpt/z;
}

double L_gauss(double X_gauss){
    return -X_gauss;
}

double scenario1_gauss(){
    return gauss();
}

double scenario1_Put(){
    double w=gauss();
    double s=S0*exp((mu-sig*sig*0.5)*tau+sig*sqrt(tau)*w);
    return s;
}

double scenario2_gauss(){
    return gauss();
}

double scenario2_Put(double s1){
    double w=gauss();
    double s=s1*exp((r-(sig*sig)/2)*(T-tau)+sig*sqrt(T-tau)*w);
    return s;
}
        
double Z_bar_gauss(double X_gauss,double Y_gauss){
    double l_gauss= L_gauss(X_gauss);
    return l_gauss+sig_inner*Y_gauss;
}

double Z_bar_Put(double ST){
    double z=X0_put-exp(-r*(T-tau))*max(K_put-ST,0.);
    return z;
}

//Fonction renvoyant l'argmin d'un tableau
int argmin(double* s,int n){
    int i=0;
    double min=s[0];
    for(int j=1;j<n;j++){
	if(min>s[j]){
	    min=s[j];
	    i=j;
	}
    }
    return i;
}

int argmin_vector(vector<double> s,int n){
    int i=0;
    double min=s[0];
    for(int j=1;j<n;j++){
	if(min>s[j]){
	    min=s[j];
	    i=j;
	}
    }
    return i;
}

double f(const double li,double c,double sig,double mi){
	return mi*fabs(li-c)/sig;
}

//Fonction renvoyant la somme des valeurs d'un tableau
double somme(double* x, int n){
    double cpt=x[0];
    for(int i=1;i<n;i++){
	cpt+=x[i];
    }
    return cpt;
}

double somme_vector(vector<double> x, int n){
    double cpt=x[0];
    for(int i=1;i<n;i++){
	cpt+=x[i];
    }
    return cpt;
}

// Approximation de la fonction de distribution cumulative normale (par des séries de Taylor)
// Fonction reprise du document de Bernt Arne Ødegaard
double FDR_Gauss(double z)
{
if (z > 6.0)
return 1.0;  // éviter les valeurs illicites
if (z < -6.0)
return 0.0;
double b1 = 0.31938153;
double b2 = -0.356563782;
double b3 = 1.781477937;
double b4 = -1.821255978;
double b5 = 1.330274429;
double p = 0.2316419;
double c2 = 0.3989423;
double a=fabs(z);
double t = 1.0/(1.0+a*p);
double b = c2*exp((-z)*(z/2.0));
double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
n = 1.0-b*n;
if ( z < 0.0 )
n = 1.0 - n;
return n;
}

int partie_entiere_exces(double c){
    return (int)c;
}

double partie_entiere(double c){
    return floor (c);
}

double Densi_Gauss(double x)
{
    return exp(-(x*x)/2)/sqrt(2*M_PI);
}

double alpha_bar(int const n,double const sig, const vector<double> M,double const c,const vector<double> L){
    double a=0;
    for (int i=0;i<n;i++){
	a+=FDR_Gauss(M[i]*(L[i]-c)/sig);
    }
    a/=n;
    return a;
}




double Biais(double a_chapeau,double a_bar){
    return a_chapeau-a_bar;
}


double Var(double a,int n){
    return (a*(1-a))/n;
}

double min_tab(vector<double> s,int n){
    double min=s[0];
    for(int j=1;j<n;j++){
	if(min>s[j]){
	    min=s[j];
	}
    }
    return min;
}

double* copie_tab(double* tab1,double* tab2, int taille1,int taille2){
    if(taille1==taille2){
	for (int i=0;i<taille2;i++){
	    tab2[i]=tab1[i];
	}
    }
    else{
	for (int i=0;i<taille2;i++){
	    if(i<taille1){
		tab2[i]=tab1[i];
	    }
	    else{
		tab2[i]=0;
	    }
	}
    }
    return tab2;
}

double* change_dim(double* M,int old_taille, int new_taille){
    if(old_taille<new_taille){
    double* M_bis=new double [new_taille];
    M_bis=copie_tab(M,M_bis,old_taille,new_taille);
    delete []M;
    double* Z=new double [new_taille];
    M=copie_tab(M_bis,Z,new_taille,new_taille);
    delete []M_bis;
    return Z;}
    else {
	return M;
    }
}

