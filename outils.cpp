#include "outils.hpp"

//Ce fichier a pour but de tester les fonctions implémenté dans le fichier outils.hpp


int main(){
    //Test de la fonction argmin
    double t[10]={51,22,34,0,-45,6,44,34,-99.52389,0};
    cout<<argmin(t,10)<<fixed<<t[argmin(t,10)]<<endl;
    cout<<min_tab(t,10)<<endl;
    
    cout<<endl;
    
   //Test de la fonction uniform
    double u[10];
    for (int i = 0 ; i<10 ; i++){
		u[i]=uniforme();
		cout<<u[i]<<endl;
	}
    
    cout<<endl;
    //Test de la fonction gauss
    
    double g[10];
    for(int j = 0; j<10 ; j++){
		g[j]=gauss();
		cout<<g[j]<<endl;
    }
    
     cout<<endl;
     //Test de la fonction indicatrice_sup
     
	 double q[2]={0.4,0.6};
	 double * l=new double[2];
     l=indicatrice_sup(q,0.5,2);
     for(int i=0; i<2; i++){
		 cout<<l[i]<<endl;
	 }
   
     //Test de la fonction moyenne
     
     double o[5]={1,1,1,1,1};
     double h;
     h=moyenne(o,5);
     cout<<h<<endl;
     
     cout<<endl;
     //Test de la fonction moyenne_indicatrice_sup
     
     double v;
     v=moyenne_indicatrice_sup(u,0.5,10);
     cout<<v<<endl;
     
     cout<<endl;
     //Test de la fonction L_gauss
     double m;
     m=L_gauss(1);
     cout<<m<<endl;
     
     cout<<endl;
     //Test de la fonction scenario1_gauss
     
     double scg[5];
     for (int i = 0 ; i<5 ; i++){
	scg[i]=scenario1_gauss();
	cout<<scg[i]<<endl;
    }
	 
    cout<<endl;
    //Test de la fonction scenario1_Put()
	 
    double scput[5];
    for (int i = 0 ; i<5 ; i++){
	scput[i]=scenario1_Put();
	cout<<scput[i]<<endl;
    }
	 
    cout<<endl;
    //Test de la fonction scenario2_gauss
	 
    double scg2[5];
    for (int i = 0 ; i<5 ; i++){
	scg2[i]=scenario2_gauss();
	cout<<scg2[i]<<endl;
    }
	
    cout<<endl;
    //Test de la fonction scenario2_Put
	 
    double scput2[5];
     for (int i = 0 ; i<5 ; i++){
	scput2[i]=scenario2_Put(scput[i]);
	cout<<scput2[i]<<endl;
    }
	
    cout<<endl;
    //Test de la fonction Z_bar_gauss
	
    double y;
    double u1=gauss();
    double u2=gauss();
    cout<<"u1="<<u1<< " " << "u2="<< u2<< endl;
    y = Z_bar_gauss(u1,u2);
    cout<< y << endl;
	
    cout<<endl;
    //Test de la fonction Z_bar_Put

    double ST1=100;
    double ST2=90;
    double test1=0;
    double test2=0;
    test1 = Z_bar_Put(ST1);
    test2 = Z_bar_Put(ST2);
    cout<<"test1="<< test1 << " " << "test2="<< test2 << endl;
	  
    cout<<endl;
    //Test de la fonction f
	  
    double test3;
    test3=f(1,2,0.5,3);
    cout<<test3<<endl;
	  
    cout<<endl;
    //Test de la fonction somme
	  
    double tab[5]={1,1,2,3,1};
    double s=0;
    s=somme(tab,5);
    cout<< s <<endl;
	  
    cout<<endl;
	   
    //Test de la fonction de répartion d'une Loi Normale Centrée Réduite.
    double x=0;
    cout<<FDR_Gauss(x)<<endl;
    double qs=5;
    cout<<FDR_Gauss(qs)<<endl;
    double z=-5;
    cout<<FDR_Gauss(z)<<endl;
    double xx=-1;
    cout<<FDR_Gauss(xx)<<endl;
    double yy=1;
    cout<<FDR_Gauss(yy)<<endl;
    
    //Test de alpha_bar
    double* MM= new double [5];
    double* LL= new double [5];
    for(int i=0;i<5;i++){
	MM[i]=0;
	LL[i]=0;
    }
    cout<<alpha_bar(5,1,MM,1,LL)<<endl;
    
    //Test de la fonction copie_tab et redim_tab en même temps.
    double *M=new double [5];
    for (int i=0;i<5;i++){
        M[i]=i;
    }
    M=change_dim(M,5,10);
    for(int i=0;i<10;i++){
        cout<<M[i]<<endl;
    }
    double *Z=new double [5];
    for (int i=0;i<5;i++){
        Z[i]=i;
    }
    
    Z=change_dim(Z,5,10);
    for(int i=0;i<10;i++){
        cout<<Z[i]<<endl;
    }
    delete []Z,M;
    
    cout<<partie_entiere(2.7)<<endl;
    cout<<partie_entiere(2)<<endl;
    cout<<partie_entiere(-2.7)<<endl;
    cout<<partie_entiere_exces(2.7)<<endl;
    cout<<partie_entiere_exces(2)<<endl;
    cout<<partie_entiere_exces(-2.7)<<endl;
    return 0;
}
