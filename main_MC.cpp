#include "MC.hpp"
#include <ctime>


const int m_1_2=159;
const int n_1_2=25199;
const int k=4000000;
const int tau_e=100000;
const int m0=2;
const int n0=500;

int main(){
	
    //SCENARIO GAUSSIEN
    scenario<m_1_2,n_1_2> s;
    scenario<889,4499> s5;
    scenario<m_1_2,n_1_2> s7;
    scenario<786,5089> s9;
    scenario<m_1_2,n_1_2> s11;
    scenario<514,7788> s13;
     
    //SCENARIO PUT
    scenario<m_1_2,n_1_2> s1;
    scenario<785,5095> s6;
    scenario<m_1_2,n_1_2> s8;
    scenario<1273,3143> s10;
    scenario<m_1_2,n_1_2> s12;
    scenario<1556,2570> s14;
    
    //Timer
    clock_t start;
    double duration;
    
    //Uniforme
    cout<<"uniform 10% 1/3 2/3 "<<endl;
    //UNIFORM 10% 1/3 2/3
    start = std::clock();
    Uniforme_gauss(s,C_G_10);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Uniforme Gauss 1/3 2/3 10% "<< duration <<'\n';
    cout<<endl;
    
    start = std::clock();
    Uniforme_Put(s1,C_P_10);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Uniforme Put 1/3 2/3 10% "<< duration <<'\n';
    cout<<endl;   
    
    cout<<"uniform 10% optimal "<<endl;
    //UNIFORM 10% OPTIMAL
    start = std::clock();
    Uniforme_gauss(s5,C_G_10);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Uniforme Gauss Optimal 10% "<< duration <<'\n';
    cout<<endl;
    
    start = std::clock();
    Uniforme_Put(s6,C_P_10);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Uniforme Put Optimal 0.1% "<< duration <<'\n';
    cout<<endl;
    
    cout<<"uniform 1% 1/3 2/3 "<<endl;
    //UNIFORM 1% 1/3 2/3
    start = std::clock();
    Uniforme_gauss(s7,C_G_1);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Uniforme Gauss 1/3 2/3 1% "<< duration <<'\n';
    cout<<endl;
    
    start = std::clock();
    Uniforme_Put(s8,C_P_1);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Uniforme Put 1/3 2/3 1% "<< duration <<'\n';
    cout<<endl;
   
    cout<<"uniform 1% optimal "<<endl;
    //UNIFORM 1% OPTIMAL
    start = std::clock();
    Uniforme_gauss(s9,C_G_1);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Uniforme Gauss Optimal 1% "<< duration <<'\n';
    cout<<endl;
    
    start = std::clock();
    Uniforme_Put(s10,C_P_1);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Uniforme Put Optimal 1% "<< duration <<'\n';
    cout<<endl;

    cout<<"uniform 0.10% 1/3 2/3 "<<endl;
    //UNIFORM 0.1% 1/3 2/3
    start = std::clock();
    Uniforme_gauss(s11,C_G_01);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Uniforme Gauss 1/3 2/3 0.1% "<< duration <<'\n';
    cout<<endl;
    
    start = std::clock();
    Uniforme_Put(s12,C_P_01);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Uniforme Put 1/3 2/3 0.1% "<< duration <<'\n';
    cout<<endl;
    
    cout<<"uniform 0.10% optimal "<<endl;
    //UNIFORM 0.1% OPTIMAL
    start = std::clock();
    Uniforme_gauss(s13,C_G_01);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Uniforme Gauss Optimal 0.1% "<< duration <<'\n';
    cout<<endl;
    
    start = std::clock();
    Uniforme_Put(s14,C_P_01);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Uniforme Put Optimal 0.1% "<< duration <<'\n';
    cout<<endl;
    
    //Sequential//
    
    start = std::clock();
    Sequential_PUT_01(151,26508,m0);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Sequential Put 0.1% "<< duration <<'\n';
    cout<<endl;
    
    start = std::clock();
    Sequential_PUT_1(205,19558,m0);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Sequential Put 1% "<< duration <<'\n';
    cout<<endl;
    
    start = std::clock();
    Sequential_PUT_10(323,12395,m0);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Sequential Put 10% "<< duration <<'\n';
    cout<<endl;
    
    start = std::clock();
    Sequential_gauss_01(71,56686,m0);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Sequential Gauss 0.1% "<< duration <<'\n';
    cout<<endl;
    
    start = std::clock();
    Sequential_gauss_1(130,30860,m0);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Sequential Gauss 1%"<< duration <<'\n';
    cout<<endl;
    
    start = std::clock();
    Sequential_gauss_10(323,12395,m0);  
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: Temps du Sequential Gauss 10%"<< duration <<'\n';
    cout<<endl;
    
    Threshold_gauss(10000,0.2,C_G_01);
    Adaptative_Gauss(m0,n0,tau_e,k,C_G_10);
    return 0;
}
