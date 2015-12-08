#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

using namespace std;

void f(double* y,double* k, int& N);							//calc f(...) with given parameters

int main (){
  int N=3;									//length of Vector y and f(y,t)
  double y[N];						//Vector y=(x,y,z)
  double k1[N];						//Vector k1 to k4
  double k2[N];
  double k3[N];
  double k4[N];
  double temp[N];
  const string file = "data.txt";
  const int tmin = 0;
  const int tmax = 100;
  //to compare with the plots: dt1=0.1, dt2=0.01, dt3=0.001
  double dt= 0.1;								//step size
  double t=0;
  
  ofstream out(file.c_str());
  
  
  for(int i=0;i<N;i++){								//initial condition
    y[i]=1;
  }
  //out<<t<<"\t"<<y[0]<<"\t"<<y[1]<<"\t"<<y[2]<<endl;
  /*
   * 
   * calculate k1,k2,k3,k4 for y_i
   * 
   * 
   */
  while(t<tmax){
    out<<t<<"\t"<<y[0]<<"\t"<<y[1]<<"\t"<<y[2]<<endl;					//write y in file data; 1st column = x, 2nd column = y, 3rd column = z;
    //calc k1
    f(y,k1,N);										//k1=f(y)
    
    //calc k2
    for(int i=0;i<N;i++){
      temp[i]=y[i]+0.5*k1[i]*dt;							//calc argument of f(...)
    }
    f(temp,k2,N);
    
    //calc k3
    for(int k=0;k<N;k++){
      temp[k]=y[k]+0.5*k2[k]*dt;
    }
    f(temp,k3,N);
    
    //calc k4
    for(int j=0;j<N;j++){
      temp[j]=y[j]+k3[j]*dt;
    }
    f(temp,k4,N);
    
    //calc y_n+1
    for(int l=0;l<N;l++){
      y[l]=y[l]+dt/6*(k1[l]+2*k2[l]+2*k3[l]+k4[l]);
    }
    t+=dt;
  }
  
  
  out.close();
  return 0; 
}


/*
 * 
 * 
 * Functions 
 * 
 * 
 * 
 */

void f(double* y, double* k, int& N){			//from given Readme.md one can see f(y,t)= Lorenz model, when y=(x,y,z)
  double* const parameter  = new double [N];		//dynamic allocation of a vector with 3 entries vor the parameter a,b and c
	      parameter[0] = 10;
	      parameter[1] = 28;
	      parameter[2] = 8.0/3;
  
  k[0] = parameter[0]*(y[1]-y[0]);		//Lorenz model 
  k[1] = y[0]*(parameter[1]-y[2])-y[1];
  k[2] = y[0]*y[1]-parameter[2]*y[2];
  delete[] parameter;				//free space
							//after function f is closed taken memory should be the same size
}