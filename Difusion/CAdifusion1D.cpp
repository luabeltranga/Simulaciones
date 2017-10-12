#include <iostream>
#include <cmath>
#include "Random64.h"

const int Lx=200;
const double p=0.5;

class LatticeGas{
private:
  int V[2];//V[0] derecha, V[1] izquierda
  int n[Lx][2],nnew[Lx][2]; //n[ix][i]
public:
  LatticeGas(void);
  void Inicie(Crandom & ran64);
  void Show(int copia);
  void Colisione(Crandom & ran64);
  void Adveccione(void);
  int DondeEstaLaBolita(void);
};
LatticeGas::LatticeGas(void){
  V[0]=1;
  V[1]=-1;
}
void LatticeGas::Inicie(Crandom &ran64){
  for(int ix = 0;ix<Lx;ix++){
    for(int ii = 0;ii<2;ii++){
      n[ix][ii]=0;nnew[ix][ii]=0;
    }
  }
  if(ran64.r()<0.5)  n[Lx/2][0]=1; else n[Lx/2][1]=1; 
}
void LatticeGas::Show(int copia){
  for(int ii = 0;ii<2;ii++){
    for(int ix = 0;ix<Lx;ix++){
      if(copia==0){
	std::cout<<n[ix][ii];     
      }
      else{
	std::cout<<nnew[ix][ii];     
      }
    }
    std::cout<<std::endl;
  }
}

void LatticeGas::Colisione(Crandom & ran64){
  for(int ix = 0;ix<Lx;ix++){
    if(ran64.r()<p){//dejelo quieto
      for(int ii =0;ii<2;ii++) nnew[ix][ii]=n[ix][ii];
    }
    else{//voltee
      for(int ii =0;ii<2;ii++) nnew[ix][ii]=n[ix][(ii+1)%2];
    }
    
  }
}

void LatticeGas::Adveccione(void){
  for(int ix = 0;ix<Lx;ix++){
    for(int ii = 0;ii<2;ii++){
      n[(ix+V[ii]+Lx)%Lx][ii]=nnew[ix][ii];
    }
 }
}

int LatticeGas::DondeEstaLaBolita(void){
  int ix=0;
  while(n[ix][0]+n[ix][1]==0){
    ix++;
  }
  return ix;
}

//----------------Funciones Globales--------

const int N=1000;

double sigma2(LatticeGas *Difusion){
  double Xprom=0,sigma2=0;
  for(int ii =0;ii<N;ii++){
    Xprom+=Difusion[ii].DondeEstaLaBolita();
  }
  Xprom/=N;

  for(int ii =0;ii<N;ii++){
    sigma2 +=std::pow(Difusion[ii].DondeEstaLaBolita()-Xprom,2);
  }
  sigma2/=(N-1);
  return sigma2;
}

int main(void){
  LatticeGas Difusion[N];
  Crandom ran64(1012);

  int tmax=400;
  //inicie todo
  for(int ii=0;ii<N;ii++) Difusion[ii].Inicie(ran64); 

  //corra
  for(int t = 0; t<tmax;t++){
    std::cout<<t<<"  "<<sigma2(Difusion)<<std::endl;
    for(int ii=0;ii<N;ii++) Difusion[ii].Colisione(ran64);
    for(int ii=0;ii<N;ii++) Difusion[ii].Adveccione(); 
  }
  
  return 0;
}
