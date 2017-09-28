// Mi Primer MOnte Carlo
#include<iostream>
#include<cmath>
#include "Random64.h"
using namespace std;

const int N=1000;

class Agente{
private:
  double dinero;
public: 
  void Inicie(double dinero0){dinero=dinero0;};
  double Getdinero(void){return dinero;};
  friend class Mercado;
};

class Mercado{
private:

public:
  void HagaTransaccionEntre(Agente & Agente1,Agente & Agente2,Crandom & ran64);
};
void Mercado::HagaTransaccionEntre(Agente & Agente1,Agente & Agente2,Crandom & ran64){
  double GananciaAgente1,GananciaAgente2,Resultado;
  double NuevoDineroAgente1,NuevoDineroAgente2;
  Resultado=2*ran64.r();
  GananciaAgente1=Resultado-1;   GananciaAgente2=(2-Resultado)-1; 
  NuevoDineroAgente1=Agente1.dinero+GananciaAgente1;
  NuevoDineroAgente2=Agente2.dinero+GananciaAgente2;
  if(NuevoDineroAgente1>0 &&NuevoDineroAgente2>0){
    Agente1.dinero=NuevoDineroAgente1;  Agente2.dinero=NuevoDineroAgente2;
  }
}


int main(void){
  Mercado Corabastos;
  Agente  Paisano[N];
  Crandom ran64(1);
  int i,j,t,Ntransacciones=N*10000;

  //Inicie todos los agentes con 10 euros
  for(i=0;i<N;i++) Paisano[i].Inicie(10);

  for(t=0;t<Ntransacciones;t++){
    //Escojo dos agentes al azar;
    i=(int) N*ran64.r();
    j=(int) (N-1)*ran64.r(); if(j>=i) j++;
    //Los pongo a interacturar;
    Corabastos.HagaTransaccionEntre(Paisano[i],Paisano[j],ran64);
  }
  //Imprimo los dineros de todos;
  for(i=0;i<N;i++) cout<<Paisano[i].Getdinero()<<endl;


  return 0;
}
