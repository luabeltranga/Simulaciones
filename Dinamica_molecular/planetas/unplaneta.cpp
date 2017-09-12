#include <iostream>
#include <cmath>

const double GM=1;

class Cuerpo;


class Cuerpo{
private:
  double x,y,Vx,Vy,Fx,Fy,m,R;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void CalculeFuerza(void);
  void Muevase(double dt);
  void Dibujese(void);
  double Getx(void){return x;};
  double Gety(void){return y;};
};

void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
  x=x0;
  y=y0;
  Vx=Vx0;
  Vy=Vy0;
  m=m0;
  R=R0;
}

void Cuerpo::CalculeFuerza(void){
  double aux=-GM*m*std::pow(x*x+y*y,-1.5);
  Fx=aux*x;
  Fy=aux*y;
}

void Cuerpo::Muevase(double dt){
  x+=Vx*dt;
  y+=Vy*dt;
  Vx+=Fx*dt/m;
  Vy+=Fy*dt/m;
}

void Cuerpo::Dibujese(void){
  std::cout<<", "<<x<<"+"<<R<<"*cos(t),"<<y<<"+"<<R<<"*sin(t)";
}
//------------------Funciones Globales---------

void InicieAnimacion(void){
  //std::cout<<"set terminal gif animate"<<std::endl;
  //std::cout<<"set output 'planeta.gif'"<<std::endl;
  std::cout<<"unset key"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set xrange [-120:120]"<<std::endl;
  std::cout<<"set yrange [-120:120]"<<std::endl;
  std::cout<<"set parametric"<<std::endl;
  std::cout<<"set trange[0:7]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;
}
void InicieCuadro(void){
  std::cout<<"plot 0,0 ";
}
void TermineCuadro(void){
  std::cout<<std::endl;
}

int main(void){
  double t, dt=0.01,tmax;
  double tdibujo;int Ndibujos;
  Cuerpo Planeta;

  double r,omega,V,T,m,R;
  R=5; m=1; r=200; omega=std::sqrt(GM/(r*r*r)); V=omega*r; T=2*M_PI/omega; tmax=1.1*T;  

  
  //-----(x0, y0, Vx0, Vy0, m0, R0);
  InicieAnimacion(); Ndibujos=50;
  Planeta.Inicie(r, 0, 0, 0.5*V, m, R);

  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    if(tdibujo<tmax/Ndibujos){
      InicieCuadro();
      Planeta.Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
    Planeta.CalculeFuerza();
    Planeta.Muevase(dt);
  }
  
  return 0;
}


