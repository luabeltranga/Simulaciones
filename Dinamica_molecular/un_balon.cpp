#include <iostream>
#include <cmath>

const double g=9.795;

class Cuerpo;

int main(void){
  double t; double dt=0.01;
  Cuerpo Balon;

  //-----(x0, y0, Vx0, Vy0, m0, R0);
  Balon.Inicie(0, 0, 12, 16, 0.457, 0.15);

  for(t=0;t<3.5;t+=dt){
    std::cout<<Balon.Getx()<<" "<<Balon.Gety()<<std::endl;
    Balon.CalculeFuerza();
    Balon.Muevase(dt);
  }
  
  return 0;
}

class Cuerpo{
private:
  double x,y,Vx,Vy,Fx,Fy,m,R;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void Muevase(double dt);
  void CalculeFuerza(void);
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
  Fx=0;
  Fy=-m*g;
}

void Cuerpo::Muevase(double dt){
  x+=Vx*dt;
  y+=Vy*dt;
  Vx+=Fx*dt/m;
  Vy+=Fy*dt/m;
}




