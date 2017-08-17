#include<iostream>
#include<vector>
#include<fstream>

const double LAMBDA=2.4;

double f(const std::vector<double> & x, double t, const int index);
void runge(std::vector<double> & x, double t, double h);


int main(void)
{
  //--------------------{R(0) ,R'(0)}
  std::vector<double> x {1.0,0.0};

  //paso de tiempo
  const double h = 0.01;

  //limite superior del radio
  double Lim= 10;

  //numero de pasos
  int N = Lim/h;
  double r;
  std::ofstream curvas;
  curvas.open("Punto_a.dat");
  
  for(int step = 1; step < N; ++step){
    r = step*h;
    curvas<< r << "   " <<x[0] << "   " << x[1]  <<  std::endl;
    runge(x, r, h);
  }

  curvas.close();
  return 0;
}
//------------ funciones que componen el sistema de ecuaciones
double f(const std::vector<double> & x, double t, const int index)
{
  if (0 == index)
    return x[1];
  else if (1 == index)
    return -(x[1]/t+LAMBDA*LAMBDA*x[0]);
  else {
    std::cerr << "Wrong index = " << index << std::endl;
    exit(1);
  }
}

//-----------metodo de runge kutta 4
void runge(std::vector<double> & x, double t, double h)
{
  const int N = x.size();
  std::vector<double> xtmp (N) ;
  std::vector<double>  k1 (N);
  std::vector<double>  k2 (N);
  std::vector<double>  k3 (N);
  std::vector<double>  k4 (N);

  //k1
  for (int ii = 0; ii < N; ++ii) {
    k1[ii] = h*f(x, t, ii);
  }

  //k2
  for (int ii = 0; ii < N; ++ii) {
    xtmp[ii] = x[ii]+k1[ii]/2;
  }
  for (int ii = 0; ii < N; ++ii) {
    k2[ii] = h*f(xtmp, t+h/2, ii);
  }

  //k3
  for (int ii = 0; ii < N; ++ii) {
    xtmp[ii] = x[ii]+k2[ii]/2;
  }
  for (int ii = 0; ii < N; ++ii) {
    k3[ii] = h*f(xtmp, t+h/2, ii);
  }

  //k4
  for (int ii = 0; ii < N; ++ii) {
    xtmp[ii] = x[ii]+k3[ii];
  }
  for (int ii = 0; ii < N; ++ii) {
    k4[ii] = h*f(xtmp, t+h, ii);
  }

  //end
  for (int ii = 0; ii < N; ++ii) {
    x[ii] = x[ii] +( k1[ii] + 2*k2[ii] +2*k3[ii] +  k4[ii])/6.0;
  }
}

