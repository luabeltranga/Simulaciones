#include<iostream>
#include<vector>
#include<fstream>

const double BETA=0.3;
const double GAMMA=0.05;

double f(const std::vector<double> & x, double t, const int index);
void runge(std::vector<double> & x, double t, double h);


int main(void)
{
  //--------------------{s(0) ,i(0)  ,r(0)}
  std::vector<double> x {0.999,0.001 , 0.0};

  //paso de tiempo
  const double h = 0.01;

  //limite superior del tiempo
  double Lim= 200;

  //numero de pasos
  int N = Lim/h;
  double t;
  std::ofstream curvas;
  curvas.open("SIR.txt");
  
  for(int step = 0; step < N; ++step){
    t = 0.0 + step*h;
    curvas<< t << "   " <<x[0] << "   " << x[1] << "  " << x[2] <<  std::endl;
    runge(x, t, h);
  }

  curvas.close();
  return 0;
}
//------------ funciones que componen el sistema de ecuaciones
double f(const std::vector<double> & x, double t, const int index)
{
  if (0 == index)
    return -BETA*x[0]*x[1];
  else if (1 == index)
    return BETA*x[0]*x[1]-GAMMA*x[1];
  else if (2 == index)
    return GAMMA*x[1];
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

