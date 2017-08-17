#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>

const double  ERR = 1.0e-7;

double f(const std::vector<double> & x,std::vector<double> & arg, double t, const int index);
void runge(std::vector<double> & x, const std::vector<double> & arg, double t, double h);

int main(void)
{
  //--------------------{s(0) ,i(0)  ,r(0)}
  std::vector<double> x {0.999,0.001 , 0.0};

  //paso de tiempo
  const double h = 0.001;
  double t;

  //Argumentos del sistema-----{BETA,GAMMA}
  std::vector<double>argumentos{0.1,0.1};
  
  //parametro R0
  double R0 = 1;

  //variable de comparacion para S_infinito
  double tmp = 0;
  

  
  std::ofstream curvas;
  curvas.open("SIR.txt");
  for(int ii = 0;ii<10000;ii++){
    // BETA=R0*GAMMA
    argumentos[0]=R0*argumentos[1];
    //----evoluciona con runge kutta el sistema hasta S_infinito
    do{
      tmp = x[0];
      t += h;
      runge(x,argumentos, t, h);
    }while(std::fabs(tmp-x[0])>ERR);
    curvas<< R0 << "   " <<x[0] << std::endl;
    R0+=0.001;
  }
  curvas.close();
  return 0;
}
//------------ funciones que componen el sistema de ecuaciones

double f(const std::vector<double> & x,const std::vector<double> & arg, double t, const int index)
{
  double BETA = arg[0];
  double GAMMA = arg[1];
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
void runge(std::vector<double> & x, const std::vector<double> & arg, double t, double h)
{
  const int N = x.size();
  std::vector<double> xtmp (N) ;
  std::vector<double>  k1 (N);
  std::vector<double>  k2 (N);
  std::vector<double>  k3 (N);
  std::vector<double>  k4 (N);

  //k1
  for (int ii = 0; ii < N; ++ii) {
    k1[ii] = h*f(x,arg, t, ii);
  }

  //k2
  for (int ii = 0; ii < N; ++ii) {
    xtmp[ii] = x[ii]+k1[ii]/2;
  }
  for (int ii = 0; ii < N; ++ii) {
    k2[ii] = h*f(xtmp,arg, t+h/2, ii);
  }

  //k3
  for (int ii = 0; ii < N; ++ii) {
    xtmp[ii] = x[ii]+k2[ii]/2;
  }
  for (int ii = 0; ii < N; ++ii) {
    k3[ii] = h*f(xtmp,arg, t+h/2, ii);
  }

  //k4
  for (int ii = 0; ii < N; ++ii) {
    xtmp[ii] = x[ii]+k3[ii];
  }
  for (int ii = 0; ii < N; ++ii) {
    k4[ii] = h*f(xtmp,arg, t+h, ii);
  }

  //end
  for (int ii = 0; ii < N; ++ii) {
    x[ii] = x[ii] +( k1[ii] + 2*k2[ii] +2*k3[ii] +  k4[ii])/6.0;
  }
}

