#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>

double f(const std::vector<double> & x,std::vector<double> & arg, double t, const int index);
void runge(std::vector<double> & x, const std::vector<double> & arg, double t, double h);
double f_lambda(double lambda);
double bisection(double a, double b, int N_MAX, double eps);




int main(void)
{
  /*
    std::cout<<bisection(2,3,100,1.0e-4)<<std::endl;
    std::cout<<bisection(5,6,100,1.0e-4)<<std::endl;
    std::cout<<bisection(8,9,100,1.0e-3)<<std::endl;
    std::cout<<bisection(11,12,100,1.0e-3)<<std::endl;
    std::cout<<bisection(14,15,100,1.0e-3)<<std::endl;
    std::cout<<bisection(17,19,100,1.0e-3)<<std::endl;
    std::cout<<bisection(20,22,100,1.0e-3)<<std::endl;
    std::cout<<bisection(23,25,100,1.0e-3)<<std::endl;
    std::cout<<bisection(27,28,100,1.0e-3)<<std::endl;
    std::cout<<bisection(30,31,100,1.0e-3)<<std::endl;
    std::cout<<bisection(33,34,100,1.0e-3)<<std::endl;
    std::cout<<bisection(36,37,100,1.0e-3)<<std::endl;
    std::cout<<bisection(39,41,100,1.0e-3)<<std::endl;
    std::cout<<bisection(42,43,100,1.0e-3)<<std::endl;
    std::cout<<bisection(45,46,100,1.0e-3)<<std::endl;
    std::cout<<bisection(49,50,100,1.0e-3)<<std::endl;
    std::cout<<bisection(52,53,100,1.0e-3)<<std::endl;
    std::cout<<bisection(55,56,100,1.0e-3)<<std::endl;
    std::cout<<bisection(58,59,100,1.0e-3)<<std::endl;
    std::cout<<bisection(61,62,100,1.0e-3)<<std::endl;
    std::cout<<bisection(64,66,100,1.0e-3)<<std::endl;
    std::cout<<bisection(67,68,100,1.0e-3)<<std::endl;
    std::cout<<bisection(70,72,100,1.0e-3)<<std::endl;
    std::cout<<bisection(74,77,100,1.0e-3)<<std::endl;
    std::cout<<bisection(77,78,100,1.0e-3)<<std::endl;    
  */

  double h = 0.1;
  double lim = 80;
  double N = lim/h;
  double lambda = 0;
  std::ofstream datos;

  datos.open("datos_b.dat");
  
  for(int ii = 0; ii < N; ii++ ){
    lambda=0.0+ii*h;
    datos<<lambda<<"  " << f_lambda(lambda) <<std::endl;
  }

  datos.close();
  return 0;
}




//------------ funciones que componen el sistema de ecuaciones

double f(const std::vector<double> & x,const std::vector<double> & arg, double t, const int index)
{
  double LAMBDA = arg[0];
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

double bisection(double a, double b, int N_MAX, double eps){
  double r;
  for(int ii = 1; ii<=N_MAX; ii++){
    r=0.5*(a+b);
    if(std::fabs(a-b)<eps){
      break;
    }
    else if(f_lambda(a)*f_lambda(r)>0){
      a=r;
    }
    else{
      b=r;
    }
  }
  return r;
}

double f_lambda(double lambda){
//--------------------{R(0) ,R'(0)}
  std::vector<double> x {1.0  , 0.0};
  //paso de r
  const double h = 0.01;
  double r;
  //limite superior del radio
  double Lim= 1;
  
  //numero de pasos del radio
  int N = Lim/h;

  //Argumentos del sistema-----{LAMBDA}
  std::vector<double>argumentos{lambda};
  
  
  for(int step = 0; step < N; ++step){
    r = 0.01 + step*h;
    runge(x,argumentos, r, h);
    
  }
  double value= x[0];
  return value;  
}
