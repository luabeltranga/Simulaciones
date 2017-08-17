#include <iostream>
#include <vector>
#include <cstdlib>

const double OMEGA = 0.2323;

double f(const std::vector<double> & x, double t, const int index);
void runge(std::vector<double> & x, double t, double h);

int main(void)
{
  const double x0 = 1.2345;
  const double h = 0.1;
  const int N = 1000;

  std::vector<double> x {x0, 0};
  double t;
  for(int step = 0; step < N; ++step){
    t = 0.0 + step*h;
    std::cout << t << "   " << x[0] << "  " << x[1] <<  std::endl;
    runge(x, t, h);
  }
  
  
  return 0;
}
double f(const std::vector<double> & x, double t, const int index)
{
  if (0 == index)
    return x[1];
  else if (1 == index)
    return -OMEGA*OMEGA*x[0];
  else {
    std::cerr << "Wrong index = " << index << std::endl;
    exit(1);
  }
}

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
