#include <iostream>
#include <cmath>
#include "Random64.h"

int main(void){
  Crandom ran64(1);
  double tau = 2;
  double mu = 3, sigma=2;
  for(int ii = 0; ii<1.0e4;ii++){
    //std::cout<<-tau*std::log(ran64.r())<<std::endl;
    std::cout<< sigma*std::sqrt(-2*std::log(ran64.r()))*std::cos(2*M_PI*ran64.r())+mu<<std::endl;
  }
  return 0;
}
