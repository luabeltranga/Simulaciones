#include <iostream>
#include "Random64.h"
#include <cmath>

int main(void){
  double test =0;
  Crandom ran64(2);

  

  for(int ii = 0; ii<10000;ii++){
    test=std::fmod(ran64.r(),0.3);    
    std::cout<<test<<std::endl;
  }
  
  return 0;
}
