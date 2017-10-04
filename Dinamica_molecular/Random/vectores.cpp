#include <iostream>
#include <cmath>
#include <fstream>
#include "Vector.h"


int main(void){
  vector3D a,b;
  a.cargue(1,2,3);
  b=a;
  b+=a;
  std::cout<<a.x()<<std::endl;
  a.show();
  b.show();
  return 0;
}


