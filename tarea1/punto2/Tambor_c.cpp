#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>

const double ERR=1.0e-7;

double f(const std::vector<double> & x,std::vector<double> & arg, double t, const int index);
void runge(std::vector<double> & x, const std::vector<double> & arg, double t, double h);
double f_lambda(double lambda);
double bisection(double a, double b);
void roots_lambda(std::vector<double> & roots, double & lambda);
void plot_functions(const std::vector<double> & roots);


int main(void)
{
  //-----------guarda los primeros 25 valores de lambda para los cuales f(r=5)=0
  std::vector<double> roots (5);
  double r = 5;
  //----------toma las raices de j_0(lambda) y las divide entre 5
  roots_lambda(roots,r);

  //dibuja las funciones de bessel para los valores de lambda en los cuales f(lambda*5)=0
  plot_functions(roots);
  
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

//metodo de biseccion

double bisection(double a,double b){
  double m;
  double fa,fm;
  //ENTRADA DE DATOS
  fa=f_lambda(a);
  //PROCESAMIENTO
  do{ //repita
  // calcule m; calcule fm;
    m=(a+b)/2; fm=f_lambda(m);
    if(fa*fm<0)// if(fa y fm tienen signos opuestos)
      {b=m;} //corra b;
  else
    {a=m;fa=fm;}//corra a;
  } while(b-a>ERR);//hasta lograr la precision deseada
  //SALIDA DE DATOS
  return (a+b)/2;
}

//funcion que calcula el valor de f(lambda*r) con r=1
double f_lambda(double lambda){
//--------------------{R(0) ,R'(0)}
  std::vector<double> x {1.0  , 0.0};
  std::vector<double> x_tmp = x;
  //paso de r
  const double h = 0.01;
  double r;
  //limite superior del radio
  double Lim= 1;
  
  //numero de pasos del radio
  int N = Lim/h;

  //limite superior de LAMBDA
  double Lima=lambda;
  
  //numero de pasos del radio
  int Na = Lima/h;
  
  //Argumentos del sistema-----{LAMBDA}
  std::vector<double>argumentos{0.01};
  
  for(int ii = 0; ii<Na;ii++){
    x_tmp = x;
    argumentos[0]=0.0+ii*h;
    for(int step = 0; step < N; ++step){
      r = 0.01 + step*h;
      runge(x_tmp,argumentos, r, h);
    }
  }
  double value= x_tmp[0];
  return value;  
}


//por medio de biseccion encuentra ceros alrededor de los numeros
//que se encontro en el punto b
void roots_lambda(std::vector<double> & roots,double & lambda){
  std::cout<<"comenze con roots"<<std::endl;
  roots[0]=bisection(2,3)/lambda;
  std::cout<<roots[0]<<std::endl;
  roots[1]=bisection(5,5.5)/lambda;
  std::cout<<roots[1]<<std::endl;
  roots[2]=bisection(8.5,9)/lambda;
  std::cout<<roots[2]<<std::endl;
  roots[3]=bisection(11.5,12)/lambda;
  std::cout<<roots[3]<<std::endl;
  roots[4]=bisection(14.5,15)/lambda;
  std::cout<<roots[4]<<std::endl;
  /*std::cout<<"comenze con roots"<<std::endl;
  roots[5]=bisection(17,19,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[6]=bisection(20,22,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[7]=bisection(23,25,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[8]=bisection(27,28,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[9]=bisection(30,31,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[10]=bisection(33,34,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[11]=bisection(36,37,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[12]=bisection(39,41,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[13]=bisection(42,43,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[14]=bisection(45,46,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[15]=bisection(49,50,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[16]=bisection(52,53,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[17]=bisection(55,56,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[18]=bisection(58,59,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[19]=bisection(61,62,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[20]=bisection(64,66,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[21]=bisection(67,68,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[22]=bisection(70,72,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[23]=bisection(74,77,100,1.0e-3)/lambda;
  std::cout<<"comenze con roots"<<std::endl;
  roots[24]=bisection(77,78,100,1.0e-3)/lambda;    
  std::cout<<"termine roots"<<std::endl;
  */
}
  
      
//dibuja para   
void plot_functions(const std::vector<double> & roots){
  std::cout<<"empeze graficas"<<std::endl;
  //paso de r
  const double h = 0.01;
  double r;
  //limite superior del radio
  double Lim= 5.0;
  
  //numero de pasos del radio
  int N = Lim/h;
  //Argumentos del sistema-----{LAMBDA}
  std::vector<double>argumentos(1);
//--------------------{R(0) ,R'(0)}
  std::vector<double> x {1.0  , 0.0};
  std::vector<double> x_tmp = x;
  std::vector<double> x_tmp2 (50);

  std::ofstream curvas;
  curvas.open("Tambor.dat");

  for(int jj=0;jj<5;jj++){
    x_tmp2[2*jj]=x[0];
    x_tmp2[2*jj+1]=x[1];
  }

  
  
  for(int step = 0; step < N; ++step){
    r = 0.01 + step*h;
    curvas<<r<<"  ";
    
    for(int ii = 0; ii<5;ii++){
      x_tmp[0]=x_tmp2[2*ii];
      x_tmp[1]=x_tmp2[2*ii+1];
    
      argumentos[0]=roots[ii];
      runge(x_tmp,argumentos, r, h);
      curvas<<x_tmp[0]<<"  ";
      x_tmp2[2*ii]=x_tmp[0];
      x_tmp2[2*ii+1]=x_tmp[1];
    }
    curvas<<std::endl;
  }
  curvas.close();
}
