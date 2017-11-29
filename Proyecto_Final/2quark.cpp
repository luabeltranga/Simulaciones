#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=1;
const int Ly=1;
const int Lz=8000;

const int Q=19;
const double W0=1.0/3;
const double DELTAT=0.8;

const double C=1.0; 

const double ETAsobreS = 0.0001;




class LatticeBoltzmann{
private:
  double w[Q];
  int V[3][Q];  //V[alpha][i]  alpha=0 es x, alpha=1 es y, alpha=2 es z
  double f[Lx][Ly][Lz][Q],fnew[Lx][Ly][Lz][Q]; // f[ix][iy][iz][i]
  double g[Lx][Ly][Lz][Q],gnew[Lx][Ly][Lz][Q]; // f[ix][iy][iz][i]
public:
  LatticeBoltzmann(void);
  double N(int ix,int iy,int iz,bool UseNew,double Ux0,double Uy0,double Uz0);
  double Ux(int ix,int iy,int iz,double P0);
  double Uy(int ix,int iy,int iz,double P0);
  double Uz(int ix,int iy,int iz,double P0);
  double Presion(int ix,int iy,int iz);
  double feq(int i,double N0,double Ux0,double Uy0,double Uz0);
  double geq(int i,double N0,double Ux0,double Uy0,double Uz0,double P0);
  void Inicie(double Ux0,double Uy0,double Uz0);
  void Colisione(void);
  void Adveccione(void);
  void Imprimase(char const * NombreArchivo,int t);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //Pesos
  w[0]=W0;
  w[1]=w[3]=w[5]=w[7]=w[13]=w[18]=1.0/18;
  w[2]=w[4]=w[6]=w[8]=w[9]=w[10]=w[11]=w[12]=w[14]=w[15]=w[16]=w[17]=1.0/36;
  
  //Vectores velocidad
  V[0][0]=0;  
  V[1][0]=0;
  V[2][0]=0;
  
  //centro del cubo
  V[0][1]=C;    V[0][2]=C;    V[0][3]=0;    V[0][4]=-C;   V[0][5]=-C;   V[0][6]=-C;   V[0][7]=0;   V[0][8]=C;  
  V[1][1]=0;    V[1][2]=C;    V[1][3]=C;    V[1][4]=C;    V[1][5]=0;    V[1][6]=-C;   V[1][7]=-C;  V[1][8]=-C;
  V[2][1]=0;    V[2][2]=0;    V[2][3]=0;    V[2][4]=0;    V[2][5]=0;    V[2][6]=0;    V[2][7]=0;   V[2][8]=0;

  //cara superior
  V[0][9]=C;    V[0][10]=0;    V[0][11]=-C;   V[0][12]=0;    V[0][13]=0;   
  V[1][9]=0;    V[1][10]=C;    V[1][11]=0;    V[1][12]=-C;   V[1][13]=0;    
  V[2][9]=C;    V[2][10]=C;    V[2][11]=C;    V[2][12]=C;    V[2][13]=C;

  //cara inferior 
  V[0][14]=C;    V[0][15]=0;    V[0][16]=-C;   V[0][17]=0;    V[0][18]=0;   
  V[1][14]=0;    V[1][15]=C;    V[1][16]=0;    V[1][17]=-C;   V[1][18]=0;    
  V[2][14]=-C;   V[2][15]=-C;   V[2][16]=-C;   V[2][17]=-C;   V[2][18]=-C;    
}
double LatticeBoltzmann::N(int ix,int iy,int iz,bool UseNew,double Ux0,double Uy0,double Uz0){
  double suma=0;
  double U2=Ux0*Ux0+Uy0*Uy0+Uz0*Uz0;
  double gamma=pow(1-U2,-0.5);
  
  for(int ii=0;ii<Q;ii++){
    if(UseNew)suma+=fnew[ix][iy][iz][ii];
    else      suma+=f[ix][iy][iz][ii];
  }
  return suma/gamma;
}
double LatticeBoltzmann::Ux(int ix,int iy,int iz,double P0){
  double suma=0,alpha=0,suma1=0;
  
  for(int ii=0;ii<Q;ii++){
    suma+=g[ix][iy][iz][ii];
    suma1+=g[ix][iy][iz][ii]*V[0][ii];
  }
  alpha=1/(3*suma+3*P0);
  
  return alpha*suma1;

}
double LatticeBoltzmann::Uy(int ix,int iy,int iz,double P0){
  double suma=0,alpha=0,suma1=0;
  
  for(int ii=0;ii<Q;ii++){
    suma+=g[ix][iy][iz][ii];
    suma1+=g[ix][iy][iz][ii]*V[1][ii];
  }
  alpha=1/(3*suma+3*P0);
  
  return alpha*suma1;
  
}
double LatticeBoltzmann::Uz(int ix,int iy,int iz,double P0){
  double suma=0,alpha=0,suma1=0;
  
  for(int ii=0;ii<Q;ii++){
    suma+=g[ix][iy][iz][ii];
    suma1+=g[ix][iy][iz][ii]*V[2][ii];
  }
  alpha=1/(3*suma+3*P0);
  
  return alpha*suma1;
}

double LatticeBoltzmann::Presion(int ix,int iy,int iz){
  double suma1=0,suma2=0;
  
  for(int ii=0;ii<Q;ii++){
    suma1+=pow(g[ix][iy][iz][ii]*(V[0][ii]*V[0][ii]+V[1][ii]*V[1][ii]+V[2][ii]*V[2][ii]),2);
    suma2+=g[ix][iy][iz][ii];
  }
  
  return -suma2/3.0+sqrt(-3*suma1*suma1+4*suma2*suma2)/3;
}
double LatticeBoltzmann::feq(int i,double N0,double Ux0,double Uy0,double Uz0){
  double U2=Ux0*Ux0+Uy0*Uy0+Uz0*Uz0;
  double gamma=pow(1-U2,-0.5);
  double CdotU=Ux0*V[0][i]+Uy0*V[1][i]+Uz0*V[2][i];
  return w[i]*N0*gamma*(1.0+3*(CdotU)/(C*C)+9/2*(CdotU*CdotU)/(C*C*C*C)-3/2*U2/(C*C));
}
double LatticeBoltzmann::geq(int i,double N0,double Ux0,double Uy0,double Uz0,double P0){
  double U2=Ux0*Ux0+Uy0*Uy0+Uz0*Uz0;
  double gamma2=1/(1-U2);
  double CdotU=Ux0*V[0][i]+Uy0*V[1][i]+Uz0*V[2][i];
  if(i==0){
    return 3*w[0]*P0*gamma2*(4-(2+C*C)/(gamma2*C*C)-2*(U2)/(C*C));
  }
  else{
    return 3*w[i]*P0*gamma2*(1/(gamma2*C*C)+4*(CdotU)/(C*C)+6*(CdotU*CdotU)/(C*C*C*C)-2*U2/(C*C));
  }
}
void LatticeBoltzmann::Inicie(double Ux0,double Uy0,double Uz0){
  double P0=0,N0=0;
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      for(int iz=0;iz<Lz;iz++){
	for(int ii=0;ii<Q;ii++){
	  P0=2.495e-7; N0=P0/0.0314;
	  if(iz<100){
	    
	    f[ix][iy][iz][ii]=feq(ii,N0,Ux0,Uy0,Uz0);
	    g[ix][iy][iz][ii]=geq(ii,N0,Ux0,Uy0,Uz0,P0);
	  }
	  else{
	    P0=1.023e-7; N0=P0/0.0314;
	    f[ix][iy][iz][ii]=feq(ii,N0,Ux0,Uy0,Uz0);
	    g[ix][iy][iz][ii]=geq(ii,N0,Ux0,Uy0,Uz0,P0);
	  }
	}
      }
    }
  }
}
void LatticeBoltzmann::Colisione(void){ //de f a fnew
  double N0,Ux0,Uy0,Uz0,P0; double sigma;
  double tau=0;
  double Utau=0;
  double UmUtau=0;
  double Lambda=0;
  double T=0;
 
  
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){ //Para cada celda
      for(int iz=0;iz<Lz;iz++){
	P0=Presion(ix,iy,iz);
       	Ux0=Ux(ix,iy,iz,P0);  Uy0=Uy(ix,iy,iz,P0); Uz0=Uz(ix,iy,iz,P0); 
	N0=N(ix,iy,iz,false,Ux0,Uy0,Uz0);  	
	T=N0/P0;
	tau=3*ETAsobreS*(4*N0-N0*log(Lambda))/(4*P0)+DELTAT/2.0;
	Utau=DELTAT/tau;
	UmUtau=1-Utau;
	Lambda=N0*M_PI*M_PI/(16*T*T*T);
	
  
	for(int ii=0;ii<Q;ii++){ //para cada dirección
	  fnew[ix][iy][iz][ii]=UmUtau*f[ix][iy][iz][ii]+Utau*feq(ii,N0,Ux0,Uy0,Uz0);
	  gnew[ix][iy][iz][ii]=UmUtau*g[ix][iy][iz][ii]+Utau*geq(ii,N0,Ux0,Uy0,Uz0,P0); 
	}
      }
    }
  }
}
void LatticeBoltzmann::Adveccione(void){ //de fnew a f
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      for(int iz=0;iz<Lz-1;iz++){
	for(int ii=0;ii<Q;ii++){
	  f[(ix+V[0][ii]+Lx)%Lx][(iy+V[1][ii]+Ly)%Ly][(iz+V[2][ii]+Lz)%Lz][ii]=fnew[ix][iy][iz][ii];
	  g[(ix+V[0][ii]+Lx)%Lx][(iy+V[1][ii]+Ly)%Ly][(iz+V[2][ii]+Lz)%Lz][ii]=gnew[ix][iy][iz][ii];
	}
      }
    }
  }
}
void LatticeBoltzmann::Imprimase(char const * NombreArchivo,int t){
  ofstream MiArchivo(NombreArchivo); double N0,Ux0,Uy0,Uz0,P0,U; double sigma;
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      for(int iz=0;iz<Lz/10;iz++){
	P0=Presion(ix,iy,iz);
       	Ux0=Ux(ix,iy,iz,P0);  Uy0=Uy(ix,iy,iz,P0); Uz0=Uz(ix,iy,iz,P0);
	N0=N(ix,iy,iz,false,Ux0,Uy0,Uz0);
	U=sqrt(Ux0*Ux0+Uy0*Uy0+Uz0*Uz0);
	//MiArchivo<<iz<<" "<<N0<<endl;
	MiArchivo<<iz<<" "<<P0/2.495e-7<<endl;
	//MiArchivo<<iz<<" "<<U<<endl;
      }
      MiArchivo<<endl;
    }
    MiArchivo.close();
  }
}
//---------------- Funciones Globales --------

int main(int argc,char ** argv){
  LatticeBoltzmann Ondas;
  int t,tmax=atoi(argv[1]);

  double Ux0=0,Uy0=0,Uz0=0.0;
  
  //Inicie
  Ondas.Inicie(Ux0,Uy0,Uz0);
  
  //Corra
  for(t=0;t<tmax;t++){
    Ondas.Colisione();
    Ondas.Adveccione();
  }
  
  //Mostrar Resultado.
  Ondas.Imprimase("dquark.dat",t);
  
  return 0;
}
