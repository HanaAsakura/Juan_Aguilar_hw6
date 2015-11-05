#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define USAGE "./a.out valor1 valor2"
#include <string.h>


void velocidad (double nE0,double nalpha, double *x, double *y, double *z, double *t,double *vx, double *vy, double *vz, int npoints);

int main(int argc, char **argv){
 
  
  FILE *fp;

  double nE0 = atof(argv[1]);
  double nalpha = atof(argv[2]);
  
  char output1[512];
  char output2[512];
  char* line = "_";
  char* trayectoria = "trayectoria";
  char* dat = ".dat";
  char name[512] = "";
  strcat(name, trayectoria);
  strcat(name,line);
  strcat(name,argv[1]);
  strcat(name,line);
  strcat(name,argv[2]);
  strcat(name,dat);
  //printf("%s\n", name);
		    
		   

  

  fp = fopen(name, "ab+");
  
  
 
  
  //if(argc!=3){
  // printf("USAGE: %s\n", USAGE);
  // exit(1);
  //}
 
  
  double *x;
  double *y;
  double *z;
  double *vx;
  double *vy;
  double *vz;
  double *t;
  double tmax = 100.0;
  double dt = 0.0001;
  int npoints = (int)(tmax/dt);

  x  = malloc(npoints * sizeof(float));
  y  = malloc(npoints * sizeof(float));
  z  = malloc(npoints * sizeof(float));
  vx = malloc(npoints * sizeof(float));
  vy = malloc(npoints * sizeof(float));
  vz = malloc(npoints * sizeof(float));
  t  = malloc(npoints * sizeof(float));

  velocidad(nE0, nalpha, x, y, z, t, vx, vy, vz, npoints);
  
  int counter;

  for (counter=0; counter<npoints; counter++){
    fprintf(fp, "%f %f %f %f\n",x[counter], y[counter], z[counter], t[counter]);
  }
  fclose(fp);
  


  return 0;
}

void velocidad (double nE0,double nalpha, double *x, double *y, double *z, double *t,double *vx, double *vy, double *vz, int npoints){
    
  double dt = 0.01;
  double pi = 3.14159;
  double B0 = 3*pow(10,5);
  double m = 1.672621777*pow(10,-27); 
  double c = 299792458; 
  double q = 1.602176565*pow(10,-19);

  double alpha = nalpha*(pi/180.0);
  double E0 = nE0*1.6021773*pow(10,-13);
  double v0 = c*pow((1-1/pow((E0/(m*pow(c,2)))+1,2)),0.5);
  double gamma =pow(1/(1-(pow(v0,2)/pow(c,2))),0.5);
  double v0x = 0.0; 
  double v0y = v0*sin(alpha); 
  double v0z = v0*cos(alpha);
  double rt = 6378100;
  double x0 = 2*rt; 
  double y0 = 0.0; 
  double z0 = 0.0;
  int i = 0;
  //printf("%f\n",gamma);
  

 
    
  x[0] = x0;
  y[0] = y0;
  z[0] = z0;
  t[0] = 0.0;
    
    
  vx[0] = v0x;
  vy[0] = v0y;
  vz[0] = v0z;
    
  
    
  double r = pow(pow(x[0],2) + pow(y[0],2) + pow(z[0],2),0.5);
    
  double Bx = (-B0*(pow(rt,3))/pow(r,5))*(3*x[0]*z[0]);
  double By = (-B0*(pow(rt,3))/pow(r,5))*(3*y[0]*z[0]);
  double Bz = (-B0*(pow(rt,3))/pow(r,5))*(2*pow(z[0],2)-pow(x[0],2)-pow(y[0],2));
  
  printf("%f %f %f\n",v0x, v0y, v0z); 
    
    
  for (i = 1; i< npoints; i++){
        
        
        
    vx[i] = (((q*dt)/(gamma*m))*(vy[i-1]*Bz-vz[i-1]*By)) + vx[i-1];
    vy[i] = -(((q*dt)/(gamma*m))*(vx[i-1]*Bz-vz[i-1]*Bx)) + vy[i-1];
    vz[i] = (((q*dt)/(gamma*m))*(vx[i-1]*By-vy[i-1]*Bx)) + vz[i-1];
    //printf("%f\n",((q*dt)/(gamma*m)));
    printf("%f %f %f\n",vx[i], vy[i], vz[i]);
    //printf("%f %f %f\n",(((q*dt)/(gamma*m))*(vy[1]*Bz-vz[1]*By)), vy[2], vz[2]);
    x[i] = x[i-1] + vx[i]*dt;
    y[i] = y[i-1] + vy[i]*dt;
    z[i] = z[i-1] + vz[i]*dt;
    t[i] = t[i-1] + dt;
        
    r = pow(pow(x[i],2) + pow(y[i],2) + pow(z[i],2),0.5);

       
    Bx = (-B0*(pow(rt,3))/pow(r,5))*(3*x[i]*z[i]);
    By = (-B0*(pow(rt,3))/pow(r,5))*(3*y[i]*z[i]);
    Bz = (-B0*(pow(rt,3))/pow(r,5))*(2*pow(z[i],2)-pow(x[i],2)-pow(y[i],2));
    printf("%f %f %f\n",Bx, By,Bz);
    
  }
  //printf("%f %f %f\n",(((q*dt)/(gamma*m))*(vy[1]*Bz-vz[1]*By)), vy[2], vz[2]);     
        
     
}
