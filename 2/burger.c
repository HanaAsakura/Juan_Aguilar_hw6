#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
 
 int nx = 41;
 int ny = 41;
 int nt = 500;
 int c = 1;
 double dx = 0.05;
 double dy = 0.05;
 double sigma = 0.0009; 
 double nu = 0.01;
 double dt = sigma*dx*dy/nu;
 double *x;
 double *y;
 int count;
 int i;
 int j;
 

 
 x = malloc(nx * sizeof(double)); //de 0 a 2
 y = malloc(ny * sizeof(double)); //de 0 a 2
 
 for (i=0;i<nx;i++){
   x[i] = (double)i/20.0;

 }
 for (i=0;i<ny;i++){
   y[i] = (double)i/20.0;

 }
 
 double **u = malloc(sizeof *u * nx);
 if (u)
   {
     for (count = 0; count < nx; count++)
       {
	 u[count] = malloc(sizeof *u[count] * ny);
	 for (j = 0; j<ny; j++){
	   u[count][j] = 1;
	   //printf ("%f\n", u[count][j]);
	 }

       }
   }

 double **v = malloc(sizeof *v * nx);
 if (v)
   {
     for (count = 0; count < nx; count++)
       {
	 v[count] = malloc(sizeof *v[count] * ny);
	 for (j = 0; j<ny; j++){
	   v[count][j] = 1;
	 }

       }
   }
 
 double **un = malloc(sizeof *un * nx);
 if (un)
   {
     for (count = 0; count < nx; count++)
       {
	 un[count] = malloc(sizeof *un[count] * ny);
	 for (j = 0; j<ny; j++){
	   un[count][j] = 1;
	 }

       }
   }

 double **vn = malloc(sizeof *vn * nx);
 if (vn)
   {
     for (count = 0; count < nx; count++)
       {
	 vn[count] = malloc(sizeof *vn[count] * ny);
	 for (j = 0; j<ny; j++){
	   vn[count][j] = 1;
		
	 }

       }
   }

 double **comb = malloc(sizeof *comb * nx);
 if (comb)
   {
     for (count = 0; count < nx; count++)
       {
	 comb[count] = malloc(sizeof *comb[count] * ny);
	 for (j = 0; j<ny; j++){
	   comb[count][j] = 1;
	 }

       }
   }
 
 //v = numpy.ones((ny,nx))
 //un = numpy.ones((ny,nx)) ##
 //vn = numpy.ones((ny,nx))
 //comb = numpy.ones((ny,nx))

 //Assign initial conditions
 for (i = 10; i<21; i++){
   for (j = 10; j<21; j++){
     u[i][j] = 2;
     v[i][j] = 2;
   }

 }
 int n;
 //FILE *fp;
 //fp = fopen("burgers_data.dat", "ab+");
 
 for (n=0; n<nt+1;n++){
   
   for (i=0; i<nx;i++){
     for (j=0; j<ny; j++){
       un[i][j] = u[i][j];
       vn[i][j] = v[i][j];
     }
   }

   for (i=0; i<nx-2; i++){
     for(j=0; j<ny-2; j++){
       
       u[i+1][j+1] = un[i+1][j+1]-(dt/dx)*un[i+1][j+1]*(un[i+1][j+1]-un[i+1][j])-(dt/dy)*vn[i+1][j+1]*(un[i+1][j+1]-un[i][j+1])+nu*((double)dt/(double)pow(dx,2))*(un[i+1][j+2]-2*un[i+1][j+1]+un[i+1][j])+((double)(nu*dt)/(double)pow(dy,2)) * (un[i+2][j+1]-2*un[i+1][j+1]+un[i][j+1]); 
	//printf ("%f\n", u[i+1][j+1]); 
    
       v[i+1][j+1] = vn[i+1][j+1]-(dt/dx)*un[i+1][j+1]*(vn[i+1][j+1]-vn[i+1][j])-(dt/dy)*vn[i+1][j+1]*(vn[i+1][j+1]-vn[i][j+1])+nu*((double)dt/(double)pow(dx,2))*(vn[i+1][j+2]-2*vn[i+1][j+1]+vn[i+1][j])+((double)(nu*dt)/(double)pow(dy,2))*(vn[i+2][j+1]-2*vn[i+1][j+1]+vn[i][j+1]); 
	//printf ("%f\n", v[i+1][j+1]);
     } 
   }
   for (i =0; i<nx; i++){
     u[0][i] = 1;
     u[nx-1][i] = 1;
     u[i][0] = 1;
     u[i][nx-1] = 1;
     v[0][i] = 1;
     v[ny-1][i] = 1; 
     v[i][0] = 1;
     v[i][ny-1] = 1;
   }
   for (i=0; i<nx; i++){
	
     for (j=0; j<ny; j++){
       printf("%f %f\n",u[i][j], v[i][j]);	
     } 
   }
   
   
   
   

   
 }
 //fclose(fp);
 



  return 0;
}
