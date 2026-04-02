#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
#include "cblas.h"
#include "time.h"

int N=20;
double *list1;
//double *list2;

int main(void)
{

/*Reading the tables*/
list1 = (double *)malloc(sizeof(double)*(N*(N+1)/2));
//list2 = (double *)malloc(sizeof(double)*(N*(N+1)/2));
FILE *FileIn;
FileIn = fopen("/home/ksardase/Work/class_2.6.3_PT/pt_matrices/M22oneline_N256.dat", "r");
int count;      
for (count=0; count < (N)*(N+1)/2; count++){
   fscanf(FileIn, "%lf", &list1[count]);
   }
/*for (count=0; count < (N)*(N+1)/2; count++){
   fscanf(FileIn, "%lf", &list2[count]);
   }*/
fclose(FileIn);

//test
/*for (count=0; count<(N)*(N+1)/2; ++count)  printf("%e\n", list1[count]);
printf("\n");*/

/*Construct vector*/
double *x;
double *y;
x = (double *)malloc(sizeof(double)*N);
y = (double *)malloc(sizeof(double)*N);
for (count=0; count < N; count++){
//for (count=(N)*(N+1)/2; count < (N)*(N+1); count++){
   x[count]=2.;
   y[count]=0.;
   printf("x_i=%e\n",x[count]);
   }

int start=clock();
double c;
cblas_dsymv(CblasRowMajor,CblasLower, N, 1.0, list1, N, x, 1, 0.0, y, 1);
for (count=0; count < N; count++){
printf("y_i=%e\n",y[count]);
}
for (count=0; count < N*(N+1)/2; count++){
printf("list1_i=%e\n",list1[count]);
}
c=cblas_ddot(N, y, 1, y, 1);
printf("output=%e\n",c);

//double m[] = {2, 3, 7, 5, 5, 8};
//double x2[] = {9, -1, 1};
//double y2[] = { 0, 0, 0};
//cblas_dsymv(CblasRowMajor,CblasLower, 3, 1.0, m, 3, x2, 1, 0.0, y2, 1);

int end=clock();
printf("Rate_total=%d\n",end-start);
printf("Rate_per_second=%ld\n",CLOCKS_PER_SEC);

//free(list1);
}
