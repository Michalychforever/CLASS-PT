//gcc-8 -o matrix_blas matrix_blas.c /opt/OpenBLAS/lib/libopenblas.a -lpthread -lm;./matrix_blas;./matrix_blas

#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
#include "time.h"

extern void dsymv_(char*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern double ddot_(int*, double*, int*, double*, int*);

int main(void) 
{

int N=256;
double* list1;
//double *list2;
char uplo = 'L';
int inc = 1;
int inc2 = 0;
double alpha = 1.;
double beta = 0.;
int count;

/*Reading the tables*/
list1 = (double*)malloc(sizeof(double) * (N*N));
//list2 = (double *)malloc(sizeof(double)*(N*(N+1)/2));
/*FILE *FileIn;
FileIn = fopen("/home/ksardase/Work/class_2.6.3_PT/pt_matrices/M22oneline_N256.dat", "r");    
for (count=0; count < (N)*(N+1)/2; count++){
   fscanf(FileIn, "%lf", &list1[count]);
   }*/
/*for (count=0; count < (N)*(N+1)/2; count++){
   fscanf(FileIn, "%lf", &list2[count]);
   }*/
//fclose(FileIn);

for (count=0; count < N*N; count++){
   list1[count]=1.;
   printf("list1_i=%e\n",list1[count]);
}

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
   //printf("x_i=%e\n",x[count]);
}

int start=clock();
double c;
dsymv_(&uplo, &N, &alpha, list1, &N, x, &inc, &beta, y, &inc);
/*for (count=0; count < N*(N+1)/2; count++){
printf("list1_i=%e\n",list1[count]);
}*/
/*for (count=0; count < N; count++){
printf("y_i=%lf\n",y[count]);
}*/
c=ddot_(&N, y, &inc, y, &inc);
printf("output=%e\n",c);

int end=clock();
printf("Rate_total=%d\n",end-start);
printf("Rate_per_second=%ld\n",CLOCKS_PER_SEC);

free(x);
free(y);
free(list1);
return 0;
}
