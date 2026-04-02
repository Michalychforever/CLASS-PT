//gcc-8 -o matrix_blas2 matrix_blas2.c /opt/OpenBLAS/lib/libopenblas.a -lpthread -lm;./matrix_blas2;./matrix_blas2

/*for real quantities*/
#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
#include "time.h"

extern void dspmv_(char*, int*, double*, double*, double*, int*, double*, double*, int*);
extern double ddot_(int*, double*, int*, double*, int*);

//extern void zspmv_(char*, int*, double complex*, double complex*, double complex*, int*, double complex*, double complex*, int*);
//extern double complex zdotu_(int*, double complex*, int*, double complex*, int*);

int main(void) 
{

int N=257;
double* list1;
//double *list2;
char uplo = 'L';
int inc = 1;
double alpha = 1.;
double beta = 0.;
int count;

/*Reading the tables*/
list1 = (double*)malloc(sizeof(double) * (N*(N+1)/2));
FILE *FileIn;
FileIn = fopen("/home/ksardase/Work/class_2.6.3_PT/pt_matrices/M22oneline_N256.dat", "r");    
for (count=0; count < N*(N+1)/2; count++){
   fscanf(FileIn, "%lf", &list1[count]);
   //printf("list1_i=%e\n",list1[count]);
}

/*Construct vector*/
double *x;
double *y;
x = (double *)malloc(sizeof(double)*N);
y = (double *)malloc(sizeof(double)*N);
for (count=0; count < N; count++){
   x[count]=2.;
   y[count]=0.;
   //printf("x_i=%e\n",x[count]);
}

int start=clock();
double c;
dspmv_(&uplo, &N, &alpha, list1, x, &inc, &beta, y, &inc);
/*for (count=0; count < N; count++){
printf("y_i=%lf\n",y[count]);
}*/
c=ddot_(&N, x, &inc, y, &inc);
printf("output=%e\n",c);

int end=clock();
printf("Rate_total=%d\n",end-start);
printf("Rate_per_second=%ld\n",CLOCKS_PER_SEC);

free(x);
free(y);
free(list1);
return 0;
}
