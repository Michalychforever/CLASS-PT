//gcc-8 -o matrix_lapack matrix_lapack.c /opt/OpenBLAS/lib/libopenblas.a -lpthread -lm;./matrix_lapack

#include "stdio.h"
#include "stdlib.h"
#include "complex.h"
#include "sys/time.h"
#include "time.h"

extern void dsymv_(char*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern double ddot_(int*, double*, int*, double*, int*);

extern void zspmv_(char*, int*, double complex*, double complex*, double complex*, int*, double complex*, double complex*, int*);
//extern struct doublecomplex zdotu_(int*, struct doublecomplex*, int*, struct doublecomplex*, int*);
extern double complex zdotu_(int*, double complex*, int*, double complex*, int*);

int main(void)
{




int N=4;
double complex *list1;
//double *list2;
char uplo = 'L';
int inc = 1;
double complex alpha = 1.;
double complex beta = 0.;
int count;

/*Reading the tables*/
list1 = (double complex*)malloc(sizeof(double complex) * (N*(N+1)/2));
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

/*for (count=0; count < N*(N+1)/2; count++){
   list1[count]=1.;
   printf("list1_real=%e list1_imag=%e\n",creal(list1[count]),cimag(list1[count]));
}*/

list1[0]=1.;
list1[1]=2.+I;
list1[2]=3. * I;
list1[3]=4.;
list1[4]=5.;
list1[5]=6.;
list1[6]=7.+2. * I;
list1[7]=8.;
list1[8]=9.;
list1[9]=0.;

for (count=0; count < N*(N+1)/2; count++){
printf("list1_real=%e list1_imag=%e\n",creal(list1[count]),cimag(list1[count]));
}

//test
/*for (count=0; count<(N)*(N+1)/2; ++count)  printf("%e\n", list1[count]);
printf("\n");*/

/*Construct vector*/
double complex *x;
double complex *y;
x = (double complex *)malloc(sizeof(double complex)*N);
y = (double complex *)malloc(sizeof(double complex)*N);
/*for (count=0; count < N; count++){
   x[count]=2.;
   y[count]=0.;
   //printf("x_i=%e\n",x[count]);
}*/

x[0]=I+1.;
x[1]=3.-I/2.;
x[2]=4.;
x[3]=5.;
printf("\n");
for (count=0; count < N; count++){
printf("x_real=%e x_imag=%e\n",creal(x[count]),cimag(x[count]));
}

int start=clock();
//double complex c;
//struct doublecomplex c;
zspmv_(&uplo, &N, &alpha, list1, x, &inc, &beta, y, &inc);

double complex c;
c=zdotu_(&N, x, &inc, y, &inc);

printf("c_real=%e c_image=%e\n",creal(c),cimag(c));

printf("\n");
for (count=0; count < N; count++){
printf("y_real=%e y_imag=%e\n",creal(y[count]),cimag(y[count]));
}

int end=clock();
printf("Rate_total=%d\n",end-start);
printf("Rate_per_second=%ld\n",CLOCKS_PER_SEC);

free(x);
free(y);
free(list1);
return 0;
}
