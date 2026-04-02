//gcc -o matrix_packed matrix_packed.c;./matrix_packed

#include <stdio.h>
#include <stdlib.h>

int N=257;
double *list1;
double *list2;

int main(void)
{

/*Reading the tables*/
list1 = (double *)malloc(sizeof(double)*(N*(N+1)/2));
list2 = (double *)malloc(sizeof(double)*(N*(N+1)/2));
FILE *FileIn;
FileIn = fopen("/home/ksardase/Work/class_2.6.3_PT/pt_matrices/M22basiconeline_N256.dat", "r");
int count;      
for (count=0; count < (N)*(N+1)/2; count++){
   fscanf(FileIn, "%lf", &list1[count]);
   //printf("%.15le\n",list1[count]);
   }
for (count=0; count < (N)*(N+1)/2; count++){
   fscanf(FileIn, "%lf", &list2[count]);
   }
fclose(FileIn);

//test
/*for (count=0; count<(N)*(N+1)/2; ++count)  printf("%e\n", list1[count]);
printf("\n");*/

/*Construct matrices*/
double m1[N][N];
double m2[N][N];
int i,j;
count=0;
for (i=0; i<N; ++i){
    for (j=0; j<=i; ++j){
        m1[i][j]=list1[count];
        m2[i][j]=list2[count];
        count++;
    }
}

//test
/*for (i=0; i<N; ++i){
    for (j=0; j<=i; ++j){
      printf("%e\n", m[i][j]);
    }
}
printf("\n");*/

/*Construct packed matrices*/
double list1p[N*(N+1)/2];
double list2p[N*(N+1)/2];
for (i=0; i<N; ++i){
    for (j=0; j<=i; ++j){
        list1p[i+(2*N-1-j)*j/2]=m1[i][j];
        list2p[i+(2*N-1-j)*j/2]=m2[i][j];
    }
}

//test
/*for (count=0; count<(N)*(N+1)/2; ++count)  printf("%e\n", list1p[count]);*/

/*Writing real output*/
FILE *FileOut1;
FileOut1 = fopen("/home/ksardase/Work/class_2.6.3_PT/pt_matrices/M22basiconeline_N256_packed.dat", "w");
for (count=0; count < (N)*(N+1)/2; count++){
   fprintf(FileOut1, "%.18le\n", list1p[count]);
}
for (count=0; count < (N)*(N+1)/2; count++){
   fprintf(FileOut1, "%.18le\n", list2p[count]);
}
fclose(FileOut1);

/*Sum of elements in matrices*/
double s1 = 0.;
double s2 = 0.;
for (i=0; i<N; ++i){
    for (j=0; j<=i; ++j){
        if (i != j) {
            s1+=2.*m1[i][j];
            s2+=2.*m2[i][j];
        }
        else{
            s1+=m1[i][j];
            s2+=m2[i][j];
        }
    }
}
printf("sumRe=%le sumIm=%le\n",s1,s2);

//free(list1);
}
