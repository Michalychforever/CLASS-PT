#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h> 

#define PI 3.14159265358979323846
#define q0 75122.6331530 
#define q1 80916.6278952
#define q2 36308.2951477
#define q3 8687.24529705
#define q4 1168.92649479
#define q5 83.8676043424
#define q6 2.50662827511
#define g05 1.772453850905516
#define g15 0.886226925452758
#define gm05 -3.544907701811032
#define reg 0.00000001
#define normJ 0.02244839026564582
#define eps2F1 10E-10
#define epsJ2F1 10E-8
#define fourPI2 39.47841760435743

double complex Pochh(double complex a, double complex b);
double complex Gamma(double complex z);
double complex GammaHelp(double complex z);
double complex Hyper2F1basic(double complex a, double complex b, double complex c, double x);
double complex Hyper2F1(double complex a, double complex b, double complex c, double x);
double complex Hyper2F1recursion(double complex a, double complex b, double complex c, double x, int n, double complex f0, double complex f1);
double complex J2F1(double complex n1, double complex n2, double complex n3, double x, double y);
double complex J2F1fast(double complex n1, double complex n2, double complex n3, double x, double y);
double complex J2F1superfast(double complex n1, double complex n2, double complex n3, double x, double y);


int main()
{
    double x,y,p, bias, k0, kmax, Delta, test1;
    double complex z1, z2, z2p3, nu1, nu2, nu3, test;
    int i, j, broj, maxbroj, i1, i2, i3, j1, j2, l1, l2, l3, Nmax;
    double h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12;
	double B222Tab[72][12];
	double etam[80], b222[50][50][50][2];
	double complex fa[30], fb[30], fa0,fa1,fa2, fb0,fb1,fb2;
	
	FILE *fp0 = fopen("Input.dat", "r"); 
	fscanf(fp0,"%lf %lf %lf %d %lf %lf\n",&x,&y,&bias,&Nmax,&k0,&kmax);
	fclose(fp0);
	
	Delta=log(kmax/k0)/(Nmax-1);
	for(i=0;i<=Nmax;i++){
	etam[i]=2*PI*(i-Nmax/2.0)/(Nmax*Delta);}
    
	/********************************************/
	/*  Importing the file for B222             */
	/********************************************/
	FILE *fp1 = fopen("B222Tab.dat", "r"); 
	for(i1=0;i1<72;i1++){
	fscanf(fp1,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&h1,&h2,&h3,&h4,&h5,&h6,&h7,&h8,&h9,&h10,&h11,&h12);
	B222Tab[i1][0]=h1; B222Tab[i1][1]=h2; B222Tab[i1][2]=h3; B222Tab[i1][3]=h4;
	B222Tab[i1][4]=h5; B222Tab[i1][5]=h6; B222Tab[i1][6]=h7; B222Tab[i1][7]=h8;
	B222Tab[i1][8]=h9; B222Tab[i1][9]=h10; B222Tab[i1][10]=h11; B222Tab[i1][11]=h12; }
	fclose(fp1);
	/********************************************/
	
	
/*	test=0;
	for(i1=1;i1<=10000;i1++){
	test=1000*J2F1superfast(-0.1+0.423*I,0.432-0.3213*I,0.144+1.231*I,0.8,0.9);
	}
	printf("%lf + I %lf\n",creal(test),cimag(test)); */
	
	double complex p1, p2, s, sold, n1, n2, n3, a, b, c;
	double complex Fmatrix1[5][5][7][30], Fmatrix2[5][5][7][30], Cmatrix1[5][5][7][30], Cmatrix2[5][5][7][30];
	double eps;
	int n, nmax, ind11, ind12, ind13, ind21, ind22, ind23, k;


		
	/********************************************/
	/*  Calculating B222                        */
	/********************************************/	
	for(i1=0;i1<Nmax+1;i1++){
	for(i2=0;i2<Nmax+1;i2++){
	for(i3=0;i3<Nmax/2.0+1;i3++){
	
	n1=-bias/2-etam[i1]*I/2+reg;
	n2=-bias/2-etam[i2]*I/2+0.9*reg;
	n3=-bias/2-etam[i3]*I/2+0.8*reg;

	p1 = g15*Gamma(1.5-n1)*Gamma(n1+n2+n3-1.5)*Gamma(1.5-n2-n3)*cpow(y,1.5-n1-n3)/(fourPI2*Gamma(n1)*Gamma(n2+n3)*Gamma(3-n1-n2-n3)); 
	p2 = g15*Gamma(n2+n3-1.5)*Gamma(1.5-n2)*Gamma(1.5-n3)*cpow(x,1.5-n2-n3)/(fourPI2*Gamma(n2)*Gamma(n3)*Gamma(3-n2-n3));

	Cmatrix1[2][2][2][0]=p1;
	Cmatrix2[2][2][2][0]=p2;

	for(l1=26;l1>=0;l1=l1-1){
	fa[l1]=0;
	fb[l1]=0;
	}

    n=25;
    	
	fa2=Hyper2F1basic(n2+n, 1.5-n1+n, n2+n3+2*n, 1-y);
	fa[n]=fa2;
	fa1=Hyper2F1basic(n2+n-1, 0.5-n1+n, n2+n3+2*n-2, 1-y);
	fa[n-1]=fa1;
	for(l1=n-2;l1>=0;l1=l1-1){
	fa0=Hyper2F1recursion(n2, 1.5-n1, n2+n3, 1-y, l1, fa1, fa2);
	fa[l1]=fa0;	
 	fa2=fa1;
 	fa1=fa0;
	}

	fb2=Hyper2F1basic(n1+n, 1.5-n2+n, 3-n2-n3+2*n, 1-y);
	fb[n]=fb2;
	fb1=Hyper2F1basic(n1+n-1, 0.5-n2+n, 1-n2-n3+2*n, 1-y);
	fb[n-1]=fb1;
	for(l1=n-2;l1>=0;l1=l1-1){
	fb0=Hyper2F1recursion(n1, 1.5-n2, 3-n2-n3, 1-y, l1, fb1, fb2);
	fb[l1]=fb0;	
 	fb2=fb1;
 	fb1=fb0;
	}

	s = 0;
	eps = 1.0;
	n = 0;
	while(eps > epsJ2F1)
	{
		sold = s;
		Fmatrix1[2][2][2][n]=fa[n]; //Hyper2F1basic(n2+n, 1.5-n1+n, n2+n3+2*n, 1-y);
		Fmatrix2[2][2][4][n]=fb[n]; //Hyper2F1basic(n1+n, 1.5-n2+n, 3-n2-n3+2*n, 1-y);
 		s = s + p1*Fmatrix1[2][2][2][n] + p2*Fmatrix2[2][2][4][n];
 		p1 = p1*((n+n2)*(n+n3)*(1.5+n-n1)*(-1.5+n+n1+n2+n3)*y)/((1+n)*(2*n+n2+n3)*(1+2*n+n2+n3)*(-0.5+n+n2+n3));
		p2 = p2*((1.5+n-n2)*(1.5+n-n3)*(3+n-n1-n2-n3)*(n+n1)*y)/((1+n)*(2.5+n-n2-n3)*(3+2*n-n2-n3)*(4+2*n-n2-n3));
		Cmatrix1[2][2][2][n+1]=p1*x/y;
 		Cmatrix2[2][2][2][n+1]=p2*x/y;
 		
		eps=sqrt(cpow(creal(s-sold),2)+cpow(cimag(s-sold),2))/sqrt(cpow(creal(s),2)+cpow(cimag(s),2));
		n=n+1;
	}			
	nmax=n+1;	
	
	n=nmax;
    	
	fa2=Hyper2F1basic(n2+n, 1.5-n1+n, n2+n3+2*n, 1-y);
	fa[n]=fa2;
	fa1=Hyper2F1basic(n2+n-1, 0.5-n1+n, n2+n3+2*n-2, 1-y);
	fa[n-1]=fa1;
	for(l1=n-2;l1>=0;l1=l1-1){
	fa0=Hyper2F1recursion(n2, 1.5-n1, n2+n3, 1-y, l1, fa1, fa2);
	fa[l1]=fa0;	
 	fa2=fa1;
 	fa1=fa0;
	}

	fb2=Hyper2F1basic(n1+n, 1.5-n2+n, 3-n2-n3+2*n, 1-y);
	fb[n]=fb2;
	fb1=Hyper2F1basic(n1+n-1, 0.5-n2+n, 1-n2-n3+2*n, 1-y);
	fb[n-1]=fb1;
	for(l1=n-2;l1>=0;l1=l1-1){
	fb0=Hyper2F1recursion(n1, 1.5-n2, 3-n2-n3, 1-y, l1, fb1, fb2);
	fb[l1]=fb0;	
 	fb2=fb1;
 	fb1=fb0;
	}

	for(l1=0;l1<nmax;l1++){
	Fmatrix1[3][2][2][l1]=Hyper2F1basic(n2+l1+1, 1.5-n1+l1, n2+n3+2*l1, 1-y);
	}
	
	for(n=0;n<nmax;n++){
		a=n2+n;
		b=1.5-n1+n;
		c=n2+n3+2*n;
		Fmatrix1[2][1][2][n]=1/(b-c)*(-y*a*Fmatrix1[3][2][2][n]-(c-a-b)*Fmatrix1[2][2][2][n]);
		Fmatrix1[4][2][2][n]=1/((a+1)*y)*((c-a-1)*Fmatrix1[2][2][2][n] + (1-c+a+b-(b-a-1)*y)*Fmatrix1[3][2][2][n]);
		Fmatrix1[3][1][2][n]=1/(b-c)*(-y*(a+1)*Fmatrix1[4][2][2][n]-(c-a-b-1)*Fmatrix1[3][2][2][n]);
		Fmatrix1[4][1][2][n]=1/((a+1)*y)*((c-a-1)*Fmatrix1[2][1][2][n] + (-c+a+b-(b-a-2)*y)*Fmatrix1[3][1][2][n]);
		Fmatrix1[1][2][2][n]=-1/(a-c)*((c-2*a+(a-b)*(1-y))*Fmatrix1[2][2][2][n] + a*y*Fmatrix1[3][2][2][n]);
		Fmatrix1[1][1][2][n]=-1/(a-c)*((c-2*a+(a-b+1)*(1-y))*Fmatrix1[2][1][2][n] + a*y*Fmatrix1[3][1][2][n]);
		Fmatrix1[0][2][2][n]=-1/(a-1-c)*((c-2*a+2+(a-b-1)*(1-y))*Fmatrix1[1][2][2][n] + (a-1)*y*Fmatrix1[2][2][2][n]);
		Fmatrix1[0][1][2][n]=-1/(a-1-c)*((c-2*a+2+(a-b)*(1-y))*Fmatrix1[1][1][2][n] + (a-1)*y*Fmatrix1[2][1][2][n]);
		for(i=-2;i<3;i++){
		Fmatrix1[2+i][0][2][n]= 1/(-1+b-c)*((1-b)*y*Fmatrix1[2+i][2][2][n] - (2-2*b+c+(-1-a+b-i)*(1-y))*Fmatrix1[2+i][1][2][n]);
		Fmatrix1[2+i][3][2][n]=1/(b*y)*((c-b)*Fmatrix1[2+i][1][2][n] - (c-2*b+(b-a-i)*(1-y))*Fmatrix1[2+i][2][2][n]);
		Fmatrix1[2+i][4][2][n]=1/((b+1)*y)*((c-b-1)*Fmatrix1[2+i][2][2][n] - (c-2*b-2+(b-a-i+1)*(1-y))*Fmatrix1[2+i][3][2][n]);
		}
		Fmatrix1[2][2][3][n]=-c/((1-y)*(a-c)*(c-b))*(((c-b)*(1-y)-a)*Fmatrix1[2][2][2][n] + y*a*Fmatrix1[3][2][2][n]);
		Fmatrix1[3][2][3][n]=-c/((1-y)*(a-c+1)*(c-b))*(((c-b)*(1-y)-a-1)*Fmatrix1[3][2][2][n] + y*(a+1)*Fmatrix1[4][2][2][n]);
		c=n2+n3+2*n+1;
		Fmatrix1[2][1][3][n]=1/(b-c)*(-y*a*Fmatrix1[3][2][3][n]-(c-a-b)*Fmatrix1[2][2][3][n]);
		Fmatrix1[4][2][3][n]=1/((a+1)*y)*((c-a-1)*Fmatrix1[2][2][3][n] + (1-c+a+b-(b-a-1)*y)*Fmatrix1[3][2][3][n]);
		Fmatrix1[3][1][3][n]=1/(b-c)*(-y*(a+1)*Fmatrix1[4][2][3][n]-(c-a-b-1)*Fmatrix1[3][2][3][n]);
		Fmatrix1[4][1][3][n]=1/((a+1)*y)*((c-a-1)*Fmatrix1[2][1][3][n] + (-c+a+b-(b-a-2)*y)*Fmatrix1[3][1][3][n]);
		Fmatrix1[1][2][3][n]=-1/(a-c)*((c-2*a+(a-b)*(1-y))*Fmatrix1[2][2][3][n] + a*y*Fmatrix1[3][2][3][n]);
		Fmatrix1[1][1][3][n]=-1/(a-c)*((c-2*a+(a-b+1)*(1-y))*Fmatrix1[2][1][3][n] + a*y*Fmatrix1[3][1][3][n]);
		Fmatrix1[0][2][3][n]=-1/(a-1-c)*((c-2*a+2+(a-b-1)*(1-y))*Fmatrix1[1][2][3][n] + (a-1)*y*Fmatrix1[2][2][3][n]);
		Fmatrix1[0][1][3][n]=-1/(a-1-c)*((c-2*a+2+(a-b)*(1-y))*Fmatrix1[1][1][3][n] + (a-1)*y*Fmatrix1[2][1][3][n]);
		for(l1=-2;l1<3;l1++){
		Fmatrix1[2+l1][0][3][n]=1/(-1+b-c)*((1-b)*y*Fmatrix1[2+l1][2][3][n] - (2-2*b+c+(-1-a+b-l1)*(1-y))*Fmatrix1[2+l1][1][3][n]);
		Fmatrix1[2+l1][3][3][n]=1/(b*y)*((c-b)*Fmatrix1[2+l1][1][3][n] - (c-2*b+(b-a-l1)*(1-y))*Fmatrix1[2+l1][2][3][n]);
		Fmatrix1[2+l1][4][3][n]=1/((b+1)*y)*((c-b-1)*Fmatrix1[2+l1][2][3][n] - (c-2*b-2+(b-a-l1+1)*(1-y))*Fmatrix1[2+l1][3][3][n]);
		}
		for(l3=0;l3<3;l3++){
		c=n2+n3+2*n+l3;
		for(l1=-2;l1<3;l1++){
		for(l2=-2;l2<3;l2++){
		Fmatrix1[2+l1][2+l2][4+l3][n]=(1+c)/((1-a-l1+c)*(1-b-l2+c)*(1-y))*(c*y*Fmatrix1[2+l1][2+l2][2+l3][n] + (-c+(1-a-b-l1-l2+2*c)*(1-y))*Fmatrix1[2+l1][2+l2][3+l3][n]); 
		}}}
		for(l3=0;l3<2;l3++){
		c=n2+n3+2*n-l3;
		for(l1=-2;l1<3;l1++){
		for(l2=-2;l2<3;l2++){
		Fmatrix1[2+l1][2+l2][1-l3][n]=1/((c-1)*c*y)*(-c*(1-c+(2*c-a-b-1-l1-l2)*(1-y))*Fmatrix1[2+l1][2+l2][2-l3][n] + (a-c+l1)*(b-c+l2)*(1-y)*Fmatrix1[2+l1][2+l2][3-l3][n]);
		}}}
		}
	
	for(l1=0;l1<nmax;l1++){
	Fmatrix2[3][2][4][l1]=Hyper2F1basic(n1+l1+1, 1.5-n2+l1, 3-n2-n3+2*l1, 1-y);
	}
	
	for(n=0;n<nmax;n++){
		a=n1+n;
		b=1.5-n2+n;
		c=3-n2-n3+2*n;
		Fmatrix2[2][1][4][n]=1/(b-c)*(-y*a*Fmatrix2[3][2][4][n]-(c-a-b)*Fmatrix2[2][2][4][n]);
		Fmatrix2[4][2][4][n]=1/((a+1)*y)*((c-a-1)*Fmatrix2[2][2][4][n] + (1-c+a+b-(b-a-1)*y)*Fmatrix2[3][2][4][n]);
		Fmatrix2[3][1][4][n]=1/(b-c)*(-y*(a+1)*Fmatrix2[4][2][4][n]-(c-a-b-1)*Fmatrix2[3][2][4][n]);
		Fmatrix2[4][1][4][n]=1/((a+1)*y)*((c-a-1)*Fmatrix2[2][1][4][n] + (-c+a+b-(b-a-2)*y)*Fmatrix2[3][1][4][n]);
		Fmatrix2[1][2][4][n]=-1/(a-c)*((c-2*a+(a-b)*(1-y))*Fmatrix2[2][2][4][n] + a*y*Fmatrix2[3][2][4][n]);
		Fmatrix2[1][1][4][n]=-1/(a-c)*((c-2*a+(a-b+1)*(1-y))*Fmatrix2[2][1][4][n] + a*y*Fmatrix2[3][1][4][n]);
		Fmatrix2[0][2][4][n]=-1/(a-1-c)*((c-2*a+2+(a-b-1)*(1-y))*Fmatrix2[1][2][4][n] + (a-1)*y*Fmatrix2[2][2][4][n]);
		Fmatrix2[0][1][4][n]=-1/(a-1-c)*((c-2*a+2+(a-b)*(1-y))*Fmatrix2[1][1][4][n] + (a-1)*y*Fmatrix2[2][1][4][n]);
		for(i=-2;i<3;i++){
		Fmatrix2[2+i][0][4][n]= 1/(-1+b-c)*((1-b)*y*Fmatrix2[2+i][2][4][n] - (2-2*b+c+(-1-a+b-i)*(1-y))*Fmatrix2[2+i][1][4][n]);
		Fmatrix2[2+i][3][4][n]=1/(b*y)*((c-b)*Fmatrix2[2+i][1][4][n] - (c-2*b+(b-a-i)*(1-y))*Fmatrix2[2+i][2][4][n]);
		Fmatrix2[2+i][4][4][n]=1/((b+1)*y)*((c-b-1)*Fmatrix2[2+i][2][4][n] - (c-2*b-2+(b-a-i+1)*(1-y))*Fmatrix2[2+i][3][4][n]);
		}
		Fmatrix2[2][2][5][n]=-c/((1-y)*(a-c)*(c-b))*(((c-b)*(1-y)-a)*Fmatrix2[2][2][4][n] + y*a*Fmatrix2[3][2][4][n]);
		Fmatrix2[3][2][5][n]=-c/((1-y)*(a-c+1)*(c-b))*(((c-b)*(1-y)-a-1)*Fmatrix2[3][2][4][n] + y*(a+1)*Fmatrix2[4][2][4][n]);
		c=3-n2-n3+2*n+1;
		Fmatrix2[2][1][5][n]=1/(b-c)*(-y*a*Fmatrix2[3][2][5][n]-(c-a-b)*Fmatrix2[2][2][5][n]);
		Fmatrix2[4][2][5][n]=1/((a+1)*y)*((c-a-1)*Fmatrix2[2][2][5][n] + (1-c+a+b-(b-a-1)*y)*Fmatrix2[3][2][5][n]);
		Fmatrix2[3][1][5][n]=1/(b-c)*(-y*(a+1)*Fmatrix2[4][2][5][n]-(c-a-b-1)*Fmatrix2[3][2][5][n]);
		Fmatrix2[4][1][5][n]=1/((a+1)*y)*((c-a-1)*Fmatrix2[2][1][5][n] + (-c+a+b-(b-a-2)*y)*Fmatrix2[3][1][5][n]);
		Fmatrix2[1][2][5][n]=-1/(a-c)*((c-2*a+(a-b)*(1-y))*Fmatrix2[2][2][5][n] + a*y*Fmatrix2[3][2][5][n]);
		Fmatrix2[1][1][5][n]=-1/(a-c)*((c-2*a+(a-b+1)*(1-y))*Fmatrix2[2][1][5][n] + a*y*Fmatrix2[3][1][5][n]);
		Fmatrix2[0][2][5][n]=-1/(a-1-c)*((c-2*a+2+(a-b-1)*(1-y))*Fmatrix2[1][2][5][n] + (a-1)*y*Fmatrix2[2][2][5][n]);
		Fmatrix2[0][1][5][n]=-1/(a-1-c)*((c-2*a+2+(a-b)*(1-y))*Fmatrix2[1][1][5][n] + (a-1)*y*Fmatrix2[2][1][5][n]);
		for(l1=-2;l1<3;l1++){
		Fmatrix2[2+l1][0][5][n]=1/(-1+b-c)*((1-b)*y*Fmatrix2[2+l1][2][5][n] - (2-2*b+c+(-1-a+b-l1)*(1-y))*Fmatrix2[2+l1][1][5][n]);
		Fmatrix2[2+l1][3][5][n]=1/(b*y)*((c-b)*Fmatrix2[2+l1][1][5][n] - (c-2*b+(b-a-l1)*(1-y))*Fmatrix2[2+l1][2][5][n]);
		Fmatrix2[2+l1][4][5][n]=1/((b+1)*y)*((c-b-1)*Fmatrix2[2+l1][2][5][n] - (c-2*b-2+(b-a-l1+1)*(1-y))*Fmatrix2[2+l1][3][5][n]);
		}
		c=3-n2-n3+2*n;
		for(l1=-2;l1<3;l1++){
		for(l2=-2;l2<3;l2++){
		Fmatrix2[2+l1][2+l2][6][n]=(1+c)/((1-a-l1+c)*(1-b-l2+c)*(1-y))*(c*y*Fmatrix2[2+l1][2+l2][4][n] + (-c+(1-a-b-l1-l2+2*c)*(1-y))*Fmatrix2[2+l1][2+l2][5][n]); 
		}}
		for(l3=0;l3<4;l3++){
		c=3-n2-n3+2*n-l3;
		for(l1=-2;l1<3;l1++){
		for(l2=-2;l2<3;l2++){
		Fmatrix2[2+l1][2+l2][3-l3][n]=1/((c-1)*c*y)*(-c*(1-c+(2*c-a-b-1-l1-l2)*(1-y))*Fmatrix2[2+l1][2+l2][4-l3][n] + (a-c+l1)*(b-c+l2)*(1-y)*Fmatrix2[2+l1][2+l2][5-l3][n]);
		}}}
		} 	

	for(i=0;i<nmax;i++){
	Cmatrix1[1][2][2][i]=((1.5-n1+i)*y*(n1-1)/((n1+n2+n3-2.5)*(3-n1-n2-n3)))*Cmatrix1[2][2][2][i];
	Cmatrix1[0][2][2][i]=((2.5-n1)*y*(n1-2)/((n1+n2+n3-3.5)*(4-n1-n2-n3)))*Cmatrix1[1][2][2][i];
	Cmatrix1[3][2][2][i]=((n1+n2+n3-1.5)*(2-n1-n2-n3)/((0.5-n1)*y*n1))*Cmatrix1[2][2][2][i];
	Cmatrix1[4][2][2][i]=((n1+n2+n3-0.5)*(1-n1-n2-n3)/((-0.5-n1)*y*(n1+1)))*Cmatrix1[3][2][2][i];
	}
	
	for(i=0;i<nmax;i++){
	for(j=0;j<5;j++){
	Cmatrix1[j][1][2][i]=((1.5-n2-n3)*(n2+n3-1)/((n1+n2+n3-2.5)*(3-n1-n2-n3)))*Cmatrix1[j][2][2][i];
	Cmatrix1[j][0][2][i]=((2.5-n2-n3)*(n2+n3-2)/((n1+n2+n3-3.5)*(4-n1-n2-n3)))*Cmatrix1[j][1][2][i];
	Cmatrix1[j][3][2][i]=((n1+n2+n3-1.5)*(2-n1-n2-n3)/((0.5-n2-n3)*(n2+n3)))*Cmatrix1[j][2][2][i];
	Cmatrix1[j][4][2][i]=((n1+n2+n3-0.5)*(1-n1-n2-n3)/((-0.5-n2-n3)*(n2+n3+1)))*Cmatrix1[j][3][2][i];
	}}
	
//	p1 = Gamma(1.5-n1)*Gamma(n1+n2+n3-1.5)*Gamma(1.5-n2-n3)*cpow(y,1.5-n1-n3)/(Gamma(n1)*Gamma(n2+n3)*Gamma(3-n1-n2-n3)); 	
	
	for(i=0;i<nmax;i++){
	for(j=0;j<5;j++){
	for(k=0;k<5;k++){
	Cmatrix1[j][k][1][i]=((1.5-n2-n3)*y*(n2+n3-1)/((n1+n2+n3-2.5)*(3-n1-n2-n3)))*Cmatrix1[j][k][2][i];
	Cmatrix1[j][k][0][i]=((2.5-n2-n3)*(n2+n3-2)/((n1+n2+n3-3.5)*(4-n1-n2-n3)))*Cmatrix1[j][k][1][i];
	Cmatrix1[j][k][3][i]=((n1+n2+n3-1.5)*(2-n1-n2-n3)/((0.5-n2-n3)*y*(n2+n3)))*Cmatrix1[j][k][2][i];
	Cmatrix1[j][k][4][i]=((n1+n2+n3-0.5)*(1-n1-n2-n3)/((-0.5-n2-n3)*y*(n2+n3+1)))*Cmatrix1[j][k][3][i];
	Cmatrix1[j][k][5][i]=((n1+n2+n3+0.5)*(-n1-n2-n3)/((-1.5-n2-n3)*y*(n2+n3+2)))*Cmatrix1[j][k][4][i];
	Cmatrix1[j][k][6][i]=((n1+n2+n3+1.5)*(-1-n1-n2-n3)/((-2.5-n2-n3)*y*(n2+n3+3)))*Cmatrix1[j][k][5][i];
	}}}
	
//	p2 = Gamma(n2+n3-1.5)*Gamma(1.5-n2)*Gamma(1.5-n3)*cpow(x,1.5-n2-n3)/(Gamma(n2)*Gamma(n3)*Gamma(3-n2-n3));  
	
	
	
	z1=0;
	for(i=0;i<72;i++){
	z2=B222Tab[i][3]+x*B222Tab[i][4]+y*B222Tab[i][5]+x*y*B222Tab[i][6]+x*x*B222Tab[i][7]+y*y*B222Tab[i][8]+x*x*y*B222Tab[i][9]+y*y*x*B222Tab[i][10]+x*x*y*y*B222Tab[i][11];
	
	n1=-bias/2-etam[i1]*I/2+B222Tab[i][0]+reg;
	n2=-bias/2-etam[i2]*I/2+B222Tab[i][1]+0.9*reg;
	n3=-bias/2-etam[i3]*I/2+B222Tab[i][2]+0.8*reg;
	
	ind11=(int) (B222Tab[i][1]+2.1);
	ind12=(int) (-B222Tab[i][0]+2.1);
	ind13=(int) (B222Tab[i][1]+B222Tab[i][2]+2.1);
	
	ind21=(int) (B222Tab[i][0]+2.1);
	ind22=(int) (-B222Tab[i][1]+2.1);
	ind23=(int) (-B222Tab[i][1]-B222Tab[i][2]+4.1);
	
	p1 = g15*Gamma(1.5-n1)*Gamma(n1+n2+n3-1.5)*Gamma(1.5-n2-n3)*cpow(y,1.5-n1-n3)/(fourPI2*Gamma(n1)*Gamma(n2+n3)*Gamma(3-n1-n2-n3)); 
	p2 = g15*Gamma(n2+n3-1.5)*Gamma(1.5-n2)*Gamma(1.5-n3)*cpow(x,1.5-n2-n3)/(fourPI2*Gamma(n2)*Gamma(n3)*Gamma(3-n2-n3));

	s = 0;
	eps = 1.0;
	n = 0;
	while(eps > epsJ2F1)
	{
		sold = s;
 		s = s + p1*Fmatrix1[ind11][ind12][ind13][n] + p2*Fmatrix2[ind21][ind22][ind23][n]; 
 		p1 = p1*((n+n2)*(n+n3)*(1.5+n-n1)*(-1.5+n+n1+n2+n3)*x)/((1+n)*(2*n+n2+n3)*(1+2*n+n2+n3)*(-0.5+n+n2+n3));
		p2 = p2*((1.5+n-n2)*(1.5+n-n3)*(3+n-n1-n2-n3)*(n+n1)*x)/((1+n)*(2.5+n-n2-n3)*(3+2*n-n2-n3)*(4+2*n-n2-n3));
 
		eps=sqrt(cpow(creal(s-sold),2)+cpow(cimag(s-sold),2))/sqrt(cpow(creal(s),2)+cpow(cimag(s),2));
		n=n+1;
	}		
		
	z1=z1+z2*s;
	};
	b222[i1][i2][i3][0]=creal(z1);
	b222[i1][i2][i3][1]=cimag(z1);   
	}};
	printf("%lf percent of B222\n",i1*100.0/Nmax);}	    
	for(i1=0;i1<Nmax+1;i1++){
	for(i2=0;i2<Nmax+1;i2++){
	for(i3=Nmax/2.0+1;i3<Nmax+1;i3++){
	b222[i1][i2][i3][0]=b222[Nmax-i1][Nmax-i2][Nmax-i3][0];
	b222[i1][i2][i3][1]=-b222[Nmax-i1][Nmax-i2][Nmax-i3][1];  
	}}}	  
	/********************************************/    
	printf("%d\n",nmax);

	
	
	/********************************************/
	/*  Calculating B222                        */
	/********************************************/	
/*	for(i1=0;i1<Nmax+1;i1++){
	for(i2=0;i2<Nmax+1;i2++){
	for(i3=0;i3<Nmax/2.0+1;i3++){
	z1=0;
	for(i=0;i<72;i++){
	z2=B222Tab[i][3]+x*B222Tab[i][4]+y*B222Tab[i][5]+x*y*B222Tab[i][6]+x*x*B222Tab[i][7]+y*y*B222Tab[i][8]+x*x*y*B222Tab[i][9]+y*y*x*B222Tab[i][10]+x*x*y*y*B222Tab[i][11];
	z1=z1+z2*J2F1fast(-bias/2-etam[i1]*I/2+B222Tab[i][0]+reg,-bias/2-etam[i2]*I/2+B222Tab[i][1]+0.9*reg,-bias/2-etam[i3]*I/2+B222Tab[i][2]+0.8*reg,x,y);
	};
	b222[i1][i2][i3][0]=creal(z1);
	b222[i1][i2][i3][1]=cimag(z1);
	}};
	printf("%lf percent of B222\n",i1*100.0/Nmax);}	    
	for(i1=0;i1<Nmax+1;i1++){
	for(i2=0;i2<Nmax+1;i2++){
	for(i3=Nmax/2.0+1;i3<Nmax+1;i3++){
	b222[i1][i2][i3][0]=b222[Nmax-i1][Nmax-i2][Nmax-i3][0];
	b222[i1][i2][i3][1]=-b222[Nmax-i1][Nmax-i2][Nmax-i3][1];
	}}}	 */	
	/********************************************/    


	/********************************************/
	/*  Output for B222                         */
	/********************************************/		
	FILE *gp1 = fopen("outb222.dat", "w");
	for(i1=0;i1<Nmax+1;i1++){
	for(i2=0;i2<Nmax+1;i2++){
	for(i3=0;i3<Nmax+1;i3++){
	fprintf(gp1, "%lf %lf\n", b222[i1][i2][i3][0],b222[i1][i2][i3][1]);
	}}}
	fclose(gp1);   
	/********************************************/

    return 0;
}

// Definitions of functions

double complex GammaHelp(double complex z)
{
	return (q0+q1*z+q2*cpow(z,2)+q3*cpow(z,3)+q4*cpow(z,4)+q5*cpow(z,5)+q6*cpow(z,6))/(z*(z+1)*(z+2)*(z+3)*(z+4)*(z+5)*(z+6));
}

double complex Gamma(double complex z)
{
	double complex result;
	
	if(creal(z)>0){
		result=GammaHelp(z) * cpow(z+5.5,z+0.5)* cexp(-z-5.5);
	}
	else
	{
		double complex help=GammaHelp(1-z) * cpow(1-z+5.5,1-z+0.5)* cexp(-1+z-5.5);
		result=PI/(csin(PI*z)*help);
	}
	return result;
}

double complex Pochh(double complex a,double complex b) 
{
    double complex result;
    result = Gamma(a+b)/Gamma(a);
    return result;                  
}

double complex Hyper2F1basic(double complex a, double complex b, double complex c, double x)
{
	double complex result, s, p;
	double eps;
	int n;
	
	if(fabs(x)<=1){
	p=1;
	s=1;
	n=1;
	eps=1.0;
	while(eps>eps2F1)
	{
		p=p*(a+n-1)*(b+n-1)*x/(n*(c+n-1));
		s=s+p;
		n=n+1;
		eps=sqrt((creal(p)*creal(p)+cimag(p)*cimag(p))/(creal(s)*creal(s)+cimag(s)*cimag(s)));
	}
	result=s;
	}
	else result=1;
	
	return result;
}

double complex Hyper2F1(double complex a, double complex b, double complex c, double x)
{
	double complex result;
		
	if(fabs(x)<=1){
	result=Hyper2F1basic(a,b,c,x);
	}
	else result=Gamma(c)*Gamma(b-a)*cpow(-x,-a)*Hyper2F1basic(a,a+1-c,a+1-b,1/x)/(Gamma(b)*Gamma(c-a))+Gamma(c)*Gamma(a-b)*cpow(-x,-b)*Hyper2F1basic(b,b+1-c,b+1-a,1/x)/(Gamma(a)*Gamma(c-b));
	
	return result;
}

double complex J2F1(double complex n1, double complex n2, double complex n3, double x, double y)
{
	double complex result, N1, N2, p1, p2, s1, s2, s1old, s2old;
	double eps;
	int n;
	
	N1 = Gamma(1.5)*Gamma(1.5-n3)*Gamma(n1+n2+n3-1.5)*Gamma(1.5-n1-n2)*cpow(y,1.5-n1-n2-n3)/(4*PI*PI*Gamma(n1)*Gamma(n2)*Gamma(n3)*Gamma(3-n1-n2-n3));
	N2 = Gamma(1.5)*Gamma(n1+n2-1.5)*cpow(y,-n3)/(4*PI*PI*Gamma(n1)*Gamma(n2));
	
	p1 = Gamma(n1)*Gamma(n2)/Gamma(n1+n2); 
	p2 = Gamma(1.5-n1)*Gamma(1.5-n2)/Gamma(3-n1-n2);
	s1 = 0;
	s2 = 0;
	eps = 1.0;
	n = 0;
	
	while(eps > epsJ2F1)
	{
		s1old = s1;
 		s1 = s1 + p1*Hyper2F1(n2+n, -1.5+n1+n2+n3+n, n1+n2+2*n, 1-x/y);
 		p1 = p1*((n+n1)*(n+n2)*(3+2*n-2*n3)*(-3+2*n+2*n1+2*n2+2*n3))/(y*2*(1+n)*(2*n+n1+n2)*(1+2*n+n1+n2)*(-1+2*n+2*n1+2*n2));
 
		s2old = s2;
		s2 = s2 + p2*Hyper2F1(1.5-n1+n, n3+n, 3-n1-n2+2*n, 1-x/y);
		p2 = p2*((3+2*n-2*n1)*(3+2*n-2*n2)*(3+n-n1-n2-n3)*(n+n3))/(y*2*(1+n)*(5+2*n-2*n1-2*n2)*(3+2*n-n1-n2)*(4+2*n-n1-n2));
 
		eps=sqrt(cpow(creal(N1*(s1-s1old)+N2*(s2-s2old)),2)+cpow(cimag(N1*(s1-s1old)+N2*(s2-s2old)),2))/sqrt(cpow(creal(N1*s1+N2*s2),2)+cpow(cimag(N1*s1+N2*s2),2));
		n=n+1;
	}
		
	result = N1*s1+N2*s2;
	return result;
}


double complex J2F1fast(double complex n1, double complex n2, double complex n3, double x, double y)
{
	double complex result, N1, N2, p1, p2, s1, s2, s1old, s2old;
	double eps;
	int n;
	
	N1 = g15*Gamma(1.5-n1)*Gamma(n1+n2+n3-1.5)*Gamma(1.5-n2-n3)/(4*PI*PI*Gamma(n2)*Gamma(n3)*Gamma(n1)*Gamma(3-n1-n2-n3));
	N2 = g15*Gamma(n2+n3-1.5)*cpow(x,1.5-n2-n3)/(4*PI*PI*Gamma(n2)*Gamma(n3));
	
	p1 = Gamma(n2)*Gamma(n3)/Gamma(n2+n3); 
	p2 = Gamma(1.5-n2)*Gamma(1.5-n3)/Gamma(3-n2-n3);
	s1 = 0;
	s2 = 0;
	eps = 1.0;
	n = 0;
	
	while(eps > epsJ2F1)
	{
		s1old = s1;
 		s1 = s1 + p1*Hyper2F1basic(n3+n, -1.5+n1+n2+n3+n, n2+n3+2*n, 1-y);
 		p1 = p1*((n+n2)*(n+n3)*(3+2*n-2*n1)*(-3+2*n+2*n1+2*n2+2*n3)*x)/(2*(1+n)*(2*n+n2+n3)*(1+2*n+n2+n3)*(-1+2*n+2*n2+2*n3));
 
		s2old = s2;
		s2 = s2 + p2*Hyper2F1basic(1.5-n2+n, n1+n, 3-n2-n3+2*n, 1-y);
		p2 = p2*((3+2*n-2*n2)*(3+2*n-2*n3)*(3+n-n1-n2-n3)*(n+n1)*x)/(2*(1+n)*(5+2*n-2*n2-2*n3)*(3+2*n-n2-n3)*(4+2*n-n2-n3));
 
		eps=sqrt(cpow(creal(N1*(s1-s1old)+N2*(s2-s2old)),2)+cpow(cimag(N1*(s1-s1old)+N2*(s2-s2old)),2))/sqrt(cpow(creal(N1*s1+N2*s2),2)+cpow(cimag(N1*s1+N2*s2),2));
		n=n+1;
	}
		
	result = N1*s1+N2*s2;
	return result;
}


double complex Hyper2F1recursion(double complex a, double complex b, double complex c, double x, int n, double complex f1, double complex f2)
{
	double complex result, N, p, q, l;
	
	q=(c+2*n+2)*(c+2*n+2);
	l=(c-a+n+1)*(c-b+n+1);
	N= ((a+n+1)*(b+n+1)*l*x*x)/(q*(q-1));
	p= (a+n)*(b+n)/(c+2*n)+l/(c+2*n+2);
	
	result=(1-x+x*p/(c+2*n+1))*f1-f2*N;
	
	return result;
}

double complex J2F1superfast(double complex n1, double complex n2, double complex n3, double x, double y)
{
	double complex result, N1, N2, p1, p2, s1, s2, s1old, s2old, fa0, fa1, fa2, fb0, fb1, fb2;
	double eps;
	int n,i1;
	double complex fa[50], fb[50];
	
	N1 = g15*Gamma(1.5-n1)*Gamma(n1+n2+n3-1.5)*Gamma(1.5-n2-n3)/(4*PI*PI*Gamma(n2)*Gamma(n3)*Gamma(n1)*Gamma(3-n1-n2-n3));
	N2 = g15*Gamma(n2+n3-1.5)*cpow(x,1.5-n2-n3)/(4*PI*PI*Gamma(n2)*Gamma(n3));
	
	p1 = Gamma(n2)*Gamma(n3)/Gamma(n2+n3); 
	p2 = Gamma(1.5-n2)*Gamma(1.5-n3)/Gamma(3-n2-n3);
	s1 = 0;
	s2 = 0;
	eps = 1.0;
	
	n=11;
	
	fa2=Hyper2F1basic(n3+n, -1.5+n1+n2+n3+n, n2+n3+2*n, 1-y);
	fa[n]=fa2;
	fa1=Hyper2F1basic(n3+n-1, -2.5+n1+n2+n3+n, n2+n3+2*n-2, 1-y);
	fa[n-1]=fa1;
	for(i1=n-2;i1>=0;i1=i1-1){
	fa0=Hyper2F1recursion(n3,-1.5+n1+n2+n3,n2+n3,1-y,i1,fa1,fa2);
	fa[i1]=fa0;	
 	fa2=fa1;
 	fa1=fa0;
	}

	fb2=Hyper2F1basic(1.5-n2+n, n1+n, 3-n2-n3+2*n, 1-y);
	fb[n]=fb2;
	fb1=Hyper2F1basic(0.5-n2+n, n1+n-1, 1-n2-n3+2*n, 1-y);
	fb[n-1]=fb1;
	for(i1=n-2;i1>=0;i1=i1-1){
	fb0=Hyper2F1recursion(1.5-n2,n1,3-n2-n3,1-y,i1,fb1,fb2);
	fb[i1]=fb0;	
 	fb2=fb1;
 	fb1=fb0;
	}
		
	n=0;
	
	while(eps > epsJ2F1)
	{
		s1old = s1;
 		s1 = s1 + p1*Hyper2F1basic(n3+n, -1.5+n1+n2+n3+n, n2+n3+2*n, 1-y);
 		p1 = p1*((n+n2)*(n+n3)*(1.5+n-n1)*(-1.5+n+n1+n2+n3)*x)/((1+n)*(2*n+n2+n3)*(1+2*n+n2+n3)*(-0.5+n+n2+n3));
 
		s2old = s2;
		s2 = s2 + p2*Hyper2F1basic(1.5-n2+n, n1+n, 3-n2-n3+2*n, 1-y);
		p2 = p2*((1.5+n-n2)*(1.5+n-n3)*(3+n-n1-n2-n3)*(n+n1)*x)/((1+n)*(2.5+n-n2-n3)*(3+2*n-n2-n3)*(4+2*n-n2-n3));

 		if(n==20){
 		eps=epsJ2F1/10;
 		}
 		else{
 		eps=sqrt(cpow(creal(N1*(s1-s1old)+N2*(s2-s2old)),2)+cpow(cimag(N1*(s1-s1old)+N2*(s2-s2old)),2))/sqrt(cpow(creal(N1*s1+N2*s2),2)+cpow(cimag(N1*s1+N2*s2),2));
 		}
 		n=n+1;
	}
		
	result = N1*s1+N2*s2;
	return result;
}


