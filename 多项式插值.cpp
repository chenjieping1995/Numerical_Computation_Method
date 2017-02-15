/*
* @ Filename: 多项式插值.cpp 
* @ Author: Jieping Chen
* @ Student Number: PB14209115
* @ Version: 1.0.0
* @ Date: 2016.10.07
* @ Function: 拉格朗日多项式插值法与牛顿多项式插值法在等距分割和切比雪夫点分割下的对比 
*/

#include <math.h>
#include <stdio.h> 
#include <iostream>
#define PI 3.14159265358979323846

using namespace std;

typedef struct tagPOINT{
	double x;
	double y;
}POINT;

//被插函数
double func(double x){
	double f = 1.0/(1.0 + pow(x,2));
	return f;
} 


//等距分割
POINT* partition(int n){
	int i;
	POINT* point = new POINT[n+1];
	for (i=0; i<=n; i++){
		point[i].x = -1.0 + 2.0/n*i;
		point[i].y = func(point[i].x);
	}
	return point;
} 


//chebyshev分割 
POINT* chebyshev(int n){
	int i;
	POINT* point = new POINT[n+1];
	for (i=0; i<=n; i++){
		point[i].x = cos( ((2.0*i+1.0)*PI) / double(2.0*n+1.0));
		point[i].y = func(point[i].x);
	}
	return point;
} 


//求取最大值函数
double max(double* e, int n){
	int i;
	double x=0;
	for (i=0; i<=n; i++){
		if (x<e[i]) x=e[i];
	}
	return x;
} 


//求和函数
double sum(double *e, int n){
	int i;
	double x=0;
	for (i=0; i<=n; i++){
		x += e[i];
	}
	return x;
} 


//将插值点带入多项式计算
double calc_func(POINT* point, double* diff, double x, int n ){
	//计算newton(x)的值
	double tmp = 1;
	double newton = diff[0];
	int i;
	for (i=0; i<n; i++){
		tmp = tmp*(x - point[i].x);
		newton = newton + tmp*diff[i+1];
	}
	return newton;
}



//newton插值多项式有关计算 
void newton1(int n){
	//计算newton插值多项式 
	int i, j;
	POINT* point = new POINT[n+1];
	double* diff = new double[n+1];
	point = partition(n);
	for (i=0; i<=n; i++){
		diff[i] = point[i].y;
	}
	
	for (i=0; i<n; i++){
		for (j=n; j>i; j--){
			//计算f(x0,...,xn)的差商 
			diff[j] = (diff[j] - diff[j-1])/(point[j].x - point[j-1-i].x);
		}
	}
	
	//计算误差
	POINT* input = new POINT[n+1];
	double* error = new double[n+1]; 
	for (j=0; j<=n; j++){
		//计算每个插值点所对应的插值函数的值 
		input[j].x = -1 + (double)(j)/(n + 1);
		input[j].y = calc_func(point, diff, input[j].x, n);
		//计算误差
		error[j] = fabs(input[j].y - func(input[j].x));
	}
	
	double L, L1;
	L = max(error, n);
	cout<<"n:"<<n<<'\t';
	cout<<"the newton L& is:"<<L<<'\t';
	
	L1 = sum(error, n)/(1.0 + n);
	cout<<"the newton L1 is:"<<L1<<'\t';
	cout<<endl;
}

//
void newton2(int n){
	//计算newton插值多项式 
	int i, j;
	POINT* point = new POINT[n+1];
	double* diff = new double[n+1];
	point = chebyshev(n);
	for (i=0; i<=n; i++){
		diff[i] = point[i].y;
	}
	
	for (i=0; i<n; i++){
		for (j=n; j>i; j--){
			//计算f(x0,...,xn)的差商 
			diff[j] = (diff[j] - diff[j-1])/(point[j].x - point[j-1-i].x);
		}
	}
	
	//计算误差
	POINT* input = new POINT[n+1];
	double* error = new double[n+1]; 
	for (j=0; j<=n; j++){
		//计算每个插值点所对应的插值函数的值 
		input[j].x = -1 + (double)(j)/(n + 1);
		input[j].y = calc_func(point, diff, input[j].x, n);
		//计算误差
		error[j] = fabs(input[j].y - func(input[j].x));
	}
	
	double L, L1;
	L = max(error, n);
	cout<<"n:"<<n<<'\t';
	cout<<"the newton L& is:"<<L<<'\t';
	
	L1 = sum(error, n)/(1.0 + n);
	cout<<"the newton L1 is:"<<L1<<'\t';
	cout<<endl;
}


//计算lagrange插值函数值 
double calc_lagrange(POINT* point, double x, int n){
	int i,j;
	double tmp;
	double fx=0;
	for (i=0; i<=n; i++){
		tmp=1.0;
		for (j=0; j<=n; j++){
			if(i != j){
				tmp *= (x-point[j].x)/(point[i].x-point[j].x);
			} 
		}
		fx += tmp*point[i].y;
	}
	return fx;
}


//lagrange插值多项式
void lagrange1(int n){
	//等距分割 
	int i, j;
	POINT* point = new POINT[n+1];
	point = partition(n);

	POINT* input = new POINT[n+1];
	double* error = new double[n+1]; 
	for (j=0; j<=n; j++){
		//计算每个插值点所对应的插值函数的值 
		input[j].x = -1 + (double)(j)/(n + 1);
		input[j].y = calc_lagrange(point, input[j].x, n);
		//计算误差
		error[j] = fabs(input[j].y - func(input[j].x));
	}
	
	double L, L1;
	L = max(error, n);
	cout<<"n:"<<n<<'\t';
	cout<<"the lagrange L& is:"<<L<<'\t';
	
	L1 = sum(error, n)/(1.0 + n);
	cout<<"the lagrange L1 is:"<<L1<<'\t';
	cout<<endl;
}

//
void lagrange2(int n){
	//等距分割 
	int i, j;
	POINT* point = new POINT[n+1];
	point = chebyshev(n);

	POINT* input = new POINT[n+1];
	double* error = new double[n+1]; 
	for (j=0; j<=n; j++){
		//计算每个插值点所对应的插值函数的值 
		input[j].x = -1 + (double)(j)/(n + 1);
		input[j].y = calc_lagrange(point, input[j].x, n);
		//计算误差
		error[j] = fabs(input[j].y - func(input[j].x));
	}
	
	double L, L1;
	L = max(error, n);
	cout<<"n:"<<n<<'\t';
	cout<<"the lagrange L& is:"<<L<<'\t';
	
	L1 = sum(error, n)/(1.0 + n);
	cout<<"the lagrange L1 is:"<<L1<<'\t';
	cout<<endl;
}
 
 
 
//主函数 
int main(){
	int n1=20, n2=40, n3=80;
	cout<<"等距分割："<<endl;
	newton1(n1);
	newton1(n2);
	newton1(n3);
	lagrange1(n1);
	lagrange1(n2);
	lagrange1(n3);
	
	cout<<endl;
	cout<<endl;
	
	cout<<"chebyshev分割："<<endl;
	newton2(n1);
	newton2(n2);
	newton2(n3);
	lagrange2(n1);
	lagrange2(n2);
	lagrange2(n3);
	return 0;
}





 
