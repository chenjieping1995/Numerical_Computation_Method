/*
* @ Filename: ����ʽ��ֵ.cpp 
* @ Author: Jieping Chen
* @ Student Number: PB14209115
* @ Version: 1.0.0
* @ Date: 2016.10.07
* @ Function: �������ն���ʽ��ֵ����ţ�ٶ���ʽ��ֵ���ڵȾ�ָ���б�ѩ���ָ��µĶԱ� 
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

//���庯��
double func(double x){
	double f = 1.0/(1.0 + pow(x,2));
	return f;
} 


//�Ⱦ�ָ�
POINT* partition(int n){
	int i;
	POINT* point = new POINT[n+1];
	for (i=0; i<=n; i++){
		point[i].x = -1.0 + 2.0/n*i;
		point[i].y = func(point[i].x);
	}
	return point;
} 


//chebyshev�ָ� 
POINT* chebyshev(int n){
	int i;
	POINT* point = new POINT[n+1];
	for (i=0; i<=n; i++){
		point[i].x = cos( ((2.0*i+1.0)*PI) / double(2.0*n+1.0));
		point[i].y = func(point[i].x);
	}
	return point;
} 


//��ȡ���ֵ����
double max(double* e, int n){
	int i;
	double x=0;
	for (i=0; i<=n; i++){
		if (x<e[i]) x=e[i];
	}
	return x;
} 


//��ͺ���
double sum(double *e, int n){
	int i;
	double x=0;
	for (i=0; i<=n; i++){
		x += e[i];
	}
	return x;
} 


//����ֵ��������ʽ����
double calc_func(POINT* point, double* diff, double x, int n ){
	//����newton(x)��ֵ
	double tmp = 1;
	double newton = diff[0];
	int i;
	for (i=0; i<n; i++){
		tmp = tmp*(x - point[i].x);
		newton = newton + tmp*diff[i+1];
	}
	return newton;
}



//newton��ֵ����ʽ�йؼ��� 
void newton1(int n){
	//����newton��ֵ����ʽ 
	int i, j;
	POINT* point = new POINT[n+1];
	double* diff = new double[n+1];
	point = partition(n);
	for (i=0; i<=n; i++){
		diff[i] = point[i].y;
	}
	
	for (i=0; i<n; i++){
		for (j=n; j>i; j--){
			//����f(x0,...,xn)�Ĳ��� 
			diff[j] = (diff[j] - diff[j-1])/(point[j].x - point[j-1-i].x);
		}
	}
	
	//�������
	POINT* input = new POINT[n+1];
	double* error = new double[n+1]; 
	for (j=0; j<=n; j++){
		//����ÿ����ֵ������Ӧ�Ĳ�ֵ������ֵ 
		input[j].x = -1 + (double)(j)/(n + 1);
		input[j].y = calc_func(point, diff, input[j].x, n);
		//�������
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
	//����newton��ֵ����ʽ 
	int i, j;
	POINT* point = new POINT[n+1];
	double* diff = new double[n+1];
	point = chebyshev(n);
	for (i=0; i<=n; i++){
		diff[i] = point[i].y;
	}
	
	for (i=0; i<n; i++){
		for (j=n; j>i; j--){
			//����f(x0,...,xn)�Ĳ��� 
			diff[j] = (diff[j] - diff[j-1])/(point[j].x - point[j-1-i].x);
		}
	}
	
	//�������
	POINT* input = new POINT[n+1];
	double* error = new double[n+1]; 
	for (j=0; j<=n; j++){
		//����ÿ����ֵ������Ӧ�Ĳ�ֵ������ֵ 
		input[j].x = -1 + (double)(j)/(n + 1);
		input[j].y = calc_func(point, diff, input[j].x, n);
		//�������
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


//����lagrange��ֵ����ֵ 
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


//lagrange��ֵ����ʽ
void lagrange1(int n){
	//�Ⱦ�ָ� 
	int i, j;
	POINT* point = new POINT[n+1];
	point = partition(n);

	POINT* input = new POINT[n+1];
	double* error = new double[n+1]; 
	for (j=0; j<=n; j++){
		//����ÿ����ֵ������Ӧ�Ĳ�ֵ������ֵ 
		input[j].x = -1 + (double)(j)/(n + 1);
		input[j].y = calc_lagrange(point, input[j].x, n);
		//�������
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
	//�Ⱦ�ָ� 
	int i, j;
	POINT* point = new POINT[n+1];
	point = chebyshev(n);

	POINT* input = new POINT[n+1];
	double* error = new double[n+1]; 
	for (j=0; j<=n; j++){
		//����ÿ����ֵ������Ӧ�Ĳ�ֵ������ֵ 
		input[j].x = -1 + (double)(j)/(n + 1);
		input[j].y = calc_lagrange(point, input[j].x, n);
		//�������
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
 
 
 
//������ 
int main(){
	int n1=20, n2=40, n3=80;
	cout<<"�Ⱦ�ָ"<<endl;
	newton1(n1);
	newton1(n2);
	newton1(n3);
	lagrange1(n1);
	lagrange1(n2);
	lagrange1(n3);
	
	cout<<endl;
	cout<<endl;
	
	cout<<"chebyshev�ָ"<<endl;
	newton2(n1);
	newton2(n2);
	newton2(n3);
	lagrange2(n1);
	lagrange2(n2);
	lagrange2(n3);
	return 0;
}





 
