/*
 * @ filename: �Ľ�Runge-Kutta��.cpp
 * @ author: �½�ƽ
 * @ student number: PB14209115
 * @ date: 2016.12.30
 * @ function: ���Ľ�Runge-Kutta����ⳣ΢�ַ��̳�ֵ���� 
 */
 
#include <stdio.h>
#include <math.h>
#include <iostream>
#define PI 3.14159265358979323846 
#define e 2.718281828

using namespace std;

double f(double x, double y){
	// f(x, y) = -50y + 50x^2 + 2x
	double f;
	f = (-50)*y + 50*pow(x,2) + 2*x;
	return f;
}

double fx(double x){
	double y;
	y = 1.0/3.0*pow(e,(-50)*x) + pow(x,2);
	return y;
}

void Runge_Kutta(double h){
	// �ⷽ�������庯��
	double x_end;
	x_end = 1;
	double x0, y0;
	double xn, yn, yn1; 
	x0 = 0;
	y0 = double(1.0/3.0);
	
	xn = x0;
	yn = y0;
	double k1, k2, k3, k4;
	int n;
	n = int(x_end/h);
	cout<<"��h="<<h<<"ʱ��n="<<n<<"����΢�ַ��̳�ֵ����Ľ�Ϊ��"<<endl;
	for (int i=1; i<=n; i++){
		k1 = f(xn, yn);
		k2 = f((xn + h/2.0), (yn + h*k1/2.0));
		k3 = f((xn + h/2.0), (yn + h*k2/2.0));
		k4 = f((xn + h), (yn + h*k3));
		yn1 = yn + h*(k1+ 2.0*k2 + 2.0*k3 +k4)/6.0;
		xn += h;
		cout<<"x"<<i<<"="<<xn<<",\t"<<"y"<<i<<"="<<yn1<<"��\t";
		double wucha = yn1 - fx(xn);
		cout<<"��ʱ�����Ϊ:"<<wucha<<endl;
		yn = yn1;
	}
	cout<<endl<<endl;
} 

int main(){
	// ������������Runge_Kutta�����������Ľ�Runge-Kutta���ⳣ΢�ַ��̳�ֵ����
	Runge_Kutta(0.1);
	Runge_Kutta(0.025);
	Runge_Kutta(0.01);
}
