/*
 * @ filename: 四阶Runge-Kutta法.cpp
 * @ author: 陈介平
 * @ student number: PB14209115
 * @ date: 2016.12.30
 * @ function: 用四阶Runge-Kutta法求解常微分方程初值问题 
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
	// 解方程组主体函数
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
	cout<<"当h="<<h<<"时，n="<<n<<"，常微分方程初值问题的解为："<<endl;
	for (int i=1; i<=n; i++){
		k1 = f(xn, yn);
		k2 = f((xn + h/2.0), (yn + h*k1/2.0));
		k3 = f((xn + h/2.0), (yn + h*k2/2.0));
		k4 = f((xn + h), (yn + h*k3));
		yn1 = yn + h*(k1+ 2.0*k2 + 2.0*k3 +k4)/6.0;
		xn += h;
		cout<<"x"<<i<<"="<<xn<<",\t"<<"y"<<i<<"="<<yn1<<"，\t";
		double wucha = yn1 - fx(xn);
		cout<<"此时的误差为:"<<wucha<<endl;
		yn = yn1;
	}
	cout<<endl<<endl;
} 

int main(){
	// 主函数，调用Runge_Kutta函数，即用四阶Runge-Kutta法解常微分方程初值问题
	Runge_Kutta(0.1);
	Runge_Kutta(0.025);
	Runge_Kutta(0.01);
}
