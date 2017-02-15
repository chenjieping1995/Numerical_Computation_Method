/*
 * @ filename: jacobi.cpp
 * @ author: 陈介平
 * @ student number: PB14209115
 * @ date: 2016.12.09
 * @ function: 实现 jacobi迭代功能 
*/
#include<stdio.h>
#include<math.h>
#define N 9

double Compare(double a[N],double b[N]){
	double c, max=0;
	int i;
	for(i=0;i<=N-1;i++){
		c = fabs(a[i]-b[i]);
		if (c>max)
			max = c;
	}
	return max;
}

void Jacobi(double A[N][N], double x[N], double b[N], double precesion){
	int i,j,k;
	double x2[N],sum;
	for(i=0;i<=N-1;i++)
		x2[i]=x[i];		//将初始迭代向量赋给x2
	k=1;	//k为迭代次数
	while(1){
		for(i=0;i<=N-1;i++){
			sum=0;
			for(j=0;j<=N-1;j++){
				if(j!=i)
					sum+=A[i][j]*x2[j];
			}
			x[i]=(b[i]-sum)/A[i][i];	//以x2为基础进行迭代求出x
		}	
		//输出每一次迭代的结果
		printf("第%d次迭代:\n",k);
		printf("x%d= ",k);
		for(i=0;i<=N-1;i++)
			printf("%lf ",x2[i]);
		printf("\n");
		printf("x%d= ",k+1);
		for(i=0;i<=N-1;i++)
			printf("%lf ",x[i]);
		printf("\n");	
		//判断是否达到迭代精度
		if(Compare(x2,x)<precesion){
			printf("达到迭代精度的方程组的解为:\n");
			for(i=0;i<=N-1;i++){
				printf("x[%d]= ",i+1);
				printf("%lf \n",x[i]);
			}
			printf("\n");
			break;
		}
		else{
			for(i=0;i<=N-1;i++)
				x2[i]=x[i];	//将第k次迭代计算得到的向量x赋给x2
			k++;
			continue;
		}
	}
}

int main(){
	double A[N][N]={{31,-13,0,0,0,-10,0,0,0},
					{-13,35,-9,0,-11,0,0,0,0},
					{0,-9,31,-10,0,0,0,0,0},
					{0,0,-10,79,-30,0,0,0,-9},
					{0,0,0,-30,57,-7,0,-5,0},
					{0,0,0,0,-7,47,-30,0,0},
					{0,0,0,0,0,-30,41,0,0},
					{0,0,0,0,-5,0,0,27,-2},
					{0,0,0,-9,0,0,0,-2,29}},
			b[N]={-15,27,-23,0,-20,12,-7,7,10};
	double x[N]={0};
	Jacobi(A,x,b,1e-6);
	return 0;
}
