// 通过编程实现数值分析中所需要的算法实例
#include<iostream>
using namespace std;

// 第一个小例子：一维多项式求值
// f(x)=A0+A1*x+A2*x^2+...+An-1*x^n-1
double YWDXS(double a[], int n, double x)
{
	double temp = a[n - 1];
	for (int i = n - 2; i >= 0; --i)
		temp = temp * x + a[i];
	return temp;
}

int main()
{
	cout << "******* 一维多项式 *********\n";
	double a[7] = { -20.0,7.0,-7.0,1.0,3.0,-5.0,2.0 };
	double x[6] = { 0.9,-0.9,1.1,-1.1,1.3,-1.3 };

	for (int i = 0; i < 6; i++)
		cout << "x(" << i <<") = "<< x[i] <<'\t'<< "f(" << i << ") = " << YWDXS(a, 7, x[i]) << endl;
	cout << "******* 一维多项式算例完成 *********\n";



}