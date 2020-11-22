#include"help.h"

int main() {
	int num = 0;// 计实特征值数
	double λ; // 实特征值
	int r[10];

	double** a = new double* [10];// 原矩阵
	for (int i = 0; i < 10; ++i)
		a[i] = new double[10];

	double** A = new double* [10];// 原矩阵
	for (int i = 0; i < 10; ++i)
		A[i] = new double[10];

	for (int i = 0; i < 10; ++i) {
		for (int j = 0; j < 10; ++j) {
			a[i][j] = (i != j) ? sin(0.5 * (i + 1.0) + 0.2 * (j + 1.0)) : 1.52 * cos(i + 1.0 + 1.2 * (j + 1.0));
		}
	}// 赋值

	double** b = new double* [10];
	for (int i = 0; i < 10; ++i)
		b[i] = new double[2];//存储特征值

	double** B = new double* [10];
	for (int i = 0; i < 10; ++i)
		B[i] = new double[2];

	cout << scientific;
	cout.precision(12);

	cout << "*********************** 原矩阵 ************************\n\n\n";
	Print_Matrix(a, 10);

	cout << "********************* 上三角化矩阵 ********************\n\n\n";
	Hessenberg(a, 10);
	Print_Matrix(a, 10);

	// 测试双步位移法
	Assign(A, a, 10);
	for (int i = 0; i < 2000; ++i)
		Step_9(A, 10);
	cout << "********************* 测试 ********************\n\n\n";
	Print_Matrix(A, 10);
	Assign(A, a, 10);

	cout << "********************* QR迭代后矩阵 ********************\n\n\n";
	QR_One_Step(a, 10, pow(10, -12));
	Print_Matrix(a, 10);

	cout << "********************* 矩阵特征值 **********************\n\n\n";
	b = Answer(a, 10);
	Print_λ(b, 10);

	cout << "********************* 双步位移矩阵 ********************\n\n\n";
	B = TZZ(A, 10, pow(10, -12), 15);
	for (int i = 0; i < 10; ++i) {
		for (int j = 0; j < 10; ++j)
			A[i][j] = (abs(A[i][j]) > pow(10, -12)) ? A[i][j] : 0.0;
	}// 消除小元素
	Print_Matrix(A, 10);
	cout << "******************* 双步位移矩阵特征值 ****************\n\n\n";
	Print_λ(B, 10);

	for (int i = 0; i < 10; ++i) {
		for (int j = 0; j < 10; ++j) {
			a[i][j] = (i != j) ? sin(0.5 * (i + 1.0) + 0.2 * (j + 1.0)) : 1.52 * cos(i + 1.0 + 1.2 * (j + 1.0));
		}
	}// 再次进行赋值用于求解特征向量
	for (int i = 0; i < 10; ++i) {
		if (b[i][1] == 0.0) {//对于实特征值情况
			++num;
		}//确定多少个实特征值
	}

	double** y = new double* [num];//存储特征向量
	for (int i = 0; i < 10; ++i)
		y[i] = new double[10];

	int j = 0;
	for (int i = 0; i < 10; ++i) {
		if (b[i][1] == 0.0) {
			r[j] = i + 1;
			λ = b[i][0];
			y[j] = Vector(Negate(a, 10, λ), 10);
			++j;
		}// 得到特征向量并保存，也可以使用vector等容器实现插尾，
		//  就不用知道有多少个实特征值再申请内存
	}

	cout << "********************* 矩阵特征向量 **********************\n\n\n";

	for (int i = 0; i < num; ++i) {
		cout << "ξ" << r[i];
		cout.width(22);
	}
	cout << endl;
	for (int i = 0; i < 10; ++i) {
		for (int j = 0; j < num; ++j) {
			cout << y[j][i] << "   ";
		}
		cout << endl << endl;
	}

	// 内存释放
	for (int i = 0; i < 10; ++i) {
		if (i < num) {
			delete[] y[i];
			delete[] b[i];
			delete[] a[i];
			delete[] B[i];
			delete[] A[i];
		}
		else {
			delete[] b[i];
			delete[] a[i];
			delete[] B[i];
			delete[] A[i];
		}
	}

	return 0;
}