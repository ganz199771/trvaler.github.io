#pragma once
#include<cmath>
#include<iostream>
using namespace std;

// 定义二维数组赋值操作
void Assign(double** a, double** b, int n) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			a[i][j] = b[i][j];
		}
	}
}

// 定义减去一定倍数的单位矩阵的函数
double** Negate(double** a, int n, double t) {
	double** c = new double* [n];
	for (int i = 0; i < n; ++i)
		c[i] = new double[n];
	Assign(c, a, n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j)
			c[i][j] = (i == j) ? (c[i][j] - t) : c[i][j];
	}
	return c;
}

// 定义一个找副对角线最大值函数
double Find_Max(double** a, int n) {
	double max = abs(a[1][0]);
	for (int i = 1; i < n; ++i) {
		max = abs(a[i][i - 1]) > max ? abs(a[i][i - 1]) : max;
	}
	return max;
}

// 定义一个找最小绝对值的函数
double Min_num(const double& a, const double& b) {
	return abs(a) < abs(b) ? abs(a) : abs(b);
}

// 定义交换函数
void  my_swap(double& a, double& b) {
	double temp = 0.0;
	temp = a;
	a = b;
	b = temp;
}

// 定义符号函数
int Sign(double& a) {
	int out;
	out = a > 0.0 ? -1 : 1;
	return out;
}

// 定义求根函数
double* Func_Resolve(const double& a, const double& b, const double& c, const double& d) {
	double* out = new double[2];
	double t = a * d - b * c;
	out[0] = (a + d)/ 2.0;
	out[1] = sqrt(abs(pow((a + d), 2.0) - 4.0 * t)) / 2.0;
	return out;
}

// 定义矩阵相乘函数
double** Multi_M(double** a, double** b, int n) {
	double total;
	double** c = new double* [n];
	for (int i = 0; i < n; ++i)
		c[i] = new double[n];
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			total = 0.0;
			for (int s = 0; s < n; ++s) {
				total+= a[i][s]*b[s][j];
			}
			c[i][j] = total;
		}
	}
	return c;
}

// 定义矩阵和数组相乘的函数
double* Multi_V(double** m, double* source,int n) {//ok
	double total;
	double* out = new double[n];
	for (int i = 0; i < n; ++i) {
		total = 0.0;
		for (int j = 0; j < n; ++j) {
			total += m[i][j] * source[j];
		}
		out[i] = total;
	}
	return out;
}

// 定义矩阵相减函数
double** Minus(double** a, double** b,int n) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			a[i][j] -= b[i][j];
		}
	}
	return a;
}

// 定义向量相乘得到矩阵的函数
double** Vec_Multi_Vec(double* a, double* b, int n) {
	double** c = new double* [n];
	for (int i = 0; i < n; ++i)
		c[i] = new double[n];
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			c[i][j] = a[i] * b[j];
		}
	}
	return c;
}

// 定义矩阵乘数的函数
double** Multi_N(double** a, const double& b,int n) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			a[i][j] *= b;
		}
	}
	return a;
}

// 定义矩阵转置函数
double** T_Matrix(double** a, int n) {
	double** out = new double* [n];
	for (int i = 0; i < n; ++i)
		out[i] = new double[n];
	Assign(out, a, n);

	for (int i = 0; i < n - 1; ++i) {
		for (int j = i + 1; j < n; ++j) {
			my_swap(out[i][j], out[j][i]);
		}
	}
	return out;
}

// 定义一个打印函数
void Print_Matrix(double** a, int n) {
	for (int i = 0; i < n; ++i) {
		std::cout << i + 1 << "行 - " << "前半行: ";
		for (int j = 0; j < n; ++j) {
			std::cout << a[i][j] << "  ";
			if (j % (n / 2) == 4 && j != n - 1) {
				std::cout << '\n';
				std::cout << i + 1 << "行 - " << "后半行: ";
			}
		}
		std::cout << std::endl;
		std::cout << "-----------------------------------------------------------\n\n";
	}
}

// 定义拟上三角函数
double** Hessenberg(double** a, int n) {
	int check, times;
	double d, c, h, t;
	double* u = new double[n];
	double* p = new double[n];
	double* q = new double[n];
	do {
		for (int r = 0; r < n - 2; ++r) {
			check = 0;
			for (int i = r + 2; i < n; ++i) {
				if (a[i][r] == 0.0)
					++check;
			}
			if (check != (n - r - 2)) {
				d = 0.0;
				for (int i = r + 1; i < n; ++i)
					d += pow(a[i][r], 2.0);
				d = sqrt(d);
				c = Sign(a[r + 1][r]) * d;
				h = c * c - c * a[r + 1][r];

				u[r + 1] = a[r + 1][r] - c;
				for (int i = 0; i < r + 1; ++i)
					u[i] = 0.0;
				for (int i = r + 2; i < n; ++i)
					u[i] = a[i][r];

				p = Multi_V(T_Matrix(a, n), u, n);
				q = Multi_V(a, u, n);
				for (int i = 0; i < n; ++i) {
					p[i] /= h;
					q[i] /= h;
				}

				t = 0.0;
				for (int i = 0; i < n; ++i)
					t += p[i] * u[i];
				t /= h;

				for (int i = 0; i < n; ++i)
					q[i] -= t * u[i];

				Minus(a, Vec_Multi_Vec(q, u, n), n);
				Minus(a, Vec_Multi_Vec(u, p, n), n);
			}
		}
		times = 0;
		for (int r = 0; r < n - 2; ++r) {
			check = 0;
			for (int i = r + 2; i < n; ++i) {
				if (a[i][r] == 0.0)
					++check;
			}
			if (check != n - r - 2) {
				times = 1;
				break;
			}
		}
	} while (times == 1);
	return a;
}

// QR迭代求特征值
// limit是精度水平
double** QR_One_Step(double** a, int n, double limit) {
	int it = 0;
	int check;
	double d, c, h,val_1,val_2;
	val_1 = val_2 = 0.0;
	double* u = new double[n];
	double* p = new double[n];
	double* w = new double[n];
	double** Q = new double* [n];
	for (int i = 0; i < n; ++i)
		Q[i] = new double[n]; // 以上是初始化
	
	do {
		val_1 = val_2;

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j)
				Q[i][j] = (i == j) ? 1.0 : 0.0;
		}// 单位化

		for (int r = 0; r < n - 1; ++r) {
			check = 0;
			for (int i = r + 1; i < n; ++i) {
				if (a[i][r] == 0.0)
					++check;
			}
			if (check != n - r - 1) {
				d = 0.0;
				for (int i = r; i < n; ++i)
					d += pow(a[i][r], 2.0);
				d = sqrt(d);
				c = Sign(a[r][r]) * d;
				h = c * c - c * a[r][r];
				u[r] = a[r][r] - c;
				for (int i = 0; i < r; ++i)
					u[i] = 0.0;
				for (int i = r + 1; i < n; ++i)
					u[i] = a[i][r];

				w = Multi_V(Q, u, n);
				Minus(Q, Multi_N(Vec_Multi_Vec(w, u, n), (1.0 / h), n), n);
				p = Multi_V(T_Matrix(a, n), u, n);
				for (int i = 0; i < n; ++i)
					p[i] /= h;
				Minus(a, Vec_Multi_Vec(u, p, n), n);
			}
		}// 每结束一次for循环就得到相应的上三角阵和正交矩阵Q

		Assign(a, Multi_M(a, Q, n), n);// A = RQ = QAQ

		val_2 = Find_Max(a, n);
		++it;
	} while (abs(val_2 - val_1) > limit && it < 450);// 验证单步迭代450步刚好

	Hessenberg(a, n);
	// 对此时的A进行上三角化去除很小的元素
	for (int i = 1; i < n; ++i) {
		if (abs(a[i][i - 1]) < limit * Min_num(a[i][i], a[i - 1][i - 1]))
			a[i][i - 1] = 0.0;
	}
	return a;
}

// 在利用QR迭代后需要求解一些复特征根及实特征根
double** Answer(double** a, int n) {
	double** b = new double* [n];
	for (int i = 0; i < n; ++i)
		b[i] = new double[2];

	for (int i = 0; i < n - 1; ++i) {
		if (a[i + 1][i] == 0.0) {// 当副对角线元素是0，以实特征根输出
			b[i][0] = a[i][i];
			b[i][1] = 0.0;
			if (i == n - 2) { // 倒数第二行检查需要另增加条件
				b[i+1][0] = a[i+1][i+1];
				b[i+1][1] = 0.0;
				break;
			}
		}
		else {// 当副对角线元素不是0，解方程以复特征根输出
			b[i][0] = Func_Resolve(a[i][i], a[i][i + 1], a[i + 1][i], a[i + 1][i + 1])[0];
			b[i][1] = Func_Resolve(a[i][i], a[i][i + 1], a[i + 1][i], a[i + 1][i + 1])[1];
			b[i + 1][0] = Func_Resolve(a[i][i], a[i][i + 1], a[i + 1][i], a[i + 1][i + 1])[0];
			b[i + 1][1] = -1.0 * Func_Resolve(a[i][i], a[i][i + 1], a[i + 1][i], a[i + 1][i + 1])[1];
			++i;
		}
	}
	return b;
}


// 声明一个打印特征值的函数 n是特征值个数
void Print_λ(double** a, int n) {
	for (int i = 0; i < n; ++i) {
		cout << "λ" << i + 1 << ":[ " << a[i][0] << " , " << a[i][1] << " ]\n\n";
	}
}

// 添加一个求实特征值对应特征向量的方法：λ是实特征值，n是向量维数
// 由于本题矩阵的性质决定了在列主元GAUSS消去法中不会出现一列全是0的情况
double* Vector(double** a, int n) {
	int ik;
	double compare, m, temp;
	double* x = new double[n];

	for (int k = 0; k < n - 1; ++k) {
		compare = abs(a[k][k]);
		ik = k;
		for (int i = k; i < n; ++i) {	
			if (abs(a[i][k]) > compare) {
				compare=abs(a[i][k]);
				ik = i;
			}
		}

		for (int j = k; j < n; ++j)
			my_swap(a[k][j], a[ik][j]);
		
		for (int i = k + 1; i < n; ++i) {
			m = a[i][k] / a[k][k];
			for (int j = k + 1; j < n; ++j)
				a[i][j] -= m * a[k][j];
		}
	}// 以上是消元下面是回代
	x[n - 1] = 1.0;
	for (int k = n - 2; k >= 0; --k) {
		temp = 0.0;
		for (int j = k + 1; j < n; ++j)
			temp += a[k][j] * x[j];
		x[k] = -temp / a[k][k];
	}
	return x;
}

// 添加一个QR分解函数,返回Q，修改A
double** QR_Split(double** a, int n) {
	int check;
	double d, c, h;
	double* u = new double[n];
	double* p = new double[n];
	double* w = new double[n];

	double** Q = new double* [n];
	for (int i = 0; i < n; ++i)
		Q[i] = new double[n];

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j)
			Q[i][j] = (i == j) ? 1.0 : 0.0;
	}
	for (int r = 0; r < n - 1; ++r) {
		check = 0;
		for (int i = r + 1; i < n; ++i) {
			if (a[i][r] == 0.0)
				++check;
		}
		if (check != n - r - 1) {
			d = 0.0;
			for (int i = r; i < n; ++i)
			    d += pow(a[i][r], 2.0);
			d = sqrt(d);
			c = Sign(a[r][r]) * d;
			h = c * c - c * a[r][r];
			u[r] = a[r][r] - c;
			for (int i = 0; i < r; ++i)
				u[i] = 0.0;
			for (int i = r + 1; i < n; ++i)
				u[i] = a[i][r];

			w = Multi_V(Q, u, n);
			Minus(Q, Multi_N(Vec_Multi_Vec(w, u, n), (1.0 / h), n), n);
			p = Multi_V(T_Matrix(a, n), u, n);
			for (int i = 0; i < n; ++i)
				p[i] /= h;
			Minus(a, Vec_Multi_Vec(u, p, n), n);
		}
	}
	return Q;
}

// 第9步函数 哪里出错了？
void Step_9(double** a, int m) {
	double s = a[m - 2][m - 2] + a[m - 1][m - 1];
	double t = a[m - 2][m - 2] * a[m - 1][m - 1] - a[m - 1][m - 2] * a[m - 2][m - 1];

	double** MK = new double* [m];
	for (int i = 0; i < m; ++i)
		MK[i] = new double[m];

	double** Lose = new double* [m];
	for (int i = 0; i < m; ++i)
		Lose[i] = new double[m];
	Assign(Lose, a, m);

	Assign(MK, Multi_M(a, a, m), m);
	Minus(MK, Multi_N(Lose, s, m), m);
	Assign(MK, Negate(MK, m, -t), m);
	
	// 对MK做QR分解
	double** Q = new double* [m];
	for (int i = 0; i < m; ++i)
		Q[i] = new double[m];
	Q = QR_Split(MK, m);
	Assign(a, Multi_M(T_Matrix(Q, m), a, m), m);
	Assign(a, Multi_M(a, Q, m), m);
}


// 双步位移法,n是矩阵维数，limit是误差允许值，it限制迭代次数
double** TZZ(double** a, int n, double limit, int it) {
	int k = 1;// k用来检查迭代次数
	int m = n;
	int check = 0;// check检查输出条件，有了n个特征值就结束循环
	double b1, b2, b3, b4;
	double** Out = new double*[n];
	for (int i = 0; i < n; ++i)
		Out[i] = new double[2];

	do {
		while (abs(a[m - 1][m - 2]) <= limit) {
			Out[m - 1][0] = a[m - 1][m - 1];
			Out[m - 1][1] = 0;
			++check;
			--m;// 下面转第4步
			if (m == 1) {
				Out[0][0] = a[0][0];
				Out[0][1] = 0;
				++check;
				return Out;// 停止计算
			}
			else if (m == 0)
				return Out;// 停止计算
			// 否则继续循环
			else {
				// 空语句只是起一个循环头尾作用
			}
		}

		// 下面开始第5步

		b1 = Func_Resolve(a[m - 2][m - 2], a[m - 2][m - 1], a[m - 1][m - 2], a[m - 1][m - 1])[0];
		b2 = Func_Resolve(a[m - 2][m - 2], a[m - 2][m - 1], a[m - 1][m - 2], a[m - 1][m - 1])[1];
		b3 = b1;
		b4 = -1.0 * b2;

		if (m == 2) {
			Out[0][0] = b1;
			Out[0][1] = b2;
			Out[1][0] = b3;
			Out[1][1] = b4;
			check += 2;
			return Out;
		}
		else {// 转第7步
			if (abs(a[m - 2][m - 3]) <= limit) {
				Out[m - 2][0] = b1;
				Out[m - 2][1] = b2;
				Out[m - 1][0] = b3;
				Out[m - 1][1] = b4;
				check += 2;
				m -= 2;
				// 转第4步
				if (m == 1) {
					Out[0][0] = a[0][0];
					Out[0][1] = 0;
					++check;
					return Out;// 停止计算
				}
				else if (m == 0)
					return Out;// 停止计算
				else {
					// 下面进入前面写好的while循环
				}
			}
			else {// 转第8步
				if (k == it)
					return Out;
				else {// 转第9步
					double** Temp_V = new double*[m];
					for (int i = 0; i < m; ++i)
						Temp_V[i] = new double[m];// 获得一个低阶矩阵
					Assign(Temp_V,a,m);// 初始化低阶矩阵
					// 执行第9步函数
					Step_9(Temp_V, m);
					// 将修改的m*m矩阵元素赋值到原始矩阵n*n
					Assign(a, Temp_V, m);
					++k;// 转第3步也就是开始的while循环
				}
			}
		}
	} while (check < n);
	return Out;
}