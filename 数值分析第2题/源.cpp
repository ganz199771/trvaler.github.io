#include"help.h"

int main() {
	int num = 0;// ��ʵ����ֵ��
	double ��; // ʵ����ֵ
	int r[10];

	double** a = new double* [10];// ԭ����
	for (int i = 0; i < 10; ++i)
		a[i] = new double[10];

	double** A = new double* [10];// ԭ����
	for (int i = 0; i < 10; ++i)
		A[i] = new double[10];

	for (int i = 0; i < 10; ++i) {
		for (int j = 0; j < 10; ++j) {
			a[i][j] = (i != j) ? sin(0.5 * (i + 1.0) + 0.2 * (j + 1.0)) : 1.52 * cos(i + 1.0 + 1.2 * (j + 1.0));
		}
	}// ��ֵ

	double** b = new double* [10];
	for (int i = 0; i < 10; ++i)
		b[i] = new double[2];//�洢����ֵ

	double** B = new double* [10];
	for (int i = 0; i < 10; ++i)
		B[i] = new double[2];

	cout << scientific;
	cout.precision(12);

	cout << "*********************** ԭ���� ************************\n\n\n";
	Print_Matrix(a, 10);

	cout << "********************* �����ǻ����� ********************\n\n\n";
	Hessenberg(a, 10);
	Print_Matrix(a, 10);

	// ����˫��λ�Ʒ�
	Assign(A, a, 10);
	for (int i = 0; i < 2000; ++i)
		Step_9(A, 10);
	cout << "********************* ���� ********************\n\n\n";
	Print_Matrix(A, 10);
	Assign(A, a, 10);

	cout << "********************* QR��������� ********************\n\n\n";
	QR_One_Step(a, 10, pow(10, -12));
	Print_Matrix(a, 10);

	cout << "********************* ��������ֵ **********************\n\n\n";
	b = Answer(a, 10);
	Print_��(b, 10);

	cout << "********************* ˫��λ�ƾ��� ********************\n\n\n";
	B = TZZ(A, 10, pow(10, -12), 15);
	for (int i = 0; i < 10; ++i) {
		for (int j = 0; j < 10; ++j)
			A[i][j] = (abs(A[i][j]) > pow(10, -12)) ? A[i][j] : 0.0;
	}// ����СԪ��
	Print_Matrix(A, 10);
	cout << "******************* ˫��λ�ƾ�������ֵ ****************\n\n\n";
	Print_��(B, 10);

	for (int i = 0; i < 10; ++i) {
		for (int j = 0; j < 10; ++j) {
			a[i][j] = (i != j) ? sin(0.5 * (i + 1.0) + 0.2 * (j + 1.0)) : 1.52 * cos(i + 1.0 + 1.2 * (j + 1.0));
		}
	}// �ٴν��и�ֵ���������������
	for (int i = 0; i < 10; ++i) {
		if (b[i][1] == 0.0) {//����ʵ����ֵ���
			++num;
		}//ȷ�����ٸ�ʵ����ֵ
	}

	double** y = new double* [num];//�洢��������
	for (int i = 0; i < 10; ++i)
		y[i] = new double[10];

	int j = 0;
	for (int i = 0; i < 10; ++i) {
		if (b[i][1] == 0.0) {
			r[j] = i + 1;
			�� = b[i][0];
			y[j] = Vector(Negate(a, 10, ��), 10);
			++j;
		}// �õ��������������棬Ҳ����ʹ��vector������ʵ�ֲ�β��
		//  �Ͳ���֪���ж��ٸ�ʵ����ֵ�������ڴ�
	}

	cout << "********************* ������������ **********************\n\n\n";

	for (int i = 0; i < num; ++i) {
		cout << "��" << r[i];
		cout.width(22);
	}
	cout << endl;
	for (int i = 0; i < 10; ++i) {
		for (int j = 0; j < num; ++j) {
			cout << y[j][i] << "   ";
		}
		cout << endl << endl;
	}

	// �ڴ��ͷ�
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