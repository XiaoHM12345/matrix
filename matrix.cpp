#include<iostream>
#include<vector>
#include<cmath>
#include "matrix.h"

using namespace std;

matrix::matrix(int i, int j)
{
	row = i;
	column = j;
	element = new double[i * j];
	for (int k = 0; k < i * j; k++)
		element[k] = 0;
}

matrix::matrix(int i, int j, initializer_list<double> input)
{
	row = i;
	column = j;
	element = new double[i * j];
	std::initializer_list<double>::iterator it = input.begin();
	for (int k = 0; k < i * j; k++)
		*(element + k) = *(it + k);
}

matrix::matrix(const matrix& A)
{
	row = A.row;
	column = A.column;
	element = new double[row * column];
	for (int i = 0; i < row * column; i++)
		element[i] = A.element[i];
}

matrix::matrix(ifstream& matrixFile)
{
	int irow, icolumn;
	double e;
	vector<double> s(0);
	matrixFile>>irow;
	matrixFile>>icolumn;
	while(matrixFile>>e)
		s.push_back(e);
	if(irow*icolumn != s.size())
		throw dimension_mismatch();
	element = new double[irow * icolumn];
	row = irow;
	column = icolumn;
	for(int i = 0; i < row * column; i++)
		element[i] = s[i];
	matrixFile.close();
}

matrix::~matrix()
{
	delete []element;
}

matrix& matrix::operator=(const matrix& A)
{
	delete []element;
	row = A.row;
	column = A.column;
	element = new double[row * column];
	for (int i = 0; i < row * column; i++)
		element[i] = A.element[i];
	return *this;
}

matrix matrix::operator+(const matrix& A)
{
	if (A.get_column() != column || A.get_row() != row)
		throw dimension_mismatch();
	matrix temp(row, column);
	for (int i = 0; i < row * column; i++)
		temp.element[i] = element[i] + A.element[i];
	return temp;
}

matrix matrix::operator-(const matrix& A)
{
	if (A.get_column() != column || A.get_row() != row)
		throw dimension_mismatch();
	matrix temp(row, column);
	for (int i = 0; i < row * column; i++)
		temp.element[i] = element[i] - A.element[i];
	return temp;
}

matrix matrix::operator*(const matrix& A)
{
	if (column != A.get_row())
		throw dimension_mismatch();
	matrix temp(row, A.column);
	int i, j, k;
	double t;
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < A.column; j++)
		{
			t = 0;
			for (k = 0; k < column; k++)
				t += element[i * column + k] * A.element[k * column + j];
			temp.element[i * temp.column + j] = t;
		}
	}
	return temp;
}

double& matrix::operator()(int i, int j) const
{
	return element[i * column + j];
}


int matrix::get_row() const
{
	return row;
}

int matrix::get_column() const
{
	return column;
}

void matrix::print()
{
	int i, j;
	for (i = 0; i < this->get_row(); i++)
	{
		for (j = 0; j < this->get_column(); j++)
			cout << element[i * column + j] << "     ";
		cout << '\n';
	}
}

void matrix::get_from_console()
{
	int i, j;
	for (i = 0; i < row; i++)
	{
		cout << "enter the " << i + 1 << "row of the matrix:" << endl;
		for (j = 0; j < column; j++)
			cin >> element[i * column + j];
	}
}

matrix inv(const matrix& A)
{
	if (A.get_column() != A.get_row())
		throw dimension_mismatch();

	matrix temp(A.get_row(), 2 * A.get_column());
	matrix C(A.get_row(), A.get_column());
	int i, j, k;
	vector<int>local = { 0,0 }, h, l;
	double max;

	//初始化增广矩阵
	for(i=0;i<A.get_row();i++)
		for (j = 0; j < 2 * A.get_column(); j++)
		{
			if (j >= A.get_column() && i == (j - A.get_column()))
				temp(i, j) = 1;
			else if (j >= A.get_column() && i != (j - A.get_column()))
				temp(i, j) = 0;
			else
				temp(i, j) = A(i, j);
		}
	C = A;

	for (i = 0; i < A.get_row(); i++)
	{
		local = find_max_of_matrix(C, h, l);
		max = C(local[0], local[1]);
		h.push_back(local[0]);
		l.push_back(local[1]);

		for (int i = 0; i < 2 * A.get_row(); i++)
			temp(local[0], i) /= max;
		for (j = 0; j < A.get_row(); j++)
		{
			double temp2 = temp(j, local[1]);
			if (j == local[0])
				continue;
			for (k = 0; k < 2 * A.get_column(); k++)
				temp(j, k) -= temp(local[0], k) * temp2;
		}

		for (int i1 = 0; i1 < A.get_row(); i1++)
			for (int j1 = 0; j1 < A.get_column(); j1++)
				C(i1, j1) = temp(i1, j1);
	}
	//排序
	for (i = 0; i < temp.get_row(); i++)//The number of loops represents that the (i, i) element should be 1
	{
		double temp3;
		k = 0;//Indicator variable
		for (j = 0; j < temp.get_column(); j++)
		{
			if (k == 1)
				continue;
			//xchange process
			if (temp(j, i) == 1)
			{
				for (int o = 0; o < temp.get_column(); o++)
				{
					temp3 = temp(i, o);
					temp(i, o) = temp(j, o);
					temp(j, o) = temp3;
					k = 1;
				}
			}
		}
	}

	for (i = 0; i < A.get_row(); i++)
		for (j = 0; j < A.get_column(); j++)
			C(i, j) = temp(i, j + A.get_column());

	return C;
}

vector<int> find_max_of_matrix(const matrix& A, vector<int>h, vector<int>l)
{
	vector<int>temp;
	temp.resize(2);
	double max = -1;
	temp[0] = temp[1] = 0;
	for (int i = 0; i < A.get_row(); i++)
	{
		int state = 0;
		for(vector<int>::iterator it1=h.begin();it1!=h.end();it1++)
			if (i == *it1)
			{
				state = 1;
				break;
			}
		if (state == 1)
			continue;
		for (int j = 0; j < A.get_column(); j++)
		{
			state = 0;
			for(vector<int>::iterator it2=l.begin();it2!=l.end();it2++)
				if (j == *it2)
				{
					state = 1;
					break;
				}
			if (state == 1)
				continue;
			if (fabs(A(i, j)) > max)
			{
				max = fabs(A(i, j));
				temp[0] = i;
				temp[1] = j;
			}
		}
	}
	return temp;
}
