#pragma once
#include<initializer_list>
#include<string>
#include<vector>

using namespace std;

class matrix
{
public:
	//构造与析构
	matrix(int i, int j);
	matrix(int i, int j, initializer_list<double> input);
	matrix(const matrix& A);
	~matrix();
	matrix& operator =(const matrix& A);
	//重载运算符
	matrix operator +(const matrix& A);
	matrix operator -(const matrix& A);
	matrix operator *(const matrix& A);
	double& operator ()(int i, int j) const;
	//获取矩阵维度
	int get_row() const;
	int get_column() const;

	void print();
	void get_from_console();
private:
	int row, column;
	double* element;
};

matrix inv(const matrix& A);
vector<int> find_max_of_matrix(const matrix& A, vector<int>h, vector<int>l);