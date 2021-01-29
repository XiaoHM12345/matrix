//
// Created by 84000 on 2021/1/28.
//

#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H
#include <initializer_list>
#include <exception>
#include <string>
#include <vector>
#include <mutex>
#include <iostream>

template<typename T> class Matrix;

class myException : public std::exception {
public:
    myException(std::string msg)  : msg_(msg){ }
    const char* what() const throw () {
        return msg_.c_str();
    }
private:
    std::string msg_;
};

Matrix<double> inv(const Matrix<double> &m);
std::vector<int> find_max_of_matrix(const Matrix<double>& m, std::vector<int> h, std::vector<int> l);

template<typename T>
class Matrix {
public:
    Matrix() : column_(0), row_(0), element_(0) {}
    Matrix(size_t row, size_t column)
    : column_(column),
        row_(row),
        element_ (row * column)
    { }
    Matrix(size_t row, size_t column, std::initializer_list<T> in)
            : column_(column), row_(row), element_(in)
    {
        if (row * column != in.size())
            throw myException("error : matrix size != initializer_list.size()");
    }

    //default copy ctor and operator= is ok;

    size_t size_of_row() const { return row_; }
    size_t size_of_column() const { return column_; }

    void print()  //not 100% thread safe!
    {
        for (int i = 0; i < row_; i++)
        {
            for (int j = 0; j < column_; j++)
                std::cout << element_[i * column_ + j] << " ";
            std::cout << "\n";  //not 100% thread safe!
        }
    }

    Matrix<T> dot_multiplication (const Matrix<T>& rhs)
    {
        if (column_ != rhs.column_ || row_ != rhs.row_)
            throw myException("error : in operator : row and column dont match");

        Matrix<T> ret(row_, column_);
        for (int i = 0; i < row_ * column_; i++)
            ret.element_[i] *= rhs.element_[i];
        return ret;
    }

    Matrix<T>& operator += (const Matrix<T>& rhs)
    {
        if (column_ != rhs.column_ || row_ != rhs.row_)
            throw myException("error : in operator : row and column dont match");
        for (int i = 0; i < row_ * column_; i++)
            element_[i] += rhs.element_[i];
        return *this;
    }

    Matrix<T>& operator -= (const Matrix<T>& rhs)
    {
        if (column_ != rhs.column_ || row_ != rhs.row_)
            throw myException("error : in operator : row and column dont match");
        for (int i = 0; i < row_ * column_; i++)
            element_[i] -= rhs.element_[i];
        return *this;
    }

    Matrix<T>& operator *= (const Matrix<T>& rhs)
    {
        *this = *this * rhs;
        return *this;
    }

    Matrix<T> operator * (const Matrix<T>& rhs) const
    {
        if (row_ != rhs.size_of_column())
            throw myException("error : row and column don't match");
        Matrix<T> ret(row_, rhs.column_);
        int i, j, k;
        double t;
        for (i = 0; i < row_; ++i) {
            for (j = 0; j < rhs.column_; ++j) {
                t = 0;
                for (k = 0; k < column_; ++k)
                    t += element_[i * column_ + k] * rhs.element_[k * column_ + j];
                ret.element_[i * ret.column_ + j] = t;
            }
        }
        return ret;
    }

    Matrix<T> operator - (const Matrix<T>& rhs) const
    {
        Matrix<T> ret(*this);
        ret -= rhs;

        return ret;
    }

    Matrix<T> operator + (/* Matrix<T> *this */const Matrix<T>& rhs) const
    {
        Matrix<T> ret(*this);
        ret += rhs;

        return ret;
    }

    T& operator () (size_t row, size_t column) const
    {
        if (row >= row_ || row < 0 || column >= column_ || column < 0)
            throw myException("error : in operator = : row must >= 0 and < row_ ,and so is column");
        return element_[row * column_ + column];
    }

    static Matrix<double> minus_one (size_t row, size_t column)
    {
        Matrix<double> ret(row, column);
        for (auto &v : ret.element_)
            v = -1;
        return ret;
    }

    static Matrix<double> zero (size_t row, size_t column)
    {
        Matrix<double> ret(row, column);
        for (auto &v : ret.element_)
            v = 0;
        return ret;
    }
    static Matrix<double>& negate(Matrix<double> &m)
    {
        for (auto &v : m.element_)
            v *= -1;
        return m;
    }
private:
    size_t column_;
    size_t row_;
    mutable std::vector<T> element_;
};

std::vector<int> find_max_of_matrix(const Matrix<double>& m, std::vector<int> h, std::vector<int> l)
{
    std::vector<int>temp;
    temp.resize(2);
    double max = -1;
    temp[0] = temp[1] = 0;
    for (int i = 0; i < m.size_of_row(); i++)
    {
        int state = 0;
        for(std::vector<int>::iterator it1=h.begin();it1!=h.end();it1++)
            if (i == *it1)
            {
                state = 1;
                break;
            }
        if (state == 1)
            continue;
        for (int j = 0; j < m.size_of_column(); j++)
        {
            state = 0;
            for(std::vector<int>::iterator it2=l.begin();it2!=l.end();it2++)
                if (j == *it2)
                {
                    state = 1;
                    break;
                }
            if (state == 1)
                continue;
            if (fabs(m(i, j)) > max)
            {
                max = fabs(m(i, j));
                temp[0] = i;
                temp[1] = j;
            }
        }
    }
    return temp;
}

Matrix<double> inv(const Matrix<double> &m)
{
    if (m.size_of_column() != m.size_of_row())
        throw myException("error : row and column don't match");

    Matrix<double> temp(m.size_of_row(), 2 * m.size_of_column());
    Matrix<double> C(m.size_of_row(), m.size_of_column());
    int i, j, k;
    std::vector<int>local = { 0,0 }, h, l;
    double max;

    //初始化增广矩阵
    for(i=0;i<m.size_of_row();i++)
        for (j = 0; j < 2 * m.size_of_column(); j++)
        {
            if (j >= m.size_of_column() && i == (j - m.size_of_column()))
                temp(i, j) = 1;
            else if (j >= m.size_of_column() && i != (j - m.size_of_column()))
                temp(i, j) = 0;
            else
                temp(i, j) = m(i, j);
        }
    C = m;

    for (i = 0; i < m.size_of_row(); i++)
    {
        local = find_max_of_matrix(C, h, l);
        max = C(local[0], local[1]);
        h.push_back(local[0]);
        l.push_back(local[1]);

        for (int s = 0; s < 2 * m.size_of_row(); s++)
            temp(local[0], s) /= max;
        for (j = 0; j < m.size_of_row(); j++)
        {
            double temp2 = temp(j, local[1]);
            if (j == local[0])
                continue;
            for (k = 0; k < 2 * m.size_of_column(); k++)
                temp(j, k) -= temp(local[0], k) * temp2;
        }

        for (int i1 = 0; i1 < m.size_of_row(); i1++)
            for (int j1 = 0; j1 < m.size_of_column(); j1++)
                C(i1, j1) = temp(i1, j1);
    }
    //排序
    for (i = 0; i < temp.size_of_row(); i++)//The number of loops represents that the (i, i) element should be 1
    {
        double temp3;
        k = 0;//Indicator variable
        for (j = 0; j < temp.size_of_column(); j++)
        {
            if (k == 1)
                continue;
            //xchange process
            if (temp(j, i) == 1)
            {
                for (int o = 0; o < temp.size_of_column(); o++)
                {
                    temp3 = temp(i, o);
                    temp(i, o) = temp(j, o);
                    temp(j, o) = temp3;
                    k = 1;
                }
            }
        }
    }

    for (i = 0; i < m.size_of_row(); i++)
        for (j = 0; j < m.size_of_column(); j++)
            C(i, j) = temp(i, j + m.size_of_column());

    return C;
}

#endif //MATRIX_MATRIX_H
