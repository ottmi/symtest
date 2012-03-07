#include "Matrix.h"
#include <iostream>
#include <iomanip>

Matrix::Matrix(unsigned int dim)
{
	_dim = dim;
	vector<double> v = vector<double>(_dim, .0);
	_m = vector<vector<double> >(_dim, v);
}

Matrix::~Matrix()
{
	// TODO Auto-generated destructor stub
}

double& Matrix::operator() (unsigned row, unsigned col)
{
	if (row >= _dim || col >= _dim)
		throw (string("Matrix subscript out of bounds"));
	return _m[row][col];
}

double Matrix::operator() (unsigned row, unsigned col) const
{
	if (row >= _dim || col >= _dim)
		throw (string("Matrix subscript out of bounds"));
	return _m[row][col];
}

void Matrix::zero()
{
	for (unsigned int i = 0; i < _dim; i++)
		for (unsigned int j = 0; j < _dim; j++)
			_m[i][j] = .0;
}

void Matrix::set(vector<vector<double> > x)
{
	for (unsigned int i = 0; i < _dim; i++)
		for (unsigned int j = 0; j < _dim; j++)
			_m[i][j] = x[i][j];
}

void Matrix::set(vector<double> x)
{
	for (unsigned int i = 0; i < _dim; i++)
		for (unsigned int j = 0; j < _dim; j++)
			_m[i][j] = x[i * _dim + j];
}

void Matrix::setDiag(double x)
{
	for (unsigned int i = 0; i < _dim; i++)
		_m[i][i] = x;
}

void Matrix::setOffDiag(double x)
{
	for (unsigned int i = 0; i < _dim; i++)
		for (unsigned int j = 0; j < _dim; j++)
			if (i != j) _m[i][j] = x;
}

void Matrix::setRow(unsigned int row, vector<double> x)
{
	for (unsigned int j = 0; j < _dim; j++)
		_m[row][j] = x[j];
}

void Matrix::setCol(unsigned int col, vector<double> x)
{
	for (unsigned int i = 0; i < _dim; i++)
		_m[i][col] = x[i];
}

void Matrix::update(Matrix& x)
{
	double sum = 0;
	for (unsigned int i = 0; i < _dim; i++)
		for (unsigned int j = 0; j < _dim; j++)
			sum += _m[i][j];

	for (unsigned int i = 0; i < _dim; i++)
		for (unsigned int j = 0; j < _dim; j++)
			_m[i][j] = _m[i][j] / sum;
}

double Matrix::getRowSum(unsigned int row)
{
	double sum = .0;
	for (unsigned int i = 0; i < _dim; i++)
		sum += _m[row][i];
	return sum;
}

double Matrix::getColSum(unsigned int col)
{
	double sum = .0;
	for (unsigned int i = 0; i < _dim; i++)
		sum += _m[i][col];
	return sum;
}

double Matrix::determinant()
{
	double det = 0;
	if (_dim == 2)
	{
		det = _m[0][0] * _m[1][1] - _m[0][1] * _m[1][0];
	} else if (_dim == 3)
	{
		det = _m[0][0] * _m[1][1] * _m[2][2] + _m[0][1] * _m[1][2] * _m[2][0] + _m[0][2] * _m[1][0] * _m[2][1] - _m[0][0] * _m[1][2] * _m[2][1]
				- _m[0][1] * _m[1][0] * _m[2][2] - _m[0][2] * _m[1][1] * _m[2][0];
	} else if (_dim == 4)
	{
		det = _m[0][3] * _m[1][2] * _m[2][1] * _m[3][0] - _m[0][2] * _m[1][3] * _m[2][1] * _m[3][0] - _m[0][3] * _m[1][1] * _m[2][2] * _m[3][0]
				+ _m[0][1] * _m[1][3] * _m[2][2] * _m[3][0] + _m[0][2] * _m[1][1] * _m[2][3] * _m[3][0] - _m[0][1] * _m[1][2] * _m[2][3] * _m[3][0]
				- _m[0][3] * _m[1][2] * _m[2][0] * _m[3][1] + _m[0][2] * _m[1][3] * _m[2][0] * _m[3][1] + _m[0][3] * _m[1][0] * _m[2][2] * _m[3][1]
				- _m[0][0] * _m[1][3] * _m[2][2] * _m[3][1] - _m[0][2] * _m[1][0] * _m[2][3] * _m[3][1] + _m[0][0] * _m[1][2] * _m[2][3] * _m[3][1]
				+ _m[0][3] * _m[1][1] * _m[2][0] * _m[3][2] - _m[0][1] * _m[1][3] * _m[2][0] * _m[3][2] - _m[0][3] * _m[1][0] * _m[2][1] * _m[3][2]
				+ _m[0][0] * _m[1][3] * _m[2][1] * _m[3][2] + _m[0][1] * _m[1][0] * _m[2][3] * _m[3][2] - _m[0][0] * _m[1][1] * _m[2][3] * _m[3][2]
				- _m[0][2] * _m[1][1] * _m[2][0] * _m[3][3] + _m[0][1] * _m[1][2] * _m[2][0] * _m[3][3] + _m[0][2] * _m[1][0] * _m[2][1] * _m[3][3]
				- _m[0][0] * _m[1][2] * _m[2][1] * _m[3][3] - _m[0][1] * _m[1][0] * _m[2][2] * _m[3][3] + _m[0][0] * _m[1][1] * _m[2][2] * _m[3][3];
	} else
	{
		/* Perform Crout's LU decomposition, the determinant will be the product
		 * of the diagonal of the lower triangular matrix L */
		det = 1;
		vector<vector<double> > lu(_dim, vector<double>(_dim));
		for (unsigned int i = 0; i < _dim; i++)
		{
			for (unsigned int j = i; j < _dim; j++)
			{
				double sum = 0.;
				for (unsigned int k = 0; k < i; k++)
					sum += lu[j][k] * lu[k][i];
				lu[j][i] = _m[j][i] - sum;
			}

			if (lu[i][i] == 0) return .0;

			for (unsigned int j = i + 1; j < _dim; j++)
			{
				double sum = 0.;
				for (unsigned int k = 0; k < i; k++)
					sum += lu[i][k] * lu[k][j];
				lu[i][j] = (_m[i][j] - sum) / lu[i][i];
			}

			det *= lu[i][i];
		}
	}

	return det;
}

void Matrix::print()
{
	cout.precision(8);
	for (unsigned int i = 0; i < _dim; i++)
	{
		for (unsigned int j = 0; j < _dim; j++)
			cout << fixed << _m[i][j] << " ";
		cout << endl;
	}
}
