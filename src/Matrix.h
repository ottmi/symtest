#ifndef MATRIX_H_
#define MATRIX_H_

#include <vector>
using namespace std;

class Matrix
{
public:
	Matrix(unsigned int dim);
	virtual ~Matrix();
	double& operator() (unsigned row, unsigned col);
	double  operator() (unsigned row, unsigned col) const;
	Matrix operator*(Matrix const &m) const;

	void zero();
	void identity();

	void set(vector< vector<double> > const &m);
	void set(vector<double> const &m);
	void setDiag(double const x);
	void setDiag(vector<double> const &x);
	void setOffDiag(double const x);
	void setRow(unsigned int const row, vector<double> const &m);
	void setCol(unsigned int const col, vector<double> const &m);

	double getRowSum(unsigned int const row) const;
	vector<double> getRowSums() const;
	double getColSum(unsigned int const col) const;
	vector<double> getColSums() const;

	void inverse();
	double determinant() const;

	void luDecomposition();
	Matrix luEvaluate(Matrix &a);
	void print() const;


private:
	unsigned int _dim;
	vector< vector<double> > _m;
};

#endif /* MATRIX_H_ */
