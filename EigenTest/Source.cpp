#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>

#include <time.h>
#include <fstream>
using namespace std;
using namespace Eigen;

ofstream g_debug("debug.txt");
/*
void ReadAFromFileX(const char* filename, MatrixXd *A)
{
	ifstream fileStream(filename, ios::in);
	int nRow, nCol, nNZero;
	if (fileStream.is_open())
	{
		fileStream >> nRow >> nCol >> nNZero;
		cout << nNZero << endl;
		*A = MatrixXd::Constant(nRow, nCol, 0);
		int row, col;
		for (int i = 0; i < nNZero; ++i)
		{
			fileStream >> row >> col;
			row -= 1;
			col -= 1;
			fileStream >> (*A)(row, col);
			//g_debug << row << " " << col << " " << (*A)(row, col) << endl;
		}
		fileStream.close();
	}
}*/
SparseMatrix<double> ReadAFromFile(const char* filename)
{
	ifstream fileStream(filename, ios::in);
	int nRow, nCol, nNZero;
	SparseMatrix<double> A;
	if (fileStream.is_open())
	{
		fileStream >> nRow >> nCol >> nNZero;
		typedef Eigen::Triplet<double> T;
		std::vector<T> tripletList;
		tripletList.reserve(nNZero);

		int row, col;
		double value;
		for (int i = 0; i < nNZero; ++i)
		{
			fileStream >> row >> col >> value;
			row -= 1;
			col -= 1;
			tripletList.push_back(T(row, col, value));
		}
		A.resize(nRow, nCol);
		A.setFromTriplets(tripletList.begin(), tripletList.end());
		//g_debug << A << endl;
		fileStream.close();
	}
	return A;
}
VectorXd ReadBFromFile(const char* filename)
{
	ifstream fileStream(filename, ios::in);
	int nRow;
	VectorXd b;
	if (fileStream.is_open())
	{
		fileStream >> nRow;
		b.resize(nRow);
		for (int i = 0; i < nRow; ++i)
		{
			fileStream >> b(i);
			//g_debug << b(i) << endl;
		}
		fileStream.close();
	}
	return b;
}
int main()
{
	SparseMatrix<double> A;
	VectorXd b, x;

	A = ReadAFromFile("A675.dat");
	b = ReadBFromFile("b675.dat");
	//g_debug << A << endl;
	//g_debug << b << endl;

	double curTime;
	curTime = (double)clock() / CLOCKS_PER_SEC;

	SimplicialLLT<SparseMatrix<double> > solver;
	solver.compute(A);
	if (solver.info() != Success) {
		// decomposition failed
		return -1;
	}
	x = solver.solve(b);
	if (solver.info() != Success) {
		// solving failed
		return -1;
	}
	g_debug << x << endl;

	double preTime = (double)clock() / CLOCKS_PER_SEC;
	cout << "Solver1 time used = " << preTime - curTime << endl;


	SimplicialLDLT<SparseMatrix<double> > solver2;
	solver2.compute(A);
	if (solver2.info() != Success) {
		// decomposition failed
		return -1;
	}
	x = solver2.solve(b);
	if (solver2.info() != Success) {
		// solving failed
		return -1;
	}
	g_debug << "-------------------------------------------------------" << endl;
	g_debug << x << endl;

	double preTime2 = (double)clock() / CLOCKS_PER_SEC;
	cout << "Solver2 time used = " << preTime2 - preTime << endl;

	ConjugateGradient<SparseMatrix<double> > solver3;
	solver3.compute(A);
	if (solver3.info() != Success) {
		// decomposition failed
		return -1;
	}
	x = solver3.solve(b);
	if (solver3.info() != Success) {
		// solving failed
		return -1;
	}
	g_debug << "-------------------------------------------------------" << endl;
	g_debug << x << endl;

	cout << "Solver2 time used = " << (double)clock() / CLOCKS_PER_SEC - preTime2 << endl;

	system("pause");
	return 0;
}