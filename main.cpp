//
//file name: main.cpp
//due date: 11.1.16
//author: Josh Morris
//course: CSCI 3180
//Program Description: Program that constructs the natural cubic spline for the function 1/(x^2 + 1)
//					using the points x = -6, -4, -2, 0, 2, 4, 6
//	input: none
//	output: data points
//			tri-diagonal system before forward elim
//			tri-diagonal system after forward elim
//			solution for z_i
//			coefficients Ai, Bi, Ci, Di 
//			Natural Cubic Spline 
//			Table containing x, f(x), S(x), and |f(x) - S(x)| for x = -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

using std::vector; using std::cin; using std::cout; using std::endl;
using std::setw; using std::endl; using std::fabs; using std::fixed;
using std::showpoint;

double Equation(double x);
void genDataPoints(vector<vector<double>>& dataPoints);
void genh_i(const vector<vector<double>>& dataPoints, vector<double> &h_i);
void genu_i(const vector<double> &h_i, vector<double> &u_i);
void genv_i(const vector<vector<double>>& dataPoints, const vector<double> &h_i, vector<double> &v_i);
void displayTDS(const vector<double>& u_i, const vector<double>& h_i, const vector<double> &v_i);
void forwardElim(vector<double>& u_i, vector<double>& h_i, vector<double>& v_i);
void displayTDSFE(const vector<double>& u_i, const vector<double>& h_i, const vector<double>& v_i);
void calcZ(vector<double>& u_i, vector<double>& h_i, vector<double>& v_i, vector<double> &z_i);
void calcABCD(const vector<double>& z_i, const vector<vector<double>> & dataPoints, const vector<double> &h_i, vector<vector<double>> &ABCD_i);
void displayNestSi(const vector<vector<double>> & dataPoints, const vector<vector<double>> &ABCD_i);
double calcSi(double x, const vector<vector<double>> & dataPoints, const vector<vector<double>> &ABCD_i);
void finalTable(const vector<vector<double>> & dataPoints, const vector<vector<double>> &ABCD_i);

int main() {
	vector<vector<double>> dataPoints(2);
	vector<vector<double>> testPoints(2);
	vector<double> h_i;
	vector<double> u_i;
	vector<double> v_i;
	vector<double> z_i;
	vector<vector<double>> ABCD_i (6, vector<double> (4, 0));

	// compute y values using equations
	genDataPoints(dataPoints);

	// calculate h_i, u_i, v_i
	genh_i(dataPoints, h_i);
	genu_i(h_i, u_i);
	genv_i(dataPoints, h_i, v_i);

	// solve for z_i
	displayTDS(u_i, h_i, v_i);
	forwardElim(u_i, h_i, v_i);
	displayTDSFE(u_i, h_i, v_i);
	calcZ(u_i, h_i, v_i, z_i);

	// compute Ai, Bi, Ci, Di
	calcABCD( z_i, dataPoints, h_i, ABCD_i);

	// display S_i(x)
	displayNestSi(dataPoints, ABCD_i);

	// print results
	finalTable(dataPoints, ABCD_i);

	return 0;
}

double Equation(double x) {
	return (1 / (x*x + 1));
}

void genDataPoints(vector<vector<double>>& dataPoints) {
	int i = 0;
	cout << "Data points (ti, yi) " << endl;
	
	for (double t = -6; t <= 6; t += 2) {
		// gen points
		dataPoints[0].push_back(t);
		dataPoints[1].push_back(Equation(t));

		//print data point
		cout << "(t" << i << ", y" << i << ") = (" << t << "," << dataPoints[1][i] << ")" << endl;

		++i;
	}

	cout << endl;

	return;
}

void genh_i(const vector<vector<double>>& dataPoints, vector<double> &h_i) {
	for (auto j = dataPoints[0].begin(); j != (dataPoints[0].end() - 1); ++j) {
		h_i.push_back(*(j + 1) - *j);
	}

	return;
}

void genu_i(const vector<double> &h_i, vector<double> &u_i) {
	for (auto j = (h_i.begin() + 1); j != (h_i.end()); ++j) {
		u_i.push_back(2*(*(j - 1) + *j));
		//cout << u_i.back() << endl;
	}

	return;
}

void genv_i(const vector<vector<double>>& dataPoints, const vector<double> &h_i, vector<double> &v_i) {
	for (int i = 1; i <= 5; ++i) {
		v_i.push_back(6*((1/h_i[i])*(dataPoints[1][i+1] - dataPoints[1][i]) - (1/h_i[i-1])*(dataPoints[1][i] - dataPoints[1][i-1])));
		//cout << v_i.back() << endl;
	}
}

void displayTDS(const vector<double>& u_i, const vector<double>& h_i, const vector<double>& v_i) {
	cout << "Tri-diagonal system for z1 through z(n-1)" << endl;
	
	cout << u_i.front() << "\t" << h_i.front() << "\t0\t0\t0\t" << v_i.front() << endl;
	
	for (int i = 1; i < 4; ++i) {
		for (int j = 0; j < (i - 1); ++j) {
			cout << "0 \t";
		}

		cout << h_i[i - 1] << "\t" << u_i[i] << "\t" << h_i[i] << "\t";

		for (int k = 0; k < (3 - i); ++k) {
			cout << "0 \t";
		}

		cout << v_i[i] << endl;
	}

	cout << "0\t0\t0\t" << h_i.back() << "\t" << u_i.back() << "\t" << v_i.back() << endl << endl;
	return;
}

void forwardElim(vector<double>& u_i, vector<double>& h_i, vector<double>& v_i) {
	double m;

	//cout << "size of h_i: " << h_i.size() << endl << endl; 
	//cout << "size of u_i: " << u_i.size() << endl << endl;

	for (int i = 0; i < 4; ++i) {
		m = h_i[i + 1] / u_i[i];
		u_i[i + 1] -= m * h_i[i + 1];
		v_i[i + 1] -= m * v_i[i];
	}

	return;
}

void displayTDSFE(const vector<double>& u_i, const vector<double>& h_i, const vector<double>& v_i) {
	cout << "Tri-diagonal system after Forward Elimination" << endl;

	cout << u_i.front() << "\t" << h_i.front() << "\t0\t0\t0\t" << v_i.front() << endl;

	for (int i = 1; i < 4; ++i) {
		for (int j = 0; j < (i - 1); ++j) {
			cout << "0 \t";
		}

		cout << "0 \t" << u_i[i] << "\t" << h_i[i] << "\t";

		for (int k = 0; k < (3 - i); ++k) {
			cout << "0 \t";
		}

		cout << v_i[i] << endl;
	}

	cout << "0\t0\t0\t0\t" << u_i.back() << "\t" << v_i.back() << endl << endl;
	return;
}

void calcZ(vector<double>& u_i, vector<double>& h_i, vector<double>& v_i, vector<double> &z_i) {
	int i = 0;

	//set size of vector to do back substitution for z values
	z_i.resize(6);
	z_i[5] = 0;

	// back substitution
	z_i[4] = v_i[4] / u_i[4];

	for (int i = 3; i >= 0; --i) {
		z_i[i] = (v_i[i] - h_i[i] * z_i[i + 1]) / u_i[i];
	}

	// add z_0 = 0 to front of z values
	z_i.insert(z_i.begin(), 0);


	// display z values
	cout << "Solutions for S''(ti) = zi" << endl;
	cout << "------------------------------------" << endl;

	for (auto iter = z_i.begin(); iter != z_i.end(); ++iter) {
		cout << "z" << i << ": " << *iter << endl;
		++i;
	}

	cout << endl;
	return;
}

void calcABCD(const vector<double>& z_i, const vector<vector<double>> & dataPoints, const vector<double>& h_i, vector<vector<double>> &ABCD_i) {
	cout << "Coefficients for Natural Cubic Splin in nested form" << endl;
	cout << "-------------------------------------------------------------" << endl;
	
	// for all S(x)
	for (int i = 0; i < 6; ++i) {
		// calculate A
		ABCD_i[i][0] = dataPoints[1][i];

		cout << "A" << i << ":\t" << setw(8) << ABCD_i[i][0] << "\t";

		// Calculate B
		ABCD_i[i][1] = (-h_i[i] / 6)*z_i[i + 1] - (h_i[i] / 3)*z_i[i] + (dataPoints[1][i + 1] - dataPoints[1][i]) / h_i[i];

		cout << "B" << i << ":\t" << setw(8) << ABCD_i[i][1] << "\t";

		// Calc C
		ABCD_i[i][2] = z_i[i] / 2;

		cout << "C" << i << ":\t" << setw(8) << ABCD_i[i][2] << "\t";

		// Calc D
		ABCD_i[i][3] = (z_i[i + 1] - z_i[i]) / (6 * h_i[i]);

		cout << "D" << i << ":\t" << setw(8) << ABCD_i[i][3] << endl;
	}

	cout << endl;

	return;
}

void displayNestSi(const vector<vector<double>> & dataPoints, const vector<vector<double>> &ABCD_i) {
	cout << "Natural Cubic Spline in Nested Form" << endl;
	cout << "-------------------------------------------------------------" << endl;
	
	for (unsigned int i = 0; i < ABCD_i.size(); ++i) {
		cout << "S" << i << " : " << ABCD_i[i][0] << " + (x - " << dataPoints[0][i] << ") * (" << ABCD_i[i][1] << " + (x - " << dataPoints[0][i] << ") * ("
			<< ABCD_i[i][2] << " + (x - " << dataPoints[0][i] << ") * (" << ABCD_i[i][3] << ")))" << endl;
	}

	cout << endl;

	return;
}

double calcSi(double x, const vector<vector<double>> & dataPoints, const vector<vector<double>> &ABCD_i) {
	int SiMax = ABCD_i.size();

	for (unsigned int i = 0; i < 6; ++i) {
		if (x <= dataPoints[0][i + 1]) {
			return ABCD_i[i][0] + (x - dataPoints[0][i]) * (ABCD_i[i][1] + (x - dataPoints[0][i]) * (ABCD_i[i][2] + (x - dataPoints[0][i]) * (ABCD_i[i][3])));
		}
	}
}

void finalTable(const vector<vector<double>> & dataPoints, const vector<vector<double>> &ABCD_i) {
	double f_xi, s_xi;
	int i = 0;

	cout << fixed << showpoint;
	cout << "Evaluation of Original and Spline functions and the absolute errors" << endl;

	cout << "---------------------------------------------------------------------------------------" << endl;
	cout << "i\t" << setw(11) << "xi\t" << setw(11) << "f(xi)\t" << setw(11) << "S(xi)\t" << setw(11) << "|f(xi) - S(xi)|" << endl;
	cout << "---------------------------------------------------------------------------------------" << endl;

	for (int x = -6; x <= 6; ++x) {
		f_xi = Equation(x);
		s_xi = calcSi(x, dataPoints, ABCD_i);

		cout << i << "\t";
		cout << setw(11) << x << "\t";
		cout << setw(11) << f_xi << "\t";
		cout << setw(11) << s_xi << "\t";
		cout << setw(11) << fabs(f_xi - s_xi) << "\t" << endl;

		++i;
	}

	return;
}