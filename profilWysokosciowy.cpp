#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <ctime>

#define Liczba_TESTOW 4
#define LICZBA_PROBEK 480
using namespace std;

struct probka
{
	double x;
	double y;
};

bool wczytaj_dane(const char*, probka*);
void LU_decomposition(double**, double*, double*, int);
double Lagrange(probka*, double, const int);
double Splines(const probka*, double, const int);
void Interpolacja(const probka*, const int, probka*, const char*, double*, int);

int main()
{
	const char* testy[Liczba_TESTOW] = { "genoa_rapallo.txt","ostrowa.txt","chelm.txt" , "ulm_lugano.txt", };
	int tests[Liczba_TESTOW] = { 16, 40, 48, 80 };
	double average[Liczba_TESTOW] = { 0 };
	probka* probki = new probka[LICZBA_PROBEK];

	for (int i = 0; i < Liczba_TESTOW; i++)
	{
		const char* filename = testy[i];

		string file_path = "..\\dane\\";
		file_path.append(filename);
		filename = file_path.c_str();

		string line;
		int count = 0;
		ifstream f(filename);
		while (getline(f, line))
			count++;

		wczytaj_dane(filename, probki);
		int j = 0;


		for (int j = 0; j < Liczba_TESTOW; j++)
		{
			int krok = tests[j];

			const int numer = LICZBA_PROBEK / krok;
			probka* wezly = new probka[numer];

			for (int k = 0, l = 0; l < numer; k += krok, l++)
			{
				wezly[l].x = probki[k].x;
				wezly[l].y = probki[k].y;

			}

			Interpolacja(probki, numer, wezly, filename, &average[j], count);
		}
	}

	for (int i = 0; i < Liczba_TESTOW; i++)
	{
		int liczba_wezlow = LICZBA_PROBEK / tests[i];
		double czas = average[i] / Liczba_TESTOW;
		cout << "liczba wezlow = " << liczba_wezlow << "   " << "czas = " << czas << "ms" << endl;
	}
	return 0;
}


bool wczytaj_dane(const char* filename, probka* samples)
{
	ifstream file;
	file.open(filename, file.in);

	if (!file.is_open())
	{
		printf_s("Failed to open the specified file.");
		return false;
	}

	string sample;

	for (int i = 0; i < LICZBA_PROBEK; i++)
	{
		getline(file, sample, ' ');
		samples[i].x = stod(sample);
		getline(file, sample);
		samples[i].y = stod(sample);
	}

	file.close();

	return true;
}

void LU_decomposition(double**M, double* b, double* x, int N)
{
	double** U = new double*[N];		//upper triangular matrix
	double** L = new double*[N];		//upper triangular matrix
	double** P = new double*[N];		//upper triangular matrix

	for (int i = 0; i < N; i++) {
		U[i] = new double[N];
		L[i] = new double[N];
		P[i] = new double[N];
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			U[i][j] = M[i][j];
			if (i == j)
				L[i][j] = P[i][j] = 1;
			else
			{
				L[i][j] = 1;
				P[i][j] = 1;
			}
		}
	}
	for (int i = 0; i < N - 1; i++)
	{
		double pivot = abs(U[i][i]);
		int pivot_index = i;

		for (int j = i + 1; j < N; j++)
		{
			if (abs(U[j][i]) > pivot)
			{
				pivot = abs(U[j][i]);
				pivot_index = j;
			}
		}

		if (pivot_index != i)
		{
			double tmp;

			for (int j = 0; j < N; j++)
			{
				if (j >= i)
				{
					tmp = U[i][j];
					U[i][j] = U[pivot_index][j];
					U[pivot_index][j] = tmp;
				}
				else
				{
					tmp = L[i][j];
					L[i][j] = L[pivot_index][j];
					L[pivot_index][j] = tmp;
				}

				tmp = P[i][j];
				P[i][j] = P[pivot_index][j];
				P[pivot_index][j] = tmp;
			}
		}

		for (int j = i + 1; j < N; j++)
		{
			L[j][i] = U[j][i] / U[i][i];

			for (int k = i; k < N; k++)
				U[j][k] = U[j][k] - L[j][i] * U[i][k];
		}
	}
	for (int i = 0; i < N; i++)
		b[i] = (**P) * b[i];

	double* y = new double[N];

	for (int i = 0; i < N; i++)
	{
		double S = 0;

		for (int j = 0; j < i; j++)
			S += L[i][j] * y[j];

		y[i] = (b[i] - S) / L[i][i];
	}

	for (int i = N - 1; i >= 0; i--)
	{
		double S = 0;

		for (int j = i + 1; j < N; j++)
			S += U[i][j] * x[j];

		x[i] = (y[i] - S) / U[i][i];
	}

	for (int i = 0; i < N; i++) {
		delete[]U[i];
		delete[]L[i];
		delete[]P[i];

	}
	delete[]U;
	delete[]L;
	delete[]P;
}

double Lagrange(probka* wezel, double distance, const int n)
{
	double elevation = 0.0;

	for (int i = 0; i < n; i++)
	{
		double a = 1;

		for (int j = 0; j < n; j++)
		{
			if (i != j)
				a *= (distance - wezel[j].x) / (wezel[i].x - wezel[j].x);
		}

		elevation += a * wezel[i].y;
	}
	return elevation;
}

double Splines(const probka* samples, double distance, const int nodes_number)
{
	int N = 4 * (nodes_number - 1);
	double** M = new double*[N];

	for (int i = 0; i < N; i++)
		M[i] = new double[N];

	double* b = new double[N];
	double* x = new double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			M[i][j] = 1;


	for (int i = 0; i < N; i++)
	{
		b[i] = 0;
		x[i] = 1;
	}
	M[0][0] = 1;
	b[0] = samples[0].y;

	double h;

	h = samples[1].x - samples[0].x;
	M[1][0] = 1;
	M[1][1] = h;
	M[1][2] = (h*h);
	M[1][3] = (h*h*h);
	b[1] = samples[1].y;

	M[2][2] = 1;
	b[2] = 0;

	h = samples[nodes_number - 1].x - samples[nodes_number - 2].x;
	M[3][4 * (nodes_number - 2) + 2] = 2;
	M[3][4 * (nodes_number - 2) + 3] = 6 * h;
	b[3] = 0;

	for (int i = 1; i < nodes_number - 1; i++)
	{
		h = samples[i].x - samples[i - 1].x;

		M[4 * i][4 * i] = 1;
		b[4 * i] = samples[i].y;

		M[4 * i + 1][4 * i] = 1;
		M[4 * i + 1][4 * i + 1] = h;
		M[4 * i + 1][4 * i + 2] = h * h;
		M[4 * i + 1][4 * i + 3] = h * h;
		b[4 * i + 1] = samples[i + 1].y;

		M[4 * i + 2][4 * (i - 1) + 1] = 1;
		M[4 * i + 2][4 * (i - 1) + 2] = 2 * h;
		M[4 * i + 2][4 * (i - 1) + 3] = 3 * h*h;
		M[4 * i + 2][4 * i + 1] = -1;
		b[4 * i + 2] = 0;

		M[4 * i + 3][4 * (i - 1) + 2] = 2;
		M[4 * i + 3][4 * (i - 1) + 3] = 6 * h;
		M[4 * i + 3][4 * i + 2] = -2;
		b[4 * i + 3] = 0;
	}

	LU_decomposition(M, b, x, N);

	double elevation = 0;
	for (int i = 0; i < nodes_number - 1; i++)
	{
		elevation = 0;

		if (distance >= samples[i].x && distance <= samples[i + 1].x)
		{
			for (int j = 0; j < 4; j++)
				elevation += x[4 * i + j] * pow(distance - samples[i].x, j);


			break;
		}
	}

	for (int i = 0; i < N; i++)
		delete[]M[i];
	delete[] x;
	delete[] b;
	delete[]M;

	return elevation;
}

void Interpolacja(const probka* samples, const int numer_wezla, probka* wezel, const char* path, double* czas, int count)
{
	string file_path = path;
	file_path.replace(3, 4, "results\\Lagrange");
	//file_path.replace(3, 4, "results\\splines");
	file_path.replace(file_path.end() - 4, file_path.end(), "_");
	file_path.append(to_string(numer_wezla));
	file_path.append(".csv");


	ofstream file;
	file.open(file_path, file.out);

	if (!file.is_open())
	{
		printf_s("Failed to save the results.");
		return;
	}

	double result;
	clock_t start = clock();

	int k = 0;
	int g = 0;
	int ile = 0;
	for (double i = wezel[0].x; i <= wezel[numer_wezla - 1].x; i += 8)
		ile++;
	int ilee = ile / 480;

	for (double i = wezel[0].x; i <= wezel[numer_wezla - 1].x; i += 8)
	{
		bool interpolate = true;

		for (int j = 0; j < numer_wezla; j++)
		{

			if ((int)wezel[j].x == i)
			{
				interpolate = false;
				result = wezel[j].y;

				break;
			}
		}

		if (interpolate) {
			result = Lagrange(wezel, i, numer_wezla);
		}
		//	result = Splines(wezel, i, numer_wezla);

		if (ilee != 0) {
			if (k % ilee == 0 && g < 480) {
				file << i << " " << result << " " << samples[g].y << endl;
				g++;
			}
			else
				file << i << " " << result << endl;
			k++;
		}
		else {
			file << i << " " << result << " " << samples[g].y << endl;

		}

	}

	file.close();
	*czas += clock() - start;

	file_path.replace(file_path.end() - 4, file_path.end(), "_nodes.csv");

	file.open(file_path, file.out);

	if (!file.is_open())
	{
		printf_s("Failed to save the results.");
		return;
	}

	for (int i = 0; i < numer_wezla; i++)
		file << wezel[i].x << " " << wezel[i].y << endl;

	file.close();
}