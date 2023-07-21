//#include <armadillo>
//using namespace arma;
#include "ThisArrays.cpp"
using namespace std;
class WaveEquationSolver {
public:
	WaveEquationSolver(int N = 80, int T = 10) : N(N), T(T) {
		h = 2.0 / (N - 1);
		dt = (0.9 * 2 / 7 * h);
		tn = round(T / dt);
		u= ThisArrays(N, tn);
		u2= ThisArrays(N, tn);
		x = linspace(-1, 1, N);
		for (int i = 0; i < N; i++) {
			u(i, 0) = sin(pi * x(i));
		}

		// Theory solution u(x,t) = sin(pi*(x-t))
		t = linspace(0, tn * dt, tn);
		for (int i = 0; i < tn; i++) {
			for (int j = 0; j < N; j++) {
				u2(j, i) = sin(pi * (x(j) - t(i)));
			}
		}
	}

	void solve() {
		for (int t = 1; t < tn; t++) {
			ThisArrays K1 = deriv_center(u.col(t - 1));
			ThisArrays K2 = deriv_center(u.col(t - 1) + 0.5 * dt * K1);
			ThisArrays K3 = deriv_center(u.col(t - 1) + 0.5 * dt * K2);
			ThisArrays K4 = deriv_center(u.col(t - 1) + dt * K3);
			u.col(t) = u.col(t - 1) + dt * (1.0 / 6 * K1 + 1.0 / 3 * K2 + 1.0 / 3 * K3 + 1.0 / 6 * K4);
			
			u(0, t) = u(N - 1, t);
		}
		u.save("u");
		u2.save("u2");
//		u.save("u.csv", csv_ascii);
//		u2.save("u2.csv", csv_ascii);

	}


	void Accuracy() {
		error = abs(u.col(tn - 1) - u2.col(tn - 1));
		std::cout << "L1_norm:" << endl << accu(error) * h << endl;
		std::cout << "Linf_norm:" << endl << max(error) << endl;
		std::cout << "L_order()" << endl << accu(error) / N << endl;
	}
private:
	const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
	int N, T;
	double dt;
	int tn;
	double h;
	ThisArrays u;
	ThisArrays u2;
	ThisArrays x;
	ThisArrays t;
	ThisArrays error;

	ThisArrays deriv_center(const ThisArrays& u) {
		ThisArrays f = u;
		ThisArrays df(N);
		for (int i = 2; i < N - 2; i++) {
			df(i) = (-f(i + 2) + 8 * f(i + 1) - 8 * f(i - 1) + f(i - 2)) / (12 * h);
		}
		df(1) = (-f(3) + 8 * f(2) - 8 * f(0) + f(N - 2)) / (12 * h);
		df(0) = (-f(2) + 8 * f(1) - 8 * f(N - 2) + f(N - 3)) / (12 * h);
		df(N - 1) = df(0);
		df(N - 2) = (-f(1) + 8 * f(N - 1) - 8 * f(N - 3) + f(N - 4)) / (12 * h);
		return -df;
	}
};
int main() {
	auto simulate = WaveEquationSolver();
	simulate.solve();

	return 0;
}