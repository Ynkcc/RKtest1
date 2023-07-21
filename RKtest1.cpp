#include <armadillo>
using namespace arma;

class WaveEquationSolver {
public:
	WaveEquationSolver(int N = 80, int T = 10) : N(N), T(T) {
		h = 2.0 / (N - 1);
		dt = (0.9 * 2 / 7 * h);//0.9
		tn = round(T / dt);
		u = mat(N, tn);
		u2 = mat(N, tn);
		for (int i = 0; i < N; i++) {
			u(i, 0) = sin(datum::pi * i*h);
		}

		// Theory solution u(x,t) = sin(datum::pi*(x-t))
		for (int i = 0; i < tn; i++) {
			for (int j = 0; j < N; j++) {
				u2(j, i) = sin(datum::pi *(j*h - i*dt));
			}
		}
	}

	void solve() {
		for (int t = 1; t < tn; t++) {
			vec K1 = deriv_center(u.col(t - 1));
			vec K2 = deriv_center(u.col(t - 1) + 0.5 * dt * K1);
			vec K3 = deriv_center(u.col(t - 1) + 0.5 * dt * K2);
			vec K4 = deriv_center(u.col(t - 1) + dt * K3);
			u.col(t) = u.col(t - 1) + dt * (1.0 / 6 * K1 + 1.0 / 3 * K2 + 1.0 / 3 * K3 + 1.0 / 6 * K4);

			u(0, t) = u(N - 1, t);
		}
		u.save("u.csv", csv_ascii);
		u2.save("u2.csv", csv_ascii);

	}


	void Accuracy() {
		error = abs(u.col(tn - 1) - u2.col(tn - 1));
		std::cout << "L1_norm:" << endl << accu(error) * h << endl;
		std::cout << "Linf_norm:" << endl << max(error) << endl;
		std::cout << "L_order()" << endl << accu(error) / N << endl;
	}
private:

	int N, T;
	double dt;
	int tn;
	double h;
	mat u;
	mat u2;
	vec x;
	vec t;
	vec error;

	vec deriv_center(const vec& u) {
		vec f = u;
		vec df(N);
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
	auto simulate = WaveEquationSolver(80,10);
	simulate.solve();
	simulate.Accuracy();
	system("pause");
	return 0;
}