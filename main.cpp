#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>

const double EPS = 1e-15;
const int MAX_ITER = 10000;
const double PI = std::acos(-1.0);

struct IterationData {
    std::vector<int> iterations;
    std::vector<double> E_values;
    std::vector<double> f_values;
};

inline double f(double E, double e, double M) {
    return E - e * std::sin(E) - M;
}

double solve_bisection(double e, double M, IterationData &data) {
    if (e < EPS) {
        return M;
    }
    double a = std::max(0.0, M - e);
    double b = std::min(2.0 * PI, M + e);
    if (f(a, e, M) * f(b, e, M) > 0) {
        a = 0.0;
        b = 2.0 * PI;
    }
    data.iterations.clear(); data.E_values.clear(); data.f_values.clear();
    double fa = f(a, e, M), fb = f(b, e, M), c = a, fc;
    for (int iter = 1; iter <= MAX_ITER; ++iter) {
        c = 0.5 * (a + b);
        fc = f(c, e, M);
        data.iterations.push_back(iter);
        data.E_values.push_back(c);
        data.f_values.push_back(fc);
        if (std::abs(fc) < EPS || (b - a) / 2.0 < EPS) break;
        if (fa * fc <= 0.0) { b = c; fb = fc; } else { a = c; fa = fc; }
    }
    return c;
}

double solve_regula_falsi(double e, double M, IterationData &data) {
    if (e < EPS) {
        return M;
    }
    double a = std::max(0.0, M - e);
    double b = std::min(2.0 * PI, M + e);
    if (f(a, e, M) * f(b, e, M) > 0) {
        a = 0.0;
        b = 2.0 * PI;
    }
    data.iterations.clear(); data.E_values.clear(); data.f_values.clear();
    double fa = f(a, e, M), fb = f(b, e, M), c = a, fc;
    for (int iter = 1; iter <= MAX_ITER; ++iter) {
        if (std::abs(fb - fa) < EPS) break;
        c = b - fb * (b - a) / (fb - fa);
        fc = f(c, e, M);
        data.iterations.push_back(iter);
        data.E_values.push_back(c);
        data.f_values.push_back(fc);
        if (std::abs(fc) < EPS) break;
        if (fa * fc <= 0.0) { b = c; fb = fc; } else { a = c; fa = fc; }
    }
    return c;
}

void save_all_iterations(const std::string &csv_filename,
                         double e, double m,
                         const std::string &method,
                         const IterationData &data,
                         std::ofstream &csv) {
    for (size_t i = 0; i < data.iterations.size(); ++i) {
        csv << e << ","
            << m << ","
            << method << ","
            << data.iterations[i] << ","
            << std::fixed << std::setprecision(15)
            << data.E_values[i] << ","
            << data.f_values[i] << "\n";
    }
}

int main() {
    std::vector< std::pair<double, double> > tests;
    tests.push_back(std::make_pair(0.0, 0.0));
    tests.push_back(std::make_pair(0.0, 1.0));
    tests.push_back(std::make_pair(0.1, 0.0));
    tests.push_back(std::make_pair(0.1, 0.5));
    tests.push_back(std::make_pair(0.1, 1.5));
    tests.push_back(std::make_pair(0.5, 0.25));
    tests.push_back(std::make_pair(0.5, 1.0));
    tests.push_back(std::make_pair(0.5, 1.75));
    tests.push_back(std::make_pair(0.9, 0.2));
    tests.push_back(std::make_pair(0.9, 1.8));
    tests.push_back(std::make_pair(0.3, 0.0));
    tests.push_back(std::make_pair(0.3, 2.0));
    tests.push_back(std::make_pair(0.8, 2.0));
    tests.push_back(std::make_pair(0.999, 1.0));
    tests.push_back(std::make_pair(0.2, 0.75));
    tests.push_back(std::make_pair(0.2, 1.25));

    std::ofstream results_out("results.txt");
    std::ofstream csv_out("iterations.csv");
    csv_out << "e,m,method,iteration,E,f(E)\n";

    IterationData data;
    for (std::vector< std::pair<double, double> >::iterator it = tests.begin(); it != tests.end(); ++it) {
        double e = it->first;
        double m = it->second;
        double M = m * PI;
        results_out << "\n===== Тест: e=" << e << ", m=" << m << " (M=" << M << ") =====\n";

        results_out << "-- Бисекция --\n";
        double E_bis = solve_bisection(e, M, data);
        save_all_iterations("iterations.csv", e, m, "bisection", data, csv_out);

        results_out << "-- Регула Фальси --\n";
        double E_false = solve_regula_falsi(e, M, data);
        save_all_iterations("iterations.csv", e, m, "regula_falsi", data, csv_out);

        results_out << "Результаты: Bisection E=" 
                    << std::fixed << std::setprecision(15) << E_bis 
                    << ", Regula Falsi E=" << E_false 
                    << ", delta=" << std::abs(E_bis - E_false) << "\n";
    }

    results_out.close();
    csv_out.close();
    std::cout << "Тестирование завершено. Итерации сохранены в iterations.csv, результаты в results.txt\n";
    return 0;
}
