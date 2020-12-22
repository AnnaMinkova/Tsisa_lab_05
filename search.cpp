#include <iostream>
#include <random>
#include <vector>
#include <algorithm>

//значения из варианта
const int N = 24;
const double a = -2.0;
const double b = 2.0;
const double stage = (double)((b - a) / (N - 1));
const double c = 1.0;
const double d = 0.0;
const double A = 2.0;

double function(const double &x) {
    return c * x + d;
}
struct Pair {
    double x;
    double y;
};


std::vector<std::pair<double, double>> f(24);

double recurssion_func(double C, double D) {
    double sum_E = 0;
    for (int i = 0; i < N; i++) {
        double x = f[i].first;
        double t = f[i].second;
        double y = C * x + D;
        sum_E = sum_E + pow((y - t), 2);
    }
    return sum_E;
}

std::vector<Pair> rand(const int N, const double noise) {
    std::vector<Pair> pair (N);
    std::random_device rand;
    std::mt19937 gen(rand());
    std::uniform_real_distribution<double> er(-0.5, 0.5);
    for (size_t i = 0; i < N; ++i) {
        pair[i].x = a + i * stage;
        pair[i].y = function( a + i * stage) + noise * er(gen);
    }
    return pair;
}

// метод дихотомии, определяющий коэффициент "C"
double Metod_Dichotomy(double a, double b, double c){
    const double Epsilon = 0.1;
    const double Delta = 0.01;
    double x_left, x_right, y_left, y_right;
    do
    {
        x_left = 0.5 * (b + a) - Delta;
        x_right = 0.5 * (b + a) + Delta;
        y_left = recurssion_func(c, x_left);
        y_right = recurssion_func(c, x_right);
        if (y_left > y_right)
        {
            a = x_left;
        }
        else{
            b = x_right;
        }
    } while ((b - a) > Epsilon);
    return (a + b) / 2;
}

// метод золотого сечения, определяющий коэффициент "D"
double Method_golden(std::vector<Pair>& p, double Dmin, double Dmax) {
    double length = std::abs(Dmax - Dmin);
    std::swap(Dmin, Dmax);
    Dmin = std::fabs(Dmin);
    Dmax = std::fabs(Dmax);
    const double e = 0.1;
    const double t = (std::sqrt(5) + 1) / 2;
    double d_k1 = Dmin + (1 - 1/t)*Dmax;
    double d_k2 = Dmin + Dmax / t;
    double f_k1 = function(-d_k1);
    double f_k2 = function(-d_k2);
    while (length > e){
        if (f_k1 < f_k2){
            Dmax = d_k2;
            d_k2 = Dmin + Dmax - d_k1;
            f_k2 = function(-d_k2);
        } else {
            Dmin = d_k1;
            d_k1 = Dmin + Dmax - d_k2;
            f_k1 = function(-d_k1);
        }
        if (d_k1 > d_k2){
            std::swap(d_k1, d_k2);
            std::swap(f_k1, f_k2);
        }
        length = std::abs(Dmax - Dmin);
    }
    return -((Dmax + Dmin) / 2);
}

void print(const double noise)
{
    std::vector<Pair> p = rand(N, noise);
    double Cmin, Cmax, Dmin, Dmax;
    std::cout << "Cmin = " << Cmin << "\nCmax = " << Cmax << "\nDmin = " << Dmin << "\nDmax = " << Dmax << std::endl;
    double w1 = Metod_Dichotomy(a,b,c); //коэффициент перед x
    double w0 = Method_golden(p, Cmin, Cmax); //свободный коэффициент
    std::cout << "w1 = " << w1 <<std::endl;
    std::cout<< "w0 = " << w0 << std::endl;
}

int main() {

    std::cout << "Function y = x + 0"<<std::endl; // c=1, d=0
    std::cout << "Without noise"<<std::endl;
    print(0.0);
    std::cout << "With noise A = " << A <<std::endl;
    print(A);
    return 0;
}
