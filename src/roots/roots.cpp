#include <roots.hpp>
#include <cmath>


static constexpr double TOL = 1e-6; //numerical tolerance ('double' since it is a floating point number)
static constexpr int MAX_ITER = 1000000; //maximum number of iterations

//static constexpr means that this file only exists in this file and the value is known at compile time

/* Bisection Method */
bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root) //this allows us to choose a location where we will store the answer
{
    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0.0)
        return false;

    for (int i = 0; i < MAX_ITER; ++i)
    {
        double c = 0.5 * (a + b);
        double fc = f(c);

        if (std::abs(fc) < TOL || (b - a) * 0.5 < TOL)
        {
            *root = c;
            return true;
        }

        if (fa * fc < 0.0)
        {
            b = c;
            fb = fc;
        }
        else
        {
            a = c;
            fa = fc;
        }
    }

    return false;
}

/* Regula Falsi (False Position) */
bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root)
{
    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0.0)
        return false;

    for (int i = 0; i < MAX_ITER; ++i)
    {
        double c = b - fb * (b - a) / (fb - fa);
        double fc = f(c);

        if (std::abs(fc) < TOL)
        {
            *root = c;
            return true;
        }

        if (fa * fc < 0.0)
        {
            b = c;
            fb = fc;
        }
        else
        {
            a = c;
            fa = fc;
        }
    }

    return false;
}

/* Newton–Raphson Method */
bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root)
{
    double x = c;

    for (int i = 0; i < MAX_ITER; ++i)
    {
        double fx = f(x);
        double gx = g(x);

        if (std::abs(gx) < 1e-12)
            return false;

        double x_new = x - fx / gx;

        if (x_new < a || x_new > b)
            return false;

        if (std::abs(x_new - x) < TOL)
        {
            *root = x_new;
            return true;
        }

        x = x_new;
    }

    return false;
}

/* Secant Method */
bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root)
{
    double x0 = a;
    double x1 = c;

    double f0 = f(x0);
    double f1 = f(x1);

    for (int i = 0; i < MAX_ITER; ++i)
    {
        if (std::abs(f1 - f0) < 1e-12)
            return false;

        double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

        if (x2 < a || x2 > b)
            return false;

        if (std::abs(x2 - x1) < TOL)
        {
            *root = x2;
            return true;
        }

        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f(x1);
    }

    return false;
}
