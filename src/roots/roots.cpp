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
bool secant(std::function<double(double)> f, double a, double b, double c, double *root)
{
    double x_old = c + TOL;  // first initial guess (slightly offset)
    double x_current = c;     // second initial guess

    // Check if the root is guaranteed to exist in [a,b]
    if (f(a) * f(b) > 0) {
        return false;
    }

    for (int i = 0; i < MAX_ITER; ++i)
    {
        double f_current = f(x_current);
        double f_old = f(x_old);

        // Avoid division by zero
        if (std::fabs(f_current - f_old) < 1e-12) {
            return false;
        }

        // Secant formula
        double x_new = x_current - f_current * ((x_current - x_old) / (f_current - f_old));

        // Check if new guess is within the interval
        if (x_new < a || x_new > b) {
            return false;
        }

        // Check if we are close enough to the root
        if (std::fabs(f(x_new)) < TOL) {
            *root = x_new;
            return true;
        }

        // Update points for next iteration
        x_old = x_current;
        x_current = x_new;
    }

    return false; // Root not found within MAX_ITER
}