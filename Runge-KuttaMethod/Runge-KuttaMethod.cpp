// Runge-KuttaMethod.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
double sigma = 1;
double testY1(double t) {
    return cos(t * t)* sqrt(1 + t);
}
double testY2(double t) {
    return sin(t * t) * sqrt(1 + t);
}
double TestF1(double t, double y1, double y2) {
    return y1 / (2 + 2 * t) - 2 * t * y2;
}
double TestF2(double t, double y1, double y2) {
    return y2 / (2 + 2 * t) + 2 * t * y1;
}
double testF1(double t, double* y) {
    return TestF1(t, y[0], y[1]);
}
double testF2(double t, double* y) {
    return TestF2(t, y[0], y[1]);
}
double V(double t, double* y) {
    return -sin(y[1]) - sigma * y[0] * y[0];
}
double Teta(double t, double* y) {
    return (y[0] * y[0] - cos(y[1])) / y[0];
}
double Z(double t, double* y) {
    return y[0] * sin(y[1]);
}
double X(double t, double* y) {
    return y[0] * cos(y[1]);
}
void RungeKutta(double (*f1)(double t, double y1, double y2, double y3, double y4), double (*f2)(double t, double y1, double y2, double y3, double y4), double (*f3)(double t, double y1, double y2, double y3, double y4), double (*f4)(double t, double y1, double y2, double y3, double y4), double* y1, double* y2, double* y3, double* y4, double h, int n, double* t, double y1_0, double y2_0, double y3_0, double y4_0) {
    y1[0] = y1_0;
    y2[0] = y2_0;
    y3[0] = y3_0;
    y4[0] = y4_0;
    double k11, k12, k13, k14, k21, k22, k23, k24;
    for (int i = 0; i < n; i++) {
        k11 = f1(t[i], y1[i], y2[i], y3[i], y4[i]);
        k12 = f2(t[i], y1[i], y2[i], y3[i], y4[i]);
        k13 = f3(t[i], y1[i], y2[i], y3[i], y4[i]);
        k14 = f4(t[i], y1[i], y2[i], y3[i], y4[i]);

        k21 = f1(t[i] + h, y1[i] + h * k11, y2[i] + h * k12, y3[i] + h * k13, y4[i] + h * k14);
        k22 = f2(t[i] + h, y1[i] + h * k11, y2[i] + h * k12, y3[i] + h * k13, y4[i] + h * k14);
        k23 = f3(t[i] + h, y1[i] + h * k11, y2[i] + h * k12, y3[i] + h * k13, y4[i] + h * k14);
        k24 = f4(t[i] + h, y1[i] + h * k11, y2[i] + h * k12, y3[i] + h * k13, y4[i] + h * k14);


        y1[i + 1] = y1[i] + h * (k11 + k21) / 2;
        y2[i + 1] = y2[i] + h * (k12 + k22) / 2;
        y3[i + 1] = y3[i] + h * (k13 + k23) / 2;
        y4[i + 1] = y4[i] + h * (k14 + k24) / 2;
    }
}

void RungeKutta(int count, double (**f)(double t, double* y), double** y, double h, int n, double* t, double* y0) {
    for (int i = 0; i < count; i++) {
        y[0][i] = y0[i];
    }
    double* k1 = new double[count]();
    double* k2 = new double[count]();
    for (int i = 0; i < n-1; i++) {
        for (int j = 0; j < count; j++) {
           k1[j] = f[j](t[i], y[i]);
        }
        double* newY = new double[count];
        for (int j = 0; j < count; j++) {
            newY[j] = y[i][j] + h * k1[j];
        }
        for (int j = 0; j < count; j++) {
            k2[j] = f[j](t[i]+h, newY);
        }
        for (int j = 0; j < count; j++) {
            y[i + 1][j] = y[i][j] + h * (k1[j] + k2[j]) / 2;
        }
    }
}
typedef double (*pFunc)(double t, double* y);
int main()
{
    
    int n = 1000;
    double h = 2.0 / n;
    double* y1 = new double[n]();
    double* y2 = new double[n]();
    double* t = new double[n];
    for (int i = 0; i < n; i++) {
        t[i] = h * i;
    }
    double (**f)(double t, double* y) = new pFunc[2];
    f[0] = testF1;
    f[1] = testF2;
    double** y = new double*[n];
    for (int i = 0; i < n; i++) {
        y[i] = new double[2]();
    }
    double* y0 = new (double[2]);
    y0[0] = testY1(t[0]);
    y0[1] = testY2(t[0]);
    std::cout << y0[0] << " " << y0[1] << "\n";
    //RungeKutta(f, y, h, n, t, y0);
    RungeKutta(2, f, y, h, n, t, y0);
    for (int i = 0; i < n; i++) {
        double testy1 = testY1(t[i]);
        double testy2 = testY2(t[i]);
        std::cout << y[i][0] << " " << testy1 << " " << y[i][1] << " " << testy2 << " " << y[i][0]-testy1<<" "<< y[i][1]-testy2 <<"\n";
    }

}