﻿// Runge-KuttaMethod.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
double sigma = 1; //глобальный параметр
#pragma  region testFunctions
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
#pragma  endregion testFunctions
#pragma  region functionsZhukovskyProblem
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
#pragma  endregion functionsZhukovskyProblem
//void RungeKutta(double (*f1)(double t, double y1, double y2, double y3, double y4), double (*f2)(double t, double y1, double y2, double y3, double y4), double (*f3)(double t, double y1, double y2, double y3, double y4), double (*f4)(double t, double y1, double y2, double y3, double y4), double* y1, double* y2, double* y3, double* y4, double h, int n, double* t, double y1_0, double y2_0, double y3_0, double y4_0) {
//    y1[0] = y1_0;
//    y2[0] = y2_0;
//    y3[0] = y3_0;
//    y4[0] = y4_0;
//    double k11, k12, k13, k14, k21, k22, k23, k24;
//    for (int i = 0; i < n; i++) {
//        k11 = f1(t[i], y1[i], y2[i], y3[i], y4[i]);
//        k12 = f2(t[i], y1[i], y2[i], y3[i], y4[i]);
//        k13 = f3(t[i], y1[i], y2[i], y3[i], y4[i]);
//        k14 = f4(t[i], y1[i], y2[i], y3[i], y4[i]);
//
//        k21 = f1(t[i] + h, y1[i] + h * k11, y2[i] + h * k12, y3[i] + h * k13, y4[i] + h * k14);
//        k22 = f2(t[i] + h, y1[i] + h * k11, y2[i] + h * k12, y3[i] + h * k13, y4[i] + h * k14);
//        k23 = f3(t[i] + h, y1[i] + h * k11, y2[i] + h * k12, y3[i] + h * k13, y4[i] + h * k14);
//        k24 = f4(t[i] + h, y1[i] + h * k11, y2[i] + h * k12, y3[i] + h * k13, y4[i] + h * k14);
//
//
//        y1[i + 1] = y1[i] + h * (k11 + k21) / 2;
//        y2[i + 1] = y2[i] + h * (k12 + k22) / 2;
//        y3[i + 1] = y3[i] + h * (k13 + k23) / 2;
//        y4[i + 1] = y4[i] + h * (k14 + k24) / 2;
//    }
//}

/// <summary>
/// Метод Рунге-Кутта принимает количество уравнений в системе, набор функций для каждого уравнения, двумерный массив для решений, где первый индекс номер состояния, второй номер переменной в системе, шаг метода, количество разбиений, набор значений, по которому производится интегрирование, набор начальных условий
/// </summary>
void RungeKutta(int count, double (**f)(double t, double* y), double** y, double h, int n, double* t, double* y0) {
    for (int i = 0; i < count; i++) {
        y[0][i] = y0[i]; //задание начальных условий
    }
    double* k1 = new double[count]();
    double* k2 = new double[count]();
    for (int i = 0; i < n-1; i++) {
        for (int j = 0; j < count; j++) { //вычисление k1 для всех функций
           k1[j] = f[j](t[i], y[i]);
        }
        double* newY = new double[count]; //определение новой точки для вычисления k2
        for (int j = 0; j < count; j++) {
            newY[j] = y[i][j] + h * k1[j];
        }
        for (int j = 0; j < count; j++) { //вычисление k2 для всех функций
            k2[j] = f[j](t[i]+h, newY);
        }
        for (int j = 0; j < count; j++) { //вычисление значений в следующей точке
            y[i + 1][j] = y[i][j] + h * (k1[j] + k2[j]) / 2;
        }
    }
}
/// <summary>
/// Вычисление погрешности
/// </summary>
double Inaccuracy(int numberDivisions, double* t, double** y) {
    double res = 0;
    double currentInaccuracy;
    for (int i = 0; i < numberDivisions; i++) {
        double testy1 = testY1(t[i]);
        double testy2 = testY2(t[i]);
        currentInaccuracy = y[i][0] - testy1;
        double otherInaccuracy = y[i][1] - testy2;
        if (otherInaccuracy > currentInaccuracy) currentInaccuracy = otherInaccuracy;
        if (currentInaccuracy > res) res = currentInaccuracy;
        //std::cout << y[i][0] << " " << testy1 << " " << y[i][1] << " " << testy2 << " " << y[i][0] - testy1 << " " << y[i][1] - testy2 << "\n";
    }
    return res;
}
typedef double (*pFunc)(double t, double* y);
/// <summary>
/// Решение тестовой задачи и возвращение пары (h,погрешность)
/// </summary>
std::pair<double, double> testTask(int numberDivisions) {
    double h = 2.0 / numberDivisions; //задание шага метода
    double* y1 = new double[numberDivisions]();
    double* y2 = new double[numberDivisions]();
    double* t = new double[numberDivisions];
    for (int i = 0; i < numberDivisions; i++) {
        t[i] = h * i; //создание оси t
    }
    double (**f)(double t, double* y) = new pFunc[2]; //объявление набора функций
    f[0] = testF1;
    f[1] = testF2;
    double** y = new double* [numberDivisions]; //выделение памяти для решения системы
    for (int i = 0; i < numberDivisions; i++) {
        y[i] = new double[2]();
    }
    double* y0 = new (double[2]); //начальные условия
    y0[0] = testY1(t[0]);
    y0[1] = testY2(t[0]);
    RungeKutta(2, f, y, h, numberDivisions, t, y0); 
    std::pair<double, double> res(h, Inaccuracy(numberDivisions, t, y));
    return res;
}
int main()
{
    std::cout << "\"X\", \"Y\"" << "\n";
    /*for (int i = 1000; i > 100; i--) {
        std::pair<double, double> inaccuracy = testTask(i);
        std::cout << inaccuracy.first<<","<< inaccuracy.second << "\n";
    }
    */
    int n = 1000;
    double h = 20.0 / n;
    double* t = new double[n];
    for (int i = 0; i < n; i++) {
        t[i] = h * i;
    }
    double (**f)(double t, double* y) = new pFunc[4];
    f[0] = V;
    f[1] = Teta;
    f[2] = Z;
    f[3] = X;
    double** y = new double* [n];
    for (int i = 0; i < n; i++) {
        y[i] = new double[4]();
    }
    double* y0 = new (double[4]);
    y0[0] =10 ;
    y0[1] =30;
    y0[2] = 100;
    y0[3] = 0; //x=0
    
    RungeKutta(4, f, y, h, n, t, y0);
    for (int i = 0; i < n; i++) 
        std::cout << y[i][3] << "," << y[i][2] << "\n";
}