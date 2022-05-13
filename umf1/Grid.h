#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
using namespace std;

class Area
{
protected:  // Доступно не только в самом классе, но и в его наследниках
    double lambda, gamma;
    double X1, X2, X3, Y1, Y2, Y3; //параметры, описывающие расчётную область в форме "Г";
    void read_area() {  // Читает параметры области и коэффициенты
        ifstream area_file;
        area_file.open("Area.txt");
        area_file >> X1 >> X2 >> X3 >> Y1 >> Y2 >> Y3 >> lambda >> gamma;
        area_file.close();
    }

public:
    double getLambda() { return lambda; }
    double getGamma() { return gamma; }
};

// Класс сетки, которые наследует класс расчётной области
class Grid : public Area
{
private:
    vector <double> X, Y;
    int nx1, nx2;
    int ny1, ny2;
    double qx1, qx2;
    double qy1, qy2;
    int N;  // Количество узлов всего

    void read_grid(string file_name, int& n1, int& n2,
        double& q1, double& q2)
    {
        ifstream GridXY;
        GridXY.open(file_name);
        GridXY >> n1 >> n2 >> q1 >> q2;
        GridXY.close();
    }

    void makeGrid(vector <double>& XY, double left, double right, int n, double q, int i0) {
        double h0;
        if (q - 1 < 1E-16)
            h0 = (right - left) / n;
        else
            h0 = (right - left) * (1 - q) / (1 - pow(q, n));
        XY[i0] = left;
        XY[n + i0] = right;
        for (int i = i0 + 1; i < n + i0; i++) {
            XY[i] = XY[i - 1] + h0;
            h0 *= q;
        }
    };

public:
    Grid() {   // Конструктор класса
        read_area();
        read_grid("GridX.txt", nx1, nx2, qx1, qx2);
        read_grid("GridY.txt", ny1, ny2, qy1, qy2);
        X.resize(nx1 + nx2 + 1);
        Y.resize(ny1 + ny2 + 1);
        makeGrid(X, X1, X2, nx1, qx1, 0);
        makeGrid(X, X2, X3, nx2, qx2, nx1);
        makeGrid(Y, Y1, Y2, ny1, qy1, 0);
        makeGrid(Y, Y2, Y3, ny2, qy2, ny1);
        for (int i = 0; i < Y.size(); i++) cout << setprecision(16) << Y[i] << endl;
        for (int i = 0; i < X.size(); i++) cout << setprecision(16) << X[i] << endl;
        nx1 += 1; nx2 += 1;
        ny1 += 1; ny2 += 1;
        N = nx1 * ny1 + nx1 * ny2 + nx2 * ny2 - nx1 - ny2;
    }

    int calc_num(int i, int j, bool flag)
    {
        if (flag)
            return i + j * nx1;
        j = j - ny1 + 1;
        return nx1 * (ny1 - 1) + i + j * (nx1 + nx2 - 1);
    }

    double getX(int i) { return X[i]; }
    double getY(int i) { return Y[i]; }
    int getNx(int i) { if (i == 1) return nx1; return nx2; }
    int getNy(int i) { if (i == 1) return ny1; return ny2; }
    int getN() { return N; }
};

