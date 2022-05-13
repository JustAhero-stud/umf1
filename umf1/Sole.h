#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include"Grid.h"
#include <iomanip>
using namespace std;



class SOLE {
private:
    vector <double> diag, ld1, ld2, ld3, rd1, rd2, rd3;   //последовательность с низу вверх -> ld2 ld1 diag rd1 rd2!!!!!
    vector <double> f;
    vector <double> u;
    int n, n_d1, n_d2, n_d3, offset1, offset2;
    int bc[6];
    double func_f(double x, double y)
    {
        //return 20 * (x + y);
        return 20 * (x * x + y * y) - 4;//квадрат
        //return -6 * (x + y) + 20 * (x * x * x + y * y * y);
        //return -12 * (x + y) + 20 * (x * x * x * x + y * y * y * y);
        //return 9* cos(2*x + 2*y);
    }

    double func_ug(int side, double x, double y)
    {
        int d = 2;
        switch (side)
        {
            /*case 1:
                return cos(2*x + 2*0.2);
            case 2:
                return cos(2 * 0.4 + 2 * y);
            case 3:
                return cos(2 * x + 2 * 0.5);
            case 4:
                return cos(2 * 0.7 + 2 * y);
            case 5:
                return cos(2 * x + 2 * 0.8);
            case 6:
                return cos(2 * 0.1 + 2 * y);*/
        case 1:
            return pow(2, d) + pow(x, d);
        case 2:
            return pow(4, d) + pow(y, d);
        case 3:
            return pow(x, d) + pow(5, d);
        case 4:
            return pow(y, d) + pow(7, d);
        case 5:
            return pow(x, d) + pow(8, d);
        case 6:
            return 1 + pow(y, d);
        }
    }

    double theta(int side, double x, double y)
    {
        switch (side)
        {
        case 1:
            return -1;
        case 2:
            return 1;
        case 3:
            return -1;
        case 4:
            return 1;
        case 5:
            return 1;
        case 6:
            return -1;
        }
    }

    //  Инициализация СЛАУ путём наложения на все узлы нулевого первого краевого условия
    void init_SOLE(Grid& grid)
    {
        n = grid.getN();
        offset1 = grid.getNx(1);
        offset2 = grid.getNx(1) + grid.getNx(2) - 1;
        n_d1 = n - 1;
        n_d2 = n - (offset1);
        n_d3 = n - (offset2);
        diag.resize(n, 1);
        ld1.resize(n_d1);
        rd1.resize(n_d1);
        ld2.resize(n_d2);
        rd2.resize(n_d2);
        ld3.resize(n_d3);
        rd3.resize(n_d3);
        f.resize(n);
        u.resize(n);
    }

    void calc_diags(Grid grid, int i, int j, bool flag)
    {
        static double lamb = grid.getLambda();
        static double gamma = grid.getGamma();
        double x0, x1, x2, y0, y1, y2, hx0, hx1, hy0, hy1;
        int m;
        x0 = grid.getX(i - 1);
        x1 = grid.getX(i);
        x2 = grid.getX(i + 1);
        hx0 = x1 - x0;
        hx1 = x2 - x1;

        y0 = grid.getY(j - 1);
        y1 = grid.getY(j);
        y2 = grid.getY(j + 1);
        hy0 = y1 - y0;
        hy1 = y2 - y1;
        m = grid.calc_num(i, j, flag);
        diag[m] = lamb * 2 * (1 / (hx1 * hx0) + 1 / (hy1 * hy0)) + gamma;//главная диагноаль





        if (m >= n_d2 && m < n_d1) {//++
            rd1[m] = -2 * lamb / (hx1 * (hx0 + hx1));
            ld1[m - 1] = -2 * lamb / (hx0 * (hx0 + hx1));
            ld3[m - (offset2 + 2)] = -2 * lamb / (hy0 * (hy0 + hy1));
        }
        else if (m == n_d1)
        {
            ld1[m - 1] = -2 * lamb / (hx0 * (hx0 + hx1));
            ld3[m - offset2] = -2 * lamb / (hy0 * (hy0 + hy1));
        }
        else if (m < offset1 && m > 0) {
            ld1[m - 1] = -2 * lamb / (hx0 * (hx0 + hx1));
            rd1[m] = -2 * lamb / (hx1 * (hx0 + hx1));
            rd2[m] = -2 * lamb / (hy1 * (hy0 + hy1));
        }
        else if (m == 0) {//++
            rd1[m] = -2 * lamb / (hx1 * (hx0 + hx1));
            rd2[m] = -2 * lamb / (hy1 * (hy0 + hy1));
        }
        else {//для 4 элементов 4 случая
            if (m < offset1 * (grid.getNy(1) - 1))
            {
                rd1[m] = -2 * lamb / (hx1 * (hx0 + hx1));
                rd2[m] = -2 * lamb / (hy1 * (hy0 + hy1));
                ld1[m - 1] = -2 * lamb / (hx0 * (hx0 + hx1));
                ld2[m - offset1] = -2 * lamb / (hy0 * (hy0 + hy1));
            }
            else
            {
                if (m >= offset1 * (grid.getNy(1) - 1) && m < offset1 * (grid.getNy(1) + 1))
                {
                    rd1[m] = -2 * lamb / (hx1 * (hx0 + hx1));
                    rd3[m] = -2 * lamb / (hy1 * (hy0 + hy1));
                    ld1[m - 1] = -2 * lamb / (hx0 * (hx0 + hx1));
                    ld2[m - offset1] = -2 * lamb / (hy0 * (hy0 + hy1));
                }
                else
                {
                    if (m >= offset1 * (grid.getNy(1) + 1) && m < offset1 * (grid.getNy(1) - 1) + offset2)//надо ли вообще???
                    {
                        rd1[m] = -2 * lamb / (hx1 * (hx0 + hx1));
                        rd3[m] = -2 * lamb / (hy1 * (hy0 + hy1));
                    }
                    else
                    {
                        rd1[m] = -2 * lamb / (hx1 * (hx0 + hx1));
                        rd3[m] = -2 * lamb / (hy1 * (hy0 + hy1));
                        ld1[m - 1] = -2 * lamb / (hx0 * (hx0 + hx1));
                        ld3[m - offset2] = -2 * lamb / (hy0 * (hy0 + hy1));
                    }
                }
            }

        }

        f[m] = func_f(x1, y1);
    }

    // Для внутренних узлов
    void make_matrix(Grid& grid)
    {
        int nx1 = grid.getNx(1);
        int nx2 = grid.getNx(2);
        int ny1 = grid.getNy(1);
        int ny2 = grid.getNy(2);

        //--------------------------------- 1 -------------------
        //--------------------------------- 1 -------------------
        for (int i = 1; i < nx1 - 1; i++)
            for (int j = 1; j < ny1 - 1; j++)
                calc_diags(grid, i, j, true);

        //--------------------------------- 2 -------------------
        //--------------------------------- 2 -------------------
        for (int i = 1; i < nx1 - 1; i++)
            for (int j = ny1 - 1; j < ny1 + ny2 - 2; j++)
                calc_diags(grid, i, j, false);

        //--------------------------------- 3 -------------------
        //--------------------------------- 3 -------------------
        for (int i = nx1 - 1; i < nx1 + nx2 - 2; i++)
            for (int j = ny1; j < ny1 + ny2 - 2; j++)
                calc_diags(grid, i, j, false);
    }

    void boundary_conditions(Grid& grid)
    {
        // Чтение краевых условий из файла
        ifstream bc_file;
        bc_file.open("bc.txt");
        for (int i = 0; i < 6; i++)
            bc_file >> bc[i];
        bc_file.close();

        // Учитывание краевых условий
        for (int i = 0; i < 6; i++)
            if (bc[i] == 1)
                bc1(i + 1, grid);
            else
                bc2(i + 1, grid);
    }

    void bc1(int side, Grid& grid)
    {
        int nx1 = grid.getNx(1);
        int nx2 = grid.getNx(2);
        int ny1 = grid.getNy(1);
        int ny2 = grid.getNy(2);
        int i, j, m;
        double x, y;
        switch (side)
        {
        case 1:
            j = 0;
            y = grid.getY(j);
            for (int i = 0; i < nx1; i++)
            {
                m = grid.calc_num(i, j, true);
                x = grid.getX(i);
                f[m] = func_ug(side, x, y);
            }
            break;

        case 2:
            i = nx1 - 1;
            x = grid.getX(i);
            for (int j = 0; j < ny1; j++)
            {
                m = grid.calc_num(i, j, true);
                y = grid.getY(j);
                f[m] = func_ug(side, x, y);
            }
            break;

        case 3:
            j = ny1 - 1;
            y = grid.getY(j);
            for (int i = nx1 - 1; i < nx1 + nx2 - 1; i++)
            {
                m = grid.calc_num(i, j, false);
                x = grid.getX(i);
                f[m] = func_ug(side, x, y);
            }
            break;

        case 4:
            i = nx1 + nx2 - 2;
            x = grid.getX(i);
            for (int j = ny1 - 1; j < ny1 + ny2 - 1; j++)
            {
                m = grid.calc_num(i, j, false);
                y = grid.getY(j);
                f[m] = func_ug(side, x, y);
            }
            break;

        case 5:
            j = ny1 + ny2 - 2;
            y = grid.getY(j);
            for (int i = 0; i < nx1 + nx2 - 1; i++)
            {
                m = grid.calc_num(i, j, false);
                x = grid.getX(i);
                f[m] = func_ug(side, x, y);
            }
            break;

        case 6:
            i = 0;
            x = grid.getX(i);
            for (int j = 0; j < ny1; j++)
            {
                m = grid.calc_num(i, j, true);
                y = grid.getY(j);
                f[m] = func_ug(side, x, y);
            }
            for (int j = ny1; j < ny1 + ny2 - 1; j++)
            {
                m = grid.calc_num(i, j, false);
                y = grid.getY(j);
                f[m] = func_ug(side, x, y);
            }
            break;
        }
    }

    void bc2(int side, Grid& grid)
    {
        int nx1 = grid.getNx(1);
        int nx2 = grid.getNx(2);
        int ny1 = grid.getNy(1);
        int ny2 = grid.getNy(2);
        double lamb = grid.getLambda();
        int i, j, m;
        double x0, y0, x1, y1;
        switch (side)
        {
        case 1:
            j = 0;
            y0 = grid.getY(j);
            for (int i = 1; i < nx1 - 1; i++)
            {
                m = grid.calc_num(i, j, true);
                x0 = grid.getX(i);
                x1 = grid.getX(i + 1);
                diag[m] = lamb / (x1 - x0);
                rd2[m] = -lamb / (x1 - x0);
                f[m] = theta(side, x0, y0);
            }
            break;

        case 2:
            i = nx1 - 1;
            x0 = grid.getX(i);
            for (int j = 1; j < ny1 - 1; j++)
            {
                m = grid.calc_num(i, j, true);
                y0 = grid.getY(j);
                y1 = grid.getY(j - 1);
                diag[m] = lamb / (y0 - y1);
                ld1[m - 1] = -lamb / (y0 - y1);
                f[m] = theta(side, x0, y0);
            }
            break;

        case 3:
            j = ny1 - 1;
            y0 = grid.getY(j);
            for (int i = nx1; i < nx1 + nx2 - 2; i++)
            {
                m = grid.calc_num(i, j, false);
                x0 = grid.getX(i);
                x1 = grid.getX(i + 1);
                diag[m] = lamb / (x1 - x0);
                rd3[m] = -lamb / (x1 - x0);
                f[m] = theta(side, x0, y0);
            }
            break;

        case 4:
            i = nx1 + nx2 - 2;
            x0 = grid.getX(i);
            for (int j = ny1; j < ny1 + ny2 - 2; j++)
            {
                m = grid.calc_num(i, j, false);
                y0 = grid.getY(j);
                y1 = grid.getY(j - 1);
                diag[m] = lamb / (y0 - y1);
                ld1[m - 1] = -lamb / (y0 - y1);
                f[m] = theta(side, x0, y0);
            }
            break;

        case 5:
            j = ny1 + ny2 - 2;
            y0 = grid.getY(j);
            for (int i = 1; i < nx1 + nx2 - 2; i++)
            {
                m = grid.calc_num(i, j, false);
                x0 = grid.getX(i);
                x1 = grid.getX(i - 1);
                diag[m] = lamb / (x1 - x0);
                ld3[m - offset2] = -lamb / (x1 - x0);
                f[m] = theta(side, x0, y0);
            }
            break;

        case 6:
            i = 0;
            x0 = grid.getX(i);
            for (int j = 1; j < ny1; j++)
            {
                m = grid.calc_num(i, j, true);
                y0 = grid.getY(j);
                y1 = grid.getY(j + 1);
                diag[m] = lamb / (y1 - y0);
                rd1[m] = -lamb / (y1 - y0);
                f[m] = theta(side, x0, y0);
            }
            for (int j = ny1; j < ny1 + ny2 - 2; j++)
            {
                m = grid.calc_num(i, j, false);
                y0 = grid.getY(j);
                y1 = grid.getY(j + 1);
                diag[m] = lamb / (y1 - y0);
                rd1[m] = -lamb / (y1 - y0);
                f[m] = theta(side, x0, y0);
            }
            break;
        }
    }

    // Вычисление нормы вектора
    double vector_norm(vector <double>& v) {
        double norm = 0;
        for (int i = 0; i < n; i++)
            norm += v[i] * v[i];
        return sqrt(norm);
    }

    // Выполнение итерационного шага
    double iterative_step(vector <double>& v, int i, double& ri2) {
        int m_offset = offset1;
        double sum = 0;
        double result;
        if (i < n_d3)
            sum += rd3[i] * v[offset2 + i];
        if (i < n_d2)
            sum += rd2[i] * v[m_offset + i];

        if (i < n_d1)
            sum += rd1[i] * v[1 + i];

        sum += diag[i] * v[i];

        if (i >= 1)
            sum += ld1[i - 1] * v[i - 1];

        if (i >= m_offset)
            sum += ld2[i - m_offset] * v[i - m_offset];
        if (i >= offset2)
            sum += ld3[i - offset2] * v[i - offset2];
        double ri = f[i] - sum;
        ri2 = ri * ri;
        result = ri / diag[i];
        return result;
    }

    // Метод Гаусса — Зейделя
    void gauss_seidel(double eps = 1.0E-14, int maxiter = 1000) {
        static double fnorm = vector_norm(f);
        double r = 2;  // относительная невязка
        double ri2;
        double ri2_sum = 0;
        int iter = 1;
        for (; iter <= maxiter && r > eps; iter++) {
            for (int i = 0; i < n; i++) {
                u[i] = u[i] + iterative_step(u, i, ri2);
                ri2_sum += ri2;
            }
            r = sqrt(ri2_sum) / fnorm;
            cout << "Current iteration: " << iter << endl;
            cout << "Approximation residual: " << r << endl;
            ri2_sum = 0;
        }
        for (int i = 0; i < n; i++) cout << setprecision(16) << u[i] << endl;
        //вывод для 12
        /*for(int i=0; i<16; i+=4)   cout << setprecision(16) <<u[13*i]<< " " << u[4+ 13 * i] << " " <<u[8+ 13 * i]<< " " <<u[12+ 13 * i]<< endl;
        cout << endl;
        for(int i=0; i<16; i+=4) cout << setprecision(16) << u[13*12+25* i] << " " << u[4 + 13 * 12+ 25 * i] << " " << u[8 + 13 * 12+ 25 * i] << " " << u[12 + 13 * 12+ 25 * i]<<" " << u[16 + 13 * 12 + 25 * i] << " " << u[20 + 13 * 12 + 25 * i] << " " << u[24 + 13 * 12 + 25 * i] << endl;*/
    }

public:
    SOLE(Grid& grid)
    {
        init_SOLE(grid);
        make_matrix(grid);
        boundary_conditions(grid);
        gauss_seidel();
    }
};
