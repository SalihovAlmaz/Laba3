
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;
const int signs = 5;

//params
//params.f - интегрируемая функция
//params.a - начало отрезка
//params.b - конец отрезка
//n - заданное разбиение
struct param
{
    double a;
    double b;
    double exactIntegral; //32/3, // Точное значение определенного интеграла x^5 на [0,2]
    double c;
    param() {

    }
    param(double a, double b, double exactIntegral = -1, double c = -1) {
        this->a = a;
        this->b = b;
        this->exactIntegral = exactIntegral;
        this->c = c;
    }
    double f(double x) {
        return 6 * (pow(x, 5));
    }
    double f1(double x) {
        return pow(x,1./15) * sqrt(1+pow(x,2));
    }
};

double numIntegral(param p, double n)
{
    const double h = (p.b - p.a) / n;


    // заполняем узловые точки
    vector<double> x;
    for (double i = 0; i <= n; i++)
    {
        x.push_back(p.a + i * h);
    }

    double sum1 = 0;
    for (double i = 1; i < n; i++)
    {
        sum1 += p.f(x[i]);
    }
    double sum2 = 0;
    for (double i = 0; i < n; i++)
    {
        sum2 += p.f(x[i] + h / 2);
    }

    // Квадратчная формула Ньютона-Котеса
    // метод парабол
    return (h / 6) * (p.f(p.a) + p.f(p.b) + 2 * sum1 + 4 * sum2);
}
double numIntegral1(param p, double n)
{
    const double h = (p.b - p.a) / n;


    // заполняем узловые точки
    vector<double> x;
    for (double i = 0; i <= n; i++)
    {
        x.push_back(p.a + i * h);
    }

    double sum1 = 0;
    for (double i = 1; i < n; i++)
    {
        sum1 += p.f1(x[i]);
    }
    double sum2 = 0;
    for (double i = 0; i < n; i++)
    {
        sum2 += p.f1(x[i] + h / 2);
    }

    // Квадратчная формула Ньютона-Котеса
    // метод парабол
    return (h / 6) * (p.f1(p.a) + p.f1(p.b) + 2 * sum1 + 4 * sum2);
}
/*
Получить результаты приблеженного счёта со всеми значениями,
начиная с n=2, и их погрешность по методу парабол
params
params.f - Заданная функция
params.a - Начало отрезка
params.b - Конец отрезка
params.exactIntegral - точное значение интеграла
params.с - с
e - Заданная точность
*/


struct row {
    double n;
    double integral;
    double deltaK;
    double deltaRunge;
    double deltaExact;
    double deltaTheory;
    row() {}
    row(double n, double integral, double deltaK, double deltaRunge) {
        this->n = n;
        this->integral = integral;
        this->deltaK = deltaK;
        this->deltaRunge = deltaRunge;

    }
    void setdeltaExact(double deltaExact) {
        this->deltaExact = deltaExact;
    }
    void setdeltaTheory(double deltaTheory) {
        this->deltaTheory = deltaTheory;
    }
};


vector<row> Process(param p, double e) {
    vector<row> result;
    double runge = -1;
    double integral, integralH2, integralH4, integralDiv2, deltaK;

    // Порядок точночти метода парабол 
    const double k = 4;

    for (double n = 2; abs(runge) > e; n *= 2)
    {
        integral = numIntegral(p, n);
        integralH2 = numIntegral(p, n * 2);
        integralH4 = numIntegral(p, n * 4);
        integralDiv2 = numIntegral(p, n / 2);
        runge = (integral - integralDiv2) / (pow(2, k) - 1);
        deltaK = (integralH2 - integral) / (integralH4 - integralH2);

        row r(n, integral, deltaK, runge);

        if (p.exactIntegral != -1) {
            r.setdeltaExact(abs(p.exactIntegral - integral));

        }
        if (p.c != -1) {
            double h = (p.b - p.a) / n;
            r.setdeltaTheory(abs(p.c * pow(h, k)));
        }
        result.push_back(r);
    }
    return result;
}

vector<row> Process1(param p, double e) {
    vector<row> result;
    double runge = -1;
    double integral, integralH2, integralH4, integralDiv2, deltaK;

    // Порядок точночти метода парабол 
    const double k = 4;

    for (double n = 2; abs(runge) > e; n*=2)
    {
        integral = numIntegral1(p,n);
        integralH2 = numIntegral1(p,n*2);
        integralH4 = numIntegral1(p, n * 4);
        integralDiv2 = numIntegral1(p, n / 2);
        runge = (integral - integralDiv2) / (pow(2, k) - 1);
        deltaK = (integralH2 - integral) / (integralH4 - integralH2);

        row r(n, integral, deltaK, runge);

        if (p.exactIntegral != -1) {
            r.setdeltaExact(abs(p.exactIntegral - integral));
       
        }
        if (p.c != -1) {
            double h = (p.b - p.a) / n;
            r.setdeltaTheory(abs(p.c * pow(h,k)));
        }
        result.push_back(r);
    }
    return result;
}



int main()
{
    // заданная точность
    const double e = 1e-7;
    //отладочный пример
    const param test(0, 1, 1, double(720) / 2880);

    const vector<row> resultTest = Process(test, e);

    const double integralTest = resultTest[resultTest.size() - 1].integral;

    cout << "index\t" << "n\t" << "integral\t" << "deltaK\t\t" << "deltRunge\t\t" << "deltaExact\t\t" << "deltaTheory\n";
   for(size_t i = 0; i < resultTest.size(); i++)
    {
            cout << i << "\t" << resultTest[i].n << "\t" << resultTest[i].integral << "\t\t" << resultTest[i].deltaK << "\t\t" << resultTest[i].deltaRunge << "\t\t" << resultTest[i].deltaExact << "\t\t" << resultTest[i].deltaTheory << "\n";
    }
   setlocale(LC_ALL, "Russian");
    cout << endl <<"Значение интеграла 6x^5 на отрезке [0, 1] = " << integralTest << endl << endl;
   const param test1(0, 1.5);

   const vector<row> resultTest1 = Process1(test1, e);

   const double integralTest1 = resultTest1[resultTest1.size() - 1].integral;

   cout << "index\t" << "n\t" << "integral\t" << "deltaK\t\t" << "deltRunge\n";
   for (size_t i = 0; i < resultTest1.size(); i++)
   {
       cout << i << "\t" << resultTest1[i].n << "\t" << resultTest1[i].integral << "\t\t" << resultTest1[i].deltaK << "\t\t" << resultTest1[i].deltaRunge << "\n";
   }
   cout << endl << "Значение интеграла x^(1/15) * sqrt(1+x^2) на отрезке [0, 1.5] = " << integralTest1 << endl << endl;

   const param test2(0.001,1.5);

   const vector<row> resultTest2 = Process1(test2, e);

   const double integralTest2 = resultTest2[resultTest2.size() - 1].integral;

   cout << "index\t" << "n\t" << "integral\t" << "deltaK\t\t" << "deltRunge\n";
   for (size_t i = 0; i < resultTest2.size(); i++)
   {
       cout << i << "\t" << resultTest2[i].n << "\t" << resultTest2[i].integral << "\t\t" << resultTest2[i].deltaK << "\t\t" << resultTest2[i].deltaRunge << "\n";
   }
   cout << endl << "Значение интеграла x^(1/15) * sqrt(1+x^2) на отрезке [0.001, 1.5] = " << integralTest2 << endl << endl;

}
 
