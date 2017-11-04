using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using static System.Math;
using static ЧМ_Лаб__5.LASSolver;

namespace ЧМ_Лаб__5
{
    class Program
    {
        public delegate double func(double x);
        const int n = 5, var = 16;
        const double a = 1, b = 2;
        static double[,] table = new double[2, n + 1];
        static StreamWriter outFile = new StreamWriter("out.txt");

        static double Func(double x)
        {
            switch (var)
            {
                case 15: return Pow(x, 2);
                case 16: return Atan(x) + 1 / x / x;
                case 17: return Exp(x) + x + 1;
                default: return 0;
            }
        }
        static void PrintFunc()
        {
            string str;
            switch (var)
            {
                case 15: str = "x^2"; break;
                case 16: str = "Atan(x) + 1 / x^2"; break;
                case 17: str = "Exp(x) + x + 1"; break;
                default: str = "0"; break;
            }
            PrintStr("Рассматриваемая функция: " + str);
            PrintStr();
        }

        static void FillTable()
        {
            double x = a;
            double delta = (b - a) / n;
            for(int i = 0; i <= n; i++)
            {
                table[0, i] = x;
                table[1, i] = Func(x);
                x = Round(x + delta, 1);
            }
        }
        static void PrintTable()
        {
            PrintStr("Таблица значений на отрезке [" + a + ";" + b + "] с шагом delta = " + ((b - a) / n));
            string xStr = "x:", yStr = "y:";
            for (int i = 0; i <= n; i++)
            {
                xStr += "\t" + table[0, i].ToString("0.0");
                yStr += "\t" + table[1, i].ToString("0.00000");
            }
            PrintStr(xStr);
            PrintStr(yStr);
            PrintStr();
        }

        static double OmegaFunc(double x, int num)
        {
            double res = 1;
            for(int i = 0; i < num; i++)
            {
                res *= x - table[0, i];
            }
            return res;
        }

        static double CalcNeigbours(int a, int b)
        {
            if (a == b)
            {
                return table[1, a];
            }
            else
            {
                return (CalcNeigbours(a + 1, b) - CalcNeigbours(a, b - 1)) / (table[0, b] - table[0, a]);
            }
        }

        static double InterPolynomial(double x)
        {
            double sum = 0;
            for(int i = 0; i <= n; i++)
            {
                sum += CalcNeigbours(0, i) * OmegaFunc(x, i);
            }
            return sum;
        }

        static void NyutonInterpolation()
        {
            double recived = 0, expected = 0;
            PrintStr("x     Ожидаемое  Полученное   Погрешность");
            for (double i = 1.1; i < 2; i += 0.2)
            {
                recived = InterPolynomial(i);
                expected = Func(i);
                PrintStr(i.ToString() + "   " + expected.ToString("0.000000") + "   " + recived.ToString("0.000000") + "     " + CalcError(expected, recived));
            }
            PrintStr();
        }

        static double CalcIntegral(func f, double a, double b, int aim = 10000)
        {
            double sum = 0, delta = (b - a) / aim;
            for (double i = a; i < b; i+= delta)
            {
                sum += f(i);
            }
            return Round(sum * delta, aim.ToString().Length - 1);
        }

        static func G(int num) => delegate (double x) { return Pow(x, num); };

        static func CoefFunc(func f, double k) => delegate (double x) { return k * f(x); };
        static func SumFunc(func f, func g) => delegate (double x) { return f(x) + g(x); };
        static func ProdFunc(func f, func g) => delegate (double x) { return f(x) * g(x); };

        static void СontinuousАpprox()
        {
            int n = 3, aim = 100000;
            double[,] matrix = new double[n, n];
            double[] vector = new double[n], c = new double[n];
            for(int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    matrix[i, j] = CalcIntegral(ProdFunc(G(i), G(j)), a, b, aim);
                }
                vector[i] = CalcIntegral(ProdFunc(Func, G(i)), a, b, aim);   
            }
            conjurateGradientsMethod(matrix, vector, n, ref c);
            string str = "Приблеженная функция имеет вид: ";
            for(int i = 0; i < n; i++)
            {
                if (i.Equals(0))
                {
                    str += (c[i].ToString("0.000000"));
                }
                else
                {
                    str += (Abs(c[i]).ToString("0.000000") + (i > 1 ? ("x^" + i) : "x"));
                }
                if (!i.Equals(n - 1))
                {
                    str += (c[i + 1] > 0 ? " + " : " - ");
                }
            }
            PrintStr(str);
            func sum = delegate (double x) { return 0; } ;
            for (int i = 0; i < n; i++)
            {
                sum = SumFunc(CoefFunc(G(i), c[i]), sum);
            }
            var f = CalcIntegral(ProdFunc(Func, Func), a, b, aim);
            var g = CalcIntegral(ProdFunc(sum, sum), a, b, aim);
            PrintStr("Погрешность приближения: " + Abs(f - g).ToString());
            PrintStr();
        }
        
        static void TableApprox()
        {
            int n = 3;
            double[,] matrix = new double[n, n];
            double[] vector = new double[n], c = new double[n];
            double sum = 0.0;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    sum = 0.0;
                    for (int k = 0; k <= Program.n; k++)
                        sum += G(i)(table[0, k]) * G(j)(table[0, k]);
                    matrix[i, j] = sum;
                }
                sum = 0.0;
                for (int k = 0; k <= Program.n; k++)
                    sum += table[1, k] * G(i)(table[0, k]);
                vector[i] = sum;
            }
            conjurateGradientsMethod(matrix, vector, n, ref c);
            string str = "Приблеженная функция имеет вид: ";
            for (int i = 0; i < n; i++)
            {
                if (i.Equals(0))
                {
                    str += (c[i].ToString("0.000000"));
                }
                else
                {
                    str += (Abs(c[i]).ToString("0.000000") + (i > 1 ? ("x^" + i) : "x"));
                }
                if (!i.Equals(n - 1))
                {
                    str += (c[i + 1] > 0 ? " + " : " - ");
                }
            }
            PrintStr(str);
            sum = 0.0;
            double sum1 = 0.0;
            func sumF = delegate (double x) { return 0; };
            for (int i = 0; i < n; i++)
            {
                sumF = SumFunc(CoefFunc(G(i), c[i]), sumF);
            }
            for (int k = 0; k <= Program.n; k++)
            {
                sum += table[1, k] * table[1, k];
                sum1 += sumF(table[0, k]) * sumF(table[0, k]);
            }
            var f = sum;
            var g = sum1;
            PrintStr("Погрешность приближения: " + Abs(f - g).ToString());
            PrintStr(); 
        }

        static void TableRevers()
        {
            for (int i = 0; i <= n; i++)
            {
                double t = table[0, i];
                table[0, i] = table[1, i];
                table[1, i] = t;
            }
        }

        private static void BackInterpolation(double y)
        {
            TableRevers();
            double x = InterPolynomial(y);
            PrintStr("Результат обратного интерполирования для y = " + y.ToString() + ": x =" + x.ToString());
            PrintStr("Значение функции в полученой точке x: y = " + Func(x));
            PrintStr("Оценка погрешности для этой точки: " + CalcError(Func(x), y).ToString("0.000000"));
            PrintStr();
            TableRevers();
        }

        static double CalcError(double expected, double recived) => Abs(expected - recived);
        
        static void PrintStr(string str)
        {
            outFile.WriteLine(str);
            Console.WriteLine(str);
        }
        static void PrintStr()
        {
            outFile.WriteLine();
            Console.WriteLine();
        }

        static void Main(string[] args)
        {
            PrintFunc();
            FillTable();
            PrintTable();
            PrintStr("Интерполяция методом Ньютона");
            NyutonInterpolation();
            PrintStr("Cреднеквадратичное приближение табличным методом");
            TableApprox();
            PrintStr("Cреднеквадратичное приближение непрерывным методом");
            СontinuousАpprox();
            PrintStr("Обратное интерполирование по формуле Ньютона");
            BackInterpolation(1.5);
            outFile.Close();
        }
    }
}
