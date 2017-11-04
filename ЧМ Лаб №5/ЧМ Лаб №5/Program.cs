using System;
using System.Collections.Generic;
using System.IO;
using static System.Math;
using static ЧМ_Лаб__5.CubicSpline;
using static ЧМ_Лаб__5.LASSolver;

namespace ЧМ_Лаб__5
{
    internal static class Program
    {
        private delegate double func(double x);

        private const int N = 5, Var = 16;
        private const double A = 1, B = 2;

        private static readonly PointD[] Table = new PointD[N + 1];
        private static readonly StreamWriter OutFile = new StreamWriter("out.txt");

        private static double Func(double x)
        {
            switch (Var)
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
            switch (Var)
            {
                case 15: str = "x^2"; break;
                case 16: str = "Atan(x) + 1 / x^2"; break;
                case 17: str = "Exp(x) + x + 1"; break;
                default: str = "0"; break;
            }
            PrintStr("Рассматриваемая функция: " + str);
            PrintStr();
        }

        private static void FillTable()
        {
            var x = A;
            const double delta = (B - A) / N;
            for(var i = 0; i <= N; i++)
            {
                Table[i].X = x;
                Table[i].Y = Func(x);
                x = Round(x + delta, 1);
            }
        }
        static void PrintTable()
        {
            PrintStr("Таблица значений на отрезке [" + A + ";" + B + "] с шагом delta = " + ((B - A) / N));
            string xStr = "x:", yStr = "y:";
            for (int i = 0; i <= N; i++)
            {
                xStr += "\t" + Table[i].X.ToString("0.0");
                yStr += "\t" + Table[i].Y.ToString("0.00000");
            }
            PrintStr(xStr);
            PrintStr(yStr);
            PrintStr();
        }

        private static double OmegaFunc(double x, int num)
        {
            double res = 1;
            for(var i = 0; i < num; i++)
            {
                res *= x - Table[i].X;
            }
            return res;
        }

        static double CalcNeigbours(int a, int b)
        {
            if (a == b)
            {
                return Table[a].Y;
            }
            else
            {
                return (CalcNeigbours(a + 1, b) - CalcNeigbours(a, b - 1)) / (Table[b].X - Table[a].X);
            }
        }

        static double InterPolynomial(double x)
        {
            double sum = 0;
            for (int i = 0; i <= N; i++)
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

        static void QuadSpline()
        {
            PrintStr("x     Ожидаемое  Полученное   Погрешность");
            for (var i = 1.1; i < 2; i += 0.2)
            {
                var recived = InterSplines(i);
                var expected = Func(i);
                PrintStr(i + "   " + expected.ToString("0.000000") + "   " + recived.ToString("0.000000") + "     " + CalcError(expected, recived));
            }
            PrintStr();
        }

        private static double CalcIntegral(func f, double a, double b, int aim = 1000)
        {
            double sum = 0, delta = (b - a) / aim;
            for (var i = a; i < b; i+= delta)
            {
                sum += f(i);
            }
            return Round(sum * delta, aim.ToString().Length - 1);
        }

        private static func G(int num) => x => Pow(x, num);
        static func CoefFunc(func f, double k) => x => k * f(x);
        private static func SumFunc(func f, func g) => x => f(x) + g(x);
        private static func ProdFunc(func f, func g) => x => f(x) * g(x);

        static void СontinuousАpprox()
        {
            int n = 3, aim = 100000;
            double[,] matrix = new double[n, n];
            double[] vector = new double[n], c = new double[n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    matrix[i, j] = CalcIntegral(ProdFunc(G(i), G(j)), A, B, aim);
                }
                vector[i] = CalcIntegral(ProdFunc(Func, G(i)), A, B, aim);
            }
            ConjurateGradientsMethod(matrix, vector, n, ref c);
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
            func sum = delegate (double x) { return 0; };
            for (int i = 0; i < n; i++)
            {
                sum = SumFunc(CoefFunc(G(i), c[i]), sum);
            }
            var f = CalcIntegral(ProdFunc(Func, Func), A, B, aim);
            var g = CalcIntegral(ProdFunc(sum, sum), A, B, aim);
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
                    for (int k = 0; k <= N; k++)
                        sum += G(i)(Table[k].X) * G(j)(Table[k].X);
                    matrix[i, j] = sum;
                }
                sum = 0.0;
                for (int k = 0; k <= N; k++)
                    sum += Table[k].Y * G(i)(Table[k].X);
                vector[i] = sum;
            }
            ConjurateGradientsMethod(matrix, vector, n, ref c);
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
            for (int k = 0; k <= N; k++)
            {
                sum += Table[k].Y * Table[k].Y;
                sum1 += sumF(Table[k].X) * sumF(Table[k].X);
            }
            var f = sum;
            var g = sum1;
            PrintStr("Погрешность приближения: " + Abs(f - g).ToString());
            PrintStr();
        }

        static void TableRevers()
        {
            for (int i = 0; i <= N; i++)
            {
                double t = Table[i].X;
                Table[i].X = Table[i].Y;
                Table[i].Y = t;
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

        private static double CalcError(double expected, double recived) => Abs(expected - recived);

        private static void PrintStr(string str)
        {
            OutFile.WriteLine(str);
            Console.WriteLine(str);
        }

        private static void PrintStr()
        {
            OutFile.WriteLine();
            Console.WriteLine();
        }

        private static void Main(string[] args)
        {
            PrintFunc();
            FillTable();
            // 1.1 формула Ньютона
            PrintStr("1. Формула Ньютона:");
            NyutonInterpolation();
            PrintStr("1. Кубический сплайн дефекта 1:");
            // 1.2 кубические сплайны дефекта 1
            BuildSpline(Table, N + 1);
            QuadSpline();
            PrintStr("Cреднеквадратичное приближение табличным методом");
            TableApprox();
            PrintStr("Cреднеквадратичное приближение непрерывным методом");
            СontinuousАpprox();
            PrintStr("Обратное интерполирование по формуле Ньютона");
            BackInterpolation(1.5);
            OutFile.Close();
        }
    }
}