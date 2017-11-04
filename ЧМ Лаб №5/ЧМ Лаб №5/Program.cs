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
            switch (x.Length)
            {
                case 1:
                    return Func(x[0]);
                case 2:
                    return (Func(x[1]) - Func(x[0])) / (x[1] - x[0]);
                default:
                    double[] x1 = new double[x.Length - 1], x2 = new double[x.Length - 1];
                    for (var i = 0; i < x.Length - 1; i ++)
                    {
                        x1[i] = x[i];
                    }
                    for (var i = 1; i < x.Length; i++)
                    {
                        x2[i - 1] = x[i];
                    }
                    return (CalcNeigbours(x2) - CalcNeigbours(x1)) / (x[x.Length - 1] - x[0]);
            }
        }

        private static double InterPolynomial(double x)
        {
            double sum = 0;
            var xi = new List<double>();
            for (var i = 0; i < N; i++)
            {
                xi.Add(Table[i].X);
                sum += CalcNeigbours(xi.ToArray()) * OmegaFunc(x, i);
            }
            return sum;
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

        private static func ProdFunc(func f, func g) => x => f(x) * g(x);
        private static void СontinuousАpprox()
        {
            const int n = 3;
            var matrix = new double[n, n];
            double[] vector = new double[n], c = new double[n];
            for(var i = 0; i < n; i++)
            {
                for (var j = 0; j < n; j++)
                {
                    matrix[i, j] = CalcIntegral(ProdFunc(G(i), G(j)), A, B);
                }
                vector[i] = CalcIntegral(ProdFunc(Func, G(i)), A, B);
            }
            ConjurateGradientsMethod(matrix, vector, n, ref c);
            double sum = 0;
            for(var i = 0; i < n; i++)
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
            PrintStr(Sqrt(CalcIntegral(ProdFunc(Func, Func), A, B) - sum).ToString());
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
            PrintStr("x     Ожидаемое  Полученное   Погрешность");
            for (var i = 1.1; i < 2; i += 0.2)
            {
                var recived = InterPolynomial(i);
                var expected = Func(i);
                PrintStr(i + "   " + expected.ToString("0.000000") + "   " + recived.ToString("0.000000") + "     " + CalcError(expected, recived));
            }
            PrintStr();
            PrintStr("1. Кубический сплайн дефекта 1:");
            // 1.2 кубические сплайны дефекта 1
            BuildSpline(Table, N + 1);
            PrintStr("x     Ожидаемое  Полученное   Погрешность");
            for (var i = 1.1; i < 2; i += 0.2)
            {
                var recived = InterSplines(i);
                var expected = Func(i);
                PrintStr(i + "   " + expected.ToString("0.000000") + "   " + recived.ToString("0.000000") + "     " + CalcError(expected, recived));
            }
            PrintStr();
            СontinuousАpprox();
            PrintStr("Обратное интерполирование по формуле Ньютона");
            BackInterpolation(1.5);
            outFile.Close();
        }
    }
}
