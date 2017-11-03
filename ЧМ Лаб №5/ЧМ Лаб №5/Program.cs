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
        const int n = 5, var = 17;
        const double a = 1, b = 2;
        static double[,] table = new double[2, n];
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

        static void FillTable()
        {
            double x = a;
            double delta = (b - a) / n;
            for(int i = 0; i < n; i++)
            {
                table[0, i] = x;
                table[1, i] = Func(x);
                x = Round(x + delta, 1);
            }
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

        static double CalcNeigbours(double[] x)
        {
            if(x.Length.Equals(1))
            {
                return Func(x[0]);
            }
            if(x.Length.Equals(2))
            {
                return (Func(x[1]) - Func(x[0])) / (x[1] - x[0]);
            }
            else
            {
                double[] x1 = new double[x.Length - 1], x2 = new double[x.Length - 1];
                for (int i = 0; i < x.Length - 1; i ++)
                {
                    x1[i] = x[i];
                }
                for (int i = 1; i < x.Length; i++)
                {
                    x2[i - 1] = x[i];
                }
                return (CalcNeigbours(x2) - CalcNeigbours(x1)) / (x[x.Length - 1] - x[0]);
            }
        }

        static double InterPolynomial(double x)
        {
            double sum = 0;
            List<double> xi = new List<double>();
            for(int i = 0; i < n; i++)
            {
                xi.Add(table[0, i]);
                sum += CalcNeigbours(xi.ToArray()) * OmegaFunc(x, i);
            }
            return sum;
        }

        static double CalcIntegral(func f, double a, double b, int aim = 1000)
        {
            double sum = 0, delta = (b - a) / aim;
            for (double i = a; i < b; i+= delta)
            {
                sum += f(i);
            }
            return Round(sum * delta, aim.ToString().Length - 1);
        }

        static func G(int num) => delegate (double x) { return Pow(x, num); };

        static func ProdFunc(func f, func g) => delegate (double x) { return f(x) * g(x); };

        static void СontinuousАpprox()
        {
            int n = 3;
            double[,] matrix = new double[n, n];
            double[] vector = new double[n], c = new double[n];
            for(int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    matrix[i, j] = CalcIntegral(ProdFunc(G(i), G(j)), a, b);
                }
                vector[i] = CalcIntegral(ProdFunc(Func, G(i)), a, b);
            }
            conjurateGradientsMethod(matrix, vector, n, ref c);
            double sum = 0;
            for(int i = 0; i < n; i++)
            {
                sum += Pow(vector[i], 2);
            }
            PrintStr(Sqrt(CalcIntegral(ProdFunc(Func, Func), a, b) - sum).ToString());
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
            FillTable();
            double recived = 0, expected = 0;
            PrintStr("x     Ожидаемое  Полученное   Погрешность");
            for (double i = 1.1; i < 2; i += 0.2)
            {
                recived = InterPolynomial(i);
                expected = Func(i);
                PrintStr(i.ToString() + "   " + expected.ToString("0.000000") + "   " + recived.ToString("0.000000") + "     " + CalcError(expected,recived));
            }
            СontinuousАpprox();
        }
    }
}
