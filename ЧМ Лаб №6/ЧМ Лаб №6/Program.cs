using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using static System.Math;

namespace ЧМ_Лаб__6
{
    class Program
    {
        private const int N = 10000, Var = 17;
        private const double A = 1, B = 2;
        private const double Delta = (B - A) / N;
        private static readonly StreamWriter OutFile = new StreamWriter("out.txt");

        private static double Func(double x)
        {
            switch (Var)
            {
                case 16: return Atan(x) + 1 / x / x;
                case 17: return Exp(x) + x + 1;
                default: return 0;
            }
        }

        // Первая производная f
        public static double FirstDerivative(double x)
        {
            switch (Var)
            {
                case 16: return 1 / (x * x + 1) - 2 / (x * x * x);
                case 17: return Exp(x) + 1;
                default: return 0;
            }
        }

        private static double RealValue
        {
            get
            {
                switch (Var)
                {
                    case 16: return 1.470753906253655163826706468653192823365752331004957604213;
                    case 17: return 7.170774270471604991870139989222345315423068476851887749120;
                    default: return 0;
                }
            }
        }

        private static void Main(string[] args)
        {
            var t = CalcTrapezeIntegral();
            PrintStr("Метод трапеции");
            PrintStr("Приближенное значение: " + DoubleConverter.ToExactString(t));
            PrintStr("Погрешность: " + Abs(t - RealValue));

            t = CalcModifiedTrapezeIntegral();
            PrintStr("Модифицированный сплайном метод трапеции");
            PrintStr("Приближенное значение: " + DoubleConverter.ToExactString(t));
            PrintStr("Погрешность: " + Abs(t - RealValue));

            t = CalcSimpsonIntegral();
            PrintStr("Метод Симпсона");
            PrintStr("Приближенное значение: " + DoubleConverter.ToExactString(t));
            PrintStr("Погрешность: " + Abs(t - RealValue));

            t = CalcGaussIntegral();
            PrintStr("Метод Гаусса(трехточечный)");
            PrintStr("Приближенное значение: " + DoubleConverter.ToExactString(t));
            PrintStr("Погрешность: " + Abs(t - RealValue));
        }

        private static double CalcTrapezeIntegral()
        {
            double sum = 0, a, b;
            for(int i = 0; i < N; i++)
            {
                a = A + Delta * i;
                b = A + Delta * (i + 1);
                sum += (b - a) * (Func(a) + Func(b)) / 2;
            }
            return sum;
        }

        private static double CalcModifiedTrapezeIntegral()
        {
            double sum = 0;
            for (var i = A + Delta; i < B; i += Delta)
            {
                sum += Func(i);
            }
            return Delta * (0.5 * (Func(A) + Func(B)) + sum) +
                   Delta * Delta / 12 * (FirstDerivative(A) - FirstDerivative(B));
        }

        private static double CalcSimpsonIntegral()
        {
            double sum = 0, a, b;
            for (int i = 0; i < N; i++)
            {
                a = A + Delta * i;
                b = A + Delta * (i + 1);
                sum += (b - a) * (Func(a) + Func(b) + 4 * Func((a + b) / 2)) / 6;
            }
            return sum;
        }

        private static double CalcGaussIntegral()
        {
            double sum = 0, a, b;
            for (int i = 0; i < N; i++)
            {
                a = A + Delta * i;
                b = A + Delta * (i + 1);
                sum += (b - a) * (5.0 / 9.0 * Func(0.5 *(a + b + (b - a) * -Sqrt(3.0 / 5.0))) + 8.0 / 9.0 * Func(0.5 * (a + b + (b - a))) + 5.0 / 9.0 * Func(0.5 * (a + b + (b - a) * Sqrt(3.0 / 5.0)))) / 2;
            }
            return sum;
        }

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
    }
}
