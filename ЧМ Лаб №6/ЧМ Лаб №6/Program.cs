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
        private const int N = 1000, Var = 16;
        private const double A = 1, B = 2;
        private static double M2 = double.MinValue, M4 = double.MinValue;
        private const double Delta = (B - A) / N;
        private static readonly StreamWriter OutFile = new StreamWriter("out.txt");

        // Функция f
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

        // Вторая производная f
        public static double SecondDerivative(double x)
        {
            switch (Var)
            {
                case 16: return 2 * (-x / Pow(x * x + 1, 2) + 3 / (x * x * x * x));
                case 17: return Exp(x);
                default: return 0;
            }
        }
        // Четвертая производная f
        public static double FourthDerivative(double x)
        {
            switch (Var)
            {
                case 16: return 24 * x * (5 * Pow(x, -7) + (1 - x * x) * Pow(x * x + 1, -4));
                case 17: return Exp(x);
                default: return 0;
            }
        }

        // Значение интеграла, полученное через Вольфрам
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

        private static void CalcM2()
        {
            for (int i = 0; i < N; i++)
            {
                double t = Abs(SecondDerivative(A + i * Delta));
                if (t > M2)
                {
                    M2 = t;
                }
            }
        }

        private static void CalcM4()
        {
            for (int i = 0; i < N; i++)
            {
                double t = Abs(FourthDerivative(A + i * Delta));
                if (t > M4)
                {
                    M4 = t;
                }
            }
        }

        // Метод трапеции
        private static void CalcTrapezeIntegral()
        {
            double sum = 0, a, b;
            for(int i = 0; i < N; i++)
            {
                a = A + Delta * i;
                b = A + Delta * (i + 1);
                sum += (b - a) * (Func(a) + Func(b)) / 2;
            }
            PrintStr("Метод трапеции");
            PrintStr("Приближенное значение: " + DoubleConverter.ToExactString(sum));
            PrintStr("Ожидаемая погрешность: " + CalcErrorForTrapeze());
            PrintStr("Реальная погрешность: " + Abs(sum - RealValue));
            PrintStr();
        }

        // Ожидаемая погрешность для метода трапеции
        private static double CalcErrorForTrapeze()
        {
            return Abs(Delta * Delta / 12 * (B - A) * M2);
        }

        /*private static void AltCalcTrapezeIntegral()
        {
            double sum = 0;
            for (int i = 1; i < N; i++)
            {
                sum += Func(A + i * Delta);
            }
            sum = Delta * ((Func(A) + Func(B)) / 2 + sum);
            PrintStr("Метод трапеции");
            PrintStr("Приближенное значение: " + DoubleConverter.ToExactString(sum));
            PrintStr("Реальная погрешность: " + Abs(sum - RealValue));
        }*/

        // Метод трапеции(модифицированный)
        private static void CalcModifiedTrapezeIntegral()
        {
            double sum = 0;
            for (int i = 1; i < N; i++)
            {
                sum += Func(A + i * Delta);
            }
            sum = Delta * (0.5 * (Func(A) + Func(B)) + sum) +
                   Delta * Delta / 12 * (FirstDerivative(A) - FirstDerivative(B));
            PrintStr("Модифицированный сплайном метод трапеции");
            PrintStr("Приближенное значение: " + DoubleConverter.ToExactString(sum));
            PrintStr("Ожидаемая погрешность: " + CalcErrorForTrapeze());
            PrintStr("Реальная погрешность: " + Abs(sum - RealValue));
            PrintStr();
        }

        private static void CalcSimpsonIntegral()
        {
            double sum = 0, a, b;
            for (int i = 0; i < N; i++)
            {
                a = A + Delta * i;
                b = A + Delta * (i + 1);
                sum += (b - a) * (Func(a) + Func(b) + 4 * Func((a + b) / 2)) / 6;
            }
            PrintStr("Метод Симпсона");
            PrintStr("Приближенное значение: " + DoubleConverter.ToExactString(sum));
            PrintStr("Ожидаемая погрешность: " + CalcErrorForSimpson());
            PrintStr("Реальная погрешность: " + Abs(sum - RealValue));
            PrintStr();
        }

        /*private static void AltCalcSimpsonIntegral()
        {
            double sum1 = 0, sum2 = 0;
            for (int i = 1; i <= N / 2; i++)
            {
                sum1 += Func(A + (2 * i - 1) * Delta);
                if(i != N / 2)
                    sum2 += Func(A + 2 * i * Delta);
            }
            double sum = Delta / 3 * ((Func(A) + Func(B)) + 4 * sum1 + 2 * sum2);
            PrintStr("Метод Симпсона");
            PrintStr("Приближенное значение: " + DoubleConverter.ToExactString(sum));
            PrintStr("Ожидаемая погрешность: " + CalcErrorForSimpson());
            PrintStr("Реальная погрешность: " + Abs(sum - RealValue));
            PrintStr();
        }*/

        private static double CalcErrorForSimpson()
        {
            return Abs(Pow(Delta, 4) * (B - A) / 2880 * M4);
        }

        private static void CalcGaussIntegral()
        {
            double sum = 0, a, b;
            for (int i = 0; i < N; i++)
            {
                a = A + Delta * i;
                b = A + Delta * (i + 1);
                sum += (b - a) * (5.0 / 9.0 * Func(0.5 *(a + b + (b - a) * -Sqrt(3.0 / 5.0))) +
                    8.0 / 9.0 * Func(0.5 * (a + b)) + 5.0 / 9.0 * Func(0.5 * (a + b + (b - a) * Sqrt(3.0 / 5.0)))) / 2;
            }
            PrintStr("Метод Гаусса(трехточечный)");
            PrintStr("Приближенное значение: " + DoubleConverter.ToExactString(sum));
            PrintStr("Реальная погрешность: " + Abs(sum - RealValue));
            PrintStr();
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

        private static void Main(string[] args)
        {

            CalcM2();
            CalcM4();

            CalcTrapezeIntegral();

            //AltCalcTrapezeIntegral();

            CalcModifiedTrapezeIntegral();

            CalcSimpsonIntegral();

            //AltCalcSimpsonIntegral();

            CalcGaussIntegral();
        }
    }
}
