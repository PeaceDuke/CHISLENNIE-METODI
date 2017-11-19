using System;
using System.IO;
using static System.Math;

namespace ЧМ_Лаб__6
{
    class Program
    {
        private const int Var = 16;
        private const double A = 1, B = 2;
        private static readonly double Eps = Pow(10, -8);
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

        // Первая производная f
        private static double FirstDerivative(double x)
        {
            switch (Var)
            {
                case 16: return 1 / (x * x + 1) - 2 / (x * x * x);
                case 17: return Exp(x) + 1;
                default: return 0;
            }
        }

        // Вторая производная f
        private static double SecondDerivative(double x)
        {
            switch (Var)
            {
                case 16: return 2 * (-x / Pow(x * x + 1, 2) + 3 / (x * x * x * x));
                case 17: return Exp(x);
                default: return 0;
            }
        }
        // Четвертая производная f
        private static double FourthDerivative(double x)
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

        private static double CalcM2(int n)
        {
            var M2 = double.MinValue;
            var delta = (B - A) / n;
            for (var i = 0; i < n; i++)
            {
                var t = Abs(SecondDerivative(A + i * delta));
                if (t > M2)
                {
                    M2 = t;
                }
            }
            return M2;
        }

        private static double CalcM4(int n)
        {
            var M4 = double.MinValue;
            var delta = (B - A) / n;
            for (var i = 0; i < n; i++)
            {
                var t = Abs(FourthDerivative(A + i * delta));
                if (t > M4)
                {
                    M4 = t;
                }
            }
            return M4;
        }

        private static double CalcApproxOrder(int n, int methodNumb)
        {
            switch (methodNumb)
            {
                case 1:
                    return 1 / Log(0.5) * Log((CalcTrapezeIntegral(4 * n, true) - CalcTrapezeIntegral(n, true)) 
                        / (CalcTrapezeIntegral(2 * n, true) - CalcTrapezeIntegral(n, true)));
                case 2:
                    return 1 / Log(0.5) * Log((CalcModifiedTrapezeIntegral(4 * n, true) - CalcModifiedTrapezeIntegral(n, true)) 
                        / (CalcModifiedTrapezeIntegral(2 * n, true) - CalcModifiedTrapezeIntegral(n, true)));
                case 3:
                    return 1 / Log(0.5) * Log((CalcSimpsonIntegral(4 * n, true) - CalcSimpsonIntegral(n, true)) 
                        / (CalcSimpsonIntegral(2 * n, true) - CalcSimpsonIntegral(n, true)));
                case 4:
                    return 1 / Log(0.5) * Log((CalcGaussIntegral(4 * n, true) - CalcGaussIntegral(n, true)) 
                        / (CalcGaussIntegral(2 * n, true) - CalcGaussIntegral(n, true)));
                default:
                    return 0;
            }
        }

        // Метод трапеции
        private static double CalcTrapezeIntegral(int n, bool nAlreadyFound = false)
        {
            double sum = 0;
            if (nAlreadyFound)
            {
                var delta = (B - A) / n;
                for (var i = 0; i < n; i++)
                {
                    var a = A + delta * i;
                    var b = A + delta * (i + 1);
                    sum += (b - a) * (Func(a) + Func(b)) / 2;
                }
                return sum;
            }
            PrintStr("=== Метод трапеции ===");
            PrintStr("N    Значение интеграла                        Ожидаемая погрешность                Реальная погрешность");
            double realError;
            do
            {
                sum = 0;
                n *= 2;
                var delta = (B - A) / n;
                for (var i = 0; i < n; i++)
                {
                    var a = A + delta * i;
                    var b = A + delta * (i + 1);
                    sum += (b - a) * (Func(a) + Func(b)) / 2;
                }
                var teoreticError = CalcErrorForTrapeze(n);
                realError = Abs(sum - RealValue);
                PrintStr(n + "    " + DoubleConverter.ToExactString(sum) + "   " + teoreticError + "   " + realError);
            } while (realError > Eps);
            PrintStr("Порядок аппроксимации: " + CalcApproxOrder(n, 1));
            PrintStr();
            return sum;
        }

        // Ожидаемая погрешность для метода трапеции
        private static double CalcErrorForTrapeze(int n)
        {
            var delta = (B - A) / n;
            return Abs(delta * delta / 12 * (B - A) * CalcM2(n));
        }

        // Метод трапеции(модифицированный)
        private static double CalcModifiedTrapezeIntegral(int n, bool nAlreadyFound = false)
        {
            double sum = 0;
            if (nAlreadyFound)
            {
                var delta = (B - A) / n;
                for (var i = 0; i < n; i++)
                {
                    sum += Func(A + i * delta);
                }
                return delta * (0.5 * (Func(A) + Func(B)) + sum) +
                      delta * delta / 12 * (FirstDerivative(A) - FirstDerivative(B));
            }
            PrintStr("=== Модифицированный метод трапеции ===");
            PrintStr("N    Значение интеграла                        Ожидаемая погрешность                Реальная погрешность");
            double realError;
            do
            {
                sum = 0;
                n *= 2;
                var delta = (B - A) / n;
                for (var i = 0; i < n; i++)
                {
                    sum += Func(A + i * delta);
                }
                sum = delta * (0.5 * (Func(A) + Func(B)) + sum) +
                      delta * delta / 12 * (FirstDerivative(A) - FirstDerivative(B));
                var teoreticError = CalcErrorForTrapeze(n);
                realError = Abs(sum - RealValue);
                PrintStr(n + "    " + DoubleConverter.ToExactString(sum) + "   " + teoreticError + "   " + realError);
            } while (realError > Eps);
            PrintStr("Порядок аппроксимации: " + CalcApproxOrder(n, 1));
            PrintStr();
            return sum;
        }

        // Метод Симпсона
        private static double CalcSimpsonIntegral(int n, bool nAlreadyFound = false)
        {
            double sum = 0;
            if (nAlreadyFound)
            {
                var delta = (B - A) / n;
                for (var i = 0; i < n; i++)
                {
                    var a = A + delta * i;
                    var b = A + delta * (i + 1);
                    sum += (b - a) * (Func(a) + Func(b) + 4 * Func((a + b) / 2)) / 6;
                }
                return sum;
            }
            PrintStr("=== Метод Симпсона ===");
            PrintStr("N    Значение интеграла                        Ожидаемая погрешность                Реальная погрешность");
            double realError;
            do
            {
                sum = 0;
                n *= 2;
                var delta = (B - A) / n;
                for (var i = 0; i < n; i++)
                {
                    var a = A + delta * i;
                    var b = A + delta * (i + 1);
                    sum += (b - a) * (Func(a) + Func(b) + 4 * Func((a + b) / 2)) / 6;
                }
                var teoreticError = CalcErrorForSimpson(n);
                realError = Abs(sum - RealValue);
                PrintStr(n + "    " + DoubleConverter.ToExactString(sum) + "   " + teoreticError + "   " + realError);
            } while (realError > Eps);
            PrintStr("Порядок аппроксимации: " + CalcApproxOrder(n, 1));
            PrintStr();
            return sum;
        }

        // Ожидаемая погрешность для метода Симпсона
        private static double CalcErrorForSimpson(int n)
        {
            var delta = (B - A) / n;
            return Abs(Pow(delta, 4) * (B - A) / 2880 * CalcM4(n));
        }

        // Метод Гаусса
        private static double CalcGaussIntegral(int n, bool nAlreadyFound = false)
        {
            double sum = 0;
            if (nAlreadyFound)
            {
                var delta = (B - A) / n;
                for (var i = 0; i < n; i++)
                {
                    var a = A + delta * i;
                    var b = A + delta * (i + 1);
                    sum += (b - a) * (5.0 / 9.0 * Func(0.5 * (a + b + (b - a) * -Sqrt(3.0 / 5.0))) +
                                      8.0 / 9.0 * Func(0.5 * (a + b)) + 5.0 / 9.0 * Func(0.5 * (a + b + (b - a) * Sqrt(3.0 / 5.0)))) / 2;
                }
                return sum;
            }
            PrintStr("=== Метод Гаусса ===");
            PrintStr("N    Значение интеграла                        Ожидаемая погрешность                Реальная погрешность");
            double realError;
            do
            {
                sum = 0;
                n *= 2;
                var delta = (B - A) / n;
                for (var i = 0; i < n; i++)
                {
                    var a = A + delta * i;
                    var b = A + delta * (i + 1);
                    sum += (b - a) * (5.0 / 9.0 * Func(0.5 *(a + b + (b - a) * -Sqrt(3.0 / 5.0))) +
                                      8.0 / 9.0 * Func(0.5 * (a + b)) + 5.0 / 9.0 * Func(0.5 * (a + b + (b - a) * Sqrt(3.0 / 5.0)))) / 2;
                }
                //todo var teoreticError = ;
                realError = Abs(sum - RealValue);
                PrintStr(n + "    " + DoubleConverter.ToExactString(sum) + "   " + realError);
            } while (realError > Eps);
            PrintStr("Порядок аппроксимации: " + CalcApproxOrder(n, 1));
            PrintStr();
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

        private static void Main(string[] args)
        {
            CalcTrapezeIntegral(1);

            CalcModifiedTrapezeIntegral(1);

            CalcSimpsonIntegral(1);

            CalcGaussIntegral(1);
        }
    }
}
