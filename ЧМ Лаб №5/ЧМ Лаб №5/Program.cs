using System;
using System.IO;
using static System.Math;
using static ЧМ_Лаб__5.CubicSplines;
using static ЧМ_Лаб__5.LASSolver;

namespace ЧМ_Лаб__5
{
    internal static class Program
    {
        private delegate double func(double x);

        public const int N = 5;
        private const int Var = 16;
        public const double A = 1, B = 2;
        public const double Delta = (B - A) / N;

        public static readonly PointD[] Table = new PointD[N + 1];
        public static readonly StreamWriter OutFile = new StreamWriter("out.txt");

        public static double Func(double x)
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

        // Печатает вид четвертой производной f
        public static void PrintFourthDerivative()
        {
            switch (Var)
            {
                case 16: PrintStr("Вид 4ой производной f: 24 * x * (5 * x^-7 + (1 - x * x) * (x * x + 1)^-4))"); break;
                case 17: PrintStr("Вид 4ой производной f: e^x"); break;
                default: PrintStr();
            }
        }

        // Пятая производная f
        public static double FifthDerivative(double x)
        {
            switch (Var)
            {
                case 16: return 24 * (- 30 / Pow(x, 7) - 12 * x * x / Pow(x * x + 1, 4) + 1 / Pow(x * x + 1, 3) + 16 * Pow(x, 4) / Pow(x * x + 1, 5));
                case 17: return Exp(x);
                default: return 0;
            }
        }

        // Печатает вид пятой производной f
        public static void PrintFifthDerivative()
        {
            switch (Var)
            {
                case 16: PrintStr("Вид 5ой производной f: 24 * (-30 / x^7 - 12 * x * x / (x * x + 1)^4 + 1 / (x * x + 1)^3 + 16 * x^4 / (x * x + 1)^5)"); break;
                case 17: PrintStr("Вид 5ой производной f: e^x"); break;
                default: PrintStr();
            }
        }

        // Шестая производная f
        public static double SixthDerivative(double x)
        {
            switch (Var)
            {
                case 16: return x * (5040 * Pow(x, -9) + Pow(x * x + 1, -6) * (2400 * x * x - 720 * Pow(x, 4) - 720));
                case 17: return Exp(x);
                default: return 0;
            }
        }

        //печать функции
        private static void PrintFunc()
        {
            string str;
            switch (Var)
            {
                case 16: str = "Atan(x) + 1 / x^2"; break;
                case 17: str = "Exp(x) + x + 1"; break;
                default: str = "0"; break;
            }
            PrintStr("Рассматриваемая функция: " + str);
            PrintStr();
        }

        //заполнение таблицы x/y
        private static void FillTable()
        {
            var x = A;
            
            for(var i = 0; i <= N; i++)
            {
                Table[i].X = x;
                Table[i].Y = Func(x);
                x = Round(x + Delta, 1);
            }
        }

        //печать таблицы x/y
        private static void PrintTable()
        {
            PrintStr("Таблица значений на отрезке [" + A + ";" + B + "] с шагом delta = " + ((B - A) / N));
            string xStr = "x:", yStr = "y:";
            for (var i = 0; i <= N; i++)
            {
                xStr += "\t" + Table[i].X.ToString("0.0");
                yStr += "\t" + Table[i].Y.ToString("0.00000");
            }
            PrintStr(xStr);
            PrintStr(yStr);
            PrintStr();
        }
        
        //составление таблицы соседних разностей
        static void NeigboursTable()
        {
            double[,] table = new double[2 * N, N];
            for (int i = 0; i < N; i++)
            {
                for (int j = i; j < 2 * N - i; j += 2)
                {
                    table[j, i] = CalcNeigbours((j - i) / 2, (j + i) / 2);
                }
            }
            PrintNeigboursTable(table);
        }

        //печать таблицы соседних разностей
        static void PrintNeigboursTable(double[,] table)
        {
            for(int i = 0; i < 2 * N; i++)
            {
                string outStr = "";
                for(int j = 0; j < N; j++)
                {
                    outStr += table[i, j] == 0 ? "\t" : table[i, j].ToString("0.000000") + "\t";
                }
                PrintStr(outStr);
            }
        }

        //Составление Омега функции
        static double OmegaFunc(double x, int num)
        {
            double res = 1;
            for(var i = 0; i < num; i++)
            {
                res *= x - Table[i].X;
            }
            return res;
        }

        //Вычисление соседних разностей
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

        //Полином интерполирования для метода Ньютона
        static double InterPolynomial(double x)
        {
            double sum = 0;
            for (int i = 0; i <= N; i++)
            {
                sum += CalcNeigbours(0, i) * OmegaFunc(x, i);
            }
            return sum;
        }

        //Вычисление факториала
        static double Factorial(double n)
        {
            if (n == 0)
                return 1;
            return n * Factorial(n - 1);
        }
      
        //Интерполяция методом Ньютона
        static void NyutonInterpolation()
        {
            var M6 = double.MinValue;
            for (var i = 1.0; i <= 2.0; i += Delta)
            {
                if (SixthDerivative(i) > M6)
                    M6 = SixthDerivative(i);
            }
            PrintStr("Заначение M6 =" + M6.ToString());
            double recived = 0, expected = 0;
            PrintStr("x     Ожидаемое  Полученное   Реальная погрешность    Оценка погрешности");
            for (double i = 1.1; i < 2; i += 0.2)
            {
                recived = InterPolynomial(i);
                expected = Func(i);
                PrintStr(i + "   " + expected.ToString("0.000000") + "   " + recived.ToString("0.000000") + "     " 
                    + CalcError(expected, recived) + "    " + Abs(M6 * OmegaFunc(i, N + 1) / Factorial(N + 1)));
            }
            PrintStr();
        }

        //Интерполяция методом кубического сплайна
        static void CubSpline()
        {
            CubicSplines.CalcSplineParams();
            var teoreticError = CalcTeoreticError();
            PrintM4M5();
            PrintStr("x     Ожидаемое  Полученное   Оценка              Погрешность");
            for (var i = 1.1; i < 2; i += Delta)
            {
                var recived = CubicSplines.CalcSplineValue(i);
                var expected = Func(i);
                PrintStr(i + "   " + expected.ToString("0.000000") + "   " + recived.ToString("0.000000") + "     " + teoreticError + "     " + CalcError(expected, recived));
            }
            PrintStr();
            teoreticError = 1.0 / 60 * Pow(Delta, 4) * M5;
            PrintStr("Погрешности для первой производной:");
            PrintStr("x     Ожидаемое   Полученное    Оценка            Погрешность");
            var j = 0;
            for (var i = 1.0; i <= 2.0; i += Delta)
            {
                var recived = _splineParams[j];
                var expected = FirstDerivative(i);
                j++;
                PrintStr(i.ToString("F1") + "   " + expected.ToString("0.000000") + "   " + recived.ToString("0.000000") + "     " + teoreticError + "     " + CalcError(expected, recived));
            }
            PrintStr();
        }

        //Вычисление интеграла
        private static double CalcIntegral(func f, int aim = 10000)
        {
            double sum = 0, a, b, delta = (B - A) / aim;
            for (int i = 0; i < aim; i++)
            {
                a = A + delta * i;
                b = A + delta * (i + 1);
                sum += (b - a) * (f(a) + f(b) + 4 * f((a + b) / 2)) / 6;
            }
            return sum;
        }

        //Функции для работы с делегатами
        private static func G(int num) => x => Pow(x, num);
        static func CoefFunc(func f, double k) => x => k * f(x);
        private static func SumFunc(func f, func g) => x => f(x) + g(x);
        private static func ProdFunc(func f, func g) => x => f(x) * g(x);

        //Метод непрерывной апроксимации
        static void СontinuousАpprox()
        {
            int n = 3;
            double[,] matrix = new double[n, n];
            double[] vector = new double[n], c = new double[n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    matrix[i, j] = CalcIntegral(ProdFunc(G(i), G(j)));
                }
                vector[i] = CalcIntegral(ProdFunc(Func, G(i)));
            }
            PrintStr("Матрица C");
            PrintMatrix(matrix, n);
            PrintStr("Вектор B");
            PrintVector(vector, n);
            ConjurateGradientsMethod(matrix, vector, n, ref c);
            string str = "Приближенная функция имеет вид: ";
            for (int i = 0; i < n; i++)
            {
                if (i.Equals(0))
                {
                    str += c[i].ToString("0.000000");
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
            var f = CalcIntegral(ProdFunc(Func, Func));
            double g = 0;
            for (int i = 0; i < n; i++)
                g += CalcIntegral(ProdFunc(Func, CoefFunc(G(i),c[i])));
            g = Abs(g);
            PrintStr("Погрешность приближения: " + Sqrt(Abs(f - g)).ToString());
            PrintStr();
        }

        //Метод табличной апроксимации
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
            PrintStr("Матрица C");
            PrintMatrix(matrix, n);
            PrintStr("Вектор B");
            PrintVector(vector, n);
            ConjurateGradientsMethod(matrix, vector, n, ref c);
            string str = "Приближенная функция имеет вид: ";
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
            double f = 0.0, g = 0.0;
            for (int i = 0; i <= N; i++)
            {
                f += Table[i].Y * Table[i].Y;
            }   
            for(int i = 0; i < n; i++)
            {
                for(int j = 0; j <= N; j++)
                {
                    g += Table[j].Y * c[i] * G(i)(Table[j].X);
                }
            }
            PrintStr("Погрешность приближения: " + Sqrt(Abs(f - g)).ToString());
            PrintStr();
        }

        //Разворот таблицы для обратного интерполирования Ньютона
        static void TableRevers()
        {
            for (int i = 0; i <= N; i++)
            {
                double t = Table[i].X;
                Table[i].X = Table[i].Y;
                Table[i].Y = t;
            }
        }

        //Обратное интерполирование Ньютона
        private static void BackInterpolation(double y)
        {
            TableRevers();
            PrintStr("Таблица обратных разностей для обратного тетода Ньютона");
            NeigboursTable();
            double x = InterPolynomial(y);
            PrintStr("С = " + y);
            PrintStr("Результат обратного интерполирования для y = " + y.ToString() + ": x = " + x.ToString());
            PrintStr("Значение функции в полученой точке x: y = " + Func(x));
            PrintStr("Погрешности для этой точки: " + CalcError(Func(x), y).ToString("0.000000"));
            PrintStr();
            TableRevers();
        }
      
        //Вычисление ошибки
        private static double CalcError(double expected, double recived) => Abs(expected - recived);

        //Печать строки в файл и консоль
        public static void PrintStr(string str)
        {
            OutFile.WriteLine(str);
            Console.WriteLine(str);
        }
        public static void PrintStr()
        {
            OutFile.WriteLine();
            Console.WriteLine();
        }

        private static void Main(string[] args)
        {
            PrintFunc();
            FillTable();
            PrintTable();
            PrintStr("Таблица соседних разностей");
            NeigboursTable();
            PrintStr("1. Формула Ньютона:");
            NyutonInterpolation();
            PrintStr("1. Кубический сплайн дефекта 1:");
            CubSpline();
            PrintStr("2. Cреднеквадратичное приближение табличным методом");
            TableApprox();
            PrintStr("2. Cреднеквадратичное приближение непрерывным методом");
            СontinuousАpprox();
            PrintStr("3. Обратное интерполирование по формуле Ньютона");
            BackInterpolation(7);
            OutFile.Close();
        }
    }
}