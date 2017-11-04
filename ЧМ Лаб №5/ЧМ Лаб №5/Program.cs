using System;
using System.Collections.Generic;
using System.Globalization;
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
        private const double Delta = (B - A) / N;

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

        // Пятая производная f
        private static double FifthDerivative(double x)
        {
            switch (Var)
            {
                //case 16: return 24 * Pow(x * x + 1, -5) * (5 * Pow(x, 4) - 10 * x * x + 1) - 720 * Pow(x, -7);
                case 16: return 24 * (- 30 / Pow(x, 7) - 12 * x * x / Pow(x * x + 1, 4) + 1 / Pow(x * x + 1, 3) + 16 * Pow(x, 4) / Pow(x * x + 1, 5));
                case 17: return Exp(x);
                default: return 0;
            }
        }

        private static void PrintFunc()
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
            
            for(var i = 0; i <= N; i++)
            {
                Table[i].X = x;
                Table[i].Y = Func(x);
                x = Round(x + Delta, 1);
            }
        }

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
            PrintTable();
            PrintStr("1. Кубический сплайн дефекта 1:");
            // 1.2 кубические сплайны дефекта 1
            BuildSpline(Table, N + 1);
            PrintStr("x     Ожидаемое  Полученное   Погрешность");
            for (var i = 1.1; i < 2; i += 0.2)
            {
                var recived = InterSplines(i, false);
                var expected = Func(i);
                PrintStr(i + "   " + expected.ToString("0.000000") + "   " + recived.ToString("0.000000") + "     " + CalcError(expected, recived));
            }
            PrintStr("Оценка погрешности аппроксимации значений первой производной в узлах интерполяции:");
            var m5 = double.MinValue;
            for (var i = 1.0; i <= 2; i += Delta)
            {
                var currentValue = Abs(FifthDerivative(i));
                if (currentValue > m5)
                    m5 = currentValue;
            }
            var max = double.MinValue;
            for (var i = 1.0; i <= 2; i += Delta)
            {
                var first = FirstDerivative(i);
                var second = InterSplines(i, true);
                var currentValue = Abs(first - second);
                if (currentValue > max)
                    max = currentValue;
            }
            PrintStr("max |f'i - mi| <= M5 / 60 * h ^ 4; i = 0,N");
            PrintStr(max + " <= " +  m5 * Pow(Delta, 4) / 60 );
        }
    }
}
