using System;
using System.Collections.Generic;
using System.IO;
using static System.Math;

namespace ЧМ_Лаб_4
{
    public delegate double D(double x, double y); 

    internal static class Program
    {
        private const double Eps = 0.0001;
        private static StreamWriter w;

        //функция f1.
        private static double F1(double x, double y, int variant)
        {
            switch (variant)
            {
                case 16:
                    return Cos(y + 0.5) + x - 0.8;
                case 17:
                    return Sin(y - 1) + x - 1.3;
            }
            return 0;
        }

        //производная f1 по x.
        private static double Dx_F1(double x, double y, int variant)
        {
            return 1;
        }

        //производная f1 по y.
        private static double Dy_F1(double x, double y, int variant)
        {
            switch (variant)
            {
                case 16:
                    return -Sin(y + 0.5);
                case 17:
                    return Cos(y - 1);
            }
            return 0;
        }

        //функция f2.
        private static double F2(double x, double y, int variant)
        {
            switch (variant)
            {
                case 16:
                    return Sin(x) - 2 * y - 1.6;
                case 17:
                    return y - Sin(x + 1) - 0.8;
            }
            return 0;
        }

        //производная f2 по x.
        private static double Dx_F2(double x, double y, int variant)
        {
            switch (variant)
            {
                case 16:
                    return Cos(x);
                case 17:
                    return Cos(x + 1);
            }
            return 0;
        }

        //производная f2 по y.
        private static double Dy_F2(double x, double y, int variant)
        {
            switch (variant)
            {
                case 16:
                    return -2;
                case 17:
                    return 1;
            }
            return 0;
        }

        static void Main(string[] args)
        {
            const int variant = 17;
            w = new StreamWriter("out.txt");
            //SimpleIterationMethod.
            double x;
            double y;
            if (variant.Equals(16))
            {
                x = -0.25;
                y = -1;
            }
            else
            {
                x = 0.5;
                y = 2;
            }
            PrintStr("Начальное приближение, определенное геометрически: x = " + x + " y = " + y);
            PrintStr();
            PrintStr("Метод простой итерации:");
            SimpleIterMethod(x, y, variant);
            PrintStr();
            //NewtonMethod.
            PrintStr("Метод Ньютона:");
            NewtonMethod(x, y, variant);
            PrintStr();
            //GradientMethod.
            PrintStr("Метод градиентного спуска:");
            GradientMethod(x, y, variant);
            w.Close();
        }

        //метод простой итерации.
        private static void SimpleIterMethod(double x, double y, int variant)
        {
            var iterCount = 0;
            PrintStr("№    x           y           Невязка     Погрешность");
            do
            {
                iterCount++;
                var newX = SIMethodFiX(x, y, variant);
                var newY = SIMethodFiY(x, y, variant);
                var pogreshnosty = Abs(newX - x);
                x = newX;
                y = newY;
                var error = Max(Abs(F1(x, y, variant)), Abs(F2(x, y, variant)));
                PrintStr(iterCount + ":   " + x.ToString("0.000000") + "    " + y.ToString("0.000000") + "    " + error.ToString("0.000000") + "    " + pogreshnosty.ToString("0.000000"));
            } while (Abs(F1(x, y, variant)) >= Eps || Abs(F2(x, y, variant)) >= Eps);
        }

        //функция fi(x) для МПИ.
        private static double SIMethodFiX(double x, double y, int variant)
        {
            switch (variant)
            {
                case 16:
                    return -Cos(y + 0.5) + 0.8;
                case 17:
                    return 1.3 - Sin(y - 1);
            }
            return 0;
        }

        //функция fi(y) для МПИ.
        private static double SIMethodFiY(double x, double y, int variant)
        {
            switch (variant)
            {
                case 16:
                    return (Sin(x) - 1.6) / 2;
                case 17:
                    return Sin(x + 1) + 0.8;
            }
            return 0;
        }

        //метод Ньютона.
        private static void NewtonMethod(double x, double y, int variant)
        {
            var iterCount = 0;
            PrintStr("№    x           y           Невязка     Погрешность");
            do
            {
                iterCount++;
                var dx = Calc2X2MatrixDeterminant(F1(x, y, variant), Dy_F1(x, y, variant),
                    F2(x, y, variant), Dy_F2(x, y, variant));
                var dy = Calc2X2MatrixDeterminant(Dx_F1(x, y, variant), F1(x, y, variant),
                    Dx_F2(x, y, variant), F2(x, y, variant));
                var d = Calc2X2MatrixDeterminant(Dx_F1(x, y, variant), Dy_F1(x, y, variant),
                    Dx_F2(x, y, variant), Dy_F2(x, y, variant));
                var newX = x - dx / d;
                var newY = y - dy / d;
                var pogreshnosty = Abs(newX - x);
                x = newX;
                y = newY;
                var error = Max(Abs(F1(x, y, variant)), Abs(F2(x, y, variant)));
                PrintStr(iterCount + ":   " + x.ToString("0.000000") + "    " + y.ToString("0.000000") + "    " + error.ToString("0.000000") + "    " + pogreshnosty.ToString("0.000000"));
            } while (Abs(F1(x, y, variant)) >= Eps || Abs(F2(x, y, variant)) >= Eps);
        }

        //считает определитель матрицы 2x2.
        private static double Calc2X2MatrixDeterminant(double el11, double el12, double el21, double el22)
        {
            return el11 * el22 - el12 * el21;
        }

        private static double V16Fsum(double x, double y)
        {
            return Pow(F1(x, y, 16), 2) + Pow(F2(x, y, 16), 2);
        }

        private static double V17Fsum(double x, double y)
        {
            return Pow(F1(x, y, 17), 2) + Pow(F2(x, y, 17), 2);
        }

        private static double Dx_V16Fsum(double x, double y)
        {
            return 2 * Cos(y + 0.5) + 2 * x - 1.6 + Sin(2 * x) - Cos(x) * (4 * y + 3.2);
        }

        private static double Dy_V16Fsum(double x, double y)
        {
            return Sin(y + 0.5) * (-2 * Cos(y + 0.5) - 2 * x  + 1.6) - 4 * Sin(x) + 8 * y + 6.4;
        }

        private static double Dx_V17Fsum(double x, double y)
        {
            return 2 * Sin(y - 1) + 2 * x - 2.6 - 2 * y * Cos(x + 1) + Sin(2 * (x + 1)) + 1.6 * Cos(x + 1);
        }

        private static double Dy_V17Fsum(double x, double y)
        {
            return Sin(2 * (y - 1)) + Cos(y - 1) * (2 * x - 2.6) + 2 * y - 2 * Sin(x + 1) - 1.6;
        }

        //метод градиентного спуска.
        private static void GradientMethod(double x, double y, int variant)
        {
            D varF = V16Fsum;
            var gradient = new List<D>();
            var gradVector = new double[2];
            switch (variant)
            {
                case 16: varF = V16Fsum; gradient.Add(Dx_V16Fsum); gradient.Add(Dy_V16Fsum); break;
                case 17: varF = V17Fsum; gradient.Add(Dx_V17Fsum); gradient.Add(Dy_V17Fsum); break;
            }
            var startX = x;
            var startY = y;
            var minIter = int.MaxValue;
            var bestLamba = 0.02;
            int iterCount = 0;
            var alpha = 2.0;
            var error = 0.0;
            for (double lambda = 0.02; lambda < 1; lambda += 0.02)
            {
                x = startX;
                y = startY;
                iterCount = 0;
                alpha = 2.0;
                do
                {
                    alpha = 2.0;
                    iterCount++;
                    while (varF((x - alpha * gradient[0](x, y)), (y - alpha * gradient[1](x, y))) >= varF(x, y))
                    {
                        alpha = alpha * lambda;
                    }
                    var newX = x - alpha * gradient[0](x, y);
                    var newY = y - alpha * gradient[1](x, y);
                    x = newX;
                    y = newY;
                    gradVector[0] = gradient[0](x, y);
                    gradVector[1] = gradient[1](x, y);
                    error = CalcVectorNorm(gradVector, 2);
                }
                while (error > Eps);
                if(iterCount < minIter)
                {
                    minIter = iterCount;
                    bestLamba = lambda;
                }
            }
            x = startX;
            y = startY;
            iterCount = 0;
            PrintStr("№    x           y           Невязка     Погрешность  lambda");
            do
            {
                alpha = 2.0;
                iterCount++;
                while (varF((x - alpha * gradient[0](x, y)), (y - alpha * gradient[1](x, y))) >= varF(x, y))
                {
                    alpha = alpha * bestLamba;
                }
                var newX = x - alpha * gradient[0](x, y);
                var newY = y - alpha * gradient[1](x, y);
                var pogreshnosty = Abs(newX - x);
                x = newX;
                y = newY;
                gradVector[0] = gradient[0](x, y);
                gradVector[1] = gradient[1](x, y);
                error = CalcVectorNorm(gradVector, 2);
                PrintStr(iterCount + ":   " + x.ToString("0.000000") + "    " + y.ToString("0.000000") + "    " + error.ToString("0.000000") + "    "  + pogreshnosty.ToString("0.000000") + "     " + bestLamba.ToString("0.00"));
            } while (error > Eps);

        }

        private static void PrintStr(string str)
        {
            w.WriteLine(str);
            Console.WriteLine(str);
        }
        private static void PrintStr()
        {
            w.WriteLine();
            Console.WriteLine();
        }

        private static double CalcVectorNorm(double[] vector, int n)
        {
            double sum = 0;
            for(int i = 0; i < n; i++)
            {
                sum += Pow(vector[i], 2);
            }
            return Sqrt(sum);
        }
    }
}
