using System;

namespace ЧМ_Лаб_4
{
    internal static class Program
    {
        private const double Eps = 0.0001;

        //функция f1.
        private static double F1(double x, double y, int variant)
        {
            switch (variant)
            {
                case 16:
                    return Math.Cos(y + 0.5) + x - 0.8;
                case 17:
                    return Math.Sin(y - 1) + x - 1.3;
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
                    return -Math.Sin(y + 0.5);
                case 17:
                    return Math.Cos(y - 1);
            }
            return 0;
        }

        //функция f2.
        private static double F2(double x, double y, int variant)
        {
            switch (variant)
            {
                case 16:
                    return Math.Sin(x) - 2 * y - 1.6;
                case 17:
                    return y - Math.Sin(x + 1) - 0.8;
            }
            return 0;
        }

        //производная f2 по x.
        private static double Dx_F2(double x, double y, int variant)
        {
            switch (variant)
            {
                case 16:
                    return Math.Cos(x);
                case 17:
                    return Math.Cos(x + 1);
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
            //SimpleIterationMethod.
            double x = -0.25;
            double y = -1;
            SimpleIterMethod(x, y, 16);
            //NewtonMethod.
            x = -0.25;
            y = -1;
            NewtonMethod(x, y, 16);
        }

        //метод простой итерации.
        private static void SimpleIterMethod(double x, double y, int variant)
        {
            var iterCount = 0;
            do
            {
                iterCount++;
                x = SIMethodFiX(x, y, 16);
                y = SIMethodFiY(x, y, 16);
            } while (Math.Abs(F1(x, y, variant)) >= Eps || Math.Abs(F2(x, y, variant)) >= Eps);
            Console.WriteLine("x: " + x + " y: " + y + " " + iterCount + " iterations");
        }

        //функция fi(x) для МПИ.
        private static double SIMethodFiX(double x, double y, int variant)
        {
            switch (variant)
            {
                case 16:
                    return -Math.Cos(y + 0.5) + 0.8;
                case 17:
                    return 1.3 - Math.Sin(y - 1);
            }
            return 0;
        }

        //функция fi(y) для МПИ.
        private static double SIMethodFiY(double x, double y, int variant)
        {
            switch (variant)
            {
                case 16:
                    return (Math.Sin(x) - 1.6) / 2;
                case 17:
                    return Math.Sin(x + 1) + 0.8;
            }
            return 0;
        }

        //метод Ньютона.
        private static void NewtonMethod(double x, double y, int variant)
        {
            var iterCount = 0;
            do
            {
                iterCount++;
                var dx = Calc2X2MatrixDeterminant(F1(x, y, variant), Dy_F1(x, y, variant),
                    F2(x, y, variant), Dy_F2(x, y, variant));
                var dy = Calc2X2MatrixDeterminant(Dx_F1(x, y, variant), F1(x, y, variant),
                    Dx_F2(x, y, variant), F2(x, y, variant));
                var d = Calc2X2MatrixDeterminant(Dx_F1(x, y, variant), Dy_F1(x, y, variant),
                    Dx_F2(x, y, variant), Dy_F2(x, y, variant));
                x = x - dx / d;
                y = y - dy / d;
            } while (Math.Abs(F1(x, y, variant)) >= Eps || Math.Abs(F2(x, y, variant)) >= Eps);
            Console.WriteLine("x: " + x + " y: " + y + " " + iterCount + " iterations");
        }

        //считает определитель матрицы 2x2.
        private static double Calc2X2MatrixDeterminant(double el11, double el12, double el21, double el22)
        {
            return el11 * el22 - el12 * el21;
        }

        //считает норму вектора (x, y).
        private static double CalcXYNorm(double x, double prevX, double y, double prevY)
        {
            return Math.Sqrt(Math.Pow(x - prevX, 2) + Math.Pow(y - prevY, 2));
        }
    }
}
