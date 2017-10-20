using System;

namespace ЧМ_Лаб_4
{
    internal static class Program
    {
        private const double Eps = 0.0001;
        //вариант 16а.
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

        private static double Dx_F1(double x, double y, int variant)
        {
            return 1;
        }

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
            //NewtonMethod.
            var x = -0.25;
            var y = -1;
            NewtonMethod(x, y, 16);
        }

        private static void SimpleIterMethod()
        {
            
        }

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
            } while (Math.Abs(F1(x, y, variant)) >= Eps && Math.Abs(F2(x, y, variant)) >= Eps);
            Console.WriteLine("x: " + x + " y: " + y + " " + iterCount + " iterations");
        }

        private static double Calc2X2MatrixDeterminant(double el11, double el12, double el21, double el22)
        {
            return el11 * el22 - el12 * el21;
        }

        private static double CalcXYNorm(double x, double prevX, double y, double prevY)
        {
            return Math.Sqrt(Math.Pow(x - prevX, 2) + Math.Pow(y - prevY, 2));
        }
    }
}
