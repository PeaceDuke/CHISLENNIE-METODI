using System;

namespace ЧМ_Лаб_4
{
    class Program
    {
        //вариант 16а.
        static double f1_16(double x, double y)
        {
            return Math.Cos(y + 0.5) + x - 0.8;
        }

        static double dx_f1_16(double x, double y)
        {
            return 1;
        }

        static double dy_f1_16(double x, double y)
        {
            return -Math.Sin(y + 0.5);
        }

        static double f2_16(double x, double y)
        {
            return Math.Sin(x) - 2 * y - 1.6;
        }

        static double dx_f2_16(double x, double y)
        {
            return Math.Cos(x);
        }

        static double dy_f2_16(double x, double y)
        {
            return -2;
        }

        //вариант 17а.
        static double f1_17(double x, double y)
        {
            return Math.Sin(y - 1) + x - 1.3;
        }

        static double dx_f1_17(double x, double y)
        {
            return 1;
        }

        static double dy_f1_17(double x, double y)
        {
            return Math.Cos(y - 1);
        }

        static double f2_17(double x, double y)
        {
            return y - Math.Sin(x + 1) - 0.8;
        }

        static double dx_f2_17(double x, double y)
        {
            return Math.Cos(x + 1);
        }

        static double dy_f2_17(double x, double y)
        {
            return 1;
        }

        static void Main(string[] args)
        {
            var variant16 = new EqationClass("cos( y + 0,5 ) + x - 0,8 = 0 // sin( x ) - 2y - 1,6 = 0", 1.0, 1.0);
            var variant17 = new EqationClass("sin( y - 1 ) + x - 1,3 = 0 // y - sin( x + 1 ) - 0,8 = 0", 1.0, 1.0);
        }

        static void SimpleIterMethod()
        {
            
        }

        static void NewtonMethod()
        {
            //дарова))
        }
    }
}
