using System;
using System.IO;
using static System.Math;

namespace ЧМ_Лаб__7
{
    class Program
    {
        public const int Var = 16;
        const double Y0 = 1, Y1 = 2;
        static StreamWriter outFile = new StreamWriter("out.txt");

        private static double TrialFunc(double x) => 1 + x + 10 * Log(Var + 1.0) * Pow(x, 3) * Pow(1 - x, 3);

        private static double FirstDerivativeTrialFunc(double x) => -30 * Pow(1 - x, 2) * Pow(x, 3) * Log(Var + 1.0) 
            - 30 * Pow(x, 2) * Pow(x - 1, 3) * Log(Var + 1.0) + 1;

        private static double SecondDerivativeTrialFunc(double x) => -60 * x * (5 * Pow(x, 3) - 10 * Pow(x, 2) + 6 * x - 1) * Log(Var + 1.0);
        
        private static double ShootingFunc(double x, double y, double z)
        {
            return RightFunc(x) - AFunc(x) * z + BFunc(x) * y - CFunc(x) * Sin(y);
        }

        private static double AFunc(double x)
        {
            switch (Var)
            {
                case 16: return 30 * (-x + 0.5);
                case 17: return 40 * (x + 1);
                default: return 0;
            }
        }

        private static double BFunc(double x)
        {
            switch (Var)
            {
                case 16: return Pow(x, 2) + 1;
                case 17: return Pow(x, 2) + 2;
                default: return 0;
            }
        }

        private static double CFunc(double x)
        {
            switch (Var)
            {
                case 16: return x + 2;
                case 17: return x + 1;
                default: return 0;
            }
        }

        private static double RightFunc(double x) => SecondDerivativeTrialFunc(x) + AFunc(x) * FirstDerivativeTrialFunc(x) - BFunc(x) 
            * TrialFunc(x) + CFunc(x) * Sin(TrialFunc(x));

        //private static double CalcGaussIntegral(double x)
        //{
        //    double sum = 0;
        //    int n = 20;
        //    var delta = x / n;
        //    for (var i = 0; i < n; i++)
        //    {
        //        var a = delta * i;
        //        var b = delta * (i + 1);
        //        sum += (b - a) * (5.0 / 9.0 * Func(0.5 * (a + b + (b - a) * -Sqrt(3.0 / 5.0))) +
        //                            8.0 / 9.0 * Func(0.5 * (a + b)) + 5.0 / 9.0 * Func(0.5 * (a + b + (b - a) * Sqrt(3.0 / 5.0)))) / 2;
        //    }
        //    return sum;
        //}

        static double maxError = double.MinValue;
        private static double CalcK(double delta, double x, double y, double z)
        {
            double[] k = new double[4], l = new double[4];
            if (flag) PrintStr("x \ty(x) полученное \ty(x) реальное \t\t Погрешность");
            do
            {
                l[0] = delta * ShootingFunc(x, y, z);
                k[0] = delta * z;
                l[1] = delta * ShootingFunc(x + delta / 2, y + k[0] / 2, z + l[0] / 2);
                k[1] = delta * (z + l[0] / 2);
                l[2] = delta * ShootingFunc(x + delta / 2, y + k[1] / 2, z + l[1] / 2);
                k[2] = delta * (z + l[1] / 2);
                l[3] = delta * ShootingFunc(x + delta, y + k[2], z + l[2]);
                k[3] = delta * (z + l[2]);
                if (flag)
                {
                    var realY = TrialFunc(x);
                    var error = Abs(realY - y);
                    if (maxError < error)
                        maxError = error;
                    PrintStr(x.ToString("0.00#####") + '\t' + y.ToString("0.0000000000000") + "\t\t" + 
                        TrialFunc(x).ToString("0.0000000000000") + '\t' + '\t' + error);
                }
                y = y + 1.0 / 6.0 * (k[0] + 2 * (k[1] + k[2]) + k[3]);
                x = Round(x + delta, 6);
                z = z + 1.0 / 6.0 * (l[0] + 2 * (l[1] + l[2]) + l[3]);
                zet = z;
            }
            while (x < 1);
            if (flag) PrintStr("Максимальная погрешность: " + maxError);
            return y;
        }

        static bool flag = false;
        static double d = 0.1;
        static double zet = 0;
        private static double KoshiSolve(double x, double y, double z)
        {
            double delta = 2 * 0.1, y1 = 0, y2 = 0, eps = 1e-5;
            int i = 0;
            do
            {
                i++;
                delta = delta / 2;
                y1 = CalcK(delta, 1, y, z);
                y2 = CalcK(delta / 2, 1, y, z);
                y2 = CalcK(delta / 2, 1, y2, zet);
            }
            while (Abs(y1 - y2) > eps);
            y1 = CalcK(delta, x, y, z);
            PrintStr("Alpha = " + z + " y(1) = " + y1 + " Погрешность "+ Abs(y1 - 2));
            PrintStr();
            d = delta;
            return y1;
        }

        private static void Solver()
        {
            double alpha1 = 0, alpha2 = 2.5, alpha3 = 0, b = 0, eps = 1e-4;
            if (Abs(KoshiSolve(0, 1, alpha1) - 2) < eps)
                return;
            if (Abs(KoshiSolve(0, 1, alpha2) - 2) < eps)
                return;
            while (Abs(b - 2) > eps)
            {
                alpha3 = (alpha1 + alpha2) / 2;
                b = KoshiSolve(0, 1, alpha3);
                if (b < 2)
                {
                    alpha1 = alpha3;
                }
                else
                {
                    alpha2 = alpha3;
                }
            }
            PrintStr("Искомая alpha = " + alpha3);
            flag = true;
            CalcK(0.0125, 0, 1, alpha3);
            PrintStr();
        }

        static void Main(string[] args)
        {
            Solver();
            HeatEquationSolver.SolveWithExplicitMesh(8);
            HeatEquationSolver.SolveWithExplicitMesh(16);
            HeatEquationSolver.SolveWithExplicitMesh(32);
            HeatEquationSolver.SolveWithImplicitMesh(8);
            HeatEquationSolver.SolveWithImplicitMesh(16);
            HeatEquationSolver.SolveWithImplicitMesh(32);
        }

        public static void PrintStr(string str)
        {
            Console.WriteLine(str);
            outFile.WriteLine(str);
        }

        public static void PrintStrNoLineBreak(string str)
        {
            Console.Write(str);
            outFile.Write(str);
        }

        public static void PrintStr()
        {
            Console.WriteLine();
            outFile.WriteLine();
        }
    }
}
