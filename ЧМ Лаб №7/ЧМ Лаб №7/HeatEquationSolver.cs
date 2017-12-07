using System;

namespace ЧМ_Лаб__7
{
    internal static class HeatEquationSolver
    {
        private static double[,] u;
        private const double X = 0.1;

        private static double Fi(double x)
        {
            return x;
        }

        private static double F(int var, double t, double x)
        {
            return 0.1 * var * Math.Sin(Math.PI * x) + 0.1 * X * Math.PI * Math.PI * t * var * Math.Sin(Math.PI * x);
        }

        private static double U(int var, double t, double x)
        {
            return x + 0.1 * t * Math.Sin(Math.PI * x) * var;
        }

        private static double CalcError(double uValue, double t, double x)
        {
            return Math.Abs(U(Program.Var, t, x) - uValue);
        }

        //Явная конечно-разностная схема.
        public static void SolveWithExplicitMesh(int N)
        {
            var h = 1.0 / N;
            var tau = h * h / 4 / X;
            var m = (int) Math.Ceiling(1 / tau) + 1;
            u = new double[m + 1, N + 1];
            double x;
            for (var j = 0; j <= N; j++)
            {
                x = j * h;
                u[0, j] = Fi(x);
            }
            double t = 0;
            for (var i = 0; t <= 1; i++)
            {
                t = i * tau;
                u[i, 0] = 0;
                u[i, N] = 1;
            }
            t = 0;
            Program.PrintStr("===== Явная разностная схема при N = " + N + " =====");
            var totalMaxError = double.MinValue;
            var n = 0;
            while(t <= 1)
            {
                t = n * tau;
                Program.PrintStrNoLineBreak("t = " + t.ToString("0.00000"));
                var maxError = double.MinValue;
                for (var j = 1; j < N; j++)
                {
                    x = j * h;
                    if (CalcError(u[n, j], t, x) > maxError)
                        maxError = CalcError(u[n, j], t, x);
                    u[n + 1, j] = u[n, j] + tau * X * (u[n, j + 1] - 2 * u[n, j] + u[n, j - 1]) / h / h +
                                  tau * F(Program.Var, t, x);
                }
                if (maxError > totalMaxError)
                    totalMaxError = maxError;
                Program.PrintStr(" Delta = " + maxError);
                n++;
            }
            Program.PrintStr("Ответ:");
            n--;
            for (var j = 0; j <= N; j++)
            {
                Program.PrintStrNoLineBreak(u[n, j].ToString("0.00000") + " ");
            }
            Program.PrintStr();
            Program.PrintStr("Max delta = " + totalMaxError);
            Program.PrintStr();
        }

        //Неявная конечно-разностная схема.
        public static void SolveWithImplicitMesh(int N)
        {
            var h = 1.0 / N;
            var tau = h;
            var m = (int)Math.Ceiling(1 / tau) + 1;
            u = new double[m + 1, N + 1];
            double x;
            for (var j = 0; j <= N; j++)
            {
                x = j * h;
                u[0, j] = Fi(x);
            }
            double t = 0;
            for (var i = 0; t <= 1; i++)
            {
                t = i * tau;
                u[i, 0] = 0;
                u[i, N] = 1;
            }
            t = 0;
            Program.PrintStr("===== Неявная разностная схема при N = " + N + " =====");
            var d = tau * X / h / h;
            var above = new double[N];
            var main = new double[N + 1];
            var below = new double[N + 1];
            var right = new double[N + 1];
            var totalMaxError = double.MinValue;
            var n = 0;
            while (t < 1)
            {
                t = n * tau;
                Program.PrintStrNoLineBreak("t = " + t.ToString("0.00000") + " ");
                var maxError = double.MinValue;
                for (var j = 0; j <= N; j++)
                {
                    x = j * h;
                    if (CalcError(u[n, j], t, x) > maxError)
                        maxError = CalcError(u[n, j], t, x);
                    var nextT = (n + 1) * tau;
                    right[j] = u[n, j] + tau * F(Program.Var, nextT, x);
                }
                InitProgonkaArrays(above, main, below, N, d);
                var uj = SolveTriDiagMatrixSystem(above, main, below, right, N);
                for (var j = 1; j < N; j++)
                {
                    u[n + 1, j] = uj[j];
                }
                if (maxError > totalMaxError)
                    totalMaxError = maxError;
                Program.PrintStr("Delta = " + maxError);
                n++;
            }
            n--;
            Program.PrintStr("Ответ:");
            for (var j = 0; j <= N; j++)
            {
                Program.PrintStrNoLineBreak(u[n, j].ToString("0.00000") + " ");
            }
            Program.PrintStr();
            Program.PrintStr("Max delta = " + totalMaxError);
            Program.PrintStr();
        }

        // Реализует алгоритм прогонки.
        private static double[] SolveTriDiagMatrixSystem(double[] above, double[] main, double[] below, double[] right, int N)
        {
            var answer = new double[N + 1];
            for (var i = 1; i <= N; i++)
            {
                var m = below[i] / main[i - 1];
                main[i] = main[i] - m * above[i - 1];
                right[i] = right[i] - m * right[i - 1];
            }
            answer[N + 1 - 1] = right[N + 1 - 1] / main[N + 1 - 1];
            for (var i = N + 1 - 2; i >= 0; i--)
                answer[i] = (right[i] - above[i] * answer[i + 1]) / main[i];
            return answer;
        }

        private static void InitProgonkaArrays(double[] above, double[] main, double[] below, int N, double d)
        {
            
            for (var j = 1; j < N; j++)
            {
                above[j] = below[j] = -d;
                main[j] = 1 + 2 * d;
            }
            //below[0] = 0;
            //above[N - 1] = 0;
            main[0] = main[N] = 1;
        }
    }
}
