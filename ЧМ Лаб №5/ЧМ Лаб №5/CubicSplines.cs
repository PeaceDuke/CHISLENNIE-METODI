using System;
using static ЧМ_Лаб__5.Program;

namespace ЧМ_Лаб__5
{
    class CubicSplines
    {
        public static double[] _splineParams;
        private static double M4 = double.MinValue;
        private static double M5 = double.MinValue;

        // Считает параметры сплайна.
        public static void CalcSplineParams()
        {
            _splineParams = new double[N + 1];
            var b = new double[N + 1];
            _splineParams[0] = b[0] = FirstDerivative(A);
            _splineParams[N] = b[N] = FirstDerivative(B);
            for (var i = 1; i < N; i++)
            {
                const double u = 1 / 2.0;
                const double lambda = 1 / 2.0;
                b[i] = 2 * (3 * lambda * (Table[i].Y - Table[i - 1].Y) / Delta + 3 * u * (Table[i + 1].Y - Table[i].Y) / Delta);
            }
            var triDiagMatrix = new double[N + 1, N + 2];
            triDiagMatrix[0, 0] = triDiagMatrix[N, N] =  1;
            var j = 1;
            for (var i = 1; i < N; i++)
            {
                triDiagMatrix[i, j] = 4;
                triDiagMatrix[i, j - 1] = triDiagMatrix[i, j + 1] = 1;
                j++;
            }
            for (var i = 0; i <= N; i++)
            {
                triDiagMatrix[i, N + 1] = b[i];
            }
            PrintTriDiagMatrix(triDiagMatrix);
            var above = new double[N];
            for (var i = 1; i < N; i++)
                above[i] = 1;
            var main = new double[N + 1];
            main[0] = main[N] = 1;
            for (var i = 1; i < N; i++)
                main[i] = 4;
            var below = new double[N + 1];
            for (var i = 1; i < N; i++)
                below[i] = 1;
            var right = b;
            _splineParams = SolveTriDiagMatrixSystem(above, main, below, right);
            PrintSplineParams();
        }

        // Печатает трехдиагональную матрицу.
        private static void PrintTriDiagMatrix(double[,] matrix)
        {
            Console.WriteLine("Трехдиагональная матрица:");
            for (var i = 0; i <= N; i++)
            {
                for (var j = 0; j <= N + 1; j++)
                {
                    Console.Write(matrix[i, j] + " ");
                    OutFile.Write(matrix[i, j] + " ");
                }
                PrintStr();
            }
            PrintStr();
        }

        // Печатает параметры сплайна.
        private static void PrintSplineParams()
        {
            PrintStr("Параметры сплайна:");
            for (var i = 0; i < N + 1; i++)
            {
                PrintStr("m" + i + ": " + _splineParams[i]);
            }
            PrintStr();
        }

        // Реализует алгоритм прогонки.
        private static double[] SolveTriDiagMatrixSystem(double[] above, double[] main, double[] below, double[] right)
        {
            var answer = new double[N + 1];
            for (var i = 1; i < N + 1; i++)
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

        // Считает значение сплайна в указанной точке на [1; 2]
        public static double CalcSplineValue(double x)
        {
            var i = 0;
            while (Table[i].X < x) i++;
            i = Table[i].X == x ? i : i - 1;
            var xi = Table[i].X;
            var t = (x - xi) / Delta;
            var fi0 = (1 - t) * (1 - t) * (1 + 2 * t);
            var fi1 = (1 - t) * (1 - t) * t;
            return (fi0 * t * Table[i].Y + fi0 * (1 - t) * Table[i + 1].Y + Delta * fi1 * t * _splineParams[i] -
                   Delta * fi1 * t * _splineParams[i + 1]) * 2;
        }

        // Считает оценку погрешности.
        public static double CalcTeoreticError()
        {
            for (var i = 1.0; i <= 2.0; i += Delta)
            {
                if (FourthDerivative(i) > M4)
                    M4 = FourthDerivative(i);
                if (FifthDerivative(i) > M5)
                    M5 = FifthDerivative(i);
            }
            
            return (M4 / 384 + M5 * Delta / 240) * Math.Pow(Delta, 4);
        }

        // Печатает значения M4 и M5.
        public static void PrintM4M5()
        {
            PrintStr("M4: " + M4);
            PrintStr("M5: " + M5);
            PrintStr();
        }
    }
}
