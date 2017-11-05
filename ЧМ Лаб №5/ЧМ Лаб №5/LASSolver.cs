using System;

namespace ЧМ_Лаб__5
{
    static class LASSolver
    {
        public static void ConjurateGradientsMethod(double[,] a, double[] b, int n, ref double[] x)
        {
            double discrepancyNorm = 0;
            var d = new double[n]; //вектор направления.
            var g = new double[n]; //вектор градиента.
            double s = 0; //скалярный шаг.
            var prevG = NegationVection(b, n);
            do
            {
                g = DiffVector(b, ProdMatrixVector(a, x, n), n);
                d = SumVectors(NegationVection(g, n), ProdVectorNumber(d, ScolarMult(g, g, n) / ScolarMult(prevG, prevG, n), n), n);
                s = ScolarMult(d, g, n) / ScolarMult(ProdMatrixVector(a, d, n), d, n);
                x = SumVectors(x, ProdVectorNumber(d, s, n), n);
                prevG = g.Clone() as double[];
                discrepancyNorm = CalcVectorNorm(DiffVector(b, ProdMatrixVector(a, x, n), n), n);
            } while (discrepancyNorm > 0.0001);
        }

        private static double[] NegationVection(double[] v, int n)
        {
            var answer = new double[n];
            for (int i = 0; i < n; i++)
            {
                answer[i] = -v[i];
            }
            return answer;
        }

        private static double CalcVectorNorm(double[] vector, int n)
        {
            double sum = 0;
            for (int i = 0; i<n; i++)
            {
                sum += vector[i] * vector[i];
            }
            return Math.Sqrt(sum);
        }

        private static double[] ProdMatrixVector(double[,] a, double[] b, int n)
        {
            double[] prodM = new double[n];
            for (int j = 0; j<n; j++)
            {
                double sum = 0;
                for (int k = 0; k<n; k++)
                {
                    sum += a[j, k] * b[k];
                }
                prodM[j] = sum;
            }
            return prodM;
        }

        private static double[] DiffVector(double[] a, double[] b, int n)
        {
            var difference = new double[n];
            for (int i = 0; i<n; i++)
            {
                difference[i] = a[i] - b[i];
            }
            return difference;
        }

        private static double[] SumVectors(double[] a, double[] b, int n)
        {
            var sum = new double[n];
            for (int i = 0; i<n; i++)
            {
                sum[i] = a[i] + b[i];
            }
            return sum;
        }

        private static double ScolarMult(double[] v1, double[] v2, int n)
        {
            double sum = 0;
            for (int i = 0; i<n; i++)
            {
                sum += v1[i] * v2[i];
            }
            return sum;
        }

        private static double[] ProdVectorNumber(double[] v, double numb, int n)
        {
            var answer = new double[n];
            for (int i = 0; i<n; i++)
            {
                answer[i] = v[i] * numb;
            }
            return answer;
        }

        public static void PrintMatrix(double[,] matrix, int n)
        {
            for (int i = 0; i<n; i++)
            {
                for (int j = 0; j<n; j++)
                    Console.Write(matrix[i, j].ToString("0.000") + "\t");
                Console.WriteLine();
            }
            Console.WriteLine();
        }

        public static void PrintVector(double[] vector, int n)
        {
            for (int i = 0; i<n; i++)
            {
                Console.Write(vector[i].ToString("0.000") + "\t");
            }
            Console.WriteLine();
            Console.WriteLine();
        }

    }
}
