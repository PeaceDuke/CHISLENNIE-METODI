using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ЧМ_Лаб__5
{
    class LASSolver
    {
        public static void conjurateGradientsMethod(double[,] a, double[] b, int n, ref double[] x)
        {
            double discrepancyNorm = 0;
            var d = new double[n]; //вектор направления.
            var g = new double[n]; //вектор градиента.
            double s = 0; //скалярный шаг.
            var prevG = negationVection(b, n);
            do
            {
                g = diffVector(b, prodMatrixVector(a, x, n), n);
                d = sumVectors(negationVection(g, n), prodVectorNumber(d, scolarMult(g, g, n) / scolarMult(prevG, prevG, n), n), n);
                s = scolarMult(d, g, n) / scolarMult(prodMatrixVector(a, d, n), d, n);
                x = sumVectors(x, prodVectorNumber(d, s, n), n);
                prevG = g.Clone() as double[];
                discrepancyNorm = calcVectorNorm(diffVector(b, prodMatrixVector(a, x, n), n), n);
            } while (discrepancyNorm > 0.0001);
        }
        public static double[] negationVection(double[] v, int n)
        {
            var answer = new double[n];
            for (int i = 0; i<n; i++)
            {
                answer[i] = -v[i];
            }
            return answer;
        }
        public static double calcVectorNorm(double[] vector, int n)
        {
            double sum = 0;
            for (int i = 0; i<n; i++)
            {
                sum += vector[i] * vector[i];
            }
            return Math.Sqrt(sum);
        }
        public static double[] prodMatrixVector(double[,] a, double[] b, int n)
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
        public static double[] diffVector(double[] a, double[] b, int n)
        {
            var difference = new double[n];
            for (int i = 0; i<n; i++)
            {
                difference[i] = a[i] - b[i];
            }
            return difference;
        }
        public static double[] sumVectors(double[] a, double[] b, int n)
        {
            var sum = new double[n];
            for (int i = 0; i<n; i++)
            {
                sum[i] = a[i] + b[i];
            }
            return sum;
        }
        public static double scolarMult(double[] v1, double[] v2, int n)
        {
            double sum = 0;
            for (int i = 0; i<n; i++)
            {
                sum += v1[i] * v2[i];
            }
            return sum;
        }
        public static double[] prodVectorNumber(double[] v, double numb, int n)
        {
            var answer = new double[n];
            for (int i = 0; i<n; i++)
            {
                answer[i] = v[i] * numb;
            }
            return answer;
        }
    }
}
