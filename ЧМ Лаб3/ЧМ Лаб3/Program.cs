using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ЧМ_Лаб3
{
    class Program
    {
        const double eps = 0.0001;
        static void Main(string[] args)
        {
            double[,] a = null;

        }

        static void iteration(double[,] a, double[] b, double[] x, int n)
        {
            int i, j;
            double norma; //чебышевская норма вектора
            double[] xn = { 0 };//вектор для текущей итерации, начальное значение
            do
            {
                norma = 0.0;
                for (i = 0; i < n; i++)
                {
                    xn[i] = -b[i];

                    for (j = 0; j < n; j++)
                    {
                        if (i != j)
                            xn[i] += a[i, j] * x[j];
                    }

                    xn[i] /= -a[i, i];
                }
                for (i = 0; i < n; i++)
                {
                    if (Math.Abs(x[i] - xn[i]) > norma)
                        norma = Math.Abs(x[i] - xn[i]); //Вычисление нормы вектора
                    x[i] = xn[i];
                }
            }
            while (norma > eps); //проверка на необходимую точность вычислений
            return;
        }

        static void readMatrix(ref double[,] a, ref double[] b, ref int n)
        {
            StreamReader r = new StreamReader("in.txt");
            n = Convert.ToInt32(r.ReadLine());
            a = new double[n, n];
            b = new double[n];
            string str;
            string[] t;
            for (int i = 0; i < n; i++)
            {
                str = r.ReadLine();
                t = str.Split(' ');
                for (int j = 0; j < n; j++)
                {
                    a[i, j] = Convert.ToDouble(t[j]);
                }
            }
            r.ReadLine();
            str = r.ReadLine();
            t = str.Split(' ');
            for (int j = 0; j < n; j++)
            {
                b[j] = Convert.ToDouble(t[j]);
            }
            r.Close();
        }
    }
}
