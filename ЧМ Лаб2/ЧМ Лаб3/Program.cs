using System;
using System.IO;
using System.Collections.Generic;
 
namespace LUGauss
{
    class Program
    {
        static void Main(string[] args)
        {
            //размерность квадратной матрицы.
            int n = 0;
            double[,] a = null;
            double[] b = null;
            //матрица А.
            readMatrix(ref a, ref b, ref n);
            //сохраним начальное состояние матрицы А и вектора b в отдельные переменные...
            //...во избежание потери этих значений после перестановок строк.
            double[,] initialA = a.Clone() as double[,];
            double[] initialB = b.Clone() as double[];
            //печать начального значения А и b.
            Console.WriteLine("Матрица А:");
            printMatrix(a, n);
            Console.WriteLine("Вектор b:");
            printVector(b, n);
            //выполняем перестановку строк в матрице А, генерируем матрицу перестановок P.
            double[,] p = new double[n, n];
            shiftA(a, p, n);
            Console.WriteLine("Матрица А после перестановок:");
            printMatrix(a, n);
            //выполняем соответствующую перестановку в векторе b через его умножение на P.
            b = prodMatrixVector(p, b, n);
            Console.WriteLine("Вектор b после перестановок:");
            printVector(b, n);
            //создаем матрицы L и U.
            double[,] u = a.Clone() as double[,];
            double[,] l = a.Clone() as double[,];
            //считаем их и выводим на экран.
            calcLU(u, l, n);
            Console.WriteLine("Матрица L:");
            printMatrix(l, n);
            Console.WriteLine("Матрица U:");
            printMatrix(u, n);
            //считаем и выводим на экран промежуточный вектор y.
            var y = calcY(l, b, n);
            Console.WriteLine("Вектор y:");
            printVector(y, n);
            //аналогично для искомого вектора x.
            var x = calcX(u, y, n);
            Console.WriteLine("Вектор x:");
            printVector(x, n);
            //считаем и выводим на экран матрицу, обратную А, используя полученное ранее LU-разложение.
            var reverseA = calcReverseA(l, u, n);
            Console.WriteLine("Обратная матрица А:");
            printMatrix(reverseA, n);
            //считаем и выводим на экран определитель А + ранг.
            var determinant = calcDeterminant(l, n);
            Console.WriteLine("Определитель матрицы А: " + determinant.ToString("0.000"));
            if (determinant >= 10e-3 || determinant <= -10e-3)
            {
                Console.WriteLine("Ранг матрицы А: " + n);
            }
            else
            {
                Console.WriteLine("Ранг матрицы А меньше чем " + n + " => Матрица А вырожденная.");
            }
            //считаем и выводим на экран число обусловленности А.
            var conditionNumber = calcNorm(a, n) * calcNorm(reverseA, n);
            Console.WriteLine("Число обусловленности матрицы А: " + conditionNumber.ToString("0.000"));
            //выполняем все необходимые проверки.
            Console.WriteLine();
            Console.WriteLine("Проверки:");
            Console.WriteLine("1. P * A = L * U");
            Console.WriteLine("Матрица P * A:");
            printMatrix(prodMatrixMatrix(p, initialA, n), n);
            Console.WriteLine("Матрица L * U:");
            printMatrix(prodMatrixMatrix(l, u, n), n);
            Console.WriteLine("2. A * x = b");
            Console.WriteLine("Начальное значение b:");
            printVector(initialB, n);
            Console.WriteLine("Посчитанное из A * x значение b:");
            printVector(prodMatrixVector(initialA, x, n), n);
            Console.WriteLine("3. A * A^(-1) = E");
            Console.WriteLine("Произведение А и А^(-1):");
            printMatrix(prodMatrixMatrix(a, reverseA, n), n);
        }
 
        //чтение матрицы
        static void readMatrix(ref double[,] a, ref double[] b, ref int n)
        {
            StreamReader r = new StreamReader("in.txt");
            n = Convert.ToInt32(r.ReadLine());
            a = new double[n, n];
            b = new double[n];
            string str;
            string[] t;
            for (int i = 0; i<n; i++)
            {
                str = r.ReadLine();
	                t = str.Split(' ');
	                for(int j = 0; j<n; j++)
	                {
	                    a[i, j] = Convert.ToDouble(t[j]);
	                }
	            }
	            r.ReadLine();
	            str = r.ReadLine();
	            t = str.Split(' ');
	            for (int j = 0; j<n; j++)
	            {
	                b[j] = Convert.ToDouble(t[j]);
	            }
	        }
	 
	        //печать матрицы.
	        static void printMatrix(double[,] matrix, int n)
	        {
	            for (int i = 0; i<n; i++)
	            {
	                for (int j = 0; j<n; j++)
	                    Console.Write(matrix[i, j].ToString("0.000") + "\t");
	                Console.WriteLine();
	            }
	            Console.WriteLine();
	        }
	 
	        //печать вектора.
	        static void printVector(double[] vector, int n)
	        {
	            for (int i = 0; i<n; i++)
	            {
	                Console.Write(vector[i].ToString("0.000") + "\t");
	            }
	            Console.WriteLine();
	            Console.WriteLine();
	        }
	 
	        //реализует перестановку строк матрицы А с выбором главного элемента по столбцу.
	        static void shiftA(double[,] a, double[,] p, int n)
	        {
	            for (int i = 0; i<n; i++)
	                p[i, i] = 1;
	            for (int j = 0; j<n; j++)
	            {
	                int indexMax = -1;
	                double numberMax = a[j, j];
	                //ищем строку с максимальным элементом
	                for (int i = j + 1; i<n; i++)
	                {
	                    if (a[i, j] > numberMax)
	                    {
	                        indexMax = i;
	                        numberMax = a[i, j];
	                    }
	                }
	                //перставляем строку, если необходимо
	                if (indexMax != -1)
	                {
	                    Console.WriteLine("Переставляем строки {0} и {1}:", j + 1, indexMax + 1);
	                    swap(a, n, j, indexMax);
	                    swap(p, n, j, indexMax);
	                    printMatrix(a, n);
	                }
	            }
	 
	        }
	 
	        //реализует перестановку строк в матрице по заданным индексам.
	        static void swap(double[,] a, int n, int first, int second)
	        {
	            for (int i = 0; i<n; i++)
	            {
	                double temp = a[first, i];
	                a[first, i] = a[second, i];
	                a[second, i] = temp;
	            }
	        }
	 
	        //раскладывает матрицу A на L * U.
	        static void calcLU(double[,] u, double[,] l, int n)
	        {
	            for (int i = 0; i<n; i++)
	            {
	                for (int j = 0; j<n; j++)
	                {
	                    double sum = 0;
	                    //считаем L.
	                    if (i >= j)
	                    {
	                        if (i == j)
	                        {
	                            u[i, j] = 1;
	                        }
	                        else
	                        {
	                            u[i, j] = 0;
	                        }
	                        for (int s = 0; s <= j - 1; s++)
	                            sum += l[i, s] * u[s, j];
	                        l[i, j] -= sum;
	                    }
	                    //считаем U.
	                    else
	                    {
	                        if (i != j)
	                        {
	                            l[i, j] = 0;
	                        }
	                        for (int s = 0; s <= i - 1; s++)
	                            sum += l[i, s] * u[s, j];
	                        u[i, j] = (u[i, j] - sum) / l[i, i];
	                    }
	                }
	            }
	        }
	 
	        //U * x = y => L * y = b. Ищет y:
	        static double[] calcY(double[,] l, double[] b, int n)
	        {
	            var y = new double[n];
	            for (int i = 0; i<n; i++)
	            {
	                double t = b[i];
	                for (int j = 0; j<n; j++)
	                {
	                    if (i == j)
	                    {
	                        continue;
	                    }
	                    t -= l[i, j] * y[j];
	                }
	                y[i] = t / l[i, i];
	            }
	            return y;
	        }
	 
	        //U * x = y. Ищет x:
	        static double[] calcX(double[,] u, double[] y, int n)
	        {
	            var x = new double[n];
	            for (int i = n - 1; i >= 0; i--)
	            {
	                double t = y[i];
	                for (int j = n - 1; j >= 0; j--)
	                {
	                    if (i == j)
	                    {
	                        continue;
	                    }
	                    t -= u[i, j] * x[j];
	                }
	                x[i] = t / u[i, i];
	            }
	            return x;
	        }
	 
	        //считает определитель матрицы A.
	        static double calcDeterminant(double[,] l, int n)
	        {
	            double determinant = 1;
	            for (int i = 0; i<n; i++)
	            {
	                determinant *= l[i, i];
	            }
	            return determinant;
	        }
	 
	        //считает норму матрицы.
	        static double calcNorm(double[,] matrix, int n)
	        {
	            double currMax = 0;
	            for (int i = 0; i<n; i++)
	            {
	                double sum = 0;
	                for (int j = 0; j<n; j++)
	                {
	                    sum += Math.Abs(matrix[i, j]);
	                }
	                if (sum > currMax)
	                {
	                    currMax = sum;
	                }
	            }
	            return currMax;
	        }
	 
	        //считает обратную матрице А.
	        //L * y = e
	        //U * x = y
	        static double[,] calcReverseA(double[,] l, double[,] u, int n)
	        {
	            List<double[]> y = new List<double[]>();
	            //вычисление вектора y
	            for (int i = 0; i<n; i++)
	            {
	                double[] t = new double[n];
	                t[i] = 1;
	                y.Add(calcY(l, t, n));
	            }
	            double[,] x = new double[n, n];
	            //вычисление вектора x
	            for (int i = 0; i<n; i++)
	            {
	                double[] t = calcX(u, y[i], n);
	                for (int j = 0; j<n; j++)
	                    x[j, i] = t[j];
	            }
	            return x;
	        }
	 
	        //считает произведение двух матриц.
	        static double[,] prodMatrixMatrix(double[,] a, double[,] b, int n)
	        {
	            double[,] prodM = new double[n, n];
	            for (int i = 0; i<n; i++)
	            {
	                for (int j = 0; j<n; j++)
	                {
	                    double sum = 0;
	                    for (int k = 0; k<n; k++)
	                    {
	                        sum += a[i, k] * b[k, j];
	                    }
	                    prodM[i, j] = sum;
	                }
	            }
	            return prodM;
	        }
	 
	        //считает произведение матрицы и вектора.
	        static double[] prodMatrixVector(double[,] a, double[] b, int n)
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
	    }
	}
