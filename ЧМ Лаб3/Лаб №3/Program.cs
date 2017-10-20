using System;
using System.IO;
using System.Collections.Generic;

namespace CMLab3
{
    class Program
    {
        const double eps = 0.0001;

        static int iterationsCount = 0;
        static double discrepancyNorm = 0;
        static double error = 0;
        static StreamWriter w;

        static void Main(string[] args)
        {
            w = new StreamWriter("out.txt");
            int n = 0;
            double[,] a = null;
            double[] b = null;
            readMatrix(ref a, ref b, ref n);
            //печать начального значения А и b.
            printStr("Матрица А:");
            printStr();
            printMatrix(a, n);
            printStr("Вектор b:");
            printStr();
            printVector(b, n);
            printStr();
            //считаем решение СЛАУ методом LU-разложения.
            var bLU = b.Clone() as double[];
            var u = a.Clone() as double[,];
            var l = a.Clone() as double[,];
            var p = new double[n, n];
            var iterationCounts = new int[4];
            int swapCount = 0;
            calcLU(u, l, p, n, ref swapCount);
            bLU = prodMatrixVector(p, bLU, n);
            var yLU = calcY(l, bLU, n);
            var xLU = calcX(u, yLU, n);
            var x = b.Clone() as double[];
            printStr("Метод простых итераций");
            printStr();
            simpleIterationMethod(a, b, n, x, xLU, out iterationCounts[0]);
            printStr("Полученный методом простой итерации ответ:");
            printStr();
            printVector(x, n);
            printStr("Полученный LU-разложением ответ:");
            printStr();
            printVector(xLU, n);
            printStr();
            x = b.Clone() as double[];
            printStr("Метод наискорейшего спуска");
            printStr();
            fastestWayDownMethod(a, b, n, ref x, xLU, out iterationCounts[1]);
            printStr("Полученный методом наискорейшего спуска ответ:");
            printStr();
            printVector(x, n);
            printStr("Полученный LU-разложением ответ:");
            printStr();
            printVector(xLU, n);
            printStr();
            x = b.Clone() as double[];
            printStr("Метод ПВР");
            printStr();
            PVRMethod(a, b, n, ref x, xLU, out iterationCounts[2]);
            printStr("Полученный методом ПВР ответ:");
            printStr();
            printVector(x, n);
            printStr("Полученный LU-разложением ответ:");
            printStr();
            printVector(xLU, n);
            printStr();
            x = b.Clone() as double[];
            printStr("Метод сопряженных градиентов");
            printStr();
            conjurateGradientsMethod(a, b, n, ref x, xLU, out iterationCounts[3]);
            printStr("Полученный методом сопряженных градиентов ответ:");
            printStr();
            printVector(x, n);
            printStr("Полученный LU-разложением ответ:");
            printStr();
            printVector(xLU, n);
            printStr();
            printStr();
            var conditionNumber = calcNorm(a, n) * calcNorm(calcReverseA(l, u, n), n);
            printStr("Число обусловленности матрицы: " + conditionNumber.ToString("0.000"));
            printStr();
            printStr("Теоретическое число итераций для метода простых итераций: " + (Math.Log(1 / eps) / 2 * conditionNumber).ToString("0.000") + "  против  " + iterationCounts[0] + " полученных");
            printStr();
            printStr("Теоретическое число итераций для метода наискорейшего спуска: " + (Math.Log(1 / eps) * conditionNumber).ToString("0.000") + "  против  " + iterationCounts[1] + " полученных");
            printStr();
            printStr("Теоретическое число итераций для метода ПВР: " + (Math.Log(1 / eps) / 4 * Math.Sqrt(conditionNumber)).ToString("0.000") + "  против  " + iterationCounts[2] + " полученных");
            printStr();
            printStr("Теоретическое число итераций для метода сопряженных градиентов: " + (Math.Log(2 / eps) / 2 * Math.Sqrt(conditionNumber)).ToString("0.000") + "  против  " + iterationCounts[3] + " полученных");
            printStr();
            w.Close();
        }

        //метод простой итерации.
        static void simpleIterationMethod(double[,] a, double[] b, int n, double[] x, double[] answer, out int iteration)
        {
            var gamma = 0.9;
            var tau = gamma * 2 / calcMatrixNorm(a, n);
            iterationsCount = 0;
            discrepancyNorm = 0;
            printStr(" №   ");
            for (int i = 0; i < n; i++)
            {
                printStr("x" + i + "        ");
            }
            printStr("tau       Невязка   q         Погрешность");
            printStr();
            do
            {
                iterationsCount++;
                var prevX = x.Clone() as double[];
                for (int i = 0; i < n; i++)
                {
                    double sum = 0;
                    for (int j = 0; j < n; j++)
                    {
                        sum += a[i, j] * prevX[j];
                    }
                    x[i] = x[i] + tau * (b[i] - sum);
                }
                //считаем следующее значение вектора x для вычисления q.
                var nextX = new double[] { 0, 0, 0, 0 };
                for (int i = 0; i < n; i++)
                {
                    double sum = 0;
                    for (int j = 0; j < n; j++)
                    {
                        sum += a[i, j] * x[j];
                    }
                    nextX[i] = x[i] + tau * (b[i] - sum);
                }
                //выводим промежуточные данные.
                discrepancyNorm = calcVectorNorm(diffVector(b, prodMatrixVector(a, x, n), n), n);
                var q = calcQ(prevX, x, nextX, n);
                error = calcVectorNorm(diffVector(x, answer, n), n);
                printStr(iterationsCount > 9 ? iterationsCount.ToString() : " " + iterationsCount);
                printStr("   ");
                for (int i = 0; i < n; i++)
                {
                    printStr(x[i].ToString("0.0000000") + (x[i] < 0 ? "  " : "   "));
                }
                printStr(tau.ToString("0.00000") + "   " + discrepancyNorm.ToString("0.00000") + (discrepancyNorm > 9 ? "  " : "   ") + q.ToString("0.00000") + "   " + error.ToString("0.00000"));
                printStr();
            } while (discrepancyNorm > eps);
            iteration = iterationsCount;
        }

        static void fastestWayDownMethod(double[,] a, double[] b, int n, ref double[] x, double[] answer, out int iteration)
        {
            error = 0;
            iterationsCount = 0;
            discrepancyNorm = 0;
            double q = 0;
            var r = new double[n]; //вектор градиента.
            double tau = 0;
            var oldX = new double[n];
            var nextX = new double[n];
            printStr(" №   ");
            for (int i = 0; i < n; i++)
            {
                printStr("x" + i + "        ");
            }
            printStr("Невязка   tau       q       Погрешность");
            printStr();
            do
            {
                iterationsCount++;
                oldX = x;
                r = diffVector(prodMatrixVector(a, x, n), b, n);
                tau = scolarMult(r, r, n) / scolarMult(prodMatrixVector(a, r, n), r, n);
                x = diffVector(x, prodVectorNumber(r, tau, n), n);
                r = diffVector(prodMatrixVector(a, x, n), b, n);
                tau = scolarMult(r, r, n) / scolarMult(prodMatrixVector(a, r, n), r, n);
                nextX = diffVector(x, prodVectorNumber(r, tau, n), n);
                discrepancyNorm = calcVectorNorm(diffVector(b, prodMatrixVector(a, x, n), n), n);
                q = calcQ(oldX, x, nextX, n);
                error = calcVectorNorm(diffVector(x, answer, n), n);
                printStr(iterationsCount > 9 ? iterationsCount.ToString() : " " + iterationsCount);
                printStr("   ");
                for (int i = 0; i < n; i++)
                {
                    printStr(x[i].ToString("0.0000000") + (x[i] > 0 ? "   " : "  "));
                }
                printStr(discrepancyNorm.ToString("0.00000") + (discrepancyNorm > 9 ? "  " : "   ") + tau.ToString("0.00000") + "   " + q.ToString("0.00000") + "   " + error.ToString("0.00000"));
                printStr();
            } while (discrepancyNorm > eps);
            iteration = iterationsCount;
        }

        //метод ПВР.
        static void PVRMethod(double[,] a, double[] b, int n, ref double[] x, double[] answer, out int iteration)
        {
            var localEps = 0.01;
            error = 0;
            iterationsCount = 0;
            var tildaX = new double[n];
            var bestOmega = 0.0;
            var bestOmegasString = string.Empty;
            double q = 0;
            var oldX = new double[n];
            var nextX = new double[n];
            var minIter = int.MaxValue;
            for (var omega = 0.1; omega < 2; omega += 0.1)
            {
                error = 0;
                tildaX = new double[n];
                x = b.Clone() as double[];
                iterationsCount = 0;
                do
                {
                    iterationsCount++;
                    for (int i = 0; i < n; i++)
                    {
                        double sum1 = 0, sum2 = 0;
                        for (int j = 0; j < i; j++)
                        {
                            sum1 += a[i, j] * x[j];
                        }
                        for (int j = i + 1; j < n; j++)
                        {
                            sum2 += a[i, j] * x[j];
                        }
                        tildaX[i] = (b[i] - sum1 - sum2) / a[i, i];
                        x[i] = x[i] + omega * (tildaX[i] - x[i]);
                    }
                    discrepancyNorm = calcVectorNorm(diffVector(b, prodMatrixVector(a, x, n), n), n); ;
                    error = calcVectorNorm(diffVector(x, answer, n), n);
                }
                while (discrepancyNorm > localEps);
                if (minIter == iterationsCount)
                {
                    bestOmegasString += " | " + omega.ToString("0.0");
                }
                if (minIter > iterationsCount)
                {
                    minIter = iterationsCount;
                    bestOmega = omega;
                    bestOmegasString = omega.ToString("0.0") + " ";
                }
                printStr("Omega = " + omega.ToString("0.0") + "; Kоличество итераций: " + iterationsCount);
                printStr();
            }
            printStr();
            printStr("Наименьшее число итераций " + minIter + " достигнуто при omega = " + bestOmegasString);
            printStr();
            iterationsCount = 0;
            error = 0;
            tildaX = new double[n];
            x = b.Clone() as double[];
            printStr(" №   ");
            for (int i = 0; i < n; i++)
            {
                printStr("x" + i + "        ");
            }
            printStr("Невязка   q         Погрешность");
            printStr();
            do
            {
                iterationsCount++;
                oldX = x.Clone() as double[];
                for (int i = 0; i < n; i++)
                {
                    double sum1 = 0, sum2 = 0;
                    for (int j = 0; j < i; j++)
                    {
                        sum1 += a[i, j] * x[j];
                    }
                    for (int j = i + 1; j < n; j++)
                    {
                        sum2 += a[i, j] * x[j];
                    }
                    tildaX[i] = (b[i] - sum1 - sum2) / a[i, i];
                    x[i] = x[i] + bestOmega * (tildaX[i] - x[i]);
                }
                nextX = x.Clone() as double[];
                for (int i = 0; i < n; i++)
                {
                    double sum1 = 0, sum2 = 0;
                    for (int j = 0; j < i; j++)
                    {
                        sum1 += a[i, j] * nextX[j];
                    }
                    for (int j = i + 1; j < n; j++)
                    {
                        sum2 += a[i, j] * nextX[j];
                    }
                    tildaX[i] = (b[i] - sum1 - sum2) / a[i, i];
                    nextX[i] = nextX[i] + bestOmega * (tildaX[i] - x[i]);
                }
                discrepancyNorm = calcVectorNorm(diffVector(b, prodMatrixVector(a, x, n), n), n);
                q = calcQ(oldX, x, nextX, n);
                error = calcVectorNorm(diffVector(x, answer, n), n);
                printStr(iterationsCount > 9 ? iterationsCount.ToString() : " " + iterationsCount);
                printStr("   ");
                for (int i = 0; i < n; i++)
                {
                    printStr(x[i].ToString("0.0000000") + "   ");
                }
                printStr(discrepancyNorm.ToString("0.0000000") + (discrepancyNorm > 9 ? "  " : "   ") + q.ToString("0.00000") + "   " + error.ToString("0.00000"));
                printStr();
            }
            while (discrepancyNorm > eps);
            iteration = iterationsCount;
        }

        //метод сопряженных градиентов.
        static void conjurateGradientsMethod(double[,] a, double[] b, int n, ref double[] x, double[] answer, out int iteration)
        {
            error = 0;
            iterationsCount = 0;
            discrepancyNorm = 0;
            var d = new double[n]; //вектор направления.
            var g = new double[n]; //вектор градиента.
            double s = 0; //скалярный шаг.
            var prevG = negationVection(b, n);
            printStr(" №   ");
            for (int i = 0; i < n; i++)
            {
                printStr("x" + i + "        ");
            }
            printStr("Невязка   q         Погрешность");
            printStr();
            do
            {
                iterationsCount++;
                g = diffVector(b, prodMatrixVector(a, x, n), n);
                d = sumVectors(negationVection(g, n), prodVectorNumber(d, scolarMult(g, g, n) / scolarMult(prevG, prevG, n), n), n);
                s = scolarMult(d, g, n) / scolarMult(prodMatrixVector(a, d, n), d, n);
                var prevX = x.Clone() as double[];
                x = sumVectors(x, prodVectorNumber(d, s, n), n);
                prevG = g.Clone() as double[];
                var nextG = diffVector(b, prodMatrixVector(a, x, n), n);
                var nextD = sumVectors(negationVection(nextG, n), prodVectorNumber(d, scolarMult(nextG, nextG, n) / scolarMult(prevG, prevG, n), n), n);
                var nextS = scolarMult(nextD, nextG, n) / scolarMult(prodMatrixVector(a, nextD, n), nextD, n);
                var nextX = sumVectors(x, prodVectorNumber(nextD, nextS, n), n);
                discrepancyNorm = calcVectorNorm(diffVector(b, prodMatrixVector(a, x, n), n), n);
                var q = calcQ(prevX, x, nextX, n);
                error = calcVectorNorm(diffVector(x, answer, n), n);
                printStr(iterationsCount > 9 ? iterationsCount.ToString() : " " + iterationsCount);
                printStr("   ");
                for (int i = 0; i < n; i++)
                {
                    printStr(x[i].ToString("0.0000000") + "   ");
                }
                printStr(discrepancyNorm.ToString("0.0000000") + (discrepancyNorm > 9 ? "  " : "   ") + q.ToString("0.0000000") + "   " + error.ToString("0.0000000"));
                printStr();
            } while (discrepancyNorm > eps);
            iteration = iterationsCount;
        }

        //отрицание вектора.
        static double[] negationVection(double[] v, int n)
        {
            var answer = new double[n];
            for (int i = 0; i < n; i++)
            {
                answer[i] = -v[i];
            }
            return answer;
        }

        //перемножение векторов.
        static double scolarMult(double[] v1, double[] v2, int n)
        {
            double sum = 0;
            for (int i = 0; i < n; i++)
            {
                sum += v1[i] * v2[i];
            }
            return sum;
        }

        //чтение матрицы.
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

        //печать строки.
        static void printStr(string outStr)
        {
            Console.Write(outStr);
            w.Write(outStr);
        }
        static void printStr()
        {
            Console.WriteLine(string.Empty);
            w.WriteLine(string.Empty);
        }

        //печать матрицы.
        static void printMatrix(double[,] matrix, int n)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write(matrix[i, j].ToString("0.000") + "\t");
                    w.Write(matrix[i, j].ToString("0.000") + "\t");
                }
                w.WriteLine();
                Console.WriteLine();
            }
        }

        //печать вектора.
        static void printVector(double[] vector, int n)
        {
            for (int i = 0; i < n; i++)
            {
                Console.Write(vector[i].ToString("0.00000") + "   ");
                w.Write(vector[i].ToString("0.00000") + "   ");
            }
            w.WriteLine();
            Console.WriteLine();
        }

        //считает норму матрицы.
        static double calcMatrixNorm(double[,] matrix, int n)
        {
            double currMax = 0;
            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int j = 0; j < n; j++)
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

        //считает норму вектора.
        static double calcVectorNorm(double[] vector, int n)
        {
            double sum = 0;
            for (int i = 0; i < n; i++)
            {
                sum += vector[i] * vector[i];
            }
            return Math.Sqrt(sum);
        }

        //считает произведение матрицы и вектора.
        static double[] prodMatrixVector(double[,] a, double[] b, int n)
        {
            double[] prodM = new double[n];
            for (int j = 0; j < n; j++)
            {
                double sum = 0;
                for (int k = 0; k < n; k++)
                {
                    sum += a[j, k] * b[k];
                }
                prodM[j] = sum;
            }
            return prodM;
        }

        //считает произведение вектора на число.
        static double[] prodVectorNumber(double[] v, double numb, int n)
        {
            var answer = new double[n];
            for (int i = 0; i < n; i++)
            {
                answer[i] = v[i] * numb;
            }
            return answer;
        }

        //реализует перестановку строк в матрице по заданным индексам.
        static void swap(double[,] a, int n, int first, int second)
        {
            for (int i = 0; i < n; i++)
            {
                double temp = a[first, i];
                a[first, i] = a[second, i];
                a[second, i] = temp;
            }
        }

        //раскладывает матрицу A на L * U.
        static void calcLU(double[,] u, double[,] l, double[,] p, int n, ref int swapCount)
        {
            for (int i = 0; i < n; i++)
                p[i, i] = 1;
            for (int i = 0; i < n; i++)
            {
                int iMax = i;
                for (int j = i + 1; j < n; ++j)
                {
                    if (Math.Abs(u[j, i]) > Math.Abs(u[iMax, i]))
                        iMax = j;
                }
                if (Math.Abs(u[iMax, i]) < 10e-6)
                    continue;
                if (i != iMax)
                {
                    swap(u, n, i, iMax);
                    swap(l, n, i, iMax);
                    swap(p, n, i, iMax);
                    swapCount++;
                }
                for (int j = 0; j < n; j++)
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
            for (int i = 0; i < n; i++)
            {
                double t = b[i];
                for (int j = 0; j < n; j++)
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

        //считает разность двух векторов.
        static double[] diffVector(double[] a, double[] b, int n)
        {
            var difference = new double[n];
            for (int i = 0; i < n; i++)
            {
                difference[i] = a[i] - b[i];
            }
            return difference;
        }

        //считает сумму двух векторов.
        static double[] sumVectors(double[] a, double[] b, int n)
        {
            var sum = new double[n];
            for (int i = 0; i < n; i++)
            {
                sum[i] = a[i] + b[i];
            }
            return sum;
        }

        static double calcQ(double[] prevX, double[] x, double[] nextX, int n)
        {
            return calcVectorNorm(diffVector(nextX, x, n), n) / calcVectorNorm(diffVector(x, prevX, n), n);
        }

        //обратная матрица.
        static double[,] calcReverseA(double[,] l, double[,] u, int n)
        {
            List<double[]> y = new List<double[]>();
            //вычисление вектора y.
            for (int i = 0; i < n; i++)
            {
                double[] t = new double[n];
                t[i] = 1;
                y.Add(calcY(l, t, n));
            }
            double[,] x = new double[n, n];
            //вычисление вектора x.
            for (int i = 0; i < n; i++)
            {
                double[] t = calcX(u, y[i], n);
                for (int j = 0; j < n; j++)
                    x[j, i] = t[j];
            }
            return x;
        }

        //норма матрицы.
        static double calcNorm(double[,] matrix, int n)
        {
            double currMax = 0;
            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int j = 0; j < n; j++)
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
    }
}