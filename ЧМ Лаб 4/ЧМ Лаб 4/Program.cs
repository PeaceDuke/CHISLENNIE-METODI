using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static System.Console;

namespace ЧМ_Лаб_4
{
    class Program
    {
        class EqationClass
        {
            private double x;
            private double y;
            private string eqaution1;
            private string eqaution2;

            public double X
            {
                get { return x; }
                set { x = value; }
            }

            public double Y
            {
                get { return y; }
                set { y = value; }
            }

            public string Eqation1
            {
                get { return eqaution1; }
            }

            public string Eqation2
            {
                get { return eqaution2; }
            }

            public EqationClass(string _eqaution, double _x, double _y)
            {
                var temp = _eqaution.Split(new string[] { " // " }, StringSplitOptions.RemoveEmptyEntries);
                eqaution1 = temp[0];
                eqaution2 = temp[1];
                x = _x;
                y = _y;
            }

            public double CalcFirstEqaution()
            {
                var eq = eqaution1.Split(' ');
                int i = 0;
                return SubCalc(ref i, eq, 1);
            }

            public double CalcSecondEqaution()
            {
                var eq = eqaution2.Split(' ');
                int i = 0;
                return SubCalc(ref i, eq, 1);
            }

            private double FuncCalc(string func, ref int startPoint, string[] eq, double multiplier)
            {
                double funcVal = SubCalc(ref startPoint, eq, 1);
                return func.Equals("sin(") ? Math.Sin(funcVal) : Math.Cos(funcVal);
            }

            private double SubCalc(ref int startPoint, string[] eq, double multiplier)
            {
                double subSumm = 0;
                bool negNext = false;
                for (int j = startPoint; j < eq.Length; j++)
                {
                    if(eq[j].Equals(")"))
                    {
                        startPoint = j;
                        return subSumm;
                    }
                    if (eq[j].Equals("sin(") || eq[j].Equals("cos("))
                    {
                        j++;
                        double funVal = FuncCalc(eq[j - 1], ref j, eq, 1);
                        subSumm += negNext ? -funVal : funVal;
                        negNext = false;
                        continue;
                    }
                    if (eq[j].Equals("x"))
                    {
                        subSumm += negNext ? -x : x;
                        negNext = false;
                        continue;
                    }
                    if (eq[j].Equals("y"))
                    {
                        subSumm += negNext ? -y : y;
                        negNext = false;
                        continue;
                    }
                    if (eq[j].Equals("2y"))
                    {
                        subSumm += 2 * (negNext ? -y : y);
                        negNext = false;
                        continue;
                    }
                    if (eq[j].Equals("2x"))
                    {
                        subSumm += 2 * (negNext ? -x : x);
                        negNext = false;
                        continue;
                    }
                    if (eq[j].Equals("-"))
                    {
                        negNext = true;
                        continue;
                    }
                    if (eq[j].Equals("+"))
                    {
                        continue;
                    }
                    if (eq[j].Equals("="))
                    {
                        break;
                    }
                    var num = Convert.ToDouble(eq[j]);
                    subSumm += negNext ? -num : num;
                    negNext = false;
                }
                return subSumm;
            }
        }

        static void Main(string[] args)
        {
            var task1 = new EqationClass("cos( y + 0,5 ) + x - 0,8 = 0 // sin( x ) - 2y - 1,6 = 0", 1.0, 1.0);
            var task2 = new EqationClass("sin( y - 1 ) + x - 1,3 = 0 // y - sin( x + 1 ) - 0,8 = 0", 1.0, 1.0);
        }

        static void SimpleIterMethod()
        {

        }
    }
}
