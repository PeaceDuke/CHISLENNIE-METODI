using System;
using System.Collections.Generic;

namespace ЧМ_Лаб_4
{
    class EqationClass
    {
        private double x;
        private double y;
        private string eqaution;

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

        public string Eqation
        {
            get { return eqaution; }
        }

        public EqationClass(string _eqaution, double _x, double _y)
        {
            eqaution = _eqaution;
            x = _x;
            y = _y;
        }

        public double CalcFirstEqaution()
        {
            var eq = eqaution.Split(' ');
            int i = 0;
            return SubCalc(ref i, eq, 1);
        }

        private void ThridPrioritySearch(ref string[] eq, int startPos = 0, int endPoint )
        {
            for (int i = startPos; i < eq.Length; i++)
            {
                if (eq[i].Equals("-"))
                {
                    eq[i - 1] = (Convert.ToDouble(eq[i - 1]) - Convert.ToDouble(eq[i + 1])).ToString();
                }
                if (eq[i].Equals("+"))
                {
                    eq[i - 1] = (Convert.ToDouble(eq[i - 1]) + Convert.ToDouble(eq[i + 1])).ToString();
                }
            }
        }

        private void SecondPrioritySearch(ref string[] eq, int startPos = 0, int endPoint)
        {
            for (int i = startPos; i < eq.Length; i++)
            {
                if (eq[i].Equals("/"))
                {
                    eq[i - 1] = (Convert.ToDouble(eq[i - 1]) / Convert.ToDouble(eq[i + 1])).ToString();
                }
                if (eq[i].Equals("*"))
                {
                    eq[i - 1] = (Convert.ToDouble(eq[i - 1]) * Convert.ToDouble(eq[i + 1])).ToString();
                }
            }
        }

        private void FirstPrioritySearch(ref string[] eq, int startPos = 0, int endPoint)
        {
            for (int i = startPos; i < eq.Length; i++)
            {
                if (eq[i].Equals("("))
                {

                }
            }
        }
    }
}
