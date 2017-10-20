using System;

namespace ЧМ_Лаб_4
{
    class Program
    {
        static void Main(string[] args)
        {
            var variant16 = new EqationClass("cos( y + 0,5 ) + x - 0,8 = 0 // sin( x ) - 2y - 1,6 = 0", 1.0, 1.0);
            var variant17 = new EqationClass("sin( y - 1 ) + x - 1,3 = 0 // y - sin( x + 1 ) - 0,8 = 0", 1.0, 1.0);
        }

        static void SimpleIterMethod()
        {
            //чонибудь
        }

        static void NewtonMethod()
        {
            //дарова))
            //дароу
        }
    }
}
