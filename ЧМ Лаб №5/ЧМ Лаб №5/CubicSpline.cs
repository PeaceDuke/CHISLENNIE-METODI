namespace ЧМ_Лаб__5
{
    internal static class CubicSpline
    {
        private static SplineCoeffs[] _splines;

        // Параметры сплайна
        private struct SplineCoeffs
        {
            public double A, B, C, D, X;
        }

        // Строит сплайна по таблице (x, y)
        public static void BuildSpline(PointD[] table, int n)
        {
            _splines = new SplineCoeffs[n];
            for (var i = 0; i < n; i++)
            {
                _splines[i].X = table[i].X;
                _splines[i].A = table[i].Y;
            }
            _splines[0].C = _splines[n - 1].C = 0.0;

            // Ищем параметры сплайна методом прогонки для трехдиагональных матриц. Прямой ход
            var alpha = new double[n - 1];
            var beta = new double[n - 1];
            alpha[0] = beta[0] = 0.0;
            for (var i = 1; i < n - 1; i++)
            {
                var hi = table[i].X - table[i - 1].X;
                var hi1 = table[i + 1].X - table[i].X;
                var a = hi;
                var c = 2.0 * (hi + hi1);
                var b = hi1;
                var f = 6.0 * ((table[i + 1].Y - table[i].Y) / hi1 - (table[i].Y - table[i - 1].Y) / hi);
                var z = a * alpha[i - 1] + c;
                alpha[i] = -b / z;
                beta[i] = (f - a * beta[i - 1]) / z;
            }

            // Обратный ход
            for (var i = n - 2; i > 0; i--)
            {
                _splines[i].C = alpha[i] * _splines[i + 1].C + beta[i];
            }

            // Находим b[i] и d[i], используя найденные ранее c[i]
            for (var i = n - 1; i > 0; --i)
            {
                var hi = table[i].X - table[i - 1].X;
                _splines[i].D = (_splines[i].C - _splines[i - 1].C) / hi;
                _splines[i].B = hi * (2.0 * _splines[i].C + _splines[i - 1].C) / 6.0 + (table[i].Y - table[i - 1].Y) / hi;
            }
        }

        // Считает значение сплайна в произвольной точке
        public static double InterSplines(double x)
        {
            var n = _splines.Length;
            SplineCoeffs s;

            if (x <= _splines[0].X) // x меньше левой границы, берем первый отрезок
            {
                s = _splines[0];
            }
            else if (x >= _splines[n - 1].X) // x больше правой границы, берем последний отрезок
            {
                s = _splines[n - 1];
            }
            else // x лежит между граничными точками сплайна, ищем нужный отрезок бинарным поиском
            {
                var i = 0;
                var j = n - 1;
                while (i + 1 < j)
                {
                    var k = i + (j - i) / 2;
                    if (x <= _splines[k].X)
                    {
                        j = k;
                    }
                    else
                    {
                        i = k;
                    }
                }
                s = _splines[j];
            }
            var dx = x - s.X;
            // Считаем значение сплайна в заданной точке
            return s.A + (s.B + (s.C / 2.0 + s.D * dx / 6.0) * dx) * dx;
        }
    }
}
