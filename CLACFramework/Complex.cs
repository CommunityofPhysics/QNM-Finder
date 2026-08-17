// File: Complex.cs

using MathNet.Numerics;
using System;
using System.Globalization;
using System.Numerics;

namespace CLACFramework;

public sealed class ScalarC
{
    public static int Mod(int i, int k)
    {
        return ((i % k) + k) % k;
    }

    public static Complex SafeSqrt(Complex z)
    {
        if (double.IsInfinity(z.Real) || double.IsInfinity(z.Imaginary))
        {
            double arg = Math.Atan2(z.Imaginary, z.Real) / 2.0;
            double mag = double.PositiveInfinity;
            return new Complex(mag * Math.Cos(arg), mag * Math.Sin(arg));
        }

        return Complex.Sqrt(z);
    }

    public static Complex IntPow(Complex z, int k)
    {
        if (k == 0) return Complex.One;
        Complex product = Complex.One;
        for (int i = 0; i < k; i++)
            product *= z;
        return product;
    }

    public static Complex[] Zero(int n)
    {
        Complex[] Zero = new Complex[n];
        Array.Fill(Zero, Complex.Zero);
        return Zero;
    }

    public static Func<Complex, Complex> One()
    {
        return z => Complex.One;
    }

    public static Func<Complex, Complex[]> One(int n)
    {
        return z =>
        {
            Complex[] One = new Complex[n];
            Array.Fill(One, Complex.One);
            return One;
        };
    }

    public static Complex[] NaN(int n)
    {
        Complex[] NaN = new Complex[n];
        Array.Fill(NaN, Complex.NaN);
        return NaN;
    }

    public static (Complex, bool) SafeEval(Func<Complex, Complex> F, Complex z)
    {
        try
        {
            Complex value = F(z);

            if (Complex.IsNaN(value) || Complex.IsInfinity(value))
                return (value, true);

            return (value, false);
        }
        catch
        {
            return (Complex.NaN, true);
        }
    }

    public static (Complex[], bool) SafeEval(Func<Complex, Complex[]> F, Complex z)
    {
        try
        {
            Complex[] value = F(z);

            int l = value.Length;

            for (int k = 0; k < l; k++)
            {
                if (Complex.IsNaN(value[k]) || Complex.IsInfinity(value[k]))
                    return (value, true);
            }

            return (value, false);
        }
        catch
        {
            try
            {
                Complex[] probe = F(ScalarC.RandomComplex(z, 1.0));
                int l = probe.Length;
                return (ScalarC.NaN(l), true);
            }
            catch
            {
                return (null!, true);
            }
        }
    }

    public static Complex RandomComplex(Complex center, double radius)
    {
        Random rng = Random.Shared;

        double theta = 2.0 * Math.PI * rng.NextDouble();
        double r = radius * Math.Sqrt(rng.NextDouble());

        return center + Complex.FromPolarCoordinates(r, theta);
    }

    public static Complex PUDifference(Complex zf, Complex zi)
    {
        Complex dz = zf - zi;
        double dzIm = dz.Imaginary;

        if (dzIm > Math.PI)
            dz -= new Complex(0.0, 2.0 * Math.PI);
        else if (dzIm < -Math.PI)
            dz += new Complex(0.0, 2.0 * Math.PI);

        return dz;
    }

    public static Complex Parse(string spec)
    {
        ArgumentNullException.ThrowIfNull(spec);
        spec = spec.Trim();

        if (!spec.StartsWith("(") || !spec.EndsWith(")"))
            throw new Exception($"Invalid ScalarC format: '{spec}'");

        string inner = spec.Substring(1, spec.Length - 2).Trim();

        string[] parts = inner.Split(',');

        if (parts.Length != 2)
            throw new Exception($"ScalarC requires two components: '{spec}'");

        string reText = string.Concat(parts[0].Where(c => !char.IsWhiteSpace(c)));
        string imText = string.Concat(parts[1].Where(c => !char.IsWhiteSpace(c)));

        if (!double.TryParse(reText, out double re))
            throw new Exception($"Invalid real part in ScalarC: '{parts[0]}'");

        if (!double.TryParse(imText, out double im))
            throw new Exception($"Invalid imaginary part in ScalarC: '{parts[1]}'");

        return new Complex(re, im);
    }
}

public sealed class DomainC
{
    private readonly Complex _center;
    private readonly double _width;
    private readonly double _height;

    // Constructor ----------------------------------------------------
    public DomainC(Complex center, double width, double height)
    {
        if (width <= 0.0) throw new ArgumentOutOfRangeException(nameof(width));
        if (height <= 0.0) throw new ArgumentOutOfRangeException(nameof(height));

        _center = center;
        _width = width;
        _height = height;
    }

    // Properties
    public Complex Center => _center;
    public double Width => _width;
    public double Height => _height;

    public double HalfWidth => _width * 0.5;
    public double HalfHeight => _height * 0.5;

    public double MinRe => _center.Real - HalfWidth;
    public double MaxRe => _center.Real + HalfWidth;
    public double MinIm => _center.Imaginary - HalfHeight;
    public double MaxIm => _center.Imaginary + HalfHeight;

    public double[] Boundary()
        => new[] { MinIm, MaxRe, MaxIm, MinRe };

    public double Perimeter => 2.0 * (_width + _height);
    public double Area => _width * _height;
    public bool IsSquare => Math.Abs(Width - Height) < 1e-14;


    // Corners --------------------------------------------------------
    public Complex BottomLeft => new(MinRe, MinIm);
    public Complex BottomRight => new(MaxRe, MinIm);
    public Complex TopLeft => new(MinRe, MaxIm);
    public Complex TopRight => new(MaxRe, MaxIm);

    public Complex[] Corners()
        => new[] { BottomLeft, BottomRight, TopRight, TopLeft };

    // Containment ----------------------------------------------------
    public bool DoesContain(Complex z)
        => (z.Real >= MinRe && z.Real <= MaxRe &&
            z.Imaginary >= MinIm && z.Imaginary <= MaxIm);

    // Transformations
    public DomainC MoveTo(Complex c)
        => new DomainC(c, Width, Height);
    public DomainC ResizeTo(double w, double h)
        => new DomainC(Center, w, h);

    // Parameterization: Starts at Bottom-Left
    public (Complex Point, Complex Tangent) Parametrize(double t)
    {
        double s = (t - Math.Floor(t)) * Perimeter;

        double s0 = Width;
        double s1 = Width + Height;
        double s2 = Width + Height + Width;

        if (s <= s0)
        {
            double x = MinRe + s;
            return (new Complex(x, MinIm), new Complex(1.0, 0.0));
        }
        else if (s <= s1)
        {
            double y = MinIm + (s - s0);
            return (new Complex(MaxRe, y), new Complex(0.0, 1.0));
        }
        else if (s <= s2)
        {
            double x = MaxRe - (s - s1);
            return (new Complex(x, MaxIm), new Complex(-1.0, 0.0));
        }
        else
        {
            double y = MaxIm - (s - s2);
            return (new Complex(MinRe, y), new Complex(0.0, -1.0));
        }

    }

    public double InverseParametrize(Complex z)
    {
        double x = z.Real;
        double y = z.Imaginary;

        double eps = 1e-12;
        double P = Perimeter;
        double s;

        if (Math.Abs(y - MinIm) < eps)
            s = x - MinRe;

        else if (Math.Abs(x - MaxRe) < eps)
            s = Width + (y - MinIm);

        else if (Math.Abs(y - MaxIm) < eps)
            s = Width + Height + (MaxRe - x);

        else if (Math.Abs(x - MinRe) < eps)
            s = Width + Height + Width + (MaxIm - y);

        else
            throw new ArgumentException("Point is not on the boundary.");

        return s / P;
    }

    // Subdivision ----------------------------------------------------
    public DomainC[] Subdivide(int mx, int my)
    {
        if (mx <= 0 || my <= 0)
            throw new ArgumentOutOfRangeException("Subdivision factors must be positive.");

        double childW = _width / mx;
        double childH = _height / my;

        double x0 = MinRe;
        double y0 = MinIm;

        DomainC[] children = new DomainC[mx * my];
        int k = 0;

        for (int i = 0; i < mx; i++)
        {
            for (int j = 0; j < my; j++)
            {
                double cx = x0 + i * childW + childW * 0.5;
                double cy = y0 + j * childH + childH * 0.5;

                children[k++] = new DomainC(new Complex(cx, cy), childW, childH);
            }
        }

        return children;
    }

    // Pretty print ----------------------------------------------------
    public string Format(string fmt)
    {
        string c = _center.Format(fmt);
        string w = _width.Format(fmt);
        string h = _height.Format(fmt);

        return $"[C: {c}, W: {w}, H: {h}]";
    }

    public override string ToString() => Format(null!);
    public string ToString(string fmt) => Format(fmt);

}

public sealed class SqDomainC
{
    private readonly DomainC _rect;

    public SqDomainC(Complex center, double edge)
        => _rect = new DomainC(center, edge, edge);

    public Complex Center => _rect.Center;
    public double Edge => _rect.Width;

    public double HalfEdge => _rect.HalfWidth;
    public double MinRe => _rect.MinRe;
    public double MaxRe => _rect.MaxRe;
    public double MinIm => _rect.MinIm;
    public double MaxIm => _rect.MaxIm;

    public double[] Boundary()
        => _rect.Boundary();

    public double Perimeter => _rect.Perimeter;
    public double Area => _rect.Area;

    public Complex BottomLeft => _rect.BottomLeft;
    public Complex BottomRight => _rect.BottomRight;
    public Complex TopLeft => _rect.TopLeft;
    public Complex TopRight => _rect.TopRight;
    public Complex[] Corners()
        => _rect.Corners();

    public bool DoesContain(Complex z)
        => _rect.DoesContain(z);

    public SqDomainC MoveTo(Complex center)
        => new SqDomainC(center, Edge);
    public SqDomainC ResizeTo(double edge)
        => new SqDomainC(Center, edge);

    public (Complex Point, Complex Tangent) Parametrize(double t)
        => _rect.Parametrize(t);

    public double InverseParametrize(Complex z)
        => _rect.InverseParametrize(z);

    public SqDomainC[] Subdivide(int m)
    {
        DomainC[] rects = _rect.Subdivide(m, m);
        SqDomainC[] squares = new SqDomainC[rects.Length];

        for (int i = 0; i < rects.Length; i++)
            squares[i] = new SqDomainC(rects[i].Center, rects[i].Width);

        return squares;
    }

    // Pretty print ----------------------------------------------------
    public string Format(string fmt)
    {
        string c = Center.Format(fmt);
        string e = Edge.Format(fmt);

        return $"[C: {c}, E: {e}]";
    }

    public override string ToString() => Format(null!);
    public string ToString(string fmt) => Format(fmt);

}

public static class NumberFormat
{
    // -------------------- double --------------------
    public static string Format(this double x)
        => x.ToString(CultureInfo.InvariantCulture);

    public static string Format(this double x, string fmt)
        => string.IsNullOrWhiteSpace(fmt)
            ? x.ToString(CultureInfo.InvariantCulture)
            : x.ToString(fmt, CultureInfo.InvariantCulture);

    // -------------------- Complex --------------------
    public static string Format(this Complex z)
        => $"({z.Real.ToString(CultureInfo.InvariantCulture)}, " +
           $"{z.Imaginary.ToString(CultureInfo.InvariantCulture)})";

    public static string Format(this Complex z, string fmt)
        => string.IsNullOrWhiteSpace(fmt)
            ? z.Format()
            : $"({z.Real.ToString(fmt, CultureInfo.InvariantCulture)}, " +
               $"{z.Imaginary.ToString(fmt, CultureInfo.InvariantCulture)})";

    // -------------------- Rounding --------------------
    public static string Round(this double x, double tol)
    {
        double v = Math.Abs(x) < tol ? 0.0 : x;
        return v.ToString(CultureInfo.InvariantCulture);
    }

    public static string Round(this Complex z, double tol)
    {
        double re = Math.Abs(z.Real) < tol ? 0.0 : z.Real;
        double im = Math.Abs(z.Imaginary) < tol ? 0.0 : z.Imaginary;

        return $"({re.ToString(CultureInfo.InvariantCulture)}, " +
               $"{im.ToString(CultureInfo.InvariantCulture)})";
    }

}

