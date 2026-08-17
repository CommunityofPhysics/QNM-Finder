// File: Vector.cs

using System.Numerics;
using MathNet.Numerics;
using System.Text;
using MathNet.Numerics.LinearAlgebra;

namespace CLACFramework;

public sealed class VectorC
{
    private readonly Matrix<Complex> _v;   // n×1 column matrix

    public int Rows => _v.RowCount;
    public int Cols => 1;
    public int Size => _v.RowCount;

    internal Matrix<Complex> Inner => _v;

    // Constructors ----------------------------------------------------

    public VectorC(Complex[] data)
        => _v = Matrix<Complex>.Build.Dense(data.Length, 1, (i, j) => data[i]);

    public VectorC((double real, double imag)[] data)
        => _v = Matrix<Complex>.Build.Dense(data.Length, 1,
               (i, j) => new Complex(data[i].real, data[i].imag));

    public VectorC(int s)
        => _v = Matrix<Complex>.Build.Dense(s, 1, Complex.Zero);

    internal VectorC(Matrix<Complex> m)
    {
        if (m.ColumnCount != 1)
            throw new ArgumentException("VectorC must wrap a column matrix.");
        _v = m;
    }

    // Indexer --------------------------------------------------------
    public Complex this[int i]
    {
        get => _v[i, 0];
        set => _v[i, 0] = value;
    }

    // Array ----------------------------------------------------------
    public Complex[] ToArray() => _v.Column(0).ToArray();

    // Array to VectorC ------------------------------------------------
    public static VectorC FromArray(double[] arr)
        => new VectorC(Matrix<Complex>.Build.Dense(arr.Length, 1, (i, j) => new Complex(arr[i], 0.0)));

    // To MatrixC ------------------------------------------------------
    public MatrixC ToMatrixC()
    {
        return new MatrixC(_v);
    }

    // Clone / Copy ---------------------------------------------------
    public VectorC Clone()
        => new VectorC(_v.Clone());

    public void CopyFrom(VectorC other)
    {
        if (other.Size != Size) throw new ArgumentException("Vector sizes must match.");

        _v.SetColumn(0, other._v.Column(0));
    }

    public void CopyFrom(VectorC other, int start, int length)
    {
        if (start < 0 || length < 0 || start + length > other.Size || length > Size)
            throw new ArgumentOutOfRangeException();

        for (int i = 0; i < length; i++)
            this[i] = other[start + i];
    }

    public VectorC GetSubvector(int start, int length)
    {
        if (start < 0 || length < 0 || start + length > Size)
            throw new ArgumentOutOfRangeException();

        VectorC v = new VectorC(length);
        v.CopyFrom(this, start, length);
        return v;
    }

    public void SetSubvector(int start, VectorC other)
    {
        int length = other.Size;

        if (start < 0 || length < 0 || start + length > Size)
            throw new ArgumentOutOfRangeException();

        for (int i = 0; i < length; i++)
            this[start + i] = other[i];
    }

    // Clear --------------------------------------------------------
    public void Clear()
    {
        _v.Clear();
    }

    // Zero -----------------------------------------------------------
    public static VectorC Zero(int n) => new VectorC(n);

    public static VectorC[] Zero(int l, int n)
    {
        VectorC[] arr = new VectorC[l];
        for (int k = 0; k < l; k++)
            arr[k] = Zero(n);
        return arr;
    }

    public VectorC ZeroLike() => Zero(Size);

    // Basis vectors --------------------------------------------------
    public static VectorC Basis(int n, int j)
    {
        VectorC e = VectorC.Zero(n);
        e[j] = Complex.One;
        return e;
    }

    // NaN ------------------------------------------------------------
    public static VectorC NaN(int n)
        => new VectorC(Matrix<Complex>.Build.Dense(n, 1, Complex.NaN));

    public static VectorC[] NaN(int l, int n)
    {
        VectorC[] arr = new VectorC[l];
        for (int k = 0; k < l; k++)
            arr[k] = NaN(n);
        return arr;
    }

    // Arithmetic -----------------------------------------------------
    public static VectorC operator +(VectorC a, VectorC b)
        => new VectorC(a._v + b._v);

    public static VectorC operator -(VectorC a)
        => new VectorC(-a._v);

    public static VectorC operator -(VectorC a, VectorC b)
        => new VectorC(a._v - b._v);

    public static VectorC operator *(Complex s, VectorC b)
        => new VectorC(s * b._v);

    public static VectorC operator *(VectorC a, Complex s)
        => new VectorC(a._v * s);

    public static VectorC operator /(VectorC a, Complex s)
        => new VectorC(a._v / s);


    // Transpose / Hermitian ------------------------------------------
    public MatrixC Transpose()
        => new MatrixC(_v.Transpose());

    public MatrixC Hermitian()
        => new MatrixC(_v.ConjugateTranspose());

    // Inner Product --------------------------------------------------
    public Complex Dot(VectorC other)
    {
        if (other.Rows != Rows)
            throw new ArgumentException("Vector sizes must match.");

        Complex sum = Complex.Zero;
        for (int i = 0; i < Rows; i++)
            sum += Complex.Conjugate(this[i]) * other[i];

        return sum;
    }

    public Complex EuclideanDot(VectorC other)
    {
        if (other.Rows != Rows)
            throw new ArgumentException("Vector sizes must match.");
        Complex sum = Complex.Zero;
        for (int i = 0; i < Rows; i++)
            sum += this[i] * other[i];
        return sum;
    }

    // Hadamard Product -----------------------------------------------
    public VectorC HadamardProduct(VectorC other)
    {
        if (other.Rows != Rows)
            throw new ArgumentException("Vector sizes must match.");

        return new VectorC(_v.PointwiseMultiply(other._v));
    }

    // Norm -----------------------------------------------------------
    public double Norm()
    {
        Complex ip = Dot(this);
        double val = Math.Max(0.0, ip.Real);
        return Math.Sqrt(val);
    }

    public void Normalize()
    {
        double norm = Norm();
        if (norm > 0.0)
            _v.Multiply(1.0 / norm, _v);
    }

    public void Normalize(VectorC other)
    {
        Complex scale = other.Dot(this);
        if (Complex.Abs(scale) > 0.0)
            _v.Multiply(1.0 / scale, _v);
    }

    // In-place operations ------------------------------------------------
    public void AddIP(VectorC other)
    {
        if (other.Rows != Rows)
            throw new ArgumentException("Vector sizes must match.");

        _v.Add(other._v, _v);
    }

    public void NegateIP()
    {
        _v.Multiply(-1.0, _v);
    }

    public void SubIP(VectorC other)
    {
        if (other.Rows != Rows)
            throw new ArgumentException("Vector sizes must match.");

        _v.Subtract(other._v, _v);
    }

    public void ScaleIP(Complex s)
        => _v.Multiply(s, _v);

    public void HadamardProductIP(VectorC other)
    {
        if (other.Rows != Rows)
            throw new ArgumentException("Vector sizes must match.");

        _v.PointwiseMultiply(other._v, _v);
    }

    // y ← y + a*x 
    public void AxPyIP(Complex a, VectorC x)
    {
        if (x.Rows != Rows) throw new ArgumentException("Vector sizes must match.");
        if (a == Complex.Zero) return;

        for (int i = 0; i < Rows; i++)
            this[i] += a * x[i];
    }

    // y ← a*this + b*other
    public void CombineIP(Complex a, VectorC other, Complex b)
    {
        if (other.Rows != Rows) throw new ArgumentException("Vector sizes must match.");

        for (int i = 0; i < Rows; i++)
            this[i] = a * this[i] + b * other[i];
    }

    public void ConjugateIP()
    {
        _v.MapInplace(c => Complex.Conjugate(c), Zeros.AllowSkip);
    }

    // Safe Eval
    public static (VectorC, bool) SafeEval(Func<Complex, VectorC> F, Complex z)
    {
        try
        {
            VectorC value = F(z);

            int n = value.Size;

            for (int i = 0; i < n; i++)
            {
                if (Complex.IsNaN(value[i]) || Complex.IsInfinity(value[i]))
                    return (value, true);
            }

            return (value, false);
        }
        catch
        {
            try
            {
                VectorC probe = F(ScalarC.RandomComplex(z, 1.0));

                int n = probe.Size;

                return (VectorC.NaN(n), true);
            }
            catch
            {
                return (null!, true);
            }
        }
    }

    public static (VectorC[], bool) SafeEval(Func<Complex, VectorC[]> F, Complex z)
    {
        try
        {
            VectorC[] value = F(z);

            int l = value.Length;
            int n = value[0].Size;

            for (int k = 0; k < l; k++)
            {
                VectorC v = value[k];

                for (int i = 0; i < n; i++)
                {
                    if (Complex.IsNaN(v[i]) || Complex.IsInfinity(v[i]))
                        return (value, true);
                }
            }

            return (value, false);
        }
        catch
        {
            try
            {
                VectorC[] probe = F(ScalarC.RandomComplex(z, 1.0));

                int l = probe.Length;
                int n = probe[0].Size;

                return (VectorC.NaN(l, n), true);
            }
            catch
            {
                return (null!, true);
            }
        }
    }

    // Pretty Print ------------------------------------------------------

    public string Format(string fmt)
    {
        StringBuilder sb = new StringBuilder(Size * 16);

        for (int i = 0; i < Size; i++)
        {
            sb.Append('[');

            sb.Append(_v[i, 0].Format(fmt));

            sb.Append(']');

            if (i < Size - 1) sb.Append('\n');
        }

        return sb.ToString();
    }

    public override string ToString() => Format(null);
    public string ToString(string fmt) => Format(fmt);

}
