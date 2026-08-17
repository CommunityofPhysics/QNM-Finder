// File: Matrix.cs

using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using System.Numerics;
using System.Text;

namespace CLACFramework;

public sealed class MatrixC
{
    private readonly Matrix<Complex> _m;

    public int Rows => _m.RowCount;
    public int Cols => _m.ColumnCount;
    public (int Rows, int Cols) Size => (Rows, Cols);
    public bool IsSquare => Rows == Cols;

    internal Matrix<Complex> Inner => _m;

    // Constructors ----------------------------------------------------
    public MatrixC(Complex[,] data)
        => _m = Matrix<Complex>.Build.DenseOfArray(data);

    public MatrixC((double real, double imag)[,] data)
    {
        int r = data.GetLength(0), c = data.GetLength(1);
        _m = Matrix<Complex>.Build.Dense(r, c,
            (i, j) => new Complex(data[i, j].real, data[i, j].imag));
    }

    public MatrixC(int rows, int cols)
    => _m = Matrix<Complex>.Build.Dense(rows, cols, Complex.Zero);

    internal MatrixC(Matrix<Complex> m)
        => _m = m;

    // Indexer --------------------------------------------------------
    public Complex this[int i, int j]
    {
        get => _m[i, j];
        set => _m[i, j] = value;
    }

    // Array ----------------------------------------------------------
    public Complex[,] ToArray() => _m.ToArray();

    public static MatrixC FromArray(double[,] arr)
    {
        int r = arr.GetLength(0), c = arr.GetLength(1);
        return new MatrixC(Matrix<Complex>.Build.Dense(r, c,
            (i, j) => new Complex(arr[i, j], 0.0)));
    }

    // Column get/set -------------------------------------------------
    public VectorC GetColumn(int col)
        => new VectorC(_m.Column(col).ToArray());

    public void SetColumn(int col, VectorC v)
    {
        if (v.Size != Rows)
            throw new ArgumentException("Column length mismatch.");
        _m.SetColumn(col, v.ToArray());
    }

    // Clone / Copy ---------------------------------------------------
    public MatrixC Clone() => new MatrixC(_m.Clone());

    public void CopyFrom(MatrixC other)
    {
        if (Rows != other.Rows || Cols != other.Cols)
            throw new ArgumentException("Matrix sizes must match.");

        _m.SetSubMatrix(0, 0, other._m);
    }

    // Clear --------------------------------------------------------
    public void Clear()
    {
        _m.Clear();
    }

    public void ZeroRow(int row)
    {
        int cols = this.Cols;
        for (int j = 0; j < cols; j++)
            this[row, j] = Complex.Zero;
    }

    public void ZeroCol(int col)
    {
        int rows = this.Rows;
        for (int i = 0; i < rows; i++)
            this[i, col] = Complex.Zero;
    }

    // Zero / Identity ------------------------------------------------
    public static MatrixC Zero(int r, int c)
        => new MatrixC(r, c);

    public static MatrixC[] Zero(int l, int r, int c)
    {
        MatrixC[] arr = new MatrixC[l];
        for (int k = 0; k < l; k++)
            arr[k] = MatrixC.Zero(r, c);
        return arr;
    }

    public MatrixC ZeroLike()
        => MatrixC.Zero(Rows, Cols);

    public static MatrixC Identity(int n)
        => new MatrixC(Matrix<Complex>.Build.DenseIdentity(n));

    // NaNs ---------------------------------------------------------
    public static MatrixC NaN(int r, int c)
        => new MatrixC(Matrix<Complex>.Build.Dense(r, c, Complex.NaN));

    public static MatrixC[] NaN(int l, int r, int c)
    {
        MatrixC[] arr = new MatrixC[l];
        for (int k = 0; k < l; k++)
            arr[k] = MatrixC.NaN(r, c);
        return arr;
    }

    // Arithmetic -----------------------------------------------------
    public static MatrixC operator +(MatrixC a, MatrixC b)
        => new MatrixC(a._m + b._m);

    public static MatrixC operator -(MatrixC a)
        => new MatrixC(-a._m);

    public static MatrixC operator -(MatrixC a, MatrixC b)
        => new MatrixC(a._m - b._m);

    public static MatrixC operator *(Complex s, MatrixC a)
        => new MatrixC(a._m * s);

    public static MatrixC operator *(MatrixC a, Complex s)
        => new MatrixC(a._m * s);

    public static MatrixC operator /(MatrixC a, Complex s)
    => new MatrixC(a._m / s);

    public static VectorC operator *(MatrixC a, VectorC v)
    {
        if (a.Cols != v.Rows)
            throw new ArgumentException("Matrix columns must match vector rows.");

        return new VectorC(a._m * v.Inner);
    }

    public static MatrixC operator *(MatrixC a, MatrixC b)
    {
        if (a.Cols != b.Rows)
            throw new ArgumentException("Matrix a columns must match matrix b rows.");

        return new MatrixC(a._m * b._m);
    }

    // Transpose / Hermitian ------------------------------------------
    public MatrixC Transpose()
        => new MatrixC(_m.Transpose());

    public MatrixC Hermitian()
        => new MatrixC(_m.ConjugateTranspose());

    // Norms / Trace / Det --------------------------------------------
    public double FrobeniusNorm() => _m.FrobeniusNorm();

    public Complex Trace()
    {
        Complex s = Complex.Zero;
        int n = Math.Min(Rows, Cols);
        for (int i = 0; i < n; i++)
            s += _m[i, i];
        return s;
    }

    public Complex Det()
    {
        if (!IsSquare)
            throw new InvalidOperationException("Determinant requires a square matrix.");

        return _m.Determinant();
    }

    // VectorC to MatrixC -------------------------------------------------------
    public static MatrixC Diagonal(Complex[] diag)
        => new MatrixC(Matrix<Complex>.Build.DiagonalOfDiagonalArray(diag));

    public static MatrixC Diagonal(VectorC v) => Diagonal(v.ToArray());

    public static MatrixC FromVectors(VectorC[] vs)
    {
        int l = vs.Length;
        int s = vs[0].Size;

        MatrixC M = new MatrixC(s, l);

        for (int i = 0; i < l; i++)
            M.SetColumn(i, vs[i]);

        return M;
    }

    public static VectorC[] ToVectors(MatrixC M)
    {
        int rows = M.Rows;
        int cols = M.Cols;

        VectorC[] vs = new VectorC[cols];

        for (int j = 0; j < cols; j++)
            vs[j] = M.GetColumn(j);

        return vs;
    }

    public VectorC ToVectorC()
    {
        if (Cols != 1)
            throw new InvalidOperationException("ToVectorC requires a single-column matrix.");
        return new VectorC(_m);
    }

    // In‑place arithmetic for MatrixC (MKL-backed)
    public void AddIP(MatrixC other)
    {
        if (Rows != other.Rows || Cols != other.Cols)
            throw new ArgumentException("Size mismatch.");
        _m.Add(other._m, _m);
    }

    public void NegateIP()
    {
        _m.Multiply(-1.0, _m);
    }

    public void SubIP(MatrixC other)
    {
        if (Rows != other.Rows || Cols != other.Cols)
            throw new ArgumentException("Size mismatch.");
        _m.Subtract(other._m, _m);
    }

    public void ScaleIP(Complex s)
    {
        _m.Multiply(s, _m);
    }

    public void DivideIP(Complex s)
    {
        _m.Multiply(1.0 / s, _m);
    }

    public void MultiplyIP(VectorC x, VectorC y)
    {
        if (Cols != x.Size || Rows != y.Size)
            throw new ArgumentException("Size mismatch.");
        _m.Multiply(x.Inner, y.Inner);
    }

    public void MultiplyIP(MatrixC B, MatrixC C)
    {
        if (Cols != B.Rows || Rows != C.Rows || B.Cols != C.Cols)
            throw new ArgumentException("Size mismatch.");
        _m.Multiply(B._m, C._m); 
    }

    public void HadamardIP(MatrixC other)
    {
        if (Rows != other.Rows || Cols != other.Cols)
            throw new ArgumentException("Size mismatch.");
        _m.PointwiseMultiply(other._m, _m);
    }

    public void HermitianIP()
    {
        if (!IsSquare)
            throw new InvalidOperationException("In-place Hermitian requires square matrix.");

        var h = _m.ConjugateTranspose();
        _m.SetSubMatrix(0, 0, h);  
    }

    public void AxPyIP(Complex s, MatrixC B)
    {
        if (Rows != B.Rows || Cols != B.Cols)
            throw new ArgumentException("Size mismatch.");

        MatrixC temp = B.Clone();
        temp.ScaleIP(s);

        AddIP(temp);
    }


    // C = alpha * A * B + beta * C -------------------------------------------------------
    public void GemmIP(MatrixC B, MatrixC C, Complex alpha, Complex beta)
    {
        if (Cols != B.Rows || Rows != C.Rows || B.Cols != C.Cols)
            throw new ArgumentException("Dimension mismatch in GemmIP.");

        C.ScaleIP(beta);

        MatrixC temp = new MatrixC(Rows, B.Cols);
        this.MultiplyIP(B, temp);

        temp.ScaleIP(alpha);
        C.AddIP(temp);
    }


    // Safe Eval -------------------------------------------------------
    public static (MatrixC, bool) SafeEval(Func<Complex, MatrixC> F, Complex z)
    {
        try
        {
            MatrixC value = F(z);

            int r = value.Rows;
            int c = value.Cols;

            for (int i = 0; i < r; i++)
            {
                for (int j = 0; j < c; j++)
                {
                    if (Complex.IsNaN(value[i, j]) || Complex.IsInfinity(value[i, j]))
                        return (value, true);
                }
            }

            return (value, false);
        }
        catch
        {
            try
            {
                MatrixC probe = F(ScalarC.RandomComplex(z, 1.0));

                int r = probe.Rows;
                int c = probe.Cols;

                return (MatrixC.NaN(r, c), true);
            }
            catch
            {
                return (null!, true);
            }
        }
    }

    public static (MatrixC[], bool) SafeEval(Func<Complex, MatrixC[]> F, Complex z)
    {
        try
        {
            MatrixC[] value = F(z);

            int l = value.Length;
            int r = value[0].Rows;
            int c = value[0].Cols;

            for (int k = 0; k < l; k++)
            {
                MatrixC M = value[k];

                for (int i = 0; i < r; i++)
                {
                    for (int j = 0; j < c; j++)
                    {
                        if (Complex.IsNaN(M[i, j]) || Complex.IsInfinity(M[i, j]))
                            return (value, true);
                    }
                }
            }

            return (value, false);
        }
        catch
        {
            try
            {
                MatrixC[] probe = F(ScalarC.RandomComplex(z, 1.0));

                int l = probe.Length;
                int r = probe[0].Rows;
                int c = probe[0].Cols;

                return (MatrixC.NaN(l, r, c), true);
            }
            catch
            {
                return (null!, true);
            }
        }
    }


    // LU Factorization ------------------------------------------------
    public LUFactorC LU()
    {
        if (!IsSquare)
            throw new InvalidOperationException("LU requires a square matrix.");

        return new LUFactorC(_m.LU());
    }

    // Solve -----------------------------------------------------------
    public VectorC Solve(VectorC b)
    {
        if (!IsSquare)
            throw new InvalidOperationException("Solve requires a square matrix.");

        if (Rows != b.Size)
            throw new ArgumentException("Matrix rows must match vector size.");

        return new VectorC(_m.Solve(b.Inner));
    }

    public void SolveIP(VectorC b, VectorC x)
    {
        if (!IsSquare)
            throw new InvalidOperationException("Solve requires a square matrix.");

        if (Rows != b.Size || Rows != x.Size)
            throw new ArgumentException("Vector sizes must match matrix rows.");

        _m.Solve(b.Inner, x.Inner);
    }

    public MatrixC Solve(MatrixC B)
    {
        if (!IsSquare)
            throw new InvalidOperationException("Solve requires a square matrix.");

        if (Rows != B.Rows)
            throw new ArgumentException("Matrix rows must match RHS rows.");

        return new MatrixC(_m.Solve(B.Inner));
    }

    public void SolveIP(MatrixC B, MatrixC X)
    {
        if (!IsSquare)
            throw new InvalidOperationException("Solve requires a square matrix.");

        if (Rows != B.Rows || Rows != X.Rows || B.Cols != X.Cols)
            throw new ArgumentException("Matrix dimensions must match.");

        _m.Solve(B.Inner, X.Inner);
    }

    // Pretty print -------------------------------------------------------
    public string Format(string fmt)
    {
        StringBuilder sb = new StringBuilder(Rows * Cols * 16);

        for (int i = 0; i < Rows; i++)
        {
            sb.Append('[');

            for (int j = 0; j < Cols; j++)
            {
                sb.Append(_m[i, j].Format(fmt));

                if (j < Cols - 1) sb.Append(", ");
            }

            sb.Append(']');
            if (i < Rows - 1) sb.Append('\n');
        }

        return sb.ToString();
    }

    public override string ToString() => Format(null!);
    public string ToString(string fmt) => Format(fmt);

}

public sealed class LUFactorC
{
    private readonly LU<Complex> _lu;

    internal LUFactorC(LU<Complex> lu)
    {
        _lu = lu;
    }

    // Solve Ax = b using the cached LU factorization
    public VectorC Solve(VectorC b)
        => new VectorC(_lu.Solve(b.Inner));

    public void SolveIP(VectorC b, VectorC x)
    {
        if (b.Size != x.Size)
            throw new ArgumentException("Vector sizes must match.");
        _lu.Solve(b.Inner, x.Inner);
    }

    // Solve AX = B for multiple RHS columns
    public MatrixC Solve(MatrixC B)
        => new MatrixC(_lu.Solve(B.Inner));

    public void SolveIP(MatrixC B, MatrixC X)
    {
        if (B.Rows != X.Rows || B.Cols != X.Cols)
            throw new ArgumentException("Matrix dimensions must match.");
        _lu.Solve(B.Inner, X.Inner);
    }

    public MatrixC U()
        => new MatrixC(_lu.U);

    public MatrixC L()
        => new MatrixC(_lu.L);

    public Permutation P() => _lu.P;

    public Complex Det()
    => _lu.Determinant;

    public int Sign
    {
        get
        {
            Permutation perm = _lu.P;
            int n = perm.Dimension;

            bool[] visited = new bool[n];
            int sign = 1;

            for (int i = 0; i < n; i++)
            {
                if (visited[i]) continue;

                int length = 0;
                int j = i;

                while (!visited[j])
                {
                    visited[j] = true;
                    j = perm[j];
                    length++;
                }

                if (length > 0 && length % 2 == 0)
                    sign *= -1;
            }

            return sign;
        }
    }
}