// File: Action.cs

using System;
using System.Numerics;

namespace CLACFramework;

public class ActionC
{
    private readonly Action<VectorC, VectorC> _op;

    public int Rows { get; }
    public int Cols { get; }
    public (int Rows, int Cols) Size => (Rows, Cols);
    public bool IsSquare => Rows == Cols;

    // Constructors ----------------------------------------------------

    public ActionC(int rows, int cols, Action<VectorC, VectorC> op)
    {
        if (rows <= 0 || cols <= 0)
            throw new ArgumentOutOfRangeException();

        _op = op ?? throw new ArgumentNullException(nameof(op));
        Rows = rows;
        Cols = cols;
    }

    public ActionC Clone() => new ActionC(Rows, Cols, _op);

    // Apply -----------------------------------------------------------

    public void Apply(VectorC x, VectorC y)
    {
        if (x.Size != Cols || y.Size != Rows)
            throw new ArgumentException("ActionC.Apply: size mismatch.");

        if (ReferenceEquals(x, y))
            throw new InvalidOperationException("ActionC.Apply: x and y must not alias.");

        y.Clear();
        _op(x, y);
    }

    public VectorC Apply(VectorC x)
    {
        if (x.Size != Cols)
            throw new ArgumentException("ActionC.Apply: size mismatch.");

        VectorC y = new VectorC(Rows);
        Apply(x, y);
        return y;
    }

    // Zero / Identity ------------------------------------------------

    public static ActionC Zero(int r, int c)
        => new ActionC(r, c, (x, y) => y.Clear());

    public ActionC ZeroLike()
        => Zero(Rows, Cols);

    public static ActionC Identity(int n)
        => new ActionC(n, n, (x, y) => y.CopyFrom(x));

    public ActionC IdentityLike()
    {
        if (!IsSquare)
            throw new InvalidOperationException("IdentityLike requires square operator.");
        return Identity(Rows);
    }

    public ActionC Clear()
    {
        return new ActionC(Rows, Cols,
            (input, output) =>
            {
                for (int i = 0; i < Rows; i++)
                    output[i] = Complex.Zero;
            });
    }

    public ActionC ZeroRow(int row)
    {
        return new ActionC(Rows, Cols,
            (input, output) =>
            {
                var prev = _op;

                prev(input, output);
                output[row] = Complex.Zero;
            });
    }

    public ActionC ZeroCol(int col)
    {
        return new ActionC(Rows, Cols,
            (input, output) =>
            {
                var prev = _op;
                prev(input, output);

                VectorC e = new VectorC(Cols);
                e[col] = Complex.One;

                VectorC Ae = new VectorC(Rows);
                prev(e, Ae);

                Complex scale = input[col];
                if (scale != Complex.Zero)
                {
                    for (int i = 0; i < Rows; i++)
                        output[i] -= Ae[i] * scale;
                }
            });
    }

    // From MatrixC ---------------------------------------------------

    public static ActionC FromMatrix(MatrixC M)
    {
        int r = M.Rows;
        int c = M.Cols;

        return new ActionC(r, c, (x, y) => M.MultiplyIP(x, y));
    }

    // Diagonal -------------------------------------------------------

    public static ActionC Diagonal(VectorC v)
    {
        int n = v.Size;
        VectorC diag = v.Clone();

        return new ActionC(n, n, (y, Dy) =>
        {
            Dy.CopyFrom(y);
            Dy.HadamardProductIP(diag);
        });
    }

    // Arithmetic -----------------------------------------------------

    public static ActionC operator +(ActionC A, ActionC B)
    {
        if (A.Rows != B.Rows || A.Cols != B.Cols)
            throw new ArgumentException("ActionC +: size mismatch.");

        int r = A.Rows;
        int c = A.Cols;

        return new ActionC(r, c, (x, y) =>
        {
            VectorC tmp = new VectorC(r);

            A.Apply(x, y);
            B.Apply(x, tmp);
            y.AddIP(tmp);
        });
    }

    public static ActionC operator -(ActionC A, ActionC B)
    {
        if (A.Rows != B.Rows || A.Cols != B.Cols)
            throw new ArgumentException("ActionC -: size mismatch.");

        int r = A.Rows;
        int c = A.Cols;

        return new ActionC(r, c, (x, y) =>
        {
            VectorC tmp = new VectorC(r);

            A.Apply(x, y);
            B.Apply(x, tmp);
            y.SubIP(tmp);
        });
    }

    public static ActionC operator -(ActionC A)
    {
        int r = A.Rows;
        int c = A.Cols;

        return new ActionC(r, c, (x, y) =>
        {
            A.Apply(x, y);
            y.NegateIP();
        });
    }

    public static ActionC operator *(Complex s, ActionC A)
    {
        int r = A.Rows;
        int c = A.Cols;

        return new ActionC(r, c, (x, y) =>
        {
            A.Apply(x, y);
            y.ScaleIP(s);
        });
    }

    public static ActionC operator *(ActionC A, Complex s)
        => s * A;

    public static ActionC operator /(ActionC A, Complex s)
        => (1.0 / s) * A;

    // Composition ----------------------------------------------------
    public static ActionC operator *(ActionC A, ActionC B)
    {
        if (A.Cols != B.Rows)
            throw new ArgumentException("ActionC *: inner size mismatch.");

        int r = A.Rows;
        int c = B.Cols;

        int n = A.Cols;


        return new ActionC(r, c, (x, y) =>
        {
            VectorC tmp = new VectorC(n);
            B.Apply(x, tmp);
            A.Apply(tmp, y);
        });
    }

    // Norm Estimate --------------------------------------------------
    public double FrobeniusNorm()
    {
        int samples = 50;
        int r = Rows;
        int c = Cols;

        double total = 0.0;
        object totalLock = new object();

        Parallel.For(0, samples, () => 0.0,
            (s, state, localSum) =>
            {
                VectorC z = new VectorC(c);
                VectorC Tz = new VectorC(r);

                Random rng = new Random(unchecked(Environment.TickCount * 31 + s));

                for (int i = 0; i < c; i++)
                    z[i] = ((rng.Next() & 1) == 0) ? Complex.One : -Complex.One;

                Apply(z, Tz);

                double nrm = Tz.Norm();
                localSum += nrm * nrm;

                return localSum;
            },

            localSum => { lock (totalLock) total += localSum; }
            
        );

        return Math.Sqrt(total / samples);
    }

    // SafeEval -------------------------------------------------------

    public static (ActionC, bool) SafeEval(Func<Complex, ActionC> F, Complex z)
    {
        try { return (F(z), false); }
        catch { return (null!, true); }
    }

    public static (ActionC[], bool) SafeEval(Func<Complex, ActionC[]> F, Complex z)
    {
        try
        {
            ActionC[] probe = F(Complex.Zero);
            int len = probe.Length;
            int r = probe[0].Rows;
            int c = probe[0].Cols;

            ActionC[] arr = F(z);

            if (arr.Length != len)
                return (new ActionC[len], true);

            for (int k = 0; k < len; k++)
                if (arr[k].Rows != r || arr[k].Cols != c)
                    return (new ActionC[len], true);

            return (arr, false);
        }
        catch
        {
            return (null!, true);
        }
    }

    // Solve -------------------------------------------------------
    public VectorC Solve(VectorC b, int depth)
    {
        MatrixC B = b.ToMatrixC();
        MatrixC X = new MatrixC(Rows, 1);
        Solver.RecursiveSolve(this, B, X, depth);

        return X.ToVectorC();
    }

    public MatrixC Solve(MatrixC B, int depth)
    {
        MatrixC X = new MatrixC(Rows, B.Cols);
        Solver.RecursiveSolve(this, B, X, depth);

        return X;
    }
    
    public void SolveIP(VectorC b, VectorC x, int depth)
    {
        MatrixC B = b.ToMatrixC();
        MatrixC X = x.ToMatrixC();
        Solver.RecursiveSolve(this, B, X, depth);
    }

    public void SolveIP(MatrixC B, MatrixC X, int depth)
    {
        Solver.RecursiveSolve(this, B, X, depth);
    }

    // Pretty print ---------------------------------------------------

    public string Format(string fmt)
    {
        BSSpec spec = new BSSpec(true, 0, (bi, bj) => BSType.InMemory);

        BlockC blkA = Blockify.Action(this, spec);
        MatrixC A = blkA.Assemble();

        return A.Format(fmt);
    }
    public override string ToString() => Format(null);
    public string ToString(string fmt) => Format(fmt);

}