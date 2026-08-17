// File: Calculus.cs

using System;
using System.Numerics;

namespace CLACFramework;

public static class Calculus
{
    // ============================================================
    // Derivatives
    // ============================================================
    public static readonly Complex Epsilon = new Complex(1e-4, 0.0);

    // C -> C
    public static Func<Complex, Complex> Derivative(Func<Complex, Complex> F)
    {
        return z =>
        {
            Complex Fb2 = F(z - 2.0 * Epsilon);
            Complex Fb1 = F(z - Epsilon);
            Complex Fc = F(z);
            Complex Ff1 = F(z + Epsilon);
            Complex Ff2 = F(z + 2.0 * Epsilon);

            Complex Fb2S = Fc + ScalarC.PUDifference(Fb2, Fc);
            Complex Fb1S = Fc + ScalarC.PUDifference(Fb1, Fc);
            Complex Ff1S = Fc + ScalarC.PUDifference(Ff1, Fc);
            Complex Ff2S = Fc + ScalarC.PUDifference(Ff2, Fc);

            Complex DF = (Fb2S - 8 * Fb1S + 8 * Ff1S - Ff2S) / (12.0 * Epsilon);

            return DF;
        };
    }

    // VectorC -> VectorC
    public static Func<Complex, VectorC> Derivative(Func<Complex, VectorC> F)
    {
        return z =>
        {
            VectorC Fb2 = F(z - 2.0 * Epsilon);
            VectorC Fb1 = F(z - Epsilon);
            VectorC Ff1 = F(z + Epsilon);
            VectorC Ff2 = F(z + 2.0 * Epsilon);

            VectorC DF = (Fb2 - 8 * Fb1 + 8 * Ff1 - Ff2) / (12.0 * Epsilon);

            return DF;
        };
    }

    // MatrixC -> MatrixC
    public static Func<Complex, MatrixC> Derivative(Func<Complex, MatrixC> F)
    {
        return z =>
        {
            MatrixC Fb2 = F(z - 2.0 * Epsilon);
            MatrixC Fb1 = F(z - Epsilon);
            MatrixC Ff1 = F(z + Epsilon);
            MatrixC Ff2 = F(z + 2.0 * Epsilon);

            MatrixC DF = (Fb2 - 8 * Fb1 + 8 * Ff1 - Ff2) / (12.0 * Epsilon);

            return DF;
        };
    }

    // ActionC -> ActionC
    public static Func<Complex, ActionC> Derivative(Func<Complex, ActionC> F)
    {
        return z =>
        {
            ActionC Fb2 = F(z - 2.0 * Epsilon);
            ActionC Fb1 = F(z - Epsilon);
            ActionC Ff1 = F(z + Epsilon);
            ActionC Ff2 = F(z + 2.0 * Epsilon);

            ActionC DF = (Fb2 - 8 * Fb1 + 8 * Ff1 - Ff2) / (12.0 * Epsilon);

            return DF;
        };
    }

    // All operator‑valued functions T(ω) used in this framework are assumed to be continuous in ω
    // and free of branch discontinuities. Therefore, phase‑unwrapping is only applied to scalar functions.

    // ============================================================
    // Loop Integrations
    // ============================================================

    // C -> C
    public static (Complex, List<Complex>) LoopIntegral(Func<Complex, Complex> Integrand, SqDomainC domain, int M)
    {
        double h = 1.0 / M;
        double tb = 0.0;

        (Complex zb, _) = domain.Parametrize(tb);
        (Complex Fb, bool IsFbNaN) = ScalarC.SafeEval(Integrand, zb);

        List<Complex> Singular = new List<Complex>();

        if (IsFbNaN)
        {
            Singular.Add(zb);
            return (Complex.NaN, Singular);
        }

        Complex sum = Complex.Zero;

        for (int i = 0; i < M; i++)
        {
            double tf = tb + h;

            (Complex zf, _) = domain.Parametrize(tf);
            (Complex Ff, bool IsFfNaN) = ScalarC.SafeEval(Integrand, zf);

            if (IsFfNaN)
            {
                Singular.Add(zf);
                return (Complex.NaN, Singular);
            }

            Complex dz = zf - zb;

            sum += 0.5 * (Fb + Ff) * dz;

            tb = tf;
            zb = zf;
            Fb = Ff;
        }

        return (sum, Singular);
    }

    public static (Complex[], List<Complex>) LoopIntegral(Func<Complex, Complex[]> Integrand, SqDomainC domain, int M)
    {
        double h = 1.0 / M;
        double tb = 0.0;

        (Complex zb, _) = domain.Parametrize(tb);
        (Complex[] Fb, bool IsFbNaN) = ScalarC.SafeEval(Integrand, zb);

        int L = Fb.Length;
        List<Complex> Singular = new List<Complex>();

        if (IsFbNaN)
        {
            Singular.Add(zb);
            return (ScalarC.NaN(L), Singular);
        }

        Complex[] sum = new Complex[L];

        for (int i = 0; i < M; i++)
        {
            double tf = tb + h;

            (Complex zf, _) = domain.Parametrize(tf);
            (Complex[] Ff, bool IsFfNaN) = ScalarC.SafeEval(Integrand, zf);

            if (IsFfNaN)
            {
                Singular.Add(zf);
                return (ScalarC.NaN(L), Singular);
            }

            Complex dz = zf - zb;

            for (int l = 0; l < L; l++)
            {
                sum[l] += 0.5 * (Fb[l] + Ff[l]) * dz;
            }

            tb = tf;
            zb = zf;
            Fb = Ff;
        }

        return (sum, Singular);
    }

    // VectorC -> VectorC
    public static (VectorC, List<Complex>) LoopIntegral(Func<Complex, VectorC> Integrand, SqDomainC domain, int M)
    {
        double h = 1.0 / M;
        double tb = 0.0;

        (Complex zb, _) = domain.Parametrize(tb);
        (VectorC Fb, bool IsFbNaN) = VectorC.SafeEval(Integrand, zb);

        int S = Fb.Size;
        List<Complex> Singular = new List<Complex>();

        if (IsFbNaN)
        {
            Singular.Add(zb);
            return (VectorC.NaN(S), Singular);
        }

        VectorC sum = VectorC.Zero(S);

        for (int i = 0; i < M; i++)
        {
            double tf = tb + h;

            (Complex zf, _) = domain.Parametrize(tf);
            (VectorC Ff, bool IsFfNaN) = VectorC.SafeEval(Integrand, zf);

            if (IsFfNaN)
            {
                Singular.Add(zf);
                return (VectorC.NaN(S), Singular);
            }

            Complex dz = zf - zb;

            sum += 0.5 * (Fb + Ff) * dz;

            tb = tf;
            zb = zf;
            Fb = Ff;
        }

        return (sum, Singular);
    }

    public static (VectorC[], List<Complex>) LoopIntegral(Func<Complex, VectorC[]> Integrand, SqDomainC domain, int M)
    {
        double h = 1.0 / M;
        double tb = 0.0;

        (Complex zb, _) = domain.Parametrize(tb);
        (VectorC[] Fb, bool IsFbNaN) = VectorC.SafeEval(Integrand, zb);

        int L = Fb.Length;
        int S = Fb[0].Size;
        List<Complex> Singular = new List<Complex>();

        if (IsFbNaN)
        {
            Singular.Add(zb);
            return (VectorC.NaN(L, S), Singular);
        }

        VectorC[] sum = VectorC.Zero(L, S);

        for (int i = 0; i < M; i++)
        {
            double tf = tb + h;

            (Complex zf, _) = domain.Parametrize(tf);
            (VectorC[] Ff, bool IsFfNaN) = VectorC.SafeEval(Integrand, zf);

            if (IsFfNaN)
            {
                Singular.Add(zf);
                return (VectorC.NaN(L, S), Singular);
            }

            Complex dz = zf - zb;

            for (int l = 0; l < L; l++)
            {
                sum[l] += 0.5 * (Fb[l] + Ff[l]) * dz;
            }

            tb = tf;
            zb = zf;
            Fb = Ff;
        }

        return (sum, Singular);
    }

    // MatrixC -> MatrixC
    public static (MatrixC, List<Complex>) LoopIntegral(Func<Complex, MatrixC> Integrand, SqDomainC domain, int M)
    {
        double h = 1.0 / M;
        double tb = 0.0;

        (Complex zb, _) = domain.Parametrize(tb);
        (MatrixC Fb, bool IsFbNaN) = MatrixC.SafeEval(Integrand, zb);

        int R = Fb.Rows;
        int C = Fb.Cols;

        List<Complex> Singular = new List<Complex>();

        if (IsFbNaN)
        {
            Singular.Add(zb);
            return (MatrixC.NaN(R, C), Singular);
        }

        MatrixC sum = MatrixC.Zero(R, C);

        for (int i = 0; i < M; i++)
        {
            double tf = tb + h;

            (Complex zf, _) = domain.Parametrize(tf);
            (MatrixC Ff, bool IsFfNaN) = MatrixC.SafeEval(Integrand, zf);

            if (IsFfNaN)
            {
                Singular.Add(zf);
                return (MatrixC.NaN(R, C), Singular);
            }

            Complex dz = zf - zb;

            sum += 0.5 * (Fb + Ff) * dz;

            tb = tf;
            zb = zf;
            Fb = Ff;
        }

        return (sum, Singular);
    }

    public static (MatrixC[], List<Complex>) LoopIntegral(Func<Complex, MatrixC[]> Integrand, SqDomainC domain, int M)
    {
        double h = 1.0 / M;
        double tb = 0.0;

        (Complex zb, _) = domain.Parametrize(tb);
        (MatrixC[] Fb, bool IsFbNaN) = MatrixC.SafeEval(Integrand, zb);

        int L = Fb.Length;
        int R = Fb[0].Rows;
        int C = Fb[0].Cols;

        List<Complex> Singular = new List<Complex>();

        if (IsFbNaN)
        {
            Singular.Add(zb);
            return (MatrixC.NaN(L, R, C), Singular);
        }

        MatrixC[] sum = MatrixC.Zero(L, R, C);

        for (int i = 0; i < M; i++)
        {
            double tf = tb + h;

            (Complex zf, _) = domain.Parametrize(tf);
            (MatrixC[] Ff, bool IsFfNaN) = MatrixC.SafeEval(Integrand, zf);

            if (IsFfNaN)
            {
                Singular.Add(zf);
                return (MatrixC.NaN(L, R, C), Singular);
            }

            Complex dz = zf - zb;

            for (int l = 0; l < L; l++)
            {
                sum[l] += 0.5 * (Fb[l] + Ff[l]) * dz;
            }

            tb = tf;
            zb = zf;
            Fb = Ff;
        }

        return (sum, Singular);
    }

    // ActionC -> ActionC
    public static (ActionC, List<Complex>) LoopIntegral(Func<Complex, ActionC> Integrand, SqDomainC domain, int M)
    {
        double h = 1.0 / M;
        double tb = 0.0;

        (Complex zb, _) = domain.Parametrize(tb);
        (ActionC Fb, bool IsFbNaN) = ActionC.SafeEval(Integrand, zb);

        int R = Fb.Rows;
        int C = Fb.Cols;

        List<Complex> Singular = new List<Complex>();

        if (IsFbNaN)
        {
            Singular.Add(zb);
            return (ActionC.NaN(R, C), Singular);
        }

        ActionC sum = ActionC.Zero(R, C);

        for (int i = 0; i < M; i++)
        {
            double tf = tb + h;

            (Complex zf, _) = domain.Parametrize(tf);
            (ActionC Ff, bool IsFfNaN) = ActionC.SafeEval(Integrand, zf);

            if (IsFfNaN)
            {
                Singular.Add(zf);
                return (ActionC.NaN(R, C), Singular);
            }

            Complex dz = zf - zb;

            sum += 0.5 * (Fb + Ff) * dz;

            tb = tf;
            zb = zf;
            Fb = Ff;
        }

        return (sum, Singular);
    }

    public static (ActionC[], List<Complex>) LoopIntegral(Func<Complex, ActionC[]> Integrand, SqDomainC domain, int M)
    {
        double h = 1.0 / M;
        double tb = 0.0;

        (Complex zb, _) = domain.Parametrize(tb);
        (ActionC[] Fb, bool IsFbNaN) = ActionC.SafeEval(Integrand, zb);

        int L = Fb.Length;
        int R = Fb[0].Rows;
        int C = Fb[0].Cols;

        List<Complex> Singular = new List<Complex>();

        if (IsFbNaN)
        {
            Singular.Add(zb);
            return (ActionC.NaN(L, R, C), Singular);
        }

        ActionC[] sum = ActionC.Zero(L, R, C);

        for (int i = 0; i < M; i++)
        {
            double tf = tb + h;

            (Complex zf, _) = domain.Parametrize(tf);
            (ActionC[] Ff, bool IsFfNaN) = ActionC.SafeEval(Integrand, zf);

            if (IsFfNaN)
            {
                Singular.Add(zf);
                return (ActionC.NaN(L, R, C), Singular);
            }

            Complex dz = zf - zb;

            for (int l = 0; l < L; l++)
            {
                sum[l] += 0.5 * (Fb[l] + Ff[l]) * dz;
            }

            tb = tf;
            zb = zf;
            Fb = Ff;
        }

        return (sum, Singular);
    }

}