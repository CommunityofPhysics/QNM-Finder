// File: Boundary.cs

using System.Numerics;
using MathNet.Numerics;
using CLACFramework;

namespace QNMFinder;

public enum BCType
{
    Dirichlet, Neumann, Robin, Sommerfeld
}

public sealed class BCSpec
{
    public BCType Type { get; }
    public Complex Alpha { get; }   // α for Robin(α) or Sommerfeld(α), or (a,b) for Robin(a,b) or Sommerfeld(a,b)
    public SigmaSpec? Sigma { get; }

    public BCSpec(BCType type, Complex alpha, SigmaSpec? sigma)
    {
        Type = type;
        Alpha = alpha;
        Sigma = sigma;
    }

    public static BCSpec Dirichlet()
        => new BCSpec(BCType.Dirichlet, Complex.Zero, null);

    public static BCSpec Neumann()
    => new BCSpec(BCType.Neumann, Complex.Zero, null);

    public static BCSpec Robin(Complex alpha)
        => new BCSpec(BCType.Robin, alpha, null);

    public static BCSpec Sommerfeld(SigmaSpec sigma)
        => new BCSpec(BCType.Sommerfeld, Complex.Zero, sigma);

    public static BCSpec Parse(string spec, SigmaLib lib)
    {
        if (string.IsNullOrWhiteSpace(spec))
            throw new ArgumentException("BC specification cannot be empty.");

        spec = spec.Trim().ToLowerInvariant();

        // -------------------------
        // 1. Dirichlet (no parentheses)
        // -------------------------
        if (spec == "dirichlet")
            return Dirichlet();

        // -------------------------
        // 2. Neumann (no parentheses)
        // -------------------------
        if (spec == "neumann")
            return Neumann();

        // -------------------------
        // 3. Robin or Sommerfeld
        // -------------------------
        bool isRobin = spec.StartsWith("robin");
        bool isSommerfeld = spec.StartsWith("sommerfeld");

        if (isRobin || isSommerfeld)
        {
            int i1 = spec.IndexOf('(');
            int i2 = spec.LastIndexOf(')');

            if (i1 < 0 || i2 < 0 || i2 <= i1 + 1)
                throw new ArgumentException($"Invalid BC: {spec}");

            string inside = spec.Substring(i1 + 1, i2 - i1 - 1).Trim();

            Complex alpha;

            // -------------------------
            // Case A: Single parameter → Robin(c)
            // -------------------------
            if (!inside.Contains(","))
            {
                if (IsInfinity(inside))
                    return Dirichlet();

                if (!double.TryParse(inside, out double c))
                {
                    if (isSommerfeld)
                    {
                        switch (inside)
                        {
                            case "custom":
                                return Sommerfeld(lib.Custom!);

                            case "incoming":
                                if (lib.IsBadSigma) return Dirichlet();
                                return Sommerfeld(lib.Incoming!);

                            case "outgoing":
                                if (lib.IsBadSigma) return Dirichlet();
                                return Sommerfeld(lib.Outgoing!);
                        }
                    }

                    throw new ArgumentException($"Invalid number in BC: {inside}");
                }

                if (c == 0.0)
                    return Neumann();

                alpha = new Complex(c, 0);
            }
            else
            {
                // -------------------------
                // Case B: Two parameters → Robin(a,b) or Sommerfeld(a,b)
                // -------------------------
                string[] parts = inside.Split(',');
                if (parts.Length != 2)
                    throw new ArgumentException($"Invalid BC parameters: {spec}");

                string sa = parts[0].Trim();
                string sb = parts[1].Trim();

                if (IsInfinity(sa) || IsInfinity(sb))
                    return Dirichlet();

                if (!double.TryParse(sa, out double a) ||
                    !double.TryParse(sb, out double b))
                    throw new ArgumentException($"Invalid numbers in BC: {inside}");

                if (a == 0.0 && b == 0.0)
                    return Neumann();

                alpha = new Complex(a, b);   // store (a,b)
            }

            // -------------------------
            // 3. Return correct BC type
            // -------------------------
            if (isRobin)
                return Robin(alpha);

            SigmaSpec sigma = new SigmaSpec(SigmaType.Custom, om => om * alpha, om => alpha);

            return Sommerfeld(sigma);
        }

        // -------------------------
        // 4. Unsupported BC
        // -------------------------
        throw new ArgumentException("Only Dirichlet, Neumann, Robin(...), or Sommerfeld(...) are supported.");
    }
    private static bool IsInfinity(string s)
    {
        s = s.Trim().ToLowerInvariant();
        return s == "inf" || s == "infinity";
    }
}

public static class Boundary
{
    public static (Func<Complex, MatrixC> R, Func<Complex, MatrixC> S, Func<Complex, MatrixC> DR, Func<Complex, MatrixC> DS)
        ImposeBC(Func<Complex, MatrixC> P, Func<Complex, MatrixC> Q, Func<Complex, MatrixC> RP, Func<Complex, MatrixC> RQ, BCSpec[] systemBC)
    {
        if (P(Complex.Zero).Size != Q(Complex.Zero).Size)
            throw new ArgumentException("Matrix function sizes must match.");

        int NP = P(Complex.Zero).Rows;
        int N = NP - 1;

        Func<Complex, MatrixC> R = P;
        Func<Complex, MatrixC> S = Q;
        Func<Complex, MatrixC> RR = RP;
        Func<Complex, MatrixC> RS = RQ;

        foreach (int i in new[] { 0, N })
        {
            BCSpec sideBC = systemBC[i / N];

            SigmaSpec? sigma = sideBC.Sigma;

            Func<Complex, MatrixC> tempR = R;
            R = omega =>
            {
                MatrixC RM = tempR(omega).Clone();
                RM.ZeroRow(i);

                switch (sideBC.Type)
                {
                    case BCType.Dirichlet:
                        RM[i, i] = Complex.One;
                        break;

                    case BCType.Neumann:
                        break;

                    case BCType.Robin:
                        RM[i, i] = sideBC.Alpha;
                        break;

                    case BCType.Sommerfeld:
                        RM[i, i] = sigma!.Map(omega);
                        break;
                }

                return RM;
            };

            Func<Complex, MatrixC> tempS = S;
            S = omega =>
            {
                MatrixC SM = tempS(omega).Clone();
                SM.ZeroRow(i);

                switch (sideBC.Type)
                {
                    case BCType.Dirichlet:
                        break;

                    case BCType.Neumann:
                        SM[i, i] = Complex.One;
                        break;

                    case BCType.Robin:
                        SM[i, i] = Complex.One;
                        break;

                    case BCType.Sommerfeld:
                        SM[i, i] = Complex.One;
                        break;
                }

                return SM;
            };

            Func<Complex, MatrixC> tempRR = RR;
            RR = omega =>
            {
                MatrixC DRM = tempRR(omega).Clone();
                DRM.ZeroRow(i);

                if (sideBC.Type == BCType.Sommerfeld)
                    DRM[i, i] = sigma!.RMap(omega);

                return DRM;
            };

            Func<Complex, MatrixC> tempRS = RS;
            RS = omega =>
            {
                MatrixC DSM = tempRS(omega).Clone();
                DSM.ZeroRow(i);

                return DSM;
            };
        }

        return (R, S, RR, RS);
    }

    public static (Func<Complex, ActionC> R, Func<Complex, ActionC> S, Func<Complex, ActionC> DR, Func<Complex, ActionC> DS)
        ImposeBC(Func<Complex, ActionC> P, Func<Complex, ActionC> Q, Func<Complex, ActionC> RP, Func<Complex, ActionC> RQ, BCSpec[] systemBC)
    {
        if (P(Complex.Zero).Size != Q(Complex.Zero).Size) throw new ArgumentException("Operator sizes must match.");

        int NP = P(Complex.Zero).Rows;
        int N = NP - 1;

        Func<Complex, ActionC> R = P;
        Func<Complex, ActionC> S = Q;
        Func<Complex, ActionC> RR = RP;
        Func<Complex, ActionC> RS = RQ;

        foreach (int i in new[] { 0, N })
        {
            BCSpec sideBC = systemBC[i / N];

            SigmaSpec? sigma = sideBC.Sigma;

            Func<Complex, ActionC> tempR = R;
            R = omega =>
            {
                ActionC RA = tempR(omega);

                return new ActionC(NP, NP, (VectorC y, VectorC Ry) =>
                {
                    RA.Apply(y, Ry);
                    Ry[i] = Complex.Zero;

                    switch (sideBC.Type)
                    {
                        case BCType.Dirichlet:
                            Ry[i] = y[i];
                            break;

                        case BCType.Neumann:
                            break;

                        case BCType.Robin:
                            Ry[i] = sideBC.Alpha * y[i];
                            break;

                        case BCType.Sommerfeld:
                            Ry[i] = sigma!.Map(omega) * y[i];
                            break;
                    }
                });
            };

            Func<Complex, ActionC> tempS = S;
            S = omega =>
            {
                ActionC SA = tempS(omega);

                return new ActionC(NP, NP, (VectorC z, VectorC Sz) =>
                {
                    SA.Apply(z, Sz);
                    Sz[i] = Complex.Zero;

                    switch (sideBC.Type)
                    {
                        case BCType.Dirichlet:
                            break;

                        case BCType.Neumann:
                            Sz[i] = z[i];
                            break;

                        case BCType.Robin:
                            Sz[i] = z[i];
                            break;

                        case BCType.Sommerfeld:
                            Sz[i] = z[i];
                            break;
                    }
                });
            };

            Func<Complex, ActionC> tempRR = RR;
            RR = omega =>
            {
                ActionC DRA = tempRR(omega);

                return new ActionC(NP, NP, (VectorC y, VectorC DRy) =>
                {
                    DRA.Apply(y, DRy);
                    DRy[i] = Complex.Zero;

                    if (sideBC.Type == BCType.Sommerfeld)
                        DRy[i] = sigma!.RMap(omega) * y[i];

                });
            };

            Func<Complex, ActionC> tempRS = RS;
            RS = omega =>
            {
                ActionC DSA = tempRS(omega);

                return new ActionC(NP, NP, (VectorC z, VectorC DSz) =>
                {
                    DSA.Apply(z, DSz);
                    DSz[i] = Complex.Zero;

                });
            };

        }

        return (R, S, RR, RS);
    }

}
