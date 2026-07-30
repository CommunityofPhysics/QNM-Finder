// File: Edge.cs

using System.Numerics;
using MathNet.Numerics;
using CLACFramework;

namespace QNMFinder;

public static class Edge
{
    public static (bool[], Func<Complex, Complex>[][], Func<Complex, Complex>[][]) SigmaMaps(Func<double, Complex>[,] rawSystemEC)
    {
        const double tol = 1e-12;
        const double huge = 1e8;

        (double leftRho, double rightRho) = RadialMap.EdgeRho;

        bool[] IsBadSigma = new bool[2];
        Func<Complex, Complex>[][] SigmaIn = new Func<Complex, Complex>[2][]; // [side][derivative]
        Func<Complex, Complex>[][] SigmaOut = new Func<Complex, Complex>[2][]; // [side][derivative]

        for (int side = 0; side < 2; side++)
        {
            double rho = (side == 0) ? leftRho : rightRho;

            Complex[,] epREC = new Complex[3, 3];

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                    epREC[i, j] = rawSystemEC[i, j](rho);
            }

            IsBadSigma[side] = false;

            // 1. A ≈ 0 AND B ≈ 0  → rawSigma returns NaN for all ω
            if (epREC[0, 2].Magnitude < tol && epREC[0, 1].Magnitude + epREC[1, 1].Magnitude < tol)
                IsBadSigma[side] = true;

            // 2. Huge C-term → σ ≈ sqrt(C/A) becomes huge for BOTH branches
            if (epREC[0, 0].Magnitude > huge || epREC[1, 0].Magnitude > huge || epREC[2, 0].Magnitude > huge)
                IsBadSigma[side] = true;

            Func<Complex, Complex>[][] rawSigmaStack = new Func<Complex, Complex>[2][]; // [branch][derivative]

            for (int branch = 0; branch < 2; branch++)
            {
                double sign = (branch == 0) ? 1.0 : -1.0;

                Func<Complex, Complex> rawSigma = omega =>
                {
                    Complex A = epREC[0, 2];
                    Complex B = epREC[0, 1] + omega * epREC[1, 1];
                    Complex C = epREC[0, 0] + omega * epREC[1, 0] + omega * omega * epREC[2, 0];

                    if (A.Magnitude < tol && B.Magnitude < tol) return Complex.NaN;

                    if (A.Magnitude < tol) return C / B;

                    Complex disc = B * B - 4.0 * A * C;
                    return (B + sign * ScalarC.SafeSqrt(disc)) / (2.0 * A);
                };

                Func<Complex, Complex> rawRSigma = omega =>
                {
                    Complex A = epREC[0, 2];
                    Complex B = epREC[0, 1] + omega * epREC[1, 1];

                    Complex RB = epREC[1, 1];
                    Complex RC = epREC[1, 0] + 2.0 * omega * epREC[2, 0];

                    if (A.Magnitude < tol && B.Magnitude < tol)
                        return Complex.NaN;

                    // Use the already-selected branch
                    Complex sigma = rawSigma(omega);
                    Complex denom = 2.0 * A * sigma - B;
                    if (denom.Magnitude < tol) return Complex.NaN; // branch point

                    return (RB * sigma - RC) / denom;
                };

                rawSigmaStack[branch] = new[] { rawSigma, rawRSigma };
            }

            // Build incoming/outgoing sigmas from branches
            Func<Complex, Complex>[] rightSigma = BranchSelector(rawSigmaStack[0], rawSigmaStack[1], true);
            Func<Complex, Complex>[] leftSigma = BranchSelector(rawSigmaStack[0], rawSigmaStack[1], false);

            SigmaIn[side] = (side == 0) ? rightSigma : leftSigma;
            SigmaOut[side] = (side == 0) ? leftSigma : rightSigma;
        }
        return (IsBadSigma, SigmaIn, SigmaOut);
    }

    private static Func<Complex, Complex>[] BranchSelector(Func<Complex, Complex>[] rawPlus, Func<Complex, Complex>[] rawMinus, bool goesRight)
    {
        Func<Complex, Complex> sigma = omega =>
        {
            Complex sp = rawPlus[0](omega);
            Complex sm = rawMinus[0](omega);

            return SelectSigmaPlus(sp, sm, omega, goesRight) ? sp : sm;
        };

        Func<Complex, Complex> rsigma = omega =>
        {
            Complex sp = rawPlus[0](omega);
            Complex sm = rawMinus[0](omega);

            return SelectSigmaPlus(sp, sm, omega, goesRight)
                ? rawPlus[1](omega) : rawMinus[1](omega);
        };

        Func<Complex, Complex>[] directionalSigma = new Func<Complex, Complex>[] { sigma, rsigma };

        return directionalSigma;
    }

    private static bool SelectSigmaPlus(Complex sigmaPlus, Complex sigmaMinus, Complex omega, bool goesRight)
    {
        double deltaPlus = sigmaPlus.Real * omega.Imaginary - sigmaPlus.Imaginary * omega.Real;
        double deltaMinus = sigmaMinus.Real * omega.Imaginary - sigmaMinus.Imaginary * omega.Real;

        bool dPP = deltaPlus > 0.0;
        bool dMP = deltaMinus > 0.0;

        if (dPP != dMP)
            return dPP == goesRight;

        if (dPP == goesRight)
            return goesRight ? deltaPlus > deltaMinus : deltaPlus < deltaMinus;

        return sigmaPlus.Real <= sigmaMinus.Real;
    }

    public static Complex[,] BranchPoints(Func<double, Complex>[,] rawSystemEC, BCSpec[] rawSystemBC)
    {
        const double tol = 1e-16;

        (double leftRho, double rightRho) = RadialMap.EdgeRho;

        Complex[,] SystemBP = new Complex[2, 2];

        for (int side = 0; side < 2; side++)
        {
            double rho = (side == 0) ? leftRho : rightRho;

            Complex[,] epREC = new Complex[3, 3];

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                    epREC[i, j] = rawSystemEC[i, j](rho);
            }

            Complex a = epREC[1, 1] * epREC[1, 1] - 4.0 * epREC[0, 2] * epREC[2, 0];
            Complex b = 2.0 * epREC[0, 1] * epREC[1, 1] - 4.0 * epREC[0, 2] * epREC[1, 0];
            Complex c = epREC[0, 1] * epREC[0, 1] - 4.0 * epREC[0, 2] * epREC[0, 0];

            if (a.Magnitude < tol && b.Magnitude < tol)
            {
                SystemBP[side, 0] = Complex.NaN;
                SystemBP[side, 1] = Complex.NaN;
            }
            else if (a.Magnitude < tol)
            {
                SystemBP[side, 0] = - c / b;
                SystemBP[side, 1] = Complex.NaN;
            }
            else
            {
                Complex disc = b * b - 4.0 * a * c;
                Complex sqrtDisc = ScalarC.SafeSqrt(disc);

                SystemBP[side, 0] = (-b + sqrtDisc) / (2.0 * a);
                SystemBP[side, 1] = (-b - sqrtDisc) / (2.0 * a);
            }

            if (rawSystemBC[side].Sigma?.Type is not
                (SigmaType.Incoming or SigmaType.Outgoing))
            {
                SystemBP[side, 0] = Complex.NaN;
                SystemBP[side, 1] = Complex.NaN;
            }
        }

        return SystemBP;
    }

}

public enum SigmaType { Custom, Incoming, Outgoing }

public sealed class SigmaSpec
{
    public SigmaType Type { get; }
    public Func<Complex, Complex> Map { get; }
    public Func<Complex, Complex> RMap { get; }

    public SigmaSpec(SigmaType type, Func<Complex, Complex> map, Func<Complex, Complex> rmap)
    {
        Type = type;
        Map = map;
        RMap = rmap;
    }

    public SigmaSpec SetType(SigmaType newType)
    {
        return new SigmaSpec(newType, Map, RMap);
    }

    public SigmaSpec SetMap(Func<Complex, Complex> newMap, Func<Complex, Complex>? newRMap)
    {
        return new SigmaSpec(Type, newMap, newRMap ?? Calculus.Derivative(newMap));
    }
}

public sealed class SigmaLib
{
    public bool IsBadSigma { get; }
    public SigmaSpec? Custom { get; }
    public SigmaSpec? Incoming { get; }
    public SigmaSpec? Outgoing { get; }

    public SigmaLib(bool isBadSigma, Func<Complex, Complex>[]? custom, Func<Complex, Complex>[]? incoming, Func<Complex, Complex>[]? outgoing)
    {
        IsBadSigma = isBadSigma;
        Custom = Store(SigmaType.Custom, custom);
        Incoming = Store(SigmaType.Incoming, incoming);
        Outgoing = Store(SigmaType.Outgoing, outgoing);

    }

    private static SigmaSpec? Store(SigmaType type, Func<Complex, Complex>[]? sigmaData)
    {
        // sigmaData[0] -> Sigma function, sigmaData[1] -> Derivative of sigma 
        if (sigmaData is null) return null;

        if (sigmaData.Length == 0)
            throw new ArgumentException("Sigma data must contain at least one mapping.");

        var sigma = sigmaData[0];
        var rsigma = (sigmaData.Length >= 2 && sigmaData[1] != null)
            ? sigmaData[1] : Calculus.Derivative(sigma);

        return new SigmaSpec(type, sigma, rsigma);
    }
}
