// File: Discretizer.cs

using IOFramework;
using MathNet.Numerics;
using System.Globalization;
using System.Numerics;

namespace CLACFramework;

public static class Discretizer
{
    // Chebyshev Nodes
    public static VectorC ChebNodesX(int N)
    {
        int NP = N + 1;

        Complex[] x = new Complex[NP];
        double factor = Math.PI / N;
        for (int j = 0; j < NP; j++)
            x[j] = - Math.Cos(factor * j);

        return new VectorC(x);
    }

    // Discretized Functions of Rho evaluated on Chebyshev Nodes
    public static VectorC VectorF(VectorC x, Func<double, Complex> F)
    {
        int NP = x.Size;

        Complex[] VF = new Complex[NP];
        for (int j = 0; j < NP; j++)
        {
            double rho = RadialMap.Rho(x[j].Real);
            VF[j] = F(rho);
        }

        return new VectorC(VF);
    }

    public static MatrixC MatrixF(VectorC x, Func<double, Complex> F)
    {
        return MatrixC.Diagonal(VectorF(x, F));
    }

    public static ActionC ActionF(VectorC x, Func<double, Complex> F)
    {
        return ActionC.Diagonal(VectorF(x, F));
    }

    public static Func<Complex, VectorC> VectorF(VectorC x, Func<double, Complex, Complex> F)
        => omega => VectorF(x, rho => F(rho, omega));

    public static Func<Complex, MatrixC> MatrixF(VectorC x, Func<double, Complex, Complex> F)
        =>omega => MatrixF(x, rho => F(rho, omega));
    
    public static Func<Complex, ActionC> ActionF(VectorC x, Func<double, Complex, Complex> F)
        => omega => ActionF(x, rho => F(rho, omega));

    // Discretized Jacobian 
    public static VectorC ChebJacobXV(VectorC x)
    {
        int NP = x.Size;

        Complex[] jacobian = new Complex[NP];
        for (int j = 0; j < NP; j++)
            jacobian[j] = RadialMap.JRhoX(x[j].Real);

        return new VectorC(jacobian);
    }

    public static MatrixC ChebJacobXM(VectorC x)
    {
        return MatrixC.Diagonal(ChebJacobXV(x));
    }

    public static ActionC ChebJacobXA(VectorC x)
    {
        return ActionC.Diagonal(ChebJacobXV(x));
    }

    // Discretized Inverse Jacobian
    public static VectorC ChebInvJacobXV(VectorC x)
    {
        int NP = x.Size;

        Complex[] invJacobian = new Complex[NP];
        for (int j = 0; j < NP; j++)
            invJacobian[j] = 1.0 / RadialMap.JRhoX(x[j].Real);

        return new VectorC(invJacobian);
    }

    public static MatrixC ChebInvJacobXM(VectorC x)
    {
        return MatrixC.Diagonal(ChebInvJacobXV(x));
    }

    public static ActionC ChebInvJacobXA(VectorC x)
    {
        return ActionC.Diagonal(ChebInvJacobXV(x));
    }

    // Discretized Differentiation Operators
    public static MatrixC ChebDiffXM(VectorC x)
    {
        int NP = x.Size;
        int N = NP - 1;

        MatrixC DX = MatrixC.Zero(NP, NP);

        double[] c = new double[NP];
        for (int j = 0; j < NP; j++)
        {
            c[j] = ((j == 0 || j == N) ? 2.0 : 1.0) * ((j % 2 == 1) ? -1.0 : 1.0);
        }

        Complex[] xArr = x.ToArray();

        for (int i = 0; i < NP; i++)
        {
            Complex xArri = xArr[i];
            double ci = c[i];

            for (int j = 0; j < NP; j++)
            {
                if (i != j) DX[i, j] = (ci / c[j]) / (xArri - xArr[j]);
            }
        }

        for (int i = 0; i < NP; i++)
        {
            Complex sum = Complex.Zero;
            for (int j = 0; j < NP; j++)
            {
                if (i != j) sum += DX[i, j];
            }
            DX[i, i] = -sum;
        }

        return DX;
    }

    public static ActionC ChebDiffXA(VectorC x)
    {
        int NP = x.Size;
        int N = NP - 1;

        double[] c = new double[NP];
        for (int j = 0; j < NP; j++)
            c[j] = ((j == 0 || j == N) ? 2.0 : 1.0) * ((j % 2 == 1) ? -1.0 : 1.0);

        Complex[] xArr = x.ToArray();

        ActionC D = new ActionC(NP, NP, (VectorC y, VectorC DXy) =>
        {
            Complex[] yArr = y.ToArray();

            for (int i = 0; i < NP; i++)
            {
                Complex xArri = xArr[i];
                Complex yArri = yArr[i];
                double ci = c[i];
                Complex sum = Complex.Zero;

                for (int j = 0; j < NP; j++)
                {
                    if (i == j) continue;
                    sum += (ci / c[j]) * (yArr[j] - yArri) / (xArri - xArr[j]);
                }

                DXy[i] = sum;
            }
        });

        return D;
    }

    // Rho Derivative on Chebyshev Nodes
    public static (VectorC x, MatrixC D) DiscretizeM(int N)
    {
        VectorC x = Discretizer.ChebNodesX(N);
        MatrixC DX = Discretizer.ChebDiffXM(x);             // d/dx
        MatrixC InvJ = Discretizer.ChebInvJacobXM(x);       // dx/dρ
        MatrixC D = InvJ * DX;                              // d/dρ

        return (x, D);
    }

    public static (VectorC x, ActionC D) DiscretizeA(int N)
    {
        VectorC x = Discretizer.ChebNodesX(N);
        ActionC DX = Discretizer.ChebDiffXA(x);             // d/dx
        ActionC InvJ = Discretizer.ChebInvJacobXA(x);       // dx/dρ
        ActionC D = InvJ * DX;                              // d/dρ

        return (x, D);
    }

    // Chebyshev Barycentric Interpolator
    public static VectorC Interpolate(VectorC VFC, int L)
    {
        int KP = VFC.Size;   // coarse nodes
        int K = KP - 1;      // coarse intervals
        int N = K * L;       // fine intervals
        int NP = N + 1;        // fine nodes

        VectorC xC = Discretizer.ChebNodesX(K);   // coarse Chebyshev nodes
        VectorC x = Discretizer.ChebNodesX(N);   // fine Chebyshev nodes

        VectorC VF = VectorC.Zero(NP);

        // Barycentric weight
        Complex[] w = new Complex[KP];
        for (int k = 0; k < KP; k++)
        {
            w[k] = (k % 2 == 0) ? Complex.One : -Complex.One;
            if (k == 0 || k == K)
                w[k] *= 0.5;
        }

        for (int i = 0; i < NP; i++)
        {
            int j = i % L;          // offset inside subinterval
            int k = (i - j) / L;    // coarse interval index

            if (j == 0)
            {
                // exactly on a coarse node
                VF[i] = VFC[k];
            }
            else
            {
                // Chebyshev–Lobatto barycentric interpolation
                Complex num = Complex.Zero;
                Complex denom = Complex.Zero;

                for (int l = 0; l < KP; l++)
                {
                    num += VFC[l] * w[l] / (x[i] - xC[l]);
                    denom += w[l] / (x[i] - xC[l]);
                }

                VF[i] = num / denom;
            }
        }

        return VF;
    }

}

public enum MapMode
{
    // Onesided: 0 ≤ r ≤ +Range // Twosided: -Range ≤ r ≤ +Range
    OneSided, TwoSided
}

public static class RadialMap
{
    public static MapMode Mode { get; set;  }
    public static (double Left, double Right) EdgeRho { get; private set; }

    public static double Scale { get; set; }
    public static double Range { get; set; }
    private static double Factor { get; set; }

    private static readonly double MinScale = 1e-16;
    private static readonly double MaxRange = 1e16;

    // Configuration
    public static void Configure(string mspec)
    {
        (Mode, Scale, Range) = Parse(mspec);

        switch (Mode)
        {
            case MapMode.OneSided:
                EdgeRho = (MinScale, Range);
                break;

            case MapMode.TwoSided:
                EdgeRho = (-Range, +Range);
                break;
        }

        Factor = Math.Asinh(Range / Scale);
    }

    // Parsing
    public static (MapMode, double, double) Parse(string spec)
    {
        if (string.IsNullOrWhiteSpace(spec)) throw new Exception("RadialMap specification is empty.");

        spec = spec.Trim();

        int p1 = spec.IndexOf('(');
        int p2 = spec.LastIndexOf(')');

        if (p1 < 0 || p2 < 0 || p2 <= p1)
            throw new Exception($"Invalid RadialMap specification: {spec}");

        string modeStr = new string(spec.Substring(0, p1).Trim().ToLowerInvariant()
            .Where(c => !char.IsWhiteSpace(c)).ToArray());

        MapMode mode = modeStr switch
        {
            "onesided" or "finite" or "infinite" => MapMode.OneSided,
            "twosided" or "bifinite" or "biinfinite" => MapMode.TwoSided,
            _ => throw new Exception($"Unsupported map mode: '{modeStr}'")
        };

        string parameters = spec.Substring(p1 + 1, p2 - p1 - 1).Trim();

        if (string.IsNullOrWhiteSpace(parameters)) throw new Exception($"{mode} requires at least a scale parameter.");

        string[] parts = parameters.Split(',', StringSplitOptions.RemoveEmptyEntries);

        if (!double.TryParse(parts[0].Trim(), NumberStyles.Float, CultureInfo.InvariantCulture, out double scale))
            throw new Exception($"Invalid scale value: '{parts[0]}'");

        double range;

        if (parts.Length == 1)
        {
            range = MaxRange;
        }
        else if (parts.Length == 2)
        {
            string rangeStr = parts[1].Trim().ToLowerInvariant();

            if (double.TryParse(rangeStr, NumberStyles.Float, CultureInfo.InvariantCulture, out double parsedRange))
                range = Math.Min(MaxRange, Math.Abs(parsedRange));

            else if (rangeStr == "inf" || rangeStr == "infinity")
                range = MaxRange;

            else throw new Exception($"Invalid range value: '{rangeStr}'");
        }
        else
        {
            throw new Exception($"{mode} accepts at most two parameters: [scale, range]");
        }

        scale = Math.Max(MinScale, Math.Abs(scale));
        range = Math.Max(scale, range);
        return (mode, scale, range);
    }

    // ρ(x)
    public static double Rho(double x) => Mode switch
    {
        MapMode.OneSided => Scale * Math.Sinh(0.5 * Factor * (1.0 + x)),
        MapMode.TwoSided => Scale * Math.Sinh(Factor * x),
        _ => throw Unsupported()
    };

    // dρ/dx
    public static double JRhoX(double x) => Mode switch
    {
        MapMode.OneSided => 0.5 * Factor * Scale * Math.Cosh(0.5 * Factor * (1.0 + x)),
        MapMode.TwoSided => Factor * Scale * Math.Cosh(Factor * x),
        _ => throw Unsupported()
    };

    // x(ρ)
    public static double X(double rho) => Mode switch
    {
        MapMode.OneSided => (2.0 / Factor) * Math.Asinh(rho / Scale) - 1.0,
        MapMode.TwoSided => (1.0 / Factor) * Math.Asinh(rho / Scale),
        _ => throw Unsupported()
    };

    // dx/dρ
    public static double JXRho(double rho) => Mode switch
    {
        MapMode.OneSided => 2.0 / (Factor * Math.Sqrt(Scale * Scale + rho * rho)),
        MapMode.TwoSided => 1.0 / (Factor * Math.Sqrt(Scale * Scale + rho * rho)),
        _ => throw Unsupported()
    };

    public static Exception Unsupported() => new Exception($"Unsupported mode '{Mode}' for RadialMap.");

}
