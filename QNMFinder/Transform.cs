// File: Transform.cs

using CLACFramework;
using System.Globalization;
using System.Numerics;
using System.Runtime;

namespace QNMFinder;

public enum TRMode
{
    Inert, Deradiate, Manual, Composite
}

public sealed class TRSpec
{
    public TRMode Mode { get; }
    public Complex Exp { get; }

    public TRSpec(TRMode mode, Complex exp)
    {  
        Mode = mode; 
        Exp = exp; 
    }
}

public enum BlendType
{
    OneSided, TwoSided
}

public static class Transformer
{
    public static double Scale { get; private set; }
    public static double Range { get; private set; }
    public static double Factor { get; private set; }

    public static TRSpec[] TRSpecs { get; private set; } = new TRSpec[2];

    public static void Configure(string trBlend, string[] trParams)
    {
        double[] trBlendParams = ParseBlend(trBlend);

        Scale = trBlendParams[0];
        Range = RadialMap.Range;
        Factor = Math.Sqrt(1 + (Scale * Scale) / (Range * Range));

        for (int i = 0; i < 2; i++)
            TRSpecs[i] = Parse(trParams[i]);
    }

    // ==============================================================
    // String Parsers
    // ==============================================================

    // Blend settings parser
    public static double[] ParseBlend(string spec)
    {
        if (string.IsNullOrWhiteSpace(spec))
            throw new Exception("Empty specification.");

        string s = spec.Trim().ToLowerInvariant();

        // Locate '(' and ')'
        int p1 = s.IndexOf('(');
        int p2 = s.LastIndexOf(')');

        if (p1 < 0 || p2 < 0 || p2 <= p1)
            throw new Exception($"Invalid format: '{spec}'");

        // Extract tag (ignored for now)
        string tag = s.Substring(0, p1).Trim();
        if (tag.Length == 0)
            throw new Exception($"Missing tag name in '{spec}'");

        if (tag != "global")
            throw new Exception($"Specification must start with 'Global': '{spec}'");

        // Extract inside parentheses
        string inner = s.Substring(p1 + 1, p2 - p1 - 1).Trim();
        if (inner.Length == 0)
            throw new Exception($"Tag '{tag}' must contain at least one numeric argument.");

        // Split by comma
        string[] parts = inner.Split(',');

        double[] values = new double[parts.Length];

        for (int i = 0; i < parts.Length; i++)
        {
            string token = parts[i].Trim();

            if (!double.TryParse(token, NumberStyles.Float, CultureInfo.InvariantCulture, out double val))
            {
                throw new Exception($"Invalid numeric value '{token}' in tag '{tag}'");
            }

            values[i] = val;
        }

        return values;
    }

    // Specification parser
    public static TRSpec Parse(string trSpec)
    {
        if (string.IsNullOrWhiteSpace(trSpec))
            return new TRSpec(TRMode.Inert, Complex.NaN);

        string lower = trSpec.Trim().ToLowerInvariant();

        if (lower == "inert")
            return new TRSpec(TRMode.Inert, Complex.NaN);

        if (lower == "deradiate")
            return new TRSpec(TRMode.Deradiate, Complex.NaN);

        bool isManual = lower.StartsWith("manual");
        bool isComposite = lower.StartsWith("composite");

        if (!isManual && !isComposite)
            throw new Exception($"Unsupported alterer specification: '{trSpec}'.");

        // Extract the parentheses content
        int p1 = lower.IndexOf('(');
        int p2 = lower.LastIndexOf(')');

        if (p1 < 0 || p2 < 0 || p2 <= p1)
            throw new Exception($"Invalid parameter specification: '{trSpec}'.");

        string inside = lower.Substring(p1 + 1, p2 - p1 - 1).Trim();
        string[] parts = inside.Split(',', StringSplitOptions.RemoveEmptyEntries);

        Complex exp = Complex.Zero;

        if (parts.Length == 1)
        {
            // Manual(c) or Composite(c)
            if (!double.TryParse(parts[0].Trim(), NumberStyles.Float, CultureInfo.InvariantCulture, out double c))
                throw new Exception($"Invalid real exponent: '{parts[0]}'.");

            exp = new Complex(c, 0.0);
        }
        else if (parts.Length == 2)
        {
            // Manual(a, b) or Composite(a, b)
            if (!double.TryParse(parts[0].Trim(), NumberStyles.Float, CultureInfo.InvariantCulture, out double a))
                throw new Exception($"Invalid real exponent: '{parts[0]}'.");

            if (!double.TryParse(parts[1].Trim(), NumberStyles.Float, CultureInfo.InvariantCulture, out double b))
                throw new Exception($"Invalid imaginary exponent: '{parts[1]}'.");

            exp = new Complex(a, b);
        }
        else throw new Exception($"Too many parameters in '{trSpec}'. Expected 1 or 2.");

        return new TRSpec(isManual ? TRMode.Manual : TRMode.Composite, exp);
    }

    // ==============================================================
    // Blend Functions
    // ==============================================================
    public static double Blend(double rho) => RadialMap.Mode switch
    {
        MapMode.OneSided => 
            rho * (Scale / Range) * (Scale / Range) * Math.Sqrt(rho * rho + Scale * Scale)
            / ((Scale * Scale * Scale * Scale) / ((2.0 * rho * rho + Scale * Scale) + 2.0 * rho * Math.Sqrt(rho * rho + Scale * Scale))
            + ((Scale / Range) * (Scale / Range) / (Factor + 1.0)) * (2.0 * rho * rho + Scale * Scale)),

        MapMode.TwoSided =>
            0.5 * (1 + Factor * rho / Math.Sqrt(rho * rho + Scale * Scale)),

        _ => throw RadialMap.Unsupported()
    };

    public static double DBlend(double rho) => RadialMap.Mode switch
    {
        MapMode.OneSided =>
            Factor * (Scale * Scale * Scale * Scale) * (Scale / Range) * (Scale / Range) / (Math.Sqrt(rho * rho + Scale * Scale)
            * ((Scale * Scale * Scale * Scale) / ((2.0 * rho * rho + Scale * Scale) + 2.0 * rho * Math.Sqrt(rho * rho + Scale * Scale))
            + ((Scale / Range) * (Scale / Range) / (Factor + 1.0)) * (2.0 * rho * rho + Scale * Scale))
            * ((Scale * Scale * Scale * Scale) / ((2.0 * rho * rho + Scale * Scale) + 2.0 * rho * Math.Sqrt(rho * rho + Scale * Scale))
            + ((Scale / Range) * (Scale / Range) / (Factor + 1.0)) * (2.0 * rho * rho + Scale * Scale))),

        MapMode.TwoSided =>
            0.5 * Factor * Scale * Scale / ((rho * rho + Scale * Scale) * Math.Sqrt(rho * rho + Scale * Scale)),

        _ => throw RadialMap.Unsupported()
    };

    public static double DDBlend(double rho) => RadialMap.Mode switch
    {
        MapMode.OneSided =>
            Factor * (Scale * Scale * Scale * Scale) * (Scale / Range) * (Scale / Range)
            * (((Scale * Scale * Scale * Scale) * (rho * rho + Scale * Scale) * (15.0 * rho * rho + 16.0 * Scale * Scale))
            / ((10.0 * rho * rho * rho * rho + 14.0 * Scale * Scale * rho * rho + 4.0 * Scale * Scale * Scale * Scale) 
            + rho * (10.0 * rho * rho + 9.0 * Scale * Scale) * Math.Sqrt(rho * rho + Scale * Scale)) 
            - ((Scale / Range) * (Scale / Range) / (Factor + 1.0)) * rho * (10.0 * rho * rho + 9.0 * Scale * Scale) * Math.Sqrt(rho * rho + Scale * Scale))
            / ((rho * rho + Scale * Scale) * (rho * rho + Scale * Scale) * ((Scale * Scale * Scale * Scale) 
            / ((2.0 * rho * rho + Scale * Scale) + 2.0 * rho * Math.Sqrt(rho * rho + Scale * Scale)) 
            + ((Scale / Range) * (Scale / Range) / (Factor + 1.0)) * (2.0 * rho * rho + Scale * Scale))
            * ((Scale * Scale * Scale * Scale) / ((2.0 * rho * rho + Scale * Scale) + 2.0 * rho * Math.Sqrt(rho * rho + Scale * Scale))
            + ((Scale / Range) * (Scale / Range) / (Factor + 1.0)) * (2.0 * rho * rho + Scale * Scale))
            * ((Scale * Scale * Scale * Scale) / ((2.0 * rho * rho + Scale * Scale) + 2.0 * rho * Math.Sqrt(rho * rho + Scale * Scale))
            + ((Scale / Range) * (Scale / Range) / (Factor + 1.0)) * (2.0 * rho * rho + Scale * Scale))),

        MapMode.TwoSided => -1.5 * Factor * rho * Scale * Scale / ((rho * rho + Scale * Scale) * (rho * rho + Scale * Scale) * Math.Sqrt(rho * rho + Scale * Scale)),

        _ => throw RadialMap.Unsupported()
    };

}

public enum TRType
{
    Inert, Constant, Variable, Mixed
}

public sealed class TRPrimit
{
    public TRType Type { get; }
    public Complex Lambda { get; }
    public Func<Complex, Complex> Gamma { get; }
    public Func<Complex, Complex> RGamma { get; }
    public TRPrimit(TRType type, Complex lambda, Func<Complex, Complex> gamma, Func<Complex, Complex> rGamma)
    {
        Type = type;
        Lambda = lambda;
        Gamma = gamma;
        RGamma = rGamma;
    }

    // ==============================================================
    // Transformation Primitives
    // ==============================================================

    public static TRPrimit[] ComputePrimitives(BCSpec[] bc, TRSpec[] tr)
    {
        TRPrimit[] primits = new TRPrimit[2];

        for (int side = 0; side < 2; side++)
        {
            BCSpec sideBC = bc[side];
            TRSpec sideTR = tr[side];

            Complex alpha = sideBC.Alpha;       // Robin α or Sommerfeld α=0
            SigmaSpec? sigma = sideBC.Sigma;    // Only Sommerfeld has SigmaSpec

            if (sideTR.Mode == TRMode.Deradiate)
            {
                if (sideBC.Type == BCType.Neumann)
                    primits[side] = new TRPrimit(TRType.Constant, Complex.Zero, null!, null!);
                else if (sideBC.Type == BCType.Robin)
                    primits[side] = new TRPrimit(TRType.Constant, Complex.ImaginaryOne * alpha.Imaginary, null!, null!);
                else if (sideBC.Type == BCType.Sommerfeld)
                    primits[side] = new TRPrimit(TRType.Variable, Complex.NaN, sigma!.Map, sigma!.RMap);
                else // Dirichlet
                    primits[side] = new TRPrimit(TRType.Constant, Complex.Zero, null!, null!);
            }
            else if (sideTR.Mode == TRMode.Manual)
            {
                // Pure, BC-independent constant factor: y = e^{-exp * rho} * u
                primits[side] = new TRPrimit(TRType.Constant, sideTR.Exp, null!, null!);
            }
            else if (sideTR.Mode == TRMode.Composite)
            {
                if (sideBC.Type == BCType.Robin)
                    primits[side] = new TRPrimit(TRType.Constant, Complex.ImaginaryOne * alpha.Imaginary + sideTR.Exp, null!, null!);
                else if (sideBC.Type == BCType.Sommerfeld)
                    primits[side] = new TRPrimit(TRType.Mixed, sideTR.Exp, sigma!.Map, sigma!.RMap);
                else // Dirichlet or Neumann
                    primits[side] = new TRPrimit(TRType.Constant, sideTR.Exp, null!, null!);
            }
            else // TRMode.Inert
            {
                primits[side] = new TRPrimit(TRType.Inert, Complex.NaN, null!, null!);
            }

        }

        return primits;
    }
}

public static class Transform
{
    // ============================================================
    // Transform System (EC and BC)
    // ============================================================

    public static (ACSpec, BCSpec[])
        TransformSystem(Func<double, Complex>[,] rawSystemEC, BCSpec[] rawSystemBC)
    {
        TRSpec[] trSpecs = Transformer.TRSpecs;

        (Func<double, Complex, Complex>[,] TC, Func<double, Complex, Complex>[,] RTC) = TransformEC(rawSystemBC, trSpecs);
        BCSpec[] systemBC = TransformBC(rawSystemBC, trSpecs);

        ACSpec systemAC = new ACSpec(rawSystemEC, TC, RTC);
        return (systemAC, systemBC);
    }

    // ==============================================================
    // Transform Interior Coefficients
    // ==============================================================

    private static (Func<double, Complex, Complex>[,] TC, Func<double, Complex, Complex>[,] RTC)
        TransformEC(BCSpec[] bc, TRSpec[] tr)
    {
        // Transformation Generating Jets
        Func<double, Complex, Complex>[,] TGJ = TransformationGeneratingJet(bc, tr);

        // Extract generating function jet entries
        Func<double, Complex, Complex> G = TGJ[0, 0];          // G
        Func<double, Complex, Complex> DG = TGJ[0, 1];         // G_ρ
        Func<double, Complex, Complex> DDG = TGJ[0, 2];        // G_ρρ

        Func<double, Complex, Complex> RG = TGJ[1, 0];         // G_ω
        Func<double, Complex, Complex> RDG = TGJ[1, 1];        // G_ωρ
        Func<double, Complex, Complex> RDDG = TGJ[1, 2];       // G_ωρρ

        // Composite quantities
        Func<double, Complex, Complex> H =
            (rho, omega) => G(rho, omega) + rho * DG(rho, omega);

        Func<double, Complex, Complex> RH =
            (rho, omega) => RG(rho, omega) + rho * RDG(rho, omega);

        Func<double, Complex, Complex> DH =
            (rho, omega) => 2.0 * DG(rho, omega) + rho * DDG(rho, omega);

        Func<double, Complex, Complex> RDH =
            (rho, omega) => 2.0 * RDG(rho, omega) + rho * RDDG(rho, omega);

        // TC matrix (transformation coefficients)
        Func<double, Complex, Complex>[,] TC = new Func<double, Complex, Complex>[3, 3];
        Func<double, Complex, Complex>[,] RTC = new Func<double, Complex, Complex>[3, 3];

        TC[0, 0] = (rho, omega) => 1.0;
        TC[0, 1] = (rho, omega) => 0.0;
        TC[0, 2] = (rho, omega) => 0.0;

        TC[1, 0] = (rho, omega) => -H(rho, omega);
        TC[1, 1] = (rho, omega) => 1.0;
        TC[1, 2] = (rho, omega) => 0.0;

        TC[2, 0] = (rho, omega) => H(rho, omega) * H(rho, omega) - DH(rho, omega);
        TC[2, 1] = (rho, omega) => -2.0 * H(rho, omega);
        TC[2, 2] = (rho, omega) => 1.0;

        // RTC matrix (omega derivative of TC)
        RTC[0, 0] = (rho, omega) => 0.0;
        RTC[0, 1] = (rho, omega) => 0.0;
        RTC[0, 2] = (rho, omega) => 0.0;

        RTC[1, 0] = (rho, omega) => -RH(rho, omega);
        RTC[1, 1] = (rho, omega) => 0.0;
        RTC[1, 2] = (rho, omega) => 0.0;

        RTC[2, 0] = (rho, omega) => 2.0 * H(rho, omega) * RH(rho, omega) - RDH(rho, omega);
        RTC[2, 1] = (rho, omega) => -2.0 * RH(rho, omega);
        RTC[2, 2] = (rho, omega) => 0.0;

        return (TC, RTC);
    }

    private static Func<double, Complex, Complex>[,]
        TransformationGeneratingJet(BCSpec[] bc, TRSpec[] tr)
    {
        TRPrimit[] prim = TRPrimit.ComputePrimitives(bc, tr);

        TRPrimit L = prim[0];
        TRPrimit R = prim[1];

        BlendType blendType = (L.Type == TRType.Inert || R.Type == TRType.Inert)
            ? BlendType.OneSided : BlendType.TwoSided;

        Func<double, Complex, Complex>[,] TGJ = new Func<double, Complex, Complex>[2, 3];

        // ------------------------------------------------------------
        // ONE-SIDED BLENDING
        // ------------------------------------------------------------

        if (blendType == BlendType.OneSided)
        {
            if (L.Type == TRType.Inert && R.Type == TRType.Inert)
            {
                TGJ[0, 0] = (rho, om) => Complex.Zero;
                TGJ[0, 1] = (rho, om) => Complex.Zero;
                TGJ[0, 2] = (rho, om) => Complex.Zero;

                TGJ[1, 0] = (rho, om) => Complex.Zero;
                TGJ[1, 1] = (rho, om) => Complex.Zero;
                TGJ[1, 2] = (rho, om) => Complex.Zero;
            }
            else
            {
                // Determine which side contributes
                TRPrimit P = (L.Type == TRType.Inert) ? R : L;

                TGJ[0, 0] = (rho, om) => GetValue(P, om);
                TGJ[0, 1] = (rho, om) => Complex.Zero;
                TGJ[0, 2] = (rho, om) => Complex.Zero;

                TGJ[1, 0] = (rho, om) => GetRValue(P, om);
                TGJ[1, 1] = (rho, om) => Complex.Zero;
                TGJ[1, 2] = (rho, om) => Complex.Zero;
            }

        }

        // ------------------------------------------------------------
        // TWO-SIDED BLENDING
        // ------------------------------------------------------------
        else // blendType == BlendType.TwoSided
        {
            TGJ[0, 0] = (rho, om) =>
            {
                double b = Transformer.Blend(rho);
                Complex GL = GetValue(L, om);
                Complex GR = GetValue(R, om);
                return (1.0 - b) * GL + b * GR;
            };

            TGJ[0, 1] = (rho, om) =>
            {
                double db = Transformer.DBlend(rho);
                Complex GL = GetValue(L, om);
                Complex GR = GetValue(R, om);
                return db * (GR - GL);
            };

            TGJ[0, 2] = (rho, om) =>
            {
                double ddb = Transformer.DDBlend(rho);
                Complex GL = GetValue(L, om);
                Complex GR = GetValue(R, om);
                return ddb * (GR - GL);
            };

            TGJ[1, 0] = (rho, om) =>
            {
                double b = Transformer.Blend(rho);
                Complex dGL = GetRValue(L, om);
                Complex dGR = GetRValue(R, om);
                return (1.0 - b) * dGL + b * dGR;
            };

            TGJ[1, 1] = (rho, om) =>
            {
                double db = Transformer.DBlend(rho);
                Complex dGL = GetRValue(L, om);
                Complex dGR = GetRValue(R, om);
                return db * (dGR - dGL);
            };

            TGJ[1, 2] = (rho, om) =>
            {
                double ddb = Transformer.DDBlend(rho);
                Complex dGL = GetRValue(L, om);
                Complex dGR = GetRValue(R, om);
                return ddb * (dGR - dGL);
            };

        }

        return TGJ;
    }

    // ==============================================================
    // Transform Boundary Conditions
    // ==============================================================
    private static BCSpec[] TransformBC(BCSpec[] bc, TRSpec[] tr)
    {
        TRPrimit[] prim = TRPrimit.ComputePrimitives(bc, tr);

        TRPrimit L = prim[0];
        TRPrimit R = prim[1];

        BCSpec[] tbc = new BCSpec[2];

        // Left side: "self" = L, "other" = R
        tbc[0] = TransformSide(bc[0], L, R);

        // Right side: "self" = R, "other" = L
        tbc[1] = TransformSide(bc[1], R, L);

        return tbc;
    }

    private static BCSpec TransformSide(BCSpec bc, TRPrimit self, TRPrimit other)
    {
        if (self.Type == TRType.Inert && other.Type == TRType.Inert)
            return bc;

        else if (self.Type == TRType.Inert && other.Type != TRType.Inert)
        {
            if (bc.Type == BCType.Dirichlet)
                return BCSpec.Dirichlet();
            else if (bc.Type == BCType.Neumann)
            {
                if (other.Type == TRType.Constant)
                    return new BCSpec(BCType.Robin, -other.Lambda, null);
                else // other.Type == TRType.Variable || other.Type == TRType.Mixed
                {
                    SigmaSpec sigma = new SigmaSpec(SigmaType.Custom, om => -GetValue(other, om), om => -GetRValue(other, om));
                    return new BCSpec(BCType.Sommerfeld, Complex.Zero, sigma);
                }
            }
            else if (bc.Type == BCType.Robin)
            {
                if (other.Type == TRType.Constant)
                    return new BCSpec(BCType.Robin, bc.Alpha - other.Lambda, null);
                else // other.Type == TRType.Variable || other.Type == TRType.Mixed
                {
                    SigmaSpec sigma = new SigmaSpec(SigmaType.Custom, om => bc.Alpha - GetValue(other, om), om => -GetRValue(other, om));
                    return new BCSpec(BCType.Sommerfeld, Complex.Zero, sigma);
                }
            }
            else // bc.Type == BCType.Sommerfeld
            {
                SigmaSpec sigma = new SigmaSpec(SigmaType.Custom, om => bc.Sigma!.Map(om) - GetValue(other, om), om => bc.Sigma!.RMap(om) - GetRValue(other, om));
                return new BCSpec(BCType.Sommerfeld, Complex.Zero, sigma);
            }
        }
        else if (self.Type == TRType.Constant)
        {
            if (bc.Type == BCType.Dirichlet)
                return BCSpec.Dirichlet();
            else if (bc.Type == BCType.Neumann)
                return new BCSpec(BCType.Robin, -self.Lambda, null);
            else if (bc.Type == BCType.Robin)
                return new BCSpec(BCType.Robin, bc.Alpha - self.Lambda, null);
            else // bc.Type == BCType.Sommerfeld
            {
                SigmaSpec sigma = new SigmaSpec(SigmaType.Custom, om => bc.Sigma!.Map(om) - GetValue(self, om), om => bc.Sigma!.RMap(om) - GetRValue(self, om));
                return new BCSpec(BCType.Sommerfeld, Complex.Zero, sigma);
            }

        }
        else if (self.Type == TRType.Variable) // Only Sommerfeld produces this
            return BCSpec.Neumann();
        else // self.Type == TRType.Mixed // Only Sommerfeld produces this
            return new BCSpec(BCType.Robin, -self.Lambda, null);
    }

    // ==============================================================
    // Value Extractor
    // ==============================================================
    static Complex GetValue(TRPrimit p, Complex om)
    {
        return p.Type switch
        {
            TRType.Inert => Complex.Zero,
            TRType.Constant => p.Lambda,
            TRType.Variable => p.Gamma(om),
            TRType.Mixed => p.Lambda + p.Gamma(om),
            _ => Complex.Zero
        };
    }

    static Complex GetRValue(TRPrimit p, Complex om)
    {
        return p.Type switch
        {
            TRType.Inert => Complex.Zero,
            TRType.Constant => Complex.Zero,
            TRType.Variable => p.RGamma(om),
            TRType.Mixed => p.RGamma(om),
            _ => Complex.Zero
        };
    }

}
