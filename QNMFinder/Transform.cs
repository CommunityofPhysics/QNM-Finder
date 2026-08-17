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

public enum TRType
{
    Inert, Constant, Variable, Mixed
}

public enum BlendType
{
    OneSided, TwoSided
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

// ==============================================================
// Transformation System
// ==============================================================

public static class Transformer
{
    public static TRSpec[] TRSpecs { get; private set; } = new TRSpec[2];

    public static void Configure(string[] trParams)
    {
        for (int i = 0; i < 2; i++)
            TRSpecs[i] = Parse(trParams[i]);
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
                double b = Blender.Blend(rho);
                Complex GL = GetValue(L, om);
                Complex GR = GetValue(R, om);
                return (1.0 - b) * GL + b * GR;
            };

            TGJ[0, 1] = (rho, om) =>
            {
                double db = Blender.DBlend(rho);
                Complex GL = GetValue(L, om);
                Complex GR = GetValue(R, om);
                return db * (GR - GL);
            };

            TGJ[0, 2] = (rho, om) =>
            {
                double ddb = Blender.DDBlend(rho);
                Complex GL = GetValue(L, om);
                Complex GR = GetValue(R, om);
                return ddb * (GR - GL);
            };

            TGJ[1, 0] = (rho, om) =>
            {
                double b = Blender.Blend(rho);
                Complex dGL = GetRValue(L, om);
                Complex dGR = GetRValue(R, om);
                return (1.0 - b) * dGL + b * dGR;
            };

            TGJ[1, 1] = (rho, om) =>
            {
                double db = Blender.DBlend(rho);
                Complex dGL = GetRValue(L, om);
                Complex dGR = GetRValue(R, om);
                return db * (dGR - dGL);
            };

            TGJ[1, 2] = (rho, om) =>
            {
                double ddb = Blender.DDBlend(rho);
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

// ==============================================================
// Blender System
// ==============================================================
public static class Blender
{
    private static double Stiff { get; set; }
    private static double Width { get; set; }
    private static double Scale { get; set; }
    private static double Range { get; set; }
    private static double Factor { get; set; }
    private static double Boost { get; set; }

    public static void Configure(string trBlend)
    {
        (Stiff, Width, Scale, Range) = Parse(trBlend);

        Factor = Math.Asinh(Range / Scale);
        Boost = 2.0 * Math.Acosh(0.5 * Factor / Math.Asinh(Width / Scale));
    }

    // Blend Specification Parser
    public static (double, double, double, double) Parse(string spec)
    {
        const double Tiny = 1e-10;

        double width = RadialMap.Width;
        double scale = RadialMap.Scale;
        double range = RadialMap.Range;

        if (string.IsNullOrWhiteSpace(spec)) throw new Exception("Blender specification is empty.");

        spec = spec.Trim();

        int p1 = spec.IndexOf('(');
        int p2 = spec.LastIndexOf(')');

        if (p1 < 0 || p2 < 0 || p2 <= p1)
            throw new Exception($"Invalid Blender specification: {spec}");

        string nameStr = new string(spec.Substring(0, p1).Trim().ToLowerInvariant()
            .Where(c => !char.IsWhiteSpace(c)).ToArray());
        if (nameStr != "blend") throw new Exception($"Blender name '{nameStr}' is unsupported. It must start with 'Blend'");

        string parameters = spec.Substring(p1 + 1, p2 - p1 - 1).Trim();

        if (string.IsNullOrWhiteSpace(parameters)) throw new Exception($"Blender requires at least a stiff parameter.");

        string[] parts = parameters.Split(',', StringSplitOptions.RemoveEmptyEntries);

        if (parts.Length < 1 || parts.Length > 3) throw new Exception($"Blender accepts parameters as (Stiff), (Stiff, Width), or (Stiff, Width, Scale).");

        if (!double.TryParse(parts[0].Trim(), NumberStyles.Float, CultureInfo.InvariantCulture, out double stiff))
            throw new Exception($"Invalid stiff value: '{parts[0]}'");

        stiff = Math.Max(Math.Abs(stiff), Tiny);

        if (parts.Length > 1)
        {
            string widthStr = parts[1].Trim().ToLowerInvariant();

            if (double.TryParse(widthStr, NumberStyles.Float, CultureInfo.InvariantCulture, out double parsedWidth))
                width = Math.Max(Math.Abs(parsedWidth), Tiny);

            else throw new Exception($"Invalid width value: '{widthStr}'");
        }

        if (parts.Length > 2)
        {
            string scaleStr = parts[2].Trim().ToLowerInvariant();

            if (double.TryParse(scaleStr, NumberStyles.Float, CultureInfo.InvariantCulture, out double parsedScale))
                scale = Math.Max(Math.Abs(parsedScale), Tiny);

            else throw new Exception($"Invalid scale value: '{scaleStr}'");
        }

        scale = Math.Min(scale, range / 2.0);
        width = Math.Min(width, Math.Sqrt(scale * range) / 2.0);

        return (stiff, width, scale, range);
    }

    // Blend Function
    public static double Blend(double rho) => RadialMap.Mode switch
    {
        MapMode.OneSided => 0.5 * (1.0 + Math.Tanh(Stiff * ((2.0 / Boost) * Math.Asinh((Math.Asinh(rho / Scale) / Factor) * Math.Sinh(Boost)) - 1.0)) / Math.Tanh(Stiff)),
        MapMode.TwoSided => 0.5 * (1.0 + Math.Tanh(Stiff * ((1.0 / Boost) * Math.Asinh((Math.Asinh(rho / Scale) / Factor) * Math.Sinh(Boost)))) / Math.Tanh(Stiff)),

        _ => throw RadialMap.Unsupported()
    };

    // Derivative of Blend Function
    public static double DBlend(double rho) => RadialMap.Mode switch
    {
        MapMode.OneSided => Stiff * Math.Sinh(Boost) / (Boost * Factor * Math.Sqrt(1.0 + (Math.Sinh(Boost) * Math.Sinh(Boost))
            * (Math.Asinh(rho / Scale) * Math.Asinh(rho / Scale)) / (Factor * Factor)) * Math.Sqrt(Scale * Scale + rho * rho)
            * Math.Pow(Math.Cosh(Stiff * (-1.0 + 2.0 * Math.Asinh(Math.Sinh(Boost) * Math.Asinh(rho / Scale) / Factor) / Boost)), 2.0) * Math.Tanh(Stiff)),

        MapMode.TwoSided => 0.5 * Stiff * Math.Sinh(Boost) / (Boost * Factor * Math.Sqrt(1.0 + (Math.Sinh(Boost) * Math.Sinh(Boost))
            * (Math.Asinh(rho / Scale) * Math.Asinh(rho / Scale)) / (Factor * Factor)) * Math.Sqrt(Scale * Scale + rho * rho)
            * Math.Pow(Math.Cosh(Stiff * Math.Asinh(Math.Sinh(Boost) * Math.Asinh(rho / Scale) / Factor) / Boost), 2.0) * Math.Tanh(Stiff)),

        _ => throw RadialMap.Unsupported()
    };

    // Second Derivative of Blend Function
    public static double DDBlend(double rho) => RadialMap.Mode switch
    {
        MapMode.OneSided => 0.5 * Stiff * (-8.0 * Stiff * (Math.Sinh(Boost) * Math.Sinh(Boost)) * Math.Tanh(Stiff * (-1.0 + 2.0 * Math.Asinh(Math.Sinh(Boost)
            * Math.Asinh(rho / Scale) / Factor) / Boost)) / (Boost * Boost * Factor * Factor * (1.0 + (Math.Sinh(Boost) * Math.Sinh(Boost))
            * (Math.Asinh(rho / Scale) * Math.Asinh(rho / Scale)) / (Factor * Factor)) * (Scale * Scale + rho * rho)) - 2.0 * rho * Math.Sinh(Boost)
            / (Boost * Factor * Math.Sqrt(1.0 + (Math.Sinh(Boost) * Math.Sinh(Boost)) * (Math.Asinh(rho / Scale) * Math.Asinh(rho / Scale)) / (Factor * Factor))
            * (Scale * Scale + rho * rho) * Math.Sqrt(Scale * Scale + rho * rho)) - 2.0 * Math.Pow(Math.Sinh(Boost), 3.0) * Math.Asinh(rho / Scale)
            / (Boost * Math.Pow(Factor, 3.0) * Scale * Math.Sqrt(1.0 + (rho * rho) / (Scale * Scale)) * (1.0 + (Math.Sinh(Boost) * Math.Sinh(Boost))
            * (Math.Asinh(rho / Scale) * Math.Asinh(rho / Scale)) / (Factor * Factor)) * Math.Sqrt(1.0 + (Math.Sinh(Boost) * Math.Sinh(Boost))
            * (Math.Asinh(rho / Scale) * Math.Asinh(rho / Scale)) / (Factor * Factor)) * Math.Sqrt(Scale * Scale + rho * rho)))
            / (Math.Pow(Math.Cosh(Stiff * (-1.0 + 2.0 * Math.Asinh(Math.Sinh(Boost) * Math.Asinh(rho / Scale) / Factor) / Boost)), 2.0) * Math.Tanh(Stiff)),

        MapMode.TwoSided => 0.5 * Stiff * (-2.0 * Stiff * (Math.Sinh(Boost) * Math.Sinh(Boost)) * Math.Tanh(Stiff * Math.Asinh(Math.Sinh(Boost)
            * Math.Asinh(rho / Scale) / Factor) / Boost) / (Boost * Boost * Factor * Factor * (1.0 + (Math.Sinh(Boost) * Math.Sinh(Boost))
            * (Math.Asinh(rho / Scale) * Math.Asinh(rho / Scale)) / (Factor * Factor)) * (Scale * Scale + rho * rho)) - rho * Math.Sinh(Boost)
            / (Boost * Factor * Math.Sqrt(1.0 + (Math.Sinh(Boost) * Math.Sinh(Boost)) * (Math.Asinh(rho / Scale) * Math.Asinh(rho / Scale))
            / (Factor * Factor)) * (Scale * Scale + rho * rho) * Math.Sqrt(Scale * Scale + rho * rho)) - Math.Pow(Math.Sinh(Boost), 3.0) * Math.Asinh(rho / Scale)
            / (Boost * Math.Pow(Factor, 3.0) * Scale * Math.Sqrt(1.0 + (rho * rho) / (Scale * Scale)) * (1.0 + (Math.Sinh(Boost) * Math.Sinh(Boost))
            * (Math.Asinh(rho / Scale) * Math.Asinh(rho / Scale)) / (Factor * Factor)) * Math.Sqrt(1.0 + (Math.Sinh(Boost) * Math.Sinh(Boost))
            * (Math.Asinh(rho / Scale) * Math.Asinh(rho / Scale)) / (Factor * Factor)) * Math.Sqrt(Scale * Scale + rho * rho)))
            / (Math.Pow(Math.Cosh(Stiff * Math.Asinh(Math.Sinh(Boost) * Math.Asinh(rho / Scale) / Factor) / Boost), 2.0) * Math.Tanh(Stiff)),

        _ => throw RadialMap.Unsupported()
    };

}
