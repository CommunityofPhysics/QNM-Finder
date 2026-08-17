// File: Kernel.cs

using CLACFramework;
using IOFramework;
using System.Numerics;

namespace QNMFinder;

public sealed class Kernel
{
    // ============================================================
    // 1. Raw System Specifications
    // ============================================================

    private readonly Func<double, Complex>[,] _RawSystemEC;
    public Func<double, Complex>[,] RawSystemEC => _RawSystemEC;

    private readonly BCSpec[] _RawSystemBC;
    public BCSpec[] RawSystemBC => _RawSystemBC;

    // ============================================================
    // 2. Transformed System Specifications
    // ============================================================

    private readonly ACSpec _SystemAC;
    public ACSpec SystemAC => _SystemAC;

    private readonly BCSpec[] _SystemBC;
    public BCSpec[] SystemBC => _SystemBC;

    private readonly Complex[,] _SystemBP;
    public Complex[,] SystemBP => _SystemBP;

    // ============================================================
    // 4. System Operators (Matrices and Actions)
    // ============================================================

    private readonly Func<VectorC, MatrixC, (Func<Complex, MatrixC> T, Func<Complex, MatrixC> RT)> _SystemOM;
    public Func<VectorC, MatrixC, (Func<Complex, MatrixC> T, Func<Complex, MatrixC> RT)> SystemOperatorMatrices => _SystemOM;

    private readonly Func<VectorC, ActionC, (Func<Complex, ActionC> T, Func<Complex, ActionC> RT)> _SystemOA;
    public Func<VectorC, ActionC, (Func<Complex, ActionC> T, Func<Complex, ActionC> RT)> SystemOperatorActions => _SystemOA;

    // ============================================================
    // 5. Constructor
    // ============================================================

    public Kernel(Func<double, Complex>[,] rawSystemEC, BCSpec[] rawSystemBC)
    {
        // Store raw system specifications
        _RawSystemEC = rawSystemEC;
        _RawSystemBC = rawSystemBC;

        (ACSpec systemAC, BCSpec[] systemBC) 
            = Transform.TransformSystem(rawSystemEC, rawSystemBC);

        // Store transformed system specifications
        _SystemAC = systemAC;
        _SystemBC = systemBC;

        // Compute branch points for the system
        Complex[,] systemBP = Edge.BranchPoints(rawSystemEC, rawSystemBC);

        _SystemBP = systemBP;

        // Store system operators
        _SystemOM = (x, D) => SystemOperators(x, D, systemAC, systemBC);
        _SystemOA = (x, D) => SystemOperators(x, D, systemAC, systemBC);
    }

    // ============================================================
    // 6. Configurer
    // ============================================================
    public static Kernel Construct(Config config)
    {
        RadialMap.Configure(config.MapSpec);

        Blender.Configure(config.TRBlend);

        Transformer.Configure(config.TRParams);

        Func<double, Complex>[,] rawSystemEC
            = Initiator.CompileSystemFunctions(config.SystemEC, config.SystemFP, config.Physics);

        Func<Complex, Complex>[][] rawCustomSigma
            = Initiator.CompileSigma(config.SigmaCM, config.SigmaFP, config.Physics);

        (bool[] isBadSigma, Func<Complex, Complex>[][] rawIncomingSigma, Func<Complex, Complex>[][] rawOutgoingSigma)
            = Edge.SigmaMaps(rawSystemEC);

        string[] rawBCStr = config.SystemBC;

        SigmaLib[] rawSigmaLibs = new SigmaLib[2];
        BCSpec[] rawSystemBC = new BCSpec[2];

        for (int side = 0; side < 2; side++)
        {
            rawSigmaLibs[side] = new SigmaLib(isBadSigma[side], rawCustomSigma[side], rawIncomingSigma[side], rawOutgoingSigma[side]);
            rawSystemBC[side] = BCSpec.Parse(rawBCStr[side], rawSigmaLibs[side]);
        }

        return new Kernel(rawSystemEC, rawSystemBC);
    }

    // ============================================================
    // 7. System Operator Builders
    // ============================================================

    private static (Func<Complex, MatrixC> T, Func<Complex, MatrixC> RT) 
        SystemOperators(VectorC x, MatrixC D, ACSpec systemAC, BCSpec[] systemBC)
    {
        (Func<double, Complex, Complex>[] F, Func<double, Complex, Complex>[] RF) = OperatorCoefficients(systemAC);

        Func<Complex, MatrixC> P = Discretizer.MatrixF(x, F[0]);

        Func<Complex, MatrixC> Q1 = Discretizer.MatrixF(x, F[1]);
        Func<Complex, MatrixC> Q2 = Discretizer.MatrixF(x, F[2]);
        Func<Complex, MatrixC> Q = omega => Q1(omega) + Q2(omega)* D;

        Func<Complex, MatrixC> RP = Discretizer.MatrixF(x, RF[0]);

        Func<Complex, MatrixC> RQ1 = Discretizer.MatrixF(x, RF[1]);
        Func<Complex, MatrixC> RQ2 = Discretizer.MatrixF(x, RF[2]);
        Func<Complex, MatrixC> RQ = omega => RQ1(omega) + RQ2(omega)* D;

        (Func<Complex, MatrixC> R, Func<Complex, MatrixC> S, Func<Complex, MatrixC> RR, Func<Complex, MatrixC> RS) 
            = Boundary.ImposeBC(P, Q, RP, RQ, systemBC);

        Func<Complex, MatrixC> T = omega => R(omega) + S(omega) * D;
        Func<Complex, MatrixC> RT = omega => RR(omega) + RS(omega) * D;

        return (T, RT);
    }

    public static (Func<Complex, ActionC> T, Func<Complex, ActionC> RT) 
        SystemOperators(VectorC x, ActionC D, ACSpec systemAC, BCSpec[] systemBC)
    {
        (Func<double, Complex, Complex>[] F, Func<double, Complex, Complex>[] RF) = OperatorCoefficients(systemAC);

        Func<Complex, ActionC> P = Discretizer.ActionF(x, F[0]);

        Func<Complex, ActionC> Q1 = Discretizer.ActionF(x, F[1]);
        Func<Complex, ActionC> Q2 = Discretizer.ActionF(x, F[2]);
        Func<Complex, ActionC> Q = omega => Q1(omega) + Q2(omega) * D;

        Func<Complex, ActionC> RP = Discretizer.ActionF(x, RF[0]);

        Func<Complex, ActionC> RQ1 = Discretizer.ActionF(x, RF[1]);
        Func<Complex, ActionC> RQ2 = Discretizer.ActionF(x, RF[2]);
        Func<Complex, ActionC> RQ = omega => RQ1(omega) + RQ2(omega) * D;

        (Func<Complex, ActionC> R, Func<Complex, ActionC> S, Func<Complex, ActionC> RR, Func<Complex, ActionC> RS)
            = Boundary.ImposeBC(P, Q, RP, RQ, systemBC);

        Func<Complex, ActionC> T = omega => R(omega) + S(omega) * D;
        Func<Complex, ActionC> RT = omega => RR(omega) + RS(omega) * D;

        return (T, RT);
    }

    private static (Func<double, Complex, Complex>[] F, Func<double, Complex, Complex>[] RF) 
        OperatorCoefficients(ACSpec systemAC)
    {
        Func<double, Complex>[,] ec = systemAC.EC;
        Func<double, Complex, Complex>[,] TC = systemAC.TC;
        Func<double, Complex, Complex>[,] RTC = systemAC.RTC;

        Func<double, Complex, Complex>[] F = new Func<double, Complex, Complex>[3];
        Func<double, Complex, Complex>[] RF = new Func<double, Complex, Complex>[3];

        // ------------------------------------------------------------
        // Computing F: F0 * y + F1 Dy + F2 DDy = 0
        // ------------------------------------------------------------

        for (int l = 0; l < 3; l++)
        {
            int k = l;

            // F_k = omega^i E_i^j T_j^k
            F[k] = (rho, omega) =>
                // i = 0
                ec[0, 0](rho) * TC[0, k](rho, omega)
              + ec[0, 1](rho) * TC[1, k](rho, omega)
              + ec[0, 2](rho) * TC[2, k](rho, omega)
              // i = 1
              + omega * (ec[1, 0](rho) * TC[0, k](rho, omega)
                  + ec[1, 1](rho) * TC[1, k](rho, omega))
              // i = 2
              + omega * omega * (ec[2, 0](rho) * TC[0, k](rho, omega));

            // RF_k = [ i omega^{i-1} E_i^j T_j^k + omega^i E_i^j RT_j^k ]
            RF[k] = (rho, omega) =>
                // i omega^{i-1} E_i^j T_j^k
                // i = 1
                (ec[1, 0](rho) * TC[0, k](rho, omega)
               + ec[1, 1](rho) * TC[1, k](rho, omega))
              // i = 2
              + 2.0 * omega * (ec[2, 0](rho) * TC[0, k](rho, omega))
              // omega^i E_i^j RT_j^k 
              // i = 0
              + (ec[0, 0](rho) * RTC[0, k](rho, omega)
                  + ec[0, 1](rho) * RTC[1, k](rho, omega)
                  + ec[0, 2](rho) * RTC[2, k](rho, omega))
              // i = 1
              + omega * (ec[1, 0](rho) * RTC[0, k](rho, omega)
                  + ec[1, 1](rho) * RTC[1, k](rho, omega))
              // i = 2
              + omega * omega * (ec[2, 0](rho) * RTC[0, k](rho, omega));
        }

        return (F, RF);
    }

}

public sealed class ACSpec
{
    // Equation coefficients EC[i,j](rho)
    // Transformation coefficients TC[j,k](rho, omega)
    // Omega derivative of transformation coefficients RTC[j,k](rho, omega)
    public Func<double, Complex>[,] EC { get; }
    public Func<double, Complex, Complex>[,] TC { get; }
    public Func<double, Complex, Complex>[,] RTC { get; }

    public ACSpec(Func<double, Complex>[,] ec, Func<double, Complex, Complex>[,] tc, Func<double, Complex, Complex>[,] rtc)
    {
        // Validate inputs are not null
        EC = ec ?? throw new ArgumentNullException(nameof(ec));
        TC = tc ?? throw new ArgumentNullException(nameof(tc));
        RTC = rtc ?? throw new ArgumentNullException(nameof(rtc));

        // Validate EC shape
        int ecRows = EC.GetLength(0);
        int ecCols = EC.GetLength(1);

        // Validate TC/RTC shapes match EC
        int tcRows = TC.GetLength(0);
        int tcCols = TC.GetLength(1);
        int rtcRows = RTC.GetLength(0);
        int rtcCols = RTC.GetLength(1);

        if (tcRows != ecRows || tcCols != ecCols)
            throw new ArgumentException("TC must have the same dimensions as EC.", nameof(tc));

        if (rtcRows != ecRows || rtcCols != ecCols)
            throw new ArgumentException("RTC must have the same dimensions as EC.", nameof(rtc));
    }

    public ACSpec(Func<double, Complex>[,] ec)
    {
        // Validate inputs are not null
        EC = ec ?? throw new ArgumentNullException(nameof(ec));

        // Validate EC shape
        int ecRows = EC.GetLength(0);
        int ecCols = EC.GetLength(1);

        (TC, RTC) = IdentityTransform(ecRows, ecCols);

    }

    public static (Func<double, Complex, Complex>[,] TC, Func<double, Complex, Complex>[,] RTC)
        IdentityTransform(int rows, int cols)
    {
        var tc = new Func<double, Complex, Complex>[rows, cols];
        var rtc = new Func<double, Complex, Complex>[rows, cols];

        for (int j = 0; j < rows; j++)
            for (int k = 0; k < cols; k++)
            {
                if (j == k)
                {
                    tc[j, k] = (rho, omega) => 1.0;
                }
                else
                {
                    tc[j, k] = (rho, omega) => 0.0;
                }

                rtc[j, k] = (rho, omega) => 0.0;
            }

        return (tc, rtc);
    }
}
