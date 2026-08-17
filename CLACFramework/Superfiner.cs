// File: Superfiner.cs

using System;
using System.Numerics;
using MathNet.Numerics;
using IOFramework;
using CLACFramework;

namespace QNMFinder;

public sealed class Superfiner
{
    private readonly Kernel _kernel;

    public Superfiner(Kernel kernel)
    {
        _kernel = kernel;
    }

    // Superfine a list of refined eigenpairs
    public (List<EigenPair>, List<EstimatedLimit>) SuperfineQNM(List<EigenPair> feeds, int I, int J)
    {
        List<EigenPair> superfined = new List<EigenPair>();
        List<EstimatedLimit> polished = new List<EstimatedLimit>();

        for (int k = 0; k < feeds.Count; k++)
        {
            EigenPair feed = feeds[k];
            EigenPair superfinedEP = feed;
            (Complex omega, VectorC y, double residual) = superfinedEP;

            Logger.WriteBoth($"\nMode {k:D2}: Superfining frequency {omega.Format("F16")}...");

            int IP = I;
            int jP = 0;

            Complex[] omegas = new Complex[J];

            for (int j = 0; j < J; j++)
            {
                IP = Math.Max(8, IP / 2);
                jP = j + 1;

                y = Discretizer.Interpolate(y, 2);
                int NP = y.Size;

                Logger.WriteBoth($"\nPass {j:D2}: Superfining {omega.Format("F16")} on {NP} nodes...");
                Logger.WriteLine("");

                feed = new EigenPair(omega, y, residual);
                superfinedEP = SuperfineSingleQNM(feed, IP, jP);

                (omega, y, residual) = superfinedEP;
                omegas[j] = omega;
            }

            Complex omegaPolished = PolishFrequency(omegas);
            double error = Complex.Abs(omegaPolished - omega);

            Logger.WriteBoth($"\nPass SE: Freq = {omegaPolished.Format("F16")}, Error = {error.Format("E6")}");

            superfined.Add(superfinedEP);
            polished.Add(new EstimatedLimit(omegaPolished, error));
        }

        return (superfined, polished);
    }

    public EigenPair SuperfineSingleQNM(EigenPair feed, int I, int K)
    {
        (Complex omega, VectorC y, double residual) = feed;

        int NP = y.Size;
        int N = NP - 1;

        // Assemble system operators and functions
        (VectorC x, ActionC D) = Discretizer.DiscretizeA(N);
        (Func<Complex, ActionC> T, Func<Complex, ActionC> RT) = _kernel.SystemOperatorActions(x, D);

        ActionC TA = T(omega);
        ActionC RTA = RT(omega);

        VectorC r = VectorC.Zero(NP);
        double OpNorm = 0.0;

        for (int i = 0; i < I; i++)
        {
            (Complex domega, VectorC dy) = NewtonUpdate(TA, RTA, y, K);

            omega += domega;
            y += dy;

            y.Normalize();

            // Update operators
            TA = T(omega);
            RTA = RT(omega);

            // Scalar residual
            r = TA.Apply(y);
            OpNorm = TA.FrobeniusNorm();
            residual = r.Norm() / OpNorm;

            Logger.WriteLine($"Iter {i:D2}: Freq = {omega.Format("F16")}, Res = {residual.Format("E6")}");

        }

        EigenPair superfinedEP = new EigenPair(omega, y, residual);
        return superfinedEP;
    }

    public static (Complex, VectorC) NewtonUpdate(ActionC TA, ActionC RTA, VectorC y, int K)
    {
        int NP = TA.Rows;
        int NPP = NP + 1;

        // Compute r = T(ω) y and s = T'(ω) y
        VectorC r = new VectorC(NP);
        VectorC s = new VectorC(NP);

        TA.Apply(y, r);
        RTA.Apply(y, s);

        // Build the bordered operator TG = [T  RT * y | y^H 0]

        ActionC TAG = new ActionC(NPP, NPP,
            (VectorC wG, VectorC TAGwG) =>
            {
                // Scratch buffers
                VectorC w = new VectorC(NP);
                VectorC TMw = new VectorC(NP);

                // Extract w = wG[0:NP], σ = wG[NP]
                w.CopyFrom(wG, 0, NP);
                Complex sigma = wG[NP];

                // top = T(ω)w + s * σ
                TA.Apply(w, TMw);
                TMw.AxPyIP(sigma, s);

                // Write top block
                TAGwG.SetSubvector(0, TMw);

                // Bottom block: y^H v
                TAGwG[NP] = y.Dot(w);
            }
        );

        // Build rG = (r, 0)
        VectorC rG = new VectorC(NPP);
        rG.SetSubvector(0, r);
        rG[NP] = Complex.Zero;

        // Solve TAG(dyG) = - rG
        VectorC dyG = - TAG.Solve(rG, K);

        // Extract dy and domega
        VectorC dy = dyG.GetSubvector(0, NP);
        Complex domega = dyG[NP];

        return (domega, dy);
    }

    public static Complex PolishFrequency(Complex[] omegas)
    {
        int n = omegas.Length;
        if (n == 0) throw new ArgumentException("Empty sequence");
        if (n == 1) return omegas[0];

        const double tiny = 1e-300;

        Complex[] eps_km1 = (Complex[])omegas.Clone();
        Complex[] eps_km2 = new Complex[n]; // ε_-1 = 0

        Complex lastEven = eps_km1[n - 1];

        for (int k = 1; k < n; k++)
        {
            int m = n - k;
            Complex[] eps_k = new Complex[m];

            for (int j = 0; j < m; j++)
            {
                Complex denom = eps_km1[j + 1] - eps_km1[j];

                if (denom.Magnitude < tiny)
                {
                    eps_k[j] = eps_km2[j + 1];
                }
                else
                {
                    eps_k[j] = eps_km2[j + 1] + Complex.One / denom;
                }
            }

            if ((k & 1) == 0)
                lastEven = eps_k[0];

            eps_km2 = eps_km1;
            eps_km1 = eps_k;
        }

        return lastEven;
    }

}
