// File: Estimator.cs

using System;
using System.Numerics;
using CLACFramework;
using IOFramework;

namespace QNMFinder;

public sealed class Estimator
{
    private readonly Kernel _kernel;

    private (double processed, double total) _area;

    public Estimator(Kernel kernel)
    {
        _kernel = kernel;
    }

    // ===================================================
    // 1. Main QNM locator
    // ===================================================

    public (List<EigenPair>, List<Complex>) EstimateQNM(MatrixC beynPencil, SqDomainC[] domains, int M, int V)
    {
        int NP = beynPencil.Rows;
        int N = NP - 1;

        (VectorC x, MatrixC D) = Discretizer.DiscretizeM(N);
        (Func<Complex, MatrixC> T, Func<Complex, MatrixC> RT) = _kernel.SystemOperatorMatrices(x, D);

        foreach (SqDomainC domain in domains)
        {
            _area.total += domain.Area;
        }

        Depth depth = new Depth(0, V);

        List<EigenPair> rootEstimates = new List<EigenPair>();
        List<Complex> rootSingularities = new List<Complex>();

        foreach (SqDomainC domain in domains)
        {
            if (DoesContainBP(domain))
            {
                Logger.WriteLine($"\nDetected branch point inside the contour. Subdividing {domain.Format("F4")}.");
                Logger.WriteLine("");

                Depth nextDepth = depth.Next();
                SqDomainC[] subdomains = domain.Subdivide(3);

                foreach (SqDomainC subdomain in subdomains)
                {
                    if (DoesContainBP(subdomain))
                    {
                        Logger.WriteLine($"\nDetected branch point inside the contour. Skipping {subdomain.Format("F4")}.");
                        Logger.WriteLine("");

                        _area.processed += subdomain.Area;
                        rootSingularities.Add(subdomain.Center);
                        continue;
                    }
                    else
                    {
                        (List<EigenPair> subEstimates, List<Complex> subSingularities) = RecursiveEstimateQNM(T, RT, beynPencil, subdomain, M, nextDepth);
                        rootEstimates.AddRange(subEstimates);
                        rootSingularities.AddRange(subSingularities);
                    }
                }

                continue;
            }

            (List<EigenPair> estimates, List<Complex> singularities) = RecursiveEstimateQNM(T, RT, beynPencil, domain, M, depth);
            rootEstimates.AddRange(estimates);
            rootSingularities.AddRange(singularities);

        }

        rootEstimates = rootEstimates.OrderBy(rEst => Complex.Abs(rEst.Value)).ThenBy(rEst => rEst.Value.Real).ToList();
        rootSingularities = rootSingularities.OrderBy(rSing => Complex.Abs(rSing)).ThenBy(rSing => rSing.Real).ToList();

        return (rootEstimates, rootSingularities);
    }

    //==================================================
    // 2. Recursive QNM locator
    //==================================================

    public (List<EigenPair>, List<Complex>) RecursiveEstimateQNM(Func<Complex, MatrixC> T, Func<Complex, MatrixC> RT, MatrixC beynPencil, SqDomainC domain, int M, Depth depth)
    {
        List<EigenPair> estimates = new List<EigenPair>();
        List<Complex> singularities = new List<Complex>();

        Func<Complex, Complex> MDZero = Moments.MomentDensity(T, RT, 0);
        (Complex MomentZero, List<Complex> localSingularities) = Calculus.LoopIntegral(MDZero, domain, M);
        int modeCount = ModeCount(MomentZero, beynPencil.Rows);
        Status(domain, modeCount, depth);

        if (modeCount < 0)
        {
            if (depth.Current >= depth.Max)
            {
                Logger.WriteLine($"\nSpurious moment {MomentZero.Format("E2")} returned. Skipping {domain.Format("F4")}.");
                Logger.WriteLine("");

                _area.processed += domain.Area;
                singularities.Add(domain.Center);
            }
            else
            {
                Logger.WriteLine($"\nSpurious moment {MomentZero.Format("E2")} returned. Subdividing {domain.Format("F4")}.");
                Logger.WriteLine("");

                _area.total -= domain.Area;
                domain = new SqDomainC(domain.Center, domain.Edge * 1.2);
                _area.total += domain.Area;

                Depth nextDepth = depth.Next();
                SqDomainC[] subdomains = domain.Subdivide(3);
                foreach (SqDomainC subdomain in subdomains)
                {
                    (List<EigenPair> subEstimates, List<Complex> subSingularities) = RecursiveEstimateQNM(T, RT, beynPencil, subdomain, M, nextDepth);
                    estimates.AddRange(subEstimates);
                    singularities.AddRange(subSingularities);
                }

            }

        }
        else if (modeCount == 0)
        {
            _area.processed += domain.Area;
        }
        else if (modeCount > 0 && modeCount <= beynPencil.Cols)
        {
            Func<Complex, MatrixC[]> BeynMDs = Moments.BeynMomentDensities(T, beynPencil, 1);
            (MatrixC[] BeynMoments, List<Complex> beynLocalSingularities) = Calculus.LoopIntegral(BeynMDs, domain, 8 * M);
            List<EigenPair> localEstimates = ExtractQNM(T, BeynMoments, modeCount);

            for (int j = 0; j < localEstimates.Count; j++)
            {
                Logger.WriteLine($"\nDiscovered {localEstimates[j].Value.Format("F12")} in {domain.Format("F4")}");
            }
            Logger.WriteLine("");

            _area.processed += domain.Area;

            estimates.AddRange(localEstimates);
            singularities.AddRange(beynLocalSingularities);

        }
        else // modeCount > beynPencil.Cols && modeCount <= beynPencil.Rows
        {
            if (depth.Current >= depth.Max)
            {
                Logger.WriteLine($"\nToo large moment {MomentZero.Format("E2")} returned. Skipping {domain.Format("F4")}.");
                Logger.WriteLine("");

                _area.processed += domain.Area;
                singularities.Add(domain.Center);
            }
            else
            {
                Logger.WriteLine($"\nToo large moment {MomentZero.Format("E2")} returned. Subdividing {domain.Format("F4")}.");
                Logger.WriteLine("");

                Depth nextDepth = depth.Next();
                SqDomainC[] subdomains = domain.Subdivide(2);

                foreach (SqDomainC subdomain in subdomains)
                {
                    (List<EigenPair> subEstimates, List<Complex> subSingularities) = RecursiveEstimateQNM(T, RT, beynPencil, subdomain, M, nextDepth);
                    estimates.AddRange(subEstimates);
                    singularities.AddRange(subSingularities);
                }
            }

        }

        singularities.AddRange(localSingularities);
        return (estimates, singularities);
    }

    //==================================================
    // 3. Extract EigenPairs from Beyn Moments
    //==================================================

    private static List<EigenPair> ExtractQNM(Func<Complex, MatrixC> T, MatrixC[] BeynMoments, int roots)
    {
        MatrixC BM0 = BeynMoments[0];
        MatrixC BM1 = BeynMoments[1];

        // SVD of BMZero
        var svd = BM0.Inner.Svd(computeVectors: true);
        MatrixC U = new MatrixC(svd.U.ToArray());
        MatrixC VT = new MatrixC(svd.VT.ToArray());
        var Sigma = svd.S; // singular values as a vector

        // Extract leading subspaces U_r, V_r, Σ_r
        int rows = BM0.Rows;
        int cols = BM0.Cols;

        MatrixC UR = MatrixC.Zero(rows, roots);
        MatrixC VR = MatrixC.Zero(cols, roots);

        MatrixC InvSigmaR = MatrixC.Zero(roots, roots);

        for (int i = 0; i < roots; i++)
        {
            for (int j = 0; j < rows; j++)
            {
                UR[j, i] = U[j, i];
            }

            for (int j = 0; j < cols; j++)
            {
                VR[j, i] = Complex.Conjugate(VT[i, j]);
            }

            InvSigmaR[i, i] = 1.0 / Sigma[i];
        }

        // Build reduced matrix B = U_r^* A1 V_r Σ_r^{-1}
        MatrixC UH = UR.Hermitian();
        MatrixC Beyn = UH * BM1 * VR * InvSigmaR;

        // Solve reduced eigenproblem B z = lambda z
        var evd = Beyn.Inner.Evd();
        var lambda = evd.EigenValues;
        var wRaw = evd.EigenVectors; // columns are eigenvectors z_k

        MatrixC RP = BM0 * VR * InvSigmaR;

        List<EigenPair> localEstimates = new List<EigenPair>();

        for (int k = 0; k < lambda.Count; k++)
        {
            VectorC w = new VectorC(wRaw.Column(k).ToArray());

            VectorC y = RP * w;
            y.Normalize();

            Complex omega = lambda[k];

            MatrixC TM = T(omega);
            VectorC r = TM * y;

            double OpNorm = TM.FrobeniusNorm();
            double residual = r.Norm() / OpNorm;

            localEstimates.Add(new EigenPair(omega, y, residual));
        }

        return localEstimates;
    }

    //==================================================
    // 4. Helper functions
    //==================================================

    // IsSommerfeld
    private bool DoesContainBP(SqDomainC domain)
    {
        Complex[,] systemBP = _kernel.SystemBP;

        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            { 
                if (domain.DoesContain(systemBP[i, j])) return true;
            }
        }

        return false;
    }

    // ModeCount from MomentZero
    private static int ModeCount(Complex momentZero, int maxModeCount)
    {
        double re = momentZero.Real;
        double im = momentZero.Imaginary;

        bool badRe = double.IsNaN(re) || double.IsInfinity(re);
        bool badIm = double.IsNaN(im) || double.IsInfinity(im);

        if (badRe || badIm) return -4;
        if (re < -1e-1) return -3;
        if (Math.Abs(im) > 1e-1) return -2;
        if (re > maxModeCount) return -1;
        return (int)Math.Round(re);
    }

    // Status update
    private void Status(SqDomainC domain, int modeCount, Depth depth)
    {
        double progress = _area.processed / _area.total;

        Console.Write($"\r{"",-40}\rProgressed {progress:P2}. ");
        Logger.WriteLine($"Found {modeCount} modes in {domain.Format("F4")} at depth {depth.Current}");

    }

}