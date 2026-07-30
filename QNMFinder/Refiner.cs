// File: Refiner.cs

using System;
using System.Numerics;
using MathNet.Numerics;
using IOFramework;
using CLACFramework;

namespace QNMFinder;

public sealed class Refiner
{
    private readonly Kernel _kernel;

    public Refiner(Kernel kernel)
    {
        _kernel = kernel;
    }
    
    // Refine a list of estimated eigenpairs
    public List<EigenPair> RefineQNM(List<EigenPair> seeds, int I)
    {
        List<EigenPair> refined = new List<EigenPair>();

        for (int k = 0; k < seeds.Count; k++)
        {
            EigenPair seed = seeds[k];

            Logger.WriteBoth($"\nMode {k:D2}: Refining frequency {seed.Value.Format("F16")}...");
            Logger.WriteLine("");
            
            EigenPair refinedEP = RefineSingleQNM(seed, I);

            refined.Add(refinedEP);
        }

        refined = refined.OrderBy(rEP => Complex.Abs(rEP.Value)).ThenBy(rEP => rEP.Value.Real).ToList();

        return refined;
    }

    public EigenPair RefineSingleQNM(EigenPair seed, int I)
    {
        (Complex omega, VectorC y, double residual) = seed;

        int NP = y.Size;
        int N = NP - 1;

        // Assemble system operators and functions
        (VectorC x, MatrixC D) = Discretizer.DiscretizeM(N);
        (Func<Complex, MatrixC> T, Func<Complex, MatrixC> RT) = _kernel.SystemOperatorMatrices(x, D);

        MatrixC TM = T(omega);
        MatrixC RTM = RT(omega);

        VectorC r = VectorC.Zero(NP);
        double OpNorm = 0.0;

        for (int i = 0; i < I; i++)
        {
            (Complex domega, VectorC dy) = NewtonUpdate(TM, RTM, y);

            omega += domega;
            y += dy;

            y.Normalize();

            // Update operators
            TM = T(omega);
            RTM = RT(omega);

            // Scalar residual
            r = TM * y;
            OpNorm = TM.FrobeniusNorm();
            residual = r.Norm() / OpNorm;

            Logger.WriteLine($"Iter {i:D2}: Freq = {omega.Format("F16")}, Res = {residual.Format("E6")}");

        }

        EigenPair refinedEP = new EigenPair(omega, y, residual);
        return refinedEP;
    }

    public static (Complex, VectorC) NewtonUpdate(MatrixC TM, MatrixC RTM, VectorC y)
    {
        int NP = TM.Rows;
        int NPP = NP + 1;

        int[] rowSizes = { NP, 1 };
        int[] colSizes = { NP, 1 };

        VectorC r = TM * y;
        VectorC s = RTM * y;

        BlockC bTMG = new BlockC(rowSizes, colSizes);

        // Top-left block: T(omega)
        bTMG.SetBlock(0, 0, TM);

        // Top-right block: s (as an n×1 matrix)
        MatrixC sG = s.ToMatrixC();
        bTMG.SetBlock(0, 1, sG);

        // Bottom-left block: y^H (as a 1×n matrix)
        MatrixC yG = y.Hermitian();
        bTMG.SetBlock(1, 0, yG);

        // Bottom-right block: 0 (1×1)
        bTMG.SetBlock(1, 1, MatrixC.Zero(1, 1));

        // Assemble the full augmented matrix
        MatrixC TMG = bTMG.Assemble();

        // Build rG = (r, 0)
        VectorC rG = new VectorC(NPP);
        rG.SetSubvector(0, r);
        rG[NP] = Complex.Zero;

        // Solve the augmented system: TMG * dyG = - rG
        VectorC dyG = - TMG.Solve(rG);

        // Extract dy and domega
        VectorC dy = dyG.GetSubvector(0, NP);
        Complex domega = dyG[NP];

        return (domega, dy);
    }

}