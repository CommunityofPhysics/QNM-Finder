// File: Guesser.cs

using System;
using System.Numerics;
using CLACFramework;
using IOFramework;

namespace QNMFinder;

public sealed class Guesser
{
    private readonly Kernel _kernel;

    public Guesser(Kernel kernel)
    {
        _kernel = kernel;
    }

    // ============================================================
    // Main QNM Guesser
    // ============================================================
    public List<EigenPair> GuessQNM(Complex W, int N, double A)
    {
        int NM = N - 1;
        int NP = N + 1;

        double wA = A * Math.Sqrt(N);

        // 1. Discretization
        (VectorC x, MatrixC D) = Discretizer.DiscretizeM(N);

        // 2. Build BC-imposed quadratic operator via Kernel
        (Func<Complex, MatrixC> T, Func<Complex, MatrixC> RT) = _kernel.SystemOperatorMatrices(x, D);
        Func<Complex, MatrixC> RRT = Calculus.Derivative(RT);

        W = FindRegularW(T, W);

        // Quadratic expansion: T(W + lambda) = P + lambda Q + lambda^2 R
        MatrixC P = T(W);
        MatrixC Q = RT(W);
        MatrixC R = 0.5 * RRT(W);

        // 3. Static condensation + linearization
        (MatrixC L, MatrixC PBI, MatrixC QBI, MatrixC RBI, MatrixC PBB, MatrixC QBB, MatrixC RBB) 
            = StaticCondenser(P, Q, R);

        // 4. Solve reduced eigenproblem L Y = ω Y
        var evd = L.Inner.Evd();
        var Lambda = evd.EigenValues;
        var bY = evd.EigenVectors;

        MatrixC Y = new MatrixC(bY);

        // 5. Reconstruct full eigenvectors and filter
        List<EigenPair> guesses = new List<EigenPair>();

        for (int i = 0; i < Lambda.Count; i++)
        {
            Complex lambda = Lambda[i];
            Complex omega = W + lambda;

            // Interior part of y is the first NM entries of Y
            VectorC yInt = Y.GetColumn(i).GetSubvector(0, NM);

            // --- Build yB ---
            MatrixC TBB = PBB + lambda * QBB + lambda * lambda * RBB;

            MatrixC TBI = PBI + lambda * QBI + lambda * lambda * RBI;
            VectorC rhs = TBI * yInt;

            VectorC yB = - TBB.LU().Solve(rhs);

            // Build full vector
            VectorC y = VectorC.Zero(NP);

            y[0] = yB[0];
            y.SetSubvector(1, yInt);
            y[N] = yB[1];

            y.Normalize();

            // Residual check with full T(ω)
            MatrixC TM = T(omega);
            VectorC r = TM * y;

            double OpNorm = TM.FrobeniusNorm();
            double residual = r.Norm() / OpNorm;

            if (residual < wA)
                guesses.Add(new EigenPair(omega, y, residual));
        }

        guesses = guesses.OrderBy(gue => gue.Value.Magnitude).ToList();

        return guesses;
    }

    // ============================================================
    // Static Condenser
    // ============================================================
    private static (MatrixC L, MatrixC PBI, MatrixC QBI, MatrixC RBI, MatrixC PBB, MatrixC QBB, MatrixC RBB)
        StaticCondenser(MatrixC P, MatrixC Q, MatrixC R)
    {
        int NP = P.Rows;
        int N = NP - 1;      // boundary indices: 0 and N
        int NM = N - 1;       // interior size

        // ------------------------------------------------------------
        // Block extraction: II, IB, BI, BB
        // ------------------------------------------------------------

        // Interior–interior blocks (NM×NM): indices 1..N-1
        MatrixC PII = MatrixC.Zero(NM, NM);
        MatrixC QII = MatrixC.Zero(NM, NM);
        MatrixC RII = MatrixC.Zero(NM, NM);

        for (int i = 1; i < N; i++)
        {
            int iM = i - 1;   // interior row index 0..NM-1

            for (int j = 1; j < N; j++)
            {
                int jM = j - 1;   // interior col index 0..NM-1

                PII[iM, jM] = P[i, j];
                QII[iM, jM] = Q[i, j];
                RII[iM, jM] = R[i, j];
            }
        }

        // Interior–boundary blocks (NM×2): cols 0 and N
        MatrixC PIB = MatrixC.Zero(NM, 2);
        MatrixC QIB = MatrixC.Zero(NM, 2);
        MatrixC RIB = MatrixC.Zero(NM, 2);

        for (int i = 1; i < N; i++)
        {
            int iM = i - 1;
            PIB[iM, 0] = P[i, 0];
            PIB[iM, 1] = P[i, N];

            QIB[iM, 0] = Q[i, 0];
            QIB[iM, 1] = Q[i, N];

            RIB[iM, 0] = R[i, 0];
            RIB[iM, 1] = R[i, N];
        }

        // Boundary–interior blocks (2×NM): rows 0 and N, cols 1..N-1
        MatrixC PBI = MatrixC.Zero(2, NM);
        MatrixC QBI = MatrixC.Zero(2, NM);
        MatrixC RBI = MatrixC.Zero(2, NM);

        for (int j = 1; j < N; j++)
        {
            int jM = j - 1;
            PBI[0, jM] = P[0, j];
            PBI[1, jM] = P[N, j];

            QBI[0, jM] = Q[0, j];
            QBI[1, jM] = Q[N, j];

            RBI[0, jM] = R[0, j];
            RBI[1, jM] = R[N, j];
        }

        // Boundary–boundary blocks (2×2): indices {0,N}×{0,N}
        MatrixC PBB = MatrixC.Zero(2, 2);
        MatrixC QBB = MatrixC.Zero(2, 2);
        MatrixC RBB = MatrixC.Zero(2, 2);

        PBB[0, 0] = P[0, 0];
        PBB[0, 1] = P[0, N];
        PBB[1, 0] = P[N, 0];
        PBB[1, 1] = P[N, N];

        QBB[0, 0] = Q[0, 0];
        QBB[0, 1] = Q[0, N];
        QBB[1, 0] = Q[N, 0];
        QBB[1, 1] = Q[N, N];

        RBB[0, 0] = R[0, 0];
        RBB[0, 1] = R[0, N];
        RBB[1, 0] = R[N, 0];
        RBB[1, 1] = R[N, N];

        // ------------------------------------------------------------
        // Quadratic Schur complement
        // ------------------------------------------------------------

        // Precompute products
        MatrixC invPBB = PBB.Solve(MatrixC.Identity(2));
        MatrixC PIB_invPBB = PIB * invPBB;
        MatrixC invPBB_PBI = invPBB * PBI;
        MatrixC QIB_invPBB = QIB * invPBB;
        MatrixC RIB_invPBB = RIB * invPBB;

        // \tilde{P} = P_II - P_IB P_BB^{-1} P_BI
        MatrixC PIIc = PII - PIB_invPBB * PBI;

        // \tilde{Q} = Q_II - Q_IB P_BB^{-1} P_BI - P_IB P_BB^{-1} Q_BI
        //            + P_IB P_BB^{-1} Q_BB P_BB^{-1} P_BI
        MatrixC QIIc = QII - QIB_invPBB * PBI - PIB_invPBB * QBI 
                       + PIB_invPBB * QBB * invPBB_PBI;

        // \tilde{R} = R_II - R_IB P_BB^{-1} P_BI - P_IB P_BB^{-1} R_BI
        //           + P_IB P_BB^{-1} R_BB P_BB^{-1} P_BI - Q_IB P_BB^{-1} Q_BI
        //           + Q_IB P_BB^{-1} Q_BB P_BB^{-1} P_BI + P_IB P_BB^{-1} Q_BB P_BB^{-1} Q_BI
        //           - P_IB P_BB^{-1} Q_BB P_BB^{-1} Q_BB P_BB^{-1} P_BI
        MatrixC RIIc = RII - RIB_invPBB * PBI - PIB_invPBB * RBI 
                       + PIB_invPBB * RBB * invPBB_PBI - QIB_invPBB * QBI 
                       + QIB_invPBB * QBB * invPBB_PBI + PIB_invPBB * QBB * (invPBB * QBI) 
                       - PIB_invPBB * QBB * invPBB * QBB * invPBB_PBI;

        // ------------------------------------------------------------
        // Linearization: L = B^{-1} A for condensed quadratic pencil
        // ------------------------------------------------------------
        LUFactorC RLU = RIIc.LU();
        MatrixC RP = RLU.Solve(PIIc);
        MatrixC RQ = RLU.Solve(QIIc);

        BlockC blkL = new BlockC(new[] { NM, NM }, new[] { NM, NM });
        blkL.SetBlock(0, 0, MatrixC.Zero(NM, NM));
        blkL.SetBlock(0, 1, MatrixC.Identity(NM));
        blkL.SetBlock(1, 0, -RP);
        blkL.SetBlock(1, 1, -RQ);

        MatrixC L = blkL.Assemble();

        return (L, PBI, QBI, RBI, PBB, QBB, RBB);
    }

    // ============================================================
    // W Finder
    // ============================================================
    private static Complex FindRegularW(Func<Complex, MatrixC> T, Complex W)
    {
        const double eps = 1e-10;
        const int maxIter = 10;

        Random rng = new Random();

        for (int iter = 0; iter < maxIter; iter++)
        {
            // Build T(W)
            MatrixC TW = T(W);
            int NP = TW.Rows;
            int N = NP - 1;

            // Extract 2×2 boundary block
            Complex a = TW[0, 0];
            Complex b = TW[0, N];
            Complex c = TW[N, 0];
            Complex d = TW[N, N];

            Complex det = a * d - b * c;

            if (det.Magnitude > eps)
                return W;

            // Otherwise perturb W slightly
            double re = (rng.NextDouble() - 0.5);
            double im = (rng.NextDouble() - 0.5);

            W += new Complex(re, im);
        }

        throw new Exception("Failed to find invertible boundary block TBB(W) after 10 perturbations.");
    }

}