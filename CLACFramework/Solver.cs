// File: Solver.cs

using System;
using System.Numerics;
using MathNet.Numerics;

namespace CLACFramework;

public class Solver
{
    // ================================================================
    // Main Recursive Solve
    // ================================================================
    public static void RecursiveSolve(ActionC A, MatrixC B, MatrixC X, int depth)
    {
        if (A == null || B == null || X == null)
            throw new ArgumentNullException();

        int ATotalRows = A.Rows;
        int ATotalCols = A.Cols;

        int BTotalRows = B.Rows;
        int BTotalCols = B.Cols;

        int XTotalRows = X.Rows;
        int XTotalCols = X.Cols;

        if (ATotalRows != ATotalCols) throw new ArgumentException("A must be quadrate.");
        if (BTotalRows != ATotalRows || XTotalRows != ATotalRows || BTotalCols != XTotalCols)
            throw new ArgumentException("Dimension mismatch between A, B, X.");

        Func<int, int, BSType> policy = (bi, bj) =>
        {
            if (bi == bj) return BSType.InMemory;
            if (Math.Abs(bi - bj) == 1) return BSType.InCache;
            return BSType.OnFile;
        };

        BSSpec specA = new BSSpec(quadrate: true, depth, policy);
        BlockC blockA = Blockify.Action(A, specA);

        BSSpec specB = new BSSpec(quadrate: false, depth, policy);
        BlockC blockB = Blockify.Matrix(B, specB);

        try
        {
            int ABlockRows = blockA.BlockRows;
            int AInteriorBlockSize = blockA.RowSizes[0];

            int ALastBlockSize = blockA.RowSizes[ABlockRows - 1];
            int ATailBlockSize = (ALastBlockSize == AInteriorBlockSize ? 0 : ALastBlockSize);

            int AInteriorBlocks = (ATailBlockSize == 0 ? ABlockRows : ABlockRows - 1);

            if (ATailBlockSize < 0)
                throw new InvalidOperationException("Blockify produced inconsistent sizes.");

            BlockFactor Interior = BlockFactor.Factor(blockA, block0: 0, blockCount: AInteriorBlocks);

            MatrixC[] bM = new MatrixC[AInteriorBlocks];
            for (int b = 0; b < AInteriorBlocks; b++)
                bM[b] = blockB.GetBlock(b, 0);

            MatrixC[] xM = MatrixC.Zero(AInteriorBlocks, AInteriorBlockSize, BTotalCols);
            Interior.SolveLeaf(bM, xM);

            if (ATailBlockSize == 0)
            {
                for (int b = 0; b < AInteriorBlocks; b++)
                {
                    MatrixC xb = xM[b];
                    int rowOffset = b * AInteriorBlockSize;
                    for (int j = 0; j < BTotalCols; j++)
                        for (int i = 0; i < AInteriorBlockSize; i++)
                            X[rowOffset + i, j] = xb[i, j];
                }
                return;
            }

            MatrixC bt = blockB.GetBlock(AInteriorBlocks, 0);

            MatrixC[] amt = new MatrixC[AInteriorBlocks];
            for (int b = 0; b < AInteriorBlocks; b++)
                amt[b] = blockA.GetBlock(b, AInteriorBlocks);

            MatrixC[] atm = new MatrixC[AInteriorBlocks];
            for (int b = 0; b < AInteriorBlocks; b++)
                atm[b] = blockA.GetBlock(AInteriorBlocks, b);

            MatrixC att = blockA.GetBlock(AInteriorBlocks, AInteriorBlocks);

            MatrixC atm_xM = new MatrixC(ATailBlockSize, BTotalCols);
            atm_xM.Clear();
            for (int b = 0; b < AInteriorBlocks; b++)
                atm[b].GemmIP(xM[b], atm_xM, Complex.One, Complex.One);

            MatrixC bTailEff = new MatrixC(ATailBlockSize, BTotalCols);
            for (int j = 0; j < BTotalCols; j++)
                for (int i = 0; i < ATailBlockSize; i++)
                    bTailEff[i, j] = bt[i, j] - atm_xM[i, j];

            MatrixC S = BuildTopLevelSchur(Interior, amt, atm, att, AInteriorBlocks, AInteriorBlockSize, ATailBlockSize);

            LUFactorC luS = S.LU();

            MatrixC xTail = new MatrixC(ATailBlockSize, BTotalCols);
            luS.SolveIP(bTailEff, xTail);

            MatrixC[] amt_xTail = new MatrixC[AInteriorBlocks];
            for (int b = 0; b < AInteriorBlocks; b++)
            {
                amt_xTail[b] = new MatrixC(AInteriorBlockSize, BTotalCols);
                amt[b].GemmIP(xTail, amt_xTail[b], Complex.One, Complex.Zero);
            }

            MatrixC[] corr = MatrixC.Zero(AInteriorBlocks, AInteriorBlockSize, BTotalCols);
            Interior.SolveLeaf(amt_xTail, corr);

            for (int b = 0; b < AInteriorBlocks; b++)
            {
                MatrixC xb = xM[b];
                MatrixC cb = corr[b];
                for (int j = 0; j < BTotalCols; j++)
                    for (int i = 0; i < AInteriorBlockSize; i++)
                        xb[i, j] -= cb[i, j];
            }

            for (int b = 0; b < AInteriorBlocks; b++)
            {
                MatrixC xb = xM[b];
                int rowOffset = b * AInteriorBlockSize;
                for (int j = 0; j < BTotalCols; j++)
                    for (int i = 0; i < AInteriorBlockSize; i++)
                        X[rowOffset + i, j] = xb[i, j];
            }

            for (int j = 0; j < BTotalCols; j++)
                for (int i = 0; i < ATailBlockSize; i++)
                    X[AInteriorBlocks * AInteriorBlockSize + i, j] = xTail[i, j];
        }
        finally
        {
            blockA.Erase();
            blockB.Erase();
        }
    }

    // ================================================================
    // Build top-level Schur complement S = A_tt - A_tm A11^{-1} A_mt
    // ================================================================
    private static MatrixC BuildTopLevelSchur(BlockFactor Interior, MatrixC[] amt, MatrixC[] atm, MatrixC att,
        int interiorBlocks, int blockRows, int tailRows)
    {
        // Compute Z = A11^{-1} A_mt
        MatrixC[] zBlocks = MatrixC.Zero(interiorBlocks, blockRows, tailRows);
        Interior.SolveLeaf(amt, zBlocks);

        // Compute A_tm * Z
        MatrixC atmZ = new MatrixC(tailRows, tailRows);
        atmZ.Clear();
        for (int b = 0; b < interiorBlocks; b++)
            atm[b].GemmIP(zBlocks[b], atmZ, Complex.One, Complex.One);

        // S = A_tt - A_tm Z
        MatrixC S = new MatrixC(tailRows, tailRows);
        for (int j = 0; j < tailRows; j++)
            for (int i = 0; i < tailRows; i++)
                S[i, j] = att[i, j] - atmZ[i, j];

        return S;
    }

}