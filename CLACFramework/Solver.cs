// File: Solver.cs

using System;
using System.Numerics;
using MathNet.Numerics;

namespace CLACFramework;

public class Solver
{
    public static void RecursiveSolve(ActionC A, MatrixC B, MatrixC X, int depth)
    {
        int ATotalRows = A.Rows;
        int ATotalCols = A.Cols;

        int BTotalRows = B.Rows;
        int BTotalCols = B.Cols;

        int XTotalRows = X.Rows;
        int XTotalCols = X.Cols;

        if (A == null || B == null || X == null) throw new ArgumentNullException();
        if (ATotalRows != ATotalCols) throw new ArgumentException("A must be quadrate.");
        if (BTotalRows != ATotalRows || XTotalRows != ATotalRows || BTotalCols != XTotalCols)
            throw new ArgumentException("Dimension mismatch between A, B, X.");

        // ------------------------------------------------------------
        // 1. Blockify A and B
        // ------------------------------------------------------------
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

        int ABlockRows = blockA.BlockRows;             // total block rows (interior + tail)
        int AInteriorBlockSize = blockA.RowSizes[0];   // interior block size

        int ALastBlockSize = blockA.RowSizes[ABlockRows - 1];
        int ATailBlockSize = (ALastBlockSize == AInteriorBlockSize ? 0 : ALastBlockSize);

        int AInteriorBlocks = (ATailBlockSize == 0 ? ABlockRows : ABlockRows - 1);

        if (ATailBlockSize < 0)
            throw new InvalidOperationException("Blockify produced inconsistent sizes.");

        // ------------------------------------------------------------
        // 2. Factor interior block operator
        // ------------------------------------------------------------
        BlockFactor Interior = BlockFactor.Factor(blockA, block0: 0, blockCount: AInteriorBlocks);

        // ------------------------------------------------------------
        // 3. Solve interior part
        // ------------------------------------------------------------
        var bM = new MatrixC[AInteriorBlocks];
        for (int b = 0; b < AInteriorBlocks; b++)
            bM[b] = blockB.GetBlock(b, 0); // thick column

        var xM = MatrixC.Zero(AInteriorBlocks, AInteriorBlockSize, BTotalCols);
        Interior.SolveLeaf(bM, xM);

        // ------------------------------------------------------------
        // 4. If no tail, assemble and return
        // ------------------------------------------------------------
        if (ATailBlockSize == 0)
        {
            for (int b = 0; b < AInteriorBlocks; b++)
            {
                var xb = xM[b];
                int rowOffset = b * AInteriorBlockSize;
                for (int j = 0; j < BTotalCols; j++)
                    for (int i = 0; i < AInteriorBlockSize; i++)
                        X[rowOffset + i, j] = xb[i, j];
            }
            return;
        }

        // ------------------------------------------------------------
        // 5. Tail exists → build top-level Schur complement
        // ------------------------------------------------------------
        MatrixC bt = blockB.GetBlock(AInteriorBlocks, 0); // tail RHS (R×nrhs)

        // A_mt: pow blocks of size Q×R
        var amt = new MatrixC[AInteriorBlocks];
        for (int b = 0; b < AInteriorBlocks; b++)
            amt[b] = blockA.GetBlock(b, AInteriorBlocks);

        // A_tm: pow blocks of size R×Q
        var atm = new MatrixC[AInteriorBlocks];
        for (int b = 0; b < AInteriorBlocks; b++)
            atm[b] = blockA.GetBlock(AInteriorBlocks, b);

        // A_tt: tail block (R×R)
        MatrixC att = blockA.GetBlock(AInteriorBlocks, AInteriorBlocks);

        // Compute A_tm * xM
        MatrixC atm_xM = new MatrixC(ATailBlockSize, BTotalCols);
        atm_xM.Clear();
        for (int b = 0; b < AInteriorBlocks; b++)
            atm[b].GemmIP(xM[b], atm_xM, Complex.One, Complex.One);

        // Effective tail RHS
        MatrixC bTailEff = new MatrixC(ATailBlockSize, BTotalCols);
        for (int j = 0; j < BTotalCols; j++)
            for (int i = 0; i < ATailBlockSize; i++)
                bTailEff[i, j] = bt[i, j] - atm_xM[i, j];

        // Build Schur complement S = A_tt - A_tm A11^{-1} A_mt
        MatrixC S = BuildTopLevelSchur(Interior, amt, atm, att, AInteriorBlocks, AInteriorBlockSize, ATailBlockSize);

        // Solve S xTail = bTailEff
        var luS = S.LU();
        MatrixC xTail = new MatrixC(ATailBlockSize, BTotalCols);
        luS.SolveIP(bTailEff, xTail);

        // ------------------------------------------------------------
        // 6. Correct interior solution: xM = xM - A11^{-1}(A_mt xTail)
        // ------------------------------------------------------------
        var amt_xTail = new MatrixC[AInteriorBlocks];
        for (int b = 0; b < AInteriorBlocks; b++)
        {
            amt_xTail[b] = new MatrixC(AInteriorBlockSize, BTotalCols);
            amt[b].GemmIP(xTail, amt_xTail[b], Complex.One, Complex.Zero);
        }

        var corr = MatrixC.Zero(AInteriorBlocks, AInteriorBlockSize, BTotalCols);
        Interior.SolveLeaf(amt_xTail, corr);

        for (int b = 0; b < AInteriorBlocks; b++)
        {
            var xb = xM[b];
            var cb = corr[b];
            for (int j = 0; j < BTotalCols; j++)
                for (int i = 0; i < AInteriorBlockSize; i++)
                    xb[i, j] -= cb[i, j];
        }

        // ------------------------------------------------------------
        // 7. Assemble final solution X
        // ------------------------------------------------------------
        for (int b = 0; b < AInteriorBlocks; b++)
        {
            var xb = xM[b];
            int rowOffset = b * AInteriorBlockSize;
            for (int j = 0; j < BTotalCols; j++)
                for (int i = 0; i < AInteriorBlockSize; i++)
                    X[rowOffset + i, j] = xb[i, j];
        }

        for (int j = 0; j < BTotalCols; j++)
            for (int i = 0; i < ATailBlockSize; i++)
                X[AInteriorBlocks * AInteriorBlockSize + i, j] = xTail[i, j];

        blockA.Erase();
        blockB.Erase();
    }

    // ----------------------------------------------------------------
    // Build top-level Schur complement S = A_tt - A_tm A11^{-1} A_mt
    // ----------------------------------------------------------------
    private static MatrixC BuildTopLevelSchur(BlockFactor Interior, MatrixC[] amt, MatrixC[] atm, MatrixC att,
        int interiorBlocks, int blockRows, int tailRows)
    {
        // Compute Z = A11^{-1} A_mt
        var zBlocks = MatrixC.Zero(interiorBlocks, blockRows, tailRows);
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