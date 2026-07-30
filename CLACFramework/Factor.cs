// File: Factor.cs

using System;
using System.Numerics;
using MathNet.Numerics;

namespace CLACFramework;

public sealed class BlockFactor
{
    private readonly BlockC blocks;        // global block matrix (shared)
    private readonly int block0;          // starting block index of this submatrix
    private readonly int blockCount;      // number of block rows/cols in this submatrix
    private readonly int blockRows;       // size of each block (assumed uniform here)

    // factorization pieces
    private readonly LUFactorC? leafLU;   // non-null only for leaf nodes
    private readonly BlockFactor? f11;    // factorization of A11 (top-left)
    private readonly BlockFactor? fS;     // factorization of Schur complement S (bottom-right)
    private readonly BlockC? schurBlocks; // explicit Schur complement blocks (for interior nodes)

    // Format: [i0, j0, i1, j1, ...] where each pair (ik, jk) is a symmetric swap of block rows/cols ik <-> jk.
    private readonly int[] swaps;

    private BlockFactor(BlockC blocks, int block0, int blockCount, int blockRows, 
        LUFactorC? leafLU, BlockFactor? f11, BlockFactor? fS, BlockC? schurBlocks, int[] swaps)
    {
        this.blocks = blocks;
        this.block0 = block0;
        this.blockCount = blockCount;
        this.blockRows = blockRows;
        this.leafLU = leafLU;
        this.f11 = f11;
        this.fS = fS;
        this.schurBlocks = schurBlocks;
        this.swaps = swaps ?? Array.Empty<int>();
    }

    // ----------------------------------------------------------------
    // Public entry: factor a block-submatrix starting at block0 with blockCount blocks
    // ----------------------------------------------------------------
    public static BlockFactor Factor(BlockC blockM, int block0, int blockCount)
    {
        if (blockCount <= 0)
            throw new ArgumentException("blockCount must be positive.");

        int blockRows = blockM.RowSizes[block0];

        // uniform block size check (square blocks)
        for (int k = block0; k < block0 + blockCount; k++)
        {
            if (blockM.RowSizes[k] != blockRows || blockM.ColSizes[k] != blockRows)
                throw new InvalidOperationException("Non-uniform block size not supported in BlockFactor.");
        }

        // Leaf case
        if (blockCount == 1)
        {
            // Try to factor a diagonal block within the submatrix range.
            // We allow searching the diagonal blocks in [block0, block0+blockCount-1]
            // (here blockCount==1 so only block0), but code supports a general scan.
            (LUFactorC lu, int[] localSwaps) = FactorPivotLeaf(blockM, block0, blockCount);
            return new BlockFactor(blockM, block0, 1, blockRows, lu, null, null, null, localSwaps);
        }

        // Recursive case
        int n1 = blockCount / 2;
        int n2 = blockCount - n1;

        // Factor A11 (top-left)
        BlockFactor f11 = Factor(blockM, block0, n1);

        // Build Schur complement S = A22 - A21 * A11^{-1} * A12
        BlockC S = BuildSchurBlock(blockM, f11, block0, n1, n2, blockRows);

        // Factor S (note: S is a fresh BlockC with local indexing 0..n2-1)
        BlockFactor fS = Factor(S, 0, n2);

        // Compose swaps: first all swaps that happened in f11 subtree, then those in fS subtree.
        // (These are global block index swaps recorded during factorization.)
        int[] composedSwaps = ConcatSwaps(f11.swaps, fS.swaps);

        return new BlockFactor(blockM, block0, blockCount, blockRows, null, f11, fS, S, composedSwaps);
    }

    // ----------------------------------------------------------------
    // Leaf factorization with block pivot search and symmetric swap
    // Returns LU factor for the chosen diagonal block and the list of swaps performed (flattened pairs).
    // ----------------------------------------------------------------
    private static (LUFactorC lu, int[] swaps) FactorPivotLeaf(BlockC blockM, int block0, int blockCount)
    {
        int start = block0;
        int end = block0 + blockCount; // exclusive
        LUFactorC? lu = null;
        int pivotIndex = -1;
        var swapsList = new System.Collections.Generic.List<int>();

        // Try diagonal blocks in the submatrix range until one LU succeeds.
        for (int k = start; k < end; k++)
        {
            MatrixC diag = blockM.GetBlock(k, k);
            try
            {
                lu = diag.LU();
                pivotIndex = k;
                break;
            }
            catch
            {
                // try next diagonal block
            }
        }

        if (lu == null)
            throw new InvalidOperationException("No LU-able diagonal block found for factorization.");

        // If pivot is not at the top of this submatrix (block0), bring it to block0 by symmetric swap.
        if (pivotIndex != block0)
        {
            SymmetricSwap(blockM, block0, pivotIndex);
            // record the global swap pair
            swapsList.Add(block0);
            swapsList.Add(pivotIndex);
        }

        // Return LU of the (now) top diagonal block and the swaps performed.
        return (lu!, swapsList.ToArray());
    }

    // ----------------------------------------------------------------
    // Symmetric swap of block rows and columns i <-> j (global indices)
    // ----------------------------------------------------------------
    private static void SymmetricSwap(BlockC blockM, int i, int j)
    {
        if (i == j) return;
        int nBlocks = blockM.BlockRows;

        // swap block rows i and j
        for (int k = 0; k < nBlocks; k++)
        {
            MatrixC a = blockM.GetBlock(i, k);
            MatrixC b = blockM.GetBlock(j, k);
            blockM.SetBlock(i, k, b);
            blockM.SetBlock(j, k, a);
        }

        // swap block cols i and j
        for (int k = 0; k < nBlocks; k++)
        {
            MatrixC a = blockM.GetBlock(k, i);
            MatrixC b = blockM.GetBlock(k, j);
            blockM.SetBlock(k, i, b);
            blockM.SetBlock(k, j, a);
        }
    }

    // ----------------------------------------------------------------
    // Concatenate two flattened swap arrays (pairs)
    // ----------------------------------------------------------------
    private static int[] ConcatSwaps(int[] a, int[] b)
    {
        if ((a == null || a.Length == 0) && (b == null || b.Length == 0))
            return Array.Empty<int>();
        if (a == null || a.Length == 0) return (int[])b.Clone();
        if (b == null || b.Length == 0) return (int[])a.Clone();

        int[] outArr = new int[a.Length + b.Length];
        Array.Copy(a, 0, outArr, 0, a.Length);
        Array.Copy(b, 0, outArr, a.Length, b.Length);
        return outArr;
    }

    // ----------------------------------------------------------------
    // Apply a sequence of swaps (flattened pairs) to a MatrixC[] array in-place.
    // The swaps are applied in the order they appear in the flattened array.
    // ----------------------------------------------------------------
    private static void ApplySwapsInPlace(MatrixC[] arr, int[] swaps)
    {
        if (swaps == null || swaps.Length == 0) return;
        for (int p = 0; p < swaps.Length; p += 2)
        {
            int i = swaps[p];
            int j = swaps[p + 1];
            // swap arr[i] and arr[j]
            var tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
        }
    }

    // ----------------------------------------------------------------
    // Undo a sequence of swaps in-place by applying them in reverse order.
    // ----------------------------------------------------------------
    private static void UndoSwapsInPlace(MatrixC[] arr, int[] swaps)
    {
        if (swaps == null || swaps.Length == 0) return;
        for (int p = swaps.Length - 2; p >= 0; p -= 2)
        {
            int i = swaps[p];
            int j = swaps[p + 1];
            var tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
        }
    }

    // ----------------------------------------------------------------
    // Build Schur complement S = A22 - A21 * A11^{-1} * A12
    // ----------------------------------------------------------------
    private static BlockC BuildSchurBlock(BlockC blockM, BlockFactor f11, int block0, int n1, int n2, int blockRows)
    {
        int[] rowSizes = new int[n2];
        int[] colSizes = new int[n2];
        for (int i = 0; i < n2; i++)
        {
            rowSizes[i] = blockRows;
            colSizes[i] = blockRows;
        }

        BlockC sBlocks = new BlockC(rowSizes, colSizes);

        for (int jb = 0; jb < n2; jb++)
        {
            // rhs = A12(:, jb)
            MatrixC[] rhs = new MatrixC[n1];
            MatrixC[] z = MatrixC.Zero(n1, blockRows, blockRows);

            for (int kb = 0; kb < n1; kb++)
                rhs[kb] = blockM.GetBlock(block0 + kb, block0 + n1 + jb);

            // apply A11^{-1} to each column block of A12
            f11.SolveLeaf(rhs, z);

            // assemble S column jb
            for (int ib = 0; ib < n2; ib++)
            {
                MatrixC s_ij = new MatrixC(blockRows, blockRows);
                s_ij.CopyFrom(blockM.GetBlock(block0 + n1 + ib, block0 + n1 + jb)); // start with A22(ib,jb)

                for (int kb = 0; kb < n1; kb++)
                {
                    MatrixC a21 = blockM.GetBlock(block0 + n1 + ib, block0 + kb);
                    MatrixC zk = z[kb];
                    a21.GemmIP(zk, s_ij, -Complex.One, Complex.One); // s_ij -= a21 * zk
                }

                sBlocks.SetBlock(ib, jb, s_ij);
            }
        }

        return sBlocks;
    }

    // ----------------------------------------------------------------
    // Solve: apply recorded swaps to RHS, perform block Schur solve, undo swaps on solution
    // ----------------------------------------------------------------
    public void SolveLeaf(MatrixC[] bBlocks, MatrixC[] xBlocks)
    {
        if (bBlocks == null || xBlocks == null)
            throw new ArgumentNullException();
        if (bBlocks.Length != blockCount || xBlocks.Length != blockCount)
            throw new ArgumentException("Block array length mismatch.");

        int pencilRows = bBlocks[0].Cols;

        // dimension checks
        for (int i = 0; i < blockCount; i++)
        {
            if (bBlocks[i].Rows != blockRows || bBlocks[i].Cols != pencilRows)
                throw new ArgumentException("bBlocks dimension mismatch.");
            if (xBlocks[i].Rows != blockRows || xBlocks[i].Cols != pencilRows)
                throw new ArgumentException("xBlocks dimension mismatch.");
        }

        // Make a local copy of RHS blocks so we can permute them in-place without affecting caller's array.
        // We will write results into xBlocks after undoing swaps.
        MatrixC[] bLocal = new MatrixC[blockCount];
        for (int i = 0; i < blockCount; i++) bLocal[i] = bBlocks[i];

        // Apply all swaps that were recorded during factorization (these are global indices).
        // However, our bLocal is indexed 0..blockCount-1 relative to this submatrix.
        // The recorded swaps are global indices; we need to map them into local indices for this submatrix.
        // Only swaps that touch indices inside [block0, block0+blockCount-1] are relevant here.
        // We'll build a local flattened swap list (pairs of local indices) in the same order.
        int[] localSwaps = ExtractLocalSwaps(swaps, block0, blockCount);

        // Apply local swaps to bLocal
        ApplySwapsInPlace(bLocal, localSwaps);

        // Base case
        if (blockCount == 1)
        {
            // Solve leaf LU in-place: leafLU solves A * x = bLocal[0] -> xBlocks[0]
            leafLU!.SolveIP(bLocal[0], xBlocks[0]);
            // No need to undo swaps because localSwaps is empty for blockCount==1 unless a swap involved this single block.
            // But we applied localSwaps already; now we must undo them on the solution before returning.
            // Build xPerm array and undo local swaps on it.
            MatrixC[] xPerm = new MatrixC[1] { xBlocks[0] };
            UndoSwapsInPlace(xPerm, localSwaps);
            xBlocks[0].CopyFrom(xPerm[0]);
            return;
        }

        // Recursive solve for blockCount > 1
        int n1 = blockCount / 2;
        int n2 = blockCount - n1;

        // Split bLocal into b1 and b2 (these are local to this submatrix)
        MatrixC[] b1 = new MatrixC[n1];
        MatrixC[] b2 = new MatrixC[n2];
        Array.Copy(bLocal, 0, b1, 0, n1);
        Array.Copy(bLocal, n1, b2, 0, n2);

        // Allocate solution blocks
        MatrixC[] x1 = MatrixC.Zero(n1, blockRows, pencilRows);
        MatrixC[] x2 = MatrixC.Zero(n2, blockRows, pencilRows);

        // Solve A11 * x1 = b1
        f11!.SolveLeaf(b1, x1);

        // a21 * x1
        MatrixC[] a21x1 = MatrixC.Zero(n2, blockRows, pencilRows);
        BlockMatrixMultiplication(blocks, block0 + n1, n2, block0, n1, x1, a21x1);

        // r2 = b2 - a21x1
        MatrixC[] r2 = MatrixC.Zero(n2, blockRows, pencilRows);
        for (int ib = 0; ib < n2; ib++)
        {
            MatrixC rb = r2[ib];
            MatrixC bb = b2[ib];
            MatrixC ax = a21x1[ib];

            for (int j = 0; j < pencilRows; j++)
                for (int i = 0; i < blockRows; i++)
                    rb[i, j] = bb[i, j] - ax[i, j];
        }

        // Solve S * x2 = r2
        fS!.SolveLeaf(r2, x2);

        // w = a12 * x2
        MatrixC[] w = MatrixC.Zero(n1, blockRows, pencilRows);
        BlockMatrixMultiplication(blocks, block0, n1, block0 + n1, n2, x2, w);

        // corr = A11^{-1} * w
        MatrixC[] corr = MatrixC.Zero(n1, blockRows, pencilRows);
        f11.SolveLeaf(w, corr);

        // x1 -= corr
        for (int ib = 0; ib < n1; ib++)
        {
            MatrixC x1b = x1[ib];
            MatrixC cb = corr[ib];
            for (int j = 0; j < pencilRows; j++)
                for (int i = 0; i < blockRows; i++)
                    x1b[i, j] -= cb[i, j];
        }

        // Compose xPerm (local ordering)
        MatrixC[] xPermLocal = new MatrixC[blockCount];
        for (int i = 0; i < n1; i++) xPermLocal[i] = x1[i];
        for (int i = 0; i < n2; i++) xPermLocal[n1 + i] = x2[i];

        // Undo local swaps to return solution in original local ordering
        UndoSwapsInPlace(xPermLocal, localSwaps);

        // Copy back into xBlocks
        for (int i = 0; i < blockCount; i++)
            xBlocks[i].CopyFrom(xPermLocal[i]);
    }

    // ----------------------------------------------------------------
    // Extract the subset of global swaps that affect this submatrix [block0, block0+blockCount-1]
    // and convert them into local index pairs (0..blockCount-1).
    // The returned flattened array contains pairs (local_i, local_j) in the same order as the original swaps.
    // ----------------------------------------------------------------
    private static int[] ExtractLocalSwaps(int[] globalSwaps, int block0, int blockCount)
    {
        if (globalSwaps == null || globalSwaps.Length == 0) return Array.Empty<int>();

        var list = new System.Collections.Generic.List<int>();

        for (int p = 0; p < globalSwaps.Length; p += 2)
        {
            int gi = globalSwaps[p];
            int gj = globalSwaps[p + 1];

            // If both swapped indices lie inside this submatrix, convert to local indices and record.
            if (gi >= block0 && gi < block0 + blockCount && gj >= block0 && gj < block0 + blockCount)
            {
                int li = gi - block0;
                int lj = gj - block0;
                list.Add(li);
                list.Add(lj);
            }
            // If only one of the indices lies inside this submatrix, then the swap moves a block
            // in/out of this submatrix. That is a global reordering that this local solver cannot
            // handle by local swaps alone. In our design we only record swaps that affect blocks
            // inside the submatrix; swaps that move blocks across submatrix boundaries must have
            // been applied to the global block matrix (they were), and the local solver will see
            // the permuted blocks in place. Therefore we do not record cross-boundary swaps here.
        }

        return list.ToArray();
    }

    // ----------------------------------------------------------------
    // Block matrix multiplication helper: yBlocks = A(row0:row0+rowCount-1, col0:col0+colCount-1) * xBlocks
    // Each yBlocks[bi] is cleared and then accumulated.
    // ----------------------------------------------------------------
    private static void BlockMatrixMultiplication(BlockC blockM, int rowBlock0, int rowBlockCount,
        int colBlock0, int colBlockCount, MatrixC[] xBlocks, MatrixC[] yBlocks)
    {
        if (rowBlockCount == 0 || colBlockCount == 0)
            throw new ArgumentException("Block dimensions are zero.");
        if (xBlocks.Length != colBlockCount || yBlocks.Length != rowBlockCount)
            throw new ArgumentException("Block array length mismatch.");

        int blockRows = blockM.RowSizes[rowBlock0];
        int pencilRows = xBlocks[0].Cols;

        for (int bi = 0; bi < rowBlockCount; bi++)
        {
            MatrixC yb = yBlocks[bi];
            if (yb.Rows != blockRows || yb.Cols != pencilRows)
                throw new ArgumentException("yBlocks dimension mismatch.");
            yb.Clear();

            int rowIndex = rowBlock0 + bi;

            for (int bj = 0; bj < colBlockCount; bj++)
            {
                MatrixC xb = xBlocks[bj];
                if (bi == 0)
                {
                    if (xb.Rows != blockRows || xb.Cols != pencilRows)
                        throw new ArgumentException("xBlocks dimension mismatch.");
                }

                MatrixC aBlock = blockM.GetBlock(rowIndex, colBlock0 + bj);
                aBlock.GemmIP(xb, yb, Complex.One, Complex.One);
            }
        }
    }
}
