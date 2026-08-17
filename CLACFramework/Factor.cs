// File: Factor.cs

using System;
using System.Numerics;
using MathNet.Numerics;

namespace CLACFramework;

public sealed class BlockFactor
{
    private readonly BlockC blocks;
    private readonly int[] externalOrder;
    private readonly int blockCount;
    private readonly int blockRows;

    private readonly LUFactorC? leafLU;
    private readonly BlockFactor? f11;
    private readonly BlockFactor? fS;

    private readonly int[]? leftOrder;
    private readonly int[]? rightOrder;

    private readonly int[] ownSwaps;

    private BlockFactor(BlockC blocks, int[] externalOrder, int blockCount, int blockRows, LUFactorC? leafLU, 
        BlockFactor? f11, BlockFactor? fS, int[]? leftOrder, int[]? rightOrder, int[] ownSwaps)
    {
        this.blocks = blocks;
        this.externalOrder = externalOrder;
        this.blockCount = blockCount;
        this.blockRows = blockRows;
        this.leafLU = leafLU;
        this.f11 = f11;
        this.fS = fS;
        this.leftOrder = leftOrder;
        this.rightOrder = rightOrder;
        this.ownSwaps = ownSwaps ?? Array.Empty<int>();
    }

    public static BlockFactor Factor(BlockC blockM, int block0, int blockCount)
    {
        if (blockM == null)
            throw new ArgumentNullException(nameof(blockM));

        if (blockCount <= 0)
            throw new ArgumentException("blockCount must be positive.", nameof(blockCount));

        if (block0 < 0 || block0 + blockCount > blockM.BlockRows)
            throw new ArgumentOutOfRangeException(nameof(block0), "Requested block range is outside the matrix.");

        int blockRows = blockM.RowSizes[block0];

        for (int k = block0; k < block0 + blockCount; k++)
        {
            if (blockM.RowSizes[k] != blockRows || blockM.ColSizes[k] != blockRows)
                throw new InvalidOperationException("Non-uniform square block size is not supported in BlockFactor.");
        }

        int[] order = new int[blockCount];

        for (int i = 0; i < blockCount; i++)
            order[i] = block0 + i;

        return FactorInternal(blockM, order);
    }

    private static BlockFactor FactorInternal(BlockC blockM, int[] externalOrder)
    {
        if (blockM == null)
            throw new ArgumentNullException(nameof(blockM));

        if (externalOrder == null)
            throw new ArgumentNullException(nameof(externalOrder));

        int blockCount = externalOrder.Length;

        if (blockCount <= 0)
            throw new ArgumentException("externalOrder must be non-empty.", nameof(externalOrder));

        int blockRows = blockM.RowSizes[externalOrder[0]];

        for (int i = 0; i < blockCount; i++)
        {
            int k = externalOrder[i];

            if (blockM.RowSizes[k] != blockRows || blockM.ColSizes[k] != blockRows)
                throw new InvalidOperationException("Non-uniform square block size is not supported in BlockFactor.");
        }

        int[] pivotedOrder = (int[])externalOrder.Clone();

        (int[] ownSwaps, LUFactorC leadingLU) = PivotLeadingBlockIfNeeded(blockM, pivotedOrder);

        if (blockCount == 1)
        {
            return new BlockFactor(blockM, (int[])externalOrder.Clone(), 1, blockRows, leadingLU, null, null, null, null, ownSwaps);
        }

        int n1 = blockCount / 2;
        int n2 = blockCount - n1;

        int[] leftOrder = new int[n1];
        int[] rightOrder = new int[n2];

        Array.Copy(pivotedOrder, 0, leftOrder, 0, n1);
        Array.Copy(pivotedOrder, n1, rightOrder, 0, n2);

        BlockFactor f11 = FactorInternal(blockM, leftOrder);
        BlockC S = BuildSchurBlock(blockM, f11, rightOrder, blockRows);

        int[] sOrder = new int[n2];

        for (int i = 0; i < n2; i++)
            sOrder[i] = i;

        BlockFactor fS = FactorInternal(S, sOrder);

        return new BlockFactor(blockM, (int[])externalOrder.Clone(), blockCount, blockRows, null, f11, fS, leftOrder, rightOrder, ownSwaps);
    }

    private static (int[] swaps, LUFactorC leadingLU) PivotLeadingBlockIfNeeded(BlockC blockM, int[] pivotedOrder)
    {
        if (TryLU(blockM.GetBlock(pivotedOrder[0], pivotedOrder[0]), out LUFactorC? lu) && lu != null)
            return (Array.Empty<int>(), lu);

        int pivot = -1;
        LUFactorC? pivotLU = null;

        for (int i = 1; i < pivotedOrder.Length; i++)
        {
            if (TryLU(blockM.GetBlock(pivotedOrder[i], pivotedOrder[i]), out LUFactorC? candidateLU) && candidateLU != null)
            {
                pivot = i;
                pivotLU = candidateLU;
                break;
            }
        }

        if (pivot < 0 || pivotLU == null)
            throw new InvalidOperationException("No LU-able diagonal block found for block-level pivoting.");

        int tmp = pivotedOrder[0];
        pivotedOrder[0] = pivotedOrder[pivot];
        pivotedOrder[pivot] = tmp;

        return (new[] { 0, pivot }, pivotLU);
    }

    private static bool TryLU(MatrixC matrix, out LUFactorC? lu)
    {
        try
        {
            lu = matrix.LU();
            return true;
        }
        catch (InvalidOperationException)
        {
            lu = null;
            return false;
        }
        catch (ArithmeticException)
        {
            lu = null;
            return false;
        }
    }

    private static void ApplySwapsInPlace(MatrixC[] arr, int[] swaps)
    {
        if (swaps == null || swaps.Length == 0)
            return;

        if ((swaps.Length & 1) != 0)
            throw new InvalidOperationException("Swap array must contain pairs of indices.");

        for (int p = 0; p < swaps.Length; p += 2)
        {
            int i = swaps[p];
            int j = swaps[p + 1];

            if (i < 0 || i >= arr.Length || j < 0 || j >= arr.Length)
                throw new InvalidOperationException("Swap index is outside the array bounds.");

            MatrixC tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
        }
    }

    private static void UndoSwapsInPlace(MatrixC[] arr, int[] swaps)
    {
        if (swaps == null || swaps.Length == 0)
            return;

        if ((swaps.Length & 1) != 0)
            throw new InvalidOperationException("Swap array must contain pairs of indices.");

        for (int p = swaps.Length - 2; p >= 0; p -= 2)
        {
            int i = swaps[p];
            int j = swaps[p + 1];

            if (i < 0 || i >= arr.Length || j < 0 || j >= arr.Length)
                throw new InvalidOperationException("Swap index is outside the array bounds.");

            MatrixC tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
        }
    }

    private static BlockC BuildSchurBlock(BlockC blockM, BlockFactor f11, int[] rightOrder, int blockRows)
    {
        if (blockM == null)
            throw new ArgumentNullException(nameof(blockM));

        if (f11 == null)
            throw new ArgumentNullException(nameof(f11));

        if (rightOrder == null)
            throw new ArgumentNullException(nameof(rightOrder));

        int n1 = f11.blockCount;
        int n2 = rightOrder.Length;

        if (n1 <= 0 || n2 <= 0)
            throw new ArgumentException("Schur complement dimensions must be positive.");

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
            MatrixC[] rhs = new MatrixC[n1];
            MatrixC[] z = MatrixC.Zero(n1, blockRows, blockRows);

            for (int kb = 0; kb < n1; kb++)
                rhs[kb] = blockM.GetBlock(f11.externalOrder[kb], rightOrder[jb]);

            f11.SolveLeaf(rhs, z);

            for (int ib = 0; ib < n2; ib++)
            {
                MatrixC sij = new MatrixC(blockRows, blockRows);

                sij.CopyFrom(blockM.GetBlock(rightOrder[ib], rightOrder[jb]));

                for (int kb = 0; kb < n1; kb++)
                {
                    MatrixC a21 = blockM.GetBlock(rightOrder[ib], f11.externalOrder[kb]);
                    MatrixC zk = z[kb];

                    a21.GemmIP(zk, sij, -Complex.One, Complex.One);
                }

                sBlocks.SetBlock(ib, jb, sij);
            }
        }

        return sBlocks;
    }

    public void SolveLeaf(MatrixC[] bBlocks, MatrixC[] xBlocks)
    {
        if (bBlocks == null)
            throw new ArgumentNullException(nameof(bBlocks));

        if (xBlocks == null)
            throw new ArgumentNullException(nameof(xBlocks));

        if (bBlocks.Length != blockCount || xBlocks.Length != blockCount)
            throw new ArgumentException("Block array length mismatch.");

        if (blockCount <= 0)
            throw new InvalidOperationException("Invalid factorization state.");

        if (bBlocks[0] == null)
            throw new ArgumentNullException(nameof(bBlocks), "bBlocks contains a null matrix.");

        int pencilRows = bBlocks[0].Cols;

        for (int i = 0; i < blockCount; i++)
        {
            if (bBlocks[i] == null)
                throw new ArgumentNullException(nameof(bBlocks), "bBlocks contains a null matrix.");

            if (xBlocks[i] == null)
                throw new ArgumentNullException(nameof(xBlocks), "xBlocks contains a null matrix.");

            if (bBlocks[i].Rows != blockRows || bBlocks[i].Cols != pencilRows)
                throw new ArgumentException("bBlocks dimension mismatch.");

            if (xBlocks[i].Rows != blockRows || xBlocks[i].Cols != pencilRows)
                throw new ArgumentException("xBlocks dimension mismatch.");
        }

        MatrixC[] bLocal = new MatrixC[blockCount];
        Array.Copy(bBlocks, bLocal, blockCount);

        ApplySwapsInPlace(bLocal, ownSwaps);

        if (blockCount == 1)
        {
            if (leafLU == null)
                throw new InvalidOperationException("Leaf factorization is missing.");

            MatrixC[] xLeaf = MatrixC.Zero(1, blockRows, pencilRows);

            leafLU.SolveIP(bLocal[0], xLeaf[0]);

            UndoSwapsInPlace(xLeaf, ownSwaps);

            xBlocks[0].CopyFrom(xLeaf[0]);
            return;
        }

        if (f11 == null || fS == null || leftOrder == null || rightOrder == null)
            throw new InvalidOperationException("Interior factorization is incomplete.");

        int n1 = f11.blockCount;
        int n2 = fS.blockCount;

        MatrixC[] b1 = new MatrixC[n1];
        MatrixC[] b2 = new MatrixC[n2];

        Array.Copy(bLocal, 0, b1, 0, n1);
        Array.Copy(bLocal, n1, b2, 0, n2);

        MatrixC[] y1 = MatrixC.Zero(n1, blockRows, pencilRows);
        f11.SolveLeaf(b1, y1);

        MatrixC[] a21y1 = MatrixC.Zero(n2, blockRows, pencilRows);
        BlockMatrixMultiplication(blocks, rightOrder, leftOrder, y1, a21y1);

        MatrixC[] r2 = MatrixC.Zero(n2, blockRows, pencilRows);

        for (int ib = 0; ib < n2; ib++)
        {
            MatrixC rb = r2[ib];
            MatrixC bb = b2[ib];
            MatrixC ax = a21y1[ib];

            for (int j = 0; j < pencilRows; j++)
            {
                for (int i = 0; i < blockRows; i++)
                    rb[i, j] = bb[i, j] - ax[i, j];
            }
        }

        MatrixC[] x2 = MatrixC.Zero(n2, blockRows, pencilRows);
        fS.SolveLeaf(r2, x2);

        MatrixC[] w = MatrixC.Zero(n1, blockRows, pencilRows);
        BlockMatrixMultiplication(blocks, leftOrder, rightOrder, x2, w);

        MatrixC[] correction = MatrixC.Zero(n1, blockRows, pencilRows);
        f11.SolveLeaf(w, correction);

        MatrixC[] x1 = MatrixC.Zero(n1, blockRows, pencilRows);

        for (int ib = 0; ib < n1; ib++)
        {
            MatrixC dst = x1[ib];
            MatrixC yy = y1[ib];
            MatrixC cc = correction[ib];

            for (int j = 0; j < pencilRows; j++)
            {
                for (int i = 0; i < blockRows; i++)
                    dst[i, j] = yy[i, j] - cc[i, j];
            }
        }

        MatrixC[] xPivoted = new MatrixC[blockCount];

        for (int i = 0; i < n1; i++)
            xPivoted[i] = x1[i];

        for (int i = 0; i < n2; i++)
            xPivoted[n1 + i] = x2[i];

        UndoSwapsInPlace(xPivoted, ownSwaps);

        for (int i = 0; i < blockCount; i++)
            xBlocks[i].CopyFrom(xPivoted[i]);
    }

    private static void BlockMatrixMultiplication(BlockC blockM, int[] rowOrder, int[] colOrder, MatrixC[] xBlocks, MatrixC[] yBlocks)
    {
        if (blockM == null)
            throw new ArgumentNullException(nameof(blockM));

        if (rowOrder == null)
            throw new ArgumentNullException(nameof(rowOrder));

        if (colOrder == null)
            throw new ArgumentNullException(nameof(colOrder));

        if (xBlocks == null)
            throw new ArgumentNullException(nameof(xBlocks));

        if (yBlocks == null)
            throw new ArgumentNullException(nameof(yBlocks));

        if (rowOrder.Length == 0 || colOrder.Length == 0)
            throw new ArgumentException("Block orders must be non-empty.");

        if (xBlocks.Length != colOrder.Length || yBlocks.Length != rowOrder.Length)
            throw new ArgumentException("Block array length mismatch.");

        if (xBlocks[0] == null)
            throw new ArgumentNullException(nameof(xBlocks), "xBlocks contains a null matrix.");

        if (yBlocks[0] == null)
            throw new ArgumentNullException(nameof(yBlocks), "yBlocks contains a null matrix.");

        int pencilRows = xBlocks[0].Cols;

        for (int bi = 0; bi < rowOrder.Length; bi++)
        {
            int row = rowOrder[bi];

            if (row < 0 || row >= blockM.BlockRows)
                throw new ArgumentOutOfRangeException(nameof(rowOrder), "rowOrder contains an invalid block row index.");
        }

        for (int bj = 0; bj < colOrder.Length; bj++)
        {
            int col = colOrder[bj];

            if (col < 0 || col >= blockM.BlockCols)
                throw new ArgumentOutOfRangeException(nameof(colOrder), "colOrder contains an invalid block column index.");
        }

        for (int bj = 0; bj < colOrder.Length; bj++)
        {
            MatrixC xb = xBlocks[bj];

            if (xb == null)
                throw new ArgumentNullException(nameof(xBlocks), "xBlocks contains a null matrix.");

            int expectedRows = blockM.ColSizes[colOrder[bj]];

            if (xb.Rows != expectedRows)
                throw new ArgumentException("xBlocks row dimension mismatch.");

            if (xb.Cols != pencilRows)
                throw new ArgumentException("xBlocks column count mismatch.");
        }

        for (int bi = 0; bi < rowOrder.Length; bi++)
        {
            MatrixC yb = yBlocks[bi];

            if (yb == null)
                throw new ArgumentNullException(nameof(yBlocks), "yBlocks contains a null matrix.");

            int expectedRows = blockM.RowSizes[rowOrder[bi]];

            if (yb.Rows != expectedRows)
                throw new ArgumentException("yBlocks row dimension mismatch.");

            if (yb.Cols != pencilRows)
                throw new ArgumentException("yBlocks column count mismatch.");

            yb.Clear();

            for (int bj = 0; bj < colOrder.Length; bj++)
            {
                MatrixC xb = xBlocks[bj];
                MatrixC aBlock = blockM.GetBlock(rowOrder[bi], colOrder[bj]);

                aBlock.GemmIP(xb, yb, Complex.One, Complex.One);
            }
        }
    }
}