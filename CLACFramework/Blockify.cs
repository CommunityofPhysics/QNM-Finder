// File: Blockify.cs

using System;
using System.Numerics;
using System.Threading.Tasks;
using System.Xml.Resolvers;

namespace CLACFramework;

public sealed class BSSpec
{
    public static bool Preload { get; set; }

    public bool Quadrate { get; }
    public int Depth { get; }

    private readonly Func<int, int, BSType> Policy;

    public BSSpec(bool quadrate, int depth, Func<int, int, BSType> policy)
    {
        Quadrate = quadrate;
        Depth = depth;

        this.Policy = policy ?? throw new ArgumentNullException(nameof(policy));
    }

    public BSType Get(int bi, int bj)
    {
        if (!Quadrate) return BSType.InMemory;
        else if (!Preload) return BSType.OnFile;
        else return Policy(bi, bj);
    }
}

public static class Blockify
{
    private static string CacheRoot()
    {
        string baseDir = Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData);
        string rootDir = Path.Combine(baseDir, "QNMFinder");
        Directory.CreateDirectory(rootDir);
        return rootDir;
    }

    private static string CreateFolder(string prefix)
    {
        string folder = Path.Combine(CacheRoot(), prefix + "_" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(folder);
        return folder;
    }

    private static BlockC SketchBlockMatrix(int rows, int cols, BSSpec spec, string prefix)
    {
        if (spec == null)
            throw new ArgumentNullException(nameof(spec));

        int depth = spec.Depth;
        bool quadrate = spec.Quadrate;

        if (depth < 0 || depth > 10)
            throw new ArgumentOutOfRangeException(nameof(spec.Depth));

        int interiorBlocks = 1 << depth;

        if (interiorBlocks > rows)
            throw new ArgumentException("Number of interior blocks cannot exceed number of rows.");

        int tailRows = rows % interiorBlocks;
        int interiorRows = rows - tailRows;
        int blockRows = interiorRows / interiorBlocks;

        int rowBlocks = interiorBlocks + (tailRows > 0 ? 1 : 0);
        int colBlocks = quadrate ? rowBlocks : 1;

        int[] rowSizes = new int[rowBlocks];
        int[] colSizes = new int[colBlocks];

        for (int i = 0; i < interiorBlocks; i++)
            rowSizes[i] = blockRows;

        if (tailRows > 0)
            rowSizes[interiorBlocks] = tailRows;

        if (quadrate)
        {
            for (int j = 0; j < interiorBlocks; j++)
                colSizes[j] = blockRows;
            if (tailRows > 0)
                colSizes[interiorBlocks] = tailRows;
        }
        else
        {
            colSizes[0] = cols;
        }

        BlockC BlockM = new BlockC(rowSizes, colSizes);
        BlockM.FilePath = CreateFolder(prefix);

        return BlockM;
    }

    private static void StreamBlockColumn(BlockC BlockM, BSSpec spec, int bj, MatrixC[] MBs)
    {
        int blockRows = BlockM.BlockRows;

        for (int bi = 0; bi < blockRows; bi++)
        {
            MatrixC MB = MBs[bi];
            BSType storage = spec.Get(bi, bj);

            switch (storage)
            {
                case BSType.InMemory:
                    BlockM.SetBlock(bi, bj, MB);
                    break;

                case BSType.InCache:
                    CacheBlock CB = BlockC.ToCacheBlock(MB);
                    BlockM.SetBlock(bi, bj, CB);
                    break;

                case BSType.OnFile:
                    string folder = BlockM.FilePath ?? throw new InvalidOperationException("Missing FilePath.");
                    string file = Path.Combine(folder, $"blk_{bi}_{bj}.bin");
                    FileBlock FB = BlockC.ToFileBlock(MB, file);
                    BlockM.SetBlock(bi, bj, FB);
                    break;
            }

            MBs[bi] = null!;
        }
    }

    private static BlockC BuildBlockColumn(int rows, int cols, BSSpec spec, string prefix,
        Action<BlockC, int, MatrixC[]> fillColumn)
    {
        BlockC BM = SketchBlockMatrix(rows, cols, spec, prefix);

        int blockRows = BM.BlockRows;
        int blockCols = BM.BlockCols;
        int[] rowSizes = BM.RowSizes;
        int[] colSizes = BM.ColSizes;

        for (int bj = 0; bj < blockCols; bj++)
        {
            MatrixC[] MBs = new MatrixC[blockRows];
            for (int bi = 0; bi < blockRows; bi++)
                MBs[bi] = new MatrixC(rowSizes[bi], colSizes[bj]);

            fillColumn(BM, bj, MBs);
            StreamBlockColumn(BM, spec, bj, MBs);
        }

        return BM;
    }

    // ------------------------------------------------------------
    // Vector
    // ------------------------------------------------------------
    public static BlockC Vector(VectorC v, BSSpec spec)
    {
        if (spec.Quadrate)
            throw new ArgumentException("VectorC cannot be quadrate.");

        int rows = v.Size;
        int cols = 1;

        return BuildBlockColumn(rows, cols, spec, "V",
            (BM, bj, MBs) =>
            {
                int blockRows = BM.BlockRows;
                int[] rowSizes = BM.RowSizes;

                for (int bi = 0; bi < blockRows; bi++)
                {
                    MatrixC MB = MBs[bi];

                    int rowStart = 0;
                    for (int k = 0; k < bi; k++)
                        rowStart += rowSizes[k];

                    for (int li = 0; li < rowSizes[bi]; li++)
                    {
                        int globalRow = rowStart + li;
                        MB[li, 0] = v[globalRow];
                    }
                }
            });
    }

    // ------------------------------------------------------------
    // Matrix
    // ------------------------------------------------------------
    public static BlockC Matrix(MatrixC B, BSSpec spec)
    {
        int rows = B.Rows;
        int cols = B.Cols;

        if (spec.Quadrate && rows != cols)
            throw new ArgumentException("Provided MatrixC is not quadrate.");

        return BuildBlockColumn(rows, cols, spec, "M",
            (BM, bj, MBs) =>
            {
                int blockRows = BM.BlockRows;
                int[] rowSizes = BM.RowSizes;
                int[] colSizes = BM.ColSizes;

                int localCols = colSizes[bj];

                Parallel.For(0, localCols, lj =>
                {
                    int globalCol = lj;
                    for (int k = 0; k < bj; k++)
                        globalCol += colSizes[k];

                    for (int bi = 0; bi < blockRows; bi++)
                    {
                        MatrixC MB = MBs[bi];

                        int rowStart = 0;
                        for (int t = 0; t < bi; t++)
                            rowStart += rowSizes[t];

                        for (int li = 0; li < rowSizes[bi]; li++)
                            MB[li, lj] = B[rowStart + li, globalCol];
                    }
                });
            });
    }

    // ------------------------------------------------------------
    // Action
    // ------------------------------------------------------------
    public static BlockC Action(ActionC A, BSSpec spec)
    {
        int rows = A.Rows;
        int cols = A.Cols;

        if (!spec.Quadrate || rows != cols)
            throw new ArgumentException("ActionC must be quadrate.");

        return BuildBlockColumn(rows, cols, spec, "A",
            (BM, bj, MBs) =>
            {
                int blockRows = BM.BlockRows;
                int[] rowSizes = BM.RowSizes;
                int[] colSizes = BM.ColSizes;

                int localCols = colSizes[bj];

                Parallel.For(0, localCols,
                    () => (e: new VectorC(cols), Ae: new VectorC(rows), prev: -1),
                    (lj, loopState, tlScratch) =>
                    {
                        int globalCol = lj;
                        for (int k = 0; k < bj; k++)
                            globalCol += colSizes[k];

                        if (tlScratch.prev >= 0)
                            tlScratch.e[tlScratch.prev] = Complex.Zero;

                        tlScratch.e[globalCol] = Complex.One;
                        tlScratch.prev = globalCol;

                        A.Apply(tlScratch.e, tlScratch.Ae);

                        for (int bi = 0; bi < blockRows; bi++)
                        {
                            MatrixC MB = MBs[bi];

                            int rowStart = 0;
                            for (int l = 0; l < bi; l++)
                                rowStart += rowSizes[l];

                            for (int li = 0; li < rowSizes[bi]; li++)
                                MB[li, lj] = tlScratch.Ae[rowStart + li];
                        }

                        return tlScratch;
                    },

                    tlScratch => { });
            });
    }

}
