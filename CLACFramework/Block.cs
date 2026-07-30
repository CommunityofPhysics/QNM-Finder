// File: Block.cs

using System;
using System.Numerics;

namespace CLACFramework;

public enum BSType
{
    InMemory, InCache, OnFile
}

public sealed class CacheBlock
{
    public readonly int Rows;
    public readonly int Cols;

    public readonly Func<MatrixC> Loader;

    public CacheBlock(int rows, int cols, Func<MatrixC> loader)
    {
        Rows = rows;
        Cols = cols;
        Loader = loader;
    }
}

public sealed class FileBlock
{
    public readonly int Rows;
    public readonly int Cols;

    public readonly string Path;

    public FileBlock(int rows, int cols, string path)
    {
        Rows = rows;
        Cols = cols;
        Path = path;
    }
}

public sealed class BlockC
{
    public int BlockRows { get; }
    public int BlockCols { get; }

    public int[] RowSizes { get; }
    public int[] ColSizes { get; }

    public int TotalRows { get; }
    public int TotalCols { get; }

    private readonly int[] rowOffset;
    private readonly int[] colOffset;

    private readonly BlockEntry[,] blocks;

    public string? FilePath { get; internal set; }

    private sealed class BlockEntry
    {
        public readonly int Rows;
        public readonly int Cols;

        public MatrixC? InMemory { get; private set; }
        public CacheBlock? InCache { get; private set; }
        public FileBlock? OnFile { get; private set; }

        public BSType StorageType
        {
            get
            {
                if (InMemory != null) return BSType.InMemory;
                if (InCache != null) return BSType.InCache;
                if (OnFile != null) return BSType.OnFile;

                throw new InvalidOperationException("BlockEntry has no data source. Failed to determine memory tier.");
            }
        }

        public BlockEntry(MatrixC MB)
        {
            Rows = MB.Rows;
            Cols = MB.Cols;

            InMemory = MB;
            InCache = null;
            OnFile = null;
        }

        public BlockEntry(CacheBlock CB)
        {
            Rows = CB.Rows;
            Cols = CB.Cols;

            InMemory = null;
            InCache = CB;
            OnFile = null;
        }

        public BlockEntry(FileBlock FB)
        {
            Rows = FB.Rows;
            Cols = FB.Cols;

            InMemory = null;
            InCache = null;
            OnFile = FB;
        }

        public MatrixC Load()
        {
            if (InMemory != null)
                return InMemory;

            if (InCache != null)
            {
                MatrixC M = InCache.Loader();
                if (M.Rows != Rows || M.Cols != Cols)
                    throw new InvalidOperationException($"Cached block size mismatch. Expected {Rows}x{Cols}, got {M.Rows}x{M.Cols}.");
                Memorize(M);
                return M;
            }

            if (OnFile != null)
            {
                MatrixC M = BlockC.Deserialize(OnFile.Path);
                if (M.Rows != Rows || M.Cols != Cols)
                    throw new InvalidOperationException($"File block size mismatch. Expected {Rows}x{Cols}, got {M.Rows}x{M.Cols}.");
                return M;
            }

            throw new InvalidOperationException("BlockEntry has no data source.");
        }

        public void Memorize(MatrixC M)
        {
            InMemory = M;
            InCache = null;
            OnFile = null;
        }
    }

    public BlockC(int[] rowSizes, int[] colSizes)
    {
        BlockRows = rowSizes.Length;
        BlockCols = colSizes.Length;

        RowSizes = rowSizes;
        ColSizes = colSizes;

        TotalRows = 0;
        TotalCols = 0;
        for (int i = 0; i < BlockRows; i++) TotalRows += RowSizes[i];
        for (int j = 0; j < BlockCols; j++) TotalCols += ColSizes[j];

        rowOffset = new int[BlockRows];
        colOffset = new int[BlockCols];

        int acc = 0;
        for (int i = 0; i < BlockRows; i++)
        {
            rowOffset[i] = acc;
            acc += RowSizes[i];
        }

        acc = 0;
        for (int j = 0; j < BlockCols; j++)
        {
            colOffset[j] = acc;
            acc += ColSizes[j];
        }

        blocks = new BlockEntry[BlockRows, BlockCols];
    }

    public void SetBlock(int bi, int bj, MatrixC MB)
    {
        ValidateBlockIndex(bi, bj);
        ValidateBlockSize(bi, bj, MB.Rows, MB.Cols);
        blocks[bi, bj] = new BlockEntry(MB);
    }

    public void SetBlock(int bi, int bj, CacheBlock CB)
    {
        ValidateBlockIndex(bi, bj);
        ValidateBlockSize(bi, bj, CB.Rows, CB.Cols);
        blocks[bi, bj] = new BlockEntry(CB);
    }

    public void SetBlock(int bi, int bj, FileBlock FB)
    {
        ValidateBlockIndex(bi, bj);
        ValidateBlockSize(bi, bj, FB.Rows, FB.Cols);
        blocks[bi, bj] = new BlockEntry(FB);
    }

    public MatrixC GetBlock(int bi, int bj)
    {
        ValidateBlockIndex(bi, bj);
        var entry = blocks[bi, bj] ?? throw new InvalidOperationException($"Block ({bi},{bj}) has not been set.");
        return entry.Load();
    }

    private void ValidateBlockIndex(int bi, int bj)
    {
        if (bi < 0 || bi >= BlockRows)
            throw new IndexOutOfRangeException($"Invalid block row index {bi}");

        if (bj < 0 || bj >= BlockCols)
            throw new IndexOutOfRangeException($"Invalid block col index {bj}");
    }

    private void ValidateBlockSize(int bi, int bj, int rows, int cols)
    {
        int expectedRows = RowSizes[bi];
        int expectedCols = ColSizes[bj];

        if (rows != expectedRows || cols != expectedCols)
            throw new ArgumentException($"Block ({bi},{bj}) must be {expectedRows}x{expectedCols}, but got {rows}x{cols}");
    }

    public void SetElement(int i, int j, Complex value)
    {
        (int bi, int li) = LocateRow(i);
        (int bj, int lj) = LocateCol(j);

        var entry = blocks[bi, bj] ?? throw new InvalidOperationException($"Block ({bi},{bj}) has not been set.");
        var block = entry.Load();

        block[li, lj] = value;

        if (entry.StorageType == BSType.OnFile)
            BlockC.Serialize(block, entry.OnFile!.Path);

        if (entry.StorageType == BSType.InCache)
            entry.Memorize(block);
    }

    public Complex GetElement(int i, int j)
    {
        (int bi, int li) = LocateRow(i);
        (int bj, int lj) = LocateCol(j);

        return GetBlock(bi, bj)[li, lj];
    }

    private (int bi, int li) LocateRow(int i)
    {
        if (i < 0 || i >= TotalRows)
            throw new IndexOutOfRangeException($"Row {i} out of range");

        for (int bi = 0; bi < BlockRows; bi++)
        {
            int start = rowOffset[bi];
            int end = start + RowSizes[bi];
            if (i < end) return (bi, i - start);
        }

        throw new Exception("LocateRow failed unexpectedly.");
    }

    private (int bj, int lj) LocateCol(int j)
    {
        if (j < 0 || j >= TotalCols)
            throw new IndexOutOfRangeException($"Col {j} out of range");

        for (int bj = 0; bj < BlockCols; bj++)
        {
            int start = colOffset[bj];
            int end = start + ColSizes[bj];
            if (j < end) return (bj, j - start);
        }

        throw new Exception("LocateCol failed unexpectedly.");
    }

    public MatrixC Assemble()
    {
        MatrixC M = new MatrixC(TotalRows, TotalCols);

        for (int bi = 0; bi < BlockRows; bi++)
        {
            for (int bj = 0; bj < BlockCols; bj++)
            {
                MatrixC B = GetBlock(bi, bj);

                int ro = rowOffset[bi];
                int co = colOffset[bj];

                for (int i = 0; i < B.Rows; i++)
                {
                    for (int j = 0; j < B.Cols; j++)
                        M[ro + i, co + j] = B[i, j];
                }
            }
        }

        return M;
    }

    public static void Serialize(MatrixC M, string file)
    {
        using var fs = new FileStream(file, FileMode.Create, FileAccess.Write);
        using var bw = new BinaryWriter(fs);

        bw.Write(M.Rows);
        bw.Write(M.Cols);

        for (int i = 0; i < M.Rows; i++)
        {
            for (int j = 0; j < M.Cols; j++)
            {
                Complex z = M[i, j];
                bw.Write(z.Real);
                bw.Write(z.Imaginary);
            }
        }
    }

    public static MatrixC Deserialize(string file)
    {
        using var fs = new FileStream(file, FileMode.Open, FileAccess.Read);
        using var br = new BinaryReader(fs);

        int r = br.ReadInt32();
        int c = br.ReadInt32();

        MatrixC M = new MatrixC(r, c);

        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                double re = br.ReadDouble();
                double im = br.ReadDouble();
                M[i, j] = new Complex(re, im);
            }
        }

        return M;
    }

    public static CacheBlock ToCacheBlock(MatrixC M)
    => new CacheBlock(M.Rows, M.Cols, () => M);

    public static FileBlock ToFileBlock(MatrixC M, string file)
    {
        Serialize(M, file);
        return new FileBlock(M.Rows, M.Cols, file);
    }

    public void Erase()
    {
        if (FilePath != null && Directory.Exists(FilePath))
        {
            try
            {
                Directory.Delete(FilePath, recursive: true);
            }
            catch { }
        }

        for (int bi = 0; bi < BlockRows; bi++)
        {
            for (int bj = 0; bj < BlockCols; bj++)
                blocks[bi, bj] = null!;
        }

        FilePath = null;
    }
}

