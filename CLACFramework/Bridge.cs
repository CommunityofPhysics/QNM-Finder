// File: Bridge.cs

using System.Numerics;
using CLACFramework;
using IOFramework;

namespace QNMFinder;

public readonly record struct Nodes(int N)
{
    public int NP => N + 1;   // nodes
}

public readonly record struct Intervals(int N, int L, int J)
{
    public Nodes Guesser => new Nodes(N);
    public Nodes Estimator => new Nodes(2 * N);
    public Nodes Refiner => new Nodes(2 * L * N);
    public Nodes SuperfinerMax => new Nodes(2 * L * N * (1 << J));
}

public sealed record EigenPair(Complex Value, VectorC Mode, double Residual)
{
    public void Deconstruct(out Complex value, out VectorC mode, out double residual)
    {
        value = Value;
        mode = Mode.Clone();
        residual = Residual;
    }
};

public sealed record EstimatedLimit(Complex Value, double Error);

public readonly record struct Depth(int Current, int Max)
{
    public Depth Next() => new Depth(Current + 1, Max);
}

public readonly record struct SDCount(int BranchPoint, int Spurious, int Oversized);

public static class Bridge
{
    // ------------------------------------------------------------
    // 1. Primer: Bridges Guesser to Estimator
    // ------------------------------------------------------------
    public static (MatrixC, SqDomainC[]) Primer(List<EigenPair> guesses, int B, Complex C, double D, double Q, double R, int S)
    {
        if (guesses.Count == 0)
            throw new Exception("No eigenpair is provided.");

        int G = guesses.Count;
        if (G < B)
        {
            Logger.WriteLine($"\nNot enough guesses to build Beyn's pencil of {B} columns. Reducing to {G}.");
            B = G;
        }

        // Double the resolution
        List<EigenPair> upscaled = new List<EigenPair>(G);
        for (int i = 0; i < G; i++)
        {
            VectorC mode = Discretizer.Interpolate(guesses[i].Mode, 2);
            upscaled.Add(new EigenPair(guesses[i].Value, mode, guesses[i].Residual));
        }

        int NP = upscaled[0].Mode.Size;

        VectorC[] W = new VectorC[B];

        for (int i = 0; i < B; i++)
        {
            VectorC Wi = new VectorC(G);
            double l = (B == 1) ? 0.0 : (double)i * (G - 1) / (B - 1);

            for (int j = 0; j < G; j++)
            {
                if (Q <= 0.0)
                {
                    int center = (int)Math.Round(l);
                    Wi[j] = (j == center ? Complex.One : Complex.Zero);
                }
                else
                {
                    double w = Math.Exp(-(j - l) * (j - l) / (2.0 * Q * Q));
                    Wi[j] = new Complex(w, 0.0);
                }
            }

            Wi.Normalize();
            W[i] = Wi;
        }

        MatrixC Weight = MatrixC.FromVectors(W);
        MatrixC guessSpace = MatrixC.FromVectors(upscaled.Select(g => g.Mode).ToArray());
        MatrixC beynPencil = guessSpace * Weight;

        if (R > 1e-10)
        {
            for (int i = 0; i < B; i++)
            {
                VectorC beynCols = beynPencil.GetColumn(i);

                for (int k = 0; k < NP; k++)
                {
                    double re = R * (2.0 * Random.Shared.NextDouble() - 1.0);
                    double im = R * (2.0 * Random.Shared.NextDouble() - 1.0);

                    beynCols[k] += new Complex(re, im);
                }

                beynCols.Normalize();
                beynPencil.SetColumn(i, beynCols);
            }
        }

        SqDomainC domain = new SqDomainC(C, D);
        SqDomainC[] domains = domain.Subdivide(S);

        return (beynPencil, domains);
    }

    // ------------------------------------------------------------
    // 2. Seeder: Bridges Estimator to Refiner
    // ------------------------------------------------------------
    public static List<EigenPair> Seeder(List<EigenPair> estimates, int E, int L, double P)
    {
        if (estimates.Count == 0)
            throw new Exception("No eigenpair is provided.");

        if (estimates.Count < E)
        {
            Logger.WriteLine($"\nNot enough estimates to build {E} seeds. Reducing to {estimates.Count}.");
            E = estimates.Count;
        }

        List<EigenPair> seeds = new List<EigenPair>(E);

        for (int i = 0; i < E; i++)
        {
            Complex omega = estimates[i].Value;

            VectorC y = Discretizer.Interpolate(estimates[i].Mode, L);
            if (P > 1e-10)
            {
                for (int k = 0; k < y.Size; k++)
                {
                    double re = P * (2.0 * Random.Shared.NextDouble() - 1.0);
                    double im = P * (2.0 * Random.Shared.NextDouble() - 1.0);

                    y[k] += new Complex(re, im);
                }
            }
            y.Normalize();

            double residual = estimates[i].Residual;

            seeds.Add(new EigenPair(omega, y, residual));
        }

        return seeds;
    }

    // ------------------------------------------------------------
    // 3. Feeder: Bridges Refiner to Superfiner
    // ------------------------------------------------------------
    public static List<EigenPair> Feeder(List<EigenPair> refined, string U)
    {
        if (refined.Count == 0)
            throw new Exception("No eigenpair is provided.");

        if (U is null)
            throw new ArgumentNullException(nameof(U), "Polish is null.");

        bool run = false;
        bool auto = false;

        switch (U.Trim().ToLowerInvariant())
        {
            case "automatic" or "auto" or "all":
                run = true;
                auto = true;
                break;

            case "manual" or "few" or "some":
                run = true;
                break;

            case "off" or "turnoff" or "no" or "none":
                break;

            default:
                throw new ArgumentException($"Invalid UltraRefine mode: '{U}'");
        }

        if (run && auto)
        {
            return refined.Select(ep => new EigenPair(ep.Value, ep.Mode, ep.Residual)).ToList();
        }
        else if (run)
        {
            Logger.WriteConsoleInv("\nConsole printing is disabled. Listing refined QNM Frequencies:");
            Logger.WriteConsoleInv("");
            for (int i = 0; i < refined.Count; i++)
            {
                EigenPair refEP = refined[i];
                Logger.WriteConsoleInv($"Mode {i:D2}: Freq = {refEP.Value.Format("F16")}, Res = {refEP.Residual.Format("E6")}");
            }

            Console.Write("\nTo superfine, enter the QNM indices (comma-separated): ");
            string? line = Console.ReadLine();

            string[] tokens = (line ?? "").Trim().ToLowerInvariant().Split(new[] { ',', ' ', '\t', ';' }, StringSplitOptions
                .RemoveEmptyEntries).Select(s => s.Trim()).Where(s => !string.IsNullOrWhiteSpace(s)).ToArray();

            if (tokens.Any(s => s == "automatic" || s == "auto" || s == "all"))
            {
                return refined.Select(ep => new EigenPair(ep.Value, ep.Mode, ep.Residual)).ToList();
            }

            int[] nums = tokens.Where(s => int.TryParse(s, out _)).Select(int.Parse).Distinct().OrderBy(i => i).ToArray();

            List<EigenPair> feeds = new List<EigenPair>(nums.Length);
            List<int> bad = new List<int>();

            foreach (int i in nums)
            {
                if (i >= 0 && i < refined.Count)
                    feeds.Add(refined[i]);
                else
                    bad.Add(i);
            }

            if (bad.Count > 0)
            {
                Console.WriteLine($"\nIndices {string.Join(", ", bad)} are out of range.\n");
            }

            return feeds;
        }

        return new List<EigenPair>(0);
    }

    // ------------------------------------------------------------
    // 4. Random EigenPair Generator
    // ------------------------------------------------------------
    public static EigenPair GenerateRandomEP(int N, string name)
    {
        Logger.WriteConsoleInv("");
        Logger.WriteBoth($"No eigenpair was returned by the {name}. Generating a random eigenpair...");

        int NP = N + 1;

        Random rng = new Random();

        // Random complex eigenvalue
        double re = rng.NextDouble() * 2.0 - 1.0;   // [-1, 1]
        double im = rng.NextDouble() * 2.0 - 1.0;   // [-1, 1]
        Complex omega = new Complex(re, im);

        // Random eigenvector of length NP
        VectorC y = new VectorC(NP);
        for (int k = 0; k < NP; k++)
        {
            double yre = rng.NextDouble() * 2.0 - 1.0;
            double yim = rng.NextDouble() * 2.0 - 1.0;
            y[k] = new Complex(yre, yim);
        }

        // Normalize eigenvector
        y.Normalize();

        // Residual is meaningless for random data → set to +∞
        return new EigenPair(omega, y, double.PositiveInfinity);
    }

}