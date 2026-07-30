// File: Program.cs

using System.Numerics;
using MathNet.Numerics;
using MathNet.Numerics.Providers.LinearAlgebra;
using IOFramework;
using CLACFramework;
using QNMFinder;

class Program
{
    static void Main()
    {
        Environment.SetEnvironmentVariable("MKL_DYNAMIC", "FALSE");
        Console.WriteLine($"\nRunning QNM Finder on {Control.MaxDegreeOfParallelism} threads with: \n{LinearAlgebraControl.Provider}");

        // ============================================================
        // Configure the Kernel
        // ============================================================

        Config config = Config.Load();
        Kernel kernel = Kernel.Construct(config);

        // ============================================================
        // Solver Parameters and Settings
        // ============================================================

        // --- Guesser ---
        Complex W = config.W;                 // Evaluation point
        int N = Math.Max(config.N, 1);        // Collocation number
        double A = config.A;                  // Accuracy

        // --- Primer ---
        int B = Math.Max(config.B, 1);        // Beyn pencil columns
        Complex C = config.C;                 // Domain center
        double D = config.D;                  // Domain edge
        double Q = config.Q;                  // Quality factor for guess weighting
        double R = config.R;                  // Randomize Beyn pencil
        int S = Math.Max(config.S, 1);        // Domain subdivision number

        // --- Estimator ---
        int M = Math.Max(config.M, 4);        // Quadrature number

        // --- Seeder ---
        int L = Math.Max(config.L, 1);        // Lapping number
        int E = Math.Max(config.E, 1);        // Estimate count
        double P = config.P;                  // Perturbation for seeding

        // --- Refiner ---
        int I = Math.Max(config.I, 1);        // Newton iteration count

        // --- Feeder ---
        string U = config.U;                  // Whether to ultra-refine parameters

        // --- Superfiner ---
        int J = Math.Max(config.J, 0);        // Maximum refinement depth

        // ---- Settings ----
        BSSpec.Preload = config.Preload;        // Toggle preloading objects to memory
        Logger.ToConsole = config.ConsolePrint; // Toggle console printing of logs

        // ============================================================
        // Build output directory structure
        // ============================================================

        string folderName = string.Join("_", config.FolderParams.Select(p => $"{p}={config.Physics[p]}"));
        string folderPath = Path.Combine("QNM_Repository", folderName);

        string subfolderName = string.Join("_", config.SubfolderParams.Select(p => $"{p}={config.Physics[p]}"));
        string subfolderPath = Path.Combine(folderPath, subfolderName);

        Directory.CreateDirectory(subfolderPath);

        // Timestamp
        string timestamp = DateTime.Now.ToString("yyyy_MM_dd_HH_mm_ss");

        // Filenames
        string logFile = $"{timestamp}_RunLog.txt";
        string configFile = $"{timestamp}_Config.xml";

        // Full paths
        string outputPath = Path.Combine(subfolderPath, logFile);
        string configCopyPath = Path.Combine(subfolderPath, configFile);

        // Initialize logger
        Logger.Init(outputPath);

        // Copy config used for this run
        File.Copy("Config.xml", configCopyPath, overwrite: true);

        // Compute interval numbers
        Intervals interval = new Intervals(N, L);

        // ============================================================
        // Run the Guesser
        // ============================================================

        Logger.WriteBoth($"\nRunning Guesser on a {RadialMap.Mode} map with {interval.Guesser.NP} collocation points...");

        Guesser guesser = new Guesser(kernel);
        List<EigenPair> guesses = guesser.GuessQNM(W, N, A);

        Logger.WriteBoth($"\nGuesser finished running with {guesses.Count} eigenpairs.");

        Logger.WriteLine("\nGuessed frequencies from Guesser:");
        Logger.WriteLine("");
        for (int i = 0; i < guesses.Count; i++)
        {
            EigenPair g = guesses[i];
            Logger.WriteLine($"Mode {i:D2}: Freq = {g.Value.Format("F16")}, Res = {g.Residual.Format("E6")}");
        }

        // ============================================================
        // Run the Primer
        // ============================================================

        if (guesses.Count == 0)
        {
            EigenPair rep = Bridge.GenerateRandomEP(interval.Guesser.N, "Guesser");
            guesses.Add(rep);
        }

        (MatrixC beynPencil, SqDomainC[] domains) = Bridge.Primer(guesses, B, C, D, Q, R, S);

        // ============================================================
        // Run the Estimator
        // ============================================================

        Logger.WriteBoth($"\nRunning Estimator on {domains.Length} domains with {interval.Estimator.NP} collocation points...");
        Logger.WriteBoth($"\nThe region is centered at {C.Format("F4")} with a total area of {(D * D).Format("F4")}.");
        Logger.WriteBoth("");

        Estimator estimator = new Estimator(kernel);
        (List<EigenPair> estimates, List<Complex> singularities) = estimator.EstimateQNM(beynPencil, domains, M);

        if (!Logger.ToConsole) Console.Write($"\r{"",-40}\r\u001b[1A");
        Logger.WriteBoth($"\nEstimator finished running with {estimates.Count} eigenpairs.");

        Logger.WriteLine("\nEstimated frequencies from Estimator:");
        Logger.WriteLine("");
        for (int i = 0; i < estimates.Count; i++)
        {
            EigenPair est = estimates[i];
            Logger.WriteLine($"Mode {i:D2}: Freq = {est.Value.Format("F16")}, Res = {est.Residual.Format("E6")}");
        }

        // ============================================================
        // Run the Seeder
        // ============================================================

        if (estimates.Count == 0)
        {
            EigenPair rep = Bridge.GenerateRandomEP(interval.Estimator.N, "Estimator");
            estimates.Add(rep);
        }

        List<EigenPair> seeds = Bridge.Seeder(estimates, E, L, P);

        // ============================================================
        // Run the Refiner
        // ============================================================

        Logger.WriteBoth($"\nRunning Refiner on {seeds.Count} seeds with {interval.Refiner.NP} collocation points...");

        Refiner refiner = new Refiner(kernel);
        List<EigenPair> refined = refiner.RefineQNM(seeds, I);

        Logger.WriteBoth($"\nRefiner finished running with {refined.Count} eigenpairs.");

        Logger.WriteLine("\nRefined frequencies from Refiner:");
        Logger.WriteLine("");
        for (int i = 0; i < refined.Count; i++)
        {
            EigenPair refEP = refined[i];
            Logger.WriteLine($"Mode {i:D2}: Freq = {refEP.Value.Format("F16")}, Res = {refEP.Residual.Format("E6")}");
        }

        // ============================================================
        // Run the Feeder
        // ============================================================

        if (refined.Count == 0)
        {
            EigenPair rep = Bridge.GenerateRandomEP(interval.Refiner.N, "Refiner");
            refined.Add(rep);
        }

        List<EigenPair> feeds = Bridge.Feeder(refined, U);

        // ============================================================
        // Run the Superfiner
        // ============================================================

        Logger.WriteBoth($"\nRunning Superfiner on {feeds.Count} feeds with {interval.SuperfinerBase.NP} base collocation points...");

        Superfiner superfiner = new Superfiner(kernel);
        (List<EigenPair> superfined, List<EstimatedLimit> polished) = superfiner.SuperfineQNM(feeds, I, J);

        Logger.WriteBoth($"\nSuperfiner finished running with {superfined.Count} eigenpairs.");

        Logger.WriteLine("\nSuperfined frequencies from Superfiner:");
        Logger.WriteLine("");
        for (int i = 0; i < superfined.Count; i++)
        {
            EigenPair supEP = superfined[i];
            Logger.WriteLine($"Mode {i:D2}: Freq = {supEP.Value.Format("F16")}, Res = {supEP.Residual.Format("E6")}");
        }

        Logger.WriteLine("\nPolished frequencies from Polisher:");
        Logger.WriteLine("");
        for (int i = 0; i < polished.Count; i++)
        {
            EstimatedLimit el = polished[i];
            Logger.WriteLine($"Mode {i:D2}: Freq = {el.Value.Format("F16")}, Error = {el.Error.Format("E6")}");
        }

        // Close logger
        Logger.Close();
    }
}