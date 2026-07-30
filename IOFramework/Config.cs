// File: Config.cs

using CLACFramework;
using System.Reflection;
using System.Xml.Linq;
using System.Numerics;

namespace IOFramework;

public class Config
{
    // ------------------------------------------------------------
    // Radial map specification
    // ------------------------------------------------------------
    public string MapSpec { get; private set; }

    // ------------------------------------------------------------
    // Physics parameters (all bound, all treated as double)
    // ------------------------------------------------------------
    public Dictionary<string, double> Physics { get; private set; }

    // ------------------------------------------------------------
    // System Specifications: Equation coefficients, boundary conditions, custom sigmas
    // ------------------------------------------------------------
    public string SystemFP { get; private set; }            // System free parameters (i.e., "double rho")
    public string[,] SystemEC { get; private set; }        // System equation coefficients [Power, Derivative] (3×3 tensor)

    public string[] SystemBC { get; private set; }         // System boundary conditions [Left, Right]
    public string SigmaFP { get; private set; }             // Sigma free parameter (i.e., "Complex omega")
    public string[] SigmaCM { get; private set; }          // Sigma custom maps [Left, Right]

    // ------------------------------------------------------------
    // Transform parameters
    // ------------------------------------------------------------
    public string TRBlend { get; private set; }             // Transformer global settings
    public string[] TRParams { get; private set; }          // Transform parameters [Left, Right]

    // ------------------------------------------------------------
    // QNMFinder parameters
    // ------------------------------------------------------------

    // --- Guesser ---
    public Complex W { get; private set; }          // Expansion origin
    public int N { get; private set; }              // Collocation number
    public double A { get; private set; }           // Accuracy

    // --- Primer ---
    public int B { get; private set; }              // Beyn pencil columns
    public Complex C { get; private set; }          // Domain center on complex plane
    public double D { get; private set; }           // Domain edge length
    public double Q { get; private set; }           // Quality factor for guess weighting
    public double R { get; private set; }           // Randomize Beyn pencil
    public int S { get; private set; }              // Domain subdivision number

    // --- Estimator ---
    public int M { get; private set; }              // Quadrature number

    // --- Seeder ---
    public int E { get; private set; }              // Estimate count
    public int L { get; private set; }              // Lapping number
    public double P { get; private set; }           // Perturbation for seed eigenvectors

    // --- Refiner ---
    public int I { get; private set; }              // Newton iteration count

    // --- Feeder ---
    public string U { get; private set; }      // Whether to autotune parameters

    // --- Superfiner ---
    public int J { get; private set; }              // Maximum refinement depth

    // ------------------------------------------------------------
    // Settings
    // ------------------------------------------------------------
    public bool Preload { get; private set; }       // Toggle preloading objects to memory
    public bool ConsolePrint { get; private set; }  // Toggle console printing of logs

    // ------------------------------------------------------------
    // Output path
    // ------------------------------------------------------------
    public List<string> FolderParams { get; private set; }
    public List<string> SubfolderParams { get; private set; }

    // ------------------------------------------------------------
    // Load configuration
    // ------------------------------------------------------------
    public static Config Load()
    {
        ValidateConfig();

        XDocument xml = XDocument.Load("Config.xml");
        XElement root = xml.Root;

        Config config = new Config();

        // ------------------------------
        // Radial map specification
        // ------------------------------
         
        config.MapSpec = root.Element("RadialMap").Value.Trim();

        // ------------------------------
        // Physics (all bound parameters)
        // ------------------------------
        config.Physics = root.Element("Parameters").Elements()
            .ToDictionary(x => x.Name.LocalName, x => double.Parse(x.Value));

        // ------------------------------
        // System equation coefficients (ECs)
        // ------------------------------
        XElement sysEC = root.Element("EquationCoefficients");

        config.SystemFP = sysEC.Attribute("free").Value.Trim();

        // Allocate 3×3 tensor
        config.SystemEC = new string[3, 3];

        // Initialize all entries to "0.0"
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
                config.SystemEC[i, j] = "0.0";
        }

        // Loop through <Power i="k"> blocks
        foreach (var coeff in sysEC.Elements())
        {
            XAttribute ai = coeff.Attribute("i");
            XAttribute aj = coeff.Attribute("j");

            if (ai == null || aj == null) continue;

            int i = int.Parse(ai.Value);
            int j = int.Parse(aj.Value);

            config.SystemEC[i, j] = coeff.Value.Trim();
        }

        XElement sysBC = root.Element("BoundaryConditions");

        // ------------------------------
        // System boundary conditions (BCs)
        // ------------------------------
        XElement boundary = sysBC.Element("Boundary");

        config.SystemBC = new string[2];
        config.SystemBC[0] = boundary.Element("Left").Value.Trim();
        config.SystemBC[1] = boundary.Element("Right").Value.Trim();

        // ------------------------------
        // Sigma custom maps (CMs)
        // ------------------------------
        XElement sigma = sysBC.Element("Sigma");

        config.SigmaFP = sigma.Attribute("free").Value.Trim();

        config.SigmaCM = new string[2];
        config.SigmaCM[0] = sigma.Element("Left").Value.Trim();
        config.SigmaCM[1] = sigma.Element("Right").Value.Trim();

        // ------------------------------
        // Transform parameters
        // ------------------------------
        XElement transform = root.Element("Transformer");

        config.TRBlend = transform.Element("Blend").Value.Trim();

        XElement trParams = transform.Element("Transform");

        config.TRParams = new string[2];
        config.TRParams[0] = trParams.Element("Left").Value.Trim();
        config.TRParams[1] = trParams.Element("Right").Value.Trim();

        // ------------------------------
        // Finder parameters
        // ------------------------------
        XElement qnmfinder = root.Element("QNMFinder");

        // --- Guesser ---
        XElement guesser = qnmfinder.Element("Guesser");
        config.W = ScalarC.Parse(guesser.Element("W").Value);
        config.N = (int)guesser.Element("N");
        config.A = (double)guesser.Element("A");

        // --- Primer ---
        XElement primer = qnmfinder.Element("Primer");
        config.B = (int)primer.Element("B");
        config.C = ScalarC.Parse(primer.Element("C").Value);
        config.D = (double)primer.Element("D");
        config.Q = (double)primer.Element("Q");
        config.R = (double)primer.Element("R");
        config.S = (int)primer.Element("S");
        // --- Estimator ---
        XElement estimator = qnmfinder.Element("Estimator");
        config.M = (int)estimator.Element("M");

        // --- Seeder ---
        XElement seeder = qnmfinder.Element("Seeder");
        config.E = (int)seeder.Element("E");
        config.L = (int)seeder.Element("L");
        config.P = (double)seeder.Element("P");

        // --- Refiner ---
        XElement refiner = qnmfinder.Element("Refiner");
        config.I = (int)refiner.Element("I");

        // --- Feeder ---
        XElement feeder = qnmfinder.Element("Feeder");
        config.U = feeder.Element("U").Value.Trim();

        // --- Superfiner ---
        XElement superfiner = qnmfinder.Element("Superfiner");
        config.J = (int)superfiner.Element("J");

        // ------------------------------
        // Settings
        // ------------------------------
        XElement settings = root.Element("Settings");
        config.Preload = bool.Parse(settings.Element("Preload").Value);
        config.ConsolePrint = bool.Parse(settings.Element("ConsolePrint").Value);

        // ------------------------------
        // Output path
        // ------------------------------
        XElement output = root.Element("OutputPath");

        config.FolderParams = output.Element("Folder")
            .Value.Split(',').Select(p => p.Trim()).ToList();

        config.SubfolderParams = output.Element("Subfolder")
            .Value.Split(',').Select(p => p.Trim()).ToList();

        return config;
    }

    // ------------------------------------------------------------
    // Validate config file existence
    // ------------------------------------------------------------
    private static void ValidateConfig()
    {
        bool configMissing = !File.Exists("Config.xml");
        bool autoRunMissing = !File.Exists("AutoRun.bat");

        if (configMissing)
        {
            Console.Write("\nConfig.xml is not found. Running with default configuration...");

            string resourceName = "QNMFinder.Resources.DefaultConfig.xml";
            string xml = LoadResource(resourceName);
            File.WriteAllText("Config.xml", xml);
        }

        if (autoRunMissing)
        {
            Console.Write("\nAutoRun.bat is not found. Creating default AutoRun.bat...");

            string resourceName = "QNMFinder.Resources.AutoRun.bat";
            string bat = LoadResource(resourceName);
            File.WriteAllText("AutoRun.bat", bat);
        }

        if(configMissing || autoRunMissing) Console.WriteLine("");
    }

    private static string LoadResource(string resourceName)
    {
        Assembly asm = Assembly.GetExecutingAssembly();
        using Stream stream = asm.GetManifestResourceStream(resourceName)
            ?? throw new Exception($"Embedded resource not found: {resourceName}");

        using StreamReader reader = new StreamReader(stream);
        return reader.ReadToEnd();
    }
}
