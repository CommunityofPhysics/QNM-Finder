// File: Initiator.cs

using MathNet.Numerics;
using Microsoft.CodeAnalysis;
using Microsoft.CodeAnalysis.CSharp;
using System.IO;
using System.Linq.Expressions;
using System.Net.NetworkInformation;
using System.Numerics;
using System.Reflection;
using System.Runtime.ConstrainedExecution;
using System.Xml.Linq;
using CLACFramework;

namespace IOFramework;

public static class Initiator
{
    // ------------------------------------------------------------
    // Simple representation of a single typed parameter
    // ------------------------------------------------------------
    private sealed record Param(string Type, string Name);

    // ------------------------------------------------------------
    // Parse a single typed parameter: "double rho" or "Complex omega"
    // ------------------------------------------------------------
    private static Param ParseSingleParam(string spec)
    {
        string[] parts = spec.Split(' ', StringSplitOptions.RemoveEmptyEntries);
        if (parts.Length != 2)
            throw new Exception($"Invalid free parameter spec: '{spec}' (expected 'type name').");

        return new Param(parts[0], parts[1]);
    }

    // ------------------------------------------------------------
    // Build full signature: free + all bound (as double)
    // e.g. "double rho, double sigma, double gamma, double mu, double m, double k"
    // ------------------------------------------------------------
    private static string BuildFullSignature(Param free, IEnumerable<string> boundNames)
    {
        string boundSig = string.Join(", ", boundNames.Select(n => $"double {n}"));
        return boundSig.Length == 0 ? $"{free.Type} {free.Name}" : $"{free.Type} {free.Name}, {boundSig}";
    }

    // ------------------------------------------------------------
    // Generic Roslyn compiler: className + methodName + signature + expr
    // ------------------------------------------------------------
    private static Delegate CompileDynamic(string className, string methodName, string signature, string expressionBody)
    {
        string code = $@"
                    using System;
                    using System.Numerics;
                    using MathNet.Numerics;

                    public static class {className}
                    {{
                        public static Complex {methodName}({signature}) => {expressionBody};
                    }}";

        var syntaxTree = CSharpSyntaxTree.ParseText(code);

        var refs = BuildMetadataReferences();

        var compilation = CSharpCompilation.Create(assemblyName: $"{className}Assembly", syntaxTrees: new[] { syntaxTree },
            references: refs, options: new CSharpCompilationOptions(OutputKind.DynamicallyLinkedLibrary));

        using var ms = new MemoryStream();
        var result = compilation.Emit(ms);

        if (!result.Success)
        {
            string errors = string.Join("\n", result.Diagnostics);
            throw new Exception($"{className} Roslyn compilation failed:\n" + errors);
        }

        ms.Seek(0, SeekOrigin.Begin);
        Assembly asm = Assembly.Load(ms.ToArray());

        Type type = asm.GetType(className)
            ?? throw new Exception($"Type {className} not found.");

        MethodInfo method = type.GetMethod(methodName)
            ?? throw new Exception($"Method {methodName} not found in {className}.");

        return DynamicDelegate(method);
    }

    // ------------------------------------------------------------
    // Compile UF, VF, WF into Func<double, Complex>
    // freeParamSpec: e.g. "double rho"
    // Physics: all bound parameters (sigma, gamma, mu, m, k, ...)
    // ------------------------------------------------------------
    public static Func<double, Complex>[,] CompileSystemFunctions(string[,] systemFunctions, string freeParam, Dictionary<string, double> physics)
    {
        Param free = ParseSingleParam(freeParam);
        if (!string.Equals(free.Type, "double", StringComparison.OrdinalIgnoreCase))
            throw new Exception($"System free parameter must be 'double', got '{free.Type}'.");

        List<string> boundNames = physics.Keys.OrderBy(n => n).ToList();
        string signature = BuildFullSignature(free, boundNames);

        int rows = systemFunctions.GetLength(0);
        int cols = systemFunctions.GetLength(1);

        Func<double, Complex>[,] systemFuncs = new Func<double, Complex>[rows, cols];

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
            {
                string expr = systemFunctions[i, j];
                if (expr == null) continue;
                expr = expr.Trim();

                string className = $"DynamicSF_{i}_{j}";
                string methodName = $"SF_{i}_{j}";

                Delegate raw = CompileDynamic(className, methodName, signature, expr);
                systemFuncs[i, j] = Wrap<double>(raw, free, boundNames, physics);
            }

        return systemFuncs;
    }

    // ------------------------------------------------------------
    // Compile Sigma into Func<Complex, Complex>
    // freeParamSpec: e.g. "Complex omega"
    // Physics: all bound parameters (sigma, gamma, mu, m, k, ...)
    // ------------------------------------------------------------
    public static Func<Complex, Complex>[][] CompileSigma(string[] customSigma, string freeParam, Dictionary<string, double> physics)
    {
        Param free = ParseSingleParam(freeParam);
        if (!string.Equals(free.Type, "Complex", StringComparison.OrdinalIgnoreCase))
            throw new Exception($"Sigma free parameter must be 'Complex', got '{free.Type}'.");

        List<string> boundNames = physics.Keys.OrderBy(n => n).ToList();
        string signature = BuildFullSignature(free, boundNames);

        // Left Sigma
        Delegate rawLeft = CompileDynamic("DynamicSigmaLeft", "SigmaLeft", signature, customSigma[0]);
        Func<Complex, Complex> SigmaLeft = Wrap<Complex>(rawLeft, free, boundNames, physics);

        // Right Sigma
        Delegate rawRight = CompileDynamic("DynamicSigmaRight", "SigmaRight", signature, customSigma[1]);
        Func<Complex, Complex> SigmaRight = Wrap<Complex>(rawRight, free, boundNames, physics);

        return new[] { new[] { SigmaLeft, Calculus.Derivative(SigmaLeft) }, new[]{ SigmaRight, Calculus.Derivative(SigmaRight) } };
    }

    // ------------------------------------------------------------
    // Generic wrapper: T is the type of the free variable
    // Roslyn method signature: (T free, double p1, double p2, ...)
    // We expose: Func<T, Complex>, binding all doubles from Physics
    // ------------------------------------------------------------
    private static Func<T, Complex> Wrap<T>(Delegate raw, Param free, List<string> boundNames, Dictionary<string, double> physics)
    {
        return (T freeValue) =>
        {
            object[] args = new object[1 + boundNames.Count];

            // First argument: free variable (rho or omega)
            args[0] = freeValue!;

            // Remaining arguments: bound parameters as doubles
            for (int i = 0; i < boundNames.Count; i++)
            {
                string name = boundNames[i];
                args[i + 1] = physics[name];
            }

            return (Complex)raw.DynamicInvoke(args);
        };
    }

    // ------------------------------------------------------------
    // Create a delegate with the correct signature dynamically
    // ------------------------------------------------------------
    private static Delegate DynamicDelegate(MethodInfo method)
    {
        List<Type> paramTypes = method.GetParameters().Select(p => p.ParameterType).ToList();
        paramTypes.Add(method.ReturnType);
        Type delegateType = Expression.GetDelegateType(paramTypes.ToArray());
        return method.CreateDelegate(delegateType);
    }

    // ------------------------------------------------------------
    // Create metadata references
    // ------------------------------------------------------------
    private static List<MetadataReference> BuildMetadataReferences()
    {
        Dictionary<string, MetadataReference> refs = new(StringComparer.OrdinalIgnoreCase);

        foreach (Assembly asm in AppDomain.CurrentDomain.GetAssemblies())
        {
            if (asm.IsDynamic || string.IsNullOrWhiteSpace(asm.Location))
                continue;

            TryAddMetadataReference(asm.Location, refs);
        }

        AddDllReferences(Path.Combine(AppContext.BaseDirectory, "Extensions"), refs);

        string currentDir = Directory.GetCurrentDirectory();

        if (!string.Equals(currentDir, AppContext.BaseDirectory, StringComparison.OrdinalIgnoreCase))
            AddDllReferences(Path.Combine(currentDir, "Extensions"), refs);

        return refs.Values.ToList();
    }

    private static void AddDllReferences(string directory, Dictionary<string, MetadataReference> refs)
    {
        if (!Directory.Exists(directory))
            return;

        foreach (string dll in Directory.EnumerateFiles(directory, "*.dll", SearchOption.AllDirectories).OrderBy(p => p, StringComparer.OrdinalIgnoreCase))
            TryAddMetadataReference(dll, refs);
    }

    private static void TryAddMetadataReference(string dll, Dictionary<string, MetadataReference> refs)
    {
        if (refs.ContainsKey(dll))
            return;

        try
        {
            AssemblyName.GetAssemblyName(dll);

            refs[dll] = MetadataReference.CreateFromFile(dll);
            Assembly.LoadFrom(dll);
        }
        catch
        {
            // Ignore native, incompatible, missing-dependency, or otherwise invalid DLLs.
        }
    }

}
