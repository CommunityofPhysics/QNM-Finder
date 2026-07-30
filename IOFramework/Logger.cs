// File: Logger.cs

using System.Numerics;
namespace IOFramework;

public static class Logger
{
    public static volatile bool ToConsole = true;

    private static StreamWriter _writer;
    private static readonly object _lock = new object();

    public static void Init(string path)
    {
        lock (_lock)
        {
            var dir = Path.GetDirectoryName(path);
            if (!string.IsNullOrEmpty(dir))
                Directory.CreateDirectory(dir);

            _writer = new StreamWriter(path, append: false);
            _writer.AutoFlush = true;

            Logger.WriteLine("");
            Logger.WriteLine("============================= Log started ==============================");
            Logger.WriteLine("");
            Logger.WriteLine($"File: {path}");
            Logger.WriteLine($"Timestamp: {DateTime.Now}");
            Logger.WriteLine("");
            Logger.WriteLine("========================================================================");
        }
    }

    public static void WriteLine(string msg)
    {
        lock (_lock)
        {
            if (_writer == null)
                return;

            _writer.WriteLine(msg);

            if (ToConsole)
                Console.WriteLine(msg);
        }
    }

    public static void WriteBoth(string msg)
    {
        lock (_lock)
        {
            if (_writer == null)
                return;

            _writer.WriteLine(msg);
            Console.WriteLine(msg);
        }
    }

    public static void WriteConsole(string msg)
    {
        lock (_lock)
        {
            if (ToConsole)
                Console.WriteLine(msg);
        }
    }

    public static void WriteConsoleInv(string msg)
    {
        lock (_lock)
        {
            if (!ToConsole)
                Console.WriteLine(msg);
        }
    }

    public static void Close()
    {
        lock (_lock)
        {
            if (_writer != null)
            {
                Logger.WriteLine("");
                Logger.WriteLine("============================== Log closed ==============================");
                Logger.WriteLine("");

                _writer.Dispose();
                _writer = null;
            }
        }
    }
}
