// File: Moments.cs

using System;
using System.Numerics;
using CLACFramework;

namespace QNMFinder;

public sealed class Moments
{
    // Moment density functions

    public static Func<Complex, Complex[]> MomentDensities(Func<Complex, MatrixC> T, Func<Complex, MatrixC> RT, int k)
    {
        int kp = k + 1;

        Complex winding = 2.0 * Math.PI * Complex.ImaginaryOne;

        return omega =>
        {
            MatrixC TM = T(omega);
            MatrixC RTM = RT(omega);

            LUFactorC luTM;
            try
            {
                luTM = TM.LU();
            }
            catch
            {
                return ScalarC.NaN(kp);
            }

            MatrixC TMInvRTM = luTM.Solve(RTM);
            Complex trace = TMInvRTM.Trace();

            Complex omegaPow = Complex.One;
            Complex[] MD = new Complex[kp];

            for (int i = 0; i < kp; i++)
            {
                MD[i] = omegaPow * trace / winding;
                omegaPow *= omega;
            }

            return MD;
        };

    }

    public static Func<Complex, Complex> MomentDensity(Func<Complex, MatrixC> T, Func<Complex, MatrixC> RT, int k)
    {
        Func<Complex, Complex[]> MDs = MomentDensities(T, RT, k);
        return omega => MDs(omega)[k];
    }

    // Beyn moment density functions

    public static Func<Complex, MatrixC[]> BeynMomentDensities(Func<Complex, MatrixC> T, MatrixC pencil, int k)
    {
        int kp = k + 1;

        Complex winding = 2 * Math.PI * Complex.ImaginaryOne;

        int Rows = pencil.Rows;
        int Cols = pencil.Cols;

        return omega =>
        {
            MatrixC TM = T(omega);

            LUFactorC luTM;
            try
            {
                luTM = TM.LU();
            }
            catch
            {
                return MatrixC.NaN(kp, Rows, Cols);
            }

            MatrixC Beyn;
            try
            {
                Beyn = luTM.Solve(pencil);
            }
            catch
            {
                return MatrixC.NaN(kp, Rows, Cols);
            }

            MatrixC[] BeynMD = new MatrixC[kp];
            Complex omegaPow = Complex.One;

            for (int i = 0; i < kp; i++)
            {
                BeynMD[i] = omegaPow * Beyn / winding;
                omegaPow *= omega;
            }

            return BeynMD;
        };
    }

    public static Func<Complex, MatrixC> BeynMomentDensity(Func<Complex, MatrixC> T, MatrixC pencil, int k)
    {
        Func<Complex, MatrixC[]> BeynMDs = BeynMomentDensities(T, pencil, k);
        return omega => BeynMDs(omega)[k];
    }

}
