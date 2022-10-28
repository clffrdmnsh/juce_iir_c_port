/*
  ==============================================================================
   This file is port of the JUCE IIR Filter Module.

    27 - 10 -2022
  ==============================================================================
*/

#include "IIRFilter.h"

double inverseRootTwo = 0.70710678118654752440;

// returns the minimum of two numbers
double min(double a, double b) {
    return a < b ? a : b;
}

// returns the maximum of two numbers
double max(double a, double b) {
    return a > b ? a : b;
}

double minus_twopi = -2*M_PI;

filt makeFirstOrderLowPass (double sampleRate, double frequency)
{
    filt filter;
    filter.order = 1;
    assert (sampleRate > 0.0);
    assert (frequency > 0 && frequency <= (sampleRate * 0.5));

    double n = tan (M_PI * frequency / sampleRate);

    filter.c0 = n;
    filter.c1 = n;
    filter.c2 = n+1;
    filter.c3 = n-1;

    return filter;
}

filt makeFirstOrderHighPass (double sampleRate, double frequency)
{
    filt filter;
    filter.order = 1;
    assert (sampleRate > 0.0);
    assert (frequency > 0 && frequency <= (sampleRate * 0.5));

    double n = tan (M_PI * frequency / sampleRate);

    filter.c0 = 1;
    filter.c1 = -1;
    filter.c2 = n+1;
    filter.c3 = n-1;

    return filter;
}

filt makeFirstOrderAllPass (double sampleRate, double frequency)
{
    filt filter;
    filter.order = 1;
    assert (sampleRate > 0.0);
    assert (frequency > 0 && frequency <= (sampleRate * 0.5));

    double n = tan (M_PI * frequency / sampleRate);

    filter.c0 = n - 1;
    filter.c1 = n + 1;
    filter.c1 = n + 1;
    filter.c2 = n - 1;

    return filter;
}

filt makeLowPass (double sampleRate, double frequency, double Q)
{
    filt filter;
    filter.order = 2;
    assert (sampleRate > 0.0);
    assert (frequency > 0 && frequency <= (sampleRate * 0.5));
    assert (Q > 0.0);

    double n = 1 / tan (M_PI * frequency / sampleRate);
    double nSquared = n * n;
    double invQ = 1 / Q;
    double c1 = 1 / (1 + invQ * n + nSquared);

    filter.c0 = c1;
    filter.c1 = c1 * 2;
    filter.c2 = c1;
    filter.c3 = 1;
    filter.c4 = c1 * 2 * (1 - nSquared);
    filter.c5 = c1 * (1 - invQ * n + nSquared);

    return filter;
}

filt makeHighPass (double sampleRate, double frequency, double Q)
{
    filt filter;
    filter.order = 2;
    assert (sampleRate > 0.0);
    assert (frequency > 0 && frequency <= (sampleRate * 0.5));
    assert (Q > 0.0);

    double n = tan (M_PI * frequency / sampleRate);
    double nSquared = n * n;
    double invQ = 1 / Q;
    double c1 = 1 / (1 + invQ * n + nSquared);

    filter.c0 = c1;
    filter.c1 = c1 * -2;
    filter.c2 = c1;
    filter.c3 = 1;
    filter.c4 = c1 * 2 * (nSquared - 1);
    filter.c5 = c1 * (1 - invQ * n + nSquared);

    return filter;
}

filt makeBandPass (double sampleRate, double frequency, double Q)
{
    filt filter;
    filter.order = 2;
    assert (sampleRate > 0.0);
    assert (frequency > 0 && frequency <= (sampleRate * 0.5));
    assert (Q > 0.0);

    double n = 1 / tan (M_PI * frequency / sampleRate);
    double nSquared = n * n;
    double invQ = 1 / Q;
    double c1 = 1 / (1 + invQ * n + nSquared);

    filter.c0 = c1 * n * invQ;
    filter.c1 = 0;
    filter.c2 = -c1 * n * invQ;
    filter.c3 = 1;
    filter.c4 = c1 * 2 * (1 - nSquared);
    filter.c5 = c1 * (1 - invQ * n + nSquared);

    return filter;
}

filt makeNotch (double sampleRate, double frequency, double Q)
{
    filt filter;
    filter.order = 2;
    assert (sampleRate > 0.0);
    assert (frequency > 0 && frequency <= (sampleRate * 0.5));
    assert (Q > 0.0);

    double n = 1 / tan (M_PI * frequency / sampleRate);
    double nSquared = n * n;
    double invQ = 1 / Q;
    double c1 = 1 / (1 + n * invQ + nSquared);
    double b0 = c1 * (1 + nSquared);
    double b1 = 2 * c1 * (1 - nSquared);

    filter.c0 = b0;
    filter.c1 = b1;
    filter.c2 = b0;
    filter.c3 = 1;
    filter.c4 = b1;
    filter.c5 = c1 * (1 - n * invQ + nSquared);

    return filter;
}

filt makeAllPass (double sampleRate, double frequency, double Q)
{
    filt filter;
    filter.order = 2;
    assert (sampleRate > 0);
    assert (frequency > 0 && frequency <= sampleRate * 0.5);
    assert (Q > 0);

    double n = 1 / tan (M_PI * frequency / sampleRate);
    double nSquared = n * n;
    double invQ = 1 / Q;
    double c1 = 1 / (1 + invQ * n + nSquared);
    double b0 = c1 * (1 - n * invQ + nSquared);
    double b1 = c1 * 2 * (1 - nSquared);

    filter.c0 = b0;
    filter.c1 = b1;
    filter.c2 = 1;
    filter.c3 = 1;
    filter.c4 = b1;
    filter.c5 = b0;

    return filter;
}

filt makeLowShelf (double sampleRate,double cutOffFrequency, double Q, double gainFactor)
{
    filt filter;
    filter.order = 2;
    assert (sampleRate > 0.0);
    assert (cutOffFrequency > 0.0 && cutOffFrequency <= sampleRate * 0.5);
    assert (Q > 0.0);
    gainFactor = pow(10, fabs(gainFactor) / 20);
    double A = max (0.0, sqrt (gainFactor));
    double aminus1 = A - 1;
    double aplus1 = A + 1;
    double omega = (2 * M_PI * max (cutOffFrequency, 2.0)) / sampleRate;
    double coso = cos (omega);
    double beta = sin (omega) * sqrt (A) / Q;
    double aminus1TimesCoso = aminus1 * coso;

    filter.c0 = A * (aplus1 - aminus1TimesCoso + beta);
    filter.c1 = A * 2 * (aminus1 - aplus1 * coso);
    filter.c2 = A * (aplus1 - aminus1TimesCoso - beta);
    filter.c3 = aplus1 + aminus1TimesCoso + beta;
    filter.c4 = -2 * (aminus1 + aplus1 * coso);
    filter.c5 = aplus1 + aminus1TimesCoso - beta;
    
    return filter;
}

filt makeHighShelf (double sampleRate, double cutOffFrequency, double Q, double gainFactor)
{
    filt filter;
    filter.order = 2;
    assert (sampleRate > 0);
    assert (cutOffFrequency > 0 && cutOffFrequency <=  (sampleRate * 0.5));
    assert (Q > 0);
    gainFactor = pow(10, fabs(gainFactor) / 20);
    double A = max (0.0, sqrt (gainFactor));
    double aminus1 = A - 1;
    double aplus1 = A + 1;
    double omega = (2 * M_PI * max (cutOffFrequency, 2.0)) / sampleRate;
    double coso = cos (omega);
    double beta = sin (omega) * sqrt (A) / Q;
    double aminus1TimesCoso = aminus1 * coso;

    filter.c0 = A * (aplus1 + aminus1TimesCoso + beta);
    filter.c1 = A * -2 * (aminus1 + aplus1 * coso);
    filter.c2 = A * (aplus1 + aminus1TimesCoso - beta);
    filter.c3 = aplus1 - aminus1TimesCoso + beta;
    filter.c4 = 2 * (aminus1 - aplus1 * coso);
    filter.c5 = aplus1 - aminus1TimesCoso - beta ;
    
    return filter;
}

filt makePeakFilter (double sampleRate, double frequency, double Q, double gainFactor)
{
    filt filter;
    filter.order = 2;
    assert (sampleRate > 0);
    assert (frequency > 0 && frequency <= (sampleRate * 0.5));
    assert (Q > 0);
    assert (gainFactor > 0);
    gainFactor = pow(10, fabs(gainFactor) / 20);
    double A = max (0.0, sqrt (gainFactor));
    double omega = (2 * M_PI * max (frequency, 2.0) / sampleRate);
    double alpha = sin (omega) / (Q * 2);
    double c2 = -2 * cos (omega);
    double alphaTimesA = alpha * A;
    double alphaOverA = alpha / A;

    filter.c0 = 1 + alphaTimesA;
    filter.c1 = c2;
    filter.c2 = 1 - alphaTimesA;
    filter.c3 = 1 + alphaOverA;
    filter.c4 = c2;
    filter.c5 = 1 - alphaOverA;
    
    return filter;
}

//pass &coeffs to the pointer argument.
double getMagnitudeForFrequency (double frequency, double sampleRate, filt Coeffs) 
{
    double complex j = 0 + 1 * I;
    coeffs coefficients;
    coefficients.filter = Coeffs;
    int order = coefficients.coeffs[6];

    assert (frequency >= 0 && frequency <= sampleRate * 0.5);

    double complex numerator = 0;
    double complex denominator = 0;
    double complex factor = 1;

    double complex jw = cexp ((minus_twopi * frequency * j) / sampleRate);

    for (int n = 0; n <= order; ++n)
    {
        numerator += coefficients.coeffs[n] * factor;
        factor *= jw;
    }

    denominator = 1.0;
    factor = jw;

    for (int n = order + 1; n <= 2 * order; ++n)
    {
        denominator += coefficients.coeffs[n] * factor;
        factor *= jw;
    }

    return cabs(numerator / denominator);
}

void getMagnitudeForFrequencyArray (double* frequencies, double* magnitudes, int numSamples, double sampleRate, filt Coeffs)
{
    double complex  j = 0 + 1*I;
    //create an instance of the coefficients union.
    coeffs coefficients;
    //assign the input argument Coeffs structure to the filter section in the coefficients union.
    coefficients.filter = Coeffs;
    
    int order = coefficients.coeffs[6];

    assert (order >= 0);
    
    for (size_t i = 0; i < numSamples; ++i)
    {
        double fq = *(frequencies+i);
        assert (fq >= 0 && fq <= sampleRate * 0.5);

        double complex numerator = 0.0, denominator = 0.0, factor = 1.0;
        double complex jw = cexp (minus_twopi * fq * j / sampleRate);

        for (int n = 0; n <= order; ++n)
        {
            numerator += coefficients.coeffs[n] * factor;
            factor *= jw;
        }

        denominator = 1.0;
        factor = jw;

        for (int n = order + 1; n <= 2 * order; ++n)
        {
            denominator += coefficients.coeffs[n] * factor;
            factor *= jw;
        }

        magnitudes[i] = cabs(numerator / denominator);
    }
}

double getPhaseForFrequency (double frequency, double sampleRate, filt Coeffs)
{
    double complex j = 0 + 1*I;
    coeffs coefficients;
    coefficients.filter = Coeffs;
    int order = coefficients.coeffs[6];

    assert (frequency >= 0 && frequency <= sampleRate * 0.5);

    double complex numerator = 0.0, denominator = 0.0, factor = 1.0;
    double complex jw = cexp ((minus_twopi * frequency * j) / sampleRate);

    for (int n = 0; n <= order; ++n)
    {
        numerator += coefficients.coeffs[n] * factor;
        factor *= jw;
    }

    denominator = 1.0;
    factor = jw;

    for (int n = order + 1; n <= 2 * order; ++n)
    {
        denominator += coefficients.coeffs[n] * factor;
        factor *= jw;
    }

    return carg(numerator / denominator);
}

void getPhaseForFrequencyArray (double* frequencies, double* phases, int numSamples, double sampleRate, filt Coeffs)
{
    assert (sampleRate > 0);
    
    double complex j = 0 + 1*I;
    coeffs coefficients;
    coefficients.filter = Coeffs;
    int order = coefficients.coeffs[6];
    
    assert (order >= 0);

    for (int i = 0; i < numSamples; ++i)
    {
        double fq = *(frequencies+i);
        
        assert (fq >= 0 && fq <= sampleRate * 0.5);

        double complex numerator = 0.0, denominator = 0.0, factor = 1.0;

        double complex jw = cexp (minus_twopi * fq * j / sampleRate);
        
        for (int n = 0; n <= order; ++n)
        {
            numerator += coefficients.coeffs[n] * factor;
            factor *= jw;
        }

        denominator = 1.0;
        factor = jw;

        for (int n = order + 1; n <= 2 * order; ++n)
        {
            denominator += coefficients.coeffs[n] * factor;
            factor *= jw;
        }
        
        phases[i] = carg(numerator / denominator);

    }
}
