/*
  ==============================================================================
   This file is port of the JUCE IIR Filter Module.

    27 - 10 -2022
  ==============================================================================
*/

#ifndef IIRFILTER_H
#define IIRFILTER_H

#include <assert.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>






/**
    Structs for IIR filter processing.
*/

    //==============================================================================
    /** A set of coefficients for use in an Filter object.
     * Order is specified in the struct, supports upto 6 coefficients, c4 and c5 are left unititialized if order is 1
     * because the are unused anyways. 
     * so accessing them migth cause errors in audio.

        @tags{DSP}
    */


    typedef struct filterCoefficients
    {
            double c0; 
            double c1; 
            double c2; 
            double c3;    
            double c4;  
            double c5;
            double order;

    } filt;


    typedef union Coefficients
    {
        filt filter;
        double coeffs[7];
    } coeffs;

    
    //==============================================================================

        /** Returns the coefficients for a first order low-pass filter. */
        filt makeFirstOrderLowPass (double sampleRate, double frequency);

        /** Returns the coefficients for a first order high-pass filter. */
        filt makeFirstOrderHighPass (double sampleRate, double frequency);

        /** Returns the coefficients for a first order all-pass filter. */
        filt makeFirstOrderAllPass (double sampleRate, double frequency);

        /** Returns the coefficients for a low-pass filter with variable Q. */
        filt makeLowPass (double sampleRate, double frequency, double Q);

        /** Returns the coefficients for a high-pass filter with variable Q. */
        filt makeHighPass (double sampleRate, double frequency, double Q);

        /** Returns the coefficients for a band-pass filter with variable Q. */
        filt makeBandPass (double sampleRate, double frequency, double Q);

        /** Returns the coefficients for a notch filter with variable Q. */
        filt makeNotch (double sampleRate, double frequency, double Q);

        /** Returns the coefficients for an all-pass filter with variable Q. */
        filt makeAllPass (double sampleRate, double frequency, double Q);

        /** Returns the coefficients for a low-pass shelf filter with variable Q and gain.

            The gain is a scale factor that the low frequencies are multiplied by, so values
            greater than 1.0 will boost the low frequencies, values less than 1.0 will
            attenuate them.
        */
        filt makeLowShelf (double sampleRate, double cutOffFrequency, double Q, double gainFactor);

        /** Returns the coefficients for a high-pass shelf filter with variable Q and gain.

            The gain is a scale factor that the high frequencies are multiplied by, so values
            greater than 1.0 will boost the high frequencies, values less than 1.0 will
            attenuate them.
        */
        filt makeHighShelf (double sampleRate, double cutOffFrequency, double Q, double gainFactor);

        /** Returns the coefficients for a peak filter centred around a
            given frequency, with a variable Q and gain.

            The gain is a scale factor that the centre frequencies are multiplied by, so
            values greater than 1.0 will boost the centre frequencies, values less than
            1.0 will attenuate them.
        */
        filt makePeakFilter (double sampleRate, double centreFrequency, double Q, double gainFactor);

   
        //==============================================================================


        /** Returns the magnitude frequency response of the filter for a given frequency
            and sample rate
        */
        double getMagnitudeForFrequency (double frequency, double sampleRate,  filt Coeffs);

        /** Returns the magnitude frequency response of the filter for a given frequency array
            and sample rate.
        */
        void getMagnitudeForFrequencyArray (double* frequencies, double* magnitudes, int numSamples, double sampleRate, filt Coeffs);

        /** Returns the phase frequency response of the filter for a given frequency and
            sample rate
        */
        double getPhaseForFrequency (double frequency, double sampleRate,  filt Coeffs);

        /** Returns the phase frequency response of the filter for a given frequency array
            and sample rate.
        */
        void getPhaseForFrequencyArray (double* frequencies, double* phases, int numSamples, double sampleRate, filt Coeffs);

    //==============================================================================


#endif
