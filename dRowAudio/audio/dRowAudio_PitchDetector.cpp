/*
  ==============================================================================

  This file is part of the dRowAudio JUCE module
  Copyright 2004-13 by dRowAudio.

  ------------------------------------------------------------------------------

  dRowAudio is provided under the terms of The MIT License (MIT):

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
  SOFTWARE.

  ==============================================================================
*/

//==============================================================================
PitchDetector::PitchDetector()
    : detectionMethod       (autoCorrelationFunction),
      sampleRate            (44100.0),
      minFrequency          (50), maxFrequency (1600),
      buffer1               (512), buffer2 (512),
      numSamplesNeededForDetection (int ((sampleRate / minFrequency) * 2)),
      currentBlockBuffer    (numSamplesNeededForDetection),
      inputFifoBuffer       (numSamplesNeededForDetection * 2),
      mostRecentPitch       (0.0)
{
    updateFiltersAndBlockSizes();
    logger = Logger::getCurrentLogger();
}

PitchDetector::PitchDetector(int bufferSize)
: detectionMethod       (autoCorrelationFunction),
sampleRate            (44100.0),
minFrequency          (50), maxFrequency (1600),
buffer1               (bufferSize), buffer2 (bufferSize),
numSamplesNeededForDetection (bufferSize),
currentBlockBuffer    (numSamplesNeededForDetection),
inputFifoBuffer       (numSamplesNeededForDetection * 2),
mostRecentPitch       (0.0)
{
    updateFilters();
    logger = Logger::getCurrentLogger();
    
}

PitchDetector::~PitchDetector()
{
}

void PitchDetector::processSamples (const float* samples, int numSamples) noexcept
{
    if (inputFifoBuffer.getNumFree() < numSamples)
        inputFifoBuffer.setSizeKeepingExisting (inputFifoBuffer.getSize() * 2);
    
    inputFifoBuffer.writeSamples (samples, numSamples);
    
    while (inputFifoBuffer.getNumAvailable() >= numSamplesNeededForDetection)
    {
        inputFifoBuffer.readSamples (currentBlockBuffer.getData(), currentBlockBuffer.getSize());
        mostRecentPitch = detectPitchForBlock (currentBlockBuffer.getData(), currentBlockBuffer.getSize());
    }
}

//==============================================================================

double PitchDetector::detectPitchAutoFFT (float* samples, int numSamples) noexcept
{
    //Array<double> pitches;
    //pitches.ensureStorageAllocated (numPitches);
    
    //while (numSamples >= numSamplesNeededForDetection)
    //{
        double pitch = detectAcFftPitchForBlock (samples, numSamples);
        
        if (pitch > 0.0) return pitch;
            //pitches.add (pitch);
            
        //numSamples -= numSamplesNeededForDetection;
        //samples += numSamplesNeededForDetection;
    //}
    /*
    if (pitches.size() == 1)
    {
        return pitches[0];
    }
    else if (pitches.size() > 1)
    {
        DefaultElementComparator<double> sorter;
        pitches.sort (sorter);
        
        const double stdDev = findStandardDeviation (pitches.getRawDataPointer(), pitches.size());
        const double medianSample = findMedian (pitches.getRawDataPointer(), pitches.size());
        const double lowerLimit = medianSample - stdDev;
        const double upperLimit = medianSample + stdDev;
        
        Array<double> correctedPitches;
        correctedPitches.ensureStorageAllocated (pitches.size());
        
        for (int i = 0; i < pitches.size(); ++i)
        {
            const double pitch = pitches.getUnchecked (i);
            
            if (pitch >= lowerLimit && pitch <= upperLimit)
                correctedPitches.add (pitch);
                }
        
        const double finalPitch = findMean (correctedPitches.getRawDataPointer(), correctedPitches.size());
        
        return finalPitch;
    }
    */
    return 0.0;
}


double PitchDetector::detectPitch (float* samples, int numSamples) noexcept
{
    Array<double> pitches;
    pitches.ensureStorageAllocated (int (numSamples / numSamplesNeededForDetection));
    
    while (numSamples >= numSamplesNeededForDetection)
    {
        double pitch = detectPitchForBlock (samples, numSamplesNeededForDetection);//0.0;
        
        if (pitch > 0.0)
            pitches.add (pitch);
        
        numSamples -= numSamplesNeededForDetection;
        samples += numSamplesNeededForDetection;
    }
    
    if (pitches.size() == 1)
    {
        return pitches[0];
    }
    else if (pitches.size() > 1)
    {
        DefaultElementComparator<double> sorter;
        pitches.sort (sorter);
        
        const double stdDev = findStandardDeviation (pitches.getRawDataPointer(), pitches.size());
        const double medianSample = findMedian (pitches.getRawDataPointer(), pitches.size());
        const double lowerLimit = medianSample - stdDev;
        const double upperLimit = medianSample + stdDev;
        
        Array<double> correctedPitches;
        correctedPitches.ensureStorageAllocated (pitches.size());
        
        for (int i = 0; i < pitches.size(); ++i)
        {
            const double pitch = pitches.getUnchecked (i);
            
            if (pitch >= lowerLimit && pitch <= upperLimit)
                correctedPitches.add (pitch);
        }
        
        const double finalPitch = findMean (correctedPitches.getRawDataPointer(), correctedPitches.size());
        
        return finalPitch;
    }
    
    return 0.0;
}

//==============================================================================
void PitchDetector::setSampleRate (double newSampleRate) noexcept
{
    sampleRate = newSampleRate;
    updateFiltersAndBlockSizes();  //updateFilters only if FFT method
}

void PitchDetector::setDetectionMethod (DetectionMethod newMethod)
{
    detectionMethod = newMethod;
}

void PitchDetector::setMinMaxFrequency (float newMinFrequency, float newMaxFrequency) noexcept
{
    minFrequency = newMinFrequency;
    maxFrequency = newMaxFrequency;

    updateFiltersAndBlockSizes(); //updateFilters only if FFT method
}

//==============================================================================
Buffer* PitchDetector::getBuffer (int stageIndex)
{
    switch (stageIndex)
    {
        case 1:     return &buffer1;    break;
        case 2:     return &buffer2;    break;
        default:    return nullptr;
    }
    
    return nullptr;
}

//==============================================================================
void PitchDetector::updateFiltersAndBlockSizes()
{
    lowFilter.setCoefficients (IIRCoefficients::makeLowPass (sampleRate, maxFrequency));
    highFilter.setCoefficients (IIRCoefficients::makeHighPass (sampleRate, minFrequency));
    
    numSamplesNeededForDetection = int (sampleRate / minFrequency) * 2;
    
    inputFifoBuffer.setSizeKeepingExisting (numSamplesNeededForDetection * 2);
    currentBlockBuffer.setSize (numSamplesNeededForDetection);
    
    buffer1.setSizeQuick (numSamplesNeededForDetection);
    buffer2.setSizeQuick (numSamplesNeededForDetection);
}

void PitchDetector::updateFilters()
{
    lowFilter.setCoefficients (IIRCoefficients::makeLowPass (sampleRate, maxFrequency));
    highFilter.setCoefficients (IIRCoefficients::makeHighPass (sampleRate, minFrequency));
}

//==============================================================================
double PitchDetector::detectPitchForBlock (float* samples, int numSamples)
{
    switch (detectionMethod)
    {
        case autoCorrelationFunction:       return detectAcfPitchForBlock (samples, numSamples);
        case squareDifferenceFunction:      return detectSdfPitchForBlock (samples, numSamples);
        //case autoCorrelationFftFunction:    return detectAcFftPitchForBlock (samples, numSamples);
        default:                            return 0.0;
    }
}



double PitchDetector::detectAcFftPitchForBlock (float* samples, int numSamples)
{
    const int minSample = int (sampleRate / maxFrequency);
    const int maxSample = int (sampleRate / minFrequency);
    
    //const int numSamples = buffer1.getSize();
    
    //logger->writeToLog ("numSamples:"+String(numSamples));
    
    const int windowSize = buffer1.getSize();
    
    Buffer magnitudes(windowSize);
    
    lowFilter.reset();
    highFilter.reset();
    lowFilter.processSamples (samples, numSamples);
    highFilter.processSamples (samples, numSamples);
    
    autocorrelateFft (samples, numSamples, buffer1.getData(), magnitudes.getData());
    normalise (buffer1.getData(), buffer1.getSize());
    
    float* bufferData = buffer1.getData();
    
    int firstNegativeZero = 0;
    
    // first peak method
    for (int i = 0; i < numSamples - 1; ++i)
    {
        if (bufferData[i] >= 0.0f && bufferData[i + 1] < 0.0f)
        {
            firstNegativeZero = i;
            break;
        }
    }
    float max = -1.0f;
    int sampleIndex = 0;
    for (int i = jmax (firstNegativeZero, minSample); i < maxSample; ++i)
    {
        if (bufferData[i] > max)
        {
            max = bufferData[i];
            sampleIndex = i;
        }
    }
    if (sampleIndex > 0)
        return sampleRate / sampleIndex;
    else
        return 0.0;
}

//==============================================================================
/** Finds the autocorrelation of a set of given samples using FFT.
 
 This will cross-correlate inputSamples with itself and put the result in
 output samples.
 Resulting FFT magnitudes are saved in magnitudes, to be used in further processing.
 
 FFT the time domain signal
 convert complex FFT output to magnitude and zero phase (i.e. power spectrum)
 take inverse FFT
 
 from http://dsp.stackexchange.com/questions/1919/efficiently-calculating-autocorrelation-using-ffts?rq=1
 
 autocorr = ifft( complex( abs(fft(inputData, n=2*N-1))**2, 0 ) )
 
 or:
 
 %% Cross correlation through a FFT
 n = 1024;
 x = randn(n,1);
 % cross correlation reference
 xref = xcorr(x,x);
 
 %FFT method based on zero padding
 // apply window first!
 fx = fft([x; zeros(n,1)]); % zero pad and FFT
 x2 = ifft(fx.*conj(fx)); % abs()^2 and IFFT
 % circulate to get the peak in the middle and drop one
 % excess zero to get to 2*n-1 samples
 x2 = [x2(n+2:end); x2(1:n)];
 % calculate the error
 d = x2-xref; % difference, this is actually zero
 fprintf('Max error = %6.2f\n',max(abs(d)));
 
 ...
 
 */
template <typename FloatingPointType> void PitchDetector::autocorrelateFft (const FloatingPointType* inputSamples, int numSamples, FloatingPointType* outputSamples, /*int fftSize,*/ FloatingPointType*  magnitudes) noexcept
{
    //int numSamples2 = numSamples*2;
    //    %FFT method based on zero padding
    //    // apply window first!
    
    //Window window(numSamples);
    
    std::copy(inputSamples, inputSamples+(numSamples*sizeof(FloatingPointType)), outputSamples);
    
    //window.applyWindow(outputSamples, numSamples);
    
    // create padding  -  might be a more efficient way of doing this
    //FloatingPointType* inputPadded = static_cast<FloatingPointType *> (malloc(numSamples2*sizeof(FloatingPointType)));
    //zeromem(inputPadded, numSamples2);
    
    //std::copy(outputSamples, outputSamples+(numSamples*sizeof(FloatingPointType)), inputPadded);
    // end create padding
    
    FFT fft(log2(numSamples));

    
    //    fx = fft([x; zeros(n,1)]); % zero pad and FFT
    fft.performFFT(outputSamples); //(inputPadded);
    
    //----------?????????
    
    fft.getMagnitudes(magnitudes);
    
    //    x2 = ifft(fx.*conj(fx)); % abs()^2 and IFFT
    
    // abs, ^2
    for (int i=0; i<numSamples; i++)
    {
        //inputPadded[i] = std::pow (std::abs (inputPadded[i]), 2);
        outputSamples[i] = std::pow (std::abs (outputSamples[i]), 2);
    }
    
    fft.performIFFT(outputSamples); //inputPadded
    
    //std::copy(inputPadded[0], inputPadded[numSamples-1], outputSamples);
    //outputSamples = inputPadded;
    
    // return ifft( complex( abs(fft(inputData, n=2*N-1))**2, 0 ) )
    
    
    //    % circulate to get the peak in the middle and drop one
    //    % excess zero to get to 2*n-1 samples
    //    x2 = [x2(n+2:end); x2(1:n)];
    //    % calculate the error
    //    d = x2-xref; % difference, this is actually zero
    //    fprintf('Max error = %6.2f\n',max(abs(d)));
    
    
    // previous code
    /*
    for (int i = 0; i < numSamples; i++)
    {
        FloatingPointType sum = 0;
        
        for (int j = 0; j < numSamples - i; j++)
            sum += inputSamples[j] * inputSamples[j + i];
            
            outputSamples[i] = sum * (static_cast<FloatingPointType> (1) / numSamples);
            }
     */
}



double PitchDetector::detectAcfPitchForBlock (float* samples, int numSamples)
{
    const int minSample = int (sampleRate / maxFrequency);
    const int maxSample = int (sampleRate / minFrequency);

    lowFilter.reset();
    highFilter.reset();
    lowFilter.processSamples (samples, numSamples);
    highFilter.processSamples (samples, numSamples);
    
    autocorrelate (samples, numSamples, buffer1.getData());
    normalise (buffer1.getData(), buffer1.getSize());

//    float max = 0.0f;
//    int sampleIndex = 0;
//    for (int i = minSample; i < maxSample; ++i)
//    {
//        const float sample = buffer1.getData()[i];
//        if (sample > max)
//        {
//            max = sample;
//            sampleIndex = i;
//        }
//    }

    float* bufferData = buffer1.getData();
//    const int bufferSize = buffer1.getSize();
    int firstNegativeZero = 0;
    
    // first peak method
    for (int i = 0; i < numSamples - 1; ++i)
    {
        if (bufferData[i] >= 0.0f && bufferData[i + 1] < 0.0f)
        {
            firstNegativeZero = i;
            break;
        }
    }
    
    // apply gain ramp
//    float rampDelta = 1.0f / numSamples;
//    float rampLevel = 1.0f;
//    for (int i = 0; i < numSamples - 1; ++i)
//    {
//        bufferData[i] *= cubeNumber (rampLevel);
//        rampLevel -= rampDelta;
//    }
    
    float max = -1.0f;
    int sampleIndex = 0;
    for (int i = jmax (firstNegativeZero, minSample); i < maxSample; ++i)
    {
        if (bufferData[i] > max)
        {
            max = bufferData[i];
            sampleIndex = i;
        }
    }
    
    
//    buffer2.setSizeQuick (numSamples);
/*    autocorrelate (buffer1.getData(), buffer1.getSize(), buffer2.getData());
    normalise (buffer2.getData(), buffer2.getSize());*/
    //buffer2.quickCopy (buffer1.getData(), buffer1.getSize());
//    differentiate (buffer1.getData(), buffer1.getSize(), buffer2.getData());
//    normalise (buffer2.getData()+2, buffer2.getSize()-2);
//    differentiate (buffer2.getData(), buffer2.getSize(), buffer2.getData());
    
/*    for (int i = minSample + 1; i < maxSample - 1; ++i)
    {
        const float previousSample = buffer2.getData()[i - 1];
        const float sample = buffer2.getData()[i];
        const float nextSample = buffer2.getData()[i + 1];

        if (sample > previousSample
            && sample > nextSample
            && sample > 0.5f)
            sampleIndex = i;
    }*/
    
    //differentiate (buffer2.getData(), buffer2.getSize(), buffer2.getData());
    //normalise (buffer2.getData() + minSample, buffer2.getSize() - minSample);
    
//    float min = 0.0f;
//    int sampleIndex = 0;
//    for (int i = minSample; i < maxSample; ++i)
//    {
//        const float sample = buffer2.getData()[i];
//        if (sample < min)
//        {
//            min = sample;
//            sampleIndex = i;
//        }
//    }


    if (sampleIndex > 0)
        return sampleRate / sampleIndex;
    else
        return 0.0;
}

double PitchDetector::detectSdfPitchForBlock (float* samples, int numSamples)
{
    const int minSample = int (sampleRate / maxFrequency);
    const int maxSample = int (sampleRate / minFrequency);
    
    lowFilter.reset();
    highFilter.reset();
    lowFilter.processSamples (samples, numSamples);
    highFilter.processSamples (samples, numSamples);
    
    sdfAutocorrelate (samples, numSamples, buffer1.getData());
    normalise (buffer1.getData(), buffer1.getSize());
    
    // find first minimum that is below a threshold
    const float threshold = 0.25f;
    const float* sdfData = buffer1.getData();
    float min = 1.0f;
    int index = 0;
    
    for (int i = minSample; i < maxSample; ++i)
    {
        const float prevSample = sdfData[i - 1];
        const float sample =  sdfData[i];
        const float nextSample =  sdfData[i + 1];
        
        if (sample < prevSample
            && sample < nextSample
            && sample < threshold)
        {
            if (sample < min)
            {
                min = sample;
                index = i;
            }
//            return sampleRate / i;
//            break;
        }
    }
    
    if (index != 0)
        return sampleRate / index;
    
    return 0.0;
}
