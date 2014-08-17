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

#if JUCE_IOS && ! DROWAUDIO_USE_FFTREAL



FFTOperation::FFTOperation (int fftSizeLog2)
    : fftProperties (fftSizeLog2)
{
	fftConfig = vDSP_create_fftsetup (fftProperties.fftSizeLog2, 0);

	fftBuffer.malloc (fftProperties.fftSize);
	fftBufferSplit.realp = fftBuffer.getData();
	fftBufferSplit.imagp = fftBufferSplit.realp + getFFTProperties().fftSizeHalved;	
}

FFTOperation::~FFTOperation()
{
	vDSP_destroy_fftsetup (fftConfig);
}

void FFTOperation::setFFTSizeLog2 (int newFFTSizeLog2)
{
	if (newFFTSizeLog2 != fftProperties.fftSizeLog2)
    {
		vDSP_destroy_fftsetup (fftConfig);
		
		fftProperties.setFFTSizeLog2 (newFFTSizeLog2);
		fftBuffer.malloc (fftProperties.fftSize);
		fftBufferSplit.realp = fftBuffer.getData();
		fftBufferSplit.imagp = fftBufferSplit.realp + getFFTProperties().fftSizeHalved;	
		
		fftConfig = vDSP_create_fftsetup (fftProperties.fftSizeLog2, 0);
	}
}

void FFTOperation::performFFT (float* samples)
{
    //convert to split complex format with evens in real and odds in imag
	vDSP_ctoz ((COMPLEX *) samples, 2, &fftBufferSplit, 1, fftProperties.fftSizeHalved);
    //calc fft
	vDSP_fft_zrip (fftConfig, &fftBufferSplit, 1, fftProperties.fftSizeLog2, FFT_FORWARD);
}

void FFTOperation::performIFFT (float* samples)
{
    /*
    float *real_p = fftBufferSplit.realp, *imag_p = fftBufferSplit.imagp;
    for (i = 0; i < fftProperties.fftSizeHalved; i++) {
        *real_p++ = magnitude[i] * cosf(phase[i]);
        *imag_p++ = magnitude[i] * sinf(phase[i]);
    }
    //from http://pkmital.com/home/2011/04/14/real-fftifft-with-the-accelerate-framework/
    */
    
    vDSP_fft_zrip(fftConfig, &fftBufferSplit, 1, fftProperties.fftSizeLog2, FFT_INVERSE);
    vDSP_ztoc(&fftBufferSplit, 1, (COMPLEX*) samples, 2, fftProperties.fftSizeHalved);
}

//============================================================================



#endif //JUCE_IOS