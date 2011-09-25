/*
 *  dRowAudio_StereoMeter.h
 *  dRowAudio
 *
 *  Created by David Rowland on 27/03/2009.
 *  Copyright 2009 dRowAudio. All rights reserved.
 *
 */

#ifndef _DROWAUDIOSTEREOMETER_H_
#define _DROWAUDIOSTEREOMETER_H_

#include "../core/dRowAudio_StandardHeader.h"

#include "dRowAudio_Meter.h"

class StereoMeter	:	public Component,
						public Timer
{
public:
	
	StereoMeter(float* leftBarValueToUse_, float* leftLineValueToUse_,
				float* rightBarValueToUse_, float* rightLineValueToUse_);
	
	~StereoMeter();
	
	void resized();
	
	void setUpdateInterval(int newUpdateInterval);
	
	void timerCallback();
	
private:
	MeterComponent *leftMeter, *rightMeter;
	float *leftBarValueToUse, *leftLineValueToUse;
	float *rightBarValueToUse, *rightLineValueToUse;
	
	JUCE_LEAK_DETECTOR (StereoMeter);
};

#endif //_DROWAUDIOSTEREOMETER_H_