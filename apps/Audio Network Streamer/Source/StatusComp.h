/*
  ==============================================================================

    StatusComp.h
    Created: 11 Jun 2011 11:47:52am
    Author:  David Rowland

  ==============================================================================
*/

#ifndef __STATUSCOMP_H_BD9DB0B2__
#define __STATUSCOMP_H_BD9DB0B2__

#include "../JuceLibraryCode/JuceHeader.h"

class StatusComponent : public Component,
						public ValueTree::Listener
{
public:
	
	enum Status {
		disconnected = 0,
		waiting,
		connectedSender,
		connectedReciever,
	};
	
	StatusComponent();
	
	~StatusComponent();
	
	void resized();
	
	void setStatus(Status newStatus);
	
	//==============================================================================
	void valueTreePropertyChanged (ValueTree& treeWhosePropertyHasChanged, const Identifier& property);
	void valueTreeChildAdded (ValueTree& parentTree, ValueTree& childWhichHasBeenAdded) {}
	void valueTreeChildRemoved (ValueTree& parentTree, ValueTree& childWhichHasBeenRemoved) {}
	void valueTreeChildOrderChanged (ValueTree& parentTreeWhoseChildrenHaveMoved) {}
	void valueTreeParentChanged (ValueTree& treeWhoseParentHasChanged) {}
	
	//==============================================================================

private:
	
	ValueTree tree;
	
	Status currentStatus;
	TextButton statusButton;
	Label statusLabel;
};


#endif  // __STATUSCOMP_H_BD9DB0B2__
