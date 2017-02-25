/*
 * CNeuron.cpp
 *
 *  Created on: 19 Feb 2013
 *      Author: kostasl
 */

#include "CNeuron.h"
using namespace SpikeN;

CNeuron::CNeuron(unsigned int ID) {
	// TODO Auto-generated constructor stub

	mbFired = false;
	miID = ID;

}

float CNeuron::stepUpdateVm(float I,double ts)
{
	ts = 0.50;
	//Update membrane Voltage According to Fitted function By Izhikevich
	//Do in two steps for numerical stability
	mfv += ts*(0.04*mfv*mfv + 5.0*mfv + 140.0-mfu + I);
	mfv += ts*(0.04*mfv*mfv + 5.0*mfv + 140.0-mfu + I);

	//Membrane Recovery
	mfu += mfa*(mfb*mfv-mfu);

	if (mfv > cFiringThres)
	{
		mbFired = true;
		mfv = mfc;
		mfu += mfd; //Update Membrane Recovery Variable
	}
	else
		mbFired = false;

	return mfv; //Return the voltage
}

bool CNeuron::hasSpiked()
{
	return mbFired;
}

CNeuron::~CNeuron() {
	// TODO Auto-generated destructor stub
}

