/*
 * CNeuron.h
 *
 *  Created on: 19 Feb 2013
 *      Author: kostasl
 */

#ifndef CNEURON_H_
#define CNEURON_H_



namespace SpikeN {



class CNeuron {
public:
	CNeuron();
	CNeuron(unsigned int ID);
	unsigned int getID();
	bool hasSpiked();

	float stepUpdateVm(float I,double ts); //Takes the injection Current and updates the Vm
										//Returns I Resets To I=0 if there is no Spike

	virtual ~CNeuron();

	static const float cFiringThres  = 30;//30mv Thres Firing Spikes

protected:
	unsigned int miID;
	double mfa,mfb,mfc,mfd;
	double mfv,mfu;
	bool mbFired;
};

};
#endif /* CNEURON_H_ */
