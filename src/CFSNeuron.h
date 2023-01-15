/*
 * CFSNeuron.h
 *
 *  Created on: 19 Feb 2013
 *      Author: kostasl
 */

#ifndef CFSNEURON_H_
#define CFSNEURON_H_

#include "CNeuron.h"
namespace SpikeN {

class CFSNeuron: public CNeuron {
public:
	CFSNeuron(unsigned int ID);
	virtual ~CFSNeuron();
};

};

#endif /* CFSNEURON_H_ */
