/*
 * CRSNeuron.h
 *
 *  Created on: 19 Feb 2013
 *      Author: kostasl
 */

#ifndef CRSNEURON_H_
#define CRSNEURON_H_

#include "CNeuron.h"
namespace SpikeN {

class CRSNeuron: public CNeuron {
public:
	CRSNeuron(unsigned int ID);
	virtual ~CRSNeuron();
};

};
#endif /* CRSNEURON_H_ */
