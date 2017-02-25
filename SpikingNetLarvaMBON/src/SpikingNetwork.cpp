//============================================================================
// Name        : SpikingNetwork.cpp
// Author      : konstantinos Lagogiannis
// Version     :
// Copyright   : Izhikevich spiking network as in 2003 Paper
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream> //File streams
#include <vector> //for Vectors
#include <math.h>
#include <unistd.h> //For Working Dir
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//Custom objs
#include "CNeuron.h"
#include "CFSNeuron.h"
#include "CRSNeuron.h"

#define PATH_MAX 255

using namespace std;

//GSL VARIABLES
const gsl_rng_type * T;
gsl_rng * r;
time_t t;


int main() {
	const float ts 						= 0.5f; //Timestep
	const unsigned int iTotalSimSteps 	= 1000;
	unsigned int itStep 				= 1;


	unsigned int inetExh				= 800;
	unsigned int inetInh				= 200;
	unsigned int inetSize 				= inetInh + inetExh;

	float* I =  new float[inetSize];//This timestep Injection current input -
	float* F =  new float[inetSize];//Next timestep Injection current input -

	memset(I,0,inetSize);
	memset(F,0,inetSize);

	float** W = new float*[inetSize]; // Vector of pointers to vectors - unitialized Weight Matrix
	float fVm; //Temp Var. Holding Membrane Voltage

	vector<SpikeN::CNeuron*> vNeurons;
	vNeurons.reserve(inetSize);

	cout << "\t SPIKING NEURAL NETWORK N: " << inetSize << endl;


	// Initialise GSL random number generator
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r,(unsigned long) time(&t));


	//Random Initialization of neuron parameters in the network
	//Weight Init Table
	for(register unsigned int i=0;i<inetSize;i++)
	{
		W[i] = new float[inetSize];
		//init Random Weights Exh
		for(register unsigned int j=0;j<inetSize;j++)
		{
			double rnd  = gsl_rng_uniform(r);
			if (j > inetExh)
				W[i][j] = -rnd; ///Inhibitory synapse
			else
				W[i][j] = 0.5*rnd; //Exhitatory

		//	cout << W[i][j];
		}
		//cout << endl;

		 //Initialize Neuron Objects
		SpikeN::CNeuron* objN;
		if (i > inetExh)
		{ //Inhibitory Should be FS cells
			objN = new SpikeN::CFSNeuron(i);
		}
		else
		{ //Exhitatory Should be RS cells
			objN = new SpikeN::CRSNeuron(i);
		}
		//cout << i << endl;
		vNeurons.push_back(objN); //Add To neuron vector

	}//End Of Init Loop

	cout << "Created N: " << vNeurons.size() << endl;

	//Open Output File
	char FilePath[PATH_MAX];
	getcwd(FilePath, PATH_MAX); // reads the current working directory into the array FilePath

	char buff[] = "SpikeRaster.dat";
	char buff2[] = "RSNeuronV.dat";
	char buff3[] = "FSNeuronV.dat";

	ofstream ofRaster(buff, ios::out );
	ofRaster << "#NeuronID\tSpikeTime" << endl;

	ofstream ofRSV(buff2, ios::out );
	ofRSV << "#t\tVm" << endl;

	ofstream ofFSV(buff3, ios::out );
	ofFSV << "#t\tVm" << endl;


	//init Thalamic Input Before 1st Timestep
	for(register unsigned int i=0;i<inetSize;i++)
	{
		double rnd = gsl_ran_gaussian(r,2.0) + 1.0; //Gaussian sigma 2 Mean 1
		if (i > inetExh)
			F[i] = 2.0*rnd; //to inhibitory neuron
		else
			F[i] = 5.0*rnd; //To Exhitatory Neuron
	}//END OF THALAMIC

	//On every timestep
	while(itStep < iTotalSimSteps)
	{
		//At Each Time Step
		///Update I with Spikes from F
		memcpy(I,F,inetSize*sizeof(float));

		//Reset F;
		memset(F,0,inetSize*sizeof(float));


		//For Each Neuron
		unsigned int i = 0;
		for (vector<SpikeN::CNeuron*>::iterator it = vNeurons.begin();it!=vNeurons.end();++it)
		{

			//Add Thalamic Input To next iteration
			double rnd = gsl_ran_gaussian(r,2.0) + 1.0; //Gaussian sigma 2 Mean 1
			if (i >= inetExh)
				F[i] = 2.0*rnd; //to inhibitory neuron
			else
				F[i] = 5.0*rnd; //To Exhitatory Neuron

			//Propagate The Injection Currents from External and Internal Spikes From the Previous Iteration
			//StepUpdateVm(I[i]) //Resets To I=0 if there is no Spike
			fVm = (*it)->stepUpdateVm(I[i],ts);

			//Record Vm RS Cell
			if (i == 0)
				ofRSV << itStep << "\t" << fVm << endl;

			if (i == (inetSize-1)) //Switched to inhibitory Cells
				ofFSV << itStep << "\t" << fVm << endl;

			if ((*it)->hasSpiked())
			{
				///Update all F injections from neuron i to neuron j
				for(register unsigned int j=0;j<inetSize;j++)
					F[j] += W[j][i];

				ofRaster << i << "\t" << itStep << endl;
			//	cout << i << "\t" << itStep << endl;
			}

			//Check if neuron Fired - Save to List of neurons have exceed threshold and fired at t
			//Update I[] +=I


			//Update membrane Voltages

				i++; //increment Neuron Index
		}//End For Each Neuron



		itStep ++;
	}



	//Close Output File
	ofRSV.close();
	ofFSV.close();
	ofRaster.close();
	//Free Memory
	vNeurons.clear();

	//CleanUp Weight Vector memory
	for (register unsigned int i=0;i<inetSize;i++)
		delete [] W[i];

	delete [] W;

	delete [] I;
	delete [] F;

	gsl_rng_free(r);


	system("gnuplot SpikeRaster.gplot");

	return 0;
}
