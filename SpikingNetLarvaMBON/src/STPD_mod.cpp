/*
# STPD_mod.cpp : Defines the entry point for the console application
# Uses the objects to setup simulations and tests.   
# An LTP and LTD test has been implemented in two seperate function. These write to a CSV file
# which can be used to draw the Exponential curves.
*/

#include "stdafx.h"
#include "synapseSW.h"
#include "IFNeuron.h"
#include "synapseEnsemble.h"
#include "PoissonSource.h"
#include "PoissonNeuron.h"
#include "synapticTransmission.h"

// ###GSL Note: For the library to work, I had to change to the Multithreaded version WinGsl_md.lib
// Also under Properties->C/C++->Code GEneration->Run Time Library Change to Multithreaded Debug
//#include <WinGsl.h >
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


using namespace std;

//Function Prototypes
void testSingleSynapse();
void testLTP(int); //Generate the LTP Exp Curve->Write to outputLTP.csv
void testLTD(int); //Generate the LTD Exp Curve->Write to outputLTD.csv
void testGSL(void); //Simple Test to check if GSL works
void testSynaptTransmission(void); //Examine the Profile and test the class by file output:synTest.csv
void testIFNeuron();
void testPoissonTrain(int nspikes); //Expectancy of single synapse Test
void testCorrelationIFNeuron(int nspikes); //Many Synapses from 1 Afferent for 2..n spikes. Mean And Var
double testS2SPoissonTrain(int nspikes); //Expectancy using Spike2Spike Jumps
//Test A number of synapses under the same input
void testSynapseCorrelation(int nSynapses);

void testSpikeTrain(float spikeInterval,int nspikes);
void BCMRule(); //Obtain the BCM like curves for expected Change
double testS2SPoissonTrainWithBursts(int nspikes,int burstLength);
double testS2SPoissonTrainWithBurstsVarPostRate(int nspikes,int burstLength);

//Simulation Time step
//const float h=0.0002f;
const float h=0.0001f;

char FilePath[PATH_MAX]; // _MAX_PATH represents the longest possible path on this OS

//Setup test variables
ISynapse* stest;
ISynapseEnsemble* synstest;
PoissonSource* Ps1,*Ps2;

	//Variables for IFTEST
    int iTestFq	= 40;
    const int iNoExSynapses = 1; //Number of synapses to test IFNeuron
    const int iNoInhSynapses = 0;//200; //Inhibitory synapses
    int IFSimulationTime = 60000000;
    int NoSyns	= 1;//Switch rule Ensemble's number of Synapses
	///According to Appleby&Elliott A~10^-3 or lower achieves stable competitive dynamics
	//THESE PARAMETERS ONY AFFECT SWITCH RULE
	float APOT	= 1.0f;///(iTestFq*200); //(iTestFq*IFSimulationTime*h); Too Slow
    float ADEP	= 0.95f;///(iTestFq*200); //(iTestFq*IFSimulationTime*h); Too Slow
    float tafPOT	= 0.013f;
    float tafDEP	= 0.020f;
    int nPOT	= 1; //Change to 1 for Poisson Neuron Test
    int nDEP	= 1;
	float StartStrength = 1.0f;


int STDP_mod_tmain(int argc, char* argv[])
{
	//Timing The process
	clock_t start, finish;
	double duration;
	
	start = clock();

	//testGSL();
	//testLTP(60);
	//testLTD(60);
	//testSingleSynapse();
	//testSynaptTransmission();
	
	

	//testSynapseCorrelation(100); //Test 2 Synapses

	//testIFNeuron();
	
	//testSpikeTrain(0.014f,4);
	//Run Expected Change Test using Poisson Source and Poisson Neuron
	//testPoissonTrain(200); //1 Afferent Mean over many synapses
	//testCorrelationIFNeuron(100); //1 Afferent with 120 synapses
	//testS2SPoissonTrain(100); //With n=1
	//BCMRule(); //For BCMTest Change Post rate to Pre-5
	//testS2SPoissonTrainWithBursts(100,15);
	testS2SPoissonTrainWithBurstsVarPostRate(100,120);
	finish = clock();

	///Print Duration of Run
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf( "\n Runtime was %2.1f seconds\n", duration );

	//Wait till Key Press
	char t;
	cin >> t;

	return 0;
}

void testSingleSynapse()
{
		stest = new synapse(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,0,false);

	//stest->SpikeArrived(0,synapse::SPIKE_POST); //Send Spike
	cout << "--- TEST the stochastic return to the OFF state " << endl;
	cout << "Off after 1ms :" << stest->StochasticReturnedtoOFF(0.001)<< endl;
	cout << "Off after 20ms :" << stest->StochasticReturnedtoOFF(0.02)<< endl;
	cout << "Off after 60ms :" << stest->StochasticReturnedtoOFF(0.06)<< endl;

	cout << "Off after 1sec :" << stest->StochasticReturnedtoOFF(1)<< endl;
	stest->Reset();
	cout << "TEST Potentiation of Single Synapse" << endl << endl;
	stest->SpikeArrived(10,synapse::SPIKE_PRE);
	cout << "Pre Spike" << endl;
	cout << "New S:" << stest->getStrength() << endl;

	stest->SpikeArrived(0.0026,synapse::SPIKE_POST);
	cout << "Post Spike 2.6ms" << endl;
	cout << "New S:" << stest->getStrength() << endl;


	stest->SpikeArrived(0.006,synapse::SPIKE_PRE);
	cout << "Pre Spike 6ms after Post" << endl;
	cout << "New S:" << stest->getStrength() << endl;

	//New Experiment pðp
	stest->Reset();
	cout << "------- DEPOT ->pðp" << endl;
	stest->SpikeArrived(10,synapse::SPIKE_POST);
	cout << "Start Post Spike 10s" << endl;
	cout << "New S:" << stest->getStrength() << endl;

	stest->SpikeArrived(0.0065,synapse::SPIKE_PRE);
	cout << "Pre Spike 6.5ms after Post" << endl;
	cout << "New S:" << stest->getStrength() << endl;
	
	stest->SpikeArrived(0.0005,synapse::SPIKE_POST);
	cout << "Post Spike after 0.5ms" << endl;
	cout << "New S:" << stest->getStrength() << endl << endl;

	stest->Reset();
	cout << "Quadriple TEST Potentiation of Single Synapse" << endl << endl;
	stest->SpikeArrived(10,synapse::SPIKE_PRE);
	cout << "Pre Spike" << endl;
	cout << "New S:" << stest->getStrength() << endl;

	stest->SpikeArrived(0.0088,synapse::SPIKE_POST);
	cout << "Post Spike 8.8ms" << endl;
	cout << "New S:" << stest->getStrength() << endl;


	stest->SpikeArrived(0.0106,synapse::SPIKE_POST);
	cout << "POst Spike 10.6ms after Post" << endl;
	cout << "New S:" << stest->getStrength() << endl;


	stest->SpikeArrived(0.0096,synapse::SPIKE_PRE);
	cout << "Pre Spike 9.6ms after Post" << endl;
	cout << "New S:" << stest->getStrength() << endl;



	stest->~synapse();

	delete stest;

}

//Assume Global vars h,NoSyns,APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP have been set
void testLTP(int trials)
{
	cout << endl << " ###  Testing LTP 1-80ms Average over:" << trials << endl;
	getcwd(FilePath, _MAX_PATH); // reads the current working directory into the array FilePath 
	ofstream ofile( strcat(FilePath,"\\..\\output\\outputLTP.csv"), ios::out );
	
	//Do Headers
	ofile << "k,t,S" << endl;

	//Ps1 = new PoissonSource(60,h,0.001);

	synstest = new synapseEnsemble(h,NoSyns,(APOT/trials),(ADEP/trials),tafPOT,tafDEP,nPOT,nDEP,0);
	//int cnt=0;
	int k=1;
	float timeWait=0.0f;
	while (k<81)
	{
		///Do 100 trials for each Spike time diff
		double Sg=0;
		synstest->Reset();

		for (int trial=0;trial<trials;trial++)
		{
				////WaiT 160ms between spike pairs so the switch returns to OFF
				//for (int w=0;w<166;w++) synstest->StepNoSpike();
				synstest->StepNoSpike(1);

				timeWait=0;
				synstest->SpikeArrived(synapse::SPIKE_PRE);

				timeWait=(double)k/1000.0;
				synstest->StepNoSpike(timeWait);

				synstest->SpikeArrived(synapse::SPIKE_POST);

			Sg=synstest->getAvgStrength();
		}//Trial

		//Save the average over 100 trials
		ofile << k << "," << timeWait << "," << Sg << endl;

		k++; //k ms tested after 60 pairs, do next one
		
	}
//Done
ofile.close();
	cout << "- Done" << endl;
}

//Assume Global vars h,NoSyns,APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP have been set
void testLTD(int trials)
{
	getcwd(FilePath, _MAX_PATH); // reads the current working directory into the array FilePath 
	ofstream ofile( strcat(FilePath,"\\..\\output\\outputLTD.csv"), ios::out );
	
	//Do Headers
	ofile << "k,t,S" << endl;
	cout << "### Test Ensemble -  LTD over 1-80ms Average over: " << trials  << endl;
//	Ps1 = new PoissonSource(60,h,0.001);
	synstest = new synapseEnsemble(h,NoSyns,(APOT/trials),(ADEP/trials),tafPOT,tafDEP,nPOT,nDEP,0);
	//int cnt=0;
	int k=81;
	float timeWait=0.0f;

	while (k>0)
	{
		///Do 100 trials for each Spike time diff
		double Sg=0;
		synstest->Reset();

		for (int trial=0;trial<trials;trial++)
		{
				////WaiT 160ms between spike pairs so the switch returns to OFF
				//for (int w=0;w<166;w++) synstest->StepNoSpike();
				synstest->StepNoSpike(1);

				timeWait=0;
				
				synstest->SpikeArrived(synapse::SPIKE_POST);
				timeWait=(double)k/1000.0;
				synstest->StepNoSpike(timeWait);

				
				synstest->SpikeArrived(synapse::SPIKE_PRE);
			
				Sg=synstest->getAvgStrength();
		}//Trial

		//Save the average over 100 trials
		ofile << k << "," << -timeWait << "," << Sg << endl;

		k--; //k ms tested after 60 pairs, do next one
		
	}

	cout << "- Done" << endl;
	ofile.close();
}


//Function to test the rates of the Poisson Class and 
// the Avg change using poisson pre and post 
void testPoisson()
{
/*
	cout << "Test Poisson" << endl;
	//Init The poisson Source
	Ps1 = new PoissonSource(60,h,0.001);
	Ps2 = new PoissonSource(60,h,0.001);

			if (Ps1->drawSpikeEvent())
			{
				synstest->SpikeArrived(synapse::SPIKE_PRE);
				cnt++;
				ofile << "1,";
			}
			else {
				ofile << "0,";
				synstest->StepNoSpike();
			}
			if (Ps1->drawSpikeEvent())
			{
				synstest->SpikeArrived(synapse::SPIKE_POST);
				ofile << "1,";
			}
			else { 
				ofile << "0,";
				synstest->StepNoSpike();
			}

*/
}


void testGSL()
{
	  //const gsl_rng_type * T;
	  gsl_rng * r;
  
	 // select random number generator 
	 r = gsl_rng_alloc (gsl_rng_mt19937);
	 gsl_rng_set(r,18211820);

	cout << gsl_ran_gaussian(r, 1) << ", " << gsl_ran_gaussian(r, 1) ;
	gsl_rng_free(r);
}

void testSynaptTransmission()
{
	double tm=0;
    getcwd(FilePath, _MAX_PATH); // reads the current working directory into the array FilePath
	ofstream ofile( strcat(FilePath,"\\synTrans.csv"), ios::out );
	
	//Do Headers
	ofile << "#t,I" << endl;
	cout << "### Test Synaptic Transmission - Spike Profile" << endl;

	synapticTransmission St(h,1);

	for (int i=0;i<550;i++) ofile << (tm+=h) << "," << St.getCurrentStep() << endl;

	ofile.close();
}

void testIFNeuron()
{	
	int	verboseperiod = 50000;
	int timetolog = verboseperiod;
	
	synapseEnsemble *synstest[iNoExSynapses+iNoInhSynapses];// = new synapseEnsemble[iNoExSynapses];
	PoissonSource *Ps[iNoExSynapses+iNoInhSynapses];//Create Seperate Poisson Sources for each afferent
	int cnt=0;
	int spikecnt = 0;
	int TotalPostSpikes = 0; //Used to get the Average Post Rate
	bool bPostEvent = false;
	double t=0;
	double Vm;
	float gamma = APOT*nPOT*tafPOT/(ADEP*nDEP*tafDEP);

	cout << "---- Test IF Neuron  with Gamma: " << gamma << endl;
	
	//Create IFNeuron
	IFNeuron ifn(h,1);
	
	//Create and register Exhitatory Synapses
	for (int i=0;i<iNoExSynapses;i++)
	{
		synstest[i] = new synapseEnsemble(h,NoSyns,APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,StartStrength,0,i);
		ifn.RegisterAfferent(synstest[i]);

		Ps[i] = new PoissonSource(iTestFq,h,0.001);
	}

	//Create and register Inhibitory Synapses
	for (int i=iNoExSynapses;i< iNoExSynapses+ iNoInhSynapses;i++)
	{
		//No Plasticity and StartStrength -1 Inhibitory
		synstest[i] = new synapseEnsemble(h,NoSyns,APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,-1,0,i,true);
		ifn.RegisterAfferent(synstest[i]);

		//According to Song Inhibitory Fixed at 10Hz
		Ps[i] = new PoissonSource(10,h,0.001);
	}


	//Open Log file for Synapse Strength
	///Record Synapse Strengths
    getcwd(FilePath, _MAX_PATH); // reads the current working directory into the array FilePath
	char *buff = strcat(FilePath,"\\synsStrength.csv");
	ofstream ofile(buff, ios::out );

	ofile << "ID Source Target S" << endl;
	//////LOG File Opened

	//Main Simulation Loop - Ends when Simulation time reached
	while (cnt < IFSimulationTime)
	{
		timetolog--;
		t+=h;
		cnt++;
		//Run Through each afferent and Poisson Source
		for (int i=0;i<iNoExSynapses+iNoInhSynapses;i++)
		{

			if (Ps[i]->drawSpikeEvent())
			{
				synstest[i]->SpikeArrived(synapse::SPIKE_PRE);


				spikecnt++;
			}
			else
			{
				synstest[i]->StepNoSpike();
			}

			
			//Log Strength of Every Exhitatory Synapse
			//if (timetolog==0 && i < iNoExSynapses) synstest[i]->logtofile(ofile);
		}
		Vm = ifn.StepRK_UpdateVm();
		if (ifn.PostSpikeOccured()) TotalPostSpikes++; 

		if (timetolog==0)
		{		
			cout<< cnt<< " t: "<< t << " Vm:" <<  Vm << " Sj:" << synstest[0]->getAvgStrength() << " F Rate:" << ifn.getFireRate() << " Avg Rate:" << TotalPostSpikes/t <<endl;
			timetolog=verboseperiod;

		}

	}
	cout << endl << "End. Spike Count: " << spikecnt << " Neuron Fire Rate: " << ifn.getFireRate() << " Avg Rate:" << TotalPostSpikes/t;
	cout << " time: " << t << " sec" << endl; 


	//Close Synapse Log File
	ofile.close();



	//delete  &synstest;
	//delete &Ps;

}

///The n synapses test from the same afferent
void testSynapseCorrelation(int nSynapses)
{	
	float StartStrength = 1;

	int	verboseperiod = 1000;
	int timetolog = 0;
	
	synapseEnsemble *synstest[150];// = new synapseEnsemble[iNoExSynapses];
	PoissonSource *Ps[1];// 1 Poisson Source for all Synapses
	int cnt=0;
	int spikecnt = 0;
	double t=0;
	bool SpikeEvent = false;
	bool PostSpikeEvent = false;
	float gamma = APOT*nPOT*tafPOT/(ADEP*nDEP*tafDEP);
	cout << "---- Test Poisson Neuron with Gamma: " << gamma << endl;
	
	//Create PNeuron - Set initial rate to TestFq*StartStrength
	PoissonNeuron pn(h,1,iTestFq-10,true); //Fixed Rate=true
	
	//Create and register Exhitatory Synapses
	for (int i=0;i<nSynapses;i++)
	{
		synstest[i] = new synapseEnsemble(h,1,APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,StartStrength,0,i,false,iTestFq);
		pn.RegisterAfferent(synstest[i]);

	}
	//A single Source for correlation Testing
		Ps[0] = new PoissonSource(iTestFq,h,0.001);


	//Open Log file for Synapse Strength
	///Record Synapse Strengths
    getcwd(FilePath, _MAX_PATH); // reads the current working directory into the array FilePath
	char *buff = strcat(FilePath,"\\..\\output\\SynCorrelationWithPoisson40-30.csv");
	ofstream ofile(buff, ios::out );

	ofile << "ID Source Target S" << endl;
	//////LOG File Opened

	//Main Simulation Loop - Ends when Simulation time reached
	//while (cnt < IFSimulationTime)
	 while (spikecnt < 10000)
	 {
		t+=h;
		cnt++;
		//Run Through each afferent and Poisson Source
		SpikeEvent = Ps[0]->drawSpikeEvent(); 
		if (SpikeEvent) spikecnt++;
			
			for (int i=0;i<nSynapses;i++)
			{
				if (SpikeEvent)
				{
					
					synstest[i]->SpikeArrived(synapse::SPIKE_PRE);
				}
				else {
					synstest[i]->StepNoSpike();
				}
				//Log Strength of Every Exhitatory Synapse
				if (PostSpikeEvent || SpikeEvent) synstest[i]->logtofile(ofile);
			}//For every Synapse
		//Step will notify synapses of any spike events
		
		PostSpikeEvent = pn.StepDrawSpike();
		if (PostSpikeEvent) spikecnt++;

		if (timetolog==verboseperiod)
		{		
			cout << cnt << " Beta:" <<  spikecnt/t << " Sj:" << synstest[0]->getAvgStrength() << " SpikesCnt:" << spikecnt <<endl;
			timetolog=0;
		}
		timetolog++;

	} //While Loop

	cout << endl << "End. Spike Count: " << spikecnt << " Neuron Fire Rate: " << pn.getFireRate();
	cout << " time: " << t << " sec" << endl; 


	//Close Synapse Log File
	ofile.close();



	//delete  &synstest;
	//delete &Ps;

}


// Do tests to obtain mean and variance of Ds under fixed spike trains pre-post pairs
void testSpikeTrain(float spikeInterval,int nspikes)
{
	int trials = 1000000;
	cout << "Testing Spike train on Single Synapse. Average over :" << trials <<endl;

	stest = new synapse(1,1,tafPOT,tafDEP,1,1,0,true);
	double Ds=0;
	double totalDs=0;


    getcwd(FilePath, _MAX_PATH); // reads the current working directory into the array FilePath
	char *buff = strcat(FilePath,"\\..\\output\\SynTrain.csv");
	ofstream ofile(buff, ios::out );

	ofile << "trial Ds" << endl;


	for (int t=0;t<trials;t++)
	{
		Ds=0;

        //cout << endl;
		for (int i=1;i<=nspikes;i++)
		{
			if ((i%2)==0)
			{
				Ds+=stest->SpikeArrived(spikeInterval,synapse::SPIKE_POST);
				//cout<< "Post-";
			}
			else
			{
				Ds+=stest->SpikeArrived(spikeInterval,synapse::SPIKE_PRE);
				//cout<< "Pre-";
			}	
		}
		//if (Ds <0) cout << "-Ve Ds:" << Ds << endl;
		ofile << t << " " <<  Ds << endl;
		totalDs+=Ds;
		stest->Reset();
	}

totalDs=totalDs/(double)trials;
delete stest;
ofile.close();

cout << endl << endl << " Expected Value of Ds:" << totalDs << " over " << trials<< " trials";
}


// Do tests to obtain mean and variance of Ds under Poisson spike trains pre-post pairs
//Use IF Neuron 1 Afferent
void testPoissonTrain(int nspikes)
{
	int trials = 1000000;
	double *Ds= new double[1000000];

	cout << "Testing Poisson Spike train on Single Synapse to an IF Neuron. Average over :" << trials <<endl;

	float StartStrength = 1.0;
	//Period of trials after which output is shown
	int	verboseperiod = 1000;
	int timetolog = 0;
	
	//Create IFNeuron
	IFNeuron ifn(h,1);

	synapseEnsemble *synstest[15];// = new synapseEnsemble[iNoExSynapses];
	PoissonSource *Ps[2];//Create Seperate Poisson Sources afferents 
	double Vm;
	int cnt=0;
	int nSynapses = 2; //Ensembles
	int NoSyns = 120; //Individual need enough to cause Post Firing
	int spikecnt = 0;
	int TotalPostRate = 0; //Used to get the Average Post Rate
	double totalDs=0;
	double lastSj=0;
	double Mean;
	double t=0;
	bool SpikeEvent = false;
	bool PostSpikeEvent = false;

	bool bPreSpike = false;
	bool bPostSpike = false;

	float gamma = APOT*nPOT*tafPOT/(ADEP*nDEP*tafDEP);
	
	
	cout << "---- Test Poisson Neuron with Gamma: " << gamma << endl;
	
	//Create PNeuron - Set initial rate to TestFq*StartStrength
	//PoissonNeuron pn(h,1,iTestFq,true);

	//Synapse Under Test -Plastic
	synstest[0] = new synapseEnsemble(h,1,APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,0,0,0,false,iTestFq);
    ifn.RegisterAfferent(synstest[0]);

	//Create and register Exhitatory Synapses
	//These Set The Post rate Non Plastic
	for (int i=1;i<nSynapses;i++)
	{
		synstest[i] = new synapseEnsemble(h,NoSyns,APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,StartStrength,0,i,true,iTestFq);
		ifn.RegisterAfferent(synstest[i]);
	}
		//Source of Synapse Under Test
		Ps[0] = new PoissonSource(iTestFq,h,0.001);
	//A single Source for 1 afferent of 150 Synapses Testing
		Ps[1] = new PoissonSource(iTestFq,h,0.001);


	//Open Log file for Synapse Strength
	///Record Synapse Strengths
    getcwd(FilePath, _MAX_PATH); // reads the current working directory into the array FilePath
	char *buff = strcat(FilePath,"\\..\\output\\SynPoissonExpIFNeuron.csv");
	ofstream ofile(buff, ios::out );

	ofile << "Spikes Mean Variance" << endl;
	//////LOG File Opened
	for (int s=13;s<=nspikes;s++)
	{
	totalDs=0;
	t=0;
	TotalPostRate=0;
	cout<< " Spikes :" << s <<endl;
		//Main Simulation Loop - Ends when Simulation time reached
	for(int tm=1;tm<=trials;tm++)
	{
		spikecnt=0; //reset count
		synstest[0]->Reset();
		//ifn.Reset();
		
		//synstest[0]->SpikeArrived(synapse::SPIKE_PRE);
		
		//Only works for
		bPreSpike = false; //Added so to ignore pi-pi po-po contributions
		bPostSpike = false;
		lastSj = synstest[0]->getAvgStrength();
		while (spikecnt < s)
		{
			t+=h;
				
				//Do the Pre Firing of the independent Afferent
				//if (nSynapses > 1)
				if (Ps[1]->drawSpikeEvent()) 
					synstest[1]->SpikeArrived(synapse::SPIKE_PRE);
				else
					synstest[1]->StepNoSpike();

				//Do the Synapse Under Test
				SpikeEvent = Ps[0]->drawSpikeEvent(); 	
				if (SpikeEvent)
				{
					
					synstest[0]->SpikeArrived(synapse::SPIKE_PRE);
					spikecnt++;
					bPreSpike = true;
					//cout << "pi-";
				}

				//No Spike Event, Step Sim time Fwd	
				if (!SpikeEvent)
				{
					synstest[0]->StepNoSpike();
				}

				Vm = ifn.StepRK_UpdateVm(); //Do Step
				if (ifn.PostSpikeOccured())
				{
					spikecnt++;
					//cout << "po-";

				}
		}//While Spikes Loop
		//cout << endl;
		//Log Strength of Every Exhitatory Synapse
//		synstest[0]->logtofile(ofile);

		timetolog++;
		Ds[tm-1] =  synstest[0]->getAvgStrength() - lastSj;

		totalDs+= Ds[tm-1];
		TotalPostRate+=ifn.getFireRate();
		Mean = (totalDs/tm); //Per Synapse 

		if (timetolog>verboseperiod)
		{		
			cout << tm << " PreRate Set:" <<  iTestFq << " Avg Post:" << (TotalPostRate/tm) << " MeanDs:" << Mean<< " Vm: " << (Vm) << endl;
			timetolog=0;
		}

	} //For trials Loop


	//Calc Mean And Var
	double Var=0;
	Mean = totalDs/trials;
	for (int tm=0;tm<trials;tm++) Var+=(Ds[tm]*Ds[tm]);
	Var = Var/trials - Mean*Mean;

	cout << endl << endl << " Mean strength s:" << (Mean) << " Var:"<< Var << " over " << trials<< " trials Total Ds:" << totalDs;

	cout << endl << "End. Spike Count: " << spikecnt << " Neuron Fire Rate: " << ifn.getFireRate();
	cout << " time: " << t << " sec" << endl; 


	ofile << s << " " << Mean << " " << Var << endl;
	} //For Spikes Loop



	//Close Synapse Log File
	ofile.close();
	

}

// Do tests to obtain mean and variance of Ds under Poisson spike trains pre-post pairs
//Use IF Neuron 1 Afferent, 1 Synapse with 1 switch, and another with 120 switches. Mean And Variance of Both is saved
void testCorrelationIFNeuron(int nspikes)
{
	int trials = 700000;
	double *Ds= new double[1000000];
	double *Ds2= new double[1000000];

	cout << "Testing Poisson Spike train of 100 and 1 Synapses to an IF Neuron. Average over :" << trials <<endl;

	float StartStrength = 1.0;
	//Period of trials after which output is shown
	int	verboseperiod = 10000;
	int timetolog = 0;
	
	//Create IFNeuron
	IFNeuron ifn(h,1);

	synapseEnsemble *synstest[15];// = new synapseEnsemble[iNoExSynapses];
	PoissonSource *Ps[2];//Create Seperate Poisson Sources afferents 
	double Vm;
	int cnt=0;
	int nSynapses = 2; //Ensembles
	int NoSyns = 100; //Individual need enough to cause Post Firing
	StartStrength = 1.0f;
	int spikecnt = 0;
	int TotalPostRate = 0; //Used to get the Average Post Rate
	double totalDs=0;
	double totalDs2=0;

	double lastSj=0;
	double lastSj2=0;

	double Mean,Mean2;
	double t=0;
	bool SpikeEvent = false;
	bool PostSpikeEvent = false;

	bool bPreSpike = false;
	bool bPostSpike = false;

	float gamma = APOT*nPOT*tafPOT/(ADEP*nDEP*tafDEP);
	
	
	cout << "---- Test IFNeuron with Gamma: " << gamma << endl;
	
	//Create PNeuron - Set initial rate to TestFq*StartStrength
	//PoissonNeuron pn(h,1,iTestFq,true);

	//Synapse Under Test -Plastic
	synstest[0] = new synapseEnsemble(h,1,APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,1,0,0,false,iTestFq);
    ifn.RegisterAfferent(synstest[0]);

	//Create and register Exhitatory Synapses
	//These Set The Post rate Non Plastic
	for (int i=1;i<nSynapses;i++)
	{
		synstest[i] = new synapseEnsemble(h,NoSyns,APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,StartStrength,0,i,false,iTestFq);
		ifn.RegisterAfferent(synstest[i]);
	}
		//Source of Synapse Under Test
		Ps[0] = new PoissonSource(iTestFq,h,0.001);
	//A single Source for 1 afferent of 150 Synapses Testing
		Ps[1] = new PoissonSource(iTestFq,h,0.001);


	//Open Log file for Synapse Strength
	///Record Synapse Strengths
    getcwd(FilePath, _MAX_PATH); // reads the current working directory into the array FilePath
	char *buff = strcat(FilePath,"\\..\\output\\ExpVarIFNeuron-1and100Syns2-100Spikes.csv");
	ofstream ofile(buff, ios::out );

	ofile << "Spikes Mean1 Variance1 Mean100 Variance100 PostRate" << endl;
	//////LOG File Opened
	for (int s=2;s<=nspikes;s++)
	{
	totalDs=0;
	totalDs2=0;
	t=0;
	TotalPostRate=0;
	cout<< " Spikes :" << s <<endl;

		//Main Simulation Loop - Ends when Simulation time reached
	for(int tm=1;tm<=trials;tm++)
	{
		spikecnt=0; //reset count
		synstest[0]->Reset();
		synstest[1]->Reset(); // New Addition
		//ifn.Reset();
		
		//synstest[0]->SpikeArrived(synapse::SPIKE_PRE);
		
		//Only works for
		bPreSpike = false; //Added so to ignore pi-pi po-po contributions
		bPostSpike = false;
		lastSj = synstest[0]->getAvgStrength();
		lastSj2 = synstest[1]->getAvgStrength();

		while (spikecnt < s)
		{
			t+=h;
				
				if (nSynapses > 1)
				{
					//Do the Pre Firing of the independent Afferent

					if (Ps[1]->drawSpikeEvent()) 
						synstest[1]->SpikeArrived(synapse::SPIKE_PRE);
					else
						synstest[1]->StepNoSpike();
				}				
				//Do the Synapse Under Test
				SpikeEvent = Ps[0]->drawSpikeEvent(); 	
				if (SpikeEvent)
				{
					
					synstest[0]->SpikeArrived(synapse::SPIKE_PRE);
					spikecnt++;
					bPreSpike = true;
					//cout << "pi-";
				}

				//No Spike Event, Step Sim time Fwd	
				if (!SpikeEvent)
				{
					synstest[0]->StepNoSpike();
				}

				Vm = ifn.StepRK_UpdateVm(); //Do Step
				if (ifn.PostSpikeOccured())
				{
					spikecnt++;
					//cout << "po-";

				}
		}//While Spikes Loop
		//cout << endl;
		//Log Strength of Every Exhitatory Synapse
//		synstest[0]->logtofile(ofile);

		timetolog++;
		Ds[tm-1] =  synstest[0]->getAvgStrength() - lastSj;
		Ds2[tm-1] = synstest[1]->getAvgStrength() - lastSj2; ///Save the ensembles Change
		
		totalDs+= Ds[tm-1];
		totalDs2+= Ds2[tm-1];

		TotalPostRate+=ifn.getFireRate();
		Mean = (totalDs/tm); //Per Synapse 
		Mean2 = (totalDs2/tm); //Per Synapse 

		if (timetolog>verboseperiod)
		{		
			cout << tm << " PreRate Set:" <<  iTestFq << " Avg Post:" << (TotalPostRate/tm) << " MeanDs:" << Mean<< " 100synsMean:" << Mean2 << endl;
			timetolog=0;
		}

	} //For trials Loop


	//Calc Mean And Var
	double Var=0;
	double Var2=0;

	Mean = totalDs/trials;
	Mean2 = totalDs2/trials;

	Var=0;
	for (int tm=0;tm<(trials);tm++)
	{
		Var+=(Ds[tm]*Ds[tm]);
		Var2+=(Ds2[tm]*Ds2[tm]);
	}

	Var = Var/trials - Mean*Mean;
	Var2 = Var2/trials - Mean2*Mean2;

	cout << endl << endl << " Mean strength of 1 S:" << (Mean) << " Var:"<< Var << " over " << trials<< " trials 100 SynMean:" << Mean2 << " Var:" << Var2;
	

	cout << endl << "End. Spike Count: " << spikecnt << " Neuron Fire Rate: " << (TotalPostRate/trials);
	cout << " time: " << t << " sec" << endl; 


	ofile << s << " " << Mean << " " << Var << " " << Mean2 << " "<< Var2 << " " << (TotalPostRate/trials) << endl;
	} //For Spikes Loop



	//Close Synapse Log File
	ofile.close();
	

}



// Do tests to obtain mean and variance of Ds under Poisson spike trains pre-post pairs
//Time-step is variable and jumps from spike to spike event -
double testS2SPoissonTrain(int nspikes)
{
	int trials = 1000000;
	cout << "Testing Poisson Spike train on Single Synapse to A poisson Neuron. Average over :" << trials <<endl;

	float StartStrength = 0.00;
	//Period of trials after which output is shown
	int	verboseperiod = 100000;
	int timetolog = 0;
	
	synapseEnsemble *synstest[1];// = new synapseEnsemble[iNoExSynapses];
	PoissonSource *Ps[2];//Create Seperate Poisson Sources for afferent and target
	int cnt=0;
	int nSynapses = 1;
	int spikecnt = 0;
	double totalDs=0;
	double *Ds= new double[10000000];
	double t=0;
	double Var=0;
	double Mean=0;
	bool SpikeEvent = false;
	bool PostSpikeEvent = false;

	bool bPreSpike = false;
	bool bPostSpike = false;

	double tpi=0; //Time of pi spike
	double tpo=0; //POst Spike time
	float gamma = APOT*nPOT*tafPOT/(ADEP*nDEP*tafDEP);
	
	
	cout << "---- Test Poisson Neuron with Gamma: " << gamma << endl;
	
	//Create PNeuron - Set initial rate to TestFq*StartStrength
	//PoissonNeuron pn(h,1,iTestFq,true);
	
	//Create and register Exhitatory Synapses
	for (int i=0;i<nSynapses;i++)
	{
		//Used to be A/nspikes/2
		synstest[i] = new synapseEnsemble(h,NoSyns,(APOT),(ADEP),tafPOT,tafDEP,nPOT,nDEP,StartStrength,0,i,false,iTestFq);
	}

	//A single Source 
		Ps[0] = new PoissonSource(iTestFq,h,0.001);
	//Target independent poisson process
		Ps[1] = new PoissonSource(iTestFq-5,h,0.001);

	//Open Log file for Synapse Strength
	///Record Synapse Strengths
    getcwd(FilePath, _MAX_PATH); // reads the current working directory into the array FilePath
	char *buff = strcat(FilePath,"\\..\\output\\SynPoissonExpVarPostRate.csv");
	ofstream ofile(buff, ios::out );

	//ofile << "ID Source Target S" << endl;
	ofile << "S" << endl;
	//////LOG File Opened

	//Main Simulation Loop - Ends when Simulation time reached
	for(int tm=1;tm<=trials;tm++)
	{
		cnt=0; //reset count
		spikecnt = 0;
		t=0; //Reset simulation time of trial

		synstest[0]->Reset();
		
		//synstest[0]->SpikeArrived(synapse::SPIKE_PRE);
		
		//Only works for
		bPreSpike = false; //Added so to ignore pi-pi po-po contributions
		bPostSpike = false;
		Ds[tm-1] = 0;
		while (spikecnt < nspikes)
		{
			
			
				tpo = Ps[1]->gsl_ran_exponential(); //Target
				tpi = Ps[0]->gsl_ran_exponential(); //Source

				if  (tpo<tpi)
				{
					synstest[0]->StepNoSpike(tpo);
					Ds[tm-1]+= synstest[0]->SpikeArrived(synapse::SPIKE_POST);
					spikecnt++;
					bPostSpike = true;
					//cout << "po-";
					cnt++;
					t+=tpo;
				}
				else
				{
					synstest[0]->StepNoSpike(tpi);
					Ds[tm-1]+= synstest[0]->SpikeArrived(synapse::SPIKE_PRE);
					spikecnt++;
					bPreSpike = true;
					cnt++;
					t+=tpi;
					//cout << "pi-";

				}

		}//While Spikes Loop
		//cout << endl;
		//Log Strength of Every Exhitatory Synapse
		//synstest[0]->logtofile(ofile);
		//ofile << Ds[tm-1] << endl; //Record the Total change after spike Train

		totalDs+=Ds[tm-1];

		Mean=totalDs/tm;
        Var+=abs((int)(Ds[tm-1]-Mean))/(trials-tm+1);
		
		timetolog++;
		if (timetolog>verboseperiod)
		{		
			cout << tm << " PreRate:" <<  iTestFq << " Real beta:" << (cnt/t) << " Sj:" << synstest[0]->getAvgStrength() << " Avg: " << (Mean) <<" Var:" << Var <<endl;
			timetolog=0;
		}

	} //For trials Loop

	cout << endl << "End. Spike Count: " << spikecnt << " Neuron Fire Rate: " << iTestFq;
	cout << " time: " << t << " sec" << endl; 


	//Close Synapse Log File
	ofile.close();
	Mean = totalDs/trials;
	Var=0;
	for (int tm=0;tm<trials;tm++) Var+=(Ds[tm]*Ds[tm]);
	Var = Var/trials - Mean*Mean;
	cout << endl << endl << " Mean strength s:" << (Mean) << " Var:"<< Var << " over " << trials<< " trials Total Ds:" << totalDs;


	delete Ds;
	delete Ps[0];
	delete Ps[1];
	return Mean; //Return the expected change
}

void BCMRule()
{
    getcwd(FilePath, _MAX_PATH); // reads the current working directory into the array FilePath
	char *buff = strcat(FilePath,"\\..\\output\\BCM50Spikes.csv");
	ofstream ofile(buff, ios::out );

	iTestFq=1;
	double Exp;
	ofile << "lpi EDs50" << endl;
	for (int i=1;i<201;i++)
	{
		iTestFq =i;
		//Test under 2 spikes for a Lpi=i
		Exp = testS2SPoissonTrain(50);
		ofile << iTestFq << " " << Exp << endl;
	}


	ofile.close();
}


////Do test for Exp And Variance Using Bursts of activity
// Do tests to obtain mean and variance of Ds under Poisson spike trains pre-post pairs
//For 2...n Spike length
double testS2SPoissonTrainWithBursts(int nspikes,int burstLength)
{
	int trials = 1000000;
	cout << "Testing Poisson Spike train on Single Synapse to A poisson Neuron with Burst. Average over :" << trials <<endl;

	float StartStrength = 0.00;
	//Period of trials after which output is shown
	int	verboseperiod = 100000;
	int timetolog = 0;
	
	synapseEnsemble *synstest[1];// = new synapseEnsemble[iNoExSynapses];
	PoissonSource *Ps[2];//Create Seperate Poisson Sources for afferent and target
	int cnt=0;
	int nSynapses = 1;
	int spikecnt = 0;
	double totalDs=0;
	double *Ds= new double[10000000];
	double t=0;
	double Var=0;
	double Mean=0;
	bool SpikeEvent = false;
	bool PostSpikeEvent = false;

	bool bPreSpike = false;
	bool bPostSpike = false;

	double tpi=0; //Time of pi spike
	double tpo=0; //POst Spike time
	float gamma = APOT*nPOT*tafPOT/(ADEP*nDEP*tafDEP);
	
	
	cout << "---- Test- Expectation variance with Bursts with Gamma: " << gamma << endl;
	
	//Create PNeuron - Set initial rate to TestFq*StartStrength
	//PoissonNeuron pn(h,1,iTestFq,true);
	
	//Create and register Exhitatory Synapses
	for (int i=0;i<nSynapses;i++)
	{
		//Used to be A/nspikes/2
		synstest[i] = new synapseEnsemble(h,NoSyns,(APOT),(ADEP),tafPOT,tafDEP,nPOT,nDEP,StartStrength,0,i,false,iTestFq);
	}

	//A single Source 
		Ps[0] = new PoissonSource(iTestFq,h,0.001);
	//Target independent poisson process
		Ps[1] = new PoissonSource(iTestFq-10,h,0.001);

	//Open Log file for Synapse Strength
	///Record Synapse Strengths
    getcwd(FilePath, _MAX_PATH); // reads the current working directory into the array FilePath
	char *buff = strcat(FilePath,"\\..\\output\\SynPoissonExpVarWithBurst15.csv");
	ofstream ofile(buff, ios::out );

	ofile << "Spikes15 Mean15 Variance15" << endl;
	//////LOG File Opened
	for (int s=2;s<=nspikes;s++)
	{
		totalDs=0;
		t=0;
		
		cout << " Spikes :" << s << " Burst Length:" << burstLength <<endl;
	//////LOG File Opened

	//Main Simulation Loop - Ends when Simulation time reached
	for(int tm=1;tm<=trials;tm++)
	{
		cnt=0; //reset count
		spikecnt = 0;
		t=0; //Reset simulation time of trial

		synstest[0]->Reset();
				
		bPreSpike = false; //Added so to ignore pi-pi po-po contributions
		bPostSpike = false;
		Ds[tm-1] = 0;
		while (spikecnt < s)
		{
			//Burst Time?
			if ((spikecnt % burstLength) == 0)
			{
				///Reset Synapse, Return to OFF and Strength 0 doesnt matter, we take the sum of changes
				synstest[0]->Reset();
				//cout << (spikecnt % burstLength)<< "Burst complete after :" << spikecnt <<endl;
			}
				tpo = Ps[1]->gsl_ran_exponential(); //Target
				tpi = Ps[0]->gsl_ran_exponential(); //Source

				if  (tpo<tpi)
				{
					synstest[0]->StepNoSpike(tpo);
					Ds[tm-1]+= synstest[0]->SpikeArrived(synapse::SPIKE_POST);
					spikecnt++;
					bPostSpike = true;
					//cout << "po-";
					cnt++;
					t+=tpo;
				}
				else
				{
					synstest[0]->StepNoSpike(tpi);
					Ds[tm-1]+= synstest[0]->SpikeArrived(synapse::SPIKE_PRE);
					spikecnt++;
					bPreSpike = true;
					cnt++;
					t+=tpi;
					//cout << "pi-";

				}

		}//While Spikes Loop
		//cout << endl;
		//Log Strength of Every Exhitatory Synapse
		//synstest[0]->logtofile(ofile);
		//ofile << Ds[tm-1] << endl; //Record the Total change after spike Train

		totalDs+=Ds[tm-1];

		Mean=totalDs/tm;
        Var+=abs((int)(Ds[tm-1]-Mean))/(trials-tm+1);
		
		timetolog++;
		if (timetolog>verboseperiod)
		{		
			cout << tm << " PreRate:" <<  iTestFq << " Real beta:" << (cnt/t) << " Sj:" << synstest[0]->getAvgStrength() << " Avg: " << (Mean) <<" Var:" << Var <<endl;
			timetolog=0;
		}

	} //For trials Loop

	cout << endl << "End. Spike Count: " << (spikecnt) << " Neuron Fire Rate: " << iTestFq;
	cout << " time: " << t << " sec" << endl; 


	Mean = totalDs/trials;
	for (int tm=0;tm<trials;tm++) Var+=(Ds[tm]*Ds[tm]);
	Var = Var/trials - Mean*Mean;
	cout << endl << endl << " Mean strength s:" << (Mean) << " Var:"<< Var << " over " << trials<< " trials Total Ds:" << totalDs;

	//Write to Log file
	ofile << s << " " << Mean << " " << Var << endl;

	}//For up to nspikes
	//Close Synapse Log File
	ofile.close();

	delete Ds;
	delete Ps[0];
	delete Ps[1];
	return Mean; //Return the expected change
}



////Do test for Exp And Variance Using Bursts of activity
// Do tests to obtain mean and variance of Ds under Poisson spike trains pre-post pairs
//For 2...n Spike length
///Used to Compare Against IFNeuron Results
double testS2SPoissonTrainWithBurstsVarPostRate(int nspikes,int burstLength)
{

int PostRate[] ={44,42,42,42,42,43,43,44,45,46,47,48,49,50,50,51,52,53,53,54,54,55,56,56,57,57,58,58,59,59,
				59,60,60,61,61,62,62,62,63,63,63,64,64,64,65,65,65,65,66,66,66,67,67,67,67,68,68,68,68,69,
				69,69,69,69,70,70,70,70,71,71,71,71,71,72,72,72,72,72,72,73,73,73,73,73,74,74,74,74,74,74,74,
				75,75,75,75,75,75,76,76};

	int trials = 1000000;
	cout << "Testing Poisson Spike train on Single Synapse to A poisson Neuron with Burst. Average over :" << trials <<endl;

	float StartStrength = 0.00;
	//Period of trials after which output is shown
	int	verboseperiod = 100000;
	int timetolog = 0;
	
	synapseEnsemble *synstest[1];// = new synapseEnsemble[iNoExSynapses];
	PoissonSource *Ps[2];//Create Seperate Poisson Sources for afferent and target
	int cnt=0;
	int nSynapses = 1;
	int spikecnt = 0;
	double totalDs=0;
	double *Ds= new double[10000000];
	double t=0;
	double Var=0;
	double Mean=0;
	bool SpikeEvent = false;
	bool PostSpikeEvent = false;
	int NoSyns = 100;
	bool bPreSpike = false;
	bool bPostSpike = false;

	double tpi=0; //Time of pi spike
	double tpo=0; //POst Spike time
	float gamma = APOT*nPOT*tafPOT/(ADEP*nDEP*tafDEP);
	
	
	cout << "---- Test- Expectation variance with Bursts with Gamma: " << gamma << endl;
	
	//Create PNeuron - Set initial rate to TestFq*StartStrength
	//PoissonNeuron pn(h,1,iTestFq,true);
	
	//Create and register Exhitatory Synapses
	for (int i=0;i<nSynapses;i++)
	{
		//Used to be A/nspikes/2
		synstest[i] = new synapseEnsemble(h,NoSyns,(APOT),(ADEP),tafPOT,tafDEP,nPOT,nDEP,StartStrength,0,i,false,iTestFq);
	}

	//A single Source 
		Ps[0] = new PoissonSource(iTestFq,h,0.001);
	//Shit Code in a hurry: Target independent poisson process
		Ps[1] = new PoissonSource(iTestFq-5,h,0.001);

	//Open Log file for Synapse Strength
	///Record Synapse Strengths
    getcwd(FilePath, _MAX_PATH); // reads the current working directory into the array FilePath
	char *buff = strcat(FilePath,"\\..\\output\\SynPoissonExpVarPostRate100Syns.csv");
	ofstream ofile(buff, ios::out );

	ofile << "Spikes Mean Variance PostRate" << endl;
	//////LOG File Opened
	for (int s=2;s<=nspikes;s++)
	{
		
		delete Ps[1];
		//Target independent poisson process
		Ps[1] = new PoissonSource(PostRate[s-2],h,0.001);

		totalDs=0;
		t=0;
		
		cout << " Spikes :" << s << " Burst Length:" << burstLength <<endl;
	//////LOG File Opened

	//Main Simulation Loop - Ends when Simulation time reached
	for(int tm=1;tm<=trials;tm++)
	{
		cnt=0; //reset count
		spikecnt = 0;
		t=0; //Reset simulation time of trial

		synstest[0]->Reset();
				
		bPreSpike = false; //Added so to ignore pi-pi po-po contributions
		bPostSpike = false;
		Ds[tm-1] = 0;
		while (spikecnt < s)
		{
			//Burst Time?
			if ((spikecnt % burstLength) == 0)
			{
				///Reset Synapse, Return to OFF and Strength 0 doesnt matter, we take the sum of changes
				synstest[0]->Reset();
				//cout << (spikecnt % burstLength)<< "Burst complete after :" << spikecnt <<endl;
			}
				tpo = Ps[1]->gsl_ran_exponential(); //Target
				tpi = Ps[0]->gsl_ran_exponential(); //Source

				if  (tpo<tpi)
				{
					synstest[0]->StepNoSpike(tpo);
					Ds[tm-1]+= synstest[0]->SpikeArrived(synapse::SPIKE_POST);
					spikecnt++;
					bPostSpike = true;
					//cout << "po-";
					cnt++;
					t+=tpo;
				}
				else
				{
					synstest[0]->StepNoSpike(tpi);
					Ds[tm-1]+= synstest[0]->SpikeArrived(synapse::SPIKE_PRE);
					spikecnt++;
					bPreSpike = true;
					cnt++;
					t+=tpi;
					//cout << "pi-";

				}

		}//While Spikes Loop
		//cout << endl;
		//Log Strength of Every Exhitatory Synapse
		//synstest[0]->logtofile(ofile);
		//ofile << Ds[tm-1] << endl; //Record the Total change after spike Train

		totalDs+=Ds[tm-1];

		Mean=totalDs/tm;
        Var+=abs((int)(Ds[tm-1]-Mean))/(trials-tm+1);
		
		timetolog++;
		if (timetolog>verboseperiod)
		{		
			cout << tm << " PreRate:" <<  iTestFq <<" Post:" << PostRate[s-2] <<" Real beta:" << (cnt/t) << " Sj:" << synstest[0]->getAvgStrength() << " Avg: " << (Mean) <<" Var:" << Var <<endl;
			timetolog=0;
		}

	} //For trials Loop

	cout << endl << "End. Spike Count: " << (spikecnt) << " Neuron Fire Rate: " << iTestFq;
	cout << " time: " << t << " sec" << endl; 


	Mean = totalDs/trials;
	for (int tm=0;tm<trials;tm++) Var+=(Ds[tm]*Ds[tm]);
	Var = Var/trials - Mean*Mean;
	cout << endl << endl << " Mean strength s:" << (Mean) << " Var:"<< Var << " over " << trials<< " trials Total Ds:" << totalDs;

	//Write to Log file
	ofile << s << " " << Mean << " " << Var << " " << PostRate[s-2] <<endl;

	}//For up to nspikes
	//Close Synapse Log File
	ofile.close();

	delete Ds;
	delete Ps[0];
	delete Ps[1];
	return Mean; //Return the expected change
}
