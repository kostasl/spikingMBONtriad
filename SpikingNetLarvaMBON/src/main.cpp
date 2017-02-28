///*
/// \brief Object oriented research code implementing spiking model of the Larval MBON triad
/// Here I use an old implementation of a Leaky integrate and fire neuron which implements a 4th Order Runge-Kutta method.
///
///
/// \todo Separate Synapse Plasticity Models, between SONG and Switch rule - Making independent models
///*





#include "stdafx.h"
#include "isynapse.h"
#include "IFNeuron.h"
#include "synapseSW.h"
#include "synapseEnsemble.h"

#include "PoissonSource.h"
#include "PoissonNeuron.h"
#include "synapticTransmission.h"

// ###GSL Note: For the library to work, I had to change to the Multithreaded version WinGsl_md.lib
// Also under Properties->C/C++->Code GEneration->Run Time Library Change to Multithreaded Debug
//#include <WinGsl.h >
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


//Function Prototypes


using namespace std;

//Simulation Time step
//const float h=0.0002f;
const float h=0.0001f;

static char FilePath[_MAX_PATH]; // _MAX_PATH represents the longest possible path on this OS



//Variables for IFTEST
static const int iTestFq	= 80;
static const int iNoExSynapses = 10; //Number of synapses to test IFNeuron
static const int iNoInhSynapses = 0;//200; //Inhibitory synapses
static const uint IFSimulationTime = 60000000;
static const int NoSyns	= 1;//Switch rule Ensemble's number of Synapses


//THESE PARAMETERS ONY AFFECT SWITCH RULE
static const float APOT	= 1.0f;///(iTestFq*200); //(iTestFq*IFSimulationTime*h); Too Slow
static const float ADEP	= 0.95f;///(iTestFq*200); //(iTestFq*IFSimulationTime*h); Too Slow
static const float tafPOT	= 0.013f;
static const float tafDEP	= 0.020f;
static const int nPOT	= 1; //Change to 1 for Poisson Neuron Test
static const int nDEP	= 1;
/// \todo
static const float StartStrength = 10.0f;


int main(int argc, char *argv[])
{

   testIFNeuron(iNoExSynapses,iNoInhSynapses,IFSimulationTime);

    return 1;
}



void testIFNeuron(int iNoExSynapses,int iNoInhSynapses,uint uiSimulationTime)
{
    int	verboseperiod = 50000;
    int timetolog = verboseperiod;

    synapseEnsemble<synapseSW,NoSyns> *synstest[iNoExSynapses+iNoInhSynapses];// = new synapseEnsemble[iNoExSynapses];
    PoissonSource *Ps[iNoExSynapses+iNoInhSynapses];//Create Separate Poisson Sources for each afferent
    int cnt=0;
    int spikecnt = 0;
    int TotalPostSpikes = 0; //Used to get the Average Post Rate
    bool bPostEvent = false;
    double t=0;
    double Vm;
    float gamma = APOT*nPOT*tafPOT/(ADEP*nDEP*tafDEP);

    cout << "---- Test IF Neuron  with Gamma: " << gamma << endl;

    //Create IFNeuron
    IFNeuron* ifn = new IFNeuron(h,1);

    //Create and register Exhitatory Synapses
    for (int i=0;i<iNoExSynapses;i++)
    {
        // Replaced synapseEnsemble(float simTimeStep,int SynapsesNumber,float A1,float A2,float tafPOT,float tafDEP,int nPOT,int nDEP,float StartStrength,short sourceID,short ID,bool bNoPlasticity,unsigned short SourceFireRate)
        //synapseSW(float A1,float A2,float tafPOT,float tafDEP,int nPOT,int nDEP,float Sreset,bool bNoPlasticity);
        synapseSW* osyn = new synapseSW(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,StartStrength,true); //This synapse is going to be copied into the ensemble

        synapseEnsemble<synapseSW,NoSyns>* posynen = new synapseEnsemble<synapseSW,NoSyns>(h,*osyn);
        synstest[i] = posynen; //new synapseEnsemble<synapseSW,1>(h,osyn); //No Plasticity
        //posynen->RegisterNeuron(ifn);
        ifn->RegisterAfferent((ISynapseEnsemble*)(posynen)); //Let the nEuron Know a bundle of synapses is connecting to it

        //Create New Efferent
        Ps[i] = new PoissonSource(iTestFq,h,0.001);
    }

    //Create and register Inhibitory Synapses
    for (int i=iNoExSynapses;i< iNoExSynapses+ iNoInhSynapses;i++)
    {
        //No Plasticity and StartStrength -1 Inhibitory
        synapseSW osynI(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,-1,true); //This synapse is going to be copied into the ensemble
        synstest[i] = new synapseEnsemble<synapseSW,NoSyns>(h,osynI); //No Plasticity
        ifn->RegisterAfferent((ISynapseEnsemble*)synstest[i]);

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
                synstest[i]->SpikeArrived(ISynapse::SPIKE_PRE);


                spikecnt++;
            }
            else
            {
                synstest[i]->StepNoSpike();
            }


            //Log Strength of Every Exhitatory Synapse
            //if (timetolog==0 && i < iNoExSynapses) synstest[i]->logtofile(ofile);
        }
        Vm = ifn->StepRK_UpdateVm();
        if (ifn->PostSpikeOccured()) TotalPostSpikes++;

        if (timetolog==0)
        {
            cout<< cnt<< " t: "<< t << " Vm:" <<  Vm << " Sj:" << synstest[0]->getAvgStrength() << " F Rate:" << ifn->getFireRate() << " Avg Rate:" << TotalPostSpikes/t <<endl;
            timetolog=verboseperiod;

        }

    }
    cout << endl << "End. Spike Count: " << spikecnt << " Neuron Fire Rate: " << ifn->getFireRate() << " Avg Rate:" << TotalPostSpikes/t;
    cout << " time: " << t << " sec" << endl;


    //Close Synapse Log File
    ofile.close();


    delete ifn;
    //delete  &synstest;
    //delete &Ps;

}
