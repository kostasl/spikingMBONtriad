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
#include "math.h"
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
typedef std::map<std::string, std::ofstream*> FileMap;
typedef FileMap::iterator Iterator;

FileMap ofiles; //declare map file object


//Variables for IFTEST
static const int iTestFq	= 10;
static const int iNoExSynapses = 10; //Number of synapses to test IFNeuron
static const int iNoInhSynapses = 0;//200; //Inhibitory synapses
static const uint IFSimulationTime = 10000000;//10000000;
static const int NoSyns	= 1;//Switch rule Ensemble's number of Synapses


//THESE PARAMETERS ONY AFFECT SWITCH RULE
static const float APOT	= 1.0f;///(iTestFq*200); //(iTestFq*IFSimulationTime*h); Too Slow
static const float ADEP	= 0.95f;///(iTestFq*200); //(iTestFq*IFSimulationTime*h); Too Slow
static const float tafPOT	= 0.013f;
static const float tafDEP	= 0.020f;
static const int nPOT	= 1; //Change to 1 for Poisson Neuron Test
static const int nDEP	= 1;
/// \todo
static const float StartStrength = 20.0f;
static const float fKCOscPeriod     = 10.0f; //Period of KC input Neuron Oscillation
static const float fKCBaselineFq    = 20.0f; //Period of KC input Neuron Oscillation

static const string strPlotCmd = "gnuplot SpikeRaster.gplot";

int main(int argc, char *argv[])
{


    ///Record Synapse Strengths

    getcwd(FilePath, _MAX_PATH); // reads the current working directory into the array FilePath
    if (chdir(DATDIR)) {
            perror("chdir to " DATDIR);
        }

    string sfilename(FilePath);
    // Set up output files
    ofiles["SynStrengthLog"] = new std::ofstream(("synsStrength.csv"),ios::out );
    //Iterator ologfile =
    (*ofiles["SynStrengthLog"])  << "#ID Source  Target  Strength" << endl;

    /// Spike Raster
    ofiles["SpikeRasterLog"] = new std::ofstream(("spikeRaster.csv"),ios::out );
    (*ofiles["SpikeRasterLog"])  << "#NeuronID\tSpikeTime" << endl;

    ///ΜΒΟΝ Neuron Membrane Voltage
    ofiles["MBONLog"] = new std::ofstream(("MBONLog.csv"),ios::out );
    (*ofiles["MBONLog"]) << "#t\tVm"  << endl;

    ///KC Neuron Membrane Voltage
    ofiles["KCLog"] = new std::ofstream(("KCLog.csv"),ios::out );
    (*ofiles["KCLog"]) << "#t\tVm"  << endl;

    testIFNeuron(iNoExSynapses,iNoInhSynapses,IFSimulationTime);


   // Close all files
    for(Iterator it = ofiles.begin(); it != ofiles.end(); ++it) {
      delete it->second;
      it->second = 0;
    }

   return 1;
}



void testIFNeuron(int iNoExSynapses,int iNoInhSynapses,uint uiSimulationTime)
{
    int	verboseperiod = 50000;
    int timetolog = verboseperiod;
    uint cnt=0;
    int spikecnt = 0;
    int TotalPostSpikes = 0; //Used to get the Average Post Rate

    double t=0;
    double Vm;
    float gamma = APOT*nPOT*tafPOT/(ADEP*nDEP*tafDEP);

    cout << "---- Test IF Neuron  with Gamma: " << gamma << endl;


    synapseEnsemble<synapseSW,NoSyns> *synstest[iNoExSynapses+iNoInhSynapses];// = new synapseEnsemble[iNoExSynapses];
    PoissonNeuron *pPsKC[iNoExSynapses+iNoInhSynapses];//Create Separate Poisson Sources for each KC afferent

    //Create IFNeuron
    IFNeuron* ifn = new IFNeuron(h,1);
    //Params synapseSW(float A1,float A2,float tafPOT,float tafDEP,int nPOT,int nDEP,float Sreset,bool bNoPlasticity);
    synapseSW osynEx(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,StartStrength,true); //This synapse is going to be copied into the ensemble
    synapseSW osynIn(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,-1,true); //This synapse is going to be copied into the ensemble

    //Create and register Exhitatory Synapses
    for (int i=0;i<iNoExSynapses;i++)
    {

        synapseEnsemble<synapseSW,NoSyns>* posynen = new synapseEnsemble<synapseSW,NoSyns>(h,osynEx);
        synstest[i] = posynen; //new synapseEnsemble<synapseSW,1>(h,osyn); //No Plasticity
        //posynen->RegisterNeuron(ifn);
        ifn->RegisterAfferent((ISynapseEnsemble*)(posynen)); //Let the nEuron Know a bundle of synapses is connecting to it

        //Create New Efferent / start IDs from 5+
        //PoissonNeuron(float timestep,short ID=0,int StartFireRate=0,bool FixedRate=false);
        pPsKC[i] = new PoissonNeuron(h,i+5,iTestFq,false);

    }


    //Open Log file for Synapse Strength
    ///Record Synapse Strengths
    getcwd(FilePath, _MAX_PATH); // reads the current working directory into the array FilePath
    char *buff = strcat(FilePath,"\\synsStrength.csv");

    ofstream ofile(buff, ios::out );

    ofile << "ID Source Target S" << endl;
    //////LOG File Opened

    //Main Simulation Loop - Ends when Simulation time reached
    while (cnt < uiSimulationTime)
    {
        timetolog--;
        t+=h;
        cnt++;
        //Run Through each afferent and Poisson Source
        for (int i=0;i<iNoExSynapses+iNoInhSynapses;i++)
        {
            pPsKC[i]->setFireRate( fKCBaselineFq + 10.0*sin(2*M_PI*t/fKCOscPeriod));
            pPsKC[i]->StepSimTime();

            if (pPsKC[i]->ActionPotentialOccured())
            {
                synstest[i]->SpikeArrived(ISynapse::SPIKE_PRE);
                (*ofiles["SpikeRasterLog"])  << pPsKC[i]->getID() << "\t" << t << endl;
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

        //Log New Membrane Voltage
        (*ofiles["MBONLog"]) << t <<"\t" << Vm  << endl;

        if (ifn->ActionPotentialOccured())
        {
            TotalPostSpikes++;
            (*ofiles["SpikeRasterLog"])  << ifn->getID() << "\t" << t << endl;
        }


        if (timetolog==0)
        {
            cout<< cnt<< " t: "<< t << " Vm:" <<  Vm << " Sj:" << synstest[0]->getAvgStrength() << " F Rate:" << ifn->getFireRate() << " Avg Rate:" << TotalPostSpikes/t <<endl;
            //Log Poisson KC neurons
            //cout << "Ps 1: " << pPsKC[1]->getFireRate() << endl;
            timetolog=verboseperiod;

        }

    }
    cout << endl << " End. Spike Count: " << spikecnt << " Neuron Fire Rate: " << ifn->getFireRate() << " Avg Rate:" << TotalPostSpikes/t;
    cout << " Total time: " << t << " sec" << endl;


    //Plot Output
    system(strPlotCmd.c_str());

    //Close Synapse Log File
    ofile.close();
    delete ifn;
    //delete  &synstest;
    //delete &Ps;


}
