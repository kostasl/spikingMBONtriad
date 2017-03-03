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

//declare map of log files - Used to Organize all output streams
FileMap ofiles;


//Variables for IFTEST
static const int iTestFq	= 30;
static const int iNoExSynapses = 10; //Number of synapses to test IFNeuron
static const int iNoInhSynapses = 0;//200; //Inhibitory synapses
static const uint IFSimulationTime = 5000000;//10000000;
static const int NoSynsWa	= 1;//Switch rule Ensemble's number of Synapses KC->DAN
static const int NoSynsWb	= 1;//Switch rule Ensemble's number of Synapses KCs->MBON
static const int NoSynsWd	= 1;//Switch rule Ensemble's number of Synapses  DAN -> MBON
static const int NoSynsWg	= 1;//Switch rule Ensemble's number of Synapses MBON -> DAN


//THESE PARAMETERS ONY AFFECT SWITCH RULE
static const float APOT	= 1.0f;///(iTestFq*200); //(iTestFq*IFSimulationTime*h); Too Slow
static const float ADEP	= 0.95f;///(iTestFq*200); //(iTestFq*IFSimulationTime*h); Too Slow
static const float tafPOT	= 0.013f;
static const float tafDEP	= 0.020f;
static const int nPOT	= 1; //Change to 1 for Poisson Neuron Test
static const int nDEP	= 1;
/// \todo
static const float StartStrength = 10.0f;
static const float fKCOscPeriod     = 150.0f; //Period of KC input Neuron Oscillation
static const float fKCBaselineFq    = 20.0f; //Baseline Spiking Rate of KC input Neuron Ontop Of Which the Oscillating one rides

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
    (*ofiles["MBONLog"]) << "#t\tVm\tSpikeRate"  << endl;

    ///DAN Neuron Membrane Voltage
    ofiles["DANLog"] = new std::ofstream(("DANLog.csv"),ios::out );
    (*ofiles["DANLog"]) << "#t\tVm\tSpikeRate"  << endl;


    ///KC Neuron Membrane Voltage
    ofiles["KCLog"] = new std::ofstream(("KCLog.csv"),ios::out );
    (*ofiles["KCLog"]) << "#t\tVm\tSpikeRate"  << endl;

    //testIFNeuron(iNoExSynapses,iNoInhSynapses,IFSimulationTime);
    testMBONTriad(iNoExSynapses,iNoInhSynapses,IFSimulationTime);

   // Close all files
    for(Iterator it = ofiles.begin(); it != ofiles.end(); ++it) {
      delete it->second;
      it->second = 0;
    }

   return 1;
}


/// \brief Instantiates N poisson excitatory neurons feeding onto an IF neuron , with no plasticity
/// \returns Exports csv files of a spike raster and a csv with the IFs Membrane voltage - Calls gnuplot and generates
/// plots in the dat subfolder
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


    //Make Synapses
    //Synapse Array
    synapseEnsemble<synapseSW,NoSynsWb> *arrKCsynsWb[iNoExSynapses+iNoInhSynapses];//KC->MBON Synapses
    PoissonNeuron *pPsKC[iNoExSynapses+iNoInhSynapses];//Create Separate Poisson Sources for each KC afferent

    //Create IFNeuron
    IFNeuron* ifn = new IFNeuron(h,1);
    //Params synapseSW(float A1,float A2,float tafPOT,float tafDEP,int nPOT,int nDEP,float Sreset,bool bNoPlasticity);
    synapseSW osynEx(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,StartStrength,true); //This synapse is going to be copied into the ensemble
    synapseSW osynIn(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,-1,true); //This synapse is going to be copied into the ensemble

    //Create and register Exhitatory Synapses
    for (int i=0;i<iNoExSynapses;i++)
    {
        //synapseEnsemble<synapseSW,NoSyns>* posynWb
        arrKCsynsWb[i] = new synapseEnsemble<synapseSW,NoSynsWb>(h,osynEx); //new synapseEnsemble<synapseSW,1>(h,osyn); //No Plasticity
        //posynen->RegisterNeuron(ifn);
        ifn->RegisterAfferent((ISynapseEnsemble*)( arrKCsynsWb[i])); //Let the nEuron Know a bundle of synapses is connecting to it

        //Create New Efferent / start IDs from 5+
        //PoissonNeuron(float timestep,short ID=0,int StartFireRate=0,bool FixedRate=false);
        pPsKC[i] = new PoissonNeuron(h,i+5,iTestFq,false);
        pPsKC[i]->RegisterEfferent((ISynapseEnsemble*)( arrKCsynsWb[i])); //Tell this neuron it has a target

    }


    //Main Simulation Loop - Ends when Simulation time reached
    //Need to Call StepSim Time on Each Neuron In the Network
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
            //Log Spike
            if (pPsKC[i]->ActionPotentialOccured())
            {
                (*ofiles["SpikeRasterLog"])  << pPsKC[i]->getID() << "\t" << t << endl;
                spikecnt++;
            }


            //Log Strength of Every Exhitatory Synapse
            //if (timetolog==0 && i < iNoExSynapses) synstest[i]->logtofile(ofile);
        }
        ifn->StepSimTime();
        Vm = ifn->getMembraneVoltage();

        //Log New Membrane Voltage
        (*ofiles["MBONLog"]) << t <<"\t" << Vm  << endl;
        //Log Spike
        if (ifn->ActionPotentialOccured())
        {
            TotalPostSpikes++;
            (*ofiles["SpikeRasterLog"])  << ifn->getID() << "\t" << t << endl;
        }


        if (timetolog==0)
        {
            cout<< cnt<< " t: "<< t << " Vm:" <<  Vm << " Sj:" << arrKCsynsWb[0]->getAvgStrength() << " F Rate:" << ifn->getFireRate() << " Avg Rate:" << TotalPostSpikes/t <<endl;
            //Log Poisson KC neurons
            //cout << "Ps 1: " << pPsKC[1]->getFireRate() << endl;
            timetolog=verboseperiod;

        }

    }
    cout << endl << " End. Spike Count: " << spikecnt << " Neuron Fire Rate: " << ifn->getFireRate() << " Avg Rate:" << TotalPostSpikes/t;
    cout << " Total time: " << t << " sec" << endl;


    //Plot Output
    int ret = system(strPlotCmd.c_str());
    if (ret >0)
        cout << "\n GnuPlots Complete, check dat subfolder!)" << endl;

    //Close Synapse Log File
    delete ifn;

    for (int i=0;i<(iNoExSynapses+iNoInhSynapses);i++)
        delete  pPsKC[i];


} // END Function IFNeuron


/// \name testMBONTriad
/// \brief Instantiates N poisson excitatory neurons representing the input KCs
/// - These feeding onto an IF MBON neurons as well as onto a IF. DAN
/// \returns spike raster (MBON ID = 1, DAN =2, KCs have >5) of whole network, and membrane voltages of the 2 IF neurons
/// There is no no plasticity activated.
void testMBONTriad(int iNoExSynapses,int iNoInhSynapses,uint uiSimulationTime)
{
    int	verboseperiod = 50000; //Report to Our every 5 secs
    int timetolog = verboseperiod;
    uint cnt     = 0; //Current Timestep
    int spikecnt = 0; //input Spikes Count
    int TotalPostSpikes = 0; //Output (MBON)Used to get the Average Post Rate

    double t            = 0; //Simulated Time in seconds

    float gamma = APOT*nPOT*tafPOT/(ADEP*nDEP*tafDEP); //This is Switch rule Plasticity Specific parameter- Not used here

    cout << "---- Test MBON Triad KCs->DAN<->MBON Neuron  with Synapse Switch rule Gamma: " << gamma << endl;


    ///Make Synapses of the Triad KC-DAN-MBON
    //Synapse Array // Holding on to pointers for Reporting Reasons
    synapseEnsemble<synapseSW,NoSynsWa> *arrKCsynsWa[iNoExSynapses+iNoInhSynapses];//KC->DAN Synapses
    synapseEnsemble<synapseSW,NoSynsWb> *arrKCsynsWb[iNoExSynapses+iNoInhSynapses];//KC->MBON Synapses
    synapseEnsemble<synapseSW,NoSynsWg> *arrMBONsynsWg[1];//MBON->DAN Synapses
    synapseEnsemble<synapseSW,NoSynsWd> *arrMBONsynsWd[1];//DAN->MBON Synapses


    //Define the synapse types and their parameters that would be used in the ensembles connecting the neurons
    //Params synapseSW(float A1,float A2,float tafPOT,float tafDEP,int nPOT,int nDEP,float Sreset,bool bNoPlasticity);
    synapseSW osynEx(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,StartStrength,true); //This synapse is going to be copied into the ensemble
    synapseSW osynIn(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,-600.0,true); //This synapse is going to be copied into the ensemble

    //Instantiate Network Neurons -KCs, DAN MBON
    PoissonNeuron *pPsKC[iNoExSynapses+iNoInhSynapses];//Create Separate Poisson Sources for each KC afferent

    //Create IFNeuron for MBON
    IFNeuron* pIfnMBON = new IFNeuron(h,1); //MBON ID -> 1
    //DAN
    IFNeuron* pIfnDAN = new IFNeuron(h,2); //DAN ID ->2

    /// Make Synaptic Connections
    // Generate KC poisson Neurons - Create and register Exhitatory Synapses to MBON and to DAN
    for (int i=0;i<iNoExSynapses;i++)
    {
        //synapseEnsemble<synapseSW,NoSyns>* posynWb
        arrKCsynsWb[i] = new synapseEnsemble<synapseSW,NoSynsWb>(h,osynEx); //new synapseEnsemble<synapseSW,1>(h,osyn); //No Plasticity
        //Connect Synapse to Neuron - Let Neuron know that this synapses is connecting to it - This will also let the synapse know of the target neuron by calling posynen->RegisterNeuron(ifn);
        pIfnMBON->RegisterAfferent((ISynapseEnsemble*)( arrKCsynsWb[i])); //Let the nEuron Know a bundle of synapses is connecting to it

        //KC -> DAN Connections
        arrKCsynsWa[i] = new synapseEnsemble<synapseSW,NoSynsWa>(h,osynEx);
        pIfnDAN->RegisterAfferent((ISynapseEnsemble*)arrKCsynsWa[i]);

        //Create New Efferent / start IDs from 5+
        //PoissonNeuron(float timestep,short ID=0,int StartFireRate=0,bool FixedRate=false);
        pPsKC[i] = new PoissonNeuron(h,i+5,iTestFq,false);
        pPsKC[i]->RegisterEfferent((ISynapseEnsemble*)( arrKCsynsWb[i])); //Tell this neuron it has a target
        pPsKC[i]->RegisterEfferent((ISynapseEnsemble*)( arrKCsynsWa[i])); //Tell this neuron it has a target


    }


    //Add DAN ->MBON Synapse
    arrMBONsynsWd[0] = new synapseEnsemble<synapseSW,NoSynsWd>(h,osynIn); //new synapseEnsemble<synapseSW,1>(h,osyn); //No Plasticity
    pIfnDAN->RegisterEfferent((ISynapseEnsemble*)arrMBONsynsWd[0]); //Register as outgoing Synapse
    pIfnMBON->RegisterAfferent((ISynapseEnsemble*)(arrMBONsynsWd[0])); //Register as Incoming Synapse

    /// Main Simulation Loop - Ends when Simulation time reached
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
                (*ofiles["SpikeRasterLog"])  << pPsKC[i]->getID() << "\t" << t << endl;
                spikecnt++;
            }


            //Log Strength of Every Exhitatory Synapse
            //if (timetolog==0 && i < iNoExSynapses) synstest[i]->logtofile(ofile);
        }
        pIfnMBON->StepSimTime();
        pIfnDAN->StepSimTime();


        ///Log - Report Output
        //Log New Membrane Voltages
        (*ofiles["MBONLog"]) << t <<"\t" << pIfnMBON->getMembraneVoltage() << "\t" << pIfnMBON->getFireRate()  << endl;
        (*ofiles["DANLog"]) << t <<"\t" << pIfnDAN->getMembraneVoltage() << "\t" << pIfnDAN->getFireRate()  << endl;


        //Log Spikes
        if (pIfnMBON->ActionPotentialOccured())
        {
            TotalPostSpikes++;
            (*ofiles["SpikeRasterLog"])  << pIfnMBON->getID() << "\t" << t << endl;
        }
        if (pIfnDAN->ActionPotentialOccured())
        {
            TotalPostSpikes++;
            (*ofiles["SpikeRasterLog"])  << pIfnDAN->getID() << "\t" << t << endl;
        }


        //Periodic Reporting to Std IO
        if (timetolog==0)
        {
            cout<< " t:"<< t << " Vm:" << (int)(pIfnMBON->getMembraneVoltage()*1000) << "mV Sj:" << arrKCsynsWb[0]->getAvgStrength() << " F Rate:" << pIfnMBON->getFireRate() << " Avg Rate:" << TotalPostSpikes/t <<endl;
            //Log Poisson KC neurons
            //cout << "Ps 1: " << pPsKC[1]->getFireRate() << endl;
            timetolog=verboseperiod;

        }

    } //End Main Sim Loop
    cout << "~~~~~~~~End Total time: " << t << " sec. ~~~~~~~~~~~~~" << endl;
    cout << endl << "~ Input Spike Count: " << spikecnt << " MBON Fire Rate: " << pIfnMBON->getFireRate() << " Avg. Rate:" << TotalPostSpikes/t;



    //Plot Output
    int iret = system(strPlotCmd.c_str());
    cout << endl << "Gnuplot Returned :" << iret << endl;

    //Close Synapse Log File

    delete pIfnMBON;
    delete pIfnDAN;

    for (int i=0;i<(iNoExSynapses+iNoInhSynapses);i++)
        delete  pPsKC[i];

} // END Function IFNeuron
