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
const double h=TIMESTEP;

static char FilePath[_MAX_PATH]; // _MAX_PATH represents the longest possible path on this OS
typedef std::map<std::string, std::ofstream*> FileMap;
typedef FileMap::iterator Iterator;

//Log Files - Declare map of log files - Used to Organize all output streams
FileMap ofiles;


/// Network Structure Variables
static const uint IFSimulationTime = 1000000;//10000000;
static const int NoSynsWa	= 0;//Switch rule Ensemble's number of Synapses KC->DAN
static const int NoSynsWb	= 1;//Switch rule Ensemble's number of Synapses KCs->MBON
static const int NoSynsWz	= 1;//Switch rule Ensemble's number of Synapses KCs->MBON
static const int NoSynsWd	= 10;//Switch rule Ensemble's number of Synapses  DAN -> MBON
static const int NoSynsWg	= 1;//Switch rule Ensemble's number of Synapses MBON -> DAN

//Kenyon Cells /Input Pattern
static const int iTestFq	= 30;
static const int iNoExSynapses = 10; //Number of synapses to test IFNeuron
static const int iNoInhSynapses = 0;//200; //Inhibitory synapses
static const float fKCOscPeriod     = 150.0f; //Period of KC input Neuron Oscillation
static const float fKCBaselineFq    = 0.0f; //Baseline Spiking Rate of KC input Neuron Ontop Of Which the Oscillating one rides


//THESE PARAMETERS ONY AFFECT SWITCH RULE
static const float APOT	= 1.0f;///(iTestFq*200); //(iTestFq*IFSimulationTime*h); Too Slow
static const float ADEP	= 0.95f;///(iTestFq*200); //(iTestFq*IFSimulationTime*h); Too Slow
static const float tafPOT	= 0.013f;
static const float tafDEP	= 0.020f;
static const int nPOT	= 1; //Change to 1 for Poisson Neuron Test
static const int nDEP	= 1;

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
    testMBONTriadConfigA(iNoExSynapses,iNoInhSynapses,IFSimulationTime);

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
    const float fExSynapseStartStrength = 50;
    const float fInSynapseStartStrength = -50;

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
    synapseSW osynEx(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,fExSynapseStartStrength,true); //This synapse is going to be copied into the ensemble
    synapseSW osynIn(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,fInSynapseStartStrength,true); //This synapse is going to be copied into the ensemble

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
            pPsKC[i]->setFireRate( fKCBaselineFq + 0.0*sin(2*M_PI*t/fKCOscPeriod));
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
/// \note Adding a new neuron to the network requires its instantiation (either Poisson or IF), then a make a synapseensemble to connect the two.
/// then call Registering Efferent at source neuron A and call RegisterEfferent on target neuron B
void testMBONTriadConfigA(int iNoExSynapses,int iNoInhSynapses,uint uiSimulationTime)
{
    int RinputFq        = 60;
    int	verboseperiod   = 50000; //Report to Our every 5 secs
    int timetolog       = verboseperiod;
    uint cnt            = 0; //Current Timestep
    int spikecnt        = 0; //input Spikes Count
    int TotalPostSpikes = 0; //Output (MBON)Used to get the Average Post Rate

    double t            = 0; //Simulated Time in seconds

    float gamma         = APOT*nPOT*tafPOT/(ADEP*nDEP*tafDEP); //This is Switch rule Plasticity Specific parameter- Not used here

    /// Synaptic Parameters  - No Plasticity - Just define Strentgh
    //Define the synapse types and their parameters that would be used in the ensembles connecting the neurons
    //Params synapseSW(float A1,float A2,float tafPOT,float tafDEP,int nPOT,int nDEP,float Sreset,bool bNoPlasticity);
    synapseSW osynWa(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,10,true); //This synapse is going to be copied into the ensemble
    synapseSW osynWb(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,20,true); //This synapse is going to be copied into the ensemble
    synapseSW osynWg(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,-10,true); //This synapse is going to be copied into the ensemble
    synapseSW osynWd(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,10,true); //This synapse is going to be copied into the ensemble
    synapseSW osynWz(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,-10,true); //This synapse is going to be copied into the ensemble
    synapseSW osynEx(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,50,true); //This synapse is going to be copied into the ensemble
    synapseSW osynIn(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,-100.0,true); // fInSynapseStartStrength//This synapse is going to be copied into the ensemble


    ///Make Synapses of the Triad KC-DAN-MBON
    //Synapse Array // Holding on to pointers for Reporting Reasons
    synapseEnsemble<synapseSW,NoSynsWa> *arrKCsynsWa[iNoExSynapses+iNoInhSynapses];//KC->DAN Synapses
    synapseEnsemble<synapseSW,NoSynsWb> *arrKCsynsWb[iNoExSynapses+iNoInhSynapses];//KC->MBON Synapses
    synapseEnsemble<synapseSW,NoSynsWz> *arrKCsynsWz[iNoExSynapses+iNoInhSynapses];//KC->MBON Synapses
    synapseEnsemble<synapseSW,NoSynsWg> *arrMBONsynsWg[1];//MBON->DAN Synapses
    synapseEnsemble<synapseSW,NoSynsWd> *arrMBONsynsWd[1];//DAN->MBON Synapses





    cout << "---- Test MBON Triad KCs->DAN<->MBON Neuron  with Synapse Switch rule Gamma: " << gamma << endl;
    cout << "Rate :" << IFFIRERATE_PERIOD << " timesteps/sec " << endl;


    //Instantiate Network Neurons -KCs, DAN MBON
    PoissonNeuron *pPsKC[iNoExSynapses+iNoInhSynapses];//Create Separate Poisson Sources for each KC afferent

    //Create IFNeuron for MBON
    IFNeuron* pIfnMBON = new IFNeuron(h,1); //MBON ID -> 1
    //DAN
    IFNeuron* pIfnDAN = new IFNeuron(h,2); //DAN ID ->2

    // (R) Add Reward Input to DAN
    PoissonNeuron* pPsR = new  PoissonNeuron(h,3,RinputFq,true);
    synapseEnsemble<synapseSW,5>* osynR = new synapseEnsemble<synapseSW,5>(h,osynEx);
    pPsR->RegisterEfferent((ISynapseEnsemble*)osynR); //Target Neuron is The DAN
    pIfnDAN->RegisterAfferent((ISynapseEnsemble*)osynR); //Register as outgoing Synapse


    /// Make Synaptic Connections
    // W betaGenerate KC poisson Neurons - Create and register Exhitatory Synapses to MBON and to DAN
    for (int i=0;i<iNoExSynapses;i++)
    {
        //W_beta synapseEnsemble<synapseSW,NoSyns>* posynWb
        arrKCsynsWb[i] = new synapseEnsemble<synapseSW,NoSynsWb>(h,osynWb); //new synapseEnsemble<synapseSW,1>(h,osyn); //No Plasticity
        //Connect Synapse to Neuron - Let Neuron know that this synapses is connecting to it - This will also let the synapse know of the target neuron by calling posynen->RegisterNeuron(ifn);
        pIfnMBON->RegisterAfferent((ISynapseEnsemble*)( arrKCsynsWb[i])); //Let the nEuron Know a bundle of synapses is connecting to it

        //W_alpha KC -> DAN Connections
        arrKCsynsWa[i] = new synapseEnsemble<synapseSW,NoSynsWa>(h,osynWa);
        pIfnDAN->RegisterAfferent((ISynapseEnsemble*)arrKCsynsWa[i]);

        //W_zeta DAN -> KC Connections
        arrKCsynsWz[i] = new synapseEnsemble<synapseSW,NoSynsWz>(h,osynWz);
        pIfnDAN->RegisterEfferent((ISynapseEnsemble*)arrKCsynsWz[i]); //Connect DAN->KC


        //Create New Efferent / start IDs from 5+
        //PoissonNeuron(float timestep,short ID=0,int StartFireRate=0,bool FixedRate=false);
        pPsKC[i] = new PoissonNeuron(h,i+5,iTestFq,false);
        pPsKC[i]->RegisterEfferent((ISynapseEnsemble*)( arrKCsynsWb[i])); //Tell this neuron it has a target
        pPsKC[i]->RegisterEfferent((ISynapseEnsemble*)( arrKCsynsWa[i])); //Tell this neuron it has a target
        pPsKC[i]->RegisterAfferent((ISynapseEnsemble*)( arrKCsynsWz[i])); //Tell this neuron it has a target
    }

    // /W_delta Add DAN ->MBON Synapse
    arrMBONsynsWd[0] = new synapseEnsemble<synapseSW,NoSynsWd>(h,osynWd); //new synapseEnsemble<synapseSW,1>(h,osyn); //No Plasticity
    pIfnDAN->RegisterEfferent((ISynapseEnsemble*)arrMBONsynsWd[0]); //Register as outgoing Synapse
    pIfnMBON->RegisterAfferent((ISynapseEnsemble*)(arrMBONsynsWd[0])); //Register as Incoming Synapse

    // W_gamma - feedback MBON ->DAN
    arrMBONsynsWg[0] = new synapseEnsemble<synapseSW,NoSynsWg>(h,osynWg); //inhibitory
    pIfnDAN->RegisterAfferent((ISynapseEnsemble*)arrMBONsynsWg[0]); //Register as Incoming Synapse
    pIfnMBON->RegisterEfferent((ISynapseEnsemble*)(arrMBONsynsWg[0])); //Register as Ougoing Synapse


    /// Main Simulation Loop - Calls stepSimTime On each neuron - Ends when Simulation time reached
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
        pPsR->StepSimTime();
        pIfnMBON->StepSimTime();
        pIfnDAN->StepSimTime();


        ///Log - Report Output
        //Log New Membrane Voltages
        (*ofiles["MBONLog"]) << t <<"\t" << pIfnMBON->getMembraneVoltage() << "\t" << pIfnMBON->getFireRate()  << endl;
        (*ofiles["DANLog"]) << t <<"\t" << pIfnDAN->getMembraneVoltage() << "\t" << pIfnDAN->getFireRate()  << endl;
        (*ofiles["KCLog"]) << t <<"\t 0 \t" << pPsKC[0]->getFireRate()  << endl;

        //Log Spikes
        if (pIfnMBON->ActionPotentialOccured())
        {
            TotalPostSpikes++;
            (*ofiles["SpikeRasterLog"])  << pIfnMBON->getID() << "\t" << t << endl;
        }
        if (pIfnDAN->ActionPotentialOccured())
        {
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


    //Free up memory
    delete pIfnMBON;
    delete pIfnDAN;
    delete pPsR;

    //Close Synapse Log File

    for (int i=0;i<(iNoExSynapses+iNoInhSynapses);i++)
        delete  pPsKC[i];

} // END Function IFNeuron
