///*
/// \brief Object oriented research code implementing spiking model of the Larval MBON triad
/// Here I use an old implementation of a Leaky integrate and fire neuron which implements a 4th Order Runge-Kutta method.
///
///
/// \todo Separate Synapse Plasticity Models, between SONG and Switch rule - Making independent models
/// \todo Change Array Handling to vector in Neuron and Synapse classes
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
static const int NoSynsWa	= 1;//Switch rule Ensemble's number of Synapses KC->DAN
static const int NoSynsWb	= 1;//Switch rule Ensemble's number of Synapses KCs->MBON
static const int NoSynsWz	= 1;//Switch rule Ensemble's number of Synapses KCs->MBON
static const int NoSynsWd	= 1;//Switch rule Ensemble's number of Synapses  DAN -> MBON
static const int NoSynsWg	= 1;//Switch rule Ensemble's number of Synapses MBON -> DAN

//Kenyon Cells /Input Pattern
static const int iTestFq	= 30;
static const int iNoExSynapses = 10; //Number of synapses to test IFNeuron
static const int iNoInhSynapses = 0;//200; //Inhibitory synapses
static const float fKCOscPeriod     = 150.0f; //Period of KC input Neuron Oscillation
static const float fKCOscAmplitude  = 20.0f;
static const float fKCBaselineFq    = 40.0f; //Baseline Spiking Rate of KC input Neuron Ontop Of Which the Oscillating one rides


//THESE PARAMETERS ONY AFFECT SWITCH RULE
static const float APOT	= 1.0f;///(iTestFq*200); //(iTestFq*IFSimulationTime*h); Too Slow
static const float ADEP	= 0.95f;///(iTestFq*200); //(iTestFq*IFSimulationTime*h); Too Slow
static const float tafPOT	= 0.013f;
static const float tafDEP	= 0.020f;
static const int nPOT	= 1; //Change to 1 for Poisson Neuron Test
static const int nDEP	= 1;

static const string strPlotCmd      = "gnuplot SpikeRaster.gplot";
static const string strShowPlot     = "evince NeuronRates.eps";

///
/// \brief main - process input params, run loop across R input strength for two circuit configurations - a baseline and a target one.
/// It compares baseline to target by producing individual plots of each but also a ratio of responses
/// \param argc
/// \param argv
/// \return exit
///
int main(int argc, char *argv[])
{
    const int nplotSteps    = 10;
    //Network Configs           Wa , Wb, Wg, Wdelta, Wzeta,Wex, Winh
    int iNaiveWeights[]     = {100, 100,   0,   100, -100 ,  20, -100};
    int iNaiveWeightssol2[] = {100, 100,-100,   100, -100 ,  20, -100};
    int iPairedWeights[]    = {100, 0  ,   0,   200, -100 ,  20, -100}; //Should show Gating

    // Wa-> 1,Wb -> 0,  Wg -> -1 ,Wd -> 2,  Wj -> -1//
    int iPairedWeights2[]   = {100, 0 , -100,  200, -100 , 20, -100}; //Should show Gating 1/2 response

    int iUnpairedWeights[]  = {10,  30,   0,   -100, -10  , 20, -100};

    int* iWeights          = iPairedWeights;
    string strTag = "Paired";

    //Kenyon Cells /Input Pattern
    float fRewardInputFq        = 80.0f; //Frequency of The R Input To The DAN
    float fKCOscPeriod          = 10.0f; //INput to KC: Period of Slow input-Neuron Oscillation
    float fKCOscAmplitude       = 20.0f; //INput to KC - Amplitude of slow input Oscillation
    float fKCBaselineFq         = 40.0f; //Baseline Spiking Rate of KC input Neuron Ontop Of Which the Oscillating one rides
    const int iInputCount       = 10;

    bool bflagPlot = true; //Produce Output plot
    bool bflagShow = true; //Show / display output file via doc. viewer

    int opt;
    //Process command line options //qmljsdebugger=port:43543,block"
    while ((opt = getopt(argc, argv, "R:P:A:B:p:s:qmljsdebugger=port:43543,block")) != -1) {
            switch (opt) {
            case 'R':
                cout << " R: " << optarg <<  endl <<  endl;
                fRewardInputFq =  atof(optarg);
                break;
            case 'P':
                fKCOscPeriod = atof(optarg);
                break;
            case 'A':
                fKCOscAmplitude = atof(optarg);
                break;
            case 'B':
                fKCBaselineFq = atof(optarg);
                break;
            case 'q':
                cout << " DEBUG MODE " << endl;
                break;
            case 'p':
                bflagPlot = optarg;
                break;
            case 's':
                bflagShow = optarg;
                break;
            case '?':
                fprintf(stderr, "\n\n Usage: %s [-R FqReward] [-P KC Input OscPeriod] [-A Amplitude Of Slow Input Osc.] [-B Input Baseline Fq]\n", argv[0]);
                exit(EXIT_SUCCESS);
            default: /* '?' */
                cout << "Uknown option '" << opt << "'" << endl;

            }
        }

       printf("Reward Input Strength (Hz) [R]=%0.2f;\nKC Input slow Osc. Period [P]=%0.2f;\nAmplitude of KC Slow Input-Oscillation [A]=%0.2f;\nBaseline on of KC input (Hz) [B]=%0.2f\n", fRewardInputFq, fKCOscPeriod, fKCOscAmplitude,fKCBaselineFq);

       if (optind > argc) {
           cout << argc << "<= " << optind << endl;
           fprintf(stderr, "Expected argument after options\n");
           exit(EXIT_FAILURE);
        }

       printf("name argument = %s\n", argv[optind]);



    ///Record Synapse Strengths to file
    getcwd(FilePath, _MAX_PATH); // reads the current working directory into the array FilePath
    if (chdir(DATDIR)) {
            perror("chdir to " DATDIR);
            exit(EXIT_FAILURE);
        }

    //string sfilename(FilePath);
    //stringstream sfilename;

    // Set up output files
    ofiles["SynStrengthLog"] = new std::ofstream(("synsStrength.csv"),ios::out );
    //Iterator ologfile =
    (*ofiles["SynStrengthLog"])  << "#ID Source  Target  Strength" << endl;

    /// Spike Raster
    ofiles["SpikeRasterLog"] = new std::ofstream(("spikeRaster.csv"),ios::out );
    (*ofiles["SpikeRasterLog"])  << "#NeuronID\tSpikeTime" << endl;

    ///ΜΒΟΝ Neuron Membrane Voltage
    ofiles["MBONLog"] = new std::ofstream(("MBONLog.csv"),ios::out );
    (*ofiles["MBONLog"]) << "#t\tR\tVm\tSpikeRate"  << endl;

    ///DAN Neuron Membrane Voltage
    ofiles["DANLog"] = new std::ofstream(("DANLog.csv"),ios::out );
    (*ofiles["DANLog"]) << "#t\tR\tVm\tSpikeRate\tInput_SpikeRate"  << endl;

    ///KC Neuron Membrane Voltage
    ofiles["KCLog"] = new std::ofstream(("KCLog.csv"),ios::out );
    (*ofiles["KCLog"]) << "#t\tER\tVm\tSpikeRate\tInput_SpikeRate"  << endl;


    const uint uiSimulationTime = 100000;
    //testIFNeuron(iNoExSynapses,iNoInhSynapses,IFSimulationTime);

    //scan R input from 0 to Max Reward Input
    for (int r=0;r<nplotSteps;r++)
    {
        testMBONTriadConfigA(iInputCount,iWeights,r*fRewardInputFq/(float)nplotSteps,fKCBaselineFq,fKCOscAmplitude,fKCOscPeriod, uiSimulationTime);

        //Add Blank lines to output / Defines a new Datablock for gnuplot
        //Log New Membrane Voltages
        (*ofiles["MBONLog"])        << endl << endl;
        (*ofiles["DANLog"])         << endl << endl;
        (*ofiles["KCLog"])          << endl << endl;
        (*ofiles["SpikeRasterLog"]) << endl << endl;

    }



    //Plot Output
    int iret;
    if (bflagPlot)
    {
        iret = system(strPlotCmd.c_str());
        cout << endl << "Gnuplot Returned :" << iret << endl;
        if (bflagShow)
        {
            iret = system(strShowPlot.c_str());
            cout << endl << "Show (evince) Returned :" << iret << endl;

        }

        string strPlotExpName= "cp NeuronRates.eps NeuronRates_";
        iret = system(strPlotExpName.append(strTag).append(".eps").c_str());

    }


    string datFilename= "cp MBONLog.csv MBONLog_";
    iret = system(datFilename.append(strTag).append(".csv").c_str());

    datFilename= "cp DANLog.csv DANLog_";
    iret = system(datFilename.append(strTag).append(".csv").c_str());

    datFilename= "cp KCLog.csv KCLog_";
    iret = system(datFilename.append(strTag).append(".csv").c_str());


   // Close all files
    for(Iterator it = ofiles.begin(); it != ofiles.end(); ++it) {
      delete it->second;
      it->second = 0;
    }

   exit(EXIT_SUCCESS);
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
/// \note Adding a new neuron to the network requires its instantiation (either Poisson or IF), then a make a synapseensemble to connect the two.
/// then call Registering Efferent at source neuron A and call RegisterEfferent on target neuron B
void testMBONTriadConfigA(int iInputCount,int * iWeights,float fRewardInputFq,float fKCBaselineFq,float fKCOscAmplitude,float fKCOscPeriod, uint uiSimulationTime)
{



    //iWeights

    const int verboseperiod   = 10000; //Report to Our every 5 secs

    int timetolog       = verboseperiod;
    uint cnt            = 0; //Current Timestep
    int spikecnt        = 0; //input Spikes Count
    int TotalPostSpikes = 0; //Output (MBON)Used to get the Average Post Rate

    double t            = 0; //Simulated Time in seconds
    const float gamma   = APOT*nPOT*tafPOT/(ADEP*nDEP*tafDEP); //This is Switch rule Plasticity Specific parameter- Not used here

    /// Synaptic Parameters  - No Plasticity - Just define Strentgh
    //Define the synapse types and their parameters that would be used in the ensembles connecting the neurons
    //Params synapseSW(float A1,float A2,float tafPOT,float tafDEP,int nPOT,int nDEP,float Sreset,bool bNoPlasticity);

    synapseSW osynWa(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,iWeights[0],true); //This synapse is going to be copied into the ensemble
    synapseSW osynWb(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,iWeights[1],true); //This synapse is going to be copied into the ensemble
    synapseSW osynWg(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,iWeights[2],true); //This synapse is going to be copied into the ensemble
    synapseSW osynWd(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,iWeights[3],true); //This synapse is going to be copied into the ensemble
    synapseSW osynWz(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,iWeights[4],true); //This synapse is going to be copied into the ensemble
    synapseSW osynEx(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,iWeights[5],true); //A general exhitatory Synapse - Used here between R->DAN, & OSNs->KC
    synapseSW osynIn(APOT,ADEP,tafPOT,tafDEP,nPOT,nDEP,iWeights[6],true); // fInSynapseStartStrength ///A general inhibitory Synapse - Used here


    ///Make Synapses of the Triad KC-DAN-MBON
    //Synapse Array // Holding on to pointers for Reporting Reasons
    synapseEnsemble<synapseSW,NoSynsWa> *arrKCsynsWa[iNoExSynapses+iNoInhSynapses];//KC->DAN Synapses
    synapseEnsemble<synapseSW,NoSynsWb> *arrKCsynsWb[iNoExSynapses+iNoInhSynapses];//KC->MBON Synapses
    synapseEnsemble<synapseSW,NoSynsWz> *arrKCsynsWz[iNoExSynapses+iNoInhSynapses];//KC->MBON Synapses
    synapseEnsemble<synapseSW,NoSynsWg> *arrMBONsynsWg[1];//MBON->DAN Synapses
    synapseEnsemble<synapseSW,NoSynsWd> *arrMBONsynsWd[1];//DAN->MBON Synapses





    cout << "---- Test MBON Triad KCs->DAN<->MBON Neuron  with Synapse Switch rule Gamma: " << gamma << endl;
    cout << "Sim Step Rate :" << IFFIRERATE_PERIOD << " timesteps/sec " << endl;
    cout << "Reward Input :" << fRewardInputFq << " Hz." << endl;

    //Instantiate Network Neurons -Inputs to KC, KC, DAN MBON
    PoissonNeuron* pPsOSN[iInputCount];//Create Separate Poisson Sources for each KC afferent

    //Create IFNeuron for MBON
    IFNeuron* pIfnMBON  = new IFNeuron(h,4); //MBON ID -> 1

    //Create IFNeuron for DAN
    IFNeuron* pIfnDAN  = new IFNeuron(h,3); //DAN ID -> 1
    //KC
    IFNeuron* pIfnKC   = new IFNeuron(h,2); //KC ID ->2

    // (R) Add Reward Input to DAN
    PoissonNeuron* pPsR = new  PoissonNeuron(h,1,fRewardInputFq,true);
    synapseEnsemble<synapseSW,1>* osynR = new synapseEnsemble<synapseSW,1>(h,osynEx);
    pPsR->RegisterEfferent((ISynapseEnsemble*)osynR); //Target Neuron is The DAN
    pIfnDAN->RegisterAfferent((ISynapseEnsemble*)osynR); //Register as outgoing Synapse


    // (R) Add Sensory Input to KC

    for (int i=0;i<iInputCount;i++)
    {
        pPsOSN[i] = new  PoissonNeuron(h,6+i,fKCBaselineFq,true);
        synapseEnsemble<synapseSW,1>* osynOSN = new synapseEnsemble<synapseSW,1>(h,osynEx);
        pPsOSN[i]->RegisterEfferent((ISynapseEnsemble*)osynOSN); //Target Neuron is The DAN
        pIfnKC->RegisterAfferent((ISynapseEnsemble*)osynOSN); //Register as outgoing Synapse

    }



    /// Make Synaptic Connections To KC
    //W_beta synapseEnsemble<synapseSW,NoSyns>* posynWb
    arrKCsynsWb[0] = new synapseEnsemble<synapseSW,NoSynsWb>(h,osynWb); //new synapseEnsemble<synapseSW,1>(h,osyn); //No Plasticity
    //Connect Synapse to Neuron - Let Neuron know that this synapses is connecting to it - This will also let the synapse know of the target neuron by calling posynen->RegisterNeuron(ifn);
    pIfnMBON->RegisterAfferent((ISynapseEnsemble*)( arrKCsynsWb[0])); //Let the nEuron Know a bundle of synapses is connecting to it

    //W_alpha KC -> DAN Connections
    arrKCsynsWa[0] = new synapseEnsemble<synapseSW,NoSynsWa>(h,osynWa);
    pIfnDAN->RegisterAfferent((ISynapseEnsemble*)arrKCsynsWa[0]);

    //W_zeta DAN -> KC Connections
    arrKCsynsWz[0] = new synapseEnsemble<synapseSW,NoSynsWz>(h,osynWz);
    pIfnDAN->RegisterEfferent((ISynapseEnsemble*)arrKCsynsWz[0]); //Connect DAN->KC


    //Connect to KC
    pIfnKC->RegisterEfferent((ISynapseEnsemble*)( arrKCsynsWb[0])); //Tell this neuron it has a target
    pIfnKC->RegisterEfferent((ISynapseEnsemble*)( arrKCsynsWa[0])); //Tell this neuron it has a target
    pIfnKC->RegisterAfferent((ISynapseEnsemble*)( arrKCsynsWz[0])); //Tell this neuron it has a target


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
        for (int i=0;i<iInputCount;i++)
        {
            pPsOSN[i]->setFireRate( fKCBaselineFq + fKCOscAmplitude*sin(2.0*M_PI*t/fKCOscPeriod));
            pPsOSN[i]->StepSimTime();

            if (pPsOSN[i]->ActionPotentialOccured())
            {
                (*ofiles["SpikeRasterLog"])  << pPsOSN[i]->getID() << "\t" << t << endl;
                spikecnt++;
            }


            //Log Strength of Every Exhitatory Synapse
            //if (timetolog==0 && i < iNoExSynapses) synstest[i]->logtofile(ofile);
        }

        pPsR->StepSimTime();
        pIfnKC->StepSimTime();
        pIfnMBON->StepSimTime();
        pIfnDAN->StepSimTime();



        ///Log - Report Output
        //Log New Membrane Voltages
        (*ofiles["MBONLog"]) << t << "\t" << fRewardInputFq <<"\t" << pIfnMBON->getMembraneVoltage() << "\t" << pIfnMBON->getFireRate()  << endl;
        (*ofiles["DANLog"])  << t << "\t" << fRewardInputFq <<"\t" << pIfnDAN->getMembraneVoltage() << "\t" << pIfnDAN->getFireRate() << "\t" <<  pPsR->getFireRate()  << endl;
        (*ofiles["KCLog"])   << t << "\t" << fRewardInputFq <<"\t" << pIfnKC->getMembraneVoltage() <<  "\t" << pIfnKC->getFireRate()  << "\t" <<  pPsOSN[0]->getFireRate() << endl;

        //Log Spikes
        if (pIfnMBON->ActionPotentialOccured())
        {
            TotalPostSpikes++;
            (*ofiles["SpikeRasterLog"])  << pIfnMBON->getID() << "\t" << t << endl;
        }
        if (pIfnDAN->ActionPotentialOccured())
            (*ofiles["SpikeRasterLog"])  << pIfnDAN->getID() << "\t" << t << endl;

        if (pIfnKC->ActionPotentialOccured())
            (*ofiles["SpikeRasterLog"])  << pIfnKC->getID() << "\t" << t << endl;


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





    //Free up memory
    delete pIfnMBON;
    delete pIfnDAN;
    delete pIfnKC;
    delete pPsR;

    //Close Synapse Log File

    for (int i=0;i<(iInputCount);i++)
        delete  pPsOSN[i];

} // END Function IFNeuron
