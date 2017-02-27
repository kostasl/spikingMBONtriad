//#include <QCoreApplication>


#include "stdafx.h"
#include "synapse.h"
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


//Function Prototypes

using namespace std;


int main(int argc, char *argv[])
{

    IFNeuron ifn;

    return 1;
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
    _getcwd(FilePath, _MAX_PATH); // reads the current working directory into the array FilePath
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
