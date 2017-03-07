///* \class SYNAPSE ENSEMBLE
/// \brief The collection of single synapses. The Ensemble propagates spike events to synapses,
/// and holds the time since the last spike event.
///  \note On every simulation time tick, either  SpikeArrived or StepNoSpike() must be called.
/// \author Kostas Lagogiannis 23/6/07
///

#include "stdafx.h"
#include "math.h"
#include "isynapse.h"
#include "synapseSW.h"
#include "IFNeuron.h"
#include "PoissonNeuron.h"
#include "INeuron.h"
#include "synapticTransmission.h"
#include "synapseEnsemble.h"
#include "isynapseensemble.h"

//Constructor Not used
template<class T,int N>
synapseEnsemble<T,N>::synapseEnsemble(float simTimeStep,short sourceID,short ID):ISynapseEnsemble()
{
    mID                     = ID;
    mbNoPlasticity          = true; //Switch Off Plasticity By Default
    mfSimTimestep           = simTimeStep;

    mdTimeSinceLastSpike    = 0;//Reset Time since last spike
    //mSourceFireRate         = 0;
//    msourceID               = sourceID;
//    mtargetID               = 0;
//    mpTargetNeuron          = 0; //Init Pointer to 0
//    bSpikeOccured           = false;
}


/// \brief Propage the spike to all synapse members, reset the time since last spike
/// \returns Change in Ds in strength
/// \note Template implementations need to be in the header file in order to be found by the calling classes
/// however this is placed here since it calls on method of INeuron, which at the time of synapseEnsemble header compilation the declaration of INeuron are not known
/// So I am forced to put this here and force particular synapseEnsemble Template instatiations for each Synapse Type at the end of this cpp file
///

template<class T,int N> float synapseEnsemble<T,N>::SpikeArrived(ISynapse::SPIKE_SITE type)
{
       double totalWeight=0.0;

       ///Only obtain the Average When using Switch Rule
#ifndef USE_SONG_LEARNING
       //NOTE: Probably need to ignore PRE spikes when in State POT
       //And POST spikes when in DEP
       mfAvgStrength = 0;
       for(int i=0;i<miSynapsesCount;i++)
       {
           //Spike Arrived Returns the Change in Weight. Need to Sum to Obtain Aggregate
           totalWeight += arrSynapses[i].SpikeArrived(mdTimeSinceLastSpike,type);
           mfAvgStrength+=arrSynapses[i].getStrength();
       }
       //If PLasticity is Off The strength is not Saved in Synapses So AvgStrength=0;
       //Save & Return Average Strength
       //if (!mbNoPlasticity) ///If running in Plasticity
       //mfAvgStrength = totalWeight; // No Averaging  ///miSynapsesCount;
#endif
       //Reset Time since last spike
       mdTimeSinceLastSpike = 0.0;
       bSpikeOccured        = true;

       //Propagate Spike to Neuron if Target exists
       if (mpTargetNeuron && type == ISynapse::SPIKE_SITE::SPIKE_PRE)
       {
           synapticTransmission* sp = new synapticTransmission((float)mfSimTimestep,mfAvgStrength);
           //Debug
           //if (mfAvgStrength < 0.0)
//           {
  //             if (mpTargetNeuron != 0) mpTargetNeuron->SpikeArrived(sp); //Spike Is passed to Neuron
    //       }
      //     else
           {
                if (mpTargetNeuron != 0) mpTargetNeuron->SpikeArrived(sp); //Spike Is passed to Neuron
           }


       }

#ifdef USE_SONG_LEARNING


       if (type == synapse::SPIKE_SITE::SPIKE_PRE)
       { //Presynaptic Event Update P
           songPot += Aspot;
           gex_a = gex_a + songMdepot*G_MAX; //Decrease Synaptic Strength

           if (gex_a <0) gex_a = 0;

       }
       else//Postsynaptic Event
       {
           songMdepot -= Asdep; //Increase potential Decrease (M(t)in synaptic Strength
           gex_a = gex_a + songPot*G_MAX; //Increase Synaptic Strength

           if (gex_a > G_MAX)
               gex_a = G_MAX; //IMpose limits 0>gex<gmax
       }

       if (!mbNoPlasticity)
           mfAvgStrength = gex_a;

       //Do time Step Decay Modification of Learning State
       UpdateSongModFunctions(); //Testing if order significant

#endif

       return (float)totalWeight;  //Return Change in s //OLD:Changed from AvgStrength (Always return New strength)
}

template<class T,int N>
double synapseEnsemble<T,N>::getAvgStrength()
{
	return mfAvgStrength;
}

template<class T,int N>
int  synapseEnsemble<T,N>::getSynapsesCount()
{
	return miSynapsesCount;
}


template<class T,int N>
short synapseEnsemble<T,N>::getsourceID()
{
	return msourceID;
}
template<class T,int N>
short synapseEnsemble<T,N>::gettargetID()
{
	return mtargetID;
}
//Logs Synapse State to open File stream
template<class T,int N>
void synapseEnsemble<T,N>::logtofile(std::ofstream &fs)
{
	fs << mID << " " << msourceID << " " <<  mtargetID << " " << mfAvgStrength << std::endl;
}

//Deprecated
//template<class T,int N>
//unsigned short synapseEnsemble<T,N>::getSourceFireRate()
//{
//	return mSourceFireRate;
//}

template<class T,int N>
ISynapse* synapseEnsemble<T,N>::getSynapseArray()
{
 return arrSynapses;
}

template<class T,int N>
synapseEnsemble<T,N>::~synapseEnsemble(void)
{
    //delete [] arrSynapses throws Seg fault!

    //for(int i=0;i<miSynapsesCount;i++)
    //    if (&arrSynapses[i])
            //delete &arrSynapses[i];
}



/// \brief Specific Template class instantiations required
template class synapseEnsemble<synapseSW,0>;
template class synapseEnsemble<synapseSW,1>;
template class synapseEnsemble<synapseSW,2>;
template class synapseEnsemble<synapseSW,5>;
template class synapseEnsemble<synapseSW,10>;
template class synapseEnsemble<synapseSW,15>;
template class synapseEnsemble<synapseSW,20>;
template class synapseEnsemble<synapseSW,30>;
template class synapseEnsemble<synapseSW,50>;
template class synapseEnsemble<synapseSW,80>;
template class synapseEnsemble<synapseSW,100>;
