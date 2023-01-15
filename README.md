# Spiking Neural Simulation of Larval Drosophila Learning    
author: Konstantinos Lagogiannis 2016 

# Description
Custom spiking neuron code used as a framework to research and validate hypothesis on neurobiology and behaviour of learning and recall in Drosophila Larvae.

It comprises of a set of C++ classes that implemend different type of Spiking Neuron Models, as well as Spike-Timing Dependent Plasticity rules for learning.
In main.cpp these are connected in a configurations that capture the anatomical details as revealed by electron microscopy in the Mushroom body output network of the larval Drosophila brain.    
The code in STPD_mod.cpp tests/validates Long-Term Potentiation / Depression synaptic learning can be repreoduced by the model.

Neurons implemented: 

    * IFNeuron : Integrate and Fire Neuron
    * CFSNeuron : Izhikevich Fast spiking  neuron
    * CRSNeuron : Izhikevich Regular Spiking  neuron 
    * PoissonSource : A Poisson type spiking neuron

Synaptic learning rules:

    * synapseSW :Implements the individual Synaptic switch as described in the Appleby & Elliot 2005 STPD paper
    * synapseEnsemble collection of individual SW synapses. The Ensemble propagates spike events to synapses,
    * synapticTransmission :  Handles the dynamics of synaptic transmission due to a pre synaptic spike event. 



<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.
