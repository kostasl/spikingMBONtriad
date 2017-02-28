#ifndef ISYNAPSE_H
#define ISYNAPSE_H



class ISynapse
{
public:
    enum SPIKE_SITE {SPIKE_PRE, SPIKE_POST};

    ISynapse();
    ISynapse(const ISynapse& obj); //Copy Constructor
    virtual float getStrength(); //Returns the Potentiation value
    virtual float SpikeArrived(double t,SPIKE_SITE type); //throws Exception
    virtual void Reset(); //Reset Strength and State
};

#endif // ISYNAPSE_H
