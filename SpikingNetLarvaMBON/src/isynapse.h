#ifndef ISYNAPSE_H
#define ISYNAPSE_H



class ISynapse
{
public:
    enum SPIKE_SITE {SPIKE_PRE, SPIKE_POST};

    ISynapse();
    ISynapse(const ISynapse& obj); //Copy Constructor

    virtual float getStrength()=0; //Returns the Potentiation value
    virtual float SpikeArrived(double t,SPIKE_SITE type)=0; //throws Exception
    virtual void Reset()=0; //Reset Strength and State
    virtual ~ISynapse();
};


#endif // ISYNAPSE_H
