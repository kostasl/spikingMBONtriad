/// \interface ISynapse
/// \brief Empty implementation of virtual functions due to a linker requirement
/// \note need to define the pure virtual functions, otherwise compiler gives "undefined reference to typeinfo for ISynapse"

#include "isynapse.h"

ISynapse::ISynapse()
{

}

 /// \brief Copy Constructor
ISynapse::ISynapse(const ISynapse& obj)
{

}

//Returns the Potentiation value
float ISynapse::getStrength()
{

}

 //throws Exception
float ISynapse::SpikeArrived(double t,SPIKE_SITE type)
{

}

 //Reset Strength and State
void ISynapse::Reset()
{

}
