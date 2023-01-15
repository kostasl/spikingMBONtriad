/// \interface ISynapse
/// \brief Empty implementation of virtual functions due to a linker requirement
/// \note need to define the pure virtual functions, otherwise compiler gives "undefined reference to typeinfo for ISynapse"

#include "isynapse.h"

ISynapse::ISynapse()
{

}


ISynapse::ISynapse(const ISynapse& obj)
{
//Need a function body Declared here
}

ISynapse::~ISynapse()
{

}
