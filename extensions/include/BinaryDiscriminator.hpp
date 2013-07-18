#ifndef BINARY_DISCRIMINATOR_ABSTRACT_BASE_EXTERNAL
#define BINARY_DISCRIMINATOR_ABSTRACT_BASE_EXTERNAL

#include <Rtypes.h>

class BinaryDiscriminator
{
	public:
		virtual Double_t operator()(Double_t const *) const = 0;
};
#endif
