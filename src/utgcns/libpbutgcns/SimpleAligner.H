#ifndef __GCON_SIMPLE_ALIGNER__
#define __GCON_SIMPLE_ALIGNER__

#include "dw.H"

class SimpleAligner {
public:
    SimpleAligner();
    void align(dagcon::Alignment& aln, double errorRate) { align(aln, errorRate, true); };
    void align(dagcon::Alignment &aln, double errorRate, bool useND); 
    void alignND(dagcon::Alignment &aln, double errorRate);
    void alignSW(dagcon::Alignment &aln, double errorRate);
};

#endif // __GCON_SIMPLE_ALIGNER__
