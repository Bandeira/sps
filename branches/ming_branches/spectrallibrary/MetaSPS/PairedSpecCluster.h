/*
 * PairedSpecCluster.h
 *
 *  Created on: Nov 30, 2011
 *      Author: aguthals
 */

#ifndef PAIREDSPECCLUSTER_H_
#define PAIREDSPECCLUSTER_H_

#include "PairedSpecSet.h"

namespace specnets
{
  class PairedSpecCluster: public PairedSpecSet
  {
  public:

    void initialize();

  protected:

    void moveSamePRMs(Spectrum::FragType fragCIDHCD);

    void moveSameSRMs(bool checkPRM, Spectrum::FragType fragCIDHCD);

    void movePRMsSRMs(Spectrum::FragType fragCIDHCD);

    void moveLeftovers(Spectrum::FragType fragCIDHCDETD);

    void moveSamePRMsAllPairs(Spectrum::FragType fragCIDHCD);

    void moveSameSRMsAllPairs(bool checkPRM, Spectrum::FragType fragCIDHCD);

    void movePRMsSRMsAllPairs(Spectrum::FragType fragCIDHCD);

  };
}

#endif /* PAIREDSPECCLUSTER_H_ */
