#!/bin/bash
#
# This is a simplified operational example. It can produce the data if unlimited computational resources are available, 
# but it is far from being realistic on a single desktop or computational node. For the paper we used a relatively
# small computational cluster with checkpoints (see the scripts in cluster/). The final figure was produced using
# OriginPro software.
#
# See also the comments inside the calc.sh. In particular, running this script out of the box, requires providing
# either scattnlay of bhfield for reference simulations.

# Cycle over RI corresponds to parts a,b,c of Fig.16
for m in 1.01 1.05 1.1; do 
  sh calc.sh 500 $m 800
done
