# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

print("WARNING: namespace 'fudge.gnds' is deprecated! Please import directly from fudge instead, e.g. 'from fudge import reactionSuite'.")

from fudge import outputChannel, documentation, externalFile, product, reactionSuite, styles, suites, sums
from fudge import channelData, covariances, productData, reactionData, reactions, resonances
