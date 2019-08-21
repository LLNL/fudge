# <<BEGIN-copyright>>
# <<END-copyright>>

"""
This module contains the distribution class and all its sub-classes.

A product outgoing distribution (i.e., outgoing angle and energy) can be represented in different ways called 
'distribution genres'.  As example, for a two-body outgoing channel only the outgoing product's angular distribution 
in the center-of-mass frame is needed, as its energy can be calculated from the angle. In this case, one genre for 
the distribution would be P(mu|E) where E is the incident energy, mu = cos( angle ), angle is the product's outgoing  
angle in the center-of-mass frame and P(mu|E) is the probability for the product to be emitted at mu when the incident 
product's energy is E. In gnd this genre is label 'angularTwoBody' (since the incident energy E is always 
present, it is ignored in the naming convention, and only the outgoing variables are used for naming).  Deterministic 
transport codes cannot use the 'angularTwoBody' directly; instead, they want the data converted to
grouped Legendre data. Storing the data in Legendre form - whether grouped, pointwise, etc. - is call 'Legendre'

For a given distribution genre, there can be several representations of the data called forms. As example, 'angularTwoBody'
data can be stored as 'pointwise' or 'equalProbableBins' data. 

Finally, some distribution genres require several 'components' to represent the data. As example, the 'angularEnergy' 
requires an angular component, P_a(mu|E), and an angularEnergy component, P_ae(E'|E,mu), where the distribution is the product 
P_a(E|mu) * P_ae(E'|E,mu).

The class distribution is the base class that should be used for each genre. That is, each genre shall be its own class which
inherits the distribution class.

The following table list the supported distribution genres as well as the components and forms.

Components                  Forms
----------------------------------------------------------
angular.component
                            pointwise
                            equalProbableBins
                            isotropic

Legendre.component
                            LLNLPointwise
                            grouped

Legendre.energyConservationComponent
                            energyConservationGrouped

CoulombElasticComponent
                            pointwise
                            grouped
                            nuclearPlusCoulombInterference

energy.component
                            equalProbableBins

angularEnergy.LLNLComponent
                            endlI3 - pointwise
                            LLNAngularEnergyEqualProbableBinsForm

LLNL_withAngularComponent
uncorrelated.component
referenceComponent
"""

import base
import angular
import energy
import energyAngular
import angularEnergy
import uncorrelated
import Legendre
import miscellaneous

from base import parseXMLNode
