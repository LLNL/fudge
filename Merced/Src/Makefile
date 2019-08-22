# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

CXX = g++
CXXFLAGS = -Wall -fopenmp

ENSCRIPT_SWITCHES = -2rjC -Ecpp --color --mark-wrapped-lines=arrow --line-numbers -H2
A2PS_SWITCHES = -r --columns=2 -j -C --tabsize=4 --pretty-print=c++ -g -b

#-------------------------------------
# Targets
#-------------------------------------
EXE = merced

OBJS = expn.o igam.o mtherr.o angle_dist.o \
  box_geom.o cm_Legendre.o coef_vector.o coherent.o Compton.o cumulative_points.o \
  data_parser.o dd_vector.o Ecm_Elab_geom.o energy_dist.o \
  energy_dist_base.o energy_function.o Energy_groups.o Eout_integrals.o \
  evaporation.o gamma_2_body_map.o gamma2body.o general_evap.o  \
  global_params.o isotropic.o joint_dist.o kalbach.o kalbach_data.o \
  Legendre_data.o Legendre2Body.o MadlandNix.o mappings.o math_util.o Maxwell.o \
  messaging.o param_base.o phase_space.o reaction.o relativistic.o  \
  standard_Legendre.o transfer.o uncorrelated.o Vcm_Vlab_Hit.o Watt.o \
  x_vector.o main.o

SOURCES = main.cpp          expn.cpp  \
        igam.cpp              mtherr.cpp \
                            protos.h  \
        angle_dist.cpp      angle_dist.hpp    \
        box_geom.cpp        box_geom.hpp    \
        cm_Legendre.cpp     cm_Legendre.hpp  \
        coef_vector.cpp     coef_vector.hpp \
        coherent.cpp        coherent.hpp \
        Compton.cpp         Compton.hpp    \
        cumulative_points.cpp   cumulative_points.hpp \
        data_parser.cpp     data_parser.hpp  \
        dd_vector.cpp       dd_vector.hpp     \
        Ecm_Elab_geom.cpp   Ecm_Elab_geom.hpp  \
        energy_dist.cpp      energy_dist.hpp    \
        energy_dist_base.cpp  energy_dist_base.hpp \
        energy_function.cpp energy_function.hpp \
        Energy_groups.cpp   Energy_groups.hpp \
        Eout_integrals.cpp  Eout_integrals.hpp \
        evaporation.cpp     evaporation.hpp \
        gamma_2_body_map.cpp  gamma_2_body_map.hpp \
        gamma2body.cpp      gamma2body.hpp \
        general_evap.cpp    general_evap.hpp \
        global_params.cpp   global_params.hpp \
        isotropic.cpp       isotropic.hpp   \
        joint_dist.cpp      joint_dist.hpp \
        kalbach.cpp         kalbach.hpp       \
        kalbach_data.cpp    kalbach_data.hpp  \
        Legendre_data.cpp   Legendre_data.hpp          \
        Legendre2Body.cpp   Legendre2Body.hpp  \
        MadlandNix.hpp      MadlandNix.cpp    \
        mappings.cpp        mappings.hpp    \
        math_util.hpp       math_util.hpp  \
        Maxwell.hpp         Maxwell.cpp    \
                            max_output_order.hpp \
        messaging.cpp       messaging.hpp  \
        param_base.cpp      param_base.hpp \
        phase_space.cpp     phase_space.hpp \
                            quadrature.hpp \
        reaction.cpp        reaction.hpp \
        relativistic.cpp    relativistic.hpp \
        standard_Legendre.cpp standard_Legendre.hpp \
        transfer.cpp        transfer.hpp \
        uncorrelated.cpp uncorrelated.hpp \
        Vcm_Vlab_Hit.cpp    Vcm_Vlab_Hit.hpp \
        Watt.hpp            Watt.cpp    \
        x_vector.cpp        x_vector.hpp

#-------------------------------------
# Action item(s)
#-------------------------------------
all: exe 

default : all

exe: $(EXE)

clean:
	@$(RM) -f *.o

cleaner: clean
	@$(RM) -f *.*~ mf*mt* yo*c*i*s* 

realclean: cleaner
	@$(RM) -f $(EXE) $(LIBRARY)

install: all
	cp $(EXE) ../bin

uninstall:

doc:

help:

#-------------------------------------
#  Suffix rules 
#-------------------------------------
.SUFFIXES : .hpp_ps .cpp_ps .cpp .hpp .c .h .o 

.hpp.hpp_ps :
#	 a2ps $(A2PS_SWITCHES) -o $@ $<
	 sed 's/	/  /g' $< | enscript $(ENSCRIPT_SWITCHES) -o $@
.cpp.cpp_ps :
#	 a2ps $(A2PS_SWITCHES) -o $@ $<
	 sed 's/	/  /g' $< | enscript $(ENSCRIPT_SWITCHES) -o $@
.c.o :
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
.cpp.o :
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
.o :
	$(CXX) -o $@ $<

#-------------------------------------
# Rules
#-------------------------------------

$(EXE): $(OBJS)
	$(CXX) $(OBJS) -g -lm  -lstdc++ -fopenmp -o $(EXE)

version.hpp: $(SOURCES)
	python create_version_hpp.py

expn.o: expn.cpp protos.h

igam.o: igam.cpp protos.h

mtherr.o: mtherr.cpp protos.h

angle_dist.o: angle_dist.hpp angle_dist.cpp \
  box_geom.hpp coef_vector.hpp data_parser.hpp dd_vector.hpp \
  Energy_groups.hpp Eout_integrals.hpp \
  global_params.hpp Legendre_data.hpp mappings.hpp math_util.hpp \
  max_output_order.hpp \
  messaging.hpp param_base.hpp quadrature.hpp relativistic.hpp transfer.hpp

box_geom.o: box_geom.hpp box_geom.cpp dd_vector.hpp global_params.hpp \
  messaging.hpp

cm_Legendre.o: cm_Legendre.hpp cm_Legendre.cpp box_geom.hpp coef_vector.hpp \
  cumulative_points.hpp data_parser.hpp dd_vector.hpp Ecm_Elab_geom.hpp \
  Energy_groups.hpp Eout_integrals.hpp \
  global_params.hpp Legendre_data.hpp mappings.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp \
  quadrature.hpp transfer.hpp Vcm_Vlab_Hit.hpp

coef_vector.o: coef_vector.hpp coef_vector.cpp dd_vector.hpp global_params.hpp \
  Legendre_data.hpp max_output_order.hpp messaging.hpp

coherent.o: coherent.hpp coherent.cpp  \
  box_geom.hpp coef_vector.hpp data_parser.hpp dd_vector.hpp \
  Energy_groups.hpp Eout_integrals.hpp \
  global_params.hpp math_util.hpp max_output_order.hpp Legendre_data.hpp \
  max_output_order.hpp messaging.hpp \
  param_base.hpp quadrature.hpp transfer.hpp x_vector.hpp

Compton.o: Compton.hpp Compton.cpp coef_vector.hpp \
  data_parser.hpp dd_vector.hpp Energy_groups.hpp Eout_integrals.hpp \
  global_params.hpp \
  Legendre_data.hpp math_util.hpp max_output_order.hpp messaging.hpp \
  param_base.hpp quadrature.hpp \
  transfer.hpp x_vector.hpp

cumulative_points.o: cumulative_points.hpp cumulative_points.cpp \
  data_parser.hpp dd_vector.hpp messaging.hpp

data_parser.o: data_parser.hpp data_parser.cpp messaging.hpp

dd_vector.o: dd_vector.cpp dd_vector.hpp data_parser.hpp \
  global_params.hpp messaging.hpp

Ecm_Elab_geom.o: Ecm_Elab_geom.hpp Ecm_Elab_geom.cpp box_geom.hpp dd_vector.hpp \
  mappings.hpp math_util.hpp messaging.hpp param_base.hpp Vcm_Vlab_Hit.hpp

energy_dist.o: energy_dist.hpp energy_dist.cpp box_geom.hpp \
  coef_vector.hpp cumulative_points.hpp data_parser.hpp dd_vector.hpp \
  energy_dist_base.hpp \
  Energy_groups.hpp Eout_integrals.hpp global_params.hpp Legendre_data.hpp \
  math_util.hpp max_output_order.hpp \
  messaging.hpp param_base.hpp quadrature.hpp standard_Legendre.hpp transfer.hpp

energy_dist_base.o: energy_dist_base.hpp energy_dist_base.cpp box_geom.hpp \
  coef_vector.hpp cumulative_points.hpp data_parser.hpp dd_vector.hpp \
  Energy_groups.hpp Eout_integrals.hpp \
  global_params.hpp Legendre_data.hpp math_util.hpp max_output_order.hpp \
  messaging.hpp \
  param_base.hpp quadrature.hpp transfer.hpp

energy_function.o: energy_function.hpp energy_function.cpp box_geom.hpp \
  coef_vector.hpp data_parser.hpp dd_vector.hpp Energy_groups.hpp \
  Eout_integrals.hpp \
  global_params.hpp Legendre_data.hpp math_util.hpp max_output_order.hpp \
  messaging.hpp \
  param_base.hpp quadrature.hpp transfer.hpp

Energy_groups.o: Energy_groups.hpp Energy_groups.cpp data_parser.hpp \
  messaging.hpp

Eout_integrals.o: Eout_integrals.cpp  Eout_integrals.hpp coef_vector.hpp \
  data_parser.hpp dd_vector.hpp Legendre_data.hpp \
  math_util.hpp messaging.hpp param_base.hpp

evaporation.o: evaporation.hpp evaporation.cpp protos.h box_geom.hpp \
  coef_vector.hpp data_parser.hpp dd_vector.hpp energy_function.hpp \
  Energy_groups.hpp Eout_integrals.hpp global_params.hpp Legendre_data.hpp \
  math_util.hpp max_output_order.hpp \
  messaging.hpp param_base.hpp quadrature.hpp transfer.hpp

gamma_2_body_map.o: gamma_2_body_map.hpp gamma_2_body_map.hpp dd_vector.hpp \
  global_params.hpp mappings.hpp math_util.hpp messaging.hpp

gamma2body.o: gamma2body.hpp gamma2body.cpp \
  coef_vector.hpp data_parser.hpp dd_vector.hpp Energy_groups.hpp \
  Eout_integrals.hpp gamma_2_body_map.hpp \
  global_params.hpp Legendre_data.hpp mappings.hpp math_util.hpp \
  max_output_order.hpp \
  messaging.hpp param_base.hpp quadrature.hpp transfer.hpp

general_evap.o: general_evap.hpp general_evap.cpp \
  coef_vector.hpp data_parser.hpp dd_vector.hpp Energy_groups.hpp \
  global_params.hpp \
  Legendre_data.hpp math_util.hpp max_output_order.hpp \
  messaging.hpp param_base.hpp quadrature.hpp \
  transfer.hpp

global_params.o: global_params.cpp global_params.hpp data_parser.hpp \
  messaging.hpp

isotropic.o: isotropic.hpp isotropic.cpp box_geom.hpp \
  coef_vector.hpp cumulative_points.hpp data_parser.hpp dd_vector.hpp \
  energy_dist.hpp energy_dist_base.hpp \
  Energy_groups.hpp Eout_integrals.hpp global_params.hpp Legendre_data.hpp \
  math_util.hpp max_output_order.hpp \
  messaging.hpp param_base.hpp quadrature.hpp standard_Legendre.hpp transfer.hpp

joint_dist.o: joint_dist.hpp joint_dist.cpp angle_dist.hpp \
  box_geom.hpp coef_vector.hpp data_parser.hpp dd_vector.hpp Energy_groups.hpp \
  Eout_integrals.hpp global_params.hpp Legendre_data.hpp mappings.hpp \
  math_util.hpp max_output_order.hpp messaging.hpp param_base.hpp \
  quadrature.hpp transfer.hpp

kalbach.o: kalbach.hpp kalbach.cpp box_geom.hpp coef_vector.hpp data_parser.hpp \
  dd_vector.hpp Ecm_Elab_geom.hpp energy_dist_base.hpp Energy_groups.hpp \
  Eout_integrals.hpp \
  global_params.hpp kalbach_data.hpp Legendre_data.hpp mappings.hpp \
  math_util.hpp max_output_order.hpp messaging.hpp param_base.hpp \
  quadrature.hpp transfer.hpp \
  Vcm_Vlab_Hit.hpp

kalbach_data.o: kalbach_data.hpp kalbach_data.cpp data_parser.hpp \
  dd_vector.hpp mappings.hpp messaging.hpp

Legendre_data.o: Legendre_data.hpp Legendre_data.cpp data_parser.hpp \
  dd_vector.hpp messaging.hpp

Legendre2Body.o: Legendre2Body.cpp Legendre2Body.hpp angle_dist.hpp  \
  box_geom.hpp coef_vector.hpp data_parser.hpp dd_vector.hpp \
  Energy_groups.hpp Eout_integrals.hpp \
  global_params.hpp Legendre_data.hpp mappings.hpp math_util.hpp \
  max_output_order.hpp \
  messaging.hpp param_base.hpp quadrature.hpp relativistic.hpp transfer.hpp

# logger.o: logger.cpp logger.hpp data_parser.hpp global_params.hpp version.hpp

mappings.o: mappings.cpp mappings.hpp dd_vector.hpp global_params.hpp \
  math_util.hpp messaging.hpp

math_util.o: math_util.cpp math_util.hpp coef_vector.hpp data_parser.hpp \
  dd_vector.hpp  global_params.hpp Legendre_data.hpp max_output_order.hpp \
  messaging.hpp param_base.hpp \
  quadrature.hpp

MadlandNix.o: MadlandNix.hpp MadlandNix.cpp protos.h box_geom.hpp coef_vector.hpp \
  data_parser.hpp dd_vector.hpp energy_function.hpp Energy_groups.hpp \
  Eout_integrals.hpp \
  evaporation.hpp global_params.hpp Legendre_data.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp \
  param_base.hpp quadrature.hpp transfer.hpp

Maxwell.o: Maxwell.hpp Maxwell.cpp protos.h box_geom.hpp coef_vector.hpp \
  data_parser.hpp dd_vector.hpp energy_function.hpp Energy_groups.hpp \
  Eout_integrals.hpp \
  evaporation.hpp global_params.hpp Legendre_data.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp \
  param_base.hpp quadrature.hpp transfer.hpp

messaging.o: messaging.cpp messaging.hpp global_params.hpp

param_base.o: param_base.cpp param_base.hpp coef_vector.hpp data_parser.hpp \
  dd_vector.hpp Energy_groups.hpp Eout_integrals.hpp Legendre_data.hpp \
  max_output_order.hpp messaging.hpp

phase_space.o: phase_space.hpp phase_space.cpp box_geom.hpp coef_vector.hpp data_parser.hpp \
  dd_vector.hpp Ecm_Elab_geom.hpp energy_dist_base.hpp Energy_groups.hpp \
  Eout_integrals.hpp Legendre_data.hpp mappings.hpp \
  math_util.hpp max_output_order.hpp messaging.hpp param_base.hpp \
  quadrature.hpp transfer.hpp Vcm_Vlab_Hit.hpp

reaction.o: reaction.cpp reaction.hpp angle_dist.hpp box_geom.hpp cm_Legendre.hpp \
  coef_vector.hpp coherent.hpp Compton.hpp cumulative_points.hpp data_parser.hpp \
  dd_vector.hpp Ecm_Elab_geom.hpp energy_dist.hpp energy_dist_base.hpp \
  energy_function.hpp Energy_groups.hpp Eout_integrals.hpp evaporation.hpp \
  gamma_2_body_map.hpp gamma2body.hpp \
  general_evap.hpp global_params.hpp isotropic.hpp joint_dist.hpp kalbach.hpp \
  kalbach_data.hpp Legendre_data.hpp Legendre2Body.hpp mappings.hpp \
  math_util.hpp MadlandNix.hpp max_output_order.hpp \
  Maxwell.hpp messaging.hpp param_base.hpp phase_space.hpp quadrature.hpp \
  relativistic.hpp standard_Legendre.hpp transfer.hpp uncorrelated.hpp \
  Vcm_Vlab_Hit.hpp Watt.hpp x_vector.hpp

relativistic.o: relativistic.cpp relativistic.hpp dd_vector.hpp global_params.hpp \
  mappings.hpp math_util.hpp messaging.hpp

standard_Legendre.o: standard_Legendre.hpp standard_Legendre.cpp box_geom.hpp \
  coef_vector.hpp cumulative_points.hpp data_parser.hpp dd_vector.hpp Energy_groups.hpp \
  Eout_integrals.hpp global_params.hpp Legendre_data.hpp math_util.hpp \
  max_output_order.hpp messaging.hpp param_base.hpp quadrature.hpp transfer.hpp

transfer.o: transfer.cpp transfer.hpp coef_vector.hpp data_parser.hpp \
  dd_vector.hpp Energy_groups.hpp global_params.hpp Legendre_data.hpp \
  max_output_order.hpp messaging.hpp quadrature.hpp 

uncorrelated.o: uncorrelated.hpp uncorrelated.cpp angle_dist.hpp \
  box_geom.hpp coef_vector.hpp cumulative_points.hpp \
  data_parser.hpp dd_vector.hpp energy_dist_base.hpp \
  energy_dist.hpp Energy_groups.hpp global_params.hpp Legendre_data.hpp \
  Legendre2Body.hpp mappings.hpp \
  math_util.hpp max_output_order.hpp messaging.hpp param_base.hpp \
  quadrature.hpp transfer.hpp

Vcm_Vlab_Hit.o: Vcm_Vlab_Hit.hpp Vcm_Vlab_Hit.cpp box_geom.hpp \
  mappings.hpp math_util.hpp messaging.hpp

Watt.o: Watt.hpp Watt.cpp protos.h box_geom.hpp coef_vector.hpp data_parser.hpp \
  dd_vector.hpp energy_function.hpp Energy_groups.hpp Eout_integrals.hpp \
  evaporation.hpp \
  global_params.hpp Legendre_data.hpp math_util.hpp max_output_order.hpp \
  messaging.hpp \
  param_base.hpp quadrature.hpp transfer.hpp

x_vector.o: x_vector.hpp x_vector.cpp data_parser.hpp dd_vector.hpp \
  global_params.hpp messaging.hpp

main.o: main.cpp angle_dist.hpp box_geom.hpp cm_Legendre.hpp coef_vector.hpp \
  coherent.hpp Compton.hpp cumulative_points.hpp data_parser.hpp dd_vector.hpp \
  Ecm_Elab_geom.hpp energy_dist.hpp energy_dist_base.hpp evaporation.hpp \
  energy_function.hpp Energy_groups.hpp Eout_integrals.hpp gamma_2_body_map.hpp  \
  gamma2body.hpp general_evap.hpp \
  global_params.hpp isotropic.hpp joint_dist.hpp kalbach.hpp kalbach_data.hpp \
  Legendre_data.hpp Legendre2Body.hpp mappings.hpp math_util.hpp \
  MadlandNix.hpp max_output_order.hpp \
  Maxwell.hpp messaging.hpp reaction.hpp param_base.hpp phase_space.hpp \
  quadrature.hpp relativistic.hpp standard_Legendre.hpp transfer.hpp \
  uncorrelated.hpp Vcm_Vlab_Hit.hpp Watt.hpp x_vector.hpp
