# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

CXX = g++
CXXFLAGS = -c -g -Wall

EXE = fig_05

SOURCES = fig_05.cpp

OBJS = fig_05.o

$(EXE): $(OBJS)
	$(CXX) $(OBJS) -g -lm  -lstdc++ -o $(EXE)	

fig_05.o: fig_05.cpp
	$(CXX) $(CXXFLAGS) -o fig_05.o fig_05.cpp

