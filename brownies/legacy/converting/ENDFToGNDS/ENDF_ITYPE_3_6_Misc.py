# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Defines dictionaries for handling ENDF ITYPE=3 and 6 data  (photo-atomic, electro-atomic and atomic relaxation)
"""

MT_AtomicConfigurations = {
    534 : '1s1/2',      535 : '2s1/2',      536 : '2p1/2',      537 : '2p3/2',      538 : '3s1/2',
    539 : '3p1/2',      540 : '3p3/2',      541 : '3d3/2',      542 : '3d5/2',      543 : '4s1/2',
    544 : '4p1/2',      545 : '4p3/2',      546 : '4d3/2',      547 : '4d5/2',      548 : '4f5/2',
    549 : '4f7/2',      550 : '5s1/2',      551 : '5p1/2',      552 : '5p3/2',      553 : '5d3/2',
    554 : '5d5/2',      555 : '5f5/2',      556 : '5f7/2',      557 : '5g7/2',      558 : '5g9/2',
    559 : '6s1/2',      560 : '6p1/2',      561 : '6p3/2',      562 : '6d3/2',      563 : '6d5/2',
    564 : '6f5/2',      565 : '6f7/2',      566 : '6g7/2',      567 : '6g9/2',      568 : '6h9/2',
    569 : '6h11/2',     570 : '7s1/2',      571 : '7p1/2',      572 : '7p3/2' }

# also create reverse lookup for translating back to ENDF-6:
AtomicConfigurations_MT = {}
for key,value in MT_AtomicConfigurations.items():
    AtomicConfigurations_MT[value] = key
