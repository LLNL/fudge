# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

try:
    from crossSectionAdjustForHeatedTarget import crossSectionAdjustForHeatedTarget as heat

except (ImportError, ModuleNotFoundError):
    from crossSectionAdjustForHeatedTarget.lib import crossSectionAdjustForHeatedTarget as heat