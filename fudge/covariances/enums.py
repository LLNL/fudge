# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import enums as enumsModule


class Type(enumsModule.Enum):
    relative = enumsModule.auto()
    absolute = enumsModule.auto()
