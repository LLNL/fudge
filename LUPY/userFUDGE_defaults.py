# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

try:
    import myFUDGE_defaults
except:
    myFUDGE_defaults = None

def add_argument(parser, name, required, *args, **kwargs):
    '''
    Adds the option *args* to *parser* by calling *parser*'s add_argument method. In addition, 
    adds a default value and more documentation to the option if a variable named *name* 
    exists in a user defined "myFUDGE_defaults.py" file. The file "myFUDGE_defaults.py" 
    must be in the user's PYTHONPATH.

    :param parser:      argparse.parse_args instance.
    :param name:        Name of the "myFYDGE_default.py" variable to use as a default.
    :param required:    If **True** and *name* is not defined in "myFUDGE_defaults.py", the add_argument **required** keyword is set to **True**.
    :param args:        Option name and its aliases.
    :param kwargs:      Keyword arguments passed to parser.add_argument.
    '''

    defaultValue = None
    if myFUDGE_defaults is not None:
        defaultValue = getattr(myFUDGE_defaults, name, None)

    if defaultValue is not None:
        kwargs['default'] = defaultValue
        if 'help' in kwargs:
            kwargs['help'] += ' Default is "%s".' % defaultValue

    if required and defaultValue is None:
        kwargs['required'] = True

    parser.add_argument(*args, **kwargs)

    if not hasattr(parser, 'FUDGE_defaultText'):
        parser.description += '\n\nUser option | default variable name:\n'
        parser.FUDGE_defaultText = True
    parser.description += '    %s | %s' % (args[-1], name)
