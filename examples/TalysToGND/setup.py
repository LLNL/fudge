import glob, os
from distutils.core import setup

setup( 
    name = 'geft',\
    version = '1.0',\
    author = 'Neil Summers',\
    author_email = 'summers21@llnl.gov',\
    packages = [ 'geft', 'geft.test', 'geft.test.driver_tests' ],\
    package_data = { \
        'data': [ os.sep.join( [ 'dicts', '*.txt' ] ), os.sep.join( [ 'data', '*.t*' ] ), os.sep.join( [ 'data', '*.pickle' ] ), os.sep.join( [ 'data', 'X4all', '*', '*.x4' ] ) ], \
        'x4i.test.driver_tests': [ '*.in', '*.out' ] \
    }, \
    url = 'http://nuclear.llnl.gov/geft',\
    license = open( 'LICENSE.txt' ).read(),\
    description = 'A tool for Generating Evaluations From Talys.',\
    long_description = open( 'README.txt' ).read()\
)
