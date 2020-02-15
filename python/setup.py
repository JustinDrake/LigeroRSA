from distutils.core import setup
from distutils.extension import Extension
import os

os.environ["CC"] = "clang++ -std=c++14"

ligero = Extension(
    'ligero',
    include_dirs=['/Users/ykashnikov/ligero/RSACeremony/include'
        ,'/usr/local/Cellar/eigen/3.3.7/include/eigen3/'
        ,'/Users/ykashnikov/ligero/RSACeremony/depends/exrandom/include'
        ,'/Users/ykashnikov/ligero/RSACeremony/depends/easyloggingpp/src'
        ],
    sources=['ligero.cpp', '../depends/easyloggingpp/src/easylogging++.cc'],
    libraries=['boost_python37-mt', 'gmp']
)

setup(
    name='ligero',
    version='0.1',
    ext_modules=[ligero])
