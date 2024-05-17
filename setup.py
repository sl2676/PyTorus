import os, sys

from distutils.core import setup, Extension
from distutils import sysconfig

#cpp_args = ['-std=c++17', '-stdlib=libc++', '-Wno-register','-mmacosx-version-min=10.7']
cpp_args = ['-std=c++17', '-stdlib=libc++', '-Wno-register', '-mmacosx-version-min=10.14', '-Wall']

ext_modules = [
	Extension(
	'PyTorus', 
	['wrap.cpp', './Torus/src/utils/WDMath.cc', './Torus/src/utils/Compress.cc'],
	include_dirs=['pybind11/include'],
	language='c++',
	extra_compile_args=cpp_args,
	),
]

setup(
	name='PyTorus',
	version='0.0.1',
	author='Sean Ly',
	author_email='seanly1101@gmail.com',
	description='torus library used for whatever',
	ext_modules=ext_modules,
)
