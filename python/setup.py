# setup.py
#from distutils.core import setup, Extension
from setuptools import setup, Extension
import distutils
from os.path import join
import numpy
import os
import subprocess
import sysconfig
import sys
from distutils.sysconfig import get_python_inc


from os.path import expanduser
home = expanduser("~")

pecos_root = join(home,'software','pecos')
pecos_build_dir = join(home,'software','pecos','build')

#boost_include = join(home,'local/boost-1.61.0/include')
boost_include = join(home,'modules-gnu-4.8.5/boost-1.54.0/include')

package_name="PyDakota"

#swig_src_include = os.getcwd()#join(pecos_root,'surrogates/python/src')
distutils.log.set_verbosity(1)
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include  = numpy.get_numpy_include()

# Generate swigpyrun.h which defines WIG_TypeQuery, SWIG_NewPointerObj, etc.
subprocess.call(['swig','-python','-external-runtime','swigpyrun.h'])

include_dirs = [
    numpy_include,
    boost_include,
    os.getcwd(),      #include swigpyrun.h generated above
#    swig_src_include, #include numpy_include.hpp
    get_python_inc(), #include Python.h
]

pecos_include_dirs=[
    join(pecos_root,'util','src'),
    join(pecos_root,'surrogates','python','src')
]

def distutils_dir_name(dname):
    """Returns the name of a distutils build directory"""
    f = "{dirname}.{platform}-{version[0]}.{version[1]}"
    return f.format(dirname=dname,
                    platform=sysconfig.get_platform(),
                    version=sys.version_info)

# The location disutils will put the compiled extension files. 
# I am not forcing this. If default location of library is not used base
# directory for build library or changed by disutils in the future then 
# this will be incorrect
distutils_build_dir = join(os.getcwd(),'build', distutils_dir_name('lib'))
if not os.path.exists(join(distutils_build_dir,package_name)):
    os.makedirs(join(distutils_build_dir,package_name))

swig_opts = ['-c++']
# distutils puts module.py of a swig extension in the current folder 
#If we want module.py files producted by swig to be in a package 
# (so can be called using package.module) use
swig_opts += ['-outdir','%s'%join(distutils_build_dir,package_name)]
#if want module.py in subpackage use
#swig_opts = ['-c++', '-outdir %s'%(join(disutils_build_dir,package_name,subspackage_name)]

options_list_srcs = ['python_helpers.cpp',"test_options_list.cpp"]

options_list = Extension(
    '_options_list',
    ['OptionsList.i']+options_list_srcs,
    include_dirs = include_dirs+[join(pecos_root,'util','src')],
    define_macros =[('COMPILE_WITH_PYTHON',None)],
    undef_macros = [],
    language='c++',
    library_dirs = [],
    libraries = [],
    extra_compile_args = ['-std=c++11'],
    swig_opts=swig_opts+['-I%s'%include_dir for include_dir in pecos_include_dirs]
)

teuchos_root = join(pecos_root,'packages','teuchos','packages','teuchos')
# Include following dir to include Teuchos_DLLExportMacro.h
teuchos_build_dir = join(
    pecos_build_dir,'packages','teuchos','packages','teuchos','src')
teuchos_include_dirs = [join(teuchos_root,'src'),teuchos_build_dir]
pecos_include_dirs=[
    join(pecos_root,'util','src'),
    join(pecos_root,'surrogates','python','src')
]
# link to teuchos library dir need to get typeinfo for Teuchos::CompObject
library_dirs=[
        join(pecos_build_dir,'util','src'),
        join(pecos_build_dir,'packages','teuchos','packages','teuchos','src'),
        join(pecos_build_dir,'surrogates','models','src')]

math_tools = Extension(
    '_math_tools',
    ['math_tools.i']+options_list_srcs,
    include_dirs = include_dirs+pecos_include_dirs+teuchos_include_dirs,
    define_macros =[('COMPILE_WITH_PYTHON',None)],
    undef_macros = [],
    language='c++',
    library_dirs = library_dirs,
    libraries = ['models','pecos_util','teuchos','blas','lapack'],
    extra_compile_args = ['-std=c++11','-Wno-unused-local-typedefs'],
    swig_opts=swig_opts+['-I%s'%include_dir for include_dir in pecos_include_dirs]+
['-I%s'%include_dir for include_dir in teuchos_include_dirs])

regression = Extension(
    '_regression',
    ['regression.i']+options_list_srcs,
    include_dirs = include_dirs+pecos_include_dirs+teuchos_include_dirs,
    define_macros =[('COMPILE_WITH_PYTHON',None)],
    undef_macros = [],
    language='c++',
    library_dirs = library_dirs,
    libraries = ['models','pecos_util','teuchos','blas','lapack'],
    extra_compile_args = ['-std=c++11','-Wno-unused-local-typedefs'],
    swig_opts=swig_opts+['-I%s'%include_dir for include_dir in pecos_include_dirs]+
['-I%s'%include_dir for include_dir in teuchos_include_dirs])

pecos_include_dirs=[
    join(pecos_root,'util','src'),
    join(pecos_root,'surrogates','models','src'),
    join(pecos_root,'surrogates','python','src')
]

approximation = Extension(
    '_approximation',
    ['approximations.i']+options_list_srcs,
    include_dirs = include_dirs+pecos_include_dirs+teuchos_include_dirs,
    define_macros =[('COMPILE_WITH_PYTHON',None)],
    undef_macros = [],
    language='c++',
    library_dirs = library_dirs,
    libraries = ['models','pecos_util','teuchos','blas','lapack'],
    extra_compile_args = ['-std=c++11','-Wno-unused-local-typedefs'],
    swig_opts=swig_opts+['-I%s'%include_dir for include_dir in pecos_include_dirs]+
['-I%s'%include_dir for include_dir in teuchos_include_dirs])

pecos_include_dirs=[join(pecos_root,'src')]
library_dirs=[
    join(pecos_build_dir,'src'),
    join(pecos_build_dir,'packages','teuchos','packages','teuchos','src'),
    join(pecos_build_dir,'packages','VPISparseGrid','src')]

univariate_polynomials = Extension(
    '_univariate_polynomials',
    ['univariate_polynomials.i'],
    include_dirs = include_dirs+pecos_include_dirs+teuchos_include_dirs,
    define_macros =[('COMPILE_WITH_PYTHON',None)],
    undef_macros = [],
    language='c++',
    library_dirs = library_dirs,
    libraries = ['pecos_src','sparsegrid','teuchos','lapack','blas'],
    extra_compile_args = ['-std=c++11','-Wno-unused-local-typedefs'],
    swig_opts=swig_opts+['-I%s'%include_dir for include_dir in pecos_include_dirs]+
['-I%s'%include_dir for include_dir in teuchos_include_dirs])

import unittest
def PyDakota_test_suite():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover(distutils_build_dir, pattern='test_*.py')
    return test_suite

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = package_name,
    version = "1.0",
    description="Python wrapper of Dakota surrogate and linear algebra utilities",
    author="John Jakeman",
    author_email='dakota-users@sandia.gov',
    license='GNU Lesser General Public License (LGPL)',
    url='https://dakota.sandia.gov',
    keywords='Approximation, Surrogates, Uncertainty Quantification, Polynomial Chaos',
    long_description=read('README'),
    packages=[package_name,package_name+'.unit'],
    package_dir = {'': 'src'},
    ext_package=package_name,
    ext_modules=[
        options_list,
        math_tools,
        regression,
        approximation,
        univariate_polynomials],
    package_data={package_name:[join('unit','data/*.gz')]},
    test_suite='setup.PyDakota_test_suite')

print join(pecos_root,'surrogates','models','unit')

print('Run python setup.py install to install to python site-packages')
print('Run python setup.py install --prefix=<dir>to install to python dir>')
string = 'export PYTHONPATH=$PYTHONPATH:%s'%distutils_build_dir
print('If you choose not to install, then you can temporarily add the build python modules to the PYTHONPATH by running in the terminal')
print (string)

# dont forget to set pythonpath
#notes:
# warning: command line option '-Wstrict-prototypes' is valid for C/ObjC but not for C++ [enabled by default] is a bug with distutils

#To create a module in a python package use
# Extension(package_name+'.'+module_name)
#To create a module in a subpackage of a python packages use
# Extension(package_name+'.'+subpackage_name+'.'+module_name)

# If get an error cannot import _modulename.so then this may be due to a
# linking error. E.g. forgot to link in a required library using the libraries
# keywork in Extension
