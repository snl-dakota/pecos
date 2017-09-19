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

include_pecos=True

from os.path import expanduser
home = expanduser("~")

pecos_root = os.path.join(os.path.split(os.path.abspath(__file__))[0],'..')
#pecos_root = os.path.join(home,'software','pecos')

# asssumes pecos is built in a subdirectory pecos/build
pecos_build_dir = join(pecos_root,'build')

# Get location of boost include from Pecos CMakeCache.txt generated when Pecos
# libraries was built
cmakecache_filename = os.path.join(pecos_build_dir,'CMakeCache.txt')
with open(cmakecache_filename, 'r') as f:
    file_string = f.read()
    index1=file_string.find('Boost_INCLUDE_DIR:PATH=')
    index2=file_string.find('=',index1)+1
    index3=file_string.find('\n',index2)
    boost_include=file_string[index2:index3]

package_name="PyDakota"

distutils.log.set_verbosity(1)
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include  = numpy.get_numpy_include()

# Generate swigpyrun.h which defines WIG_TypeQuery, SWIG_NewPointerObj, etc.
subprocess.call(['swig','-python','-external-runtime','swigpyrun.h'])

base_include_dirs = [
    numpy_include,
    boost_include,
    join(os.getcwd(),'cpp_src'),
    os.getcwd(),   #include swigpyrun.h generated above
    get_python_inc(), #include Python.h
]

surrogates_include_dirs=[
    join(pecos_root,'util','src'),
    join(pecos_root,'surrogates','models','src'),
    ]

pecos_include_dirs=[
    join(pecos_root,'src')
]

teuchos_root = join(pecos_root,'packages','teuchos','packages','teuchos')
# Include following dir to include Teuchos_DLLExportMacro.h
teuchos_build_dir = join(
    pecos_build_dir,'packages','teuchos','packages','teuchos','src')
teuchos_include_dirs = [
    join(teuchos_root,'src'),
    teuchos_build_dir]

pydakota_srcs = [join('cpp_src','python_helpers.cpp')]

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

base_swig_opts = ['-c++']
# distutils puts module.py of a swig extension in the current folder
# If we want module.py files producted by swig to be in a package
# (so can be called using package.module) use
package_swig_opts = base_swig_opts+[
    '-outdir',
    '%s'%join(distutils_build_dir,package_name)]

options_list_include_dirs = base_include_dirs+[join(pecos_root,'util','src')]+\
  surrogates_include_dirs
options_list = Extension(
    '_options_list',
    [join('cpp_src','OptionsList.i')]+pydakota_srcs,
    include_dirs = options_list_include_dirs,
    define_macros =[('COMPILE_WITH_PYTHON',None)],
    undef_macros = [],
    language='c++',
    library_dirs = [],
    libraries = [],
    extra_compile_args = ['-std=c++11'],
    swig_opts=package_swig_opts+['-I%s'%include_dir for include_dir in options_list_include_dirs]
)

# link to teuchos library dir need to get typeinfo for Teuchos::CompObject
math_tools_library_dirs=[
    join(pecos_build_dir,'util','src'),
    join(pecos_build_dir,'packages','teuchos','packages','teuchos','src'),
    join(pecos_build_dir,'surrogates','models','src')]

math_tools_include_dirs=base_include_dirs+surrogates_include_dirs+\
  teuchos_include_dirs
math_tools = Extension(
    '_math_tools',
    [join('cpp_src','math_tools.i')]+pydakota_srcs,
    include_dirs = math_tools_include_dirs,
    define_macros =[('COMPILE_WITH_PYTHON',None)],
    undef_macros = [],
    language='c++',
    library_dirs = math_tools_library_dirs,
    libraries = ['models','pecos_util','teuchos','blas','lapack'],
    extra_compile_args = ['-std=c++11','-Wno-unused-local-typedefs'],
    swig_opts=package_swig_opts+['-I%s'%include_dir for include_dir in math_tools_include_dirs])

regression_include_dirs = base_include_dirs+surrogates_include_dirs+\
  teuchos_include_dirs
regression_library_dirs = [
    join(pecos_build_dir,'util','src'),
    join(pecos_build_dir,'packages','teuchos','packages','teuchos','src'),
    join(pecos_build_dir,'surrogates','models','src')]
regression = Extension(
    '_regression',
    [join('cpp_src','regression.i')]+pydakota_srcs,
    include_dirs = regression_include_dirs,
    define_macros =[('COMPILE_WITH_PYTHON',None)],
    undef_macros = [],
    language='c++',
    library_dirs = regression_library_dirs,
    libraries = ['models','pecos_util','teuchos','blas','lapack'],
    extra_compile_args = ['-std=c++11','-Wno-unused-local-typedefs'],
    swig_opts=package_swig_opts+['-I%s'%include_dir for include_dir in regression_include_dirs])

approximation_include_dirs=base_include_dirs+surrogates_include_dirs+\
  teuchos_include_dirs
approximation_library_dirs = regression_library_dirs
approximation_libraries = ['models','pecos_util','teuchos','blas','lapack']
if include_pecos:
    approximation_include_dirs+=[join(pecos_root,'surrogates','pecos_wrapper','src')]+pecos_include_dirs+[join(pecos_root,'packages','teuchos','packages','teuchos','src'),join(pecos_root,'packages','VPISparseGrid','src')]
    approximation_library_dirs += [
        join(pecos_build_dir,'src'),
        join(pecos_build_dir,'packages','VPISparseGrid','src'),
        join(pecos_build_dir,'surrogates','pecos_wrapper','src')]
approximation_libraries = ['pecos_wrapper','pecos_src','sparsegrid']+approximation_libraries#order of libraries is important
approximation_srcs=[join('cpp_src','approximations.i')]+pydakota_srcs
print approximation_srcs
approximation = Extension(
    '_approximation',
    approximation_srcs,
    include_dirs = approximation_include_dirs,
    define_macros =[('COMPILE_WITH_PYTHON',None)],
    undef_macros = [],
    language='c++',
    library_dirs = approximation_library_dirs,
    libraries = approximation_libraries,
    extra_compile_args = ['-std=c++11','-Wno-unused-local-typedefs'],
    swig_opts=package_swig_opts+['-I%s'%include_dir for include_dir in approximation_include_dirs])

univariate_polynomials_include_dirs=base_include_dirs+pecos_include_dirs+teuchos_include_dirs+surrogates_include_dirs
univariate_polynomials_library_dirs=[
    join(pecos_build_dir,'src'),
    join(pecos_build_dir,'packages','teuchos','packages','teuchos','src'),
    join(pecos_build_dir,'packages','VPISparseGrid','src')]
univariate_polynomials = Extension(
    '_univariate_polynomials',
    [join('cpp_src','univariate_polynomials.i')]+pydakota_srcs,
    include_dirs = univariate_polynomials_include_dirs,
    define_macros =[('COMPILE_WITH_PYTHON',None)],
    undef_macros = [],
    language='c++',
    library_dirs = univariate_polynomials_library_dirs,
    libraries = ['pecos_src','sparsegrid','teuchos','lapack','blas'],
    extra_compile_args = ['-std=c++11','-Wno-unused-local-typedefs'],
    swig_opts=package_swig_opts+['-I%s'%include_dir for include_dir in univariate_polynomials_include_dirs])

subpackage_name = 'swig_examples'
# The location disutils will put the compiled extension files. 
# I am not forcing this. If default location of library is not used base
# directory for build library or changed by disutils in the future then 
# this will be incorrect
if not os.path.exists(join(distutils_build_dir,package_name, subpackage_name)):
    os.makedirs(join(distutils_build_dir,package_name,subpackage_name))

std_vector_example_include_dirs=base_include_dirs+['cpp_src/swig_examples_src']
swig_opts = ['-c++','-outdir','%s'%join(
    distutils_build_dir,package_name,subpackage_name)]
std_vector_example = Extension(
    '%s._std_vector_example'%subpackage_name,
    [join('cpp_src','swig_examples_src','std_vector_example.i')],
    include_dirs = std_vector_example_include_dirs,
    define_macros =[('COMPILE_WITH_PYTHON',None)],
    undef_macros = [],
    language='c++',
    library_dirs = [],
    libraries = [],
    extra_compile_args = ['-std=c++11'],
    swig_opts=swig_opts+['-I%s'%include_dir for include_dir in std_vector_example_include_dirs])

options_list_interface_include_dirs=['cpp_src/swig_examples_src']+options_list_include_dirs
swig_opts = ['-c++','-outdir','%s'%join(
    distutils_build_dir,package_name,subpackage_name)]
options_list_interface = Extension(
    '%s._options_list_interface'%subpackage_name,
    [join('cpp_src','swig_examples_src','options_list_interface.i'),
     join('cpp_src','swig_examples_src','options_list_interface.cpp')]+pydakota_srcs,
    include_dirs = options_list_interface_include_dirs,
    define_macros =[('COMPILE_WITH_PYTHON',None)],
    undef_macros = [],
    language='c++',
    library_dirs = [],
    libraries = [],
    extra_compile_args = ['-std=c++11'],
    swig_opts=swig_opts+['-I%s'%include_dir for include_dir in options_list_interface_include_dirs])

enum_example_include_dirs=['cpp_src/swig_examples_src']
swig_opts = ['-c++','-outdir','%s'%join(
    distutils_build_dir,package_name,subpackage_name)]
enum_example = Extension(
    '%s._enum_example'%subpackage_name,
    [join('cpp_src','swig_examples_src','enum_example.i')],
    include_dirs = enum_example_include_dirs,
    define_macros =[('COMPILE_WITH_PYTHON',None)],
    undef_macros = [],
    language='c++',
    library_dirs = [],
    libraries = [],
    extra_compile_args = ['-std=c++11'],
    swig_opts=swig_opts+['-I%s'%include_dir for include_dir in enum_example_include_dirs])

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
    packages=[package_name,package_name+'.unit',package_name+'.swig_examples'],
    package_dir = {'': 'python_src'},
    ext_package=package_name,
    ext_modules=[
        options_list,
        math_tools,
        regression,
        approximation,
        std_vector_example,
        options_list_interface,
        enum_example,
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

# on osx use
# conda install numpy scipy matplotlib cmake boost=1.61.0 openblas=0.2.19 gcc swig
# need to do export DYLD_FALLBACK_LIBRARY_PATH=~/miniconda2/envs/fenics/lib # this gets gfotran library on  library path
