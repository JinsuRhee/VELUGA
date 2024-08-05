import sys
from setuptools import setup, find_packages

try:
    from numpy.distutils.core import Extension, setup
except: ImportError:
    print("Numpy is not installed. Installing...")
    import subprocess
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'numpy'])
    from numpy.distutils.core import Extension, setup



# Fortran extension
fortran_extension = Extension(
	name=[	'src.fortran.get_ptcl_py', 
			'src.get_cell_py',
			'src.fortran.get_amr_py',
			'src.fortran.find_domain_py'], 
	sources=[	'src/fortran/get_ptcl_py.f90',
				'src/fortran/get_cell_py.f90',
				'src/fortran/get_amr_py.f90',
				'src/fortran/find_domain_py.f90'], )

setup(
	name='VELUGA', 
	version='0.01',
	packages=find_packages(),
	ext_modules=[fortran_extension],
	install_requires=['numpy', 'scipy', 'pickle5', 'h5py'],
	)
