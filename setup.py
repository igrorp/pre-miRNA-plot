
from setuptools import setup, find_packages

setup(
	name='premirnaplot',
	version=1.16,
	description='pre-miRNA secondary structure prediction image generator',
	author='Igor Paim',
	author_email='igorpaim8@gmail.com',
	url='https://github.com/igrorp/pre-miRNA-plot',
	license='GPL',
	classifiers=[
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'Programming Language :: Python :: 3',
	],
	keywords='pre-miRNA',
	packages=find_packages(),
	install_requires=[
		'matplotlib',
		'scikit-learn',
		'numpy',
		'svglib',
		'svgwrite',
		'pandas',
	],
	scripts=[
		'premirnaplot/premirnaplot',
	],
	long_description=''
)
