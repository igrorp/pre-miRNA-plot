
from setuptools import setup, find_packages
from help import desctxt

setup(name='premiRNAplot',
      version=0.2,
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
          'numpy'
      ],
      long_description=desctxt
)
