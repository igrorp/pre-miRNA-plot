
from setuptools import setup, find_packages
from help import desctxt

setup(name='premiRNAplot',
      version=1.1,
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
    #   packages=find_packages(),
    #     install_requires=[
    #     'matplotlib',
    #     'scikit-learn',
    #     'numpy',
    #     'svglib',
    #     'svgwrite',
    #     'pandas',
    #   ],
    #   install_requires=[
    #       'matplotlib==3.2.1',
    #       'scikit-learn==0.22.2',
    #       'numpy>=1.14.6',
    #       'svglib==1.0.0',
    #       'svgwrite==1.4',
    #       'pandas==1.0.3',
    #   ],
      scripts=[
          'src/premirnaplot.py'
      ],
      long_description=desctxt
)
