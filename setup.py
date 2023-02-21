#!/usr/bin/env python
# coding:utf-8

from setuptools import find_packages, setup

setup(
    name='CASCADE',
    version='0.0.1',
    keywords=('pip','CASCADE')
    description='CASCADE: a scCAS cell type annotation method dedicated to differentiating and imbalanced types',
    long_description="CASCADE provides an effective and efficient way to automatically annotate celltypes in scCAS datasets. All CASCADE wheels distributed on PyPI are MIT licensed.",
    license='MIT License',
    url='https://github.com/BioX-NKU/CASCADE',
    author='Yuhang Jia, Siyu Li',
    packages=find_packages(),
    python_requires='>3.6.0',    
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.10',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics'    ],
     install_requires=[
        'numpy>=1.21.6,<1.22',
        'pandas>=1.4.2',
        'scipy>=1.8.1',
        'scikit-learn>=1.1.1',
        'numba>=0.55.1',
        'scanpy>=1.9.1',
        'matplotlib==3.5.2'
        'anndata>=0.8.0',
        'episcanpy>=0.3.2',
        'torch>=1.11.0',
         
      
         
    ]
)
