#!/usr/bin/env python
from setuptools import setup, find_packages


install_requires = [
    # 'collections', 
    # 'itertools', ### default pack in python3
    # 'subprocess', 
    # 'multiprocessing',  ### only python2
    # 'regex',
    'numpy>=1.12.1', 
    'pyfaidx>=0.5.3.1', 
    'pysam>=0.11.2.2',
    'pysamstats>=1.0.1', 
    'scipy>=0.19.0',
    'pandas>=0.20.1'
]

tests_require = [
    'nose',
    'mock'
]

extras_require = {
    'docs': [
        'Sphinx>=1.1', 
        'numpydoc>=0.5'
    ]
}


setup(name='mosaicpc',
    version='0.1.2',
    url='https://github.com/parklab/MosaicPC-Calico',
    license='MIT',
    author='Yanmei Dou, Minseok Kwon',
    author_email='daniel.minseok.kwon@gmail.com',
    description='Detection of mosaic variants',
    keywords=['genomics', 'bioinformatics'],
    classifiers=[
        'Operating System :: OS Independent',
        'Topic :: Software Development :: Libraries',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    packages=find_packages(exclude=['tests']),
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    zip_safe=False,
    install_requires=install_requires,
    setup_requires=['nose>=1.0'],
    test_suite='nose.collector',
    entry_points={
        'console_scripts': [
            'mosaicpc=mosaicpc.cli:cli',
        ]
    })