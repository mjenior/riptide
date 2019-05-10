#!/usr/bin/python

from distutils.core import setup

setup(name='riptide',
	version='1.7.1', 
	description='Reaction Inclusion by Parsimony and Transcript Distribution (RIPTiDe)',
	author='Matthew Jenior',
	author_email='mattjenior@gmail.com',
	url='https://github.com/mjenior/riptide',
	packages=['riptide'],
    install_requires=[
          'cobra',
          'symengine',
    ],
    license='MIT',
    long_description=open('README.rst').read()
	)
