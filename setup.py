#!/usr/bin/python

from distutils.core import setup

setup(name='riptide',
	version='1.5', 
	description='Reaction Inclusion by Parsimony and Transcript Distribution (RIPTiDe)',
	author='Matt Jenior',
	author_email='mattjenior@gmail.com',
	url='https://github.com/mjenior/riptide',
	packages=['riptide'],
    license='MIT',
    long_description=open('README.md').read()
	)
