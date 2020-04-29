"""setup.py: python package setup for Guidefinder

"""

from setuptools import setup


setup(
    name='guidefinder',
    version='0.1',
    packages=['guidefinder'],
    license='????',
    description='GuideFinder: globally design gRNAs for any CRISPR-Cas system in any small genome',
    long_description=open('README.rst').read(),
    classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Programming Language :: Python :: 3.7',
                 'Programming Language :: Python :: 3.8',
                 'Development Status :: 3 - Alpha'],
    keywords='CRISPR-Cas',
    url='http://example.com',
    test_suite='pytest',
    author='Adam Rivers',
    author_email='adam.rivers@usda.gov',
    install_requires=['biopython>=1.74', 'pybedtools>=0.8.0', 'nmslib>=2.0.6', 'pandas>=1.0.3'],
    python_requires='>=3.7',
    tests_require=['pytest'],
    include_package_data=True,
    entry_points={'console_scripts':['guidefinder=guidefinder.cli:main']},
    zip_safe=False)
