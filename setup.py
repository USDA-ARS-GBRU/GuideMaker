"""setup.py: python package setup for predictPAM

"""

from setuptools import setup


setup(
    name='predictPAM',
    version='0.1',
    packages=['predictPAM'],
    license='????',
    description='Python program to predict target sequence for CRISPR-Cas',
    long_description=open('README.rst').read(),
    classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Programming Language :: Python :: 3.6',
                 'Programming Language :: Python :: 3.5',
                 'Development Status :: 3 - Alpha'],
    keywords='CRISPR-Cas',
    url='http://example.com',
    test_suite='pytest',
    author='Ravin Poudel',
    author_email='ravin.poudel@usda.gov',
    install_requires=['biopython>=1.70'],
    python_requires='>3.7',
    tests_require=['pytest'],
    include_package_data=True,
    entry_points={'console_scripts':['predictPAM=predictPAM.main:main']},
    zip_safe=False)
