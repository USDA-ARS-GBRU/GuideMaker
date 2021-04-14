"""setup.py: python package setup for GuideMaker

"""

from setuptools import setup
import versioneer

requirements = [
    # package requirements go here
    'biopython==1.76',
    'numpy >=1.11',
    'pybedtools',
    'nmslib>=2.0.4',
    'pandas>=1.0.0',
    'pyyaml==5.4',
    'regex==2020.11.13',
    'altair',
    'streamlit',
    'pytest>=4.6',
    'pytest-cov',
    'streamlit-tags',
    'pdoc3'
]


setup(
    name='guidemaker',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=['guidemaker'],
    license='CC0 1.0 Universal (CC0 1.0) Public Domain Dedication',
    description='GuideMAker: globally design gRNAs for any CRISPR-Cas system in any small genome',
    long_description=open('README.md').read(),
    classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Programming Language :: Python :: 3.7',
                 'Programming Language :: Python :: 3.8',
                 'Development Status :: 3 - Alpha'],
    keywords='CRISPR-Cas',
    url='http://tinyecology.com',
    test_suite='pytest',
    author='Adam Rivers',
    author_email='adam.rivers@usda.gov',
    install_requires=requirements,
    test_suite='pytest',
    tests_require=['pytest'],
    include_package_data=True,
    entry_points={'console_scripts':['guidemaker = guidemaker.cli:main',]},
    zip_safe=False)
