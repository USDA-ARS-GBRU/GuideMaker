"""setup.py: python package setup for GuideMaker

"""

from setuptools import setup
import versioneer

requirements = [
    # package requirements go here
    'biopython==1.76',
    'numpy >=1.13',
    'pybedtools',
    'nmslib>=2.1.1',
    'pandas>=1.0.0',
    'pyyaml>=5.4.1',
    'regex==2020.11.13',
    'altair',
    'streamlit=0.80.0',
    'pytest>=4.6',
    'pytest-cov',
    'streamlit-tags==1.1.9',
    'pdoc3'
]


setup(
    name='guidemaker',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='GuideMAker: globally design gRNAs for any CRISPR-Cas system in any small genome',
    license='CC0 1.0 Universal (CC0 1.0) Public Domain Dedication',
    author='Adam Rivers',
    author_email='adam.rivers@usda.gov',
    url='http://tinyecology.com',
    long_description=open('README.md').read(),
    packages=['guidemaker'],
    entry_points={
        'console_scripts': [
            'guidemaker=guidemaker.cli:main'
        ]
    },
    install_requires=requirements,
    python_requires='>=3.6',
    test_suite='pytest',
    tests_require=['pytest'],
    keywords='CRISPR-Cas',
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Development Status :: 3 - Alpha'
    ]
)



    
    
    
    
    
    
  
    
