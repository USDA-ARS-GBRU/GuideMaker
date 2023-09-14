"""setup.py: python package setup for GuideMaker

"""

from setuptools import setup
import versioneer

requirements = [
    # package requirements go here
    'biopython>=1.81',
    'numpy >=1.11',
    'pybedtools>=0.9.1',
    'nmslib>=2.1.1',
    'pandas>=1.5',
    'pyyaml>=6.0.1',
    'regex==2020.11.13',
    'altair==4.1.0',
    'jsonschema==3.2.0',
    'streamlit>=1.26.0',
    'streamlit_tags>=1.2.8',
    'pytest>=7.4',
    'pytest-cov>=4.1',
    'pdoc3>=0.10.0',
    'onnxruntime>=1.15.1',
]



setup(
    name='guidemaker',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='GuideMaker: Software to design gRNAs pools in non-model genomes and CRISPR-Cas',
    license='CC0 1.0 Universal (CC0 1.0) Public Domain Dedication',
    author='Adam Rivers',
    author_email='adam.rivers@usda.gov',
    url='https://guidemaker.org',
    long_description=open('README.md').read(),
    packages=['guidemaker'],
    entry_points={
        'console_scripts': [
            'guidemaker=guidemaker.cli:main'
        ]
    },
    install_requires=requirements,
    python_requires='>=3.8',
    test_suite='pytest',
    tests_require=['pytest'],
    keywords='CRISPR-Cas',
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.8',
        'Development Status :: 3 - Alpha'
    ]
)



    
    
    
    
    
    
  
    
