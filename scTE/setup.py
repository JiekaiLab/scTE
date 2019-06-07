from setuptools import setup
import glob,os

def readme():
      with open('README.rst') as f:
	       return f.read()

setup(name='scTE',
        version='1.0',
        description='Tools for estimating differential enrichment of Transposable Elements and other highly repetitive regions',
        long_description=readme(),
        classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        ],
        python_requires=">=3.6",
        keywords='..',
        url='..',
        author='..',
        author_email='he_jiangping@gibh.ac.cn',
        license='GPLv3',
        packages=[
          'scTE',
        ],
        platforms=[
          'Linux',
          'MacOS'
        ],
        install_requires=[
          'argparse',
          'pysam>=0.9',
          'pybedtools',
        ],
        include_package_data=True,
        zip_safe=False,
        scripts=[
          'bin/scTE',
        ]
        )
