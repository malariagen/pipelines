import glob
import os

from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name='vector_genotype_concordance',
    version='1.0.0',
    description='Computes genotype concordance for malaria vetors',
    long_description=read('README.md'),
    long_description_content_type="text/markdown",
    packages=[],
    author='Pathogen Informatics',
    author_email='path-help@sanger.ac.uk',
    url='https://github.com/malariagen/pipelines',
    scripts=glob.glob('scripts/*'),
    setup_requires=['nose>=1.3'],
    install_requires=['wheel', 'intake', 'zarr', 'numpy', 'pandas', 'scikit-allel', 'requests', 'aiohttp', 'gcsfs'],
    test_suite='nose.collector',
    tests_require=['nose >= 1.3', 'ddt'],
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience  :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'License :: OSI Approved :: MIT License'
    ],
)
