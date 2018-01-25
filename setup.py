# -*- coding: utf-8 -*-

from distutils.core import setup

from setuptools import find_packages

from refseq_masher import __version__, program_name, program_desc

classifiers = """
Development Status :: 3 - Alpha
Environment :: Console
License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)
Intended Audience :: Science/Research
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python :: 3.5
Programming Language :: Python :: 3.6
Programming Language :: Python :: Implementation :: CPython
Operating System :: POSIX :: Linux
""".strip().split('\n')


def get_long_description():
    import codecs
    with codecs.open('README.md', encoding='UTF-8') as f:
        return f.read()


setup(
    name=program_name,
    version=__version__,
    packages=find_packages(exclude=['tests']),
    url='https://github.com/phac-nml/{}'.format(program_name),
    license='Apache v2.0',
    author='Peter Kruczkiewicz',
    author_email='peter.kruczkiewicz@gmail.com',
    description=program_desc,
    long_description=get_long_description(),
    keywords='Mash MinHash RefSeq Taxonomic Classification Containment Sequencing',
    classifiers=classifiers,
    package_dir={program_name: program_name},
    package_data={program_name: ['data/*.msh', 'data/*.csv']},
    install_requires=[
        'numpy>=1.12.1',
        'pandas>=0.20.1',
        'click',
    ],
    extras_require={
        'test': ['pytest>=3.0.7',],
    },
    entry_points={
        'console_scripts': [
            '{0}={0}.cli:cli'.format(program_name),
        ],
    }
)
