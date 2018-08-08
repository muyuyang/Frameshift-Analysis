import os
from setuptools import setup


setup(
    name = 'Frameshift_Analysis',
    description = 'Analyze frameshift mutations between genomes',
    packages = ['Frameshift_Analysis'],
    author = 'Muyu Yang',
    author_email = 'muyu.yang1127@gmail.com',
    url = 'https://github.com/muyuyang/Frameshift-Analysis',
    classifiers = [
        'Operating system :: OS independent',
        'Programming language :: Python'
        ],
    install_requires = [
        'ete3',
        'Biopython'
        ],
    )