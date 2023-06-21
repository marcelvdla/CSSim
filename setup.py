"""For package installation.

References:
https://godatadriven.com/blog/a-practical-guide-to-using-setup-py/
"""

from setuptools import find_packages, setup

def readme():
    with open('README.md', encoding='utf-8') as f:
        content = f.read()
    return content

setup(
    name="forest-ecoystem",
    packages=find_packages(),
    version='0.1.0',
    description="Implmentation of" 
        " `Networks of forest ecosystems: Mathematical modeling of their" 
        " biotic pump mechanism and resilience to certain patch deforestation`",
    long_description=readme(),
    author="Jared Frazier, Dominique Weltevreden," 
            " Eva Lampret, Marcel van de Lagemaat",
    author_email="cscidev001@gmail.com",
    license='MIT',
)
