from setuptools import setup, find_packages


REQUIREMENTS = [i.strip() for i in open("requirements.txt").readlines()]

setup(
    name = "methylation_analysis",
    version = "0.2",
    packages = find_packages(),
    install_requires = REQUIREMENTS,
)