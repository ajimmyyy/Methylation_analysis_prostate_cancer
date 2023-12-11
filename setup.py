from setuptools import setup, find_packages

setup(
    name = "make_file",
    version = "0.1",
    packages = find_packages(),
    install_requires=[
        'pandas',
        'matplotlib',
    ],
)