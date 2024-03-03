from setuptools import setup, find_packages

setup(
    name = "methylation_analysis",
    version = "0.2",
    packages = find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'matplotlib',
        'tqdm',
        'scikit-learn',
        'seaborn',
        'scipy',
        'multipledispatch',
        "biopython",
        "scikit-learn-extra",
        "biomart",
        "goatools",
        "imblearn",
        "graphviz",
        "joblib"
    ],
)