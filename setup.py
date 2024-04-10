from setuptools import setup, find_packages

setup(
    name = "methylation_analysis",
    version = "0.2",
    packages = find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'matplotlib',
        'seaborn',
        'tqdm',
        'scikit-learn',
        "scikit-learn-extra",
        "xgboost",
        "imblearn",
        'scipy',
        'multipledispatch',
        "biopython",
        "biomart",
        "goatools",
        "joblib",
        "shap"
    ],
)