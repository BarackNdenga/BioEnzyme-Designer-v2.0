#!/usr/bin/env python3
"""
BioEnzyme Designer v2.0 — Package setup.

Install:
    pip install -e .

Usage:
    python -m bioenzyme_v2 --pdb enzyme.pdb --improve activity
    python -m bioenzyme_v2 --gui
    python -m bioenzyme_v2 --api
"""

from setuptools import setup, find_packages

setup(
    name="bioenzyme_v2",
    version="2.0.0",
    description="Computational enzyme engineering toolkit with ML, "
                "multimeric analysis, ΔΔG validation, GUI, REST API, "
                "and cloud deployment.",
    author="BioEnzyme Team",
    python_requires=">=3.10",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.83",
        "scikit-learn>=1.4.0",
        "numpy>=1.26.0",
        "pandas>=2.2.0",
        "matplotlib>=3.8.0",
        "plotly>=5.20.0",
        "py3Dmol>=2.1.0",
        "requests>=2.31.0",
    ],
    extras_require={
        "api": [
            "fastapi>=0.110.0",
            "uvicorn>=0.29.0",
            "python-multipart>=0.0.9",
        ],
        "gui": [],  # tkinter is bundled with Python
        "cloud": [
            "boto3>=1.34.0",
        ],
        "all": [
            "fastapi>=0.110.0",
            "uvicorn>=0.29.0",
            "python-multipart>=0.0.9",
            "boto3>=1.34.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "bioenzyme=bioenzyme_v2.cli:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "License :: OSI Approved :: MIT License",
    ],
)
