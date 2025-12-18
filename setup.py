from setuptools import setup, find_packages

setup(
    name="proteomics-explorer",
    version="1.0.0",
    author="Arina Atom",
    description="Easy access to proteomics datasets",
    url="https://github.com/arinaatom-cyber/multiomics-platform",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[
        "pandas",
        "requests",
    ],
    extras_require={
        "ui": ["ipywidgets"],
    },
)
