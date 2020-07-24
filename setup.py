from setuptools import setup

setup(
    name="mgxsim",
    version="0.0.1",
    description="metagenomics utils.",
    packages=["mgxsim"],
    install_requires = [
        "biopython",
        "pandas",
    ],
)

