from setuptools import setup, find_packages

setup(
    name="oset",
    version="0.1.0",
    description="An implementation of OSET in Python.",
    author="Reza Sameni",
    author_email="reza.sameni@gmail.com",
    license="BSD-3-Clause",
    packages=find_packages("src"),
    package_dir={"": "src"},
    python_requires=">=3.6",
)
