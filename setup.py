from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf8") as fh:
    long_description = fh.read()

with open("requirements.txt") as f:
    install_requires = f.read().splitlines()

setup(
    name="chemical_space_explorer",
    version="0.1.0",
    author="Jason Hlozek",
    author_email="jason.hlozek@uct.ac.za",
    url="https://github.com/JHlozek/chemical_space_checker",
    description="Interactive plotting package for UMAPs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    python_requires=">=3.10",
    install_requires=install_requires,
    packages=find_packages(exclude=("utilities")),
    classifiers=[
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ],
    keywords="qsar machine-learning chemistry computer-aided-drug-design",
    project_urls={"Source Code": "https://github.com/JHlozek/chemical_space_explorer"},
    include_package_data=True,
)
