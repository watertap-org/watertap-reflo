from setuptools import setup, find_namespace_packages

setup(
    name="watertap-reflo",
    version="0.2.0rc0",
    packages=find_namespace_packages(where="src"),
    package_dir={"": "src"},
    author="WaterTAP-REFLO contributors",
    python_requires=">=3.9",
    install_requires=[
        # "watertap @ https://github.com/watertap-org/watertap/archive/main.zip", # uncomment if we need to point to main mid release cycle
        "watertap==1.4.0",
        "idaes_pse==2.8.0",
        "pyomo>=6.6.1,<6.9.3",
        "nrel-pysam>=7.0.0",
        "pint<0.25",
        "requests>=2.32",
    ],
    extras_require={
        "dev": [
            "nbsphinx",  # jupyter notebook support for sphinx
            "jinja2<3.1.0",  # see watertap-org/watertap#449
            "Sphinx",  # docs
            "sphinx_rtd_theme >=0.30",  # docs
        ]
    },
    package_data={
        "": [
            "*.yaml",
            "*.json",
            "*.pkl",
            "*.csv",
        ],
    },
)
