from setuptools import setup, find_namespace_packages

setup(
    name="watertap-reflo",
    version="0.1.0rc0",
    packages=find_namespace_packages(where="src"),
    package_dir={"": "src"},
    author="WaterTAP-REFLO contributors",
    python_requires=">=3.8",
    install_requires=[
        "watertap == 0.11",
        "pyomo==6.7.1",
        "nrel-pysam == 5.1.0",
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
