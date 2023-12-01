from setuptools import setup, find_namespace_packages

setup(
    name="watertap-reflo",
    version="0.1.0.dev0",
    packages=find_namespace_packages(where="src"),
    package_dir={"": "src"},
    author="WaterTAP-REFLO contributors",
    python_requires=">=3.8",
    install_requires=[
        "watertap @ https://github.com/watertap-org/watertap/archive/refs/tags/pr967.zip",
        # "watertap @ https://github.com/watertap-org/watertap/archive/main.zip",
        # "watertap <= 0.8",
        "pyomo <= 6.5",
        "pytest >= 7",
        "nrel-pysam == 3.0.2",
    ],
    extras_require={
        "dev": [
            "nbsphinx",  # jupyter notebook support for sphinx
            "jinja2<3.1.0",  # see watertap-org/watertap#449
            "Sphinx",  # docs
            "sphinx_rtd_theme >=0.30",  # docs
        ]
    },
)
