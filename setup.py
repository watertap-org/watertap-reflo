from setuptools import setup, find_namespace_packages

setup(
    name="watertap-seto",
    version="0.1.0.dev0",
    packages=find_namespace_packages(where="src"),
    package_dir={"": "src"},
    author="WaterTAP-SETO contributors",
    python_requires=">=3.8",
    install_requires=[
        "watertap @ https://github.com/watertap-org/watertap/archive/main.zip",
        "pytest >= 7",
        "nrel-pysam == 3.0.2"
    ],
)
