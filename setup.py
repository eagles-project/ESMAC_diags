from setuptools import setup, find_packages

setup(
        name="esmac_diags",
        version = "2.0.0",
        license="BSD",
        packages=find_packages(
            where='src'
            ),
        package_dir={"": "src"},
        install_requires=[
            "numpy",
            "pytest",
            "matplotlib",
            "netCDF4",
            "pip"]
        )

