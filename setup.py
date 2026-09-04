from setuptools import setup, find_packages
setup(
    name="juxtascope",
    version="0.1.0",
    description="Post-annotation doublet/mixed-population detection for "
                "segmentation-based spatial transcriptomics",
    packages=find_packages(),
    install_requires=["scanpy>=1.9","anndata","numpy","pandas","scipy"],
    entry_points={"console_scripts":["juxtascope=juxtascope.cli:main"]},
    python_requires=">=3.9",
)
