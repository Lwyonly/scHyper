from setuptools import setup, find_packages
import glob

files = glob.glob("scHyper/**/*.*", recursive=True)
files = [path.replace("\\", "/") for path in files]
files = [path.replace("scHyper/", "") for path in files]

with open("README.rst", "rt", encoding="utf-8") as fh:
    long_description_readme_md = fh.read()

setup(
    name='scHyper',
    author='Li',
    author_email='contact@author.com',
    version='1.0.13',
    description='scHyper: Reconstructing cell-cell communication through hypergraph neural networks',
    long_description=long_description_readme_md,
    long_description_content_type='text/markdown',
    packages=find_packages(),
    url='https://github.com/Lwyonly/scHyper',
    license='MIT',
    install_requires=[
        'numpy>=1.24.3',
        'scipy==1.11.2',
        'pandas==2.1.4',
        'torch>=2.0.1',
        'scikit-learn>=1.3.0',
        'tqdm',
        'h5py',
    ],
    Requires_Python='>=3.9',
)
