from setuptools import setup, find_packages
import os.path
import codecs

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
              delim = '"' if '"' in line else "'"
              return line.split(delim)[1]


with open("README.md", "r") as fh:
    long_description = fh.read()
  
setup(
    name='lsdo_motor',
    version=get_version("lsdo_motor/__init__.py"),
    author='suspensemax',
    author_email='swaminathanrajashekar@gmail.com',
    license='MIT',
    keywords='python motor PMSM',
    url='http://github.com/LSDOlab/suspensemax/lsdo_motor',
    description='The model is valid for PMSM',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(),
    python_requires='>=3.7',
    platforms=['any'],
    install_requires=[
        'numpy',
        'csdl',
        'modopt',
        'python_csdl_backend',
        'scipy',
        'seaborn',
        'matplotlib',
    ],
    
)
