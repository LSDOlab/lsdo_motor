from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

    
setup(
    name='lsdo_motor',
    version='0.1.0',
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
