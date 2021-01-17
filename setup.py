from setuptools import setup

setup(
    name='bodas',
    version='0.0.0',
    author='Francesco Urbani',
    author_email='francescourbanidue@gmail.com',
    packages=['bodas'],
    scripts=[],
    url='https://github.com/urbanij/bodas',
    license='LICENSE.txt',
    description="The Asymptotic Bode Diagram in Python",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    install_requires=[
        "numpy>=1.19.2",
        "matplotlib>=3.3.2",
   ],
   keywords="bode, bode Diagram, bode plot, ",
   python_requires=">=3.7",
)
