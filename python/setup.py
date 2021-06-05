from setuptools import setup, find_packages

setup(
    name="eos",
    version="0.1.0",
    install_requires=[
        "scipy",
        "numpy"
    ],
    author="Sho Hirose",
    author_email="sho.hirose@gmail.com",
    description="Package for cubic equation of states",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
    python_requires='>=3.6'
    )