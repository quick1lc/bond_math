import setuptools

with open("README.txt", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bond_math",
    version="0.1.0",
    author="LC Yarnelle",
    author_email="larry.c.yarnelle@gmail.com",
    description="Help with implied horizon rates via the curve class",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/quick1lc/bond_math",
    packages=setuptools.find_packages(),
    install_requires=['pandas'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
