import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Chemistry_NewtonRaphson",
    version="0.1.3",
    author="Leticia",
    author_email="leticiapequeno30@gmail.com",
    description="Equation-State Solver with Newton-Raphson method",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Leticia-maria/Chemistry_NewtonRaphson/archive/refs/tags/v_01.tar.gz",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'sympy',
        'numpy'
    ],
    keywords='Newton-Raphson',
)
