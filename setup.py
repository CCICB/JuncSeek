from setuptools import setup, find_packages

setup(
    name="juncseq",                   # Package name
    version="0.1.0",                     # Initial version
    author="Patricia Sullivan",                  # Your name
    author_email="psullivan@ccia.org.au",
    description="Detect and analyze aberrant splicing events using splice junction data",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/CCICB/juncseq",  # Repository URL
    packages=find_packages(),            # Automatically find package directories
    classifiers=[                        # Metadata for PyPI
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',             # Minimum Python version
    install_requires=[                   # Dependencies for your project
        # Add required libraries here
    ],
    entry_points={
        "console_scripts": [
            "juncseq=juncseq:main",    # Maps 'juncseq' command to the main() function in juncseq.py
        ],
    },
)