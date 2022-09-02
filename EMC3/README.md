This package provides classes that represents the EMC3-EIRENE input and 
output files for an easy access with python.
Furthermore, tools for analysing simulations obtained by EMC3-EIRENE are 
provided.

Typically this package is installed as "EMC3".

The package requires additional python packages:
- FortranReader (provided via the GitLab https://src.ipp.kfa-juelich.de)
- fortranformat (https://pypi.python.org/pypi/fortranformat)

The following lines may help during the installation:

    git clone git@src.ipp.kfa-juelich.de:emc3/emc3-eirene-utilities.git EMC3
    git clone git@src.ipp.kfa-juelich.de:RACK/FortranReader.git
    wget https://pypi.python.org/packages/7f/8a/5c2361b8a45238f593ef4824b24f2dd122dd3294297424aa37486de08209/fortranformat-0.2.5.tar.gz#md5=3b8cd134f1c2cb02a7a1119a086ea7f6
    tar -xzf fortranformat-0.2.5.tar.gz
    mv fortranformat-0.2.5/fortranformat/ .
    rm -r fortranformat-0.2.5*