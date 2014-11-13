#Population Genetics project

Analyzing our dataset using dadi. Two scripts so far, one to convert
Mike's datafile to the format that dadi wants, and one to start doing our
analysis.

Usage:

    python makeSnpFile.py inputFile.out --dadi > outputFile.txt

Or, for help:

    python makeSnpFile.py --help

This will give you an `outputFile.txt` which you can use to generate an allele
frequency spectrum using dadi.
