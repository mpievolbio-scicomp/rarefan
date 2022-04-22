#! /bin/bash

# Build executable via pyinstaller

cp -ur app/utilities ./

pyinstaller -y \
            --clean \
            --workpath pkg-build \
            --distpath pkg-dist \
            --name rarefan \
            --add-data REPIN_ecology/REPIN_ecology/build/libs/REPIN_ecology.jar:./ \
            --add-data utilities/checkers.py:./ \
            utilities/rarefan_cli.py

rm -r utilities
