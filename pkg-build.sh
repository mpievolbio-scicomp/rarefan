#! /bin/bash

# Build executable via pyinstaller

cp -ur app/utilities ./
ls utilities

pyinstaller -y \
            --clean \
            --workpath pkg-build \
            --distpath pkg-dist \
            --name rarefan \
            --add-data REPIN_ecology/REPIN_ecology/build/libs/REPIN_ecology.jar:./ \
            --collect-data utilities.checkers \
            -p app/utilities \
            app/utilities/rarefan_cli.py

rm -r utilities
