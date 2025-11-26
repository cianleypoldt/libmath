# Libmath

A small C library providing raw vector and matrix operations for simulation and rendering.

### Builds and link

''' bash
set -e

( cd libmath && ./build.sh )

gcc -o proj main.c libmath/libmath.a -lm
'''

### Syntax
