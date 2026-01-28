#!/usr/bin/env bash
set -e

gcc -Wall -Wextra -c linalg.c -lm -o libmath.o

ar rcs libmath.a libmath.o
