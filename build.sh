#!/usr/bin/env bash
set -e

gcc -c linalg.c -lm -o libmath.o

ar rcs libmath.a libmath.o
