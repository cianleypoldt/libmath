#!/usr/bin/env bash
set -e

gcc -c linalg.c -o linalg.o

ar rcs libmath.a libmath.o
