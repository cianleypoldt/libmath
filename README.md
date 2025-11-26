# Libmath

A small C library providing raw vector and matrix operations for simulation and rendering.

## Builds and link

```bash
set -e

( cd libmath && ./build.sh )

gcc -o proj main.c libmath/libmath.a
```

## API

```C
typedef struct {
    double X, Y, Z;
} dv3;

dv3 v1 = {4, 2, 0};
dv3 v2 = {7, 1, 3};

vec3_cross(v1, v1, v2);
```

```C
double axis[3] = {4, 9, 5};
double mat3[9];
double Res[9];

mat3_assign_identity(mat3);
mat3_rotate_axis(Res, mat4, axis, 0.5 * M_PI);

```