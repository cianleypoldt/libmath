# libmath

Small, self contained math library for simulation and graphics.

- No allocations inside the math library.
- does not define or enforce composite types; Everything operates on raw `f64*`.
- Assumes column major layout. Define `ROW_MAJOR` before including if you want the opposite.
- Most functions take a target pointer as first argument.

## Examples

```C
typedef struct {
    f64 X, Y, Z;
} dv3;

dv3 v1 = {4, 2, 0};
dv3 v2 = {7, 1, 3};

vec3_cross(&v1.X, &v1.X, &v2.X); // result stored in v1
```

```C
f64 axis[3] = {4, 9, 5};

f64 mat3[9];
f64 Res[9];

mat3_assign_identity(mat3);
mat3_rotate_axis(Res, mat3, axis, 0.5 * M_PI);

```

## Builds and Link

```bash
set -e
( cd libmath && ./build.sh )
gcc -o proj main.c libmath/libmath.a
```
