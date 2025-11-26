# libmath

A minimal C math library intended for personal use in simulation/graphics code.

- Column-major layout by default. Define `ROW_MAJOR` before including the header to use row-major
- No dynamic allocations anywhere in the library
- All operations work on raw `f64*` buffers. No custom structs or composite types are provided or enforced
- Most functions follow the pattern: `fn(target, sourceA, sourceB, ...)`.  
  `target` may alias one of the sources for in-place updates.
  
## Examples

```c
f64 axis[3] = {4, 9, 5};

f64 mat3[9];
f64 result[9];

mat3_assign_identity(mat3);
mat3_rotate_axis(result, mat3, axis, 0.5 * M_PI);
```

```c
typedef struct {
    f64 X, Y, Z;
} dv3;

dv3 v1 = {4, 2, 0};
dv3 v2 = {7, 1, 3};

// overwrite v1 with the cross product
vec3_cross(&v1.X, &v1.X, &v2.X);
```
All vector and matrix operations permit in-place assignment.

## Build / Link

The library is platform and compiler agnostic.

```bash
set -e
( cd libmath && ./build.sh )
gcc -o proj main.c libmath/libmath.a
```
