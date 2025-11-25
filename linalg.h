#ifndef LIN_ALG_H
#define LIN_ALG_H

#include <stdint.h>

typedef double   f64;
typedef int32_t  i32;
typedef int64_t  i64;
typedef uint32_t u32;
typedef uint64_t u64;

#ifndef M_PI
#        define M_PI 3.14159265359
#endif

#ifndef COL_MAJOR
#        define MATNIDX(row, col, rowN, colN) ((col) * rowN + (row))
#else
#        define MATNIDX(row, col, rowN, colN) ((row) * colN + (row))
#endif
#define MAT4IDX(row, col) MATNIDX((row), (col), 4, 4)
#define MAT3IDX(row, col) MATNIDX((row), (col), 3, 3)

void clearN(f64* res, u32 N);
void assignN(f64* res, const f64* v, u32 N);
void fillN(f64* res, f64 val, u32 N);

f64  len3(const f64* v);
f64  dot3(const f64* a, const f64* b);
void add3(f64* res, const f64* a, const f64* b);
void scale3(f64* res, const f64* v, f64 factor);
void norm3(f64* res, const f64* v);
void cross3(f64* res, const f64* left, const f64* right);

void mul3x3vec3(f64* res, const f64* LHS, const f64* v);
void outer_product3(f64* Res, const f64* left, const f64* right);
void v3to_cross_product3x3(f64* Res, const f64* v);

void scale3x3_diagonal(f64* Res, const f64* M, f64 factor);
void add3x3(f64* Res, const f64* A, const f64* B);
void mul3x3scalar(f64* Res, const f64* M, f64 factor);
void mul3x3(f64* Res, const f64* LHS, const f64* RHS);
void rotate3x3X(f64* Res, const f64* M, f64 angle);
void rotate3x3Y(f64* Res, const f64* M, f64 angle);
void rotate3x3Z(f64* Res, const f64* M, f64 angle);
void rotate3x3_axis(f64* Res, const f64* M, const f64* axis, f64 angle);

void mul4x4vec4(f64* res, const f64* LHS, const f64* v);
void demote4x4to3x3(f64* Res, const f64* M);
void promote3x3to4x4(f64* Res, const f64* M);

void mul4x4(f64* Res, const f64* LHS, const f64* RHS);
void translate4x4(f64* Res, const f64* M, const f64* v);

void assign_identityNxN(f64* M, u32 N);
void transposeNxN(f64* Res, const f64* M, u32 N);

// debug
void printN(const f64* v, u32 N);
void printNxN(const f64* M, u32 rowN, u32 colN);
void printN_named(const f64* v, u32 N, const char* name);
void printNxN_named(const f64* M, u32 rowN, u32 colN, const char* name);

#endif
