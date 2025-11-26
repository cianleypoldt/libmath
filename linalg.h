#ifndef LIN_ALG_H
#define LIN_ALG_H

#include <stdint.h>

typedef double   f64;
typedef int32_t  i32;
typedef int64_t  i64;
typedef uint32_t u32;
typedef uint64_t u64;

#ifndef M_PI
#        define M_PI 3.14159265358979323846
#endif

#ifndef COL_MAJOR
#        define MATNIDX(row, col, rowN, colN) ((col) * rowN + (row))
#else
#        define MATNIDX(row, col, rowN, colN) ((row) * colN + (col))
#endif
#define MAT4IDX(row, col) MATNIDX((row), (col), 4, 4)
#define MAT3IDX(row, col) MATNIDX((row), (col), 3, 3)

f64  vec3_len(const f64* v);
f64  vec3_dot(const f64* a, const f64* b);
void vec3_add(f64* res, const f64* a, const f64* b);
void vec3_scale(f64* res, const f64* v, f64 factor);
void vec3_norm(f64* res, const f64* v);
void vec3_cross(f64* res, const f64* left, const f64* right);

void mat3_mul_vec3(f64* res, const f64* LHS, const f64* v);
void vec3_outer_product(f64* Res, const f64* left, const f64* right);
void vec3_to_skew_symmetric_mat3(f64* Res, const f64* v);

void mat3_scale_diagonal(f64* Res, const f64* M, f64 factor);
void mat3_add(f64* Res, const f64* A, const f64* B);
void mat3_mul_scalar(f64* Res, const f64* M, f64 factor);
void mat3_mul(f64* Res, const f64* LHS, const f64* RHS);
void mat3_rotate_X(f64* Res, const f64* M, f64 angle);
void mat3_rotate_Y(f64* Res, const f64* M, f64 angle);
void mat3_rotate_Z(f64* Res, const f64* M, f64 angle);
void mat3_rotate_axis(f64* Res, const f64* M, const f64* axis, f64 angle);

void mat4_mul_vec4(f64* res, const f64* LHS, const f64* v);
void mat3_to_mat4(f64* Res, const f64* M);
void mat4_to_mat3(f64* Res, const f64* M);

void mat4_mul(f64* Res, const f64* LHS, const f64* RHS);
void mat4_translate(f64* Res, const f64* M, const f64* v);

void matN_assign_identity(f64* M, u32 N);
void matN_transpose(f64* Res, const f64* M, u32 N);

void zeroN(f64* res, u32 N);
void assignN(f64* res, const f64* v, u32 N);
void fillN(f64* res, f64 val, u32 N);

// debug
void vecN_print(const f64* v, u32 N);
void matN_print(const f64* M, u32 rowN, u32 colN);
void vecN_print_named(const f64* v, u32 N, const char* name);
void matN_print_named(const f64* M, u32 rowN, u32 colN, const char* name);

#endif
