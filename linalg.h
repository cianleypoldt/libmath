#ifndef LIN_ALG_H
#define LIN_ALG_H

#include <stdint.h>

#ifdef LIBMATH_SINGLE_PRECISION
typedef float real_t;
#else
typedef double real_t;
#endif

typedef float    f32;
typedef double   f64;
typedef int32_t  i32;
typedef int64_t  i64;
typedef uint32_t u32;
typedef uint64_t u64;

#define PI 3.14159265358979323846

#ifdef COL_MAJOR
#        define MATNIDX(row, col, rowN, colN) ((col) * (rowN) + (row))
#else
#        define MATNIDX(row, col, rowN, colN) ((row) * (colN) + (col))
#endif
#define MAT4IDX(row, col) MATNIDX((row), (col), 4, 4)
#define MAT3IDX(row, col) MATNIDX((row), (col), 3, 3)

real_t vec3_len(const real_t* v);
real_t vec3_dot(const real_t* a, const real_t* b);
void   vec3_add(real_t* res, const real_t* a, const real_t* b);
void   vec3_scale(real_t* res, const real_t* v, real_t factor);
void   vec3_norm(real_t* res, const real_t* v);
void   vec3_cross(real_t* res, const real_t* left, const real_t* right);

void mat3_mul_vec3(real_t* res, const real_t* LHS, const real_t* v);
void vec3_outer_product(real_t* Res, const real_t* left, const real_t* right);
void vec3_to_skew_symmetric_mat3(real_t* Res, const real_t* v);

void mat3_scale_diagonal(real_t* Res, const real_t* M, real_t factor);
void mat3_add(real_t* Res, const real_t* A, const real_t* B);
void mat3_mul_scalar(real_t* Res, const real_t* M, real_t factor);
void mat3_mul(real_t* Res, const real_t* LHS, const real_t* RHS);
void mat3_rotate_X(real_t* Res, const real_t* M, real_t angle);
void mat3_rotate_Y(real_t* Res, const real_t* M, real_t angle);
void mat3_rotate_Z(real_t* Res, const real_t* M, real_t angle);
void mat3_rotate_axis(real_t*       Res,
                      const real_t* M,
                      const real_t* axis,
                      real_t        angle);

void mat4_mul_vec4(real_t* res, const real_t* LHS, const real_t* v);
void mat3_to_mat4(real_t* Res, const real_t* M);
void mat4_to_mat3(real_t* Res, const real_t* M);

void mat4_mul(real_t* Res, const real_t* LHS, const real_t* RHS);
void mat4_translate(real_t* Res, const real_t* M, const real_t* v);

void matN_assign_identity(real_t* M, u32 N);
void matN_transpose(real_t* Res, const real_t* M, u32 N);

void zeroN(real_t* res, u32 N);
void assignN(real_t* res, const real_t* v, u32 N);
void fillN(real_t* res, real_t val, u32 N);

// debug
void vecN_print(const real_t* v, u32 N);
void matN_print(const real_t* M, u32 rowN, u32 colN);
void vecN_print_named(const real_t* v, u32 N, const char* name);
void matN_print_named(const real_t* M, u32 rowN, u32 colN, const char* name);

#endif
