#ifndef C_UTILS
#define C_UTILS

#include <stdint.h>

typedef float    f32;
typedef double   f64;
typedef int32_t  i32;
typedef int64_t  i64;
typedef uint32_t u32;
typedef uint64_t u64;

typedef f64 real_t;

#define PI 3.14159265358979323846

#ifdef COL_MAJOR
#        define MAT_IDX(row, col, rows, cols) ((col) * (rows) + (row))
#else
#        define MAT_IDX(row, col, rows, cols) ((row) * (cols) + (col))
#endif
#define MAT3IDX(row, col) MATNIDX((row), (col), 3, 3)
#define MAT4IDX(row, col) MATNIDX((row), (col), 4, 4)

void fillN(real_t* Res, real_t value, u32 n);
void zeroN(real_t* Res, u32 n);
void copyN(real_t* Res, const real_t* v, u32 n);

real_t vec_dot(const real_t* a, const real_t* b, u32 n);
real_t vec_length(const real_t* v, u32 n);
real_t vec_length2(const real_t* v, u32 n);
void   vec_add(real_t* Res, const real_t* a, const real_t* b, u32 n);
void   vec_sub(real_t* Res, const real_t* a, const real_t* b, u32 n);
void   vec_scale(real_t* Res, const real_t* v, real_t factor, u32 n);
void   vec_normalize(real_t* Res, const real_t* v, u32 n);

void vec_outer_product(real_t* Res, const real_t* a, const real_t* b, u32 n_b, u32 n_a);

void mat_identity(real_t* Res, u32 n);
void mat_transpose(real_t* Res, const real_t* M, u32 rowsA, u32 colsA);
void mat_resize(real_t* Res, const real_t* M, u32 rowsM, u32 colsM, u32 rowsRes, u32 colsRes);
void mat_add(real_t* Res, const real_t* A, const real_t* B, u32 rows, u32 cols);
void mat_sub(real_t* Res, const real_t* A, const real_t* B, u32 rows, u32 cols);
void mat_mul(real_t* Res, const real_t* A, const real_t* B, u32 rowsA, u32 colsA, u32 colsB);
void mat_vec_mul(real_t* Res, const real_t* M, const real_t* v, u32 n, u32 cols);
void mat_scalar_mul(real_t* Res, const real_t* M, real_t factor, u32 rows, u32 cols);
void mat_plane_rotation(real_t* Res, const real_t* u, const real_t* v, real_t angle, u32 n);  // |u| = |v| = 1  ;   u * v = 0

void cross_vec3(real_t* Res, const real_t* a, const real_t* b);
void mat3_line_proj(real_t* Res, const real_t* a);
void mat3_plane_proj(real_t* Res, const real_t* a);
void mat3_X_rotation(real_t* Res, real_t angle);
void mat3_Y_rotation(real_t* Res, real_t angle);
void mat3_Z_rotation(real_t* Res, real_t angle);
void mat3_axis_rotation(real_t* Res, const real_t* v, real_t angle);  // |v| = 1

void vec3_to_cross_mat3(real_t* Res, const real_t* v);

void mat4_identity(real_t* Res);
void mat4_mul(real_t* Res, real_t* A, real_t* B);
void mat4_perspective_proj(real_t* Res, real_t fov, real_t aspect, real_t near_plane, real_t far_plane);

// debug
void vec_print(const real_t* v, u32 n);
void vec_print_named(const real_t* v, u32 n, const char* name);
void mat_print(const real_t* M, u32 rows, u32 cols);
void mat_print_named(const real_t* M, u32 rows, u32 cols, const char* name);

#endif
