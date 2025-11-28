#include "linalg.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

void fillN(real_t* res, real_t val, u32 N) {
        for (int i = 0; i < N; i++) {
                res[i] = val;
        }
}

void zeroN(real_t* res, u32 N) {
        fillN(res, 0, N);
}

void assignN(real_t* res, const real_t* v, u32 N) {
        if (res == v) {
                return;
        }
        memcpy(res, v, N * sizeof(real_t));
}

real_t vec3_len(const real_t* v) {
        return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

real_t vec3_dot(const real_t* a, const real_t* b) {
        return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

void vec3_add(real_t* res, const real_t* a, const real_t* b) {
        res[0] = a[0] + b[0];
        res[1] = a[1] + b[1];
        res[2] = a[2] + b[2];
}

void vec3_scale(real_t* res, const real_t* v, const real_t factor) {
        res[0] = v[0] * factor;
        res[1] = v[1] * factor;
        res[2] = v[2] * factor;
}

void vec3_norm(real_t* res, const real_t* v) {
        real_t len = vec3_len(v);
        if (len != 0) {
                vec3_scale(res, v, 1 / len);
        }
}

void vec3_cross(real_t* res, const real_t* left, const real_t* right) {
        real_t tmp[3];
        tmp[0] = left[1] * right[2] - left[2] * right[1];
        tmp[1] = left[2] * right[0] - left[0] * right[2];
        tmp[2] = left[0] * right[1] - left[1] * right[0];
        assignN(res, tmp, 3);
}

void mat3_mul_vec3(real_t* res, const real_t* LHS, const real_t* v) {
        for (int i = 0; i < 3; i++) {
                res[i] = v[0] * LHS[MAT3IDX(i, 0)] + v[1] * LHS[MAT3IDX(i, 1)] +
                         v[2] * LHS[MAT3IDX(i, 2)];
        }
}

void vec3_outer_product(real_t* Res, const real_t* left, const real_t* right) {
        for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                        Res[MAT3IDX(row, col)] = left[row] * right[col];
                }
        }
}

void vec3_to_skew_symmetric_mat3(real_t* Res, const real_t* v) {
        zeroN(Res, 9);
        Res[MAT3IDX(0, 1)] = -v[2];  // Res_{1,2} = -v_z
        Res[MAT3IDX(0, 2)] = v[1];   // Res_{1,3} = v_y
        Res[MAT3IDX(1, 2)] = -v[0];  // Res_{2,3} = -v_x

        Res[MAT3IDX(1, 0)] = -Res[MAT3IDX(0, 1)];
        Res[MAT3IDX(2, 0)] = -Res[MAT3IDX(0, 2)];
        Res[MAT3IDX(2, 1)] = -Res[MAT3IDX(1, 2)];
}

void mat3_scale_diagonal(real_t* Res, const real_t* M, const real_t factor) {
        assignN(Res, M, 9);
        Res[MAT3IDX(0, 0)] = M[MAT3IDX(0, 0)] * factor;
        Res[MAT3IDX(1, 1)] = M[MAT3IDX(1, 1)] * factor;
        Res[MAT3IDX(2, 2)] = M[MAT3IDX(2, 2)] * factor;
}

void mat3_add(real_t* Res, const real_t* A, const real_t* B) {
        for (int i = 0; i < 9; i++) {
                Res[i] = A[i] + B[i];
        }
}

void mat3_mul_scalar(real_t* Res, const real_t* M, real_t factor) {
        for (int i = 0; i < 9; i++) {
                Res[i] = M[i] * factor;
        }
}

void mat3_mul(real_t* Res, const real_t* LHS, const real_t* RHS) {
        real_t tmp[9] = { 0 };
        for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                        for (int idx = 0; idx < 3; idx++) {
                                tmp[MAT3IDX(row, col)] +=
                                        LHS[MAT3IDX(row, idx)] *
                                        RHS[MAT3IDX(idx, col)];
                        }
                }
        }
        assignN(Res, tmp, 9);
}

void mat3_rotate_X(real_t* Res, const real_t* M, real_t angle) {
        real_t s = sin(angle);
        real_t c = cos(angle);
        real_t m3_rot[9];
        matN_assign_identity(m3_rot, 3);
        m3_rot[MAT3IDX(1, 1)] = c;
        m3_rot[MAT3IDX(1, 2)] = -s;
        m3_rot[MAT3IDX(2, 1)] = s;
        m3_rot[MAT3IDX(2, 2)] = c;
        mat3_mul(Res, m3_rot, M);
}

void mat3_rotate_Y(real_t* Res, const real_t* M, real_t angle) {
        real_t s = sin(angle);
        real_t c = cos(angle);
        real_t m3_rot[9];
        matN_assign_identity(m3_rot, 3);
        m3_rot[MAT3IDX(0, 0)] = c;
        m3_rot[MAT3IDX(0, 2)] = s;
        m3_rot[MAT3IDX(2, 0)] = -s;
        m3_rot[MAT3IDX(2, 2)] = c;
        mat3_mul(Res, m3_rot, M);
}

void mat3_rotate_Z(real_t* Res, const real_t* M, real_t angle) {
        real_t s = sin(angle);
        real_t c = cos(angle);
        real_t m3_rot[9];
        matN_assign_identity(m3_rot, 3);
        m3_rot[MAT3IDX(0, 0)] = c;
        m3_rot[MAT3IDX(0, 1)] = -s;
        m3_rot[MAT3IDX(1, 0)] = s;
        m3_rot[MAT3IDX(1, 1)] = c;
        mat3_mul(Res, m3_rot, M);
}

void mat3_rotate_axis(real_t*       Res,
                      const real_t* M,
                      const real_t* axis,
                      real_t        angle) {
        // v_new = v * cos(angle) + (axis x v) sin(angle) + axis(axis dot v)(1 - cos(angle))
        // to matrix form:
        // aaM = Identity * cos(angle) + (v3_to_cross_product4x4(axis)) * sin(angle) + axis*transpose(axis)(1 - cos(angle))

        real_t u[3];
        vec3_norm(u, axis);

        real_t s = sin(angle);
        real_t c = cos(angle);

        real_t identity[9];
        matN_assign_identity(identity, 3);
        mat3_mul_scalar(identity, identity, c);

        real_t cross[9];
        vec3_to_skew_symmetric_mat3(cross, u);
        mat3_mul_scalar(cross, cross, s);

        real_t dot[9];
        matN_assign_identity(dot, 3);
        vec3_outer_product(dot, u, u);
        mat3_mul_scalar(dot, dot, (1 - c));

        real_t aaM[9];
        mat3_add(aaM, identity, cross);
        mat3_add(aaM, aaM, dot);

        mat3_mul(Res, M, aaM);
}

void mat3_to_mat4(real_t* Res, const real_t* M) {
        zeroN(Res, 16);
        for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                        Res[MAT4IDX(row, col)] = M[MAT3IDX(row, col)];
                }
        }
        Res[MAT4IDX(3, 3)] = 1.0;
}

void mat4_to_mat3(real_t* Res, const real_t* M) {
        for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                        Res[MAT3IDX(row, col)] = M[MAT4IDX(row, col)];
                }
        }
}

void mat4_mul_vec4(real_t* res, const real_t* LHS, const real_t* v) {
        for (int i = 0; i < 4; i++) {
                res[i] = v[0] * LHS[MAT4IDX(i, 0)] + v[1] * LHS[MAT4IDX(i, 1)] +
                         v[2] * LHS[MAT4IDX(i, 2)] + v[3] * LHS[MAT4IDX(i, 3)];
        }
}

void mat4_mul(real_t* Res, const real_t* LHS, const real_t* RHS) {
        real_t tmp[16] = { 0 };
        for (int row = 0; row < 4; row++) {
                for (int col = 0; col < 4; col++) {
                        for (int idx = 0; idx < 4; idx++) {
                                tmp[MAT4IDX(row, col)] +=
                                        LHS[MAT4IDX(row, idx)] *
                                        RHS[MAT4IDX(idx, col)];
                        }
                }
        }
        assignN(Res, tmp, 16);
}

void mat4_translate(real_t* Res, const real_t* M, const real_t* v) {
        assignN(Res, M, 16);
        Res[MAT4IDX(0, 3)] = M[MAT4IDX(0, 3)] + v[0];
        Res[MAT4IDX(1, 3)] = M[MAT4IDX(1, 3)] + v[1];
        Res[MAT4IDX(2, 3)] = M[MAT4IDX(2, 3)] + v[2];
}

void matN_transpose(real_t* Res, const real_t* M, u32 N) {
        for (int row = 0; row < N; row++) {
                for (int col = 0; col < N; col++) {
                        Res[col * N + row] = M[row * N + col];
                }
        }
}

void matN_assign_identity(real_t* M, u32 N) {
        zeroN(M, N * N);
        for (int i = 0; i < N; i++) {
                M[MATNIDX(i, i, N, N)] = 1;
        }
}

void vecN_print(const real_t* v, u32 N) {
        printf("{ ");
        for (int i = 0; i < N - 1; i++) {
                printf("%f, ", v[i]);
        }
        printf("%f }\n", v[N - 1]);
}

void matN_print(const real_t* M, u32 rowN, u32 colN) {
        for (int row = 0; row < rowN; row++) {
                printf("{ ");
                for (int col = 0; col < colN - 1; col++) {
                        printf("%f, ", M[MATNIDX(row, col, rowN, colN)]);
                }
                printf("%f }\n", M[MATNIDX(row, colN - 1, rowN, colN)]);
        }
}

void vecN_print_named(const real_t* v, u32 N, const char* name) {
        printf("vec%u \"%s\":\n", N, name);
        vecN_print(v, N);
}

void matN_print_named(const real_t* M, u32 rowN, u32 colN, const char* name) {
        printf("mat%ux%u \"%s\":\n", rowN, colN, name);
        matN_print(M, rowN, colN);
}
