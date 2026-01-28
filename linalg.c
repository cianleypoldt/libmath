#include "linalg.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void zeroN(real_t* Res, u32 n) {
        memset(Res, 0, ((size_t) n) * sizeof(real_t));
}

void fillN(real_t* Res, real_t value, u32 n) {
        for (u32 i = 0; i < n; i++) {
                Res[i] = value;
        }
}

void copyN(real_t* Res, const real_t* v, u32 n) {
        if (Res == v) {
                return;
        }
        memcpy(Res, v, ((size_t) n) * sizeof(real_t));
}

real_t vec_dot(const real_t* a, const real_t* b, u32 n) {
        real_t sum = 0;
        for (u32 i = 0; i < n; i++) {
                sum += a[i] * b[i];
        }
        return sum;
}

real_t vec_length(const real_t* v, u32 n) {
        return (real_t) sqrt((double) vec_length2(v, n));
}

real_t vec_length2(const real_t* v, u32 n) {
        real_t len2 = 0;
        for (u32 i = 0; i < n; i++) {
                len2 += v[i] * v[i];
        }
        return len2;
}

void vec_normalize(real_t* Res, const real_t* v, u32 n) {
        real_t len = vec_length(v, n);
        if (len != 0) {
                vec_scale(Res, v, 1 / len, n);
                return;
        }
        zeroN(Res, n);
}

void vec_add(real_t* Res, const real_t* a, const real_t* b, u32 n) {
        for (u32 i = 0; i < n; i++) {
                Res[i] = a[i] + b[i];
        }
}

void vec_sub(real_t* Res, const real_t* a, const real_t* b, u32 n) {
        for (u32 i = 0; i < n; i++) {
                Res[i] = a[i] - b[i];
        }
}

void vec_scale(real_t* Res, const real_t* v, real_t factor, u32 n) {
        for (u32 i = 0; i < n; i++) {
                Res[i] = v[i] * factor;
        }
}

void vec_outer_product(real_t* Res, const real_t* a, const real_t* b, u32 len_a, u32 len_b) {
        for (u32 row = 0; row < len_a; row++) {
                for (u32 col = 0; col < len_b; col++) {
                        Res[MAT_IDX(row, col, len_a, len_b)] = a[row] * b[col];
                }
        }
}

real_t mat_trace(real_t* M, u32 n) {
        real_t trace = 0;
        for (u32 i = 0; i < n; i++) {
                trace += M[MAT_IDX(i, i, n, n)];
        }
        return trace;
}

real_t mat_determinant(real_t* M, u32 n) {  // upgrade to LU decomposition in future
        if (n == 1) {
                return M[0];
        } else if (n == 2) {
                return M[MAT_IDX(0, 0, n, n)] * M[MAT_IDX(1, 1, n, n)] -
                       M[MAT_IDX(0, 1, n, n)] * M[MAT_IDX(1, 0, n, n)];
        }

        real_t det = 0;

        for (u32 col = 0; col < n; col++) {
                // build minor matrix
                real_t* minor = (real_t*) malloc((n - 1) * (n - 1) * sizeof(real_t));
                for (u32 i = 1; i < n; i++) {
                        u32 minor_col = 0;
                        for (u32 j = 0; j < n; j++) {
                                if (j == col) {
                                        continue;
                                }
                                minor[MAT_IDX(i - 1, minor_col, n - 1, n - 1)] = M[MAT_IDX(i, j, n, n)];
                                minor_col++;
                        }
                }

                // Recursive expansion
                real_t sign = (col % 2 == 0) ? 1.0 : -1.0;
                det += sign * M[MAT_IDX(0, col, n, n)] * mat_determinant(minor, n - 1);

                free(minor);
        }

        return det;
}

void mat_transpose(real_t* Res, const real_t* A, u32 rowsA, u32 colsA) {
        real_t* tmp = Res == A ? malloc((size_t) rowsA * (size_t) colsA * sizeof(real_t)) : Res;

        for (u32 col = 0; col < rowsA; col++) {
                for (u32 row = 0; row < colsA; row++) {
                        tmp[MAT_IDX(row, col, colsA, rowsA)] = A[MAT_IDX(col, row, rowsA, colsA)];
                }
        }
        if (tmp != Res) {
                copyN(Res, tmp, rowsA * colsA);
                free(tmp);
        }
}

void mat_identity(real_t* Res, u32 n) {
        zeroN(Res, n * n);
        for (u32 i = 0; i < n; i++) {
                Res[MAT_IDX(i, i, n, n)] = 1;
        }
}

void mat_add(real_t* Res, const real_t* A, const real_t* B, u32 rows, u32 cols) {
        vec_add(Res, A, B, rows * cols);
}

void mat_sub(real_t* Res, const real_t* A, const real_t* B, u32 rows, u32 cols) {
        vec_sub(Res, A, B, rows * cols);
}

void mat_mul(real_t* Res, const real_t* A, const real_t* B, u32 rowsA, u32 colsA, u32 colsB) {
        real_t* tmp = (Res == A || Res == B) ? malloc(rowsA * colsB * sizeof(real_t)) : Res;

        for (u32 row = 0; row < rowsA; row++) {
                for (u32 col = 0; col < colsB; col++) {
                        real_t sum = 0;
                        for (u32 i = 0; i < colsA; i++) {
                                sum += A[MAT_IDX(row, i, rowsA, colsA)] * B[MAT_IDX(i, col, colsA, colsB)];
                        }
                        tmp[MAT_IDX(row, col, rowsA, colsB)] = sum;
                }
        }
        if (tmp != Res) {
                copyN(Res, tmp, rowsA * colsB);
                free(tmp);
        }
}

void mat_vec_mul(real_t* Res, const real_t* A, const real_t* v, u32 rowsA, u32 colsA) {
        zeroN(Res, rowsA);
        for (u32 row = 0; row < rowsA; row++) {
                for (u32 col = 0; col < colsA; col++) {
                        Res[row] += A[MAT_IDX(row, col, rowsA, colsA)] * v[col];
                }
        }
}

void mat_scalar_mul(real_t* Res, const real_t* M, real_t factor, u32 rows, u32 cols) {
        vec_scale(Res, M, factor, rows * cols);
}

void mat_plane_proj(real_t* Res, const real_t* u, const real_t* v, u32 n) {
        real_t* M = (real_t*) malloc((size_t) n * (size_t) n * sizeof(real_t));
        vec_length(u, n);
}

void mat_plane_rotation(real_t* Res, const real_t* u, const real_t* v, real_t angle, u32 n) {
        size_t  mat_size = (size_t) n * (size_t) n;
        real_t* buffer   = (real_t*) malloc(2 * mat_size * sizeof(real_t));
        real_t* M1       = buffer + 0 * mat_size;
        real_t* M2       = buffer + 1 * mat_size;
        real_t* M3       = Res;

        real_t* uv = M1;
        real_t* vu = M2;
        vec_outer_product(uv, u, v, n, n);
        vec_outer_product(vu, v, u, n, n);

        real_t* K = M1;
        mat_sub(K, vu, uv, n, n);

        real_t* vv = M2;
        real_t* uu = M3;
        vec_outer_product(vv, v, v, n, n);
        vec_outer_product(uu, u, u, n, n);

        real_t* P = M2;
        mat_add(P, vv, uu, n, n);

        mat_identity(Res, n);

        mat_scalar_mul(K, K, sin(angle), n, n);
        mat_scalar_mul(P, P, cos(angle) - 1.0, n, n);

        mat_add(Res, Res, K, n, n);
        mat_add(Res, Res, P, n, n);

        free(buffer);
}

void mat3_X_rotation(real_t* Res, real_t angle) {
        real_t s = sin(angle);
        real_t c = cos(angle);
        mat_identity(Res, 3);
        Res[MAT_IDX(1, 1, 3, 3)] = c;
        Res[MAT_IDX(1, 2, 3, 3)] = -s;
        Res[MAT_IDX(2, 1, 3, 3)] = s;
        Res[MAT_IDX(2, 2, 3, 3)] = c;
}

void mat3_Y_rotation(real_t* Res, real_t angle) {
        real_t s = sin(angle);
        real_t c = cos(angle);
        mat_identity(Res, 3);
        Res[MAT_IDX(0, 0, 3, 3)] = c;
        Res[MAT_IDX(0, 2, 3, 3)] = s;
        Res[MAT_IDX(2, 0, 3, 3)] = -s;
        Res[MAT_IDX(2, 2, 3, 3)] = c;
}

void mat3_Z_rotation(real_t* Res, real_t angle) {
        real_t s = sin(angle);
        real_t c = cos(angle);
        mat_identity(Res, 3);
        Res[MAT_IDX(0, 0, 3, 3)] = c;
        Res[MAT_IDX(0, 1, 3, 3)] = -s;
        Res[MAT_IDX(1, 0, 3, 3)] = s;
        Res[MAT_IDX(1, 1, 3, 3)] = c;
}

void mat3_axis_rotation(real_t* Res, const real_t* v, real_t angle) {
        real_t s = sin(angle);
        real_t c = cos(angle);

        real_t u_cross[9];
        vec3_to_cross_mat3(u_cross, v);
        mat_scalar_mul(u_cross, u_cross, s, 3, 3);

        real_t uu[9];
        vec_outer_product(uu, v, v, 3, 3);
        mat_scalar_mul(uu, uu, (1 - c), 3, 3);

        mat_identity(Res, 3);
        mat_scalar_mul(Res, Res, c, 3, 3);

        mat_add(Res, Res, u_cross, 3, 3);
        mat_add(Res, Res, uu, 3, 3);
}

void vec3_to_cross_mat3(real_t* Res, const real_t* v) {
        memset(Res, 0, 9 * sizeof(real_t));
        Res[MAT_IDX(0, 1, 3, 3)] = -v[2];
        Res[MAT_IDX(0, 2, 3, 3)] = v[1];
        Res[MAT_IDX(1, 2, 3, 3)] = -v[0];

        Res[MAT_IDX(1, 0, 3, 3)] = v[2];
        Res[MAT_IDX(2, 0, 3, 3)] = -v[1];
        Res[MAT_IDX(2, 1, 3, 3)] = v[0];
}

// debug
void vec_print(const real_t* v, u32 n) {
        printf("{ ");
        for (u32 i = 0; i < n - 1; i++) {
                printf("%f, ", v[i]);
        }
        printf("%f }\n", v[n - 1]);
}

void vec_print_named(const real_t* v, u32 n, const char* name) {
        printf("vec%i \"%s\":\n", n, name);
        vec_print(v, n);
}

void mat_print(const real_t* M, u32 rows, u32 cols) {
        for (u32 row = 0; row < rows; row++) {
                printf("{ ");
                for (u32 col = 0; col < cols - 1; col++) {
                        printf("%f, ", M[MAT_IDX(row, col, rows, cols)]);
                }
                printf("%f }\n", M[MAT_IDX(row, cols - 1, rows, cols)]);
        }
}

void mat_print_named(const real_t* M, u32 rows, u32 cols, const char* name) {
        printf("mat%ix%i \"%s\":\n", rows, cols, name);
        mat_print(M, rows, cols);
}
