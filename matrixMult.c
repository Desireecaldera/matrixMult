#include "matrixMult.h"

int main(int argc, char *argv[]) {
    //if (freopen(argv[1], "r", stdin) == 0)
    if (freopen("in.txt", "r", stdin) == 0) oops("Cannot open the input file.\n", -1);

    int **a1, **b1, **c1, **a2, **b2, **c2; // matrices
    int m1, k1, n1, m2, k2, n2; // dimensions of the matices m x k and k x n

    allocateAndLoadMatrices(&a1, &b1, &c1, &m1, &k1, &n1);
    allocateAndLoadMatrices(&a2, &b2, &c2, &m2, &k2, &n2);

    // the real magic happens in here
    pthread_t **tids1 = multiply(a1, b1, c1, m1, k1, n1);
    pthread_t **tids2 = multiply(a2, b2, c2, m2, k2, n2);

    free_tids(tids1, m1);
    free_tids(tids2, m2);

    // dispaly results of matrix multiplication
    printf("\nMATRIX A1\n");
    displayMatrix(a1, m1, k1);
    freeMatrix(a1, m1);
    printf("\nMATRIX B1\n");
    displayMatrix(b1, k1, n1);
    freeMatrix(b1, k1);
    printf("\nMATRIX A1 x B1\n");
    displayMatrix(c1, m1, n1);
    freeMatrix(c1, m1);

    printf("\nMATRIX A2\n");
    displayMatrix(a2, m2, k2);
    freeMatrix(a2, m2);
    printf("\nMATRIX B2\n");
    displayMatrix(b2, k2, n2);
    freeMatrix(b2, k2);
    printf("\nMATRIX A2 x B2\n");
    displayMatrix(c2, m2, n2);
    freeMatrix(c2, m2);

    return 0;
}

void *matrixThread(void *param) {
    // map the parameter onto the structure
    MATRIX_CELL *cell = (MATRIX_CELL *) param;

    //multiply the assigned col and row
    int sum = 0;
    for (int index = 0; index < cell->k; index++) {
        sum += cell->a[cell->i][index] * cell->b[index][cell->j];
    }
    //assiging result to a cell in C matrix
    // *(*(matrix +0) + 0)
    // *( *(cell->c+cell->i) + (cell->j) ) = sum;
    cell->c[cell->i][cell->j] = sum;
    //printf("ROW: %d COL: %d SUM: %d\n", cell->i, cell->j, sum);

    //TODO: implement

    free(cell);

    return NULL;
}

void allocateAndLoadMatrices(int ***a, int ***b, int ***c, int *m, int *k, int *n)
// takes pointers to two-dimensional matrices, so they can be allocated in here
// and used by the caller
// should call the function loadMatrix(int ***matrix, int m, int n)
{
    //reading in dimensions
    if (scanf("%d %d %d", m, k, n) == 0) {
        oops("Cannot read matrix sizes.\n", -2);
    }

    //allocating space for matrix A
    *a = malloc((*m * sizeof(int *)));
    for (int i = 0; i < *m; i++) {
        *((*a) + i) = malloc(*k * sizeof(int));
    }
    loadMatrix(a, *m, *k);

    //allocating space for matrix B
    *b = malloc((*k * sizeof(int *)));
    for (int j = 0; j < *k; j++) {
        *((*b) + j) = malloc(*n * sizeof(int));
    }
    loadMatrix(b, *k, *n);

    //allocating space for matrix C
    *c = malloc((*m * sizeof(int *)));
    for (int l = 0; l < *m; l++) {
        *((*c) + l) = malloc(*n * sizeof(int));
    }

    // TODO: implement(double check)
}

void loadMatrix(int ***matrix, int m, int n) {
    int test = 0;
    for (int row = 0; row < m; row++) {
        for (int col = 0; col < n; col++) {
            scanf("%d", &test);
            *(*((*matrix) + row) + col) = test;
        }
    }

    // TODO: implement(double check)
}

void freeMatrix(int **matrix, int m) {
    for (int row = 0; row < m; row++) {
        free(*(matrix + row));
    }

    // TODO: implement
}

pthread_t **multiply(int **a, int **b, int **c, int m, int k, int n) {
    pthread_t **tids = alloc_tids(m, n);
    for (int row = 0; row < m; row++) {
        for (int col = 0; col < n; col++) {
            MATRIX_CELL *p = malloc(sizeof(MATRIX_CELL));
            p->k = k;
            p->a = a;
            p->b = b;
            p->c = c;

            p->i = row;
            p->j = col;
            pthread_create((*(tids + row) + col), NULL, matrixThread, (void *) p);

        }

    }

    join(tids, m, n);

    //for **matrix
    //*(matrix + 0)  is first row
    //*(*(matrix +0) + 0) is first element of the row
    //*(*(c +i) + j) +=    *(*(a +i) + k)  times  *(*(b +k) + j)


    //TODO: implement (work on)

    return tids;
}

pthread_t **alloc_tids(int m, int n) {
    pthread_t **tids;
    //create m x n threads
    tids = malloc(m * sizeof(pthread_t *));
    for (int i = 0; i < m; i++) {
        *(tids + i) = malloc(n * sizeof(pthread_t));
    }

    // TODO: implement (double check)

    return tids;
}

void free_tids(pthread_t **tids, int m) {
    for (int row = 0; row < m; row++) {
        free(*(tids + row));
    }
    // TODO: implement
}

void join(pthread_t **tids, int m, int n) {
    for (int row = 0; row < m; row++) {
        for (int col = 0; col < n; col++) {
            pthread_join(*(*(tids + row) + col), NULL);
        }
    }
    // TODO: implement
}

void displayMatrix(int **matrix, int m, int n) {
    int test = 0;
    for (int row = 0; row < m; row++) {
        for (int col = 0; col < n; col++) {

            //for **matrix 
            //*(matrix + 0)  is first row
            //*(*(matrix +0) + 0) is first element of the row 
            test = *(*(matrix + row) + col);
            printf("%d       ", test);

        }
        printf("\n");
    }
    // TODO: implement
}
