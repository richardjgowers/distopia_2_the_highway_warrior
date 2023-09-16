#include <benchmark/benchmark.h>
#include <iostream>
#include <random>

#include "qcprot.h"

namespace QCProt {
    double **MatInit(const int rows, const int cols)
    {
        int             i;
        double        **matrix = NULL;
        double         *matspace = NULL;

        matspace = (double *) calloc((rows * cols), sizeof(double));
        if (matspace == NULL)
        {
            perror("\n ERROR");
            printf("\n ERROR: Failure to allocate matrix space in MatInit(): (%d x %d)\n", rows, cols);
            exit(EXIT_FAILURE);
        }

        /* allocate room for the pointers to the rows */
        matrix = (double **) malloc(rows * sizeof(double *));
        if (matrix == NULL)
        {
            perror("\n ERROR");
            printf("\n ERROR: Failure to allocate room for row pointers in MatInit(): (%d)\n", rows);
            exit(EXIT_FAILURE);
        }

        /*  now 'point' the pointers */
        for (i = 0; i < rows; i++)
            matrix[i] = matspace + (i * cols);

        return(matrix);
    }
}


void BM_qcprot(benchmark::State &state) {
    double **frag_a, **frag_b;
    double rotmat[9];
    double weight[10];

    frag_a = QCProt::MatInit(3, 10);
    frag_b = QCProt::MatInit(3, 10);

    for (auto _ : state) {
        CalcRMSDRotationalMatrix(frag_a, frag_b, 10, rotmat, weight);
    }
}


BENCHMARK(BM_qcprot);

BENCHMARK_MAIN();
