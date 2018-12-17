
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <immintrin.h>

float gm(float *U, float *W, float a, int k, int n)
{
    float res = 0;
    float sum_of_weights = 0;
    float weighted_sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum_of_weights += W[i];
        weighted_sum += powf(W[i] * U[i] - a, k);
    }
    return weighted_sum / sum_of_weights;
}

float vect_gm(float *U, float *W, float a, int k, int n)
{
}

float gen_gm(float *U, float *W, float a, int k, int n, int mode)
{
    float res;
    if (mode == 0)
    {
        res = gm(U, W, a, k, n);
    }
    else
    {
        res = vect_gm(U, W, a, k, n);
    }
    return res;
}

float parallel_gm(float *U, float *W, float a, int k, int n, int mode, int nb_threads)
{
    if (mode == 0)
    {
    }
}
void main(int argc, char const *argv[])
{
    char *end;
    int NUM_THREADS = strtol(argv[1], &end, 10);

    float res = 0;
    int n = 100000;

    float U[n];
    float W[n];
    for (int i = 0; i < n; i++)
    {
        W[i] = 1;
        U[i] = 10;
    }

    if (NUM_THREADS == 1)
    {
        res = gm(U, W, 0, 1, n);
    }
    printf("%10g\n", res);
}
