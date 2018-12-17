
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

void main(int argc, char const *argv[])
{
    char *end;
    int NUM_THREADS = strtol(argv[1], &end, 10);

    float res = 0;
    int n = 1000;

    float U[n];
    float W[n];
    for (int i = 0; i < n; i++)
    {
        W[i] = 1;
        U[i] = 12414424.268 * i;
    }

    if (NUM_THREADS == 1)
    {
        res = gm(U, W, 0, 1, 100);
    }
    printf("%10g\n", res);
}
