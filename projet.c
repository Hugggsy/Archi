
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <immintrin.h>

struct thread_args
{
    float *U;
    float *V;
    float a;
    int k;
    int n;
    int mode;
};

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

float gen_gm(void *thread_args)
{
    float res;
    struct thread_args *arg = (struct thread_args *)thread_args;
    float *U = arg->U;
    float *V = arg->V;
    float a = arg->a;
    float k = arg->k;
    int n = arg->n;
    int mode = arg->mode;
    if (mode == 0)
    {
        res = gm(U, V, a, k, n);
    }
    else
    {
        res = vect_gm(U, V, a, k, n);
    }
    return res;
}

float parallel_gm(float *U, float *W, float a, int k, int n, int mode, int nb_threads)
{
    pthread_t thread[nb_threads];
    pthread_attr_t attr;
    int error_code;
    long t;
    int size_of_partition = n / nb_threads;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for (t = 0; t < nb_threads; t++)
    {
        // Create static partition of U and W
        float partition_of_u[size_of_partition];
        float partition_of_v[size_of_partition];
        for (int i = size_of_partition * t; i < size_of_partition * (t + 1); i++)
        {
            partition_of_u[i] = U[i];
            partition_of_v[i] = W[i];
        }
        // Create struct for passing args to compute function
        struct thread_args args;
        args.U = partition_of_u;
        args.V = partition_of_v;
        args.a = a;
        args.n = size_of_partition;
        args.mode = mode;
        printf("Main: creating thread %ld\n", t);
        error_code = pthread_create(&thread[t], &attr, gen_gm, (void *)&args);
        if (error_code)
        {
            printf("ERROR; return code from pthread_create() is %d\n", error_code);
            exit(-1);
        }
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
        res = parallel_gm(U, W, 0, 1, n, 0, 1);
    }
    printf("%10g\n", res);
}
