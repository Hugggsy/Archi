/***********************
 * Projet ES Architecture matérielle et logiciel des ordinateurs
 * Auteurs: H. Depretz et Y. Pradat
 * Décembre 2018
 * **********************/

#define N_FLOAT 4
#define N_ALIGN 16
#include <sys/time.h> // for timing
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <immintrin.h>
#include <xmmintrin.h>
#include <pmmintrin.h>

struct thread_res
{
    float sum_of_weights;
    float sum_of_weighted_diff;
};

struct thread_args
{
    float *U;
    float *W;
    float a;
    int k;
    int n;
    int mode;
    struct thread_res res;
};

double now()
{
    // Retourne l'heure actuelle en secondes
    struct timeval t;
    double f_t;
    gettimeofday(&t, NULL);
    f_t = t.tv_usec;
    f_t = f_t / 1000000.0;
    f_t += t.tv_sec;
    return f_t;
}
//INITIALISATION DES VECTEURS U ET W
void init(float *U, float *W, unsigned int N)
{
    unsigned int i;
    for (i = 0; i < N; i++)
    {
        U[i] = (float)rand() / RAND_MAX / N;
        W[i] = (float)rand() / RAND_MAX / N;
    }
}

//CALCULE LA SOMME D'UN VECTEUR DE TAILLE N
float sum(float *W, int N)
{
    unsigned int i;
    float s = 0;
    for (i = 0; i < N; i++)
        s += W[i];
    return s;
}

//REMPLIT LE VECTEUR A PAR LA VALEUR a
void init_A(float a, float *A)
{
    unsigned int i;
    for (i = 0; i < N_FLOAT; i++)
        A[i] = a;
}
//VERSION VECTORIELLE
void run_sse(float *U, float *W, float *A, float *RV, int k, int N)
{
    unsigned int i;
    unsigned int j;
    __m128 mm_U, mm_W, mm_A, mm_t, mm_p; //On déclare cinq registres vectoriels
    mm_A = _mm_load_ps(&A[0]);

    for (i = 0; i < N; i += 4)
    {
        mm_p = _mm_set1_ps(1.0f); // Pour calculer la puissance k
        mm_U = _mm_load_ps(&U[i]); //On charge (U[4i], U[4i+1], U[4i+2], U[4i+3])
        mm_W = _mm_load_ps(&W[i]); //On charge (W[4i], W[4i+1], W[4i+2], W[4i+3])
        mm_t = _mm_sub_ps(_mm_mul_ps(mm_U, mm_W), mm_A);
        for (j = 0; j < k; j++)
        { //On calcule la puissance k
            mm_p = _mm_mul_ps(mm_p, mm_t);
        }
        _mm_store_ps(&RV[i], mm_p); //On stocke  (p[4i], p[4i+1], p[4i+2], p[4i+3])
    }
}

//EXECUTE CALCUL SCALAIRE ET RETOURNE RESULTAT
float gm(float *U, float *W, float a, int k, int N)
{
    unsigned int i;
    float sum_W, sum_R = 0;
    for (i = 0; i < N; i++)
    {
        sum_W += W[i];
        sum_R += powf(U[i] * W[i] - a, k);
    }
    return sum_R / sum_W;
}

float sum_R(float *U, float *W, float a, int k, int N)
{
}

//EXECUTE CALCUL VECTORIEL ET RETOURNE RESULTAT
float vect_gm(float *U, float *W, float a, int k, int N)
{

    float A[N_FLOAT] __attribute__((aligned(N_ALIGN)));
    float RV[N] __attribute__((aligned(N_ALIGN))); //résultat de la version vectorielle

    init_A(a, A);
    run_sse(U, W, A, RV, k, N);
    return sum(RV, N) / sum(W, N);
}

float gen_gm(void *thread_args)
{
    struct thread_args *arg = (struct thread_args *)thread_args;
    float *U = arg->U;
    float *W = arg->W;
    float a = arg->a;
    float k = arg->k;
    int n = arg->n;
    int mode = arg->mode;
    if (mode == 0)
    {
        arg->res.sum_of_weights = sum(W, n);
        arg->res.sum_of_weights = sum(W, n);
    }
    else
    {
        res = vect_gm(U, W, a, k, n);
    }
    pthread_exit(NULL);
}

float parallel_gm(float *U, float *W, float a, int k, int n, int mode, int nb_threads)
{
    pthread_t thread[nb_threads];
    pthread_attr_t attr;
    int error_code;
    void *status;
    long t;
    int size_of_partition = n / nb_threads;
    // Store result in struct table
    struct thread_args args[nb_threads];

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
        struct thread_args arg = args[t];
        arg.U = partition_of_u;
        arg.W = partition_of_v;
        arg.a = a;
        arg.n = size_of_partition;
        arg.mode = mode;
        printf("Main: creating thread %ld\n", t);
        error_code = pthread_create(&thread[t], &attr, gen_gm, (void *)&args);
        if (error_code)
        {
            printf("ERROR; return code from pthread_create() is %d\n", error_code);
            exit(-1);
        }
        pthread_attr_destroy(&attr);
        for (t = 0; t < nb_threads; t++)
        {
            error_code = pthread_join(thread[t], &status);
            if (error_code)
            {
                printf("ERROR: pthread_join() is %d\n", error_code);
                exit(-1);
            }
            printf("Join with thread %ld with status %ld\n", t, (long)status);
        }
        pthread_exit(NULL);
    }
}

int main(int argc, char const *argv[])
{
    // On lit les arguments passé avec l'appel de main
    char *end;
    int NUM_THREADS = strtol(argv[1], &end, 10);
    int N = strtol(argv[2], &end, 10);

    //Ici on déclare tous nos vecteurs en prenant soin de les aligner
    float U[N] __attribute__((aligned(N_ALIGN)));
    float W[N] __attribute__((aligned(N_ALIGN)));
    init(U, W, N);

    float rs, rv, rp;
    float a = 0;
    int k = 1;
    double t;

    // Calcul scalaire
    t = now();
    rs = gm(U, W, a, k, N);
    t = now() - t;
    printf("MP = %10.3g Temps du code scalaire : %f seconde(s)\n", rs, t);

    // Calcul vectoriel
    t = now();
    rv = vect_gm(U, W, a, k, N);
    t = now() - t;
    printf("MP = %10.3g Temps du code vectoriel : %f seconde(s)\n", rv, t);

    // Calcul parallèle
    t = now();
    rp = parallel_gm(U, W, 0, 1, N, 0, 1);
    t = now() - t;
    printf("MP = %10.3g Temps du code parallele : %f seconde(s)\n", rp, t);
}
