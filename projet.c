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
    srand(time(NULL));
    unsigned int i;
    for (i = 0; i < N; i++)
    {
        U[i] = (float)rand() / RAND_MAX;
        W[i] = (float)rand() / RAND_MAX;
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

//CALCULE LA SOMME PONDERE ET MISE A LA PUISSANCE D'UN VECTEUR
double sum_diff(float *U, float *W, float a, int k, int N)
{
    unsigned int i;
    float sum_R = 0;
    for (i = 0; i < N; i++)
    {
        sum_R += powf(U[i] * W[i] - a, k);
    }
    return sum_R;
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
        mm_p = _mm_set1_ps(1.0f);  // Pour calculer la puissance k
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
    float sum_W = sum(W, N);
    float sum_R = sum_diff(U, W, a, k, N);
    return sum_R / sum_W;
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

//FONCTION APPELE PAR CHAQUE THREAD
void *gen_gm(void *thread_args)
{
    struct thread_args *arg = (struct thread_args *)thread_args;
    float *U = arg->U;
    float *W = arg->W;
    float a = arg->a;
    int k = arg->k;
    int n = arg->n;
    int mode = arg->mode;
    if (mode == 0)
    {
        //Execute le calcul scalairement
        arg->res.sum_of_weights = sum(W, n);
        arg->res.sum_of_weighted_diff = sum_diff(U, W, a, k, n);
    }
    else
    {
        //Execute le calcul vectoriel
        arg->res.sum_of_weights = sum(W, n);
        float A[N_FLOAT] __attribute__((aligned(N_ALIGN)));
        float RV[n] __attribute__((aligned(N_ALIGN))); //résultat de la version vectorielle
        init_A(a, A);
        run_sse(U, W, A, RV, k, n);
        arg->res.sum_of_weighted_diff = sum(RV, n);
    }
    pthread_exit(NULL);
}

//EXECUTE LE CALCUL DE FACON MULTITHREADE, SOIT SCALAIRE, SOIT VECTORIEL
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
    double sum_of_weights = 0;
    double sum_of_weighted_diff = 0;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for (t = 0; t < nb_threads; t++)
    {

        // Create partition pointers of U and W
        float *partition_of_u = U + size_of_partition * t;
        float *partition_of_v = W + size_of_partition * t;
        // Create struct for passing args to compute function
        struct thread_args *arg = &args[t];
        arg->U = partition_of_u;
        arg->W = partition_of_v;
        arg->a = a;
        arg->k = k;
        arg->n = size_of_partition;
        arg->mode = mode;
        error_code = pthread_create(&thread[t], &attr, gen_gm, (void *)arg);
        if (error_code)
        {
            printf("ERROR; return code from pthread_create() is %d\n", error_code);
            exit(-1);
        }
    }
    for (t = 0; t < nb_threads; t++)
    {
        error_code = pthread_join(thread[t], &status);
        if (error_code)
        {
            printf("ERROR: pthread_join() is %d\n", error_code);
            exit(-1);
        }
        struct thread_args arg = args[t];
        sum_of_weights += arg.res.sum_of_weights;
        sum_of_weighted_diff += arg.res.sum_of_weighted_diff;
    }
    pthread_attr_destroy(&attr);
    return sum_of_weighted_diff / sum_of_weights;
}

int main(int argc, char const *argv[])
{
    int NUM_THREADS = 4;
    int N = 100000;

    //Ici on déclare tous nos vecteurs en prenant soin de les aligner
    float U[N] __attribute__((aligned));
    float W[N] __attribute__((aligned));
    init(U, W, N);

    float rs, rv, rps, rpv = 0;
    // Variable pour le calcul de variance
    float a = sum(U, N) / N;
    int k = 2;
    double t, tbase, tperf;

    // Calcul scalaire
    t = now();
    rs = gm(U, W, a, k, N);
    tbase = now() - t;
    printf("Variance = %10.5g Temps du code scalaire : %f seconde(s)\n", rs, tbase);

    // Calcul vectoriel
    t = now();
    rv = vect_gm(U, W, a, k, N);
    t = now() - t;
    printf("Variance = %10.5g Temps du code vectoriel : %f seconde(s)\n", rv, t);

    // Calcul parallèle scalaire
    t = now();
    rps = parallel_gm(U, W, a, k, N, 0, NUM_THREADS);
    t = now() - t;
    printf("Variance = %10.5g Temps du code parallele scalaire : %f seconde(s)\n", rps, t);

    // Calcul parallèle vectorielle
    t = now();
    rpv = parallel_gm(U, W, a, k, N, 1, NUM_THREADS);
    tperf = now() - t;
    printf("Variance = %10.5g Temps du code parallele vectorielle : %f seconde(s)\n", rpv, tperf);

    printf("La version multithreade vectorielle est %f fois plus rapide que la version de base\n", tbase / tperf);
}
