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
#include <xmmintrin.h>
#include <pmmintrin.h>

double now(){
    // Retourne l'heure actuelle en secondes
    struct timeval t; double f_t;
    gettimeofday(&t, NULL);
    f_t = t.tv_usec; f_t = f_t/1000000.0; f_t +=t.tv_sec;
    return f_t;
}

//INITIALISATION DES VECTEURS U ET W
void init(float *U, float *W, unsigned int N){
    unsigned int i;
    for(i=0; i<N; i++){
        U[i] = (float)rand () / RAND_MAX / N; 
        W[i] = (float)rand () / RAND_MAX / N;  
    }
}

//CALCULE LA SOMME D'UN VECTEUR DE TAILLE N
float sum(float *v, int N){
    unsigned int i;
    float s = 0;
    for(i=0; i<N; i++) s += v[i];
    return s;
}

//REMPLIT LE VECTEUR A PAR LA VALEUR a
void init_A(float a, float *A){
    unsigned int i;
    for(i=0; i<N_FLOAT; i++)
        A[i] = a;
}

//VERSION VECTORIELLE
void run_sse(float *U, float *W, float *A, float *RV, int k, int N){
    unsigned int i;
    unsigned int j;
    __m128 mm_U, mm_W, mm_A, mm_t, mm_p; //On déclare cinq registres vectoriels
    mm_A = _mm_load_ps(&A[0]);
    mm_p = _mm_set1_ps(1.0f); // Pour calculer la puissance k

    for(i=0; i<N; i+=4){
        mm_U = _mm_load_ps(&U[i]);  //On charge (U[4i], U[4i+1], U[4i+2], U[4i+3])
        mm_W = _mm_load_ps(&W[i]);  //On charge (W[4i], W[4i+1], W[4i+2], W[4i+3])
        mm_t = _mm_mul_ps(mm_U, mm_W);
        mm_t = _mm_sub_ps(mm_t, mm_A);
        for(j=0; j<k; j++){ //On calcule la puissance k
            mm_p = _mm_mul_ps(mm_p, mm_t);
        }
        _mm_store_ps(&RV[i], mm_p);  //On stocke  (p[4i], p[4i+1], p[4i+2], p[4i+3])
      }
   }

//EXECUTE CALCUL SCALAIRE ET RETOURNE RESULTAT
float gm(float *U, float *W, float a, int k, int N)
{
    unsigned int i;
    float sum_W, sum_R = 0;
    for (i=0; i<N; i++)
    {
        sum_W += W[i];
        sum_R += powf(W[i] * U[i] - a, k);
    }
    return sum_R / sum_W;
}

//EXECUTE CALCUL VECTORIEL ET RETOURNE RESULTAT
float vect_gm(float *U, float *W, float a, int k, int N)
{
    float A[N_FLOAT] __attribute__((aligned(N_ALIGN)));
    float RV[N] __attribute__((aligned(N_ALIGN))); //résultat de la version vectorielle

    init_A(a, A);
    run_sse(U, W, A, RV, k, N);
    return sum(RV, N)/sum(W, N);
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

    float rs, rv;
    float a=0;
    int k=1;
    double t;
    
    // Calcul scalaire
    t = now();
    rs = gm(U, W, a, k, N);
    t = now()-t;
    printf("S = %10.3g Temps du code scalaire : %f seconde(s)\n", rs, t);

    // Calcul vectoriel
    t = now();
    rv = vect_gm(U, W, a, k, N);
    t = now()-t;
    printf("S = %10.3g Temps du code vectoriel : %f seconde(s)\n", rv, t);

    return 0;
}
