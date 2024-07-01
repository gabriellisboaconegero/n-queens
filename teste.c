#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "rainhas.h"
#define uint unsigned int

static clock_t inicio, fim;

//------------------------------------------------------------------------------
#define CRONOMETRA(call,t) inicio = clock(); (call); fim = clock(); t=fim-inicio

//------------------------------------------------------------------------------
static void mostra_resposta(unsigned int n, unsigned int *r, casa *c, unsigned int l) {
    for (unsigned int i=0; i<n; i++) {
        printf("%u, ", r[i]);
    }
#ifdef DEBUG
    printf("\n");
    for (unsigned int i=0; i<n; i++) {
        for (unsigned int j=0; j<2*n; j++) {
            for (unsigned int k=0; k < l; k++){
                if (i == (c[k].linha-1) && j/2 == (c[k].coluna-1)){
                    printf("\x1b[1;31m█\x1b[m");
                    goto defer;              
                }

            }
            if (j/2 == r[i]-1)
                printf("\x1b[1;33m█\x1b[m");
            else if ((j/2 + i) % 2)
                printf("\x1b[1;30m█\x1b[m");
            else
                printf("█");
defer:
        }
        printf("\n");
    }
#endif

    printf("\n");
}
//------------------------------------------------------------------------------
// preenche proibido[0..2n-1] com as posições das diagonais do tabuleiro n x n
//
// devolve &(proibido[2n]): o endereço a partir do qual proibir novas posições

static casa *proibe_diagonais(unsigned int n, casa *proibido) {
    // proíbe todas as casas nas diagonais
    for(unsigned int i = 0, p = 1; i < 2 * n; i+=2, p++) {
        // diagonal "principal"
        proibido[i].linha = proibido[i].coluna = p;

        // "outra" diagonal
        proibido[i+1].linha = p;
        proibido[i+1].coluna = n - p + 1;
    }



    return proibido + 2*n;
}

static casa *proibe_linha(unsigned int n, casa *proibido, unsigned int linha) {
    for(unsigned int i = 0; i < n; i++) {
        proibido[i].linha = linha;
        proibido[i].coluna = i+1;
    }

    return proibido + n;
}
static casa *proibe_coluna(unsigned int n, casa *proibido, unsigned int coluna) {
    for(unsigned int i = 0; i < n; i++) {
        proibido[i].linha = i+1;
        proibido[i].coluna = coluna;
    }

    return proibido + n;
}

static casa *proibe_random(unsigned int n, casa *proibido, int max) {
    for(int i = 0; i < max; i++) {
        proibido[i].linha = rand() % n + 1;
        proibido[i].coluna = rand() % n + 1;
    }

    return proibido + max;
}

//------------------------------------------------------------------------------
/* int main1 (int argc, char **argv) { */
/*     int adj_mat[] = { */
/*         1, 0, 0, 0, 1, 1, 1, 1, */
/*         0, 1, 0, 0, 0, 0, 1, 0, */
/*         0, 0, 1, 0, 0, 1, 0, 0, */
/*         0, 0, 0, 1, 0, 1, 1, 0, */
/*         1, 0, 0, 0, 1, 1, 0, 0, */
/*         1, 0, 1, 1, 1, 1, 0, 0, */
/*         1, 1, 0, 1, 0, 0, 1, 0, */
/*         1, 0, 0, 0, 0, 0, 0, 1, */
/*     }; */
/*     uint adj_mat_sz = 8; */
/*     int *C = calloc(adj_mat_sz, sizeof(int)); */
/*     for (uint i = 0; i < adj_mat_sz; i++) */
/*         C[i] = 1; */
/*     uint c_uns_count = adj_mat_sz; */
/*     uint *I = calloc(adj_mat_sz, sizeof(uint)); */
/*     uint I_sz = 0; */
/*     rainhas_ci_(5, adj_mat, adj_mat_sz, I_sz, I, C, c_uns_count); */
/*     for (uint i = 0; i < adj_mat_sz; i++){ */
/*         printf("%d ", I[i]); */
/*     } */
/*     printf("\n"); */
/* } */
int main (int argc, char **argv) {

    if (argc < 2){
        fprintf(stderr, "[ERRO]: %s <tam>\n", argv[0]);
        exit(1);
    }
    unsigned int n = atoi(argv[1]);
    unsigned int *resposta = malloc(n*sizeof(unsigned int));

    casa *proibido = malloc(n*n*3*sizeof(casa));

    unsigned int k = 5*n;   //deve ser alterado conforme o tipo de proibição
    // proibe_diagonais(n, proibido); 
    // k = 2*n;
    // proibe_coluna(n, 
    // proibe_coluna(n, 
    // proibe_coluna(n, 
    // proibe_linha(n, 
    // proibe_linha(n, 
    // proibe_linha(n, proibido, 
    // 1), n), n/2+1), 1), n), n/2+1); 
    // k = 6*n;
    srand(1234);
    proibe_random(n, proibido, k);

    printf("backtracking: ");
    long int tempo_bt;
    CRONOMETRA(rainhas_bt(n, k, proibido, resposta), tempo_bt);
    printf("%ld\n", tempo_bt);
    mostra_resposta(n, resposta, proibido, k);

   /* memset(resposta, 0, n*sizeof(unsigned int));

    printf("grafo: ");
    long int tempo_ci;
    CRONOMETRA(rainhas_ci(n, k, proibido, resposta), tempo_ci);
    printf("%ld\n", tempo_ci);
    mostra_resposta(n, resposta, proibido, k);

    /* printf("%.2f\n", (double)tempo_ci/(double)tempo_bt); */

    free(proibido);
    free(resposta);

    return 0;
}
