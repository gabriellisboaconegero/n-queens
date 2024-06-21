#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "rainhas.h"

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

//------------------------------------------------------------------------------
int main (int argc, char **argv) {

    unsigned int n = atoi(argv[1]);
    unsigned int *resposta = malloc(n*sizeof(unsigned int));

    unsigned int k = 2 * n;
    casa *proibido = malloc(k*sizeof(casa));

    proibe_diagonais(n, proibido);

    printf("backtracking: ");
    long int tempo_bt;
    CRONOMETRA(rainhas_bt(n, k, proibido, resposta), tempo_bt);
    printf("%ld\n", tempo_bt);
    mostra_resposta(n, resposta, proibido, k);

//     printf("grafo: ");
//     long int tempo_ci;
//     CRONOMETRA(rainhas_ci(n, k, proibido, resposta), tempo_ci);
//     printf("%ld\n", tempo_ci);
//     mostra_resposta(n, resposta, proibido, k);

//     printf("%.2f\n", (double)tempo_ci/(double)tempo_bt);

    free(proibido);
    free(resposta);

    return 0;
}
