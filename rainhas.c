#include "rainhas.h"
#include <stdlib.h>
#include <string.h>
#define uint unsigned int

//------------------------------------------------------------------------------
// computa uma resposta para a instância (n,c) do problema das n rainhas 
// com casas proibidas usando backtracking
//
//    n é o tamanho (número de linhas/colunas) do tabuleiro
//
//    c é um vetor de k 'struct casa' indicando as casas proibidas
//
//    r é um vetor de n posições (já alocado) a ser preenchido com a resposta:
//      r[i] = j > 0 indica que a rainha da linha i+1 fica na coluna j;
//      r[i] = 0     indica que não há rainha nenhuma na linha i+1
//
// devolve r
static int rainhas_bt_(uint n, uint *cols, uint *diags1, uint *diags2, uint *mat, uint row, uint *r);
static int rainhas_bt_(uint n, uint *cols, uint *diags1, uint *diags2, uint *mat, uint row, uint *r){
    // Caso base
    if (row == n)
        return 1;

    for (uint col = 0; col < n; col++){
        // 1. Se for uma casa proibida
        // 2. Se tiver uma rainha na coluna
        // 3. Se tiver uma rainha na diagonal principal
        // 4. Se tiver uma rainha na diagonal secundaria
        if (mat[row*n + col] == 1   ||
            cols[col] == 1          ||
            diags1[col + row] == 1  ||
            diags2[row-col+n-1] == 1)
            continue;

        // Coloca rainha na coluna e diagonais atual
        r[row] = col+1;
        cols[col] = diags1[row+col] = diags2[row-col+n-1] = 1;

        // Se achou retorna
        if (rainhas_bt_(n, cols, diags1, diags2, mat, row+1, r) == 1)
            return 1;

        r[row] = 0;
        // Se não faz o backtaking
        cols[col] = diags1[row+col] = diags2[row-col+n-1] = 0;
    }
    // Não achou nenhum lugar volta tudo
    return 0;
}

unsigned int *rainhas_bt(unsigned int n, unsigned int k, casa *c, unsigned int *r) {
    uint *cols = calloc(n, (sizeof (uint)));
    memset(cols, 0, n*(sizeof (uint)));

    uint *mat  = calloc(n*n, (sizeof (uint)));
    memset(mat, 0, n*n*(sizeof (uint)));
    for (uint i = 0; i < k; i++)
        mat[(c[i].linha-1)*n + c[i].coluna-1] = 1;

    uint *diags1  = calloc(2*n, (sizeof (uint)));
    memset(diags1, 0, 2*n*(sizeof (uint)));

    uint *diags2  = calloc(2*n, (sizeof (uint)));
    memset(diags2, 0, 2*n*(sizeof (uint)));
    
#ifdef DEBUG
    for (unsigned int i=0; i<n; i++) {
        for (unsigned int j=0; j<n; j++) {
            if (mat[i*n + j])
                printf("*");
            else
                printf(".");
        }
        printf("\n");
    }
#endif

    rainhas_bt_(n, cols, diags1, diags2, mat, 0, r);
    free(cols);
    free(diags1);
    free(diags2);
    free(mat);

    return r;
}
//------------------------------------------------------------------------------
// computa uma resposta para a instância (n,c) do problema das n
// rainhas com casas proibidas usando a modelagem do problema como
// conjunto independente de um grafo
//
// n, c e r são como em rainhas_bt()

unsigned int *rainhas_ci(unsigned int n, unsigned int k, casa *c, unsigned int *r) {

    n = n;
    k = k;
    c = c;

    return r;
}
