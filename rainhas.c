#include "rainhas.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#define uint unsigned int
#define MAX(a, b) ((a) > (b) ? (a) : (b))

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
static uint maior_sol_sz;
static uint *maior_sol;
static int rainhas_bt_(uint n, uint *cols, uint *diags2, uint *diags1, uint *mat, uint row, uint *r);
static int rainhas_bt_(uint n, uint *cols, uint *diags2, uint *diags1, uint *mat, uint row, uint *r){
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
            diags2[col + row] == 1  ||
            diags1[row-col+n-1] == 1)
            continue;

        // Coloca rainha na coluna e diagonais atual
        r[row] = col+1;
        cols[col] = diags2[row+col] = diags1[row-col+n-1] = 1;

        // Se achou retorna
        if (rainhas_bt_(n, cols, diags2, diags1, mat, row+1, r) == 1)
            return 1;

        r[row] = 0;
        // Se não faz o backtaking
        cols[col] = diags2[row+col] = diags1[row-col+n-1] = 0;
    }
    // Salva maior sequencia até agora
    if (row+1 > maior_sol_sz){
	maior_sol_sz = row;
	memcpy(maior_sol, r, n*sizeof(uint));
    }
    // Não achou nenhum lugar volta tudo
    return 0;
}

unsigned int *rainhas_bt(unsigned int n, unsigned int k, casa *c, unsigned int *r) {
    // Vetor de colunas:
    //      cols[i] == 1: se tem rainha na coluna i
    //      cols[i] == 0: c.c
    uint *cols = calloc(n, (sizeof (uint)));

    // Matriz de casas proibidas (indexada por i*n + j):
    //      mat[i*n+j] == 1: se casa é proibida
    //      mat[i*n+j] == 0: c.c
    uint *mat  = calloc(n*n, (sizeof (uint)));
    for (uint i = 0; i < k; i++)
        mat[(c[i].linha-1)*n + c[i].coluna-1] = 1;

    // Vetor de diagonal secundaria:
    //      diags2[i] == 1: Se tem rainha na diagonal secundaria
    //      diags2[i] == 0: c.c
    // Fórmula de indexação (dado linha i e coluna j):
    //      diag = i + j
    uint *diags2  = calloc(2*n, (sizeof (uint)));

    // Vetor de diagonal principal:
    //      diags1[i] == 1: Se tem rainha na diagonal principal
    //      diags1[i] == 0: c.c
    // Fórmula de indexação (dado linha i e coluna j):
    //      diag = i - j + n - 1
    uint *diags1  = calloc(2*n, (sizeof (uint)));

    // Vetor de maior sequencia, guarda a maior sequencia
    // em caso de não ter solução o problema
    maior_sol = calloc(n, (sizeof (uint)));
    
    if (!rainhas_bt_(n, cols, diags2, diags1, mat, 0, r))
        memcpy(r, maior_sol, n*sizeof(uint));

    free(cols);
    free(diags1);
    free(diags2);
    free(mat);
    free(maior_sol);

    return r;
}
//------------------------------------------------------------------------------
// computa uma resposta para a instância (n,c) do problema das n
// rainhas com casas proibidas usando a modelagem do problema como
// conjunto independente de um grafo
//
// n, c e r são como em rainhas_bt()

// TODO: Mudar c_uns_count para int
// TODO: Implementar função para retirar vértice v de C

// Para todo vértice em C verifica se é vizinho de v e decrementa o grau
// se grau < 0 então i não está no grafo
uint retira_vizinhos(uint t, int *mat_adj, uint v, int *c, uint c_uns_count){
    for (uint i = 0; i < t ; i++){
	if (v != i && mat_adj[v*t + i] == 1){
	    c[i]--;
	    if (c[i] == 0)
		c_uns_count--;
	}
    }
    return c_uns_count;
}

// Para todo vértice em C verifica se é vizinho de v e decrementa o grau
// se grau > 0 então i está no grafo
uint adiciona_vizinhos(uint t, int *mat_adj, uint v, int *c, uint c_uns_count){
    for (uint i = 0; i < t; i++){
	if (v != i && mat_adj[v*t + i] == 1){
	    c[i]++;
	    if (c[i] == 1)
		c_uns_count++;
	}
    }
    return c_uns_count;
}

static int rainhas_ci_(uint n, int *mat_adj, uint t, uint I_sz, uint *I, int *c, uint c_uns_count){
    if (I_sz == n)
        return 1;
    if (I_sz + c_uns_count < n)
        return 0;
	// Remove vertice aleatório v de C
	c[v] = 0;
	c_uns_count--;
	// Adiciona v em I
	I[I_sz++] = v;
	// Remove vizinhos de v de C
	c_uns_count = retira_vizinhos(t, mat_adj, v, c, c_uns_count);
	// Chama recursivamente
	if (rainhas_ci_(n, mat_adj, t, I_sz, I, c, c_uns_count))
		return 1;
	
	// ------------------- Backtracking -------------------
	// Se não achou remove v de I
	I_sz--;
	// Adiciona vizinhos de v em C
	c_uns_count = adiciona_vizinhos(t, mat_adj, v, c, c_uns_count);
	// Chama recursivamente
	return rainhas_ci_(n, mat_adj, t, I_sz, I, c, c_uns_count);
	// ------------------- Backtracking -------------------
}

unsigned int *rainhas_ci(unsigned int n, unsigned int k, casa *c, unsigned int *r) {
    int *mat_adj = calloc(n*n*n*n, sizeof(int));
    uint t = n*n;

    for (uint i = 0; i < t; i++){
        for (uint j = i; j < t; j++){
            uint row_i = i/n;
            uint col_i = i%n;
            uint row_j = j/n;
            uint col_j = j%n;

            // Verifica se casa i é vizinha de casa j
            // 1. Verifica linha
            // 2. Verifica coluna
            // 3. Verifica diagonal secundaria
            // 4. Verifica diagonal principal
            mat_adj[i*t + j] =  (row_i == row_j) ||
                                (col_i == col_j) ||
                                (row_i + col_i == row_j + col_j) ||
                                (row_i - col_j == row_j - col_j);
        }
    }
	// Vetor de situação dos vértices
	// C[i] < 1: se vértice i não está no grafo
	// C[i] == 1: se vértice i está no grafo
	int *C = calloc(t, sizeof(int));
	memset(C, 1, t*sizeof(int));
	// Retira as casa proibidas do grafo
	for (uint i = 0; i < k; i++)
		C[(c[i].linha-1)*n + c[i].coluna-1] = 0;
	// Após retirar os vértices proibidos, a quantidade de vértices no grafo é t - k
	uint c_uns_count = t - k;
	
	// Vetor de vértices independentes
	uint *I = calloc(n, sizeof(uint));
	uint I_sz = 0;

	rainhas_ci_(n, mat_adj, t, I_sz, I, C, c_uns_count);
	free(C);
	free(mat_adj);
	free(I);
    return r;
}
