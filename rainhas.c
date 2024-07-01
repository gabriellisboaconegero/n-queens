#include "rainhas.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#define uint unsigned int
#define MAX(a, b) ((a) > (b) ? (a) : (b))


// Estrutura de lista de adjacência
typedef struct Vertice {
	uint grau;	// Grau do vértice
	uint *vizinhos;	// Lista de vizinhos
} Vertice;

typedef struct Grafo {
	uint size;			// Tamanho do grafo (em vertices)
	Vertice *vertices;	// Lista de vertices
} Grafo;

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

// Tamanho da maior solução
static uint maior_sol_sz;
// Vetor que guarda maior solução
static uint *maior_sol;

// n: Tamanho do lado do tabuleiro
// cols: vetor booleano que indica
//      cols[i] == 1: Tem rainha nessa coluna
//      cols[i] == 0: c.c
// diags2: vetor booleano que indica
//      diags2[i] == 1: Tem rainha nessa diagona secundaria
//      diags2[i] == 0: c.c
// diags1: vetor booleano que indica
//      diags1[i] == 1: Tem rainha nessa diagona principal
//      diags1[i] == 0: c.c
// mat: Matriz booleana que indica
//      mat[i][j] == 1: casa (i,j) é proibida
//      mat[i][j] == 0: c.c
// row: Linha do tabuleiro atual da recursão
// r: Vetor de solução atual da recusão
// sol_sz: tamanho da solução atual da recursão
// TODO: Implementar otimização de achar primeira linha livre
//     (evita ficar verificando a linha inteira até perceber que ela não tem espaço livre)
//     usar vetor rows_free_count que indica
//          rows_free_count[i] == k: Quer dizer que a linha i tem k casas não livres (proibida ou atacada por rainha)
// TODO: Refatorar código (utilizar strucst e mais funções de auxilio)

static uint maior_sol_bt_sz = 0;
static uint *maior_sol_bt;


static int countRemainingValues(uint n, uint row, uint *cols, uint *diags2, uint *diags1, uint *mat) {
    int count = 0;
    for (uint col = 0; col < n; col++) {
        if (mat[row*n + col] == 0 &&
            cols[col] == 0 &&
            diags2[col + row] == 0 &&
            diags1[row - col + n - 1] == 0) {
            count++;
        }
    }
    return count;
}

static int findNextRow(uint n, uint *cols, uint *diags2, uint *diags1, uint *mat, int *rows_used) {
    uint min_remaining = n + 1;
    int best_row = -1;
    for (uint row = 0; row < n; row++) {
        if (rows_used[row])
            continue;
        uint remaining = countRemainingValues(n, row, cols, diags2, diags1, mat);
        if (remaining < min_remaining) {
            min_remaining = remaining;
            best_row = row;
        }
    }
    return best_row;
}

static int rainhas_bt_(uint n, uint *cols, uint *diags2, uint *diags1, uint *mat, uint *r, uint sol_sz, int *rows_used) {
    // Caso base: Se sol_sz é igual a n, encontramos uma solução completa
    if (sol_sz > maior_sol_bt_sz) {
        maior_sol_bt_sz = sol_sz;
        memcpy(maior_sol_bt, r, n * sizeof(uint));
    }

    if (sol_sz == n)
        return 1;

    int row = findNextRow(n, cols, diags2, diags1, mat, rows_used);
    if (row == -1)
        return 0;

    for (uint col = 0; col < n; col++) {
        if (mat[row*n + col] == 1 || cols[col] == 1 || diags2[col + row] == 1 || diags1[row - col + n - 1] == 1)
            continue;

        r[row] = col + 1;
        cols[col] = diags2[row + col] = diags1[row - col + n - 1] = 1;
        rows_used[row] = 1;

        if (rainhas_bt_(n, cols, diags2, diags1, mat, r, sol_sz + 1, rows_used) == 1)
            return 1;

        r[row] = 0;
        cols[col] = diags2[row + col] = diags1[row - col + n - 1] = 0;
        rows_used[row] = 0;
    }

    return 0;
}

unsigned int *rainhas_bt(unsigned int n, unsigned int k, casa *c, unsigned int *r) {
    uint *cols = calloc(n, sizeof(uint));
    uint *mat = calloc(n * n, sizeof(uint));
    for (uint i = 0; i < k; i++)
        mat[(c[i].linha - 1) * n + c[i].coluna - 1] = 1;

    uint *diags2 = calloc(2 * n, sizeof(uint));
    uint *diags1 = calloc(2 * n, sizeof(uint));
    maior_sol_bt = calloc(n, sizeof(uint));
    int *rows_used = calloc(n, sizeof(int));

    rainhas_bt_(n, cols, diags2, diags1, mat, r, 0, rows_used);
    memcpy(r, maior_sol_bt, n * sizeof(uint));

    free(cols);
    free(diags1);
    free(diags2);
    free(mat);
    free(maior_sol_bt);
    free(rows_used);

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
// TODO: Fazer rainhas_ci_ funcionar para respostas menores que n (ou seja respostas onde tem rainhas faltando)
// TODO: Ver problema de quando resposta tem tamanho 1, pode ser que pegue uma casa proibida (evitar isso)
// TODO: altera para lista de adjacência - não insere vértices proibidos

/*static int retira_random_avalible_vert(uint n, Grafo *c, uint *res){
    for (uint i = 0; i < n; i++){
        if (c[i].stack == 1){
            *res = i;
            return 1;
        }
    }
    return 0;
}

static void libera_grafo(Grafo* g){
    for(uint i = 0; i < g->size; i++)
        if (g->vertices[i].vizinhos != NULL)
            free(g->vertices[i].vizinhos);
    free(g->vertices);
    free(g);
}

inline Vertice *get_vertice(Grafo *g, uint i){
	return &(g->vertices[i]);
}

static copia_grafo(Grafo *g){
    Grafo *res = calloc(g->size, sizeof(Grafo));
    for (uint i = 0; i < g->size; i++){
        Vertice *res_u = get_vertice(res, i);
        Vertice *g_u = get_vertice(g, i);
        res_u->vizinhos = calloc(4*g->size, sizeof(uint));
        for (uint j = 0; j < g_u->grau; j++)
            res_u->vizinhos[j] = g->vertices[i].vizinhos[j];
    }
    res->size = g->size;
    return res;
}

static uint *remove_vertice(Grafo *g, uint v_){
    Vertice *v = get_vertice(g, v_);
    g->size--;
    for (uint i = 0; i < v->grau; i++){
        Vertice *u = get_vertice(g, v->vizinhos[i]);
        for (uint j = 0; j < u->grau-1; j++)
            if (u->vizinhos[j] >= v){
                u->vizinhos[j] = u->vizinhos[u->grau+1];
                break;
            }
        u->grau--;
    }

    Vertice *v_vizinhanca = calloc(v->grau, sizeof(uint));
    v_vizinhanca->grau = v->grau;
    v_vizinhanca->vizinhos = v->vizinhos;

    v->vizinhos = NULL;
    v->grau = 0;

    return v_vizinhanca;
}

static void remove_vizinhanca(Grafo *g, Vertice *v){
    for (int i = 0; i < v->grau; i++){
        remove_vertice(g, v->vizinhos[i]);
    }

    free(v->vizinhos);
}

static int rainhas_ci_(uint n, uint I_sz, uint *I, Grafo *c){
    if (I_sz == n)
        return 1;
    if (I_sz + c->size < n)
        return 0;
    // Remove vertice aleatório v de C
    uint v;
    if (!retira_random_avalible_vert(c, &v))
        return 0;
    // Marca vertice como retirado
    Grafo *c_copy = copia_grafo(c);
    Vertice *v_vizinhanca = remove_vertice(c_copy, v);

    // Chama recursivamente
    int res = rainhas_ci_(n, I_sz, I, c_copy);
    if (res){
        libera_grafo(c_copy);
        return 1;
    }

    // Adiciona v em I
    I[I_sz++] = v;
    remove_vizinhanca(c_copy, v_vizinhanca);

    // ------------------- Backtracking -------------------
    // Chama recursivamente
    res = rainhas_ci_(n, I_sz, I, c_copy);
    libera_grafos(c_copy);
    return res;
    // ------------------- Backtracking -------------------
}

inline static uint get_col(uint v_id, uint n){
    return v_id%n;
}
inline static uint get_row(uint v_id, uint n){
    return v_id/n;
}

static lista_adjacencia *adiciona_vizinhos_lista(lista_adjacencia* lista, uint i, uint j){
    lista[i].neighbors[lista[i].size++] = j;
    lista[j].neighbors[lista[j].size++] = i;

    return lista;
}

static lista_adjacencia *cria_lista_adjacecia(uint n, uint k, casa *c, uint *c_uns_count){
    uint size_list = n*n;

    // Matriz para dizer se casa é proibida ou não
    uint *mat = calloc(size_list, sizeof(uint));
    for (uint i = 0; i < k; i++)
        // c[i].linha-1 pois é indexado começando em 1
        mat[(c[i].linha-1)*n + c[i].coluna-1] = 1;

    
    struct lista_adjacencia *lista = calloc(size_list, sizeof(struct lista_adjacencia));
    
    // preenche a lista de adjacencia, ignorando casas proibidas
    for (uint i = 0; i < size_list; i++){
        if (mat[i]) // Se casa é proibida, ignora
            continue;
        lista[i].stack = 1;
        *c_uns_count += 1;

        if (lista[i].neighbors == NULL)
            lista[i].neighbors = calloc(4*n, sizeof(uint));

        uint row_i = get_row(i, n); 
        uint col_i = get_col(i, n);
        for (uint j = i + 1; j < n*n; j++){
            uint row_j = get_row(j, n);
            uint col_j = get_col(j, n);
            if (mat[j])
                continue;
            
            if (lista[j].neighbors == NULL)
                lista[j].neighbors = calloc(4*n, sizeof(uint));

            // verifica se casa i é vizinha de casa j
            // 1. Verifica linha
            // 2. Verifica coluna
            // 3. Verifica diagonal secundaria
            // 4. Verifica diagonal principal
            if (row_i == row_j || col_i == col_j || row_i + col_i == row_j + col_j || row_i - col_i == row_j - col_j){
                adiciona_vizinhos_lista(lista, i, j);
            }
        }
        
    }

    free(mat);
    return lista;
}

unsigned int *rainhas_ci(unsigned int n, unsigned int k, casa *c, unsigned int *r) {
    //int *mat_adj = cria_matriz_adjacencia(n, k, c);
    uint c_uns_count = 0;
    lista_adjacencia *lista = cria_lista_adjacecia(n, k, c, &c_uns_count);
    //imrpime a lista
    // for (uint i = 0; i < n*n; i++){
    //     printf("%d: ", i);
    //     for (uint j = 0; j < lista[i].size; j++){
    //         printf("%d ", lista[i].neighbors[j]);
    //     }
    //     printf("\n");
    // }

    uint t = n*n;

    // Vetor de vértices independentes
    uint *I = calloc(n, sizeof(uint));
    uint I_sz = 0;

    rainhas_ci_(n, t, I_sz, I, lista, c_uns_count);
    for (uint i = 0; i < n; i++){
        printf("%d ", I[i]);
        r[get_row(I[i], n)] = get_col(I[i], n)+1;
    }
    printf("\n");
    libera_lista(lista, t);
    free(I);
    return r;
}
*/