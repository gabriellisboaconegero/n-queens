#include "rainhas.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#define uint unsigned int
#define MAX(a, b) ((a) > (b) ? (a) : (b))


// Estrutura de lista de adjacência
typedef struct lista_adjacencia{
    uint size;   // tamanho da lista de vizinhos
    uint *neighbors;     // lista de vizinhos
    int stack;   // pilha para backtracking
} lista_adjacencia;

// Tamanho da maior solução bt
static uint maior_sol_bt_sz;
// Vetor que guarda maior solução bt
static uint *maior_sol_bt;

// Tamanho da maior solução ci
static uint maior_sol_ci_sz;
// Vetor que guarda maior solução ci
static uint *maior_sol_ci;


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
static int rainhas_bt_(uint n, uint *cols, uint *diags2, uint *diags1, uint *mat, uint row, uint *r, uint sol_sz){
    // Caso base
    if (row == n){
        // Verificar se a solução até agora é a maior pois pode ser que seja
        // e quando voltar a chamada recursiva ela não vai ser salva
        // vide caso 5 com diagonais cortadas para entender
        if (sol_sz > maior_sol_bt_sz){
            maior_sol_bt_sz = sol_sz;
            memcpy(maior_sol_bt, r, n*sizeof(uint));
        }
        return sol_sz == n;
    }
    // Não desce se não tiver como melhorar
    // Nesse caso se não tiver mais linhas que somam
    // mais que a maior solução
    // OMG: A diferença que essa bomba faz, subiu o topo de uns 14 para 17
    if (sol_sz + n-row <= maior_sol_bt_sz)
        return 0;

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
        if (rainhas_bt_(n, cols, diags2, diags1, mat, row+1, r, sol_sz+1) == 1)
            return 1;

        r[row] = 0;
        // Se não faz o backtaking
        cols[col] = diags2[row+col] = diags1[row-col+n-1] = 0;
    }
    /* for (int i = 0; i < n; i++) */
    /*     printf("%d ", r[i]); */
    /* printf("\n"); */
    // Salva maior sequencia até agora
    if (sol_sz > maior_sol_bt_sz){
        maior_sol_bt_sz = sol_sz;
        memcpy(maior_sol_bt, r, n*sizeof(uint));
    }
    // Não conseguiu colocar nessa linha, tenta colocar na proxima
    // É nessario para casos de uma linha inteira coberta ou em que a recursão
    // para sem chegar no caso base (row == n)
    return rainhas_bt_(n, cols, diags2, diags1, mat, row+1, r, sol_sz);
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
    maior_sol_bt = calloc(n, (sizeof (uint)));

    if (!rainhas_bt_(n, cols, diags2, diags1, mat, 0, r, 0))
        memcpy(r, maior_sol_bt, n*sizeof(uint));

    free(cols);
    free(diags1);
    free(diags2);
    free(mat);
    free(maior_sol_bt);

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

// Para todo vértice em C verifica se é vizinho de v e decrementa o grau
// se grau < 0 então i não está no grafo
static uint retira_vizinhos(uint v, lista_adjacencia *c, uint c_uns_count){
    for (uint i = 0; i < c[v].size ; i++){
        c[c[v].neighbors[i]].stack--;
        if (c[i].stack == 1)
            c_uns_count--;
    }
    return c_uns_count;
}

// Para todo vértice em C verifica se é vizinho de v e decrementa o grau
// se grau > 0 então i está no grafo
static uint adiciona_vizinhos(uint v, lista_adjacencia *c, uint c_uns_count){
    for (uint i = 0; i < c[v].size; i++){
        c[c[v].neighbors[i]].stack++;
        if (c[i].stack == 1)
            c_uns_count++;
    }
    return c_uns_count;
}

// Retorna primeiro vertice que está no grafo
static int get_avalible_vert(uint n, lista_adjacencia *c, uint *res){
    for (uint i = 0; i < n; i++){
        if (c[i].stack == 1){
            *res = i;
            return 1;
        }
    }
    return 0;
}

static void libera_lista(lista_adjacencia* lista, uint t){
    for(uint i = 0; i < t; i++)
        if (lista[i].neighbors != NULL)
            free(lista[i].neighbors);
    free(lista);
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
static int rainhas_ci_(uint n, uint t, uint I_sz, uint *I, lista_adjacencia *c, uint c_uns_count){
    if (I_sz == n){
        memcpy(maior_sol_ci, I, n*sizeof(uint));
        maior_sol_ci_sz = I_sz;
        return 1;
    }
    if (I_sz + c_uns_count < maior_sol_ci_sz)
        return 0;
    // Remove vertice v de C
    uint v;
    if (!get_avalible_vert(t, c, &v)){
        // Se I for maior solução até agora
        if (I_sz > maior_sol_ci_sz){
            memcpy(maior_sol_ci, I, n*sizeof(uint));
            maior_sol_ci_sz = I_sz;
        }
        return 0;
    }
    // Marca vertice como retirado
    c[v].stack = 0;
    c_uns_count--;

    // Adiciona v em I
    I[I_sz++] = v;

    // Remove vizinhos de v de C
    c_uns_count = retira_vizinhos(v, c, c_uns_count);

    // Chama recursivamente
    if (rainhas_ci_(n, t, I_sz, I, c, c_uns_count))
        return 1;

    // ------------------- Backtracking -------------------
    // Se não achou remove v de I
    I[I_sz--] = 0;

    // Adiciona vizinhos de v em C
    c_uns_count = adiciona_vizinhos(v, c, c_uns_count);

    // Chama recursivamente
    int res = rainhas_ci_(n, t, I_sz, I, c, c_uns_count);
    // Volta o vertice retirado, pois na chamda anterior ele não está retirado
    c[v].stack = 1;
    return res;
    // ------------------- Backtracking -------------------
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
    maior_sol_ci = calloc(n, (sizeof (uint)));
    maior_sol_ci_sz = 0;

    if (!rainhas_ci_(n, t, I_sz, I, lista, c_uns_count))
        memcpy(I, maior_sol_ci, n*sizeof(uint));

    for (uint i = 0; i < n; i++){
        printf("%d ", I[i]);
        if (I[i] != 0)
            r[get_row(I[i], n)] = get_col(I[i], n)+1;
    }
    printf("\n");
    libera_lista(lista, t);
    free(I);
    free(maior_sol_ci);
    return r;
}
