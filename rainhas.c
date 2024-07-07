#include "rainhas.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#define uint unsigned int
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

// Vetor para guardar solução da maior resposta
static uint maior_sol_sz;
static uint *maior_sol;

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

typedef struct backtraking_params {
    // n: Tamanho do lado do tabuleiro
    // cols: vetor booleano que indica
    //      cols[i] == 1: Tem rainha nessa coluna
    //      cols[i] == 0: c.c
    // diags2: vetor booleano que indica
    //      diags2[i] == 1: Tem rainha nessa diagona secundaria
    //      diags2[i] == 0: c.c
    //      Fórmula de indexação (dado linha i e coluna j):
    //          diag = i+j
    // diags1: vetor booleano que indica
    //      diags1[i] == 1: Tem rainha nessa diagona principal
    //      diags1[i] == 0: c.c
    //      Fórmula de indexação (dado linha i e coluna j):
    //          diag = i-j+n-1
    // mat: Matriz booleana que indica
    //      mat[i*n+j] == 1: casa (i,j) é proibida
    //      mat[i*n+j] == 0: c.c
    // row: Linha do tabuleiro atual da recursão
    // sol: Vetor de solução atual da recursão
    // sol_sz: tamanho da solução atual da recursão
    uint n;
    uint *cols;
    uint *diags2;
    uint *diags1;
    uint *mat;
    uint *sol;
    uint sol_sz;
    int *linhas_usadas;
} backtraking_params;

// Calcula quantas colunas estão livres para colocar uma rainha
static uint col_restantes(backtraking_params *bt, uint row){
    uint n = bt->n;
    uint count = 0;
    for (uint col = 0; col < n; col++){
        if (bt->mat[row*n+col] == 0 &&
            bt->cols[col] == 0 &&
            bt->diags2[col+row] == 0 &&
            bt->diags1[row-col+n-1] == 0){
            count++;
        }
    }
    return count;
}

// Pega próxima linha com menor número de casas livres para colocar rainhas
// OBS: Pode pegar linhas totalmente proibidas, está correto
static uint prox_linha(backtraking_params *bt){
    uint n = bt->n;
    uint menor = n+1;
    uint melhor = n+1;
    for (uint row = 0; row < n; row++){
        if (bt->linhas_usadas[row])
            continue;
        uint num_col_res = col_restantes(bt, row);
        if (num_col_res < menor){
            menor = num_col_res;
            melhor = row;
        }
    }
    return melhor;
}

// Calcula quantas linhas ainda já foram usadas
static uint quantidade_linhas_usadas(uint n, int *linhas_usadas){
    uint count = 0;
    for (uint i = 0; i < n; i++)
        if (linhas_usadas[i])
            count++;
    return count;
}

// Retorna se casa não está livre para colocar rainha
//     1: Se casa (row, col) não está livre
//     0: c.c
static inline int eh_casa_livre(backtraking_params *bt, uint row, uint col){
    // 1. Se for uma casa proibida
    // 2. Se tiver uma rainha na coluna
    // 3. Se tiver uma rainha na diagonal principal
    // 4. Se tiver uma rainha na diagonal secundaria
    uint n = bt->n;
    return  !(bt->mat[row*n+col] == 1     ||
            bt->cols[col] == 1          ||
            bt->diags2[col+row] == 1    ||
            bt->diags1[row-col+n-1] == 1);
}

// Marca que rainha foi colocada na casa (row, col)
static inline void coloca_rainha(backtraking_params *bt, uint row, uint col){
    uint n = bt->n;
    bt->sol[row] = col+1;
    bt->cols[col] = bt->diags2[row+col] = bt->diags1[row-col+n-1] = 1;
    bt->linhas_usadas[row] = 1;
    bt->sol_sz++;
}

// Desmarca que rainha foi colocada na casa (row, col)
static inline void retira_rainha(backtraking_params *bt, uint row, uint col){
    uint n = bt->n;
    bt->sol_sz--;
    bt->linhas_usadas[row] = 0;
    bt->cols[col] = bt->diags2[row+col] = bt->diags1[row-col+n-1] = 0;
    bt->sol[row] = 0;
}

static int rainhas_bt_recursivo(backtraking_params *bt){
    // verifica se eh a melhor solucao ate o momento (sol_sz se refere a solucao criada na recursao anterior, 
    // por isso a comparacao feita logo no inicio).
    uint n = bt->n;
    if (bt->sol_sz > maior_sol_sz){
        maior_sol_sz = bt->sol_sz;
        memcpy(maior_sol, bt->sol, n * sizeof(uint));
    }

    // verifica se a solucao atual eh a completa
    if (bt->sol_sz == n)
        return 1;
    // Se a solução atual mais a quantidade de linhas não usadas não passar
    // da maior solução então corta busca
    if (bt->sol_sz + n-quantidade_linhas_usadas(n, bt->linhas_usadas) <= maior_sol_sz)
        return 0;

    // Pega próxima linha para tentar colocar rainha
    uint row = prox_linha(bt);
    if (row == n+1)
        return 0;

    // Itera sobre as colunas da linha
    for (uint col = 0; col < n; col++){
        if (!eh_casa_livre(bt, row, col))
            continue;

        // Coloca rainha na coluna e diagonais atual
        coloca_rainha(bt, row, col);

        // Se achou solucao retorna 1
        if (rainhas_bt_recursivo(bt) == 1)
            return 1;

        // Se não faz o backtracking
        retira_rainha(bt, row, col);
    }

    // Se não conseguiu solução na linha marca como usada
    bt->linhas_usadas[row] = 1;
    int res = rainhas_bt_recursivo(bt);
    // Volta linha, pois está voltando do backtracking
    bt->linhas_usadas[row] = 0;
    return res;
}

unsigned int *rainhas_bt(unsigned int n, unsigned int k, casa *c, unsigned int *r){
    backtraking_params bt;

    bt.n      = n;
    bt.sol_sz = 0;

    // Aloca memória
    bt.cols          = calloc(n,   sizeof (uint));
    bt.mat           = calloc(n*n, sizeof (uint));
    bt.diags2        = calloc(2*n, sizeof (uint));
    bt.diags1        = calloc(2*n, sizeof (uint));
    bt.sol           = calloc(n,   sizeof (uint));
    bt.linhas_usadas = calloc(n,   sizeof (int));
    maior_sol        = calloc(n,   sizeof (uint));
    for (uint i = 0; i < k; i++)
        bt.mat[(c[i].linha-1)*n+c[i].coluna-1] = 1;

    // Chama função recursiva
    rainhas_bt_recursivo(&bt);
    memcpy(r, maior_sol, n * sizeof(uint));

    // Libera memória
    free(bt.cols);
    free(bt.mat);
    free(bt.diags2);
    free(bt.diags1);
    free(bt.sol);
    free(bt.linhas_usadas);
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
// TODO: Fazer rainhas_ci_recursivo funcionar para respostas menores que n (ou seja respostas onde tem rainhas faltando)
// TODO: Ver problema de quando resposta tem tamanho 1, pode ser que pegue uma casa proibida (evitar isso)
// TODO: altera para lista de adjacência - não insere vértices proibidos

// Para todo vértice em C verifica se é vizinho de v e decrementa o grau
// se grau < 0 então i não está no grafo
// TODO: Refatorar código (utilizar strucst e mais funções de auxilio)
// Estrutura de lista de adjacência
typedef struct lista_adjacencia{
    uint tam;   // tamanho da lista de vizinhos
    uint *vizinhos;     // lista de vizinhos
    int ref_count;   // pilha para backtracking
    uint grau;
} lista_adjacencia;
static void retira_vizinhos(uint v, lista_adjacencia *c){
    for (uint i = 0; i < c[v].tam ; i++){
        uint u = c[v].vizinhos[i];
        c[u].ref_count--;
        c[u].grau--;
    }
}

// Para todo vértice em C verifica se é vizinho de v e decrementa o grau
// se grau > 0 então i está no grafo
static void adiciona_vizinhos(uint v, lista_adjacencia *c){
    for (uint i = 0; i < c[v].tam; i++){
        uint u = c[v].vizinhos[i];
        c[u].ref_count++;
        c[u].grau++;
    }
}

// Retorna primeiro vertice que está no grafo
static int get_avalible_vert_primeiro(uint n, lista_adjacencia *c, uint *res){
    for (uint i = 0; i < n; i++){
        if (c[i].ref_count == 1){
            *res = i;
            return 1;
        }
    }
    return 0;
}
// Retorna vertice de menor grau que está no grafo
static int get_avalible_vert_menor_grau(uint n, lista_adjacencia *c, uint *res){
    // Nenhum vertice tem como ter grau maior que n (numero de vertices)
    uint menor_grau = n+1;
    for (uint i = 0; i < n; i++){
        if (c[i].ref_count == 1 && c[i].grau < menor_grau){
            menor_grau = c[i].grau;
            *res = i;
        }
    }
    return menor_grau != n+1;
}

static uint grau2(lista_adjacencia *c, uint v){
    uint res = 0;
    for (uint i = 0; i < c[v].tam; i++)
        res += c[c[v].vizinhos[i]].grau;
   return res; 
}

static uint grau3(lista_adjacencia *c, uint n, uint v){
    uint res = 0;
    // Vetor que diz se essa casa pode ou não ser considerada na contagem
    int *proibido = calloc(n, sizeof(int));
    // Preenche com os vizinhos de v
    for (uint i = 0; i < c[v].tam; i++)
        proibido[c[v].vizinhos[i]] = 1;
    for (uint i = 0; i < c[v].tam; i++){
        uint u = c[v].vizinhos[i];
        for (uint j = 0; j < c[u].tam; j++){
            if (!proibido[c[u].vizinhos[j]])
                res += c[c[u].vizinhos[j]].grau;
        }
    }
    free(proibido);
    return res; 
}
// Retorna vertice que tem menor grau e maior grau2 e menor grau3
//  grau2 é a soma dos graus do vizinhos de v
//  grau3 é a soma dos graus do vizinhos dos vizinhos de v
//      que não incluem os vizinhos de v
static int get_avalible_vert_melhor_deg_3(uint n, lista_adjacencia *c, uint *res){
    // Nenhum vertice tem como ter grau maior que n (numero de vertices)
    uint menor_grau = n+1;
    uint grau2_a, grau2_b;
    for (uint i = 0; i < n; i++){
        if (c[i].ref_count == 1){
            // Se for grau menor ou
            // se grau for igual escolhe o que tem mairo grau2
            if ((c[i].grau < menor_grau) ||
                (c[i].grau == menor_grau && (grau2_a=grau2(c, i)) > (grau2_b=grau2(c, *res)))||
                (grau2_a == grau2_b && grau3(c, n, i) < grau3(c, n,*res))){
                menor_grau = c[i].grau;
                *res = i;
            }
        }
    }
    return menor_grau != n+1;
}

// Retorna vertice que tem menor grau e maior grau2
//  grau2 é a soma dos graus do vizinhos de v
static int get_avalible_vert_melhor_deg_2(uint n, lista_adjacencia *c, uint *res){
    // Nenhum vertice tem como ter grau maior que n (numero de vertices)
    uint menor_grau = n+1;
    for (uint i = 0; i < n; i++){
        if (c[i].ref_count == 1){
            // Se for grau menor ou
            // se grau for igual escolhe o que tem mairo grau2
            if ((c[i].grau < menor_grau) ||
                (c[i].grau == menor_grau && grau2(c, i) > grau2(c, *res))){
                menor_grau = c[i].grau;
                *res = i;
            }
        }
    }
    return menor_grau != n+1;
}

static void libera_lista(lista_adjacencia* lista, uint t){
    for(uint i = 0; i < t; i++)
        if (lista[i].vizinhos != NULL)
            free(lista[i].vizinhos);
    free(lista);
}

inline static uint get_col(uint v_id, uint n){
    return v_id%n;
}
inline static uint get_row(uint v_id, uint n){
    return v_id/n;
}

static lista_adjacencia *adiciona_vizinhos_lista(lista_adjacencia* lista, uint i, uint j){
    lista[i].vizinhos[lista[i].tam++] = j;
    lista[j].vizinhos[lista[j].tam++] = i;

    return lista;
}

// função que calcula o número de linhas restantes através do grafo
static uint linhas_restantes(lista_adjacencia *c, uint n){
    uint *linha_ocupada = calloc(n, sizeof(uint));
    uint count = 0;
    for (uint i = 0; i < n*n; i++){
        if (c[i].ref_count == 1){
            uint row = get_row(i, n);
            if (!linha_ocupada[row]){
                linha_ocupada[row] = 1;
                count++;
            }
        }
    }

    free(linha_ocupada);
    return count;
}
// função que calcula o número de linhas restantes através do grafo
static uint colunas_restantes(lista_adjacencia *c, uint n){
    uint *coluna_ocupada = calloc(n, sizeof(uint));
    uint count = 0;
    for (uint i = 0; i < n*n; i++){
        if (c[i].ref_count == 1){
            uint col = get_col(i, n);
            if (!coluna_ocupada[col]){
                coluna_ocupada[col] = 1;
                count++;
            }
        }
    }

    free(coluna_ocupada);
    return count;
}

// Variavel para testar heuristicas, ALTERAR para testar novas
int (*get_avalible_vert)(uint, lista_adjacencia *, uint *) = get_avalible_vert_melhor_deg_2;
static int rainhas_ci_recursivo(uint n, uint t, uint I_sz, uint *I, lista_adjacencia *c){
    if (I_sz == n){
        memcpy(maior_sol, I, n*sizeof(uint));
        maior_sol_sz = I_sz;
        return 1;
    }
    // Se não tiver linhas ou colunas suficientes para passar da maior solução
    if (I_sz + MIN(colunas_restantes(c, n), linhas_restantes(c, n)) <= maior_sol_sz){
        return 0;
    }
    // Remove vertice v de C
    uint v;
    if (!get_avalible_vert(t, c, &v)){
        // Se I for maior solução até agora
        if (I_sz > maior_sol_sz){
            memcpy(maior_sol, I, n*sizeof(uint));
            maior_sol_sz = I_sz;
        }
        return 0;
    }
    // Marca vertice como retirado
    c[v].ref_count = 0;

    // Adiciona v em I
    I[I_sz++] = v;

    // Remove vizinhos de v de C
    retira_vizinhos(v, c);

    // Chama recursivamente
    if (rainhas_ci_recursivo(n, t, I_sz, I, c))
        return 1;

    // ------------------- Backtracking -------------------
    // Se não achou remove v de I
    I[I_sz--] = 0;

    // Adiciona vizinhos de v em C
    adiciona_vizinhos(v, c);

    // Chama recursivamente
    int res = rainhas_ci_recursivo(n, t, I_sz, I, c);
    // Volta o vertice retirado, pois na chamda anterior ele não está retirado
    c[v].ref_count = 1;
    return res;
    // ------------------- Backtracking -------------------
}

static lista_adjacencia *cria_lista_adjacecia(uint n, uint k, casa *c){
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
        lista[i].ref_count = 1;

        if (lista[i].vizinhos == NULL)
            lista[i].vizinhos = calloc(4*n, sizeof(uint));

        uint row_i = get_row(i, n); 
        uint col_i = get_col(i, n);
        for (uint j = i + 1; j < n*n; j++){
            uint row_j = get_row(j, n);
            uint col_j = get_col(j, n);
            if (mat[j])
                continue;
            
            if (lista[j].vizinhos == NULL)
                lista[j].vizinhos = calloc(4*n, sizeof(uint));

            // verifica se casa i é vizinha de casa j
            // 1. Verifica linha
            // 2. Verifica coluna
            // 3. Verifica diagonal secundaria
            // 4. Verifica diagonal principal
            if (row_i == row_j || col_i == col_j || row_i + col_i == row_j + col_j || row_i - col_i == row_j - col_j){
                adiciona_vizinhos_lista(lista, i, j);
            }
        }
        // Grau é quantos vizinhos tem
        lista[i].grau = lista[i].tam;
    }

    free(mat);
    return lista;
}

unsigned int *rainhas_ci(unsigned int n, unsigned int k, casa *c, unsigned int *r) {
    lista_adjacencia *lista = cria_lista_adjacecia(n, k, c);
    //imrpime a lista
    // for (uint i = 0; i < n*n; i++){
    //     printf("%d: ", i);
    //     for (uint j = 0; j < lista[i].tam; j++){
    //         printf("%d ", lista[i].vizinhos[j]);
    //     }
    //     printf("\n");
    // }
    /* for (uint i = 0; i < n*n; i++){ */
    /*     printf("%d: %d\n", i, lista[i].tam); */
    /* } */

    uint t = n*n;

    // Vetor de vértices independentes
    uint *I = calloc(n, sizeof(uint));
    uint I_sz = 0;
    maior_sol = calloc(n, (sizeof (uint)));
    maior_sol_sz = 0;

    if (!rainhas_ci_recursivo(n, t, I_sz, I, lista))
        memcpy(I, maior_sol, n*sizeof(uint));

    for (uint i = 0; i < n; i++){
#ifdef DEBUG
        printf("%d ", I[i]);
#endif
        if (I[i] != 0)
            r[get_row(I[i], n)] = get_col(I[i], n)+1;
    }
#ifdef DEBUG
    printf("\n");
#endif
    libera_lista(lista, t);
    free(I);
    free(maior_sol);
    return r;
}
