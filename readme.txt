# BACKTRACKING
## Estruturas de Dados

struct backtraking_state {
    n: Tamanho do lado do tabuleiro
    cols: vetor booleano que indica
         cols[i] == 1: Tem rainha nessa coluna
         cols[i] == 0: c.c
    diags2: vetor booleano que indica
         diags2[i] == 1: Tem rainha nessa diagona secundaria
         diags2[i] == 0: c.c
         Fórmula de indexação (dado linha i e coluna j):
             diag = i+j
    diags1: vetor booleano que indica
         diags1[i] == 1: Tem rainha nessa diagona principal
         diags1[i] == 0: c.c
         Fórmula de indexação (dado linha i e coluna j):
             diag = i-j+n-1
    mat: Matriz booleana que indica
         mat[i*n+j] == 1: casa (i,j) é proibida
         mat[i*n+j] == 0: c.c
    row: Linha do tabuleiro atual da recursão
    sol: Vetor de solução atual da recursão
    sol_sz: tamanho da solução atual da recursão
}

## Funções principais
- prox_linha: Seleciona a próxima linha a ser explorada com base na quantidade de colunas disponíveis, 
  buscando sempre a linha com menos opções para otimizar a busca.
- rainhas_bt_recursivo: Função recursiva de backtracking que tenta posicionar rainhas no tabuleiro,
  respeitando as restrições e casas proibidas.
- rainhas_bt: Função principal que inicializa as estruturas necessárias e chama a função de backtracking.

## O algoritmo
O algoritmo de backtracking consiste em selecionar a linha do tabuleiro que possui menor quantidade de
casas disponíveis e inserir uma rainha, verificando se a solução atual é a de maior tamanho encontrada até
então e guardando tal informação. As podas utilizadas para melhorar o desempenho do código são: verificar se
a quantidade de linhas usadas juntamente com o número de linhas disponíveis poderiam ultrapassar a maior solução já
encontrada e verificar se a solução já é n. O controle de colunas, diagonais e linhas utilizadas por alguma rainha é
feito através de vetores auxiliares.

# CONJUNTO INDEPENDENTE

## Estruturas de Dados
typedef struct lista_adjacencia{
    tam: tamanho da lista de vizinhos
    vizinhos: vetor de vizinhos
    ref_count: pilha para backtracking
    grau: grau do vertice
} lista_adjacencia;
	
## Funções principais
- rainhas_ci: Função principal que inicializa as estruturas necessárias e chama a função de conjunto independente.
- rainhas_ci_recursivo: Função recursiva que tenta posicionar rainhas no tabuleiro,
  respeitando as restrições e casas proibidas.
- get_avalible_vert_melhor_deg_2: Função que retorna o vertice que possui menor grau e maior segundo grau.

O grafo utilizado para a implementação do algoritmo se baseia em considerar cada casa do tabuleiro como um vértice,
com arestas que incidem entre casas onde duas rainhas conseguiriam se atacar. A estrutura utilizada para representar
tal grafo é uma lista de adjacência, na qual as casas proibidas não possuem vizinhos. Dessa maneira, se busca o maior
conjunto independente possível do grafo, adicionando a cada recursão um vértice ao grupo independente e retirando o mesmo
e seus vizinhos do grafo. Para tornar o algoritmo mais eficiente, é selecionado o vértice com menor primeiro grau e maior
segundo grau (https://www.gcsu.edu/sites/files/page-assets/node-808/attachments/ballardmyer.pdf), visando manter o maior
número de vértices no grafo ao retirar os vizinhos do vértice escolhido. Além disso, é utilizada uma poda para verificar
se o número de colunas ou linhas ainda sem rainhas são o suficiente para atingir a maior solução.

A poda é feita verificando escolhendo o menor valor entre número de linhas livres e colunas livres (uma linha/coluna é
livre se existe alguma casa nela que se pode colocar uma rainha). Se o tamanho da solução atual mais o tamanho desse
valor é menor(ou igual) que o tamanho da maior solução então é feita a podagem, pois não existe forma de colocar mais
rainhas do que a maior solução.

O algoritmo de conjunto independente (chamado de CI) funciona utilizando o conceito de Refrence Counting
(Contagem de Referência). O campo `ref_count` na struct lista_adjacencia é a variável que armazena uma valor para saber
quantas vezes esse vértice já foi retirado do grafo (porém o valor é negativo e começa em 1 para fins de simplicidade).
Então se um vértice possui `ref_count` igual a 1 ele não foi retirado nenhuma vez do grafo então pertence ao grafo, se tiver
valor menor que 1 então ele já foi retirado -`ref_count`+1 vezes do grafo e não pertence ao grafo. Esse coceito é utilizado para
escolher um vértice, calcular grau2 e grau3, verificar se casa é valida, etc. Isso pois é usando `ref_count` que é determinado se
o grafo contém aquele vértice ou não.

Utilizar o `ref_count` remove a cópia constante dos conjuntos e reestruturação do grafo, pois a estrutura do grafo em si não é
alterada, apenas a forma que os vértices são vistos é.

# CASOS DE TESTE
Foram utilizados os seguintes casos de teste para garantir um bom desempenho do algoritmo para diferentes casas proibidas:
1. Diagonais: Proíbe as diagonais do tabuleiro. Segue uma tabela para esse teste:

		entrada	backtracking	CI
		10	141		4411
		15	1020		9333
		20	1482		20167
		25	4738		1164999
		30	1084		95507
		35	904		80401
		40	9318		1775138
		45	2539		2860773
		50	4227		-
		55	15992		66615
		60	21128		-
		65	14158		-
		70	4517		-
		75	122575		1025336
		OBS: O algoritmo de backtracking funcionou ainda para entradas 100-170 em tempo razoável (menos de um minuto) para esse teste.
		OBS: '-' significa que demorou mais de um minuto para rodar.

2. Linhas: Proíbe linhas inteiras do tabuleiro.
3. Colunas: Proíbe colunas inteiras.
4. Rand: Proíbe aleatoriamente 5*n casas do tabuleiro.
5. Janela: Proíbe uma linha e uma coluna do centro do tabuleiro, além das primeiras e últimas linhas e colunas.

Vale ressaltar o comportamento do teste 3. Onde os valores resultantes do backtracking são muito maiores que o
de conjunto independente para tamanho de tabuleiros maiores que 7. Isso acontece pois o algoritmo de CI consegue fazer a poda
de acordo com a quantidade de linhas e colunas que ainda tem casa livres para colocar, enquanto o backtracking faz apenas para
linhas livres.