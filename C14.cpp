#include<stdio.h>
#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<windows.h>
#include<omp.h>

#define LINK 15     //tamanho da cadeia de ponteiros de cada nó da árvore
#define VETMAX 14   //tamanho do vetor do vetor de subciclos auxiliar

//Definição de cada ramo no vetor linear
struct Node{
    int barra_1, barra_2, chave, id; //declaração dos vértices de cada ramo, e configuração da chave - ligada/desligada/sempre ligada
    float resist, reat; //valores da impedância do ramo, real e imaginário
    float I_real, I_imag; //valores de corrente, real e imaginário, do ramo
};

//Definição do elemento de Potência de Carga e Tensão para cada barra
struct Tensao_Carga{
    float pot_at, pot_reat; //valores de potência ativa e reativa
    float V_real, V_imag; //valores de tensão real e imaginário
};

//Estrutura dinâmica que armazena as barras ativas para construir a árvore do sistema elétrico
struct V_Node{
        int v, cont_p, id;// índice da barra, contador de conexões da barra (filhos)
        struct V_Node* p[LINK];  //vetor de ponteiro para os filhos
    };
typedef V_Node* VN;

//Estrutura dinâmica auxiliar para guardar as barras do sistema para futura comparação e construção da árvore
struct V_aux{
    int v;
    struct V_aux* prox;
};
typedef V_aux* v_aux;

void PRIM(struct Node *N, float FC_ini[], int ramos, int subst, int barras);
void GRASP(struct Node *N, struct Node *N_best, float FC_ini[], int ramos, int subst, int barras);
void VTree(struct Node *N, struct V_aux *vaux, int ramos, int subst, int barras, struct V_Node *T_node);
void Tree(struct Node *N, int ramos, int barras, struct V_Node *T_node, struct V_aux *vaux);
void printTree_pre(struct V_Node *T_node);
void printTree_pos(struct V_Node *T_node);
void printTree_pont(struct V_Node *T_node);
void Apagar_Tree(struct V_Node *T_node);
float Tree_pos_corrente_real(struct V_Node *T_node, struct Node *N, struct Tensao_Carga *TC, int ramos);
float Tree_pos_corrente_imag(struct V_Node *T_node, struct Node *N, struct Tensao_Carga *TC, int ramos);
void Tree_pre_tensao_real(struct V_Node *T_node, struct Node *N, struct Tensao_Carga *TC, int ramos);
void Tree_pre_tensao_imag(struct V_Node *T_node, struct Node *N, struct Tensao_Carga *TC, int ramos);
float Perdas_Ativas(struct V_Node *T_node, struct Node *N, int ramos, int barras);
float Perdas_Reativas(struct V_Node *T_node, struct Node *N, int ramos, int barras);
void FCR(struct V_Node *T_node, struct Node *N, struct Tensao_Carga *TC, int barras, int ramos, float *p_at_2, float *p_reat_2,
    float *delta_P, float *P_per_1, float *P_per_2, float lim, int v_ini);
void Tree_pre_tensao(struct V_Node *T_node);
//void Busca_Ciclos(struct V_Node *T_node, struct Node *N, int ramos, int vet_ciclos[], int v1, int v2, int subst, int vet_ciclos_aux[][14], int *i_aux, int *j_aux);
void Ciclo_v1(struct V_Node *T_node, struct Node *N, int ramos, int v1, int v2, int subst, int vet_ciclos_aux[][VETMAX],
    int *i_aux, int *j_aux);
void Ciclo_v2(struct V_Node *T_node, struct Node *N, int ramos, int v1, int v2, int subst, int vet_ciclos_aux[][VETMAX],
    int *i_aux, int *j_aux);
void Busca_Ciclos(struct V_Node *T_node, struct Node *N, int ramos, int barras, int subst, int **M_ciclos, int vet_ciclos[]);
void PR_1(struct V_Node *T_node, struct Node *N, struct Tensao_Carga *TC, struct V_aux *vaux, int ramos, int barras, int subst,
    int **M_ciclos, int vet_ciclos[], float *p_at_2, float *p_reat_2, float *delta_P, float *P_per_1,float *P_per_2, float lim,
    int v_ini, float p_base, float *p_at_top, float *p_reat_top, struct Node M[][16], int n_pr, float FC_ini[]);
void PR_2(struct V_Node *T_node, struct Node *N, struct Tensao_Carga *TC, struct V_aux *vaux, int ramos, int barras, int subst,
    int **M_ciclos, int vet_ciclos[], float *p_at_2, float *p_reat_2, float *delta_P, float *P_per_1,float *P_per_2, float lim,
    int v_ini, float p_base, float *p_at_top, float *p_reat_top, struct Node M[][16],int n_pr);
void Viz_1(struct V_Node *T_node, struct Node *N, struct Node *N_best, struct Node *N_worst,struct Tensao_Carga *TC, struct V_aux *vaux,
    int ramos, int barras, int subst, int **M_ciclos, int vet_ciclos[], float *p_at_2, float *p_reat_2, float *delta_P, float *P_per_1,
    float *P_per_2, float lim, int v_ini, float p_base, float *p_at_top, float *p_reat_top, struct Node M_pr[][16], int *n_pr,
    int *inicio_pr, int vet_block[], float FC_ini[], bool *achei_1, float *soma_fluxo_princ, FILE *arq3, FILE *arq5, FILE *arq6,
    int *bt_cont, float *perda_worst, float *p_at_aux, float *p_reat_aux, float *perdas_pr, int *p1, int *p2);
void Viz_2(struct V_Node *T_node, struct Node *N, struct Node *N_best, struct Node *N_worst,struct Tensao_Carga *TC, struct V_aux *vaux,
    int ramos, int barras, int subst, int **M_ciclos, int vet_ciclos[], float *p_at_2, float *p_reat_2, float *delta_P, float *P_per_1,
    float *P_per_2, float lim, int v_ini, float p_base, float *p_at_top, float *p_reat_top, struct Node M_pr[][16], int *n_pr,
    int *inicio_pr, int vet_block[], float FC_ini[], bool *achei_2, float *soma_fluxo_princ, FILE *arq3, FILE *arq5, FILE *arq6,
    int *bt_cont, float *perda_worst, float *p_at_aux, float *p_reat_aux, float *perdas_pr, int *p1, int *p2);
//void Reiniciar(struct Node *N, struct Tensao_Carga *TC, int ramos, int barras, int v_ini);

int main(int argc, char**argv){

    clock_t T_sol_ini_1, T_sol_ini_2, T_total_1, T_total_2, T_BT_1, T_BT_2, T_BT_proc1, T_BT_proc2, T_BT_proc3, T_PR_1, T_PR_2;
    double Tempo_sol_ini, Tempo_total, Tempo_BT, Tempo_BT_proc1, Tempo_BT_proc2, Tempo_PR;

    int barras, ramos, subst, i, j, k, k1, k2, m;
    //qtde de barras, qtde de ramos, nó da subestação

    float v_init_real, v_init_imag, v_base, p_base, z_base, P_per_1 = 0.0, P_per_2, delta_P = 1.0, p_at, p_reat;
    float p_at_aux, p_reat_aux, p_at_aux_2, p_reat_aux_2, p_at_top, p_reat_top;
    double lim;
    //tensao inicial real e imaginaria, tensao base do sistema em kva
    //potencia base em kva, tolerância de perdas e perdas ativas do sistema

    float p_at_2, p_reat_2;
    //potências ativas e reativas referentes a cada resultado de uma solução vizinha

    float aux1, aux2, aux3; //variáveis de auxílio a leitura e tratamento de dados
    int aux4, aux5, aux6, bt_cont;   //variáveis de auxílio a leitura e tratamento de dados

    //Escolha do formato de trabalho da Impedância base
    //-----------------------------------------------------------------------------------------------------------//
    int op = 0;
    int cont_ramos;

    if( argc < 2) return -1;
    op = atoi(argv[1]);
    switch(op)
    {
        case 1:
            printf("\n***** op 1 - Formato de Porcentagem.\n");
            break;

        case 2:
            printf("\n***** 2 - Formato Comumgem.\n");
            break;

        default:
                return -2; // so aceita 1 ou 2
    };  // switch()

//    while((op!=1)&&(op!=2)){
//        system("cls");
//        printf("\nEscolha a forma de tratamento para a Impedancia Base - Z_base.");
//        printf("\n1 - Formato de Porcentagem.");
//        printf("\n2 - Formato comum.");
//        printf("\n\nDigite sua opcao: ");
//        scanf("%d",&op);
//    }

    //leitura doa dados dos arquivos e preparação das variáveis
    //-----------------------------------------------------------------------------------------------------------//
    FILE *arq1;
    FILE *arq2;

    system("cls");

    arq1 = fopen("14_barras.dat","r");
    arq2 = fopen("FC_14_barras.dat","r");

	fscanf(arq1, "%d", &barras);
	fscanf(arq1, "%d", &ramos);
	fscanf(arq1, "%d", &subst);
	fscanf(arq1, "%f", &v_init_real);
	v_init_imag = 0.0;
	fscanf(arq1, "%f", &v_base);
	fscanf(arq1, "%f", &p_base);
	fscanf(arq1, "%f", &lim);

	//Solução corrente, cada solução vizinha por iteração, melhor solução encontrada na vizinhança
	struct Node N[ramos], N_aux[ramos], N_best[ramos], N_worst[ramos], N_top[ramos];

	struct Tensao_Carga TC[barras];
	float FC_ini[ramos];//vetor de fluxo de carga inicial de cada ramo
	//int M_ciclos[ramos-(barras-1)][ramos];
	int vet_ciclos[ramos], vet_block[ramos];
	v_aux vaux = NULL;

	//Matriz que armazena os ciclos da solução corrente - geradora - de cada iteração da BT
	int **M_ciclos;
    M_ciclos = (int**)malloc ((ramos-(barras-1)) * sizeof(int *));
    for (int i = 0; i < (ramos-(barras-1)); ++i){
      M_ciclos[i] = (int*)malloc (ramos * sizeof (int));
    }

    //Matriz que armazena as 10 melhores soluções encontradas no procedimento GRASP/BT para utilizar no Path Relinking
    //A inserção de cada solução é na forma de uma fila, mas de forma estática
    /*struct Node **M_pr;
    M_pr = (struct Node**)malloc (10 * sizeof(struct Node *));
    for (int i = 0; i < ramos; ++i){
      M_pr[i] = (struct Node*)malloc (ramos * sizeof (struct Node));
    }*/

    //struct Node *M_pr;
    //M_pr = (struct Node *) malloc(10 * ramos *sizeof(struct Node));
    struct Node M_pr[10][16];
    int inicio_pr, n_pr, it_pr; //Variáveis que armazenam a posição de início da fila, quantos elementos estão presentes na mesma, e o iterador
    float perdas_pr; //Armazena o valor de perdas elétricas da última solução a entrar na matriz/file
    inicio_pr = 0;
    n_pr = 0;

	//teste de leitura dos dados do arquivo
	/*printf("\nNro barras: %d",barras);
	printf("\nNro ramos: %d",ramos);
	printf("\nNo da subestação: %d",subst);
	printf("\nTensao de referencia: %f",v_init_real);
	printf("\nTensao base do sistema: %f",v_base);
	printf("\nPotencia base do sistema: %f",p_base);
	printf("\nTolerancia: %f",lim);

	system("pause");*/

    //preenchendo o vetor de potências e tensões de cada barra, vetor 1 a n (barras), com seus valores iniciais
	for (i=0;i<barras;i++){
	    fscanf(arq1, "%d", &j);
	    //salvando na posição j (barra j) seus valores de Pot. At., Reat., e tensões iniciais (real e imaginária)
        fscanf(arq1, "%f", &TC[j-1].pot_at);
        TC[j-1].pot_at = (float)(TC[j-1].pot_at / p_base);
        fscanf(arq1, "%f", &aux1);
        fscanf(arq1, "%f", &aux2);
        //TC[j-1].pot_reat = (float)((aux1 - aux2) / p_base); //Q_ind - Q_cap
        TC[j-1].pot_reat = (float)((aux1-aux2) / p_base);
        TC[j-1].V_real = v_init_real;
        //TC[j-1].V_real = v_base;
        TC[j-1].V_imag = 0.0;

        //printf("\nBarra %d -- Pot. At.: %f --- Pot. Reat.: %f",i+1,TC[i].pot_at,TC[i].pot_reat);
	}
	//system("pause");


	//preenchendo o vetor de informações de ramos existentes, condição de chave de cada ramo
	//impedâncias de cada ramo, e suas correntes
	for(i=0;i<ramos;i++){
        fscanf(arq1, "%d", &N[i].barra_1);
        fscanf(arq1, "%d", &N[i].barra_2);

        //mantém em N[i].barra_1 o menor índice da aresta
        if(N[i].barra_2 < N[i].barra_1){
            aux4 = N[i].barra_1;
            N[i].barra_1 = N[i].barra_2;
            N[i].barra_2 = aux4;
        }

        //Setando a configuração do sistema igual à do arquivo
        /*if(i<=barras-2){
            N[i].chave = 1;
            //printf("\nChave FECHADA -- %d - %d",N[i].barra_1,N[i].barra_2);
        }else{
            N[i].chave = 0;
            //printf("\nChave ABERTA -- %d - %d",N[i].barra_1,N[i].barra_2);
        }*/

        //system("pause");

        //N[i].chave = 1; //Situação para a solução inicial ser utilizada
        N[i].chave = 0; //Situação para PRIM e GRASP
        fscanf(arq1, "%f", &N[i].resist);
        fscanf(arq1, "%f", &N[i].reat);
        N[i].I_real = 0.0;
        N[i].I_imag = 0.0;

        //Ajustando a Resistencia e Reatancia com a Impedancia Base - z_base
        if(op == 1){
            z_base = 100.0;
            N[i].resist = (float)(N[i].resist / z_base);
            N[i].reat = (float)(N[i].reat / z_base);
        }else{
            z_base = (float)(1000 * (pow(v_base,2) / p_base));
            N[i].resist = (float)(N[i].resist / z_base);
            N[i].reat = (float)(N[i].reat / z_base);
        }

        //printf("\nRamo %d - %d \tResist-Reat \t%f - %f",N[i].barra_1,N[i].barra_2,N[i].resist,N[i].reat);
	}

	//system("pause");

    //estruturando o vetor de fluxo de carga inicial
	for(i=0;i<ramos;i++){
        fscanf(arq2, "%d", &j);
        //setando o valor de fluxo de carga inicial do ramo da posição j
        //referente ao ramo do vetor N, de posição j-1
        fscanf(arq2, "%f", &FC_ini[j-1]);
	}

	//Zerando o vetor de proibição de movimentos
	//Iniciam com 1 para que no início da Busca Tabu este vetor seja percorrido e retorne a 0 todos seus valores
	for(i=0;i<ramos;i++){
        vet_block[i] = 1;
	}

    fclose(arq1);
    fclose(arq2);

    //for (i=0;i<barras;i++){
	    //printf("\nBarra %d -- Pot. At.: %f --- Pot. Reat.: %f",i+1,TC[i].pot_at,TC[i].pot_reat);
	//}
	//-----------------------------------------------------------------------------------------------------------//

	//Início do cálculo de tempo total da proposta
	T_total_1 = clock();

	//Execução do algoritmo de PRIM para a obtenção da primeira solução
	//-----------------------------------------------------------------------------------------------------------//

	T_sol_ini_1 = clock();

    //PRIM(N,FC_ini,ramos,subst,barras);

	GRASP(N,N_best,FC_ini,ramos,subst,barras);


	T_sol_ini_2 = clock();
	Tempo_sol_ini = ((T_sol_ini_2 - T_sol_ini_1) * 1000.0) / CLOCKS_PER_SEC;
    printf("\nTempo de execucao da Solução Inicial: %g ms",Tempo_sol_ini);

    system("pause");

	//Setando informações da raiz da árvore do sistema elétrico com a barra da subestação
    //-----------------------------------------------------------------------------------------------------------//


    VN T_node = NULL, T_node_aux = NULL, T_aux = NULL;
    VN vt1; //ponteiro da árvore do sistema elétrico, variável de auxílio para inserção de novos nós
    vt1 = (VN)malloc(sizeof(V_Node));
    vt1->v = subst; //Setando a raiz da árvore do sistema elétrico com o índice da barra da subestação
    vt1->cont_p = 0;
    for(i=0;i<LINK;i++){
        vt1->p[i] = NULL;
    }
    T_node = vt1;
    T_node_aux = vt1; //Raiz da árvore referente a cada solução vizinha
    //-----------------------------------------------------------------------------------------------------------//

    //Gerando a solução inicial
    //-----------------------------------------------------------------------------------------------------------//
    //T_node->cont_p = 0;
    VTree(N,vaux,ramos,subst,barras,T_node);

    //printf("\n");
    printTree_pre(T_node);

    //printf("\n");
    //printTree_pont(T_node);

    //system("pause");

    //printf("\n");
    //printTree_pos(T_node);
    FCR(T_node,N,TC,barras,ramos,&p_at_2,&p_reat_2,&delta_P,&P_per_1,&P_per_2,lim,v_init_real);

    p_at = p_at_2;
    p_reat = p_reat_2;
    p_at_aux = p_at_2;
    p_reat_aux = p_reat_2;
    printf("\n\nSolucao Inicial !!!");
    printf("\nRamos Desligados: ");
    for(i=0;i<ramos;i++){
        if(N[i].chave == 0){
            //printf(" (%d-%d) ",N[i].barra_1,N[i].barra_2);
            printf(" (%d) ",i);
        }
    }
    printf("\nPerdas Ativas da Sol. Inicial: %f",p_at*p_base);
    printf("\nPerdas Reativas da Sol. Inicial: %f",p_reat*p_base);

    /*cont_ramos = 0;
    for(int c_ramos=0;c_ramos<ramos;c_ramos++){
        if(N[c_ramos].chave == 0){
            cont_ramos++;
        }
    }*/

    //Armazenando a primeira solução na Matriz do PR
    for(i=0;i<ramos;i++){
        M_pr[0][i].barra_1 = N[i].barra_1;
        M_pr[0][i].barra_2 = N[i].barra_2;
        M_pr[0][i].chave = N[i].chave;
        M_pr[0][i].id = N[i].id;
        M_pr[0][i].resist = N[i].resist;
        M_pr[0][i].reat = N[i].reat;
        M_pr[0][i].I_real = N[i].I_real;
        M_pr[0][i].I_imag = N[i].I_imag;
    }
    n_pr++;
    perdas_pr = p_at;

    /*printf("\nRamos Desligados: ");
    for(i=0;i<ramos;i++){
        if(M_pr[0][i].chave == 0){
            //printf(" (%d-%d) ",N[i].barra_1,N[i].barra_2);
            printf(" (%d) ",i);
        }
    }

    system("pause");*/

    //Apagar_Tree(T_node);
    //printTree_pre(T_node);
    //system("pause");
    //printTree_pre(T_node);
    //-----------------------------------------------------------------------------------------------------------//


    FILE *arq3, *arq4, *arq5, *arq6;
    //FILE *arq9, *arq10;
    arq3 = fopen("todas_sols.txt","w");
    arq4 = fopen("melhores_sols.txt","w");
    arq5 = fopen("media_fluxos.txt","w");
    arq6 = fopen("iteracao.txt","w");
    float soma_fluxo, soma_fluxo_princ, soma_worst, perda_worst;
    //soma_worst = 100000; //valor estritamente maior do que qualquer somatório de fluxo

    T_BT_1 = clock();

    //Busca Tabu
    //-----------------------------------------------------------------------------------------------------------//
    int v1, v2, p1, p2;
    bool achei_1, achei_2, parada, worst, infactivel;
    worst = false;

    soma_fluxo_princ = 0.0;
    for(int i_ciclo=0;i_ciclo<ramos;i_ciclo++){
        if(N[i_ciclo].chave == 1){
            soma_fluxo_princ = (float) soma_fluxo_princ + FC_ini[i_ciclo];
        }
    }
    printf("\nMedia de fluxo de carga da solucao inicial: %f\n", (float) soma_fluxo_princ/(barras-1));

    /*for(int j=0;j<barras;j++){
        printf("\nBarra %d --- Tensão: %f",j,TC[j].V_real);
    }
    system("pause");*/

    fprintf (arq3,"%f\n",p_at_2*p_base);
    fprintf (arq4,"%f\n",p_at_2*p_base);
    fprintf (arq5,"%f\n",(float) soma_fluxo_princ/(barras-1));
    fprintf (arq6,"%d\n",1);

    //soma_worst = 100000.0;
    perda_worst = 0.0;
    //system("cls");

    T_BT_proc1 = clock();

    for(bt_cont=0;bt_cont<1000;bt_cont++){

        if(!worst){

        //soma_worst = 100000.0;
        //perda_worst = 0.0;
        //printf("\nIteracao BT %d: \n",bt_cont);
        //printTree_pre(T_node);
        //printf("\n");
        achei_1 = false;
        achei_2 = false;
        //"Zerando" os vetores de ciclos
        //-----------------------------------------------------------------------------------------------------------//

        Busca_Ciclos(T_node,N,ramos,barras,subst,M_ciclos,vet_ciclos);


        //-----------------------------------------------------------------------------------------------------------//

        //Percorre a Matriz de ciclos, linha por linha e coluna por coluna modificando os índices dos ramos da solução
        //geradora com os índices da Matriz de ciclos, identificando os diferentes vizinhos.
        //Para cada solução vizinha se cria sua árvore dinâmica e se avalia seu fluxo de carga radial.
        //Caso as perdas sejam menores que a solução com menor perdas, esta solução é guardada no vetor N_best.
        //Ao final de cada iteração da Busca Tabu volta-se a solução armazenada em N_best para N e se reinicia o processo.
        //A melhor solução então é ilustrada ao usuário.
        //-----------------------------------------------------------------------------------------------------------//
        for(int proib=0;proib<ramos;proib++){
            if(vet_block[proib] != 0){
                vet_block[proib]--;
            }
        }

        /*int contm1 = 1, stotal = 0, contm2 = 0;
        for(int im1=0;im1<(ramos-(barras-1))-1;im1++){
            for(int il1=0;il1<ramos;il1++){
                if((M_ciclos[im1][il1] == 1)||(M_ciclos[im1][il1] == 2)){
                    contm1++;
                }
            }
            for(int im2=im1+1;im2<(ramos-(barras-1));im2++){
                for(int il2=0;il2<ramos;il2++){
                    if((M_ciclos[im2][il2] == 1)||(M_ciclos[im2][il2] == 2)){
                        contm2++;
                    }
                }
                stotal = stotal + (contm1 * contm2);
                contm2 = 0;
            }
            contm1 = 0;
        }
        printf("\nTotal possivel de sols: %d",stotal);
        system("pause");*/


        Viz_1(T_node,N,N_best,N_worst,TC,vaux,ramos,barras,subst,M_ciclos,vet_ciclos,&p_at_2,&p_reat_2,&delta_P,&P_per_1,
            &P_per_2,lim,v_init_real,p_base,&p_at_top,&p_reat_top,M_pr,&n_pr,&inicio_pr,vet_block,FC_ini,&achei_1,&soma_fluxo_princ,
            arq3,arq5,arq6,&bt_cont,&perda_worst,&p_at_aux,&p_reat_aux,&perdas_pr,&p1,&p2);

        /*Viz_2(T_node,N,N_best,N_worst,TC,vaux,ramos,barras,subst,M_ciclos,vet_ciclos,&p_at_2,&p_reat_2,&delta_P,&P_per_1,
            &P_per_2,lim,v_init_real,p_base,&p_at_top,&p_reat_top,M_pr,&n_pr,&inicio_pr,vet_block,FC_ini,&achei_2,&soma_fluxo_princ,
            arq3,arq5,arq6,&bt_cont,&perda_worst,&p_at_aux,&p_reat_aux,&perdas_pr,&p1,&p2);*/

        //system("pause");


        //Atualizando o vetor N com a melhor solução encontrada em uma iteração
        //O melhor vizinho do Vetor N em cada iteração passa a ser a próxima solução inicial
        if(achei_1){

            vet_block[p1] = 3;
            vet_block[p2] = 3;
            //vet_block[ramos]++;
            /*printf("\nIteracao %d\n",bt_cont+1);
            printf("Ramos bloqueados: ");
            for(int proib=0;proib<ramos;proib++){
                if(vet_block[proib] != 0){
                    printf("%d ",proib);
                }
            }
            printf("\n");*/

            Apagar_Tree(T_node);
            for(m=0;m<ramos;m++){
                N[m].barra_1 = N_best[m].barra_1;
                N[m].barra_2 = N_best[m].barra_2;
                N[m].chave = N_best[m].chave;
                N[m].id = N_best[m].id;
                N[m].resist = N_best[m].resist;
                N[m].reat = N_best[m].reat;
                N[m].I_real = N_best[m].I_real;
                N[m].I_imag = N_best[m].I_imag;
                /*if(N[m].chave == 0){
                    printf("\nRamos abertos da solução de bt_cont %d: %d", bt_cont, m);
                }*/
            }
            T_node->cont_p = 0;
            VTree(N,vaux,ramos,subst,barras,T_node);
            delta_P = 1.0;
            FCR(T_node,N,TC,barras,ramos,&p_at_2,&p_reat_2,&delta_P,&P_per_1,&P_per_2,lim,v_init_real);
            fprintf (arq4,"%f\n",p_at_2*p_base);
            p_at_aux = p_at_2;
            p_reat_aux = p_reat_2;

            soma_fluxo_princ = 0.0;
            for(int i_ciclo=0;i_ciclo<ramos;i_ciclo++){
                if(N[i_ciclo].chave == 1){
                    soma_fluxo_princ = (float) soma_fluxo_princ + FC_ini[i_ciclo];
                }
            }

            //printf("\nTerminando BT sem piorar a solução.");
            //printf("\nBreak: %d --- Melhor Perda Ativa: %f",bt_cont+1,p_at_aux*p_base);
            //printf("\nBreak: %d --- Perda Reativa: %f",bt_cont+1,p_reat_aux*p_base);
            //printf("\nBreak: %d --- Media de fluxo: %f\n",bt_cont+1,(float) (soma_fluxo_princ/(barras-1)));
            //t2=clock();
            //Tempo = ((t2 - t1) * 1000.0) / CLOCKS_PER_SEC;
            //printf("\nTempo de execucao: %g ms",Tempo);

            /*cont_ramos = 0;
            for(int c_ramos=0;c_ramos<ramos;c_ramos++){
                if(N[c_ramos].chave == 0){
                    cont_ramos++;
                }
            }
            printf("\nNro de ramos desligados: %d",cont_ramos);*/

            //system("pause");
        }else{
            worst = true;

            //Salvando a melhor solução encontrada ate o momento
            for(m=0;m<ramos;m++){
                N_top[m].barra_1 = N[m].barra_1;
                N_top[m].barra_2 = N[m].barra_2;
                N_top[m].chave = N[m].chave;
                N_top[m].id = N[m].id;
                N_top[m].resist = N[m].resist;
                N_top[m].reat = N[m].reat;
                N_top[m].I_real = N[m].I_real;
                N_top[m].I_imag = N[m].I_imag;
            }
            Apagar_Tree(T_node);
            T_node->cont_p = 0;
            VTree(N,vaux,ramos,subst,barras,T_node);
            delta_P = 1.0;
            FCR(T_node,N_top,TC,barras,ramos,&p_at_2,&p_reat_2,&delta_P,&P_per_1,&P_per_2,lim,v_init_real);
            p_at_top = p_at_2;
            p_reat_top = p_reat_2;

            /*Apagar_Tree(T_node);
            T_node->cont_p = 0;
            VTree(N,vaux,ramos,subst,barras,T_node);
            delta_P = 1.0;
            FCR(T_node,N,TC,barras,ramos,&p_at_2,&p_reat_2,&delta_P,&P_per_1,&P_per_2,lim,v_init_real);
            fprintf (arq4,"%f\n",p_at_2*p_base);*/
            p_at_aux_2 = p_at_aux;
            p_reat_aux_2 = p_reat_aux;

            soma_fluxo = 0.0;
            for(int i_ciclo=0;i_ciclo<ramos;i_ciclo++){
                if(N[i_ciclo].chave == 1){
                    soma_fluxo = (float) soma_fluxo + FC_ini[i_ciclo];
                }
            }

            printf("\nPerdas Ativas antes de piorar a solucao: %f",p_at_aux_2*p_base);
            printf("\nPerdas Reativas antes de piorar a solucao: %f",p_reat_aux_2*p_base);
            printf("\nBreak: %d --- Media de fluxo: %f\n",bt_cont+1,(float) (soma_fluxo/(barras-1)));
            printf("\nRamos desligados: ");
            for(int iter=0;iter<ramos;iter++){
                if(N[iter].chave == 0){
                    printf(" (%d) ",iter);
                }
            }
            printf("\n");

            for(m=0;m<ramos;m++){
                N[m].barra_1 = N_worst[m].barra_1;
                N[m].barra_2 = N_worst[m].barra_2;
                N[m].chave = N_worst[m].chave;
                N[m].id = N_worst[m].id;
                N[m].resist = N_worst[m].resist;
                N[m].reat = N_worst[m].reat;
                N[m].I_real = N_worst[m].I_real;
                N[m].I_imag = N_worst[m].I_imag;
            }
            Apagar_Tree(T_node);
            T_node->cont_p = 0;
            VTree(N,vaux,ramos,subst,barras,T_node);
            delta_P = 1.0;
            FCR(T_node,N,TC,barras,ramos,&p_at_2,&p_reat_2,&delta_P,&P_per_1,&P_per_2,lim,v_init_real);
            fprintf (arq4,"%f\n",p_at_2*p_base);
            p_at_aux = p_at_2;
            p_reat_aux = p_reat_2;

            soma_fluxo_princ = 0.0;
            for(int i_ciclo=0;i_ciclo<ramos;i_ciclo++){
                if(N[i_ciclo].chave == 1){
                    soma_fluxo_princ = (float) soma_fluxo_princ + FC_ini[i_ciclo];
                }
            }

            for(int proib=0;proib<ramos;proib++){
                vet_block[proib] = 1;
            }

            printf("\nPiorando a solucao.");
            printf("\nBreak: %d --- Pior Perda Ativa: %f",bt_cont+1,p_at_aux*p_base);
            printf("\nBreak: %d --- Perda Reativa: %f",bt_cont+1,p_reat_aux*p_base);
            printf("\nBreak: %d --- Media de fluxo: %f\n",bt_cont+1,(float) (soma_fluxo_princ/(barras-1)));
            printf("\nRamos desligados: ");
            for(int iter=0;iter<ramos;iter++){
                if(N[iter].chave == 0){
                    printf(" (%d) ",iter);
                }
            }
            T_BT_proc2 = clock();
            printf("\n");

            //printf("\nENTROU NO BREAK em bt_cont = %d\n",bt_cont+1);
            //break;
        }
    }else{

        //printf("\naqui 1");
        //printf("\nIteracao BT %d: \n",bt_cont);
        //printTree_pre(T_node);
        //printf("\n");
        achei_1 = false;
        //"Zerando" os vetores de ciclos
        //-----------------------------------------------------------------------------------------------------------//

        Busca_Ciclos(T_node,N,ramos,barras,subst,M_ciclos,vet_ciclos);

        //-----------------------------------------------------------------------------------------------------------//

        //Percorre a Matriz de ciclos, linha por linha e coluna por coluna modificando os índices dos ramos da solução
        //geradora com os índices da Matriz de ciclos, identificando os diferentes vizinhos.
        //Para cada solução vizinha se cria sua árvore dinâmica e se avalia seu fluxo de carga radial.
        //Caso as perdas sejam menores que a solução com menor perdas, esta solução é guardada no vetor N_best.
        //Ao final de cada iteração da Busca Tabu volta-se a solução armazenada em N_best para N e se reinicia o processo.
        //A melhor solução então é ilustrada ao usuário.
        //-----------------------------------------------------------------------------------------------------------//
        for(int proib=0;proib<ramos;proib++){
            if(vet_block[proib] != 0){
                vet_block[proib]--;
            }
        }

        Viz_1(T_node,N,N_best,N_worst,TC,vaux,ramos,barras,subst,M_ciclos,vet_ciclos,&p_at_2,&p_reat_2,&delta_P,&P_per_1,
            &P_per_2,lim,v_init_real,p_base,&p_at_top,&p_reat_top,M_pr,&n_pr,&inicio_pr,vet_block,FC_ini,&achei_1,&soma_fluxo_princ,
            arq3,arq5,arq6,&bt_cont,&perda_worst,&p_at_aux,&p_reat_aux,&perdas_pr,&p1,&p2);

        //Atualizando o vetor N com a melhor solução encontrada em uma iteração
        //O melhor vizinho do Vetor N em cada iteração passa a ser a próxima solução inicial
        if(achei_1){

            vet_block[p1] = 3;
            vet_block[p2] = 3;
            //vet_block[ramos]++;
            /*printf("\nIteracao %d\n",bt_cont+1);
            printf("Ramos bloqueados: ");
            for(int proib=0;proib<ramos;proib++){
                if(vet_block[proib] != 0){
                    printf("%d ",proib);
                }
            }
            printf("\n");*/

            Apagar_Tree(T_node);
            for(m=0;m<ramos;m++){
                N[m].barra_1 = N_best[m].barra_1;
                N[m].barra_2 = N_best[m].barra_2;
                N[m].chave = N_best[m].chave;
                N[m].id = N_best[m].id;
                N[m].resist = N_best[m].resist;
                N[m].reat = N_best[m].reat;
                N[m].I_real = N_best[m].I_real;
                N[m].I_imag = N_best[m].I_imag;
                /*if(N[m].chave == 0){
                    printf("\nRamos abertos da solução de bt_cont %d: %d", bt_cont, m);
                }*/
            }
            T_node->cont_p = 0;
            VTree(N,vaux,ramos,subst,barras,T_node);
            delta_P = 1.0;
            FCR(T_node,N,TC,barras,ramos,&p_at_2,&p_reat_2,&delta_P,&P_per_1,&P_per_2,lim,v_init_real);
            fprintf (arq4,"%f\n",p_at_2*p_base);
            //p_at_aux_2 = p_at_2;
            //p_reat_aux_2 = p_reat_2;
            p_at_aux = p_at_2;
            p_reat_aux = p_reat_2;

            soma_fluxo_princ = 0.0;
            for(int i_ciclo=0;i_ciclo<ramos;i_ciclo++){
                if(N[i_ciclo].chave == 1){
                    soma_fluxo_princ = (float) soma_fluxo_princ + FC_ini[i_ciclo];
                }
            }

            //printf("\nBreak: %d --- Melhor Perda Ativa: %f",bt_cont+1,p_at_aux*p_base);
            //printf("\nBreak: %d --- Perda Reativa: %f",bt_cont+1,p_reat_aux*p_base);
            //printf("\nBreak: %d --- Media de fluxo: %f\n",bt_cont+1,(float) (soma_fluxo_princ/(barras-1)));
            /*cont_ramos = 0;
            for(int c_ramos=0;c_ramos<ramos;c_ramos++){
                if(N[c_ramos].chave == 0){
                    cont_ramos++;
                }
            }
            printf("\nNro de ramos desligados: %d",cont_ramos);*/

            //Inserindo a melhor solução da iteração corrente na matriz/fila das melhores soluções


            //system("pause");
        }else{
            printf("\nENTROU NO BREAK em bt_cont = %d\n",bt_cont+1);
            break;
        }
    }
    }
    fclose(arq3);
    fclose(arq4);
    fclose(arq5);
    fclose(arq6);

    T_BT_proc3 = clock();
    Tempo_BT_proc1 = ((T_BT_proc2 - T_BT_proc1) * 1000.0) / CLOCKS_PER_SEC;
    Tempo_BT_proc2 = ((T_BT_proc3 - T_BT_proc2) * 1000.0) / CLOCKS_PER_SEC;
    printf("\nTempo de execucao da BT - Procedimento 1: %g ms",Tempo_BT_proc1);
    printf("\nTempo de execucao da BT - Procedimento 2: %g ms",Tempo_BT_proc2);


    T_BT_2 = clock();
    Tempo_BT = ((T_BT_2 - T_BT_1) * 1000.0) / CLOCKS_PER_SEC;
    printf("\nTempo de execucao da BT: %g ms",Tempo_BT);


    Apagar_Tree(T_node);
    T_node->cont_p = 0;
    VTree(N,vaux,ramos,subst,barras,T_node);
    delta_P = 1.0;
    FCR(T_node,N,TC,barras,ramos,&p_at_2,&p_reat_2,&delta_P,&P_per_1,&P_per_2,lim,v_init_real);
    p_at = p_at_2;
    p_reat = p_reat_2;
    //Ticks[1] = clock();
    //t2 = clock();

    //Escolhendo entre a melhor solução encontrada ANTES de piorar a solução e ao final DEPOIS de piorar a solução
    if(p_at_top < p_at){

        Apagar_Tree(T_node);
        T_node->cont_p = 0;
        VTree(N_top,vaux,ramos,subst,barras,T_node);
        delta_P = 1.0;
        FCR(T_node,N_top,TC,barras,ramos,&p_at_2,&p_reat_2,&delta_P,&P_per_1,&P_per_2,lim,v_init_real);

        soma_fluxo_princ = 0.0;
        for(int i_ciclo=0;i_ciclo<ramos;i_ciclo++){
            if(N_top[i_ciclo].chave == 1){
                soma_fluxo_princ = (float) soma_fluxo_princ + FC_ini[i_ciclo];
            }
        }

        printf("\nPerdas Ativas da Sol. Final: %f",p_at_top*p_base);
        printf("\nPerdas Reativas da Sol. Final: %f",p_reat_top*p_base);
        printf("\nBreak: %d --- Media de fluxo: %f",bt_cont+1,(float) (soma_fluxo_princ/(barras-1)));
        printf("\nRamos deligados: ");
        for(i=0;i<ramos;i++){
            if(N_top[i].chave == 0){
                //printf(" (%d-%d) ",N[i].barra_1,N[i].barra_2);
                printf(" (%d) ",i);
            }
        }
        for(int j=0;j<barras;j++){
            printf("\nBarra %d --- Tensão: %f",j,TC[j].V_real);
        }
    }else{

        soma_fluxo_princ = 0.0;
        for(int i_ciclo=0;i_ciclo<ramos;i_ciclo++){
            if(N[i_ciclo].chave == 1){
                soma_fluxo_princ = (float) soma_fluxo_princ + FC_ini[i_ciclo];
            }
        }

        printf("\nPerdas Ativas da Sol. Final: %f",p_at*p_base);
        printf("\nPerdas Reativas da Sol. Final: %f",p_reat*p_base);
        printf("\nBreak: %d --- Media de fluxo: %f",bt_cont+1,(float) (soma_fluxo_princ/(barras-1)));
        printf("\nRamos deligados: ");
        for(i=0;i<ramos;i++){
            if(N[i].chave == 0){
                //printf(" (%d-%d) ",N[i].barra_1,N[i].barra_2);
                printf(" (%d) ",i);
            }
        }
        /*for(int j=0;j<barras;j++){
            printf("\nBarra %d --- Tensão: %f",j,TC[j].V_real);
        }*/

        //Passando a solução para N_top, para N_top guardar sempre a melhor solução
        for(m=0;m<ramos;m++){
            N_top[m].barra_1 = N[m].barra_1;
            N_top[m].barra_1 = N[m].barra_2;
            N_top[m].barra_1 = N[m].chave;
            N_top[m].barra_1 = N[m].id;
            N_top[m].barra_1 = N[m].resist;
            N_top[m].barra_1 = N[m].reat;
            N_top[m].barra_1 = N[m].I_real;
            N_top[m].barra_1 = N[m].I_imag;
        }
        p_at_top = p_at;
        p_reat_top = p_reat;
    }
    //system("pause");

    //system("cls");
    printf("\n\nPartindo para PR!!!!!!");
    //system("pause");

    T_PR_1 = clock();

    PR_1(T_node,N_top,TC,vaux,ramos,barras,subst,M_ciclos,vet_ciclos,&p_at_2,&p_reat_2,&delta_P,&P_per_1,
       &P_per_2,lim,v_init_real,p_base,&p_at_top,&p_reat_top,M_pr,n_pr,FC_ini);

    T_PR_2 = clock();
    Tempo_PR = ((T_PR_2 - T_PR_1) * 1000.0) / CLOCKS_PER_SEC;
    printf("\nTempo de execucao do PR: %g ms",Tempo_PR);

    Tempo_total = ((T_PR_2 - T_total_1) * 1000.0) / CLOCKS_PER_SEC;
    printf("\nTempo de execucao TOTAL: %g ms",Tempo_total);

    //system("pause");
    //system("cls");
    //printf("\n\nPartindo para PR_2!!!!!!");
    //system("pause");

    /*PR_2(T_node,N_top,TC,vaux,ramos,barras,subst,M_ciclos,vet_ciclos,&p_at_2,&p_reat_2,&delta_P,&P_per_1,
       &P_per_2,lim,v_init_real,p_base,&p_at_top,&p_reat_top,M_pr,n_pr);*/

    //printf("\nPerdas Ativas da Melhor Sol.: %f",p_at_top*p_base);


    printf("\nPerdas Ativas da Sol. Final: %f",p_at_top*p_base);
    printf("\nPerdas Reativas da Sol. Final: %f",p_reat_top*p_base);
    printf("\nRamos Desligados: ");
    for(i=0;i<ramos;i++){
        if(N_top[i].chave == 0){
            //printf(" (%d-%d) ",N[i].barra_1,N[i].barra_2);
            printf(" (%d) ",i);
        }
    }
    //printf("\nBreak: %d --- Media de fluxo: %f",bt_cont+1,(float) (soma_fluxo_princ/(barras-1)));
    //printf("\nRamos deligados: ");

    //Tempo = ((t2 - t1) * 1000.0) / CLOCKS_PER_SEC;
    //printf("\n\nTempo de execucao: %g ms\n",Tempo);
    system("pause");

}

void PRIM(struct Node *N, float FC_ini[], int ramos, int subst,int barras){
    int i, j, k, m, n, b1, b2, v, cont = 0;
    int cont2 = 0 ; //indica a qtde de ramos iniciais que irão possuir a chave ligada
    float dmin;
    bool Found;
    dmin = 1.0; //valor estritamente mais baixo do que qualquer valor no vetor de fluxo de carga

    //Escolhendo o ramo inicial, de menor caminho de fluxo de carga inicial, a partir da subestação
    //-----------------------------------------------------------------------------------------------------------//
    for(i=0;i<ramos;i++){
        if((N[i].barra_1 == subst)||(N[i].barra_2 == subst)){
            if(FC_ini[i] > dmin){
                dmin = FC_ini[i];
                j = i;
            }
        }
    }
    //Setando ramo ativo inicial, de menor caminho, à partir da subestação
    N[j].chave = 1;
    //printf("\nRamo %d a entrar na solucao: %d-%d --- FC_ini. %f",j+1,N[j].barra_1,N[j].barra_2,FC_ini[j]);
    cont2 = 1;
    //-----------------------------------------------------------------------------------------------------------//

    //Setando com chave ligada todos os ramos que possuem a subestação
    //-----------------------------------------------------------------------------------------------------------//
    /*for(i=0;i<ramos;i++){
        if((N[i].barra_1 == subst)||(N[i].barra_2 == subst)){
            N[i].chave = 1;
            cont2++;
        }
    }*/
    //-----------------------------------------------------------------------------------------------------------//

    //Executa a qtde de ramos/arestas restantes a entrar na árvore
    //A cada execução se avaliam todos os possíveis ramos e se escolhe somente o que possui o menor custo
    for(n=0;n<barras-cont2-1;n++){
        //Percorre e compara ramo a ramo se existe uma ligação entre eles
        //e se um ramo conectado esta se conectando com um ramo desligado
        dmin = 0.0;
        for(i=0;i<ramos;i++){
            for(j=0;j<ramos;j++){
                if(i != j){ //Evita-se tratar ramos idênticos
                    if(N[i].chave == 1){
                        //Avaliando se o ramo "j" encontrado possui chave desligada
                        //e se está conectado ao ramo "i" com chave ligada
                        if(((N[i].barra_1 == N[j].barra_1)||(N[i].barra_1 == N[j].barra_2)||
                            (N[i].barra_2 == N[j].barra_1)||(N[i].barra_2 == N[j].barra_2))
                            &&(N[j].chave == 0)){
                                //printf("\nAvaliando ramo %d -- FC: %f",j+1,FC_ini[j]);
                                if(FC_ini[j] >= dmin){
                                    cont = 0;
                                    //Para o ramo "j" com chave desligada encontrado, e válido (FC_ini[j]>=dmin)
                                    //identifica-se o índice da barra que não é comum aos 2 ramos "i" e "j"
                                    if((N[i].barra_1 == N[j].barra_1)||(N[i].barra_2 == N[j].barra_1)){
                                        v = N[j].barra_2;
                                    }else{
                                        v = N[j].barra_1;
                                    }
                                    //Percorre-se todos os ramos se verificando, a partir da barra de índice "v" que não é
                                    //comum aos 2 ramos "i" e "j", se existem ramos "k", ligados ao ramo "j" com chave ligada
                                    //com a barra de índice "v", no intuito de se evitar ciclos
                                    for(k=0;k<ramos;k++){
                                        if(k != j){
                                            if(((N[k].barra_1 == v)||(N[k].barra_2 == v))&&(N[k].chave == 1)){
                                                cont++;
                                                break;
                                                //Caso um ramo "k" esteja com chave ligada o programa força a saída do laço
                                                //pois o ramo "j" se conecta a um ramo "k" que possui chave ligada
                                                //e portanto cria um ciclo
                                            }
                                        }
                                    }
                                    //Caso o contador "cont" seja maior que zero, indica que se obteve um ciclo, e que o ramo "j"
                                    //não é mais válido para escolha, então força-se novamente o programa a quebrar o laço
                                    //Caso o contador tenha continuado NULO, se atualiza dmin, e guarda-se o índice do ramo "j"
                                    if(cont>0){
                                        continue;
                                    }else{
                                        dmin = FC_ini[j];
                                        m = j;
                                    }
                                }
                        }
                    }
                }
            }
        }
        //printf("\nRamo %d a entrar na solucao: %d-%d --- FC_ini. %f",m+1,N[m].barra_1,N[m].barra_2,FC_ini[m]);
        N[m].chave = 1;
    }
}

void GRASP(struct Node *N, struct Node *N_best, float FC_ini[], int ramos, int subst,int barras){
    int i, j, j1, k, m, n, b1, b2, v, cont = 0, i_max, i_min;
    float RCL[ramos], S_fc_1 = 0, S_fc_2; //vetor da lista restrita, somatórios de fluxo de carga
    float m_Scf_1 = 0 , m_Scf_2 = 0, m_Scf_aux = 0; //médias de fluxos de carga
    int i_rcl, n_rcl, index_rcl[ramos], index, r;
    int it, cont2; //indica a qtde de ramos iniciais que irão possuir a chave ligada
    float dmin, alpha = 0.4;
    bool Found;
    int parada;

    FILE *arq5;

    arq5 = fopen("GRASP_sols.txt","w");

    for(i=0;i<ramos;i++){
        S_fc_1 = S_fc_1 + FC_ini[i];
    }
    //printf("\nSomatorio total de fluxos de carga: %f\n",S_fc_1);
    m_Scf_1 = (float) S_fc_1 / ramos; //Guardando a média da soma dos fluxos de carga de TODOS os ramos
    //printf("\nMedia de total de fluxos de carga: %f\n",m_Scf_1);

    //while(m_Scf_1 != 677.583031){
    //do{
    for(it=0;it<10;it++){

        //Desligando todas as chaves dos ramos
        for(i=0;i<ramos;i++){
            N[i].chave = 0;
        }

        //Buscando o primeiro ramo a possuir chave fechada
        //-----------------------------------------------------------------------------------------------------------//
        cont2 = 0;
        dmin = 0.0; //valor estritamente mais baixo do que qualquer valor no vetor de fluxo de carga
        Found = false;
        index = 0;
        n_rcl = 0;
        m_Scf_aux = m_Scf_1;
        for(i=0;i<ramos;i++){
            if((N[i].barra_1 == subst)||(N[i].barra_2 == subst)){
                if(n_rcl < ramos){
                    i_rcl = n_rcl;
                    //Percorrendo o vetor de índices de ramos escolhidos a entrar na lista restrita
                    //a fim de que não existam ramos repetidos
                    for(j1=0;j1<index;j1++){
                        if(index_rcl[index] == i)
                        Found = true;
                    }
                    //Não encontrando índices repetidos de ramos no vetor de índices de ramos
                    //possíveis para a atual iteração do GRASP, insere-se todos os fluxos de
                    //carga dos ramos possíveis da solução em ordem decrescente
                    if(Found == false){
                        while((i_rcl > 0) && (RCL[i_rcl-1] < FC_ini[i])){
                            RCL[i_rcl] = RCL[i_rcl-1];
                            index_rcl[i_rcl] = index_rcl[i_rcl-1];
                            i_rcl--;
                        }
                        RCL[i_rcl] = FC_ini[i]; //Insere na lista restrita os fluxos de carga
                        n_rcl++;
                        index_rcl[i_rcl] = i; //Insere na lista os índices dos ramos escolhidos
                        index++;
                    }
                }
            }
        }

        //srand(time(0));
        srand((unsigned)time(NULL));
        //printf("\nGRASP solucao %d.",it+1);

        /*for(i=0;i<index;i++){
            printf("\nRamo Candidato - Fluxo -- na lista: %d - %f",index_rcl[i]+1,RCL[i]);
        }

        printf("\n");*/
        i_max = 0;
        i_min = index-1;
        i=0;
        while(i!=index){
            if((RCL[i] <= RCL[i_max]) && (RCL[i] >= (float)(RCL[i_min] + alpha*(RCL[i_max]-RCL[i_min])))){
                i++;
            }else{
                for(j=i;j<=index;j++){
                    RCL[j]=RCL[j+1];
                    index_rcl[j]=index_rcl[j+1];
                    index--;
                }
            }
        }

        /*for(i=0;i<index;i++){
            printf("\nRamos possiveis - Fluxo -- na lista: %d - %f",index_rcl[i]+1,RCL[i]);
        }

        printf("\n");*/
        i = rand()%index;
        N[index_rcl[i]].chave = 1;
        //printf("\nRamo Inicial/Fluxo escolhido: %d / %f\n",index_rcl[i]+1,RCL[i]);

        //system("pause");

        /*index = 0;
        for(i=0;i<n_rcl;i++){
            printf("\nFluxo do ramo %d, %d-%d escolhido: %f",index_rcl[index]+1,N[index_rcl[index]].barra_1,N[index_rcl[index]].barra_2,RCL[i]);
            index++;
        }*/

        cont2++;
        //-----------------------------------------------------------------------------------------------------------//

        //Executa a qtde de ramos/arestas restantes a entrar na árvore
        //A cada execução se avaliam todos os possíveis ramos e se escolhe uma percela dos ramos de maior fluxo de carga
        for(n=0;n<barras-cont2-1;n++){
            //Percorre e compara ramo a ramo se existe uma ligação entre eles
            //e se um ramo conectado esta se conectando com um ramo desligado
            dmin = 1.0;
            Found = false;
            index = 0;
            n_rcl = 0;
            m_Scf_aux = m_Scf_1;
            for(i=0;i<ramos;i++){
                for(j=0;j<ramos;j++){
                    if(i != j){ //Evita-se tratar ramos idênticos
                        if(N[i].chave == 1){
                            //Avaliando se o ramo "j" encontrado possui chave desligada
                            //e se está conectado ao ramo "i" com chave ligada
                            if(((N[i].barra_1 == N[j].barra_1)||(N[i].barra_1 == N[j].barra_2)||
                                (N[i].barra_2 == N[j].barra_1)||(N[i].barra_2 == N[j].barra_2))
                                &&(N[j].chave == 0)){
                                cont = 0;
                                //Para o ramo "j" com chave desligada encontrado, e válido (FC_ini[j]>=dmin)
                                //identifica-se o índice da barra que não é comum aos 2 ramos "i" e "j"
                                if((N[i].barra_1 == N[j].barra_1)||(N[i].barra_2 == N[j].barra_1)){
                                    v = N[j].barra_2;
                                }else{
                                    v = N[j].barra_1;
                                }
                                //Percorre-se todos os ramos se verificando, a partir da barra de índice "v" que não é
                                //comum aos 2 ramos "i" e "j", se existem ramos "k", ligados ao ramo "j" com chave ligada
                                //com a barra de índice "v", no intuito de se evitar ciclos
                                for(k=0;k<ramos;k++){
                                    if(k != j){
                                        if(((N[k].barra_1 == v)||(N[k].barra_2 == v))&&(N[k].chave == 1)){
                                            cont++;
                                            break;
                                            //Caso um ramo "k" esteja com chave ligada o programa força a saída do laço
                                            //pois o ramo "j" se conecta a um ramo "k" que possui chave ligada
                                            //e portanto cria um ciclo
                                        }
                                    }
                                }
                                //Caso o contador "cont" seja maior que zero, indica que se obteve um ciclo, e que o ramo "j"
                                //não é mais válido para escolha, então força-se novamente o programa a quebrar o laço
                                //Caso o contador tenha continuado NULO, se atualiza dmin, e guarda-se o índice do ramo "j"
                                if(cont>0){
                                    continue;
                                }else{
                                    if(n_rcl < ramos){
                                        i_rcl = n_rcl;
                                        //Percorrendo o vetor de índices de ramos escolhidos a entrar na lista restrita
                                        //a fim de que não existam ramos repetidos
                                        Found = false;
                                        for(j1=0;j1<n_rcl;j1++){
                                            if(index_rcl[j1] == j){
                                                Found = true;
                                            }
                                        }
                                        //Não encontrando índices repetidos de ramos no vetor de índices de ramos
                                        //possíveis para a atual iteração do GRASP, insere-se todos os fluxos de
                                        //carga dos ramos possíveis da solução em ordem decrescente
                                        if(Found == false){
                                            while((i_rcl > 0) && (RCL[i_rcl-1] < FC_ini[j])){
                                                RCL[i_rcl] = RCL[i_rcl-1];
                                                index_rcl[i_rcl] = index_rcl[i_rcl-1];
                                                i_rcl--;
                                            }
                                            RCL[i_rcl] = FC_ini[j]; //Insere na lista restrita os fluxos de carga
                                            n_rcl++;
                                            index_rcl[i_rcl] = j; //Insere na lista os índices dos ramos escolhidos
                                            index++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            //printf("\nRamo %d a entrar na solucao: %d-%d --- FC_ini. %f",m+1,N[m].barra_1,N[m].barra_2,FC_ini[m]);

            //srand(time(0));
            srand((unsigned)time(NULL));

            //Ajustando os vetores RCL e index_rcl pela fórmula do GRASP
            /*for(i=0;i<index;i++){
            printf("\nRamo Candidato - Fluxo -- na lista: %d - %f",index_rcl[i]+1,RCL[i]);
            }

            printf("\n");*/
            i_max = 0;
            i_min = index-1;
            i=0;
            while(i!=index){
                if((RCL[i] <= RCL[i_max]) && (RCL[i] >= (float)(RCL[i_min] + alpha*(RCL[i_max]-RCL[i_min])))){
                    i++;
                }else{
                    for(j=i;j<=index;j++){
                        RCL[j]=RCL[j+1];
                        index_rcl[j]=index_rcl[j+1];
                        index--;
                    }
                }
            }

            /*for(i=0;i<index;i++){
                printf("\nRamos possiveis - Fluxo -- na lista: %d - %f",index_rcl[i]+1,RCL[i]);
            }

            printf("\n");*/
            i = rand()%index;
            N[index_rcl[i]].chave = 1;
            //printf("\nRamo Inicial/Fluxo escolhido: %d / %f\n",index_rcl[i]+1,RCL[i]);
            //system("pause");
        }

        //Aqui faz-se o calculo para ver qual solução entra em Nbest
        //-----------------------------------------------------------------------------------------------------------//
        //Somatório e média dos fluxos de carga dos ramos ativos
        S_fc_2 = 0;
        printf("\nGRASP solucao %d.",it+1);
        printf("\nRamos fechados: ");
        for(j1=0;j1<ramos;j1++){
            if(N[j1].chave == 1){
                S_fc_2 = S_fc_2 + FC_ini[j1];
            }else{
                //printf("(%d-%d) ",N[j1].barra_1,N[j1].barra_2);
                printf("(%d) ",j1);
            }
        }
        m_Scf_2 = (float) S_fc_2 / (barras-1);
        //m_Scf_1 = (float) S_fc_2 / (barras-1);
        printf("\nMedia FC: %f\n",m_Scf_2);
        //printf("\nMedia FC: %f\n",m_Scf_1);

        fprintf (arq5,"%f\n",m_Scf_2);

        //Sleep(300);
        //system("pause");

        //Caso a média dos fluxos de carga dos ramos ativos seja MAIOR que a média do fluxo de carga
        //de TODOS os ramos, a solução é reservada
        if(m_Scf_2 > m_Scf_1){
            m_Scf_1 = m_Scf_2;
            for(m=0;m<ramos;m++){
                N_best[m].barra_1 = N[m].barra_1;
                N_best[m].barra_2 = N[m].barra_2;
                N_best[m].chave = N[m].chave;
                N_best[m].id = N[m].id;
                N_best[m].resist = N[m].resist;
                N_best[m].reat = N[m].reat;
                N_best[m].I_real = N[m].I_real;
                N_best[m].I_imag = N[m].I_imag;
            }
        }

        /*printf("\nDigite a parada 1-Continua -- 0-Para: ");
        scanf("%d",&parada);
        if(parada == 0){
            break;
        }*/
    }
    //}while(m_Scf_2 != 677.583031);

    //Retorna ao vetor de solução, a melhor solução encontrada no GRASP
    for(m=0;m<ramos;m++){
        N[m].barra_1 = N_best[m].barra_1;
        N[m].barra_2 = N_best[m].barra_2;
        N[m].chave = N_best[m].chave;
        N[m].id = N_best[m].id;
        N[m].resist = N_best[m].resist;
        N[m].reat = N_best[m].reat;
        N[m].I_real = N_best[m].I_real;
        N[m].I_imag = N_best[m].I_imag;
    }

    fclose(arq5);
}

void VTree(struct Node *N, struct V_aux *vaux, int ramos, int subst, int barras, struct V_Node *T_node){

    int i, j, aux=0;
    bool Found;
    v_aux v1, v2;

    //-----------------------------------------------------------------------------------------------------------//
    //Construindo o vetor dinâmico que irá conter os índices das barras, menos o da subestação
    for(i=1;i<=barras;i++){
        if(i == subst){
            continue;
        }else{
            if(vaux == NULL){
                v1 = (v_aux)malloc(sizeof(V_aux));
                v1->v = i;
                v1->prox = NULL;
                vaux = v1;
            }else{
                v2 = (v_aux)malloc(sizeof(V_aux));
                v2->v = i;
                v2->prox = NULL;
                v1->prox = v2;
                v1 = v1->prox;
            }
        }
    }
    //-----------------------------------------------------------------------------------------------------------//

    /*v1 = vaux;
    printf("\nAntes de Tree.\n");
    while(v1!=NULL){
        printf(" %d",v1->v);
        v1 = v1->prox;
    }*/


    //Constrói a árvore do sistema elétrico
    Tree(N,ramos,barras,T_node,vaux);

    //Apagando o vetor dinâmico vaux - retornando seu ponteiro inicial para NULL
    v1 = vaux;
    while(vaux != NULL){
        if(v1->prox == NULL){
            vaux = NULL;
            free(v1);
        }else{
            v2 = v1;
            vaux->prox = v1->prox;
            v1 = v1->prox;
            free(v2);
        }
    }

    /*v1 = vaux;
    printf("\nDepois de Tree.\n");
    while(v1!=NULL){
        printf(" %d",v1->v);
        v1 = v1->prox;
    }*/
}

void Tree(struct Node *N, int ramos, int barras, struct V_Node *T_node, struct V_aux *vaux){
    int i, j, v;
    bool Found = false;
    v_aux v1, v2;

    /*v1 = vaux;
    printf("\nAntes da Recursao de Tree.\n");
    while(v1!=NULL){
        printf(" %d",v1->v);
        v1 = v1->prox;
    }
    system("pause");*/

    v1 = vaux;
    v2 = vaux;

    while((vaux!=NULL)&&(v1!=NULL)){
        Found = false;
        for(i=0;i<ramos;i++){
            if((((T_node->v == N[i].barra_1)&&(v1->v == N[i].barra_2))||((T_node->v == N[i].barra_2)&&(v1->v == N[i].barra_1)))&&(N[i].chave == 1)){
                Found =  true;

                //Adicionando uma posição ao contador de ponteiros do elemento corrente da árvore
                T_node->cont_p = T_node->cont_p + 1;

                //Verificando em que N.barra se encontra o índice encontrado
                if(N[i].barra_1 == T_node->v){
                    v = N[i].barra_2;
                }else{
                    v = N[i].barra_1;
                }

                //Apagando o elemento do vetor de índices de barra, cuja barra foi encontrada
                //Lê-se o vetor auxiliar de índides de barras SEMPRE do início
                while((v2!=v1)&&(v2->prox!=v1)){
                    v2=v2->prox;
                }
                 if(v2 == v1){
                    vaux = v1->prox;
                }else{
                    v2->prox = v1->prox;
                }
                free(v1);
                break;
            }
        }
        if(Found){
            VN v_novo;
            v_novo = (VN)malloc(sizeof(V_Node));
            v_novo->v = v;
            v_novo->cont_p = 0;
            for(j=0;j<LINK;j++){
                v_novo->p[j] = NULL;
            }
            T_node->p[(T_node->cont_p)-1] = v_novo;
            Tree(N,ramos,barras,T_node->p[(T_node->cont_p)-1],vaux);
            v1 = vaux;
            v2 = vaux;
            //return;
        }else{
            v1 = v1->prox;
        }

    }
}

//Percurso Pré-Ordem
void printTree_pre(struct V_Node *T_node){
    int i;
    if (T_node != NULL) {
        printf("<%d",T_node->v);
        for (i = 0; i < LINK; i++) {
            printTree_pre(T_node->p[i]);
        }
        printf(">");
    }
}

void printTree_pont(struct V_Node *T_node){
    int i;
    for (i = 0; i < LINK; i++) {
        if(T_node->p[i] != NULL){
            printf("\nBarra %d --- P[%d] -> %d",T_node->v,i,T_node->p[i]->v);
            printTree_pont(T_node->p[i]);
        }else{
            printf("\nBarra %d --- P[%d] -> NULL",T_node->v,i);
            break;
        }
    }
}

//Percurso Pós-Ordem
void printTree_pos(struct V_Node *T_node){
    int i;
    printf("<");
    for (i = LINK-1; i >= 0; i--) {
         if (T_node->p[i] == NULL) {
             continue;
         }else{

             printTree_pos(T_node->p[i]);

         }
    }
    printf("%d>",T_node->v);
}

void Apagar_Tree(struct V_Node *T_node){
    int i, v1, v2;
    for (i = 0; i < LINK; i++){
        if(T_node->p[i] != NULL){
            Apagar_Tree(T_node->p[i]);
            free(T_node->p[i]);
            T_node->p[i] = NULL;
        }else{
            return;
        }
    }
}

float Tree_pos_corrente_real(struct V_Node *T_node, struct Node *N, struct Tensao_Carga *TC, int ramos){
    int i, j, v1;
    float i_real, soma = 0.0;

    for(i=0;i<LINK;i++){
        if(T_node->p[i] == NULL){
            //printf("\n%d - %d - pt %d",T_node->v,v1,i);
            return 0.0;
        }else{
            v1 = T_node->p[i]->v;
            //printf("\n%d - %d - pt %d",T_node->v,v1,i);
            for(j=0;j<ramos;j++){
                if((N[j].barra_1 == T_node->v)&&(N[j].barra_2 == v1)||(N[j].barra_1 == v1)&&(N[j].barra_2 == T_node->v)){
                    break;
                }
            }
            i_real = (float) (((TC[v1-1].pot_at * TC[v1-1].V_real) + (TC[v1-1].pot_reat * TC[v1-1].V_imag)) / (pow(TC[v1-1].V_real,2) + pow(TC[v1-1].V_imag,2)));
            //printf("\n\nBarra %d - Pts: %f %f",v1,TC[v1-1].pot_at,TC[v1-1].pot_reat);
            //printf("\nRamo %d - %d --- i_real: %f",T_node->v,v1,i_real);
            N[j].I_real = (float) (N[j].I_real + i_real);
            N[j].I_real = (float) (N[j].I_real + Tree_pos_corrente_real(T_node->p[i],N,TC,ramos));
            soma = soma + N[j].I_real;
            if(T_node->p[i+1] == NULL){
                return soma;
            }
        }
    }
}

float Tree_pos_corrente_imag(struct V_Node *T_node, struct Node *N, struct Tensao_Carga *TC, int ramos){
    int i, j, v1;
    float i_imag, soma = 0.0;

    for(i=0;i<LINK;i++){
        if(T_node->p[i] == NULL){
            //printf("\n%d - %d - pt %d",T_node->v,v1,i);
            return 0.0;
        }else{
            v1 = T_node->p[i]->v;
            //printf("\n%d - %d - pt %d",T_node->v,v1,i);
            for(j=0;j<ramos;j++){
                if((N[j].barra_1 == T_node->v)&&(N[j].barra_2 == v1)||(N[j].barra_1 == v1)&&(N[j].barra_2 == T_node->v)){
                    break;
                }
            }
            i_imag = (float) (((TC[v1-1].pot_at * TC[v1-1].V_imag) - (TC[v1-1].pot_reat * TC[v1-1].V_real)) / (pow(TC[v1-1].V_real,2) + pow(TC[v1-1].V_imag,2)));
            N[j].I_imag = (float) (N[j].I_imag + i_imag);
            N[j].I_imag = (float) (N[j].I_imag + Tree_pos_corrente_imag(T_node->p[i],N,TC,ramos));
            soma = soma + N[j].I_imag;
            if(T_node->p[i+1] == NULL){
                return soma;
            }
        }
    }
}

float Perdas_Ativas(struct V_Node *T_node, struct Node *N, int ramos, int barras){
    int i, j, v1;
    float p_at = 0.0;

    for(i=0;i<ramos;i++){
        if(N[i].chave == 1){
            p_at = p_at + (N[i].resist * (pow(N[i].I_real,2) + pow(N[i].I_imag,2)));
        }
    }
    //printf("\nPerdas ativas: %f",p_at*10000);
    return(p_at);
}

float Perdas_Reativas(struct V_Node *T_node, struct Node *N, int ramos, int barras){
    int i, j, v1;
    float p_reat = 0.0;

    for(i=0;i<ramos;i++){
        if(N[i].chave == 1){
            p_reat = p_reat + (N[i].reat * (pow(N[i].I_real,2) + pow(N[i].I_imag,2)));
        }
    }
    //printf("\nPerdas ativas: %f",p_reat*10000);
    return(p_reat);
}

void Tree_pre_tensao_real(struct V_Node *T_node, struct Node *N, struct Tensao_Carga *TC, int ramos){
    int i, j, k, v1;

    for(i=0;i<LINK;i++){
        if(T_node->p[i] == NULL){
            return;
        }else{
            v1 = T_node->p[i]->v;
            for(j=0;j<ramos;j++){
                if((N[j].barra_1 == T_node->v)&&(N[j].barra_2 == v1)||(N[j].barra_1 == v1)&&(N[j].barra_2 == T_node->v)){
                    break;
                }
            }
            TC[v1-1].V_real = TC[T_node->v - 1].V_real -(N[j].resist * N[j].I_real) +(N[j].reat * N[j].I_imag);
            Tree_pre_tensao_real(T_node->p[i],N,TC,ramos);
        }
    }
}

void Tree_pre_tensao_imag(struct V_Node *T_node, struct Node *N, struct Tensao_Carga *TC, int ramos){
    int i, j, k, v1;

    for(i=0;i<LINK;i++){
        if(T_node->p[i] == NULL){
            return;
        }else{
            v1 = T_node->p[i]->v;
            for(j=0;j<ramos;j++){
                if((N[j].barra_1 == T_node->v)&&(N[j].barra_2 == v1)||(N[j].barra_1 == v1)&&(N[j].barra_2 == T_node->v)){
                    break;
                }
            }
            TC[v1-1].V_imag = TC[T_node->v - 1].V_imag -(N[j].resist * N[j].I_imag) -(N[j].reat * N[j].I_real);
            Tree_pre_tensao_imag(T_node->p[i],N,TC,ramos);
        }
    }
}

void Ciclo_v1(struct V_Node *T_node, struct Node *N, int ramos, int v1, int v2, int subst, int vet_ciclos_aux[][VETMAX], int *i_aux, int *j_aux){
    //Buscando a lista de ramos do ciclo v1-v2, partindo de v1
    int i;
    for(i=0;i<LINK;i++){
        if(T_node->v != v1){
            if(*i_aux == 0){
                if(T_node->p[i] == NULL){
                    return;
                }else{
                    Ciclo_v1(T_node->p[i],N,ramos,v1,v2,subst,vet_ciclos_aux,i_aux,j_aux);
                }
            }else{
                vet_ciclos_aux[0][*i_aux] = T_node->v;
                //printf("%d ",T_node->v);
                *i_aux = *i_aux + 1;
                return;
            }
        }else{
            vet_ciclos_aux[0][*i_aux] = T_node->v;
            //printf("%d ",T_node->v);
            *i_aux = *i_aux + 1;
            return;
        }
    }
}

void Ciclo_v2(struct V_Node *T_node, struct Node *N, int ramos, int v1, int v2, int subst, int vet_ciclos_aux[][VETMAX], int *i_aux, int *j_aux){
    //Buscando a lista de ramos do ciclo v1-v2, partindo de v2
    int i;
    for(i=0;i<LINK;i++){
        if(T_node->v != v2){
            if(*j_aux == 0){
                if(T_node->p[i] == NULL){
                    return;
                }else{
                    Ciclo_v2(T_node->p[i],N,ramos,v1,v2,subst,vet_ciclos_aux,i_aux,j_aux);
                }
            }else{
                vet_ciclos_aux[1][*j_aux] = T_node->v;
                //printf("%d ",T_node->v);
                *j_aux = *j_aux + 1;
                return;
            }
        }else{
            vet_ciclos_aux[1][*j_aux] = T_node->v;
            //printf("%d ",T_node->v);
            *j_aux = *j_aux + 1;
            return;
        }
    }
}


void FCR(struct V_Node *T_node, struct Node *N, struct Tensao_Carga *TC, int barras, int ramos,
    float *p_at_2, float *p_reat_2, float *delta_P, float *P_per_1, float *P_per_2, float lim, int v_ini){
    int cont=0, i;

	for(i=0;i<barras;i++){
        TC[i].V_real = v_ini;
        TC[i].V_imag = 0.0;
    }

	while((*delta_P > lim)&&(cont < 20)){
	//while(*delta_P > lim){
        cont++;

        //printf("\nITERACAO %d\n",cont);

        //Começa a calcular as tensões APÓS a primeira iteração
        if(cont > 1){
            Tree_pre_tensao_real(T_node,N,TC,ramos);
            Tree_pre_tensao_imag(T_node,N,TC,ramos);
        }

        /*for(int j=0;j<barras;j++){
            printf("\nBarra %d --- Tensão: %f",j,TC[j].V_real);
        }
        system("pause");*/

        for(i=0;i<ramos;i++){
            N[i].I_real = 0.0;
            N[i].I_imag = 0.0;
        }

        Tree_pos_corrente_real(T_node,N,TC,ramos);
        Tree_pos_corrente_imag(T_node,N,TC,ramos);

        /*printf("\nIt. %d",cont);
        for(i=0;i<ramos;i++){
            if(N[i].chave == 1)
            printf("\nRamo %d-%d --- Cor: %f %f",N[i].barra_1,N[i].barra_2,N[i].I_real,N[i].I_imag);
        }
        system("pause");*/

        /*for(i=0;i<barras;i++){
            TC[i].V_real = v_ini;
            TC[i].V_imag = 0.0;
        }*/

        *p_at_2 = Perdas_Ativas(T_node,N,ramos,barras);
        *p_reat_2 = Perdas_Reativas(T_node,N,ramos,barras);
        //printf("\nPerdas Ativas: %f\n",p_at*p_base);
        //printf("\nPerdas Reativas: %f\n",p_reat*p_base);

        *P_per_2 = *p_at_2;
        *delta_P = fabs((*P_per_2) - (*P_per_1));
        if((*delta_P) > lim){
            (*P_per_1) = (*P_per_2);
        }
        //printf("\nIteracao %d.",cont);
        //printf("\nDelta %f.",*delta_P);
        //printf("\nPerdas Ativas: %f",(*p_at_2)*1000);
        //printf("\nPerdas Reativas: %f\n",(*p_reat_2)*1000);
        //system("pause");
        //system("cls");
	}
	//printf("\nIteracoes necessarias: %d",cont);
	//printf("\nPerdas Ativas: %f",(*p_at_2)*1000);
    //printf("\nPerdas Reativas: %f\n",(*p_reat_2)*1000);
	//system("pause");
}

void Busca_Ciclos(struct V_Node *T_node, struct Node *N, int ramos, int barras, int subst, int **M_ciclos, int vet_ciclos[]){
    int i, j, aux4, v1, v2, k;
    bool parada;

    for(i=0;i<(ramos-(barras-1));i++){
        for(j=0;j<ramos;j++){
            M_ciclos[i][j] = 0;
            vet_ciclos[j] = 0;
        }
    }
    //-----------------------------------------------------------------------------------------------------------//

    //Busca os subciclos da solução armazenada em N para a árvore em T_node
    //Preenche a Matriz de subciclos, e os imprime na tela
    //Cada linha da Matriz se subciclos possui os índices das arestas pertencentes à cada subciclo encontrado
    //-----------------------------------------------------------------------------------------------------------//
    aux4 = 0;
    int vet_ciclos_aux[2][VETMAX], i_aux, j_aux;
    while(aux4 < ramos-(barras-1)){

        /*for (int teste=0;teste<VETMAX;teste++){
            vet_ciclos_aux[0][teste] = -1;
            vet_ciclos_aux[1][teste] = -1;
        }*/

        for(j=0;j<ramos;j++){
            i_aux = 0;
            j_aux = 0;
            if(N[j].chave == 0){

                //system("cls");
                v1 = N[j].barra_1;
                v2 = N[j].barra_2;
                //printf("\nInfo:");
                //printf("\nRamo: %d - Vertices: %d - %d\n",j,N[j].barra_1,N[j].barra_2);

                //Construção da submatriz que contém o ciclo correspondente à v1-v2
                //printf("v1-v2: %d-%d\n",v1,v2);
                int i_teste, j_teste;
                Ciclo_v1(T_node,N,ramos,v1,v2,subst,vet_ciclos_aux,&i_aux,&j_aux);
                Ciclo_v2(T_node,N,ramos,v1,v2,subst,vet_ciclos_aux,&i_aux,&j_aux);


                /*printf("\nÍndices do ramo desligado: %d - %d",v1,v2);
                printf("\n");
                for (int teste=0;teste<VETMAX;teste++){
                    printf("- %d -",vet_ciclos_aux[0][teste]);
                }
                printf("\n");
                for (int teste=0;teste<VETMAX;teste++){
                    printf("- %d -",vet_ciclos_aux[1][teste]);
                }
                system("pause");*/

                /*printf("\nCiclos_v1\n");
                for(int i1=0;i1<i_aux;i1++){
                    printf("%d ",vet_ciclos_aux[0][i1]);
                }
                printf("\nCiclos_v2\n");
                for(int i1=0;i1<j_aux;i1++){
                    printf("%d ",vet_ciclos_aux[1][i1]);
                }*/

                /*printf("i_aux: %d --- j_aux: %d\n",i_aux,j_aux);
                for(i_teste=0;i_teste<i_aux;i_teste++){
                    printf("%d ",vet_ciclos_aux[0][i_teste]);
                }
                printf("\n");
                for(j_teste=0;j_teste<j_aux;j_teste++){
                    printf("%d ",vet_ciclos_aux[1][j_teste]);
                }*/
                //system("pause");

                //Caminha 2 indicadores no submatriz de ciclo a fim de encontrar o primeiro elemento igual às duas linhas da submatriz
                parada = false;
                int i_vet = 0, j_vet = 0;
                for(i_vet=0;i_vet<i_aux;i_vet++){
                    for(j_vet=0;j_vet<j_aux;j_vet++){
                        if(vet_ciclos_aux[0][i_vet] == vet_ciclos_aux[1][j_vet]){
                            i_teste = i_vet;
                            j_teste = j_vet;
                            parada = true;
                            break;
                        }
                    }
                    if(parada) break;
                }
                //printf("\ni_teste - j_teste: %d - %d",i_teste,j_teste);

                /*printf("\nCiclos_v1 - Corrigido\n");
                for(int i1=0;i1<=i_teste;i1++){
                    printf("%d ",vet_ciclos_aux[0][i1]);
                }
                printf("\nCiclos_v2 - Corrigido\n");
                for(int i1=0;i1<=j_teste;i1++){
                    printf("%d ",vet_ciclos_aux[1][i1]);
                }
                printf("\n");*/

                //Inserindo Ciclo_v1, linha 0 de vet_ciclos_aux, na M_ciclos
                int loop;
                int i_vet_1 = 0;
                if(i_teste != 0){ //Verificando linha com somente 1 elemento - raiz por exemplo
                    loop=0;
                    while (loop < i_teste){
                        for(k=0;k<ramos;k++){
                            if((N[k].barra_1 == vet_ciclos_aux[0][i_vet_1])&&(N[k].barra_2 == vet_ciclos_aux[0][i_vet_1 + 1])
                            ||(N[k].barra_1 == vet_ciclos_aux[0][i_vet_1 + 1])&&(N[k].barra_2 == vet_ciclos_aux[0][i_vet_1])){
                                M_ciclos[aux4][k] = 1;
                                i_vet_1++;
                                break;
                            }
                        }
                        loop++;
                    }
                }

                //Inserindo Ciclo_v2, linha 1 de vet_ciclos_aux, na M_ciclos
                int j_vet_1 = 0;
                if(j_teste != 0){ //Verificando linha com somente 1 elemento - raiz por exemplo
                    loop=0;
                    while (loop < j_teste){
                        for(k=0;k<ramos;k++){
                            if((N[k].barra_1 == vet_ciclos_aux[1][j_vet_1])&&(N[k].barra_2 == vet_ciclos_aux[1][j_vet_1 + 1])
                            ||(N[k].barra_1 == vet_ciclos_aux[1][j_vet_1 + 1])&&(N[k].barra_2 == vet_ciclos_aux[1][j_vet_1])){
                                M_ciclos[aux4][k] = 2;
                                j_vet_1++;
                                break;
                            }
                        }
                        loop++;
                    }
                }

                //Inserindo o ramo desligado na M_ciclos
                for(k=0;k<ramos;k++){
                    if((N[k].barra_1 == vet_ciclos_aux[0][0])&&(N[k].barra_2 == vet_ciclos_aux[1][0])
                    ||(N[k].barra_1 == vet_ciclos_aux[1][0])&&(N[k].barra_2 == vet_ciclos_aux[0][0])){
                        M_ciclos[aux4][k] = 3;
                    }
                }

                /*printf("\n");
                /*for(k=0;k<ramos;k++){
                    //system("cls");
                    //printf("\n%d ",vet_ciclos[k]);
                    M_ciclos[aux4][k] = vet_ciclos[k];
                    if(vet_ciclos[k] == 1){
                        //printf(" -- %d - %d",N[k].barra_1,N[k].barra_2);
                    }
                }
                M_ciclos[aux4][j] = 2;
                aux4++;
                //system("pause");
                for(k=0;k<ramos;k++){
                    vet_ciclos[k] = 0;
                }*/
                aux4++;
            }
        }
    }

    //Imprime na tela a Matriz de subciclos
    /*printf("\n");
    for(i=0;i<(ramos-(barras-1));i++){
        printf("\nBT iteracao: %d\n",bt_cont);
        for(j=0;j<ramos;j++){
            printf(" %d",M_ciclos[i][j]);
        }
        for(j=0;j<ramos;j++){
            if(M_ciclos[i][j] == 1){
                printf("\nRamo: %d -- Barras: %d - %d",j,N[j].barra_1,N[j].barra_2);
            }
            if(M_ciclos[i][j] == 2){
                printf("\nCiclo do ramo: %d -- Barras: %d - %d",j,N[j].barra_1,N[j].barra_2);
            }
        }
        printf("\n");
        system("pause");
        system("cls");
    }*/
}

void PR_1(struct V_Node *T_node, struct Node *N, struct Tensao_Carga *TC, struct V_aux *vaux, int ramos, int barras, int subst, int **M_ciclos, int vet_ciclos[],
    float *p_at_2, float *p_reat_2, float *delta_P, float *P_per_1,float *P_per_2, float lim, int v_ini, float p_base, float *p_at_top, float *p_reat_top,
    struct Node M[][16], int n_pr, float FC_ini[]){

    FILE *arq7, *arq8;
    arq7 = fopen("sols_pr.txt","w");
    arq8 = fopen("med_fluxo_pr.txt","w");

    bool factivel, similar;
    int inicio, fim, i, j, j1, k, elem, vet_pr[ramos], i_pr, i_dif, j_dif, i_desativo, i_ativo;
    int chave_guia, chave_int, iter_pr=1, cont, sols, sols_total, i_ciclo;
    sols_total = 0;
    //inicio = 0;
    //fim = 9;

    //Soluções guia e intermediária
    struct Node N_guia[ramos], N_int[ramos];
    float p_at_guia, p_at_int, p_reat_guia, p_reat_int, p_at, p_reat, soma_fluxo_pr_int;

    for(i=0;i<9;i++){
        for(j=0;j<ramos;j++){
            N_guia[j].barra_1 = M[i][j].barra_1;
            N_guia[j].barra_2 = M[i][j].barra_2;
            N_guia[j].chave = M[i][j].chave;
            N_guia[j].id = M[i][j].id;
            N_guia[j].resist = M[i][j].resist;
            N_guia[j].reat = M[i][j].reat;
            N_guia[j].I_real = M[i][j].I_real;
            N_guia[j].I_imag = M[i][j].I_imag;
        }
        /*Apagar_Tree(T_node);
        T_node->cont_p = 0;
        VTree(N_guia,vaux,ramos,subst,barras,T_node);
        *delta_P = 1.0;
        FCR(T_node,N_guia,TC,barras,ramos,p_at_2,p_reat_2,delta_P,P_per_1,P_per_2,lim,v_ini);
        printf("\nSolução %d de M_pr - Perdas Ativas: %f",i,(float) (*p_at_2)*p_base);*/
    }
    //system("pause");

    //Aplicação do PR em todos os pares de soluções encontradas
    //system("cls");
    //while(inicio < fim){
    for(inicio=0;inicio<n_pr-2;inicio++){
        for(fim=inicio+1;fim<n_pr-1;fim++){
        //system("cls");
        sols = 0;
        //printf("\nIteracao %d do PR - Sols. %d - %d",iter_pr,inicio,fim);

        //Calcula as perdas das funções para decidir qual solução é a guia e qual a intermediária
        //A solução escolhida como guia sempre terá as menores perdas elétricas

        //Calculando as perdas da solução da extremidade INICIO da matriz
        for(i=0;i<ramos;i++){
            N_guia[i].barra_1 = M[inicio][i].barra_1;
            N_guia[i].barra_2 = M[inicio][i].barra_2;
            N_guia[i].chave = M[inicio][i].chave;
            N_guia[i].id = M[inicio][i].id;
            N_guia[i].resist = M[inicio][i].resist;
            N_guia[i].reat = M[inicio][i].reat;
            N_guia[i].I_real = M[inicio][i].I_real;
            N_guia[i].I_imag = M[inicio][i].I_imag;
        }
        Apagar_Tree(T_node);
        T_node->cont_p = 0;
        VTree(N_guia,vaux,ramos,subst,barras,T_node);
        *delta_P = 1.0;
        FCR(T_node,N_guia,TC,barras,ramos,p_at_2,p_reat_2,delta_P,P_per_1,P_per_2,lim,v_ini);
        p_at = *p_at_2;

        //Calculando as perdas da solução da extremidade FIM da matriz
        for(i=0;i<ramos;i++){
            N_guia[i].barra_1 = M[fim][i].barra_1;
            N_guia[i].barra_2 = M[fim][i].barra_2;
            N_guia[i].chave = M[fim][i].chave;
            N_guia[i].id = M[fim][i].id;
            N_guia[i].resist = M[fim][i].resist;
            N_guia[i].reat = M[fim][i].reat;
            N_guia[i].I_real = M[fim][i].I_real;
            N_guia[i].I_imag = M[fim][i].I_imag;
        }
        Apagar_Tree(T_node);
        T_node->cont_p = 0;
        VTree(N_guia,vaux,ramos,subst,barras,T_node);
        *delta_P = 1.0;
        FCR(T_node,N_guia,TC,barras,ramos,p_at_2,p_reat_2,delta_P,P_per_1,P_per_2,lim,v_ini);

        //Escolhendo como guia a solução de menor perdas elétricas
        if((float) ((*p_at_2)*p_base) < (float) (p_at*p_base)){
            //Guia->fim ----- Intermediária->inicio
            for(i=0;i<ramos;i++){
                N_int[i].barra_1 = M[inicio][i].barra_1;
                N_int[i].barra_2 = M[inicio][i].barra_2;
                N_int[i].chave = M[inicio][i].chave;
                N_int[i].id = M[inicio][i].id;
                N_int[i].resist = M[inicio][i].resist;
                N_int[i].reat = M[inicio][i].reat;
                N_int[i].I_real = M[inicio][i].I_real;
                N_int[i].I_imag = M[inicio][i].I_imag;
            }
            p_at_guia = *p_at_2;
        }else{
            //Guia->inicio ----- Intermediária->fim
            for(i=0;i<ramos;i++){
                N_guia[i].barra_1 = M[inicio][i].barra_1;
                N_guia[i].barra_2 = M[inicio][i].barra_2;
                N_guia[i].chave = M[inicio][i].chave;
                N_guia[i].id = M[inicio][i].id;
                N_guia[i].resist = M[inicio][i].resist;
                N_guia[i].reat = M[inicio][i].reat;
                N_guia[i].I_real = M[inicio][i].I_real;
                N_guia[i].I_imag = M[inicio][i].I_imag;

                N_int[i].barra_1 = M[fim][i].barra_1;
                N_int[i].barra_2 = M[fim][i].barra_2;
                N_int[i].chave = M[fim][i].chave;
                N_int[i].id = M[fim][i].id;
                N_int[i].resist = M[fim][i].resist;
                N_int[i].reat = M[fim][i].reat;
                N_int[i].I_real = M[fim][i].I_real;
                N_int[i].I_imag = M[fim][i].I_imag;
            }
        }

        Apagar_Tree(T_node);
        T_node->cont_p = 0;
        VTree(N_guia,vaux,ramos,subst,barras,T_node);
        *delta_P = 1.0;
        FCR(T_node,N_guia,TC,barras,ramos,p_at_2,p_reat_2,delta_P,P_per_1,P_per_2,lim,v_ini);
        p_at_guia = *p_at_2;
        p_reat_guia = *p_reat_2;
        //printf("\nN_guia - Perdas: %f",p_at_guia*p_base);

        Apagar_Tree(T_node);
        T_node->cont_p = 0;
        VTree(N_int,vaux,ramos,subst,barras,T_node);
        *delta_P = 1.0;
        FCR(T_node,N_int,TC,barras,ramos,p_at_2,p_reat_2,delta_P,P_per_1,P_per_2,lim,v_ini);
        p_at_int = *p_at_2;
        p_reat_int = *p_reat_2;
        //printf("\nN_int - Perdas: %f",p_at_int*p_base);

        //printf("\n");
        //Iniciando o procedimento de Path Relinking
        do{
            similar = true; //Variável lógica que testa a igualdade entre as 2 soluções
            //i = 0;
            i_pr = 0;
            //Percorre ambos vetores de ramos das soluções guia e intermediária
            //Para cada característica diferente de chave, salva-se tal índice de ramo em um vetor e se altera a variável de similaridade
            //vet_pr é um vetor que armazena os índices dos ramos da solução intermediária cujas chaves são diferentes da guia
            for(i=0;i<ramos;i++){
                if(N_guia[i].chave != N_int[i].chave){
                    similar = false;
                    vet_pr[i_pr] = i;
                    i_pr++;
                }
            }

            if(similar == false){  //Testa de se a posição i parou dentro de N como posição válida ou se passou do limite e as soluções são idênticas

                //Quando as soluções são diferentes, se escolhe aleatoriamente um índice do vetor de índices vet_pr
                //para se gerar a primeira solução intermediária...o índice do ramo escolhido é armazenado em i
                srand((unsigned)time(NULL));
                Sleep(300);
                i_dif = rand()%i_pr;
                i = vet_pr[i_dif];

                factivel = true;
                chave_guia = N_guia[i].chave;
                chave_int = N_int[i].chave;
                elem = i; //Guardando a posição da característica de chave diferente entre as soluções guia/intermediaria

                /*cont=0;
                for(int p=0;p<ramos;p++){
                    if(N_int[p].chave == 0){
                        cont++;
                    }
                }
                printf("\nQtde de ramos abertos: %d",cont);*/

                Apagar_Tree(T_node);
                T_node->cont_p = 0;
                VTree(N_int,vaux,ramos,subst,barras,T_node);
                *delta_P = 1.0;
                FCR(T_node,N_int,TC,barras,ramos,p_at_2,p_reat_2,delta_P,P_per_1,P_per_2,lim,v_ini);
                p_at_int = *p_at_2;
                p_reat_int = *p_reat_2;
                //printf("\nN_int produzida - Perdas: %f",p_at_int*p_base);
                Busca_Ciclos(T_node,N_int,ramos,barras,subst,M_ciclos,vet_ciclos);


                //Avalia a factibilidade da solução para posterior comparação com a melhor solução encontrada até o momento
                for(int w=0;w<barras;w++){
                    if((TC[w].V_real < 0.93)||(TC[w].V_real > 1.05)){
                        factivel = false;
                        //printf("\tN_int produzida INF - Perdas: %f",p_at_int*p_base);
                        break;
                    }
                }

                soma_fluxo_pr_int = 0.0;
                for(int i_ciclo=0;i_ciclo<ramos;i_ciclo++){
                    if(N_int[i_ciclo].chave == 1){
                        soma_fluxo_pr_int = (float) soma_fluxo_pr_int + FC_ini[i_ciclo];
                    }
                }
                //printf("\nMedia de fluxo de carga da solucao inicial: %f\n", (float) soma_fluxo_pr_int/(barras-1));

                //Salvando no arquivo as perdas das soluções intermediarias
                fprintf (arq7,"%f\n",p_at_int*p_base);
                fprintf (arq8,"%f\n",(float) soma_fluxo_pr_int/(barras-1));

                //Caso a solução intermediaria encontrada seja factivel, verifica-se se esta é melhor do que
                //a melhor solução encontrada até o momento, armazenada em N_top, parâmetro passado para N nessa função
                if(factivel == true){
                    if((float) (p_at_int*p_base) < (float) ((*p_at_top)*p_base)){
                        for(int m=0;m<ramos;m++){
                            N[m].barra_1 = N_int[m].barra_1;
                            N[m].barra_2 = N_int[m].barra_2;
                            N[m].chave = N_int[m].chave;
                            N[m].id = N_int[m].id;
                            N[m].resist = N_int[m].resist;
                            N[m].reat = N_int[m].reat;
                            N[m].I_real = N_int[m].I_real;
                            N[m].I_imag = N_int[m].I_imag;
                        }
                        *p_at_top = p_at_int;
                        *p_reat_top = p_reat_int;
                        //printf("\nMelhor perdas em PR: %f",*p_at_top);
                        //printf("\tN_int produzida FAC - Perdas: %f",p_at_int*p_base);
                    }
                }

                //A geração de uma nova solução intermediaria segue normalmente não se atendendo à factibilidade da solução
                //Procura-se em cada linha da matriz de ciclos, na posição do ramo escolhido, a solução que possui a característica
                //de chave diferente entre as soluções -- pega-se a primeira linha que atende a este comportamento e se armazena
                //o índice da linha em j

                //Encontrar a que ciclo pertence o ramo da chave escolhida
                i_pr = 0;
                for(j=0;j<(ramos-(barras-1));j++){
                    if(chave_int == 1){ //Ramo aberto em N_guia procura um ramo fechado na linha de ciclo j da matriz de ciclos
                        if ((M_ciclos[j][elem] == 1)||(M_ciclos[j][elem] == 2)){ //Ramos fechados
                            vet_pr[i_pr] = j;
                            i_pr++;
                        }
                    }else{ //Ramo fechado em N_guia procura o ramo aberto na linha de ciclo j da matriz de ciclos
                        if (M_ciclos[j][elem] == 3){ //Ramo aberto
                            break;
                        }
                    }
                }

                if(chave_int == 1){
                    srand((unsigned)time(NULL));
                    Sleep(300);
                    j_dif = rand()%i_pr;
                    j = vet_pr[j_dif];
                }

                i_pr = 0;
                //Percorrendo a linha de índice j encontrado anteriormente a fim de encontrar os índices de ramos
                //restantes pertencentes ao ramo de índice elem, e completar tais índices em um vetor
                if(chave_int == 0){
                    for(k=0;k<ramos;k++){
                        if(k!=elem){ //PULA o índice do ramo encontrado
                            //Salva em vet_pr somente os ramos com chave fechada
                            if ((M_ciclos[j][k] == 1)||(M_ciclos[j][k] == 2)){
                                vet_pr[i_pr] = k;
                                i_pr++;
                            }
                        }
                    }
                }else{
                    for(k=0;k<ramos;k++){
                        if(k!=elem){ //PULA o índice do ramo encontrado
                            //Salva em vet_pr somente os ramos com chave fechada
                            if (M_ciclos[j][k] == 3){
                                i_desativo = k;
                                break;
                            }
                        }
                    }
                }

                //Trocando a chave do ramo de N_int que diferenciou de N_guia
                if(N_int[elem].chave == 1){
                //if(chave_int == 1){
                    //Quando desligo um novo ramo do ciclo, o ramo desligado ORIGINAL se ativa
                    N_int[elem].chave = 0;
                    N_int[i_desativo].chave = 1;
                }else{
                    //Quando se ativa o ramo desligado, se escolhe aleatoriamente um ramo dentre os ramos ativos, e o desativa
                    N_int[elem].chave = 1;
                    srand((unsigned)time(NULL));
                    Sleep(300);
                    i_ativo = rand()%i_pr;
                    N_int[vet_pr[i_ativo]].chave = 0;
                }
                sols++;
            }
            //system("pause");
        }while(similar == false);
        //printf("\n%d solucoes geradas na iteracao %d do PR_1!\n",sols,iter_pr);
        sols_total = sols_total + sols;
        //inicio++;
        //fim--;
        iter_pr++;
        //system("pause");
        }
    }
    printf("\nTotal de solucoes geradas no PR_1: %d!\n",sols_total);
    fclose(arq7);
    fclose(arq8);
}

void PR_2(struct V_Node *T_node, struct Node *N, struct Tensao_Carga *TC, struct V_aux *vaux, int ramos, int barras, int subst, int **M_ciclos, int vet_ciclos[],
    float *p_at_2, float *p_reat_2, float *delta_P, float *P_per_1,float *P_per_2, float lim, int v_ini, float p_base, float *p_at_top, float *p_reat_top,
    struct Node M[][16], int n_pr){

    bool factivel, similar, escolha_1, escolha_2, similar_2,similar_3;
    int inicio, fim, i, j, j1, k, elem, vet_pr[ramos], i_pr, i_dif, j_dif, i_desativo, i_ativo, viz;
    int chave_guia, chave_int, iter_pr=1, cont, sols, it1, it2, ramos_fechados, status, ramos_fechados_aux, sols_total;
    sols_total = 0;
    //inicio = 0;
    //fim = 7;

    //Soluções guia e intermediária
    struct Node N_guia[ramos], N_int[ramos], N_aux[ramos];
    float p_at_guia, p_at_int, p_reat_guia, p_reat_int, p_at, p_reat, p_at_base, p_reat_base;

    for(i=0;i<9;i++){
        for(j=0;j<ramos;j++){
            N_guia[j].barra_1 = M[i][j].barra_1;
            N_guia[j].barra_2 = M[i][j].barra_2;
            N_guia[j].chave = M[i][j].chave;
            N_guia[j].id = M[i][j].id;
            N_guia[j].resist = M[i][j].resist;
            N_guia[j].reat = M[i][j].reat;
            N_guia[j].I_real = M[i][j].I_real;
            N_guia[j].I_imag = M[i][j].I_imag;
        }
        Apagar_Tree(T_node);
        T_node->cont_p = 0;
        VTree(N_guia,vaux,ramos,subst,barras,T_node);
        *delta_P = 1.0;
        FCR(T_node,N_guia,TC,barras,ramos,p_at_2,p_reat_2,delta_P,P_per_1,P_per_2,lim,v_ini);
        printf("\nSolução %d de M_pr - Perdas Ativas: %f",i,(float) (*p_at_2)*p_base);
    }
    system("pause");

    //Utiliza as soluções em pares, partindo das extremidades da matriz para o centro
    system("cls");
    //while(inicio < fim){
    for(it1=0;it1<n_pr-2;it1++){
        for(it2=it1+1;it2<n_pr-1;it2++){
            //system("cls");
            sols = 0;
            //printf("\nIteracao %d do PR",iter_pr);
            printf("\nIteracao %d do PR - Sols. %d - %d",iter_pr,it1,it2);

            //N_guia com as melhores perdas -- N_int com as piores perdas
            for(i=0;i<ramos;i++){
                N_guia[i].barra_1 = M[it2][i].barra_1;
                N_guia[i].barra_2 = M[it2][i].barra_2;
                N_guia[i].chave = M[it2][i].chave;
                N_guia[i].id = M[it2][i].id;
                N_guia[i].resist = M[it2][i].resist;
                N_guia[i].reat = M[it2][i].reat;
                N_guia[i].I_real = M[it2][i].I_real;
                N_guia[i].I_imag = M[it2][i].I_imag;

                N_int[i].barra_1 = M[it1][i].barra_1;
                N_int[i].barra_2 = M[it1][i].barra_2;
                N_int[i].chave = M[it1][i].chave;
                N_int[i].id = M[it1][i].id;
                N_int[i].resist = M[it1][i].resist;
                N_int[i].reat = M[it1][i].reat;
                N_int[i].I_real = M[it1][i].I_real;
                N_int[i].I_imag = M[it1][i].I_imag;
            }

            Apagar_Tree(T_node);
            T_node->cont_p = 0;
            VTree(N_guia,vaux,ramos,subst,barras,T_node);
            *delta_P = 1.0;
            FCR(T_node,N_guia,TC,barras,ramos,p_at_2,p_reat_2,delta_P,P_per_1,P_per_2,lim,v_ini);
            p_at_guia = *p_at_2;
            p_reat_guia = *p_reat_2;
            //printf("\nP_at_guia: %f",*p_at_2);
            printf("\nN_guia - Perdas: %f",(float) p_at_guia*p_base);

            Apagar_Tree(T_node);
            T_node->cont_p = 0;
            VTree(N_int,vaux,ramos,subst,barras,T_node);
            *delta_P = 1.0;
            FCR(T_node,N_int,TC,barras,ramos,p_at_2,p_reat_2,delta_P,P_per_1,P_per_2,lim,v_ini);
            p_at_base = *p_at_2;
            p_reat_base = *p_reat_2;
            //printf("\nP_at_int: %f",*p_at_2);
            printf("\nN_base - Perdas: %f",(float) p_at_base*p_base);

            printf("\n");
            //Iniciando o procedimento de Path Relinking
            do{
                //printf("\nEntrou DO\n");
            //system("pause");
                similar = true; //Variável lógica que testa a igualdade entre as 2 soluções
                //i = 0;
                i_pr = 0;
                ramos_fechados = 0;
                //Percorre ambos vetores de ramos das soluções guia e intermediária
                //Para cada característica diferente de chave, salva-se tal índice de ramo em um vetor e se altera a variável de similaridade
                //vet_pr é um vetor que armazena os índices dos ramos da solução intermediária cujas chaves são diferentes da guia
                for(i=0;i<ramos;i++){
                    if((N_guia[i].chave == 1)&&(N_int[i].chave == 1)){
                        ramos_fechados++;
                    }
                    if((N_guia[i].chave == 1)&&(N_int[i].chave == 0)){
                        similar = false;
                        vet_pr[i_pr] = i;
                        i_pr++;
                    }
                }

                /*printf("\n\nN_guia e N_int\n");
                for(int m=0;m<ramos;m++){
                    printf("%d ",N_guia[m].chave);
                }
                printf("\n");
                for(int m=0;m<ramos;m++){
                    printf("%d ",N_int[m].chave);
                }
                system("pause");*/

                 similar_3 = true;
                            for(i=0;i<ramos;i++){
                                if((N_guia[i].chave != N_int[i].chave)){
                                    similar_3 = false;
                                    break;
                                }
                            }

                if((similar == false) && (similar_3!=true)) {
                        //printf("\nEntrou IF\n");
                        //system("pause");//Testa de se a posição i parou dentro de N como posição válida ou se passou do limite e as soluções são idênticas

                    //Quando as soluções são diferentes, se escolhe aleatoriamente um índice do vetor de índices vet_pr
                    //para se gerar a primeira solução intermediária...o índice do ramo escolhido é armazenado em i
                    srand((unsigned)time(NULL));
                    Sleep(300);
                    i_dif = rand()%i_pr;
                    i = vet_pr[i_dif];

                    chave_guia = N_guia[i].chave;
                    chave_int = N_int[i].chave;
                    elem = i; //Guardando a posição da característica de chave diferente entre as soluções guia/intermediaria

                    Apagar_Tree(T_node);
                    T_node->cont_p = 0;
                    VTree(N_int,vaux,ramos,subst,barras,T_node);
                    *delta_P = 1.0;
                    FCR(T_node,N_int,TC,barras,ramos,p_at_2,p_reat_2,delta_P,P_per_1,P_per_2,lim,v_ini);
                    //p_at_int = *p_at_2;
                    //p_reat_int = *p_reat_2;
                    Busca_Ciclos(T_node,N_int,ramos,barras,subst,M_ciclos,vet_ciclos);

                    //Encontrar a que ciclo pertence o ramo da chave escolhida
                    for(j=0;j<(ramos-(barras-1));j++){
                        /*if(chave_int == 1){ //Ramo aberto em N_guia procura um ramo fechado na linha de ciclo j da matriz de ciclos
                            if ((M_ciclos[j][elem] == 1)||(M_ciclos[j][elem] == 2)){ //Ramos fechados
                                vet_pr[i_pr] = j;
                                i_pr++;
                            }
                        }else{*/ //Ramo fechado em N_guia procura o ramo aberto na linha de ciclo j da matriz de ciclos
                            if (M_ciclos[j][elem] == 3){ //Ramo aberto
                                break;
                            }
                        //}
                    }
                              escolha_1 = false;
                            escolha_2 = false;
                    //Avaliação da vizinhança de PR do ramo escolhido
                    for(k=0;k<ramos;k++){

                        if((M_ciclos[j][k] == 1)||(M_ciclos[j][k] == 2)){
                            N_int[k].chave = 0;
                            N_int[elem].chave = 1;

                            Apagar_Tree(T_node);
                            T_node->cont_p = 0;
                            VTree(N_int,vaux,ramos,subst,barras,T_node);
                            *delta_P = 1.0;
                            FCR(T_node,N_int,TC,barras,ramos,p_at_2,p_reat_2,delta_P,P_per_1,P_per_2,lim,v_ini);
                            p_at_int = *p_at_2;
                            p_reat_int = *p_reat_2;
                            factivel = true;
                            sols++;

                            //Avalia a factibilidade da solução para posterior comparação com a melhor solução encontrada até o momento
                            for(int w=0;w<barras;w++){
                                if((TC[w].V_real < 0.93)||(TC[w].V_real > 1.05)){
                                    factivel = false;
                                    //printf("\tN_int produzida INF - Perdas: %f",p_at_int*p_base);
                                    break;
                                }
                            }

                            //printf("\nSolucao %d da iteracao %d -- Factibilidade (%d) -- Perdas: %f",sols,iter_pr,factivel,(float) p_at_int*p_base);
                            /*if (sols == 6){

                                printf("\n\nN_guia e N_int\n");
                                for(int m=0;m<ramos;m++){
                                    printf("%d ",N_guia[m].chave);
                                }
                                printf("\n");
                                for(int m=0;m<ramos;m++){
                                    printf("%d ",N_int[m].chave);
                                }
                                system("pause");

                            }*/

                            //Caso a solução intermediaria encontrada seja factivel, verifica-se se esta é melhor do que
                            //a melhor solução encontrada até o momento, armazenada em N_top, parâmetro passado para N nessa função
                            if(factivel == true){
                                if((float) (p_at_int*p_base) < (float) ((*p_at_top)*p_base)){
                                    for(int m=0;m<ramos;m++){
                                        N[m].barra_1 = N_int[m].barra_1;
                                        N[m].barra_2 = N_int[m].barra_2;
                                        N[m].chave = N_int[m].chave;
                                        N[m].id = N_int[m].id;
                                        N[m].resist = N_int[m].resist;
                                        N[m].reat = N_int[m].reat;
                                        N[m].I_real = N_int[m].I_real;
                                        N[m].I_imag = N_int[m].I_imag;
                                    }
                                    *p_at_top = p_at_int;
                                    *p_reat_top = p_reat_int;
                                    //printf("\nMelhor perdas em PR: %f",*p_at_top);
                                    //printf("\tN_int produzida FAC - Perdas: %f",p_at_int*p_base);
                                }
                            }

                            similar_2 = true;
                            for(int i2=0;i2<ramos;i2++){
                                if((N_guia[i2].chave != N_int[i2].chave)){
                                    similar_2 = false;
                                    break;
                                }
                            }

                            if(!similar_2){
                                ramos_fechados_aux = 0;
                                for(int i1=0;i1<ramos;i1++){
                                    if((N_guia[i1].chave == 1)&&(N_int[i1].chave == 1)){
                                        ramos_fechados_aux ++;
                                    }
                                }

                               // escolha_1 = false;
                              //  escolha_2 = false;
                                if(ramos_fechados_aux > ramos_fechados){
                                    escolha_1 = true;
                                    //if(ramos_fechados_aux == (ramos_fechados+1)){
                                    for(int m=0;m<ramos;m++){
                                        N_aux[m].barra_1 = N_int[m].barra_1;
                                        N_aux[m].barra_2 = N_int[m].barra_2;
                                        N_aux[m].chave = N_int[m].chave;
                                        N_aux[m].id = N_int[m].id;
                                        N_aux[m].resist = N_int[m].resist;
                                        N_aux[m].reat = N_int[m].reat;
                                        N_aux[m].I_real = N_int[m].I_real;
                                        N_aux[m].I_imag = N_int[m].I_imag;
                                    }
                                    //}
                                    ramos_fechados = ramos_fechados_aux;
                                }

                                if(!escolha_1){
                                    if(ramos_fechados_aux == ramos_fechados){
                                        escolha_2 = true;
                                        //if(ramos_fechados_aux == (ramos_fechados+1)){
                                            for(int m=0;m<ramos;m++){
                                                N_aux[m].barra_1 = N_int[m].barra_1;
                                                N_aux[m].barra_2 = N_int[m].barra_2;
                                                N_aux[m].chave = N_int[m].chave;
                                                N_aux[m].id = N_int[m].id;
                                                N_aux[m].resist = N_int[m].resist;
                                                N_aux[m].reat = N_int[m].reat;
                                                N_aux[m].I_real = N_int[m].I_real;
                                                N_aux[m].I_imag = N_int[m].I_imag;
                                            }
                                        //}
                                        //ramos_fechados = ramos_fechados_aux;
                                    }
                                }

                                if(!escolha_1 && !escolha_2){
                                    if(ramos_fechados_aux < ramos_fechados){
                                        //if(ramos_fechados_aux == (ramos_fechados+1)){
                                            for(int m=0;m<ramos;m++){
                                                N_aux[m].barra_1 = N_int[m].barra_1;
                                                N_aux[m].barra_2 = N_int[m].barra_2;
                                                N_aux[m].chave = N_int[m].chave;
                                                N_aux[m].id = N_int[m].id;
                                                N_aux[m].resist = N_int[m].resist;
                                                N_aux[m].reat = N_int[m].reat;
                                                N_aux[m].I_real = N_int[m].I_real;
                                                N_aux[m].I_imag = N_int[m].I_imag;
                                            }
                                        //}
                                        ramos_fechados = ramos_fechados_aux;
                                    }
                                }
                            } else{
                                for(int m=0;m<ramos;m++){
                                                N_aux[m].barra_1 = N_int[m].barra_1;
                                                N_aux[m].barra_2 = N_int[m].barra_2;
                                                N_aux[m].chave = N_int[m].chave;
                                                N_aux[m].id = N_int[m].id;
                                                N_aux[m].resist = N_int[m].resist;
                                                N_aux[m].reat = N_int[m].reat;
                                                N_aux[m].I_real = N_int[m].I_real;
                                                N_aux[m].I_imag = N_int[m].I_imag;
                                            }
                            break;}

                            //Retornando os status de chaves iniciais
                            N_int[k].chave = 1;
                            N_int[elem].chave = 0;
                        }


}                   // }
                    for(int m=0;m<ramos;m++){
                        N_int[m].barra_1 = N_aux[m].barra_1;
                        N_int[m].barra_2 = N_aux[m].barra_2;
                        N_int[m].chave = N_aux[m].chave;
                        N_int[m].id = N_aux[m].id;
                        N_int[m].resist = N_aux[m].resist;
                        N_int[m].reat = N_aux[m].reat;
                        N_int[m].I_real = N_aux[m].I_real;
                        N_int[m].I_imag = N_aux[m].I_imag;
                    }


               }

                //system("pause");
            }while(!similar_3);

            /*for(int m=0;m<ramos;m++){
                N_int[m].barra_1 = N_aux[m].barra_1;
                N_int[m].barra_2 = N_aux[m].barra_2;
                N_int[m].chave = N_aux[m].chave;
                N_int[m].id = N_aux[m].id;
                N_int[m].resist = N_aux[m].resist;
                N_int[m].reat = N_aux[m].reat;
                N_int[m].I_real = N_aux[m].I_real;
                N_int[m].I_imag = N_aux[m].I_imag;
            }*/

            printf("\n%d solucoes geradas na iteracao %d do PR!\n",sols,iter_pr);
            sols_total = sols_total + sols;
            //inicio++;
            //fim--;
            iter_pr++;
            //system("pause");
            //system("cls");
        }
    }
    printf("\nTotal de solucoes geradas no PR_1: %d\n",sols_total);
}

void Viz_1(struct V_Node *T_node, struct Node *N, struct Node *N_best, struct Node *N_worst,struct Tensao_Carga *TC, struct V_aux *vaux, int ramos, int barras, int subst, int **M_ciclos, int vet_ciclos[],
    float *p_at_2, float *p_reat_2, float *delta_P, float *P_per_1,float *P_per_2, float lim, int v_ini, float p_base, float *p_at_top, float *p_reat_top,
    struct Node M_pr[][16], int *n_pr, int *inicio_pr, int vet_block[], float FC_ini[], bool *achei_1, float *soma_fluxo_princ,
    FILE *arq3, FILE *arq5, FILE *arq6, int *bt_cont, float *perda_worst, float *p_at_aux, float *p_reat_aux, float *perdas_pr,
    int *p1, int *p2){

    int i, j, k1, k2, viz, w, m, it_pr;
    float soma_fluxo;
    bool infactivel;
    //viz = 0;

    for(i=0;i<(ramos-(barras-1));i++){
        //Salvando o índice do ramo inicialmente desligado do subciclo referente à ele
        for(j=0;j<ramos;j++){
            if(M_ciclos[i][j] == 3){
                k1 = j;
            }
        }
        //Ativando a chave do ramo desligado na linha i da Matriz de subciclos
        N[k1].chave = 1;
        for(j=0;j<ramos;j++){
            //achei = false;
            if((M_ciclos[i][j] == 1)||(M_ciclos[i][j] == 2)){
                k2 = j;
                //Desativando a chave do índice j da Matriz de subciclo
                N[k2].chave = 0;

                //if((vet_block[k1] == 0)||(vet_block[k2] == 0)){
                if((vet_block[k1] == 0)&&(vet_block[k2] == 0)){

                    soma_fluxo = 0.0;
                    for(int i_ciclo=0;i_ciclo<ramos;i_ciclo++){
                        if(N[i_ciclo].chave == 1){
                            soma_fluxo = (float) soma_fluxo + FC_ini[i_ciclo];
                        }
                    }

                    //if(soma_fluxo > (*soma_fluxo_princ)){
                    if((float) (soma_fluxo/(barras-1)) > (float) (0.99 * (((*soma_fluxo_princ))/(barras-1)))){
                    //if((float) (soma_fluxo/(barras-1)) > (float) (0.995 * ((soma_fluxo_princ)/(barras-1)))){

                        Apagar_Tree(T_node);
                        //viz++;
                        T_node->cont_p = 0;
                        VTree(N,vaux,ramos,subst,barras,T_node);
                        //printf("\nRamo desligado: %d - %d",N[k1].barra_1,N[k1].barra_2);
                        //printf("\nRamo ligado: %d - %d\n",N[k2].barra_1,N[k2].barra_2);
                        //printTree_pre(T_node);
                        *delta_P = 1.0;
                        //printTree_pont(T_node);
                        infactivel = false;
                        FCR(T_node,N,TC,barras,ramos,p_at_2,p_reat_2,delta_P,P_per_1,P_per_2,lim,v_ini);

                        for(int w=0;w<barras;w++){
                            if((TC[w].V_real < 0.93)||(TC[w].V_real > 1.05)){
                                infactivel = true;
                                break;
                            }
                        }

                        soma_fluxo = 0.0;
                        for(int i_ciclo=0;i_ciclo<ramos;i_ciclo++){
                            if(N[i_ciclo].chave == 1){
                                soma_fluxo = (float) soma_fluxo + FC_ini[i_ciclo];
                            }
                        }

                        //printf("aqui");
                        fprintf (arq3,"%f\n",(*p_at_2)*p_base);
                        fprintf (arq5,"%f\n",(float) soma_fluxo/(barras-1));
                        fprintf (arq6,"%d\n",(*bt_cont)+1);
                        //printf("\n\nVizinho %d --- Iteracao %d.",viz,bt_cont+1);
                        //printf("\nRamo desligado: %d - %d",N[k1].barra_1,N[k1].barra_2);
                        //printf("\nRamo da vez: %d - %d",N[k2].barra_1,N[k2].barra_2);
                        //printf("\nPerdas Ativas da Sol.: %f",p_at_2*p_base);
                        //printf("\nPerdas Reativas da Sol.: %f\n",p_reat_2*p_base);

                        //Guarda em N_worst a solução com a maior perda
                        if(!infactivel){
                            if((float) ((*p_at_2)*p_base) > (*perda_worst)){
                                (*perda_worst) = (float) ((*p_at_2)*p_base);
                                for(m=0;m<ramos;m++){
                                    N_worst[m].barra_1 = N[m].barra_1;
                                    N_worst[m].barra_2 = N[m].barra_2;
                                    N_worst[m].chave = N[m].chave;
                                    N_worst[m].id = N[m].id;
                                    N_worst[m].resist = N[m].resist;
                                    N_worst[m].reat = N[m].reat;
                                    N_worst[m].I_real = N[m].I_real;
                                    N_worst[m].I_imag = N[m].I_imag;
                                }
                            }

                            //system("pause");
                            if((float) ((*p_at_2)*p_base) < (float) ((*p_at_aux)*p_base)){

                                *p1 = k1;
                                *p2 = k2;

                                (*achei_1) = true;
                                *p_at_aux = (*p_at_2);
                                *p_reat_aux = (*p_reat_2);
                                for(m=0;m<ramos;m++){
                                    N_best[m].barra_1 = N[m].barra_1;
                                    N_best[m].barra_2 = N[m].barra_2;
                                    N_best[m].chave = N[m].chave;
                                    N_best[m].id = N[m].id;
                                    N_best[m].resist = N[m].resist;
                                    N_best[m].reat = N[m].reat;
                                    N_best[m].I_real = N[m].I_real;
                                    N_best[m].I_imag = N[m].I_imag;
                                }

                                //Inserindo a solução da iteração corrente na matriz/fila das melhores soluções
                                if((*n_pr) < 10){ //Insere elemento quando ainda há espaço para inserção
                                    for(it_pr=0;it_pr<ramos;it_pr++){
                                        M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].barra_1 = N[it_pr].barra_1;
                                        M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].barra_2 = N[it_pr].barra_2;
                                        M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].chave = N[it_pr].chave;
                                        M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].id = N[it_pr].id;
                                        M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].resist = N[it_pr].resist;
                                        M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].reat = N[it_pr].reat;
                                        M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].I_real = N[it_pr].I_real;
                                        M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].I_imag = N[it_pr].I_imag;
                                    }
                                    (*n_pr)++;
                                }else{ //Aqui a matriz se encontra cheia, então se remove um elemento para a entrada de um novo
                                    (*inicio_pr) = ((*inicio_pr) + 1) % 10;
                                    (*n_pr)--;
                                    for(it_pr=0;it_pr<ramos;it_pr++){
                                        M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].barra_1 = N[it_pr].barra_1;
                                        M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].barra_2 = N[it_pr].barra_2;
                                        M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].chave = N[it_pr].chave;
                                        M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].id = N[it_pr].id;
                                        M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].resist = N[it_pr].resist;
                                        M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].reat = N[it_pr].reat;
                                        M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].I_real = N[it_pr].I_real;
                                        M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].I_imag = N[it_pr].I_imag;
                                    }
                                    (*n_pr)++;
                                }
                                *perdas_pr = *p_at_aux;
                            }
                        }
                        //Retornando o valor da chave do índice j da Matriz de subciclos
                    }
                }
                N[k2].chave = 1;
            }else{
                continue;
            }
        }
        //Retornando o valor da chave do índice i da Matriz de subciclos
        N[k1].chave = 0;
    }
    //printf("\nQtde total de vizinhos em Viz_1: %d",viz);
}

void Viz_2(struct V_Node *T_node, struct Node *N, struct Node *N_best, struct Node *N_worst,struct Tensao_Carga *TC, struct V_aux *vaux, int ramos, int barras, int subst, int **M_ciclos, int vet_ciclos[],
    float *p_at_2, float *p_reat_2, float *delta_P, float *P_per_1,float *P_per_2, float lim, int v_ini, float p_base, float *p_at_top, float *p_reat_top,
    struct Node M_pr[][16], int *n_pr, int *inicio_pr, int vet_block[], float FC_ini[], bool *achei_2, float *soma_fluxo_princ,
    FILE *arq3, FILE *arq5, FILE *arq6, int *bt_cont, float *perda_worst, float *p_at_aux, float *p_reat_aux, float *perdas_pr,
    int *p1, int *p2){


    int i1, i2, j, j1, j2, k1, k2, k3, k4, viz, w, m, it_pr;
    float soma_fluxo;
    bool infactivel;
    viz = 0;

    for(i1=0;i1<(ramos-(barras-1))-1;i1++){
        for(j=0;j<ramos;j++){
            if(M_ciclos[i1][j] == 3){
                k1 = j;
                break;
            }
        }
        N[k1].chave = 1;
        for(i2=i1+1;i2<(ramos-(barras-1));i2++){
        //Salvando o índice do ramo inicialmente desligado do subciclo referente à ele
            for(j=0;j<ramos;j++){
                if(M_ciclos[i2][j] == 3){
                    k2 = j;
                    break;
                }
            }
            N[k2].chave = 1;
            for(j1=0;j1<ramos;j1++){
                //achei = false;
                if((M_ciclos[i1][j1] == 1)||(M_ciclos[i1][j1] == 2)){
                    k3 = j1;
                    //Desativando a chave do índice j da Matriz de subciclo
                    N[k3].chave = 0;

                    for(j2=0;j2<ramos;j2++){
                        if((M_ciclos[i2][j2] == 1)||(M_ciclos[i2][j2] == 2)){
                            k4 = j2;
                            N[k4].chave = 0;

                            //if((vet_block[k1] == 0)||(vet_block[k2] == 0)){
                            if((vet_block[k1] == 0)&&(vet_block[k2] == 0)&&(vet_block[k3] == 0)&&(vet_block[k4] == 0)){

                                soma_fluxo = 0.0;
                                for(int i_ciclo=0;i_ciclo<ramos;i_ciclo++){
                                    if(N[i_ciclo].chave == 1){
                                        soma_fluxo = (float) soma_fluxo + FC_ini[i_ciclo];
                                    }
                                }

                                //if(soma_fluxo > soma_fluxo_princ){
                                if((float) (soma_fluxo/(barras-1)) > (float) (0.995 * ((*soma_fluxo_princ)/(barras-1)))){
                                //if((float) (soma_fluxo/(barras-1)) > (float) (0.995 * ((soma_fluxo_princ)/(barras-1)))){

                                    Apagar_Tree(T_node);
                                    viz++;
                                    printf("\nVizinho %d de Viz_2!",viz);
                                    T_node->cont_p = 0;
                                    VTree(N,vaux,ramos,subst,barras,T_node);
                                    //printf("\nRamo desligado: %d - %d",N[k1].barra_1,N[k1].barra_2);
                                    //printf("\nRamo ligado: %d - %d\n",N[k2].barra_1,N[k2].barra_2);
                                    //printTree_pre(T_node);
                                    *delta_P = 1.0;
                                    //printTree_pont(T_node);
                                    infactivel = false;
                                    FCR(T_node,N,TC,barras,ramos,p_at_2,p_reat_2,delta_P,P_per_1,P_per_2,lim,v_ini);

                                    for(int w=0;w<barras;w++){
                                        if((TC[w].V_real < 0.93)||(TC[w].V_real > 1.05)){
                                            infactivel = true;
                                            break;
                                        }
                                    }

                                    soma_fluxo = 0.0;
                                    for(int i_ciclo=0;i_ciclo<ramos;i_ciclo++){
                                        if(N[i_ciclo].chave == 1){
                                            soma_fluxo = (float) soma_fluxo + FC_ini[i_ciclo];
                                        }
                                    }

                                    //printf("aqui");
                                    fprintf (arq3,"%f\n",(*p_at_2)*p_base);
                                    fprintf (arq5,"%f\n",(float) soma_fluxo/(barras-1));
                                    fprintf (arq6,"%d\n",(*bt_cont)+1);
                                    printf("\n\nVizinho %d --- Iteracao %d.",viz,(*bt_cont)+1);
                                    //printf("\nRamo desligado: %d - %d",N[k1].barra_1,N[k1].barra_2);
                                    //printf("\nRamo da vez: %d - %d",N[k2].barra_1,N[k2].barra_2);
                                    printf("\nPerdas Ativas da Sol.: %f",(*p_at_2)*p_base);
                                    //printf("\nPerdas Reativas da Sol.: %f\n",p_reat_2*p_base);

                                    //Guarda em N_worst a solução com a maior perda
                                    if(!infactivel){
                                        if((float) ((*p_at_2)*p_base) > (*perda_worst)){
                                            (*perda_worst) = (float) ((*p_at_2)*p_base);
                                            for(m=0;m<ramos;m++){
                                                N_worst[m].barra_1 = N[m].barra_1;
                                                N_worst[m].barra_2 = N[m].barra_2;
                                                N_worst[m].chave = N[m].chave;
                                                N_worst[m].id = N[m].id;
                                                N_worst[m].resist = N[m].resist;
                                                N_worst[m].reat = N[m].reat;
                                                N_worst[m].I_real = N[m].I_real;
                                                N_worst[m].I_imag = N[m].I_imag;
                                            }
                                        }

                                        //system("pause");
                                        if((float) ((*p_at_2)*p_base) < (float) ((*p_at_aux)*p_base)){

                                            *p1 = k3;
                                            *p2 = k4;

                                            (*achei_2) = true;
                                            *p_at_aux = (*p_at_2);
                                            *p_reat_aux = (*p_reat_2);
                                            for(m=0;m<ramos;m++){
                                                N_best[m].barra_1 = N[m].barra_1;
                                                N_best[m].barra_2 = N[m].barra_2;
                                                N_best[m].chave = N[m].chave;
                                                N_best[m].id = N[m].id;
                                                N_best[m].resist = N[m].resist;
                                                N_best[m].reat = N[m].reat;
                                                N_best[m].I_real = N[m].I_real;
                                                N_best[m].I_imag = N[m].I_imag;
                                            }

                                            //Inserindo a solução da iteração corrente na matriz/fila das melhores soluções
                                            if((*n_pr) < 10){ //Insere elemento quando ainda há espaço para inserção
                                                for(it_pr=0;it_pr<ramos;it_pr++){
                                                    M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].barra_1 = N[it_pr].barra_1;
                                                    M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].barra_2 = N[it_pr].barra_2;
                                                    M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].chave = N[it_pr].chave;
                                                    M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].id = N[it_pr].id;
                                                    M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].resist = N[it_pr].resist;
                                                    M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].reat = N[it_pr].reat;
                                                    M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].I_real = N[it_pr].I_real;
                                                    M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].I_imag = N[it_pr].I_imag;
                                                }
                                                (*n_pr)++;
                                            }else{ //Aqui a matriz se encontra cheia, então se remove um elemento para a entrada de um novo
                                                (*inicio_pr) = ((*inicio_pr) + 1) % 10;
                                                (*n_pr)--;
                                                for(it_pr=0;it_pr<ramos;it_pr++){
                                                    M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].barra_1 = N[it_pr].barra_1;
                                                    M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].barra_2 = N[it_pr].barra_2;
                                                    M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].chave = N[it_pr].chave;
                                                    M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].id = N[it_pr].id;
                                                    M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].resist = N[it_pr].resist;
                                                    M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].reat = N[it_pr].reat;
                                                    M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].I_real = N[it_pr].I_real;
                                                    M_pr[((*inicio_pr) + (*n_pr))%10][it_pr].I_imag = N[it_pr].I_imag;
                                                }
                                                (*n_pr)++;
                                            }
                                            *perdas_pr = *p_at_aux;
                                        }
                                    }
                                    //Retornando o valor da chave do índice j da Matriz de subciclos
                                }
                            }
                            if(k3 != k4){
                                N[k4].chave = 1;
                            }
                        }
                    }
                    N[k3].chave = 1;
                }
            }
            //Retornando o valor da chave do índice i da Matriz de subciclos
            N[k2].chave = 0;
        }
        N[k1].chave = 0;
    }
    printf("\nQtde total de vizinhos em Viz_2: %d",viz);
}
