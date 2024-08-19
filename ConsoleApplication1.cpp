//#include "stdafx.h"
#include <ilcplex/ilocplex.h>
#include <time.h>
#include <fstream>
#include <math.h>
#include <vector>


ILOSTLBEGIN

typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumVarArray2> IloNumVarArray3;

using namespace std;

time_t t_ini, t_fim; //referente a time.h
double tempo;
double tempototal;
const float lowvalue = 0.001;

int main() {
	IloEnv env;
	try
	{
		//Numero de periodos t (T), numero de biodigestores i (N), periodo j do biodigestor J

		IloInt indiceT, indiceN, indiceJ, indiceT_aux;
		ifstream Arq2;
		Arq2.open("revistaB_sazonal_5_8_11.txt", ios::app);
		Arq2 >> indiceN; //numero de biodigestores
		Arq2 >> indiceT; //numero de periodos
		Arq2 >> indiceT_aux; //numero de periodos auxiliares (pré e pós demanda)



		IloNumArray DemT(env, indiceT); //tamanho T períodos
		IloNumArray TiN(env, indiceN); //Tamanho N biodigestores

		//Ler a demanda

		for (int t = 0; t < indiceT; t++)
		{
			Arq2 >> DemT[t];
		}




		//Ler a ti (vetor com tamanho N periodos)
		cout << "tamanho de cada" << endl;
		for (int i = 0; i < indiceN; i++)
		{

			Arq2 >> TiN[i];

		}

		//Producao do biodigestor i no periodo j
		IloArray<IloNumArray> Pij(env, indiceN);
		for (int i = 0; i < indiceN; i++)
			Pij[i] = IloNumArray(env, TiN[i]);
		for (int i = 0; i < indiceN; i++)
		{
			for (int j = 0; j < TiN[i]; j++)
			{//preencher o mes
				Arq2 >> Pij[i][j];
				cout << "prod vale" << Pij[i][j] << endl;
			}
			//			cout << "tamo junto" << endl;

		}


		time(&t_ini);
		IloModel mod(env);

		//Variaveis de decisao (Kijt)

		IloNumVarArray3 Kijt(env, indiceN);

		IloNumArray3 resp_Kijt(env, indiceN);

		for (int i = 0; i < indiceN; i++) {
			Kijt[i] = IloNumVarArray2(env, TiN[i]);
			resp_Kijt[i] = IloNumArray2(env, TiN[i]);

			for (int j = 0; j < TiN[i]; j++) {
				Kijt[i][j] = IloNumVarArray(env, indiceT);
				resp_Kijt[i][j] = IloNumArray(env, indiceT);

				for (int t = 0; t < indiceT; t++) {
					Kijt[i][j][t] = IloNumVar(env, 0.0, 1.0, ILOINT);
				}
			}
		}

	
		//Variavel Prodt
		IloNumVarArray ProdT(env);
		IloNumArray resp_prod(env);

		for (int t = 0; t < indiceT; t++) {
			ProdT.add(IloNumVar(env, 0, IloInfinity, ILOINT));

		}


		//Variavel de decisao binária

		IloNumVarArray2 Yit(env, indiceN);
		IloNumArray2 resp_Yi(env, indiceN);


		for (int i = 0; i < indiceN; i++) {
			Yit[i] = IloNumVarArray(env, indiceT);
			resp_Yi[i] = IloNumArray(env, indiceT);
			for (int t = 0; t < indiceT; t++) {
				Yit[i][t] = IloNumVar(env, 0.0, 1, ILOINT);
			}
		}





		//FO
		//preciso excluir períodos pré e pós demanda (são auxiliares)
		//NÃO posso incluir desvio de períodos auxiliar na FO
		

		IloExpr fo(env);
		
		for (int t = 0; t < indiceT; t++) {
			fo += ProdT[t]; //somando todos os desvios
		}


		mod.add(IloMinimize(env, fo));
		fo.end();

		//Rest2

		//IloExpr expr2(env);

		IloInt inicioProducao = 0; //esse é o Prod[0]
		IloInt inicioDemanda = 0; //esse é o DemT[0]

		//Restricao 2: sem estoque anterior
		IloConstraint c3(inicioProducao == ProdT[0]);
		c3.setName("periodo zero");
		mod.add(c3);
		
	

		for (int t = 1; t < indiceT; t++)
		{
			IloExpr expr3(env);
			for (int i = 0; i < indiceN; i++) {
				for (int j = 0; j < TiN[i]; j++) {
					expr3 += Kijt[i][j][t] * Pij[i][j]; //quantidade produzida
				}
			}
			
			expr3 = expr3 + ProdT[t - 1] - DemT[t];
			
			IloConstraint c(expr3 == ProdT[t]);
			char O[50];
			int nome;
			nome = sprintf_s(O, "IT%d", t);
			c.setName(O);
			mod.add(c);
			expr3.end();
		}


		
		
		//Rest 6: biodigestor só pode ligar em um dia do periodo
		for (int i = 0; i < indiceN; i++)  //para cada biodigestor
		{
			for (int t = 0; t < indiceT; t++)
			{
				IloExpr expr6(env);
				for (int j = 0; j < TiN[i]; j++)
				{
					expr6 = expr6 + Kijt[i][j][t];
				}
				char O[50];
				int nome;
				nome = sprintf_s(O, "O%d%d", i, t);
				IloConstraint c6(expr6 <= 1);
				c6.setName(O);
				mod.add(c6);
				expr6.end();
			}
		}

		//	cout << "agora restricao de ligar" << endl;

			//cada bio tem seu proprio ti. Entao preciso limitar de acordo com o bio
		cout << "\nKKK" << endl;
		
		//Limitando os valores de z (iniciais)
		//cout << "Limitando os valores de K predemanda)" << endl;
		for (int i = 0; i < indiceN; i++)
		{
			for (int t = 0; t < TiN[i]; t++)
			{
				
				for (int j = t + 1; j < TiN[i]; j++)
				{
					IloExpr expr5(env);
					expr5 = Kijt[i][j][t];
					IloConstraint c7(expr5 == 0);
					mod.add(c7);
					expr5.end();
					
				}
			}
		}
	
		int aux = 0;
		//Limitando os valores de z (finais)
		for (int i = 0; i < indiceN; i++) { //todos os biodigestores	
			aux = 0;
			for (int t = (indiceT - TiN[i]); t < indiceT; t++) { //apenas os ultimos periodos
				aux++;

				for (int j = 0; j < aux; j++) //nao precisa ir ate Tin[i] pois os periodos adicionais sao criados baseados no tin[i]
				{
					Kijt[i][j][t].setUB(0);
				}
			}
		}
		

		//Restricao para que uma sequencia inteira de dias seja ativada junta, ou não seja nenhum dia

		IloInt TTT = 0;
		for (int i = 0; i < indiceN; i++) //para cada biodigestor
		{
			for (int t = 0; t < indiceT - indiceT_aux; t++) //0 até 5-3+1 ou 20-4 
			{

				
				//getchar();
				//getchar();
				IloExpr expr6(env);
				TTT = t;
				for (int j = 0; j < TiN[i]; j++)
				{
					expr6 += Kijt[i][j][TTT];
					TTT = TTT + 1;
				}

				IloConstraint c8(expr6 == Yit[i][t] * TiN[i]);
				char O[50];
				int nome;
				nome = sprintf_s(O, "IT%d%d", i, t);
				c8.setName(O);
				mod.add(c8);
				expr6.end();
			}
		}


		//n pode iniciar producao no pos demanda
		//nenhum y pode ser positivo no pos demanda
		for (int i = 0; i < indiceN; i++) //para cada biodigestor
		{
			
			IloExpr expr10(env);
			for (int t = indiceT - indiceT_aux; t < indiceT; t++) //
			{
				cout << "rest: " << t << endl;
				expr10 = Yit[i][t];
				IloConstraint c10(expr10 == 0);
				c10.setName("Restricao de forçar");
				mod.add(c10);
				expr10.end();
			}
		}
		
		IloCplex cplex(mod);
		cplex.setParam(IloCplex::Param::TimeLimit, 7200); //limite de 2 horas

		/*
		//0 = meio termo, 1 = feasibility, 2 = optimalitty, 3 mais optimality, 4 mais feasibility
		cplex.setParam(IloCplex::Param::Emphasis::MIP, 3);

		//-1 sem probing, 0 automatico, 1,2,3 moderado ate mt agressivo
		cplex.setParam(IloCplex::Param::MIP::Strategy::Probe, 3);

		//cortes agressivos para deixar o modelo mais restrito
		cplex.setParam(IloCplex::Param::MIP::Cuts::Covers,3);
		cplex.setParam(IloCplex::Param::MIP::Cuts::Cliques, 3);
		//cplex.setParam(IloCplex::Param::MIP::Cuts::Disjunctive, 3);
		cplex.setParam(IloCplex::Param::MIP::Cuts::LiftProj, 3);
		cplex.setParam(IloCplex::Param::MIP::Cuts::LocalImplied, 3);
	*/

		if (!cplex.solve()) {
			cplex.exportModel("artigo.lp");
			env.error() << "Erro ao obter resposta." << endl;
			env.out() << "Status = " << cplex.getStatus() << endl;
			env.out() << "FO = " << cplex.getObjValue() << endl;
			throw(-1);
		}


		time(&t_fim);
		double tempototal = 0;
		tempototal = (t_fim - t_ini) / double(CLOCKS_PER_SEC) * 1000;

		env.out() << "Status = " << cplex.getStatus() << endl;
		env.out() << "FO = " << cplex.getObjValue() << endl;
		cplex.exportModel("artigo.lp");
		
		for (int i = 0; i < indiceN; i++)
		{
			for (int j = 0; j < TiN[i]; j++)
				cplex.getValues(resp_Kijt[i][j], Kijt[i][j]);
		}
		
		for (int i = 0; i < indiceN; i++)
		{
			cplex.getValues(resp_Yi[i], Yit[i]);
		}
		


		//TRATAMENTO VARIAVEIS DE DECISAO
		for (int i = 0; i < indiceN; i++)
		{
			for (int t = 0; t < indiceT; t++)
			{
				if (resp_Yi[i][t] < 1)
				{
					resp_Yi[i][t] = 0;
				}
				else
					resp_Yi[i][t] = 1;
			}
		}

		for (int i = 0; i < indiceN; i++)
		{
			for (int t = 0; t < indiceT; t++)
			{

				for (int j = 0; j < TiN[i]; j++)
				{

					if (resp_Kijt[i][j][t] < 1)
					{
						resp_Kijt[i][j][t] = 0;
					}
					else
						resp_Kijt[i][j][t] = 1;
				}
			}
		}



		cplex.exportModel("artigo.lp");
		IloNum media;
		media = cplex.getObjValue() / (indiceT);

		cplex.getValues(resp_prod, ProdT);
		
		ofstream myfile;
		myfile.open("saida_revistaB_sazonal_5_8_11_TESTE4.txt");
		
		

		for (int i = 0; i < indiceN; i++)
		{
			int primeiravez = 1;
			int primeirodia = 1000;
			int numerodevezes = 0;
			myfile << "\nBIO: " << i << ": " << endl;
			for (int t = 0; t <indiceT; t++)
			{
				myfile << resp_Yi[i][t] << " ";
				if (resp_Yi[i][t] == 1) 
				{
					numerodevezes++;
					if (primeiravez == 1)
					{
						primeirodia = t + 1;
						primeiravez = 2;
					}

				}
				
			}
			myfile << "\nNumber of used times: " << numerodevezes << endl;

			if(primeirodia == 1000)
				myfile << "BIO NOT USED: " << endl;
			else			
				myfile << "First Usage day: " << primeirodia << endl;
		}
		
		
		vector <int> resp_dtotal;
		IloNum varsum = 0;
		IloNum soma2 = 0;
		

		for (int t = 0; t <indiceT; t++)
			soma2 = soma2 + resp_prod[t];

		IloNum media2 = soma2 / indiceT;

		for (int t = 0; t < indiceT; t++)
			varsum = varsum + pow((resp_prod[t] - media2), 2);

		IloNum Variance = varsum / indiceT;
		IloNum desviopadrao = sqrt(Variance);

		
		myfile << "\n\nMin =" << cplex.getObjValue() << endl;
		myfile << "Time =" << tempototal << endl;
		myfile << "GAP =" << cplex.getMIPRelativeGap() << endl;
		myfile << "AABS = " << media << endl;
		myfile << "AABS (manually obtained) = " << media2 << endl;
		myfile << "Standard deviation = " << desviopadrao << endl;
		
		for (int i = 0; i < indiceN; i++)
		{
			myfile << "\n\nBiodigester " << i << ": ";
			for (int t = 0; t < indiceT; t++)
			{
				myfile << "\nDay (horizon)" << t << ":";
				for (int j = 0; j < TiN[i]; j++)
				{
					myfile << "\nday (cycle) J" << j << ": ";

					myfile << resp_Kijt[i][j][t] << " ";
				}
			}
		}
		myfile << "\n\nProduction" << endl;
		for (int t = 0; t < indiceT; t++)
		{
			myfile << resp_prod[t] << " ";
		}

		myfile << "\n\nDemand \n\n" << endl;
		for (int t = 0; t < indiceT; t++)
		{
			myfile << DemT[t] << " ";
		}

	
		

		//int total[3][423]; //total[i][t]
		int soma = 0;
		for (int i = 0; i < indiceN; i++)
		{
			myfile << "\nBiodigester: " << i << ": " << endl;
			for (int t = 0; t < indiceT; t++)
			{
				soma = 0;
				for (int j = 0; j < TiN[i]; j++)
				{
					soma = soma + Pij[i][j] * resp_Kijt[i][j][t];
				}
				myfile << soma << " ";
			}
		}
		

	}

	catch (IloException& e) {
		cerr << "Concert exception caught: " << e << endl;
		ofstream myfile;
		myfile.open("saida_revista_constante_5.txt");
		myfile << "Concert exception caught: " << e << endl;
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
		ofstream myfile;
		myfile.open("saida_revista_constante_5.txt");
		myfile << "Unknown exception caught" << endl; 
	}

	env.end();
	getchar();

}