#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <gmp.h>
#include <time.h>
#include <string.h>
// ---------
#include "mpi.h"

int aPath[100][2];

int LookHorizontaly(int u, int v, int u1, int v1, int n, int m, int aRoute[n][m]){
	int i =0;
	int j =0;
    for (i=0; i < m; i++){
        if(i!=v && aRoute[u][i] != 0){
            if(i == v1){
				for(j=0; j<100; j++){
					if(aPath[j][0] == -12345 && aPath[j][1] == -12345){
						aPath[j][0]=u;
						aPath[j][1]=i;
						return 1;
					}
				}
                return 1;
			}
            if(LookVerticaly(u, i, u1, v1, n, m, aRoute) == 1){
				for(j=0; j<100; j++){
					if(aPath[j][0] == -12345 && aPath[j][1] == -12345){
						aPath[j][0]=u;
						aPath[j][1]=i;
						return 1;
					}
				}
                return 1;
			}
		}
	}
	return 0;
}
int LookVerticaly(int u, int v, int u1, int v1, int n, int m, int aRoute[n][m]){
	int i=0;
	int j =0;
    for (i=0; i < n; i++){
        if(i != u && aRoute[i][v] != 0){
            if(LookHorizontaly(i, v, u1, v1, n, m, aRoute) == 1){
				for(j=0; j < 100; j++){
					if(	aPath[j][0] == -12345 && aPath[j][1] == -12345){
						aPath[j][0]=i;
						aPath[j][1]=v;
						return 1;
					}
				}
                return 1;
			}
		}
	}
    return 0;
}
void PrintOut(int m, int n, int aSupply[m], int aDemand[n], int aCost[n][m], int aDual[n][m], int aRoute[n][m]){
	int nCost = 0;
    int x,y,i;
    for (i=0; i < m; i++){
        printf("%i		", aDemand[i]);
    }
    printf("\n");

    for (x = 0; x < n; x++){
        for (y = 0; y < m; y++){
            nCost += aCost[x][y] * aRoute[x][y];
            if(aRoute[x][y] == 0){
                printf("[<%i>  %i]  	", aCost[x][y], aDual[x][y]);
            }else{
                printf("[<%i>(%i)]  	", aCost[x][y], aRoute[x][y]);
            }
        }
       printf(": %i\n", aSupply[x]);
    }
	printf("Cost: %i\n\n", nCost);
}

int main(int argc, char *argv[]){

	//Reading inputs.txt to get a list of input files
	int NumInputs=0;
	FILE * fi;
	char * line = NULL;
	size_t flen = 0;
	ssize_t read;
	fi = fopen("inputs.txt", "r");
	if (fi == NULL){
		printf("[ERROR] inputs.txt Not Found\n");
		return 0;
	}
	while ((read = getline(&line, &flen, fi)) != -1) {
		NumInputs++;
	}
	rewind(fi);
	char inputFiles[NumInputs][20];
	int i=0;
	int j=0;

	while ((read = getline(&line, &flen, fi)) != -1) {
		strtok(line, "\r\n");
		strcpy(inputFiles[i], line);
		i++;
	}
	fclose(fi);
	if(NumInputs==0){
		printf("[Error] inputs.txt is empty\n");
		return 0;
	}
	// -- finished reading inputs.txt

	int rank,size,len;
	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;
    MPI_Request request;

	int remainder;
	int rank_start[size];
	int rank_end[size];

	// Divide Chunks
	if(NumInputs>=size){
		int remainder = NumInputs % size;
		if (rank==0){
			if (remainder==0){
				for(i=0; i<size; i++){
					rank_start[i]=i*(NumInputs/size);
					rank_end[i]=(i+1)*(NumInputs/size);
				}
			}else{
				for(i=0; i<size; i++){
					rank_start[i]=i*(NumInputs/size);
					rank_end[i]=(i+1)*(NumInputs/size);
				}
				rank_start[size-1]=(size-1)*(NumInputs/size);
				rank_end[size-1]=NumInputs;
			}
			for(i=1; i<size; i++){
				MPI_Send(&rank_start[i],1,MPI_INT,i,1,MPI_COMM_WORLD);
				MPI_Send(&rank_end[i],1,MPI_INT,i,1,MPI_COMM_WORLD);
			}
		}else{
			MPI_Recv(&rank_start[rank],1,MPI_INT,0,1,MPI_COMM_WORLD,&status);
			MPI_Recv(&rank_end[rank],1,MPI_INT,0,1,MPI_COMM_WORLD,&status);
		}
	}else{
		for(i=0; i<size; i++){
			rank_start[i]=i;
			rank_end[i]=i+1;
		}
	}
	
	int	irank=0;

	//printf("Rank %i: Start %i, End %i\n", rank, rank_start[rank], rank_end[rank]);
	// Read file(s) from specific rank then calculate result
	for (irank=rank_start[rank]; irank < rank_end[rank] && irank < NumInputs; irank++){

		FILE * fr;
		char * line = NULL;
		char * pch;
		size_t flen = 0;
		ssize_t read;

		fr = fopen(inputFiles[irank], "r"); // Read Input File
		if (fr == NULL){
			printf("[ERROR] Rank %i %s Not Found\n", rank, inputFiles[irank]);
			MPI_Finalize();
			return 0;
		}else{
			printf("-- Rank %i Reading inputs from %s\n", rank, inputFiles[irank]);
		}

		int lineCount, colCount, rowCount, NS, ND;
		lineCount=0;
		// Get Supply Count and Demand Count
		while ((read = getline(&line, &flen, fr)) != -1) {
			colCount=0;
			pch = strtok (line,", ");
			while (pch != NULL){
				colCount++;
				pch = strtok(NULL, ", ");
			}
			if(lineCount==0){
				NS = colCount;
			}else if(lineCount==1){
				ND = colCount;
			}else{
				break;
			}
			lineCount++;
		}
		fclose(fr);

		int aCost[NS][ND];
		int aSupply[NS];
		int aDemand[ND];

		
		fr = fopen(inputFiles[irank], "r");
		lineCount=0;
		rowCount=0;

		// Get Supply, Demand, and Cost
		while ((read = getline(&line, &flen, fr)) != -1) {
			pch = strtok (line,", ");
			colCount=0;
			while (pch != NULL){
				if(lineCount==0){
					aSupply[colCount]=atoi(pch);
				}else if(lineCount==1){
					aDemand[colCount]=atoi(pch);
				}else{
					aCost[rowCount][colCount]=atoi(pch);
					//printf("aCost[%i][%i]=%i\n",rowCount,colCount,aCost[rowCount][colCount]);
				}
				colCount++;
				pch = strtok(NULL, ", ");
			}
			if(lineCount>1){
				rowCount++;
			}
			lineCount++;
		}

		fclose(fr);
		if (line){
			free(line);
		}

		// Create Output file
		char fileWrite[20];
		sprintf(fileWrite, "output%i.txt", irank+1);
		FILE *fw = fopen(fileWrite, "w");
		if (fw == NULL){
			printf("Error opening file!\n");
			exit(1);
		}

		fprintf(fw, "#Input:\n\n");

		for(i=0;i<ND;i++){
			fprintf(fw, "%4d", aDemand[i]);
		}
		fprintf(fw, "\n----------------------------------\n");
		for(i=0;i<NS;i++){
			for(j=0;j<ND;j++){
				fprintf(fw, "%4d", aCost[i][j]);
			}
			fprintf(fw, "	:	%i", aSupply[i]);
			fprintf(fw, "\n");
		}
		fprintf(fw, "\n");

		fprintf(fw, "#Output:\n\n", fileWrite);

		int n = sizeof(aSupply)/sizeof(aSupply[0]);
		int m = sizeof(aDemand)/sizeof(aDemand[0]);
		long nVeryLargeNumber = 99999999999;
		double elipsis = 0.001;
		int i,j;
		int aRoute[n][m];
		int aDual[n][m];
		int u=0;
		int v=0;
		int z=0;
		int aS[m];
		int aD[n];
		int PivotN = -1;
		int PivotM = -1;
		int aPathMax = 1;
		int x;
		int chunk,tid;
		chunk=4;
		omp_set_num_threads(4);
		#pragma omp parallel
    	{
	        tid = omp_get_thread_num();
	        #pragma omp for schedule (static, chunk)
			for (i=0; i < n; i++){
				for (j=0; j< m; j++){
					aRoute[i][j]=0;
					
				}
			}
			#pragma omp for schedule (static, chunk)
			for (i=0; i<n; i++){
				for (j=0; j<m; j++){
					aDual[i][j]=-1;

				}
			}
			
			//-- NorthWest 
			#pragma omp for schedule (static, chunk)
			for (i=0; i<m; i++){
				aS[i]=0;
			}
			#pragma omp for schedule (static, chunk)
			for (i=0; i<n; i++){
				aD[i]=0;
			}

			#pragma omp single
			{
				while (u < n && v < m){
					if( (aDemand[v] - aS[v]) < (aSupply[u] - aD[u])){
						z = aDemand[v] - aS[v];
						aRoute[u][v] = z;
						aS[v] +=z;
						aD[u] +=z;
						v++;
					}else{
						z = aSupply[u] - aD[u];
						aRoute[u][v] = (double)z;
						aS[v] +=z;
						aD[u] +=z;
						u++;
					}
				}
			}//-- End: NorthWest


    		tid = omp_get_thread_num();
	        #pragma omp for schedule (static, chunk)
			for(i=0; i < 100; i++){
				aPath[i][0] = -12345;
				aPath[i][1] = -12345;
				
			}
			
			//-- GetDual
			#pragma omp for schedule (static, chunk)
			for (u=0; u < n; u++){
				for (v = 0; v < m; v++){
					aDual[u][v] = -1;//Null Value
					if(aRoute[u][v] == 0){
						for(i=0; i < 100; i++){//Reset aPath
							aPath[i][0] = -12345;
							aPath[i][1] = -12345;
						}
						aPath[0][0] = u;
						aPath[0][1] = v;
						if(LookHorizontaly(u, v, u, v, n, m, aRoute) == 0){
							//printf("rank %d Path Error\n",rank);
						}
						z = -1;
						x = 0;
						aPathMax=1;
						for(i=0; i < 100 && aPathMax; i++){
							if(aPath[i][0] != -12345 && aPath[i][1] != -12345){
								x += z * aCost[aPath[i][0]][aPath[i][1]];
								z *= -1;
							}else{
								aPathMax=0;
							}
						}
						aDual[u][v] = x;
					}
				}
			}
		}
		//-- End: GetDual

		// PrintOut NorthWest to Output file
		// PrintOut(m,n, aSupply, aDemand, aCost, aDual, aRoute);
		fprintf(fw, "North West:\n", fileWrite);
		int nCost = 0;
		int y;
		for (i=0; i < m; i++){
			fprintf(fw, "%11d   ", aDemand[i]);
		}
		fprintf(fw, "\n");
		for (x = 0; x < n; x++){
			for (y = 0; y < m; y++){
				nCost += aCost[x][y] * aRoute[x][y];
				if(aRoute[x][y] == 0){
					fprintf(fw, "[<%2d>  %4d ] ", aCost[x][y], aDual[x][y]);
				}else{
					fprintf(fw, "[<%2d> (%4d)] ", aCost[x][y], aRoute[x][y]);
				}
			}
			fprintf(fw, ": %i\n", aSupply[x]);
		}
		fprintf(fw, "Cost: %i\n\n", nCost);
		//End: PrintOut

		int IsOptimal=0;
		int w, t, nMin, nMax, aPathLen;
		
		nCost = 0;
		while(IsOptimal==0){
			nMax = -1;

			for (u=0; u < n; u++){
				for (v = 0; v < m; v++){
					aDual[u][v] = -1;//Null Value
					if(aRoute[u][v] == 0){
						for(i=0; i < 100; i++){//Reset aPath
							aPath[i][0] = -12345;
							aPath[i][1] = -12345;
						}
						aPath[0][0] = u;
						aPath[0][1] = v;
						if(LookHorizontaly(u, v, u, v, n, m, aRoute) == 0){
							//printf("Rank %d Path Error\n", rank);
						}
						z = -1;
						x = 0;
						aPathMax=1;
						for(i=0; i < 100 && aPathMax; i++){
							if(aPath[i][0] != -12345 && aPath[i][1] != -12345){
								x += z * aCost[aPath[i][0]][aPath[i][1]];
								z *= -1;
							}else{
								aPathMax=0;
							}
						}
						aDual[u][v] = x;
					}
				}
			}
			//End: GetDual

			for (u=0; u < n; u++){
				for (v=0; v < m; v++){
					x = aDual[u][v];
					if (x > nMax){
						nMax = x;
						PivotN = u;
						PivotM = v;
					}
				}
			}

			if(nMax <= 0){
				IsOptimal=1;
				fprintf(fw, "Finished! Final Cost: %i\n", nCost);
			}else{
				fprintf(fw, "Better Optimized:\n");
				//BetterOptimal
				//-- FindPath
				for(i=0; i < 100; i++){//Reset aPath
					aPath[i][0] = -12345;
					aPath[i][1] = -12345;
				}
				aPath[0][0] = PivotN;
				aPath[0][1] = PivotM;
				aPathLen=1;
				if(LookHorizontaly(PivotN, PivotM, PivotN, PivotM, n, m, aRoute) == 0){
					printf("Path Error\n");
				}
				for(i=0; i < 100; i++){
					if(aPath[i][0]!= -12345 && aPath[i][1]!= -12345){
						aPathLen++;
					}
				}
				//-- End: FindPath

				nMin  = 999999;
				for (w=1; w < aPathLen; w=w+2){
					t = aRoute[aPath[w][0]][aPath[w][1]];
					if(t < nMin){
						nMin = t;
					}
				}
				for (w=1; w < aPathLen; w=w+2){
					aRoute[aPath[w][0]][aPath[w][1]]         -= nMin;
					aRoute[aPath[w - 1][0]][aPath[w - 1][1]] += nMin;
				}				
				//End: BetterOptimal	

				//-- PrintOut BetterOptimal to Output File
				//PrintOut(m,n, aSupply, aDemand, aCost, aDual, aRoute);
				nCost = 0;
				for (i=0; i < m; i++){
					fprintf(fw, "%11d   ", aDemand[i]);
				}
				fprintf(fw, "\n");
				for (x = 0; x < n; x++){
					for (y = 0; y < m; y++){
						nCost += aCost[x][y] * aRoute[x][y];
						if(aRoute[x][y] == 0){
							fprintf(fw, "[<%2d>  %4d ] ", aCost[x][y], aDual[x][y]);
						}else{
							fprintf(fw, "[<%2d> (%4d)] ", aCost[x][y], aRoute[x][y]);
						}
					}
					fprintf(fw, ": %i\n", aSupply[x]);
				}
				fprintf(fw, "Cost: %i\n\n", nCost);
				//End: PrintOut
			}
		}
		fclose(fw);
		//PrintOut(m,n, aSupply, aDemand, aCost, aDual, aRoute);
		printf("\nProduct %s\n", inputFiles[irank]);
		for (x = 0; x < n; x++){
			for (y = 0; y < m; y++){
				if(aRoute[x][y] > 0){
					printf("src:%i --> dest:%i  quantity:%i\n", x, y, aRoute[x][y]);
				}
			}
		}
		printf("Printed result in [output%i.txt]\n", irank+1);
	}
	MPI_Finalize();
	return 0;
}