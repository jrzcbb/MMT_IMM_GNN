
/*
function [assignment, cost] = assignmentoptimal(distMatrix)
*/

//#include <mex.h>
//#include <matrix.h>

//���������ȫ��Ϊ����Ӧmatlabд��
//#include "StdAfx.h"
#include <stdlib.h>
#include <string.h>
#include "GNNTrakerCom.h"
#include <stdio.h>

//#define CHECK_FOR_INF
#define ONE_INDEXING

//#define mxMalloc malloc

void assignmentoptimal(double *assignment, double *cost, double *distMatrix, int nOfRows, int nOfColumns);
void buildassignmentvector(double *assignment, bool *starMatrix, int nOfRows, int nOfColumns);
void computeassignmentcost(double *assignment, double *cost, double *distMatrix, int nOfRows);
void step2a(double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
void step2b(double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
void step3 (double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
void step4 (double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim, int row, int col);
void step5 (double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);

//void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
//{
//	double *assignment, *cost, *distMatrix;
//	int nOfRows, nOfColumns;
//	
//	/* Input arguments */
//	nOfRows    = mxGetM(prhs[0]);
//	nOfColumns = mxGetN(prhs[0]);
//	distMatrix = mxGetPr(prhs[0]);
//	
//	/* Output arguments */
//	plhs[0]    = mxCreateDoubleMatrix(nOfRows, 1, mxREAL);
//	plhs[1]    = mxCreateDoubleScalar(0);
//	assignment = mxGetPr(plhs[0]);
//	cost       = mxGetPr(plhs[1]);
//	
//	/* Call C-function */
//	assignmentoptimal(assignment, cost, distMatrix, nOfRows, nOfColumns);	
//}

//nOfRows �� columnsLen , nOfColumns �� RowsLen
void assignmentoptimal(double *assignment, double *cost, double *distMatrixIn, int nOfRows, int nOfColumns)
{
	double *distMatrix, *distMatrixTemp, *distMatrixEnd, *columnEnd, value, minValue;
	bool *coveredColumns, *coveredRows, *starMatrix, *newStarMatrix, *primeMatrix;
	int nOfElements, minDim, row, col;
#ifdef CHECK_FOR_INF
	bool infiniteValueFound;
	double maxFiniteValue, infValue;
#endif
	
	/* initialization */
	*cost = 0;
	for(row=0; row<nOfRows; row++)
#ifdef ONE_INDEXING
		assignment[row] =  0.0;
#else
		assignment[row] = -1.0;
#endif
	
	/* generate working copy of distance Matrix */
	/* check if all matrix elements are positive */
	nOfElements   = nOfRows * nOfColumns;

	//distMatrix    = (double *)mxMalloc(nOfElements * sizeof(double));
	distMatrix    = (double *)malloc(nOfElements * sizeof(double));
	
	distMatrixEnd = distMatrix + nOfElements;
	for(row=0; row<nOfElements; row++)
	{
		value = distMatrixIn[row];
		//if(mxIsFinite(value) && (value < 0))
		//	mexErrMsgTxt("All matrix elements have to be non-negative.");
		if (value < 0)
		{
			printf("All matrix elements have to be non-negative.");
		}
		distMatrix[row] = value;
	}

#ifdef CHECK_FOR_INF
	/* check for infinite values */
	maxFiniteValue     = -1;
	infiniteValueFound = false;
	
	distMatrixTemp = distMatrix;
	while(distMatrixTemp < distMatrixEnd)
	{
		value = *distMatrixTemp++;

		//if(mxIsFinite(value))
		if (MAX_INF == value)
		{
			if(value > maxFiniteValue)
				maxFiniteValue = value;
		}
		else
			infiniteValueFound = true;
	}
	if(infiniteValueFound)
	{
		if(maxFiniteValue == -1) /* all elements are infinite */
			return;
		
		/* set all infinite elements to big finite value */
		if(maxFiniteValue > 0)
			infValue = 10 * maxFiniteValue * nOfElements;
		else
			infValue = 10;
		distMatrixTemp = distMatrix;
		while(distMatrixTemp < distMatrixEnd)
			//if(mxIsInf(*distMatrixTemp++))
			if (MAX_INF != (*distMatrixTemp++))
			{
				*(distMatrixTemp-1) = infValue;
			}
	}
#endif
				
	/* memory allocation */
	//coveredColumns = (bool *)mxCalloc(nOfColumns,  sizeof(bool));
	//coveredRows    = (bool *)mxCalloc(nOfRows,     sizeof(bool));
	//starMatrix     = (bool *)mxCalloc(nOfElements, sizeof(bool));
	//primeMatrix    = (bool *)mxCalloc(nOfElements, sizeof(bool));
	//newStarMatrix  = (bool *)mxCalloc(nOfElements, sizeof(bool)); /* used in step4 */


	coveredColumns = (bool *)malloc(nOfColumns * sizeof(bool));
	memset(coveredColumns, 0, nOfColumns * sizeof(bool));
	coveredRows    = (bool *)malloc(nOfRows * sizeof(bool));
	memset(coveredRows, 0, nOfRows * sizeof(bool));
	starMatrix     = (bool *)malloc(nOfElements * sizeof(bool));
	memset(starMatrix, 0, nOfElements * sizeof(bool));
	primeMatrix    = (bool *)malloc(nOfElements * sizeof(bool));
	memset(primeMatrix, 0, nOfElements * sizeof(bool));
	newStarMatrix  = (bool *)malloc(nOfElements * sizeof(bool)); /* used in step4 */
	memset(newStarMatrix, 0, nOfElements * sizeof(bool));
	/* preliminary steps */
	//�켣���㼣�ľ������ С�� �㼣���켣�ľ������ 
	if(nOfRows <= nOfColumns)
	{
		minDim = nOfRows;
		//�����켣��,����㼣���켣�е���Сֵ��������С���������㼣���Ϊtrue
		for(row=0; row<nOfRows; row++)
		{
			/* find the smallest element in the row */
			//Ѱ��ÿ����Сֵ
			distMatrixTemp = distMatrix + row;
			minValue = *distMatrixTemp;
			//ÿ��ƫ����������Ϊ��nOfRows
			distMatrixTemp += nOfRows;	
			//���ÿ����Сֵ������ <= �У�
			while(distMatrixTemp < distMatrixEnd)
			{
				value = *distMatrixTemp;
				if(value < minValue)
					minValue = value;
				distMatrixTemp += nOfRows;
			}
			
			/* subtract the smallest element from each element of the row */
			//��ȥÿ����Сֵ����ø�������λ��
			distMatrixTemp = distMatrix + row;
			
			while(distMatrixTemp < distMatrixEnd /*&& MAX_INF != minValue*/)
			{
				//if (MAX_INF != (*distMatrixTemp))
				//{
					*distMatrixTemp -= minValue;
				//}
				distMatrixTemp += nOfRows;
			}
		}
		
		/* Steps 1 and 2a */
		//���ÿ�м�ȥ��Сֵ���±���Ϣ
		for(row=0; row<nOfRows; row++)
		{
			for(col=0; col<nOfColumns; col++)
			{
				//printf("%0.3f      ", distMatrix[row + nOfRows*col]);
				if(distMatrix[row + nOfRows*col] == 0)
					if(!coveredColumns[col])
					{
						starMatrix[row + nOfRows*col] = true; //��ʶ��row�е�col����������
						coveredColumns[col]           = true; //��ʶcol�������Ž⣬ע��û�м����е����Ž�
						break;
					}
			}
		//	printf("\n");
		}
			
				
	}
	else /* if(nOfRows > nOfColumns) */ //�켣���㼣�ľ������ ���� �㼣���켣�ľ������ 
	{
		minDim = nOfColumns;
		//�����㼣��������㼣���켣�е���Сֵ��������С��
		for(col=0; col<nOfColumns; col++)
		{
			/* find the smallest element in the column */
			//Ѱ��ÿ����Сֵ
			//��ȡÿ�е��׸�����
			distMatrixTemp = distMatrix     + nOfRows*col;
			//�е����ݳ���
			columnEnd      = distMatrixTemp + nOfRows;
			
			minValue = *distMatrixTemp++;	
			//��ȡ�������Ž��ֵ
			while(distMatrixTemp < columnEnd)
			{
				value = *distMatrixTemp++;
				if(value < minValue)
					minValue = value;
			}
			
			/* subtract the smallest element from each element of the column */
			//ÿ�м�ȥ���Ž��ֵ����ʶ�˸������Ž��λ��
			distMatrixTemp = distMatrix + nOfRows*col;
			while(distMatrixTemp < columnEnd)
			{
				//if (MAX_INF != (*distMatrixTemp))
				//{
					*distMatrixTemp++ -= minValue;
				//}
				//distMatrixTemp++;
			}
		}
		
		/* Steps 1 and 2a */
		//���ÿ�м�ȥ��Сֵ�����Ž�λ��
		for(col=0; col<nOfColumns; col++)
			for(row=0; row<nOfRows; row++)
				if(distMatrix[row + nOfRows*col] == 0)
					if(!coveredRows[row])
					{
						starMatrix[row + nOfRows*col] = true; //row��col������
						coveredColumns[col]           = true; //��ʶcol�������Ž�
						coveredRows[row]              = true; //��ʶrow�������Ž�
						break;
					}
		for(row=0; row<nOfRows; row++)
			coveredRows[row] = false;  //�������о������Ž�
		
	}	
	
	/* move to step 2b */
	step2b(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);

	/* compute cost and remove invalid assignments */
	computeassignmentcost(assignment, cost, distMatrixIn, nOfRows);
	
	/* free allocated memory */
	//mxFree(distMatrix);
	//mxFree(coveredColumns);
	//mxFree(coveredRows);
	//mxFree(starMatrix);
	//mxFree(primeMatrix);
	//mxFree(newStarMatrix);

	free(distMatrix);
	free(coveredColumns);
	free(coveredRows);
	free(starMatrix);
	free(primeMatrix);
	free(newStarMatrix);


	return;
}

/********************************************************/
//����켣����С�ڻ���ڵ㼣���������
void buildassignmentvector(double *assignment, bool *starMatrix, int nOfRows, int nOfColumns)
{
	int row, col;
	//�ҳ��켣���㼣�������±�
	for(row=0; row<nOfRows; row++) //�켣
		for(col=0; col<nOfColumns; col++) //�㼣
			if(starMatrix[row + nOfRows*col]) //�켣���㼣������ֵ
			{
#ifdef ONE_INDEXING
				assignment[row] = col + 1; /* MATLAB-Indexing */
#else
				assignment[row] = col;
#endif
				break;
			}
}

/********************************************************/
void computeassignmentcost(double *assignment, double *cost, double *distMatrix, int nOfRows)
{
	int row, col;
#ifdef CHECK_FOR_INF
	double value;
#endif
	
	for(row=0; row<nOfRows; row++)
	{
#ifdef ONE_INDEXING
		col = assignment[row]-1; /* MATLAB-Indexing */
#else
		col = assignment[row];
#endif

		if(col >= 0)
		{
#ifdef CHECK_FOR_INF
			value = distMatrix[row + nOfRows*col];
			//if(mxIsFinite(value))
			if (MAX_INF == value)
			{
				*cost += value;
			}
			else
#ifdef ONE_INDEXING
				assignment[row] =  0.0;
#else
				assignment[row] = -1.0;
#endif

#else
			*cost += distMatrix[row + nOfRows*col];
#endif
		}
	}
}

/********************************************************/
void step2a(double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{
	bool *starMatrixTemp, *columnEnd;
	int col;
	
	/* cover every column containing a starred zero */
	for(col=0; col<nOfColumns; col++)
	{
		starMatrixTemp = starMatrix     + nOfRows*col;
		columnEnd      = starMatrixTemp + nOfRows;
		while(starMatrixTemp < columnEnd){
			if(*starMatrixTemp++)
			{
				coveredColumns[col] = true;
				break;
			}
		}	
	}

	/* move to step 3 */
	step2b(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

/********************************************************/
void step2b(double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{
	int col, nOfCoveredColumns;
	
	/* count covered columns */
	nOfCoveredColumns = 0;
	//�켣���㼣����ֵ�ĸ���
	for(col=0; col<nOfColumns; col++)
		if(coveredColumns[col])		//�����Ƿ������Ž�
			nOfCoveredColumns++;
	//�켣����С�ڻ���ڵ㼣����	
	if(nOfCoveredColumns == minDim)  
	{
		/* algorithm finished */
		buildassignmentvector(assignment, starMatrix, nOfRows, nOfColumns);
	}
	else
	{
		/* move to step 3 */
		step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
	}
	
}

/********************************************************/
//����켣������㼣����������±��ʶ��
//1���ҵ��㼣��û�����Ž�ĵ㼣
//2���ҵ��õ㼣���켣��������Ž�
//3���Ƿ�Ϊ��һ�����Ž�
//�ǣ�step4����
//���ǣ�
//   a�����������ű�ʶ�����ҵ�״̬��
//   b������һ�����Ž�ĵ㼣����Ϊû���ҵ���
//	 c�����±�����֪�����õ㼣��Ϊ��һ���Ž⣻����step4��
void step3(double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{
	bool zerosFound;
	int row, col, starCol;

	zerosFound = true;
	while(zerosFound)
	{
		zerosFound = false;	
		//�����㼣��Ϣ
		for(col=0; col<nOfColumns; col++)
		{
			if(!coveredColumns[col]) //�ҳ��ĸ���û�����Ž�
			{
				for(row=0; row<nOfRows; row++) //�����켣��Ϣ
				{
					if((!coveredRows[row]) && (distMatrix[row + nOfRows*col] == 0)) //����켣���㼣�����Ž�
					{
						/* prime zero */
						primeMatrix[row + nOfRows*col] = true; //row�켣���㼣�����Ž�

						/* find starred zero in current row */
						for(starCol=0; starCol<nOfColumns; starCol++) //�����㼣��Ϣ,�ҵ���һ�����Ž��±꣺���Ǹõ㼣�������켣�е����Ž�
						{
							if(starMatrix[row + nOfRows*starCol]) //��һ�����Ž�
								break;
						}
							

						if(starCol == nOfColumns) /* no starred zero found */
						{
							/* move to step 4 */
							step4(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim, row, col);
							return;
						}
						else
						{
							coveredRows[row]        = true;  //�켣�ҵ����Ž�
							coveredColumns[starCol] = false; //���õ㼣�������켣�����Ž�����Ϊfalse
							zerosFound              = true;
							break;
						}
					}
				}
			}
		}
			
	}
	
	/* move to step 5 */
	step5(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

/********************************************************/

void step4(double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim, int row, int col)
{	
	int n, starRow, starCol, primeRow, primeCol;
	int nOfElements = nOfRows*nOfColumns;
	
	/* generate temporary copy of starMatrix */
	//newStarMatrix : �㼣���켣�����Ž�
	for(n=0; n<nOfElements; n++)
		newStarMatrix[n] = starMatrix[n];
	
	/* star current zero */
	//���������������
	newStarMatrix[row + nOfRows*col] = true;

	/* find starred zero in current column */
	starCol = col;
	//�����켣���ҵ��õ�㼣���켣�����Ž⣬��ʷ����
	for(starRow=0; starRow<nOfRows; starRow++)
		if(starMatrix[starRow + nOfRows*starCol])
			break;
	//�ҵ���ʷ���Ź켣starRow
	while(starRow<nOfRows)
	{
		/* unstar the starred zero */
		newStarMatrix[starRow + nOfRows*starCol] = false;
	
		/* find primed zero in current row */
		//�켣���ŵĵ㼣��primeCol
		primeRow = starRow;
		for(primeCol=0; primeCol<nOfColumns; primeCol++)
			if(primeMatrix[primeRow + nOfRows*primeCol])
				break;
								
		/* star the primed zero */
		newStarMatrix[primeRow + nOfRows*primeCol] = true;
	
		/* find starred zero in current column */
		starCol = primeCol;
		for(starRow=0; starRow<nOfRows; starRow++)
			if(starMatrix[starRow + nOfRows*starCol])
				break;
	}	

	/* use temporary copy as new starMatrix */
	/* delete all primes, uncover all rows */
	for(n=0; n<nOfElements; n++)
	{
		primeMatrix[n] = false;
		starMatrix[n]  = newStarMatrix[n];
	}
	for(n=0; n<nOfRows; n++)
		coveredRows[n] = false;
	
	/* move to step 2a */
	step2a(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

/********************************************************/
//����켣���㼣���Ž���Ϣ:
//1����ȡδ�ҵ����Ž�Ĺ켣������Щ�켣���ҵ��㼣���켣��Сֵh�� �켣�е����Ž�
//2����ȡ�ҵ����Ž�Ĺ켣�����켣����Щ��ľ����� + h�� �������Ž� + h;����δ�ҵ��켣�ĵ㼣��ȫ��+h��
//3����ȡδƥ�䵽���Ź켣�ĵ㼣�����õ㼣�����й켣�ľ����� - h�� 
//��ȡ���켣���㼣�Լ��㼣���켣�����Ž��������ʱ������Ĺ켣rowΪδ�ҵ�״̬����step3�д���
void step5(double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{
	double h, value;
	int row, col;
	
	/* find smallest uncovered element h */
	h = /*mxGetInf();	*/MAX_INF;
	//�����켣���㼣�ľ�����ҳ�û�����Ž�켣������ȡ�켣���㼣����Сֵ
	for(row=0; row<nOfRows; row++)
	{
		//�Ƿ�����·��ƥ���ϸù켣
		if(!coveredRows[row]) 
		{//��
			for(col=0; col<nOfColumns; col++) //�����㼣��Ϣ
			{
				//�õ㼣�Ƿ��ѱ�·��ƥ����
				if(!coveredColumns[col])
				{//��
					//��ȡû��ƥ��·���ĵ㼣����
					value = distMatrix[row + nOfRows*col];
					//�ҳ���С·������h
					if(value < h)
					{
						h = value;
					}
				}
			}
				
		}
			
	}
		
	
	/* add h to each covered row */
	//�����켣���㼣�ľ���������ҵ��켣���㼣������ֵ���У�������Сֵ
	for(row=0; row<nOfRows; row++)
	{
		//�ҳ��ѱ�ƥ��㼣�Ĺ켣
		if(coveredRows[row])
		{
			//�㼣��Ϣ����
			for(col=0; col<nOfColumns; col++)
			{
				//����Сֵ�ӵ���ƥ��·���Ĺ켣��
				distMatrix[row + nOfRows*col] += h;
			}
		}
			
	}
		
	
	/* subtract h from each uncovered column */
	//�����㼣���켣�ľ�����ҳ��㼣δƥ�����ŵ�����������й켣���õ㼣����ȥ��Сֵ
	for(col=0; col<nOfColumns; col++)
	{
		if(!coveredColumns[col])
		{
			for(row=0; row<nOfRows; row++)
			{
				distMatrix[row + nOfRows*col] -= h;
			}
		}
	}
		
	
	/* move to step 3 */
	step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

