
/*
function [assignment, cost] = assignmentoptimal(distMatrix)
*/

//#include <mex.h>
//#include <matrix.h>

//这个代码完全是为了适应matlab写的
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

//nOfRows 是 columnsLen , nOfColumns 是 RowsLen
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
	//轨迹到点迹的距离个数 小于 点迹到轨迹的距离个数 
	if(nOfRows <= nOfColumns)
	{
		minDim = nOfRows;
		//遍历轨迹数,求出点迹到轨迹中的最小值（列数最小），并将点迹标记为true
		for(row=0; row<nOfRows; row++)
		{
			/* find the smallest element in the row */
			//寻找每列最小值
			distMatrixTemp = distMatrix + row;
			minValue = *distMatrixTemp;
			//每列偏移量即行数为：nOfRows
			distMatrixTemp += nOfRows;	
			//求出每列最小值，（列 <= 行）
			while(distMatrixTemp < distMatrixEnd)
			{
				value = *distMatrixTemp;
				if(value < minValue)
					minValue = value;
				distMatrixTemp += nOfRows;
			}
			
			/* subtract the smallest element from each element of the row */
			//减去每列最小值，获得该列最优位置
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
		//标记每列减去最小值的下标信息
		for(row=0; row<nOfRows; row++)
		{
			for(col=0; col<nOfColumns; col++)
			{
				//printf("%0.3f      ", distMatrix[row + nOfRows*col]);
				if(distMatrix[row + nOfRows*col] == 0)
					if(!coveredColumns[col])
					{
						starMatrix[row + nOfRows*col] = true; //标识第row列第col行数据最优
						coveredColumns[col]           = true; //标识col行有最优解，注：没有计算行的最优解
						break;
					}
			}
		//	printf("\n");
		}
			
				
	}
	else /* if(nOfRows > nOfColumns) */ //轨迹到点迹的距离个数 大于 点迹到轨迹的距离个数 
	{
		minDim = nOfColumns;
		//遍历点迹数，求出点迹到轨迹中的最小值（行数最小）
		for(col=0; col<nOfColumns; col++)
		{
			/* find the smallest element in the column */
			//寻找每行最小值
			//获取每行的首个数据
			distMatrixTemp = distMatrix     + nOfRows*col;
			//行的数据长度
			columnEnd      = distMatrixTemp + nOfRows;
			
			minValue = *distMatrixTemp++;	
			//获取行中最优解的值
			while(distMatrixTemp < columnEnd)
			{
				value = *distMatrixTemp++;
				if(value < minValue)
					minValue = value;
			}
			
			/* subtract the smallest element from each element of the column */
			//每行减去最优解的值，标识了该行最优解的位置
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
		//标记每行减去最小值的最优解位置
		for(col=0; col<nOfColumns; col++)
			for(row=0; row<nOfRows; row++)
				if(distMatrix[row + nOfRows*col] == 0)
					if(!coveredRows[row])
					{
						starMatrix[row + nOfRows*col] = true; //row列col行最优
						coveredColumns[col]           = true; //标识col行有最优解
						coveredRows[row]              = true; //标识row列有最优解
						break;
					}
		for(row=0; row<nOfRows; row++)
			coveredRows[row] = false;  //所有列中均无最优解
		
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
//处理轨迹个数小于或等于点迹个数的情况
void buildassignmentvector(double *assignment, bool *starMatrix, int nOfRows, int nOfColumns)
{
	int row, col;
	//找出轨迹到点迹的最优下标
	for(row=0; row<nOfRows; row++) //轨迹
		for(col=0; col<nOfColumns; col++) //点迹
			if(starMatrix[row + nOfRows*col]) //轨迹到点迹的最优值
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
	//轨迹到点迹最优值的个数
	for(col=0; col<nOfColumns; col++)
		if(coveredColumns[col])		//该行是否有最优解
			nOfCoveredColumns++;
	//轨迹个数小于或等于点迹个数	
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
//处理轨迹最优与点迹最优情况，下标标识：
//1、找到点迹中没有最优解的点迹
//2、找到该点迹到轨迹距离的最优解
//3、是否为第一个最优解
//是：step4处理
//不是：
//   a、将该行最优标识设置找到状态；
//   b、将第一个最优解的点迹设置为没有找到；
//	 c、重新遍历，知道将该点迹变为第一最优解；进入step4；
void step3(double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{
	bool zerosFound;
	int row, col, starCol;

	zerosFound = true;
	while(zerosFound)
	{
		zerosFound = false;	
		//遍历点迹信息
		for(col=0; col<nOfColumns; col++)
		{
			if(!coveredColumns[col]) //找出哪个点没有最优解
			{
				for(row=0; row<nOfRows; row++) //遍历轨迹信息
				{
					if((!coveredRows[row]) && (distMatrix[row + nOfRows*col] == 0)) //处理轨迹到点迹的最优解
					{
						/* prime zero */
						primeMatrix[row + nOfRows*col] = true; //row轨迹到点迹的最优解

						/* find starred zero in current row */
						for(starCol=0; starCol<nOfColumns; starCol++) //遍历点迹信息,找到第一个最优解下标：就是该点迹到其他轨迹中的最优解
						{
							if(starMatrix[row + nOfRows*starCol]) //第一个最优解
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
							coveredRows[row]        = true;  //轨迹找到最优解
							coveredColumns[starCol] = false; //将该点迹到其他轨迹的最优解设置为false
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
	//newStarMatrix : 点迹到轨迹的最优解
	for(n=0; n<nOfElements; n++)
		newStarMatrix[n] = starMatrix[n];
	
	/* star current zero */
	//将处理点设置最优
	newStarMatrix[row + nOfRows*col] = true;

	/* find starred zero in current column */
	starCol = col;
	//遍历轨迹，找到该点点迹到轨迹的最优解，历史最优
	for(starRow=0; starRow<nOfRows; starRow++)
		if(starMatrix[starRow + nOfRows*starCol])
			break;
	//找到历史最优轨迹starRow
	while(starRow<nOfRows)
	{
		/* unstar the starred zero */
		newStarMatrix[starRow + nOfRows*starCol] = false;
	
		/* find primed zero in current row */
		//轨迹最优的点迹：primeCol
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
//处理轨迹到点迹最优解信息:
//1、获取未找到最优解的轨迹，在这些轨迹中找到点迹到轨迹最小值h； 轨迹中的最优解
//2、获取找到最优解的轨迹，将轨迹到这些点的距离项 + h； 已有最优解 + h;（除未找到轨迹的点迹，全部+h）
//3、获取未匹配到最优轨迹的点迹，将该点迹到所有轨迹的距离项 - h； 
//获取，轨迹到点迹以及点迹到轨迹的最优解情况，此时待处理的轨迹row为未找到状态，到step3中处理
void step5(double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{
	double h, value;
	int row, col;
	
	/* find smallest uncovered element h */
	h = /*mxGetInf();	*/MAX_INF;
	//遍历轨迹到点迹的距离项，找出没有最优解轨迹，并获取轨迹到点迹的最小值
	for(row=0; row<nOfRows; row++)
	{
		//是否已有路径匹配上该轨迹
		if(!coveredRows[row]) 
		{//无
			for(col=0; col<nOfColumns; col++) //遍历点迹信息
			{
				//该点迹是否已被路径匹配上
				if(!coveredColumns[col])
				{//否
					//获取没有匹配路径的点迹距离
					value = distMatrix[row + nOfRows*col];
					//找出最小路劲距离h
					if(value < h)
					{
						h = value;
					}
				}
			}
				
		}
			
	}
		
	
	/* add h to each covered row */
	//遍历轨迹到点迹的距离项，将已找到轨迹到点迹的最优值的行，加上最小值
	for(row=0; row<nOfRows; row++)
	{
		//找出已被匹配点迹的轨迹
		if(coveredRows[row])
		{
			//点迹信息遍历
			for(col=0; col<nOfColumns; col++)
			{
				//将最小值加到以匹配路径的轨迹上
				distMatrix[row + nOfRows*col] += h;
			}
		}
			
	}
		
	
	/* subtract h from each uncovered column */
	//遍历点迹到轨迹的距离项，找出点迹未匹配最优的数据项，将所有轨迹到该点迹均减去最小值
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

