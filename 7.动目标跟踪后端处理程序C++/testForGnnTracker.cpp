// testForGnnTracker.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
//#include <stdlib.h>
#include <malloc.h>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include "GNNTrakerCom.h"

using namespace std;
using namespace Eigen;

int main()
{

	//Matrix4d aa;
	//EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	int nMaxLostLen = 5;
	trackerData *pLostData = new trackerData();

	fstream fileReador;
	
	fileReador.open(".\\trackerPoint.txt", ios::in);

	if (!fileReador.is_open())
	{
		printf("\nfile open file, please check the file is exist or the name is trackerPoint.txt?\n");
		return 0;
	}
	//char readBuf[10240] = {0};
	char *readBuf = (char *)malloc(sizeof(char) * 10240);
	trackerData *trakcerPoins = (trackerData *)malloc(sizeof(trackerData) * 1024);
	int nPointsLen = 0;
	int nReadLen = 0;
	std::string strData = "";
	int nLoopDataLen = 0;
	int nTotalLoop = 0;
	//GNN_Traker *pTraker = createGNNTracker();

	GNN_Traker *pTraker = (GNN_Traker *)malloc(sizeof(GNN_Traker));
	pTraker = createGNNTracker();
	((GNN_Traker *)pTraker)->setTrackerInfo(4, 3);
	int nLostLen = 0;
	while(!fileReador.eof())
	{
		fileReador.getline(readBuf, 1000);
		nReadLen = strlen(readBuf);
		nPointsLen = 0;
		nLoopDataLen = 0;
		for (int i = 0; i < nReadLen; i++)
		{
			//数据起始位置
			if ('[' == readBuf[i])
			{
				strData.clear();

				continue;
			}
			else if (']' == readBuf[i])
			{
				trakcerPoins[nLoopDataLen].dTargetY = atof(strData.c_str());
				trakcerPoins[nLoopDataLen].nOutId = TARGET_NULL;
				nLoopDataLen++;
				strData.clear();
				break;
			}
			//当个数据结束位置
			else if (',' == readBuf[i])
			{
				trakcerPoins[nLoopDataLen].dTragetX = atof(strData.c_str());
				strData.clear();
			}
			else if (';' == readBuf[i])
			{
				trakcerPoins[nLoopDataLen].dTargetY = atof(strData.c_str());
				trakcerPoins[nLoopDataLen].nOutId = TARGET_NULL;
				nLoopDataLen++;
				strData.clear();
			}
			else 
			{
				strData.append(1, readBuf[i]);
			}
		}
		//printf("\b\b\b\b\b\b\b\b\b\b");

		pTraker->GNNTracker_Sequt(nTotalLoop, trakcerPoins, nLoopDataLen, nLostLen, pLostData, nMaxLostLen);
		nTotalLoop++;
		//printf("Runing for % 2d frame\n", nTotalLoop);
		

	}

	free(pTraker);
	delete(pLostData);
	free(readBuf);
	free(trakcerPoins);
	fileReador.close();

	return 0;
}

