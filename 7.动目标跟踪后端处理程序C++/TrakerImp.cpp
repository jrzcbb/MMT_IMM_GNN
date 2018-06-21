//#include "stdafx.h"
#include "stdafx.h"
#include "TrakerImp.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <Eigen/Dense>
#include <iostream>
using namespace std;
using namespace Eigen;
#define NO_PRINT

TrakerImp::TrakerImp(void)
{
	m_nFirsLoopId = -1;
	m_nMax_linking_distance = 8;
	m_nMax_gap_closing = 4;
	m_nMaxTargetId = 1;
	m_pMapTraker = NULL;
	m_nThreshold = m_nMax_linking_distance * m_nMax_linking_distance;
	//Alpha Beta 
	m_dAlpha = 0.3;
	m_dBeta  = 0.2;
	m_nHisDataLen = 0;
}

TrakerImp::~TrakerImp(void)
{
	//内存释放。。。
	TarcksMap *pDel = NULL;
	while(NULL != m_pMapTraker)
	{
		pDel = m_pMapTraker;
		m_pMapTraker = m_pMapTraker->pNext;
		delete pDel->pTaget;
		delete pDel;
	}
}

void TrakerImp::clearTrakers()
{
	TarcksMap *pDel = NULL;
	while (NULL != m_pMapTraker)
	{
		pDel = m_pMapTraker;
		m_pMapTraker = m_pMapTraker->pNext;
		delete pDel->pTaget;
		delete pDel;
	}
	m_nFirsLoopId = -1;
}

void TrakerImp::insertNewTraker(int nLoopId, trackerData *pData, int nLengTh)
{
	if (NULL == pData || 0 >= nLengTh)
	{
		return;
	}
	targetTracker *pNewTrackers = NULL;
	for (int nTarckerIndex = 0; nTarckerIndex < nLengTh; nTarckerIndex++)
	{
		pNewTrackers = createTracker(nLoopId, pData + nTarckerIndex);
		if (NULL == pNewTrackers)
		{
			printf("error: createTracker fail!\n");
			continue;
		}
		insertToMapTrackerToTail(pNewTrackers);
	}
}

void TrakerImp::setTrackerInfo(int nMax_linking_distance, int nMax_gap_closing)
{
	m_nMax_gap_closing = nMax_gap_closing;
	m_nMax_linking_distance = nMax_linking_distance;
	m_nThreshold = m_nMax_linking_distance * m_nMax_linking_distance;
}

void TrakerImp::GNNTracker_Sequt(int nTotalLoop, trackerData *pTrackData, int nDataLen, int &nLoseTrackerLen, trackerData *pLostData, int nMaxLostLen)
{
	if (0 >= nDataLen)
	{
		return;
	}
	if(0 > m_nFirsLoopId)
	{
		m_nFirsLoopId = nTotalLoop;
	}
	//
	if (m_nFirsLoopId == nTotalLoop)
	{
		insertNewTraker(nTotalLoop, pTrackData, nDataLen);
		
#ifndef NO_PRINT
		//测试打印
		printf("\n******************** %d the end tracker  First == *******************\n", nTotalLoop);
		TarcksMap *pMap = m_pMapTraker;
		targetTracker *pHisTracker = NULL;
		while (NULL != pMap)
		{
			pHisTracker = pMap->pTaget;
			printf("nLoopId = %d, dx = %0.3f, dy = %0.3f, speed = %0.3f, ange = %0.3f, Id = %d \n", pHisTracker->nLoopId,
				pHisTracker->dTragetX, pHisTracker->dTargetY,
				pHisTracker->dSpeed, pHisTracker->dDirection, pHisTracker->nOutId);
			pMap = pMap->pNext;
		}
		printf("********************end*******************\n");
#endif //

		return;
	}
	calYuceTracker(nTotalLoop);
	int nRowLen = 0;
	trackerData *pHisGjData = calTrackerPointByLoop(nTotalLoop - 1, nRowLen);

	double *pTrackMap = NULL;
	double dcost = 0;
	targetTrackerEx *pTracker = NULL;
	trackerData *pMatchingData = NULL;
	if (NULL != pHisGjData)
	{
		pTrackMap = createMunkresData(pHisGjData, nRowLen, pTrackData, nDataLen);
		if (NULL != pTrackMap)
		{
			memset(m_dTrackerRows, 0, sizeof(double) * nDataLen);
			assignmentoptimal(m_dTrackerRows, &dcost, pTrackMap, nDataLen, nRowLen);
#ifndef NO_PRINT
			printf("the xyl data is : ");
			for (int i = 0; i < nDataLen; i++)
			{
				printf("%d  ", (int)m_dTrackerRows[i]);
			}
			printf("\n");
#endif
			//m_dTrackerRows为i列(点迹) m_dTrackerrows[i]行（轨迹），
			for (int i = 0; i < nDataLen; i++)
			{
				if (0 == m_dTrackerRows[i])
				{
					continue;
				}
				if (/*1 == m_dTrackerRows[i] &&*/ MAX_INF == (*(pTrackMap + i +  nDataLen * ((int)m_dTrackerRows[i] - 1))))
				{
					m_dTrackerRows[i] = 0;
					continue;
				}
				//pTracker = getTrackerById(nTotalLoop - 1, pHisGjData[(int)m_dTrackerRows[i] - 1].nOutId);
				pTracker = getTrackerById(nTotalLoop - 1, pHisGjData[(int)m_dTrackerRows[i] - 1].nOutId);
				assert(NULL != pTracker && "the cal x y l fail!\n");
				pMatchingData = pTrackData + i;
				pMatchingData->nOutId = pTracker->nOutId;
				if (1 == pTracker->nPtLen)
				{
					calSpeedAndDirect(pTracker, pMatchingData, nTotalLoop);
				}
				else
				{
					calAlphaAndBeta(nTotalLoop, pTracker, pMatchingData);
				}
				
				updateMapData(nTotalLoop, pTracker, pMatchingData);
			}
		}
	}

	//计算最短路径：每个数据最多匹配计算m_nMax_linking_distance次；
	//历史最小loop
	int nHisMinLoopId = nTotalLoop - m_nMax_gap_closing;
	if (0 > nHisMinLoopId)
	{
		nHisMinLoopId = 0;
	}
	int nNextTrackId = nTotalLoop - 2;

	int nMatchingId = 0;
	bool bMatching = false;


	double dMinRange = 100000;
	trackerData *pTrackerToPoint = NULL;
	//遍历待匹配目标,计算点迹到轨迹的最优的遍历	
	for (int nTrackerId = 0; nTrackerId < nDataLen; nTrackerId++)
	{
		pMatchingData = (pTrackData + nTrackerId);
		//该点已经匹配最优路径，过滤在轨迹到点迹中已经匹配了的点
		if (TARGET_NULL != pMatchingData->nOutId)
		{
			continue;
		}
		pTracker = NULL;
		dMinRange = 100000;
		bMatching = false;
		nNextTrackId = nTotalLoop - 2;
		//最多匹配计算m_nMax_linking_distance次
		for (; nNextTrackId >= nHisMinLoopId; nNextTrackId--)
		{
			//相邻数据路径匹配,新点迹到轨迹的距离
			pTracker = (targetTrackerEx *)matchingPointToTracker(nNextTrackId, pMatchingData, nMatchingId, dMinRange);
			if(NULL != pTracker)
			{
				//遍历占用
				pMatchingData->nOutId = TARGET_OCCUPY;
				bMatching = true;
				//遍历轨迹目标匹配最优距离
				pTrackerToPoint = matchingTrackerToPoint(pTracker, dMinRange, pTrackData, nDataLen);
				pMatchingData->nOutId = TARGET_NULL;
				if (NULL != pTrackerToPoint)
				{
					pMatchingData = pTrackerToPoint;
					pMatchingData->nOutId = TARGET_OCCUPY;
					nTrackerId--;
				}
				break;
			}
		}
		//未找到匹配路径
		if (false == bMatching)
		{
			insertNewTraker(nTotalLoop, pMatchingData, 1);
		}
		else
		{
			pMatchingData->nOutId = nMatchingId;
			pTracker = getTrackerById(pTracker->nLoopId, pTracker->nOutId);
			assert(NULL != pTracker);
			//更新数据
			if (1 == pTracker->nPtLen)
			{
				calSpeedAndDirect(pTracker, pMatchingData, nTotalLoop);
			}
			else
			{
				calAlphaAndBeta(nTotalLoop, pTracker, pMatchingData);
			}
			updateMapData(nTotalLoop, pTracker, pMatchingData);
		}
	}
	nLoseTrackerLen = 0;
	deleteTracke(nHisMinLoopId, pLostData, nMaxLostLen, nLoseTrackerLen);

#ifndef NO_PRINT
	//测试打印
	printf("\n******************** %d the end tracker*******************\n", nTotalLoop);
	TarcksMap *pMap = m_pMapTraker;
	targetTracker *pHisTracker = NULL;
	while (NULL != pMap)
	{
		pHisTracker = pMap->pTaget;
		printf("nLoopId = %d, dx = %0.3f, dy = %0.3f, speed = %0.3f, ange = %0.3f, Id = %d \n", pHisTracker->nLoopId,
			pHisTracker->dTragetX, pHisTracker->dTargetY,
			pHisTracker->dSpeed, pHisTracker->dDirection, pHisTracker->nOutId);
		pMap = pMap->pNext;
	}
	printf("********************end*******************\n");
#endif
}

void TrakerImp::releaseTraker()
{
	delete this;
}

void TrakerImp::insertToMapTrakerToHeador(targetTracker *pData)
{
	if (NULL == pData)
	{
		return;
	}
	//printf("**********insert new gj data ****************\n");
	//printf("nLoopId = %d, x = %0.3f		y = %0.3f nOutId = %d\n", pData->nLoopId,
	//	pData->dTragetX, pData->dTargetY, pData->nOutId);

	if (NULL == m_pMapTraker)
	{
		m_pMapTraker = new TarcksMap;
		m_pMapTraker->pNext = NULL;
		m_pMapTraker->pTaget = pData;
		m_pTailTraker = m_pMapTraker;
		return;
	}
	TarcksMap *pNewMap = new TarcksMap;
	pNewMap->pNext = m_pMapTraker;
	pNewMap->pTaget = pData;
	m_pMapTraker = pNewMap;
}

void TrakerImp::insertToMapTrackerToTail(targetTracker *pData)
{
	if (NULL == pData)
	{
		return;
	}
	//printf("**********insert new gj data ****************\n");
	//printf("nLoopId = %d, x = %0.3f		y = %0.3f nOutId = %d\n", pData->nLoopId,
	//	pData->dTragetX, pData->dTargetY, pData->nOutId);

	if (NULL == m_pMapTraker)
	{
		printf("m_pMapTraker is NULL.....\n");
		m_pMapTraker = new TarcksMap;
		m_pMapTraker->pNext = NULL;
		m_pMapTraker->pTaget = pData;
		m_pTailTraker = m_pMapTraker;
		m_pTailTraker->nId = pData->nOutId;
		return;
	}
	TarcksMap *pNewMap = new TarcksMap;
	pNewMap->pNext = NULL;
	pNewMap->pTaget = pData;
	pNewMap->nId = pData->nOutId;
	m_pTailTraker->pNext = pNewMap;
	m_pTailTraker= pNewMap;
}

targetTrackerEx * TrakerImp::createTracker(int &nLoopId, trackerData *pData)
{
	targetTrackerEx *pNewTracker = NULL;
	if (NULL == pData)
	{
		return NULL;
	}
	pNewTracker = new targetTrackerEx;
	pNewTracker->nPtLen = 1;
	pNewTracker->dTargetY = pData->dTargetY;
	pNewTracker->dTragetX = pData->dTragetX;
	pNewTracker->dYcTargetY = -1;
	pNewTracker->dYcTragetX = -1;
	pNewTracker->fAng = pData->fAng;
	pNewTracker->nOutId = ++m_nMaxTargetId;
	pNewTracker->dDirection = 0;
	pNewTracker->dSpeed = 0;
	pData->nOutId = pNewTracker->nOutId;
	pNewTracker->nLoopId = nLoopId;
	pNewTracker->pNextTrack = NULL;

	//IMM 与 Kalman专区
	pNewTracker->xPred << pData->dTragetX, 0, pData->dTargetY, 0;
	pNewTracker->xEstm << pData->dTragetX, 0, pData->dTargetY, 0;
	pNewTracker->PosState << pNewTracker->xPred,pNewTracker->xPred,pNewTracker->xPred;

	return pNewTracker;
}

targetTrackerEx * TrakerImp::getTrackerById(int nLoop, int nId)
{
	TarcksMap *pMap = m_pMapTraker;
	targetTrackerEx *pRstTracker = NULL;
	while(NULL != pMap)
	{
		if (nLoop == pMap->pTaget->nLoopId && nId == pMap->pTaget->nOutId)
		{
			pRstTracker = (targetTrackerEx *)pMap->pTaget;
			break;
		}
		pMap->pTaget;
		pMap = pMap->pNext;
	}
	return pRstTracker;
}

targetTracker *TrakerImp::matchingPointToTracker(int nMatchLoopId, trackerData *pData, int &nMatchId, double &dMinRange)
{
	double dRstMinRange = 100000;
	double dTempRange = 0;
	double dX = 0.0;
	double dY = 0.0;
	targetTracker *pRstTracker = NULL;
	TarcksMap *pMap = m_pMapTraker;
	targetTrackerEx *pHisTracker = NULL;
	while (NULL != pMap)
	{
		pHisTracker = (targetTrackerEx *)pMap->pTaget;
		if (nMatchLoopId != pHisTracker->nLoopId)
		{
			pMap = pMap->pNext;
			continue;
		}
		dX = pHisTracker->dYcTragetX - pData->dTragetX;
		dY = pHisTracker->dYcTargetY - pData->dTargetY;
		dTempRange = dX * dX + dY * dY;
		if (m_nThreshold < dTempRange)
		{
			pMap = pMap->pNext;
			continue;
		}
		dTempRange = sqrt(dTempRange);
		if (dRstMinRange > dTempRange)
		{
			dRstMinRange = dTempRange;
			nMatchId = pHisTracker->nOutId;
			pRstTracker = pHisTracker;
		}
		pMap = pMap->pNext;
	}
	dMinRange = dRstMinRange;
	return pRstTracker;
}

//点迹到轨迹存在最优距离，求出该轨迹点到处点迹点之外的所有未匹配点的轨迹
trackerData *TrakerImp::matchingTrackerToPoint(targetTracker *pTracker, double dMinRange, trackerData *pTrackData, int nDataLen)
{
	trackerData *pRstData = NULL;
	trackerData *pMatchingData = NULL;
	double dX = ((targetTrackerEx *)pTracker)->dYcTragetX;
	double dY = ((targetTrackerEx *)pTracker)->dYcTargetY;
	double dRstX = 0.0;
	double dRstY = 0.0;
	double dTempRange = 0;
	//遍历点迹
	for (int nPointOffset = 0; nPointOffset < nDataLen; nPointOffset++)
	{
		pMatchingData = pTrackData + nPointOffset;
		//过滤已分配轨迹的点迹，以及占用的点迹
		if (TARGET_NULL != pMatchingData->nOutId || TARGET_OCCUPY == pMatchingData->nOutId)
		{
			continue;
		}
		dRstX = dX - pMatchingData->dTragetX;
		dRstY = dY - pMatchingData->dTargetY;
		dTempRange = dRstX * dRstX + dRstY * dRstY;
		if (m_nThreshold < dTempRange)
		{
			continue;
		}
		dTempRange = sqrt(dTempRange);
		if (dMinRange > dTempRange)
		{
			pRstData = pMatchingData;
			dMinRange = dTempRange;
		}
	}
	return pRstData;
}

void TrakerImp::updateMapData(int nLoopId, targetTracker *pTracker, trackerData *pUpdateData)
{
	//printf("\n************change latest tracke***************\n");
	//printf("history to current data : \n");
	//printf("nHisLoop = %d	nLoop = %d\n", pTracker->nLoopId, nLoopId);
	//
	//printf("pHisOutId = %d	OutId = %d \n", pTracker->nOutId, pUpdateData->nOutId);
	//
	//printf("(x = %0.3f	  y = %0.3f)=>(x = %0.3f	  y = %0.3f)\n", pTracker->dTragetX, pTracker->dTargetY, 
	//	pUpdateData->dTragetX, pUpdateData->dTargetY);
	
	pTracker->nLoopId = nLoopId;
	pTracker->nPtLen++;
	pTracker->nOutId = pUpdateData->nOutId;
	pTracker->dTragetX = pUpdateData->dTragetX;
	pTracker->dTargetY = pUpdateData->dTargetY;
	pTracker->dSpeed = pUpdateData->dSpeed;
	pTracker->dDirection = pUpdateData->dDirection;
}

void TrakerImp::deleteTracke(int nMinId, trackerData *pLostData, int nMaxLen, int &nRealLostLen)
{
	//m_vcLoseIds.clear();
	TarcksMap *pMap = m_pMapTraker;
	TarcksMap *pPreMap = pMap;
	TarcksMap *pDeleteMap = NULL;
	int nDeletLen = 0;
	int nTalLen = 0;
	nRealLostLen = 0;
	while(NULL != pMap)
	{
		nTalLen++;
		if (nMinId > pMap->pTaget->nLoopId)
		{
			//删除Map
			//m_vcLoseIds.push_back(pMap->pTaget->nOutId);
			pDeleteMap = pMap;
			pMap = pMap->pNext;
			if (m_pMapTraker == pDeleteMap)
			{
				m_pMapTraker = pMap;
				pPreMap = m_pMapTraker;
			}
			else
			{
				pPreMap->pNext = pDeleteMap->pNext;
			}
			if (nDeletLen < nMaxLen)
			{
				cpyTrackerData(pDeleteMap->pTaget, pLostData + nRealLostLen);
				nRealLostLen++;
			}
			delete pDeleteMap->pTaget;
			delete pDeleteMap;
			nDeletLen++;
			continue;
		}
		pPreMap = pMap;
		pMap = pMap->pNext;
	}
	if (NULL == m_pMapTraker)
	{
		m_nFirsLoopId = -1;
	}
	m_pTailTraker = pPreMap;
}

void TrakerImp::getLoseTrackers(int *pTargetIds, int nLoseTrackerLen)
{
	nLoseTrackerLen = -1;
	if (NULL == pTargetIds)
	{
		return;
	}
	int nLen = m_vcLoseIds.size();
	for (int i = 0; i < nLoseTrackerLen; i++)
	{
		if (i < nLen)
		{
			pTargetIds[i] = m_vcLoseIds.at(i);
		}
		else
			break;
	}
}


trackerData * TrakerImp::calTrackerPointByLoop(int nLoopId, int &nLen)
{
	nLen = 0;
	
	TarcksMap *pMap = m_pMapTraker;
	targetTrackerEx *pDealTracker = NULL;
	while (NULL != pMap)
	{
		pDealTracker = (targetTrackerEx *)pMap->pTaget;
		if (nLoopId == pDealTracker->nLoopId)
		{
			m_trackersPoint[nLen].dTargetY = pDealTracker->dYcTargetY;
			m_trackersPoint[nLen].dTragetX = pDealTracker->dYcTragetX;
			m_trackersPoint[nLen].nOutId = pDealTracker->nOutId;
			nLen++;
		}
		pMap = pMap->pNext;
	}
	if (0 >= nLen)
	{
		return NULL;
	}
	return m_trackersPoint;
}

double * TrakerImp::createMunkresData(trackerData *pRowData, int nRowLen, trackerData *pColData, int nColLen)
{
	//printf("\n*****************gj****************\n");
	//for (int i = 0; i < nRowLen; i++)
	//{
	//	printf("%0.3f  %0.3f\n", (pRowData + i)->dTragetX, (pRowData + i)->dTargetY );
	//}
	//printf("\n*****************dj****************\n");
	//for (int i = 0; i < nColLen; i++)
	//{
	//	printf("%0.3f  %0.3f\n", (pColData + i)->dTragetX, (pColData + i)->dTargetY );
	//}

	double *pRstData = NULL;
	if (NULL == pRowData || NULL == pColData || 0 >= nRowLen || 0 >= nColLen)
	{
		return pRstData;
	}
	//轨迹到点迹的距离差

	double dX = 0.0;
	double dY = 0.0;
	double dTempRange = 0;
	int nTotalLen = sizeof(double) * nColLen * nRowLen;
	memset(m_dTrackerMap, 0, nTotalLen);
	//pRstData = (double *)malloc(sizeof(double) * nColLen * nRowLen);
	//assert(NULL != pRstData && "malloc mem fail !");
	trackerData *pPoint = NULL;
	bool bFind = false;
	//
	for (int i = 0; i < nColLen; i++)
	{
		pPoint = pColData + i;
		for (int j = 0; j < nRowLen; j++)
		{
			dX = (pRowData + j)->dTragetX - pPoint->dTragetX;
			dY = (pRowData + j)->dTargetY - pPoint->dTargetY;
			dTempRange = dX * dX + dY * dY;
			if (m_nThreshold < dTempRange)
				//pRstData[i + nColLen * j] = MAX_INF;
				m_dTrackerMap[i + nColLen *j] = MAX_INF;
			else
			{
				m_dTrackerMap[i + nColLen * j] = dTempRange;
				//pRstData[i + nColLen * j] = /*sqrt*/(dTempRange);
				bFind = true;
			}
		}
	}
#ifndef NO_PRINT
	printf("\n******************xyl****************************\n");
	for (int i = 0; i < nRowLen; i++)
	{
		for (int j = 0; j < nColLen; j++)
		{
			printf("%0.3f  ", *(m_dTrackerMap + j + i * nColLen));
		}
		printf("\n");
	}
	printf("\n\n");
#endif

	if (false == bFind)
	{
		//free(pRstData);
		pRstData = NULL;
		memset(m_dTrackerMap, 0, nTotalLen);
		return pRstData;
	}
	pRstData = m_dTrackerMap;
	return pRstData;
}

void TrakerImp::calSpeedAndDirect(targetTrackerEx *pGjData, trackerData *pPointData, int nLoopId)//首次匹配
{
	if (NULL == pGjData || NULL == pPointData)
	{
		return;
	}
	double dX = pGjData->dYcTragetX - pPointData->dTragetX;
	double dY = pGjData->dYcTargetY - pPointData->dTargetY;
	double dSpeed = 0;
	double dRange = sqrt(dX * dX + dY * dY);
	dSpeed = dRange / 2 / (nLoopId - pGjData->nLoopId);
	double dAng = 0;
	dAng = acos( (pPointData->dTragetX - pGjData->dYcTragetX) / (dRange + 1e-3) ) / GNN_pi * 180;
	if (0 > (pPointData->dTargetY - pGjData->dYcTargetY))
	{
		dAng = 360 - dAng;
	}
	pPointData->dDirection = dAng;
	pPointData->dSpeed = dSpeed;
	pGjData->xEstm << pPointData->dTragetX, (pPointData->dTragetX - pGjData->dYcTragetX) / 2 / (nLoopId - pGjData->nLoopId), pPointData->dTargetY, (pPointData->dTargetY - pGjData->dYcTargetY) / 2 / (nLoopId - pGjData->nLoopId);

}

void TrakerImp::cpyTrackerData(trackerData *srcData, trackerData *disData)
{
	if (NULL == srcData || NULL == disData)
	{
		return;
	}
	memcpy(disData, srcData, sizeof(trackerData));
}

void TrakerImp::calYuceTracker(int nCurLoopId)
{
	int nMaxLoop = nCurLoopId - 1;
	int nMinLoop = nMaxLoop - m_nMax_gap_closing;
	int nHisLoopId = 0;
	if (0 > nMinLoop)
	{
		nMinLoop = 0;
	}
	int nBetwnLoop = 0;
	double dAngDir = 0;
	double dVDir = 0;
	double dDirFx = 0;

	TarcksMap *pMap = m_pMapTraker;
	m_nHisDataLen = 0;
	targetTrackerEx *pDealTracker = NULL;
	while (NULL != pMap)
	{
		nHisLoopId = pMap->pTaget->nLoopId;
		pDealTracker = (targetTrackerEx *)pMap->pTaget;
		if (nMinLoop <= nHisLoopId && nMaxLoop >= nHisLoopId)
		{
			//IMM混合预测
			Vector4d xPred1, xPred2,xPred3;
			Matrix4d Cov1, Cov2, Cov3;
			Cov1 = pDealTracker->F1 * pDealTracker->mErrorEstm * pDealTracker->F1.transpose() + pDealTracker->Gamma;
			Cov2 = pDealTracker->F2 * pDealTracker->mErrorEstm * pDealTracker->F2.transpose() + pDealTracker->Gamma;
			Cov3 = pDealTracker->F3 * pDealTracker->mErrorEstm * pDealTracker->F3.transpose() + pDealTracker->Gamma;

			xPred1 = pDealTracker->F1 * pDealTracker->xEstm;
			xPred2 = pDealTracker->F2 * pDealTracker->xEstm;
			xPred3 = pDealTracker->F3 * pDealTracker->xEstm;

			if ((nCurLoopId - nHisLoopId)>=2) //如果间隔大于一帧，则需要预测远一点
			{
				for (int i = 2; i <= (nCurLoopId - nHisLoopId);i++)
				{
					xPred1 = pDealTracker->F1 * xPred1;
					xPred2 = pDealTracker->F2 * xPred2;
					xPred3 = pDealTracker->F3 * xPred3;
					Cov1 = pDealTracker->F1 * Cov1 * pDealTracker->F1.transpose() + pDealTracker->Gamma;
					Cov2 = pDealTracker->F2 * Cov2 * pDealTracker->F2.transpose() + pDealTracker->Gamma;
					Cov3 = pDealTracker->F3 * Cov3 * pDealTracker->F3.transpose() + pDealTracker->Gamma;

				}
			}
			// Mixup
			//printf("%f", pDealTracker->Mp(1));
			//cout <<  pDealTracker->Mp(2) << endl;
			pDealTracker->mErrorPred = pDealTracker->Mp(0)*Cov1 + pDealTracker->Mp(1)*Cov2 + pDealTracker->Mp(2)*Cov3;
			pDealTracker->xPred = pDealTracker->Mp(0)*xPred1 + pDealTracker->Mp(1)*xPred2 + pDealTracker->Mp(2)*xPred3;
			

			//dAngDir = pDealTracker->dDirection / 180 * GNN_pi;
			//dVDir = pDealTracker->dSpeed;
			//nBetwnLoop = 2 * (nCurLoopId - nHisLoopId);
			pDealTracker->dYcTargetY = pDealTracker->xPred(3);// pDealTracker->dTargetY + dVDir * sin(dAngDir) * nBetwnLoop;
			pDealTracker->dYcTragetX = pDealTracker->xPred(1);// pDealTracker->dTragetX + dVDir * cos(dAngDir) * nBetwnLoop;
			m_nHisDataLen++;
		}
		pMap = pMap->pNext;
	}
}

void TrakerImp::calAlphaAndBeta(int nCurLoopId, targetTrackerEx *pOldGj, trackerData *pNewGj)
{
	double dXVelPrdct = 0;
	double dYVelPrdct = 0;
	double dXPosPrdct = 0;
	double dYPosPrdct = 0;
	double dXInnov = 0;
	double dYInnov = 0;
	double fAngDir = 0;
	double dVelAmpSmooth = 0;
	double dVelAngSmooth = 0;

	int nBetwnLoop = 2 * (nCurLoopId - pOldGj->nLoopId);
	
	//IMM 与 Kalman 专区
	Vector2d innov0, innov1, innov2;
	Vector3d cj,ui0,ui1,ui2;
	Vector4d x0, x1, x2, Pos0, Pos1, Pos2, xPred0, xPred1, xPred2;
	Matrix4d P0, P1, P2, I,M2; I << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
	Matrix4d pPred0, pPred1, pPred2;
	Matrix2d R,M1,S0,S1,S2; R << 10, 0, 0, 10;
	Matrix <double, 4, 2> K;
	cj = pOldGj->pij * pOldGj->Mp;
	ui0 = (1 / cj(0)) * pOldGj->pij.col(0).cross(pOldGj->Mp);
	ui1 = (1 / cj(1)) * pOldGj->pij.col(1).cross(pOldGj->Mp);
	ui2 = (1 / cj(2)) * pOldGj->pij.col(2).cross(pOldGj->Mp);

	x0 = pOldGj->PosState.col(0) * ui0(0) + pOldGj->PosState.col(1) * ui0(1) + pOldGj->PosState.col(2) * ui0(2);
	x1 = pOldGj->PosState.col(0) * ui1(0) + pOldGj->PosState.col(1) * ui1(1) + pOldGj->PosState.col(2) * ui1(2);
	x2 = pOldGj->PosState.col(0) * ui2(0) + pOldGj->PosState.col(1) * ui2(1) + pOldGj->PosState.col(2) * ui2(2);
	Pos0 = pOldGj->PosState.col(0); Pos1 = pOldGj->PosState.col(1); Pos2 = pOldGj->PosState.col(2);
	Vector4d Dif0, Dif1, Dif2;
	Dif0 = Pos0 - x0; Dif1 = Pos1 - x0; Dif2 = Pos2 - x0;
	P0 = (pOldGj->CovInit1 + (Dif0)*(Pos0.transpose() - x0.transpose())) * ui0(0) +
		 (pOldGj->CovInit2 + (Dif1)*(Pos1.transpose() - x0.transpose())) * ui0(1) +
		 (pOldGj->CovInit3 + (Dif2)*(Pos2.transpose() - x0.transpose())) * ui0(2);
	Dif0 = Pos0 - x1; Dif1 = Pos1 - x1; Dif2 = Pos2 - x1;
	P1 = (pOldGj->CovInit1 + (Dif0)*(Pos0.transpose() - x1.transpose())) * ui1(0) +
		(pOldGj->CovInit2 + (Dif1)*(Pos1.transpose() - x1.transpose())) * ui1(1) +
		(pOldGj->CovInit3 + (Dif2)*(Pos2.transpose() - x1.transpose())) * ui1(2);
	Dif0 = Pos0 - x2; Dif1 = Pos1 - x2; Dif2 = Pos2 - x2;
	P2 = (pOldGj->CovInit1 + (Dif0)*(Pos0.transpose() - x2.transpose())) * ui2(0) +
		(pOldGj->CovInit2 + (Dif1)*(Pos1.transpose() - x2.transpose())) * ui2(1) +
		(pOldGj->CovInit3 + (Dif2)*(Pos2.transpose() - x2.transpose())) * ui2(2);
	//Kalman过程
	//F1 过程
	for (int i = 1; i <= (nCurLoopId - pOldGj->nLoopId);i++)
	{
		if (i != 1)
		{
			xPred0 = pOldGj->F1 * xPred0;
			pPred0 = pOldGj->F1 * pPred0 + pOldGj->Gamma;
		}
		xPred0 = pOldGj->F1 * pOldGj->xEstm;
		pPred0 = pOldGj->F1 * pOldGj->mErrorEstm + pOldGj->Gamma;
	}
	M1 = pOldGj->H * pPred0 * pOldGj->H.transpose() + R;
	S0 = M1;
	M1.transposeInPlace(); 
	M1 = M1.inverse().eval();
	innov0(0) = pNewGj->dTragetX - xPred0(0);
	innov0(1) = pNewGj->dTargetY - xPred0(2);
	K = pPred0*pOldGj->H.transpose()*M1;
	M2 = I - K*pOldGj->H;
	pOldGj->PosState.col(0) = xPred0 + K * innov0;
	pOldGj->CovInit1 = (I - K * pOldGj->H) * pPred0 * M2.transpose() + K * R * K.transpose();
	
	//F2
	for (int i = 1; i <= (nCurLoopId - pOldGj->nLoopId); i++)
	{
		if (i != 1)
		{
			xPred1 = pOldGj->F2 * xPred1;
			pPred1 = pOldGj->F2 * pPred1 + pOldGj->Gamma;
		}
		xPred1 = pOldGj->F2 * pOldGj->xEstm;
		pPred1 = pOldGj->F2 * pOldGj->mErrorEstm + pOldGj->Gamma;
	}
	M1 = pOldGj->H * pPred1 * pOldGj->H.transpose() + R;
	S1 = M1;
	M1.transposeInPlace(); M1 = M1.inverse().eval();
	innov1(0) = pNewGj->dTragetX - xPred1(0);
	innov1(1) = pNewGj->dTargetY - xPred1(2);
	K = pPred1*pOldGj->H.transpose()*M1;
	M2 = I - K*pOldGj->H;
	pOldGj->PosState.col(1) = xPred1 + K * innov1;
	pOldGj->CovInit2 = (I - K * pOldGj->H) * pPred1 * M2.transpose() + K * R * K.transpose();
	
	//F3
	for (int i = 1; i <= (nCurLoopId - pOldGj->nLoopId); i++)
	{
		if (i != 1)
		{
			xPred2 = pOldGj->F3 * xPred2;
			pPred2 = pOldGj->F3 * pPred2 + pOldGj->Gamma;
		}
		xPred2 = pOldGj->F1 * pOldGj->xEstm;
		pPred2 = pOldGj->F1 * pOldGj->mErrorEstm + pOldGj->Gamma;
	}
	M1 = pOldGj->H * pPred2 * pOldGj->H.transpose() + R;
	S2 = M1;
	M1.transposeInPlace(); M1 = M1.inverse().eval();
	innov2(0) = pNewGj->dTragetX - xPred2(0);
	innov2(1) = pNewGj->dTargetY - xPred2(2);
	K = pPred2*pOldGj->H.transpose()*M1;
	M2 = I - K*pOldGj->H;
	pOldGj->PosState.col(2) = xPred2 + K * innov2;
	pOldGj->CovInit3 = (I - K * pOldGj->H) * pPred2 * M2.transpose() + K * R * K.transpose();

	//计算概率
	double Lfun1, Lfun2, Lfun3,LfunSum;
	double cTemp;
	cTemp = (innov0.transpose() * S0.inverse() * innov0)  ;
	Lfun1 = (1 / sqrt(fabs(2 * GNN_pi * S0.determinant()))) * exp((-1 / 2) * cTemp);
	cTemp = (innov1.transpose() * S1.inverse() * innov1);
	Lfun2 = (1 / sqrt(fabs(2 * GNN_pi * S1.determinant()))) * exp((-1 / 2) * cTemp);
	cTemp = (innov2.transpose() * S2.inverse() * innov2);
	Lfun3 = (1 / sqrt(fabs(2 * GNN_pi * S2.determinant()))) * exp((-1 / 2) * cTemp);
	LfunSum = Lfun1 + Lfun2 + Lfun3;
	Lfun1 = Lfun1 / LfunSum;
	Lfun2 = Lfun2 / LfunSum;
	Lfun3 = Lfun3 / LfunSum;
	Vector3d c(Lfun1, Lfun2, Lfun3);
	cTemp = c.transpose() * cj;
	pOldGj->Mp(0) = c(0)*cj(0) * (1 / cTemp);
	pOldGj->Mp(1) = c(1)*cj(1) * (1 / cTemp);
	pOldGj->Mp(2) = c(2)*cj(2) * (1 / cTemp);

	//状态混合
	Pos0 = pOldGj->PosState.col(0); Pos1 = pOldGj->PosState.col(1); Pos2 = pOldGj->PosState.col(2);
	pOldGj->xEstm = pOldGj->PosState.col(0)*pOldGj->Mp(0) + pOldGj->PosState.col(1)*pOldGj->Mp(1) + pOldGj->PosState.col(2)*pOldGj->Mp(2);
	pOldGj->mErrorEstm = (pOldGj->CovInit1 + (Pos0 - pOldGj->xEstm) * (Pos0.transpose() - pOldGj->xEstm.transpose())) * pOldGj->Mp(0) +
						(pOldGj->CovInit2 + (Pos1 - pOldGj->xEstm) * (Pos1.transpose() - pOldGj->xEstm.transpose())) * pOldGj->Mp(1) +
						(pOldGj->CovInit3 + (Pos2 - pOldGj->xEstm) * (Pos2.transpose() - pOldGj->xEstm.transpose())) * pOldGj->Mp(2);
	
	//更新


	/*
	fAngDir = pOldGj->dDirection * GNN_pi / 180;

	dXVelPrdct = pOldGj->dSpeed * cos(fAngDir);
	dYVelPrdct = pOldGj->dSpeed * sin(fAngDir);
	dXPosPrdct = pOldGj->dTragetX + nBetwnLoop * dXVelPrdct;
	dYPosPrdct = pOldGj->dTargetY + nBetwnLoop * dYVelPrdct;

	dXInnov = pNewGj->dTragetX - dXPosPrdct;
	dYInnov = pNewGj->dTargetY - dYPosPrdct;

	dXPosPrdct += dXInnov * m_dAlpha;
	dYPosPrdct += dYInnov * m_dAlpha;

	pNewGj->dTragetX = dXPosPrdct;
	pNewGj->dTargetY = dYPosPrdct;

	pNewGj->fRange = (float)sqrt(dXPosPrdct * dXPosPrdct + dYPosPrdct * dYPosPrdct);

	dXVelPrdct = dXVelPrdct + dXInnov * m_dBeta / nBetwnLoop;
	dYVelPrdct = dYVelPrdct + dYInnov * m_dBeta / nBetwnLoop;

	dVelAmpSmooth = sqrt(dXVelPrdct * dXVelPrdct + dYVelPrdct * dYVelPrdct);
	pNewGj->dSpeed = dVelAmpSmooth;
	dVelAngSmooth = acos(dXVelPrdct / (dVelAmpSmooth + 1e-5)) / GNN_pi * 180;
	*/
	pNewGj = pOldGj;
	pNewGj->dTragetX = pNewGj->xEstm(0);
	pNewGj->dTargetY = pNewGj->xEstm(2);
	pNewGj->dSpeed = (float)sqrt(pNewGj->xEstm(1)*pNewGj->xEstm(1) + pNewGj->xEstm(3)*pNewGj->xEstm(3));
	dVelAngSmooth = acos(pNewGj->xEstm(1) / (pNewGj->dSpeed + 1e-5)) / GNN_pi * 180;
	if (0 > dYVelPrdct)
	{
		dVelAngSmooth = 360 - dVelAngSmooth;
	}
	pNewGj->dDirection = dVelAngSmooth;
}

