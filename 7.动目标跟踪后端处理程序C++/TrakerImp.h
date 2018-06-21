#pragma once

#include "GNNTrakerCom.h"

#include <vector>

#define TARGET_OCCUPY -456


struct targetTrackerEx : public targetTracker
{
	double dYcTragetX; 
	double dYcTargetY;
};
class TrakerImp : public GNN_Traker
{
public:
	TrakerImp(void);
	~TrakerImp(void);



	virtual void clearTrakers() override;

public:

	//新增航迹数据
	void insertNewTraker(int nLoopId, trackerData *pData, int nLengTh);


protected:

	virtual void setTrackerInfo(int nMax_linking_distance, int nMax_gap_closing);

/*	virtual void GNNTracker_Sequt(int nTotalLoop, trackerData *pTrackData, int nDataLen, int &nLoseTrackerLen);*/

	virtual void GNNTracker_Sequt(int nTotalLoop, trackerData *pTrackData, int nDataLen, int &nLoseTrackerLen, trackerData *pLostData, int nMaxLostLen);


	virtual void releaseTraker();

	void insertToMapTrakerToHeador(targetTracker *pData);

	void insertToMapTrackerToTail(targetTracker *pData);

	targetTrackerEx *createTracker(int &nLoopId, trackerData *pData);

	targetTrackerEx *getTrackerById(int nLoop, int nId);

	//点迹到轨迹的最优距离
	targetTracker *matchingPointToTracker(int nMatchLoopId, trackerData *pData, int &nMatchId, double &dMinRange);

	//轨迹到点迹的最优距离, nLoopId : 匹配历史圈数ID;pTarget : 匹配目标点迹，dMinRange : 点迹匹配的最小距离，
	trackerData *matchingTrackerToPoint(targetTracker *pTracker, double dMinRange, trackerData *pTrackData, int nDataLen);

	//更新目标列表
	void updateMapData(int nLoopId, targetTracker *pTracker, trackerData  *pUpdateData);

	//多余数据删除
	//int deleteTracke(int nMinId);

	void deleteTracke(int nMinId, trackerData *pLostData, int nMaxLen, int &nRealLostLen);

	virtual void getLoseTrackers(int *pTargetIds, int nLoseTrackerLen);

	//通过loopId号获取轨迹信息
	trackerData *calTrackerPointByLoop(int nLoopId, int &nLen);

	//通过轨迹信息以及点迹信息生成匈牙利算法数据；
	double *createMunkresData(trackerData *pRowData, int nRowLen, trackerData *pColData, int nColLen);

	void calSpeedAndDirect(targetTrackerEx *pGjData, trackerData *pPointData, int nLoopId);

	void cpyTrackerData(trackerData *srcData, trackerData *disData);

	void calYuceTracker(int nCurLoopId);


	void calAlphaAndBeta(int nCurLoopId, targetTrackerEx *pOldGj, trackerData *pNewGj);


private:
	//跟踪列表
	TarcksMap *m_pMapTraker;
	//最后一个跟踪数据，用于链表运算
	TarcksMap *m_pTailTraker;
	//第一圈数据累积ID号
	int m_nFirsLoopId;
	//门限
	int m_nThreshold;
	//最大对象ID号： 目前累加方式维护有越界风险
	int m_nMaxTargetId;
	//最大关联距离
	int m_nMax_linking_distance;
	//计算圈数
	int m_nMax_gap_closing;
	//丢失的数据
	std::vector<int> m_vcLoseIds;

	trackerData m_trackersPoint[1024];

	//轨迹位置的预测
	//targetTracker m_trackerCal[1024];

	double m_dTrackerRows[1024];

	double m_dTrackerMap[1024 * 1024];

	int m_nHisDataLen;



	//Alpha Beta滤波
	double m_dAlpha;

	double m_dBeta;

};

