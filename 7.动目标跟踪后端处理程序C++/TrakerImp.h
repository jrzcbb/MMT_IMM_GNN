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

	//������������
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

	//�㼣���켣�����ž���
	targetTracker *matchingPointToTracker(int nMatchLoopId, trackerData *pData, int &nMatchId, double &dMinRange);

	//�켣���㼣�����ž���, nLoopId : ƥ����ʷȦ��ID;pTarget : ƥ��Ŀ��㼣��dMinRange : �㼣ƥ�����С���룬
	trackerData *matchingTrackerToPoint(targetTracker *pTracker, double dMinRange, trackerData *pTrackData, int nDataLen);

	//����Ŀ���б�
	void updateMapData(int nLoopId, targetTracker *pTracker, trackerData  *pUpdateData);

	//��������ɾ��
	//int deleteTracke(int nMinId);

	void deleteTracke(int nMinId, trackerData *pLostData, int nMaxLen, int &nRealLostLen);

	virtual void getLoseTrackers(int *pTargetIds, int nLoseTrackerLen);

	//ͨ��loopId�Ż�ȡ�켣��Ϣ
	trackerData *calTrackerPointByLoop(int nLoopId, int &nLen);

	//ͨ���켣��Ϣ�Լ��㼣��Ϣ�����������㷨���ݣ�
	double *createMunkresData(trackerData *pRowData, int nRowLen, trackerData *pColData, int nColLen);

	void calSpeedAndDirect(targetTrackerEx *pGjData, trackerData *pPointData, int nLoopId);

	void cpyTrackerData(trackerData *srcData, trackerData *disData);

	void calYuceTracker(int nCurLoopId);


	void calAlphaAndBeta(int nCurLoopId, targetTrackerEx *pOldGj, trackerData *pNewGj);


private:
	//�����б�
	TarcksMap *m_pMapTraker;
	//���һ���������ݣ�������������
	TarcksMap *m_pTailTraker;
	//��һȦ�����ۻ�ID��
	int m_nFirsLoopId;
	//����
	int m_nThreshold;
	//������ID�ţ� Ŀǰ�ۼӷ�ʽά����Խ�����
	int m_nMaxTargetId;
	//����������
	int m_nMax_linking_distance;
	//����Ȧ��
	int m_nMax_gap_closing;
	//��ʧ������
	std::vector<int> m_vcLoseIds;

	trackerData m_trackersPoint[1024];

	//�켣λ�õ�Ԥ��
	//targetTracker m_trackerCal[1024];

	double m_dTrackerRows[1024];

	double m_dTrackerMap[1024 * 1024];

	int m_nHisDataLen;



	//Alpha Beta�˲�
	double m_dAlpha;

	double m_dBeta;

};

