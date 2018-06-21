#include <Eigen/Dense>

#ifndef GNNTRACKERCOM_H
#define GNNTRACKERCOM_H


#define TARGET_NULL -111

#define MAX_INF 9999999

#define GNN_pi 3.14159265358979323846
using namespace Eigen;
	//��������

	struct trackerData
	{
		
		unsigned int nOutId;			//ID�ţ��������λ��	
		double dTragetX;    // range * cos(Ang)
		double dTargetY;	//range * sin(Ang)
		double dSpeed;		//�ٶ�
		double dDirection;  //����
		float fRange;
		float fAng;
		//Kalman �� IMM 
		//EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		Vector4d xPred;     //Vector4f ����4ά������
		Vector4d xEstm;
		Matrix4d mErrorPred;
		Matrix4d mErrorEstm;
		Matrix4d Gamma;
		Matrix3d pij;
		Vector3d Mp;
		Matrix <double, 4, 3> PosState;
		Matrix <double, 2, 4> H;
		Matrix4d CovInit1;
		Matrix4d CovInit2;
		Matrix4d CovInit3;
		Matrix4d F1;
		Matrix4d F2;
		Matrix4d F3;
		//struct IMMPar;
		trackerData()
		{
			nOutId = TARGET_NULL;
			dTragetX = -1;
			dTargetY = -1;
			dSpeed = 0;
			dDirection = 0;
			H << 1, 0, 0, 0,
				0, 0, 1, 0;
			xPred << dTragetX, 0, dTargetY, 0;
			xEstm << dTragetX, 0, dTargetY, 0;
			mErrorPred << 100, 50, 0, 0,
				50, 50.2500000000000, 0, 0,
				0, 0, 100, 50,
				0, 0, 50, 50.2500000000000;
			mErrorEstm << 100, 50, 0, 0,
				50, 50.2500000000000, 0, 0,
				0, 0, 100, 50,
				0, 0, 50, 50.2500000000000;
			Gamma << 2, 2, 0, 0,
				2, 2, 0, 0,
				0, 0, 2, 2,
				0, 0, 2, 2;
			PosState << xPred, xPred, xPred;
			pij << 0.9, 0.05, 0.05,
				0.1, 0.8, 0.1,
				0.05, 0.15, 0.8;
			Mp << 0.3, 0.3, 0.4;
			F1 << 1, 2, 0, 0,
				0, 1, 0, 0,
				0, 0, 1, 2,
				0, 0, 0, 1;
			F2 << 1, 1.99959386818404, 0, -0.0349030407992163,
				0, 0.999390827019096, 0, 0.0348994967025010,
				0, 0.0349030407992163, 1, 1.99959386818404,
				0, -0.0348994967025010, 0, 0.999390827019096;
			F3 << 1, 1.99634659474160, 0, 0.104624091709670,
				0, 0.994521895368273, 0, -0.104528463267653,
				0, -0.104624091709670, 1, 1.99634659474160,
				0, 0.104528463267653, 0, 0.994521895368273;
		}
	};
	/*/IMM����
	struct IMMPar
	{
		Matrix3f pij;
		Vector3f Mp;
		Matrix <float, 4, 3> PosState;
		Matrix4f CovInit1;
		Matrix4f CovInit2;
		Matrix4f CovInit3;
		Matrix4f F1;
		Matrix4f F2;
		Matrix4f F3;
		IMMPar()
		{

		}
	};
	*/
	//����·��
	struct targetTracker : public trackerData
	{
		int nLoopId;		//֡���
		int nPtLen;			//����֡��	
		targetTracker *pNextTrack;  //�Ƿ���Ҫά������·�������Ǽ����·�����ڹ�������ĿǰԤ��
	};

	//Ŀ������б�
	struct TarcksMap
	{
		int nId;				//�켣ID��
		int nObjType;			//Ŀ������
		double dThreatDegree;   //��в�̶�
		double begTime;			//��ʼ����ʱ��
		targetTracker *pTaget;  //����Ŀ��
		TarcksMap *pNext;		//�¸�Ŀ��
	};

	//��ô���ԣ�

//�򵥵��״�����㷨,���㷨�������⣺����Ŀ��ID��Ϊ�ۼ�ģʽ������Խ�磻loop����Խ��
class GNN_Traker
{
public:
	//��������,nMax_linking_distance : ���������룻nMax_gap_closing �� ƥ��Ȧ��
	virtual void setTrackerInfo(int nMax_linking_distance, int nMax_gap_closing) = 0;

	//�����㷨�ӿڣ�nTotalLoop ����ǰȦ����pTrackData ���״�Ŀ�����ݣ� nDataLen �� �״����ݳ��ȣ�nLoseTrackerLen �� ��ʧ����
	virtual void GNNTracker_Sequt(int nTotalLoop, trackerData *pTrackData, int nDataLen, int &nLoseTrackerLen, trackerData *pLostData, int nMaxLostLen) = 0;

	//��ʧ���ݵ�����Ϣ,����ʽ������GNNTracker_Sequt��Ḳ����һ�εĶ�ʧ���ݣ�pTargetIds[out] : ��ʧID�ţ�nLoseTrackerLen ����ʧ���ݳ���
	virtual void getLoseTrackers(int *pTargetIds, int nLoseTrackerLen) = 0;

	//�ͷ��ڴ�
	virtual void releaseTraker() = 0;

	virtual void clearTrakers() = 0;

	//Ԥ���ˡ���
	//virtual void clearLoop() = 0;

	//ͨ��Ȧ��ID��ȡMap
	//virtual TarcksMap *getMapByLoop(int nLoop) = 0;
};

extern "C" GNN_Traker *createGNNTracker();

extern "C" void assignmentoptimal(double *assignment, double *cost, double *distMatrix, int nOfRows, int nOfColumns);

#endif