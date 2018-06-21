#include <Eigen/Dense>

#ifndef GNNTRACKERCOM_H
#define GNNTRACKERCOM_H


#define TARGET_NULL -111

#define MAX_INF 9999999

#define GNN_pi 3.14159265358979323846
using namespace Eigen;
	//跟踪数据

	struct trackerData
	{
		
		unsigned int nOutId;			//ID号，是链表的位置	
		double dTragetX;    // range * cos(Ang)
		double dTargetY;	//range * sin(Ang)
		double dSpeed;		//速度
		double dDirection;  //方向
		float fRange;
		float fAng;
		//Kalman 与 IMM 
		//EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		Vector4d xPred;     //Vector4f 定义4维列向量
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
	/*/IMM参数
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
	//跟踪路径
	struct targetTracker : public trackerData
	{
		int nLoopId;		//帧序号
		int nPtLen;			//跟上帧数	
		targetTracker *pNextTrack;  //是否需要维护跟踪路径，还是计算出路径后不在管理？？？目前预留
	};

	//目标跟踪列表
	struct TarcksMap
	{
		int nId;				//轨迹ID号
		int nObjType;			//目标类型
		double dThreatDegree;   //威胁程度
		double begTime;			//初始跟踪时间
		targetTracker *pTaget;  //跟踪目标
		TarcksMap *pNext;		//下个目标
	};

	//怎么测试？

//简单的雷达跟踪算法,该算法存在问题：跟踪目标ID号为累加模式，可能越界；loop可能越界
class GNN_Traker
{
public:
	//参数设置,nMax_linking_distance : 最大关联距离；nMax_gap_closing ： 匹配圈数
	virtual void setTrackerInfo(int nMax_linking_distance, int nMax_gap_closing) = 0;

	//跟踪算法接口，nTotalLoop ：当前圈数；pTrackData ：雷达目标数据； nDataLen ： 雷达数据长度；nLoseTrackerLen ： 丢失数据
	virtual void GNNTracker_Sequt(int nTotalLoop, trackerData *pTrackData, int nDataLen, int &nLoseTrackerLen, trackerData *pLostData, int nMaxLostLen) = 0;

	//丢失数据点数信息,覆盖式：调用GNNTracker_Sequt后会覆盖上一次的丢失数据；pTargetIds[out] : 丢失ID号；nLoseTrackerLen ：丢失数据长度
	virtual void getLoseTrackers(int *pTargetIds, int nLoseTrackerLen) = 0;

	//释放内存
	virtual void releaseTraker() = 0;

	virtual void clearTrakers() = 0;

	//预留了。。
	//virtual void clearLoop() = 0;

	//通过圈数ID获取Map
	//virtual TarcksMap *getMapByLoop(int nLoop) = 0;
};

extern "C" GNN_Traker *createGNNTracker();

extern "C" void assignmentoptimal(double *assignment, double *cost, double *distMatrix, int nOfRows, int nOfColumns);

#endif