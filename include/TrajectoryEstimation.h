#ifndef _TRAJECTORYESTIMATION
#define _TRAJECTORYESTIMATION

#include <stdio.h>
#include "MathLib.h"
#include <boost/circular_buffer.hpp>
#include "eigen3/Eigen/Dense"

// 1 cm nosiy
#define MAX_POSITION_NOISE 0.01

enum PREDICTION_METHOD { PREDICTION_2NDORDER, PREDICTION_BALLISTIC, PREDICTION_SVRRBF, PREDICTION_GMR };

using namespace MathLib;

struct sObServation {
	double time;
	MathLib::Vector3 pos;
	MathLib::Vector3 vel;
	MathLib::Matrix3 R;
	MathLib::Vector3 omega;
	MathLib::Vector3 angularMomentum;
};


class TrajectoryEstimation
{
protected:
	SpatialInertia mSpatialInertia;

	double mDt;
	int    mNbReservoir;
	int    mStartPredictionFrame;

	double mAirDrag;
	MathLib::Vector3 mGravity;

	double mCurrentTime;
	MathLib::Vector3 mCurrentPos;
	MathLib::Vector3 mCurrentVel;
	MathLib::Vector3 mCurrentAngularVel;
	MathLib::Matrix3 mCurrentR;
	MathLib::Vector3 mCurrentAngularMomentum;
	MathLib::Vector3 mCurrentAccel;

	MathLib::Matrix mRinv;
	MathLib::Vector mAngularVelocity;
	MathLib::Matrix3 mW, mRdot;

	sObServation mObservationOld;
	sObServation mObservation;

	deque<MathLib::Matrix> mlD;
	deque<MathLib::Vector> mlC;

	int *mReservoirIndexK;

public:
	boost::circular_buffer<sObServation> mObservationReservoirX;
	bool mReadyToPredict;


	TrajectoryEstimation(double gravity[3],double integrationStep=(1./240.), int nbReservoir=21, int StartPredictionFrame=9);
	~TrajectoryEstimation();


	void setAirdrag(double lAirdrag);
	// 
	void setMomentOfInertia(MathLib::Matrix3 Iobj);
	void setMomentOfInertia(Eigen::Matrix3d Iobj);
	
	void resetTrajectoryEstimation(void);

	// observation
	void setNewObservationset(double t, MathLib::Vector3 position, MathLib::Matrix3 orientation);
	void setNewObservationset(double t, Eigen::Vector3d position);
	
	//
	void InitFitting(void);
	int estimateCurrentState(bool DelNoisyData=false);

	// prediction
	void PredictNextPos(double dt, MathLib::Vector3& pos);
	int PredictNextPos(double duration, MathLib::Matrix& trajectory);
	int PredictNextPos(double duration, double dt, MathLib::Matrix& trajectory);

	int PredictNextPosVel(double duration, MathLib::Matrix& postrj, MathLib::Matrix& veltrj);
	int PredictNextPosVel(double duration, double dt, MathLib::Matrix& postrj, MathLib::Matrix& veltrj);
	int PredictNextPosVel(double duration, double dt, Eigen::MatrixXd& postrj, bool constraint=false, int Direction=0, double value=0);

	void PredictEndPos(double duration, MathLib::Vector3& pos);
	void PredictEndPos(double duration, double dt, MathLib::Vector3& pos);
	void PredictEndPosAndVel(double duration, double dt, MathLib::Vector3& pos, MathLib::Vector3& vel);
	void PredictEndPosAndVel_Predicted(double duration, double dt, MathLib::Vector3& pos, MathLib::Vector3& vel);

	void DistanceCatchingObject(double duration, double dt, MathLib::Vector3& pos, MathLib::Vector3& vel);

	void PredictNextOrient(double dt, MathLib::Matrix3& orient);
	int PredictNextOrient(double duration, MathLib::Matrix& dirx, MathLib::Matrix& diry);
	int PredictNextOrient(double duration, double dt, MathLib::Matrix& dirx, MathLib::Matrix& diry);

	void PredictEndOrient(double duration, MathLib::Vector3& dirx, MathLib::Vector3& diry);
	void PredictEndOrient(double duration, double dt, MathLib::Vector3& dirx, MathLib::Vector3& diry);

	double getCurrentTime();

};

#endif // _TRAJECTORYESTIMATION
