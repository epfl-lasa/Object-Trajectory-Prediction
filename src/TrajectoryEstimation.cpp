#include "TrajectoryEstimation.h"

TrajectoryEstimation::TrajectoryEstimation(double gravity[3], double integrationStep, int nbReservoir, int StartPredictionFrame)
{
	mGravity(0) =  gravity[0];
	mGravity(1) =  gravity[1];
	mGravity(2) = gravity[2];
	mAirDrag=0.00;

	cout<<"mAirDrag "<<mAirDrag <<endl;
	mDt = integrationStep;
	mNbReservoir = nbReservoir;
	mStartPredictionFrame = StartPredictionFrame;

	mRinv.Resize(3,3);
	mAngularVelocity.Resize(3);

	mObservationReservoirX.set_capacity(mNbReservoir);
	mObservationReservoirX.clear();

	mReservoirIndexK = (int *)malloc(sizeof(int)*mNbReservoir);
	InitFitting();
	resetTrajectoryEstimation();
}


TrajectoryEstimation::~TrajectoryEstimation()
{

}
void TrajectoryEstimation::setAirdrag(double lAirdrag)
{
	mAirDrag = lAirdrag;
	InitFitting();
	resetTrajectoryEstimation();
}
void TrajectoryEstimation::setMomentOfInertia(Eigen::Matrix3d Iobj)
{
	MathLib::Matrix3 handle;
	for (int i=0; i<3;i++)
		for (int j=0; j<3;j++)
		{
			{
				handle(i,j)=Iobj(i,j);
			}
		}
	mSpatialInertia.mInertiaMoment.Set(handle);
}

void TrajectoryEstimation::setMomentOfInertia(MathLib::Matrix3 Iobj)
{
	mSpatialInertia.mInertiaMoment.Set(Iobj);
}


void TrajectoryEstimation::resetTrajectoryEstimation(void)
{
	mObservation.time = -1.0;
	mObservationOld.time = -1.0;

	mObservationReservoirX.clear();

	mReadyToPredict = false;
}

void TrajectoryEstimation::setNewObservationset(double t, Eigen::Vector3d position)
{
	MathLib::Vector3 handle;handle(0)=position(0);handle(1)=position(1);handle(2)=position(2);
	MathLib::Matrix3 matrix;matrix.Zero();

	setNewObservationset(t, handle,matrix);
}

// observation
void TrajectoryEstimation::setNewObservationset(double t, MathLib::Vector3 position, MathLib::Matrix3 orientation)
{
	MathLib::Matrix3 dR;
	int acnt;
	unsigned int lReservoirSize;

	mObservation.time = t;

	/*	if ((position(0)<0.0)&&(position(2)<mObservationOld.pos(2)))
	{
	//	cout<<"Position "<<position(0)<<" "<<position(1)<<" "<<position(2)<<endl;
		double a=mObservationOld.pos(2)-position(2);
		position(2)=position(2)+a;
	//	cout<<"Position "<<position(0)<<" "<<position(1)<<" "<<position(2)<<endl;
	}
	else
	{
	//	cout<<"Everything is awesome "<<endl;
	}*/
	mObservation.pos.Set( position );
	mObservation.R.Set( orientation );

	// calculate velocity
	mObservation.vel = ( mObservation.pos - mObservationOld.pos)/(mObservation.time - mObservationOld.time );
	//mObservation.vel.Print("mObservation.vel");

	// calculate angular velocity
	dR = (mObservation.R-mObservationOld.R)/(mObservation.time - mObservationOld.time );
	mW = dR*mObservation.R.Transpose();
	mObservation.omega(0) = (-mW(1,2)+ mW(2,1))/2;
	mObservation.omega(1) = ( mW(0,2)- mW(2,0))/2;
	mObservation.omega(2) = (-mW(0,1)+ mW(1,0))/2;

	// calculate angular momentum
	mObservation.angularMomentum = mObservation.R * mSpatialInertia.mInertiaMoment * mObservation.R.Transpose() * mObservation.omega;

	// replace old value
	if( mObservationOld.time < 0 )
	{
		mObservationReservoirX.clear();
		mObservationReservoirX.push_back(mObservation);
		mObservationOld = mObservation;
		return;
	}
	else if( (mObservation.time - mObservationReservoirX[0].time) >= 1.0 )
	{
		mObservationReservoirX.clear();
		mObservationReservoirX.push_back(mObservation);
		mObservationOld = mObservation;
		return;
	}
	else{
		if( (mObservation.time - mObservationOld.time) <= mDt ){
			return;
		}

		mObservationOld = mObservation;
	}

	mObservationReservoirX.push_back(mObservation);
	lReservoirSize = mObservationReservoirX.size();

	if( lReservoirSize >= mStartPredictionFrame ){
		mReadyToPredict = TRUE;
	}

	// calculate angular momentum
	if( mReadyToPredict ){
		mCurrentTime = mObservationReservoirX[lReservoirSize-1].time;

		/*
		if( estimateCurrentState(true) == -1)
		{
			estimateCurrentState(false);
			lReservoirSize -= 1;
		}
		 */
		estimateCurrentState(false);

		mCurrentAngularMomentum.Zero();
		acnt = 0;
		for(int i=0; i<lReservoirSize; i++ ){
			if( (mObservationReservoirX[i].R*mSpatialInertia.mInertiaMoment*mObservationReservoirX[i].R.Transpose()*mObservationReservoirX[i].omega).Norm() < 50.0*6.0){ //15.0
				mCurrentAngularMomentum += mObservationReservoirX[i].R*mSpatialInertia.mInertiaMoment*mObservationReservoirX[i].R.Transpose()*mObservationReservoirX[i].omega;
				acnt++;
			}
		}

		mCurrentAngularMomentum /= (double)acnt*(1.2);
		mCurrentR.Set(mObservation.R);
	}
}

void TrajectoryEstimation::InitFitting(void)
{
	MathLib::Matrix lA(6,6), lAi(6,6), lB(6,6);
	MathLib::Vector lvb(6);

	lA.Identity();
	lA(0,3+0) = mDt;
	lA(1,3+1) = mDt;
	lA(2,3+2) = mDt;
	lA(3,3+0) = (1.0 - 4.0*mAirDrag*mDt);
	lA(4,3+1) = (1.0 - 4.0*mAirDrag*mDt);
	lA(5,3+2) = (1.0 - 0.0*mAirDrag*mDt);

	lB.Identity();
	lvb.Zero();
	lvb(2) = 0.5*mGravity(2)*mDt*mDt;
	lvb(5) = mGravity(2)*mDt;

	mlD.clear();
	mlC.clear();
	lAi.Identity();
	for(int i=0; i<(1.0/mDt); i++)
	{
		lAi *= lA;
		mlD.push_back(lAi);

		if(i>0){
			lB += mlD.at(i-1);
		}
		mlC.push_back(lB*lvb);
	}
}

//void TrajectoryEstimation::estimateCurrentState(int kindex[], int nbObservation)
//{
//
//}


// if return is -1 : recall this one
int TrajectoryEstimation::estimateCurrentState(bool DelNoisyData)
{
	MathLib::Matrix lA(6,6), lB(6,6);
	MathLib::Vector lvb(6);

	MathLib::Vector Xx(2), Xy(2), Xz(2);
	MathLib::Vector X0(6), Xn(6);
	MathLib::Vector3 Xpos;

	double maxerror=0.0;
	double error;
	int errorIndex = 0;

	int nbObservation = mObservationReservoirX.size();

	MathLib::Matrix lCx(nbObservation-1,2);
	MathLib::Matrix lCy(nbObservation-1,2);
	MathLib::Matrix lCz(nbObservation-1,2);
	MathLib::Vector lYx(nbObservation-1), lYy(nbObservation-1), lYz(nbObservation-1);

	int k;
	for(int i=1; i<nbObservation; i++)
	{
		k = (int)( (mObservationReservoirX[i].time - mObservationReservoirX[0].time)/mDt );
		if(k<(1.0/mDt) && (k>=1))
		{
			lYx(i-1) = mObservationReservoirX[i].pos(0) - (mlC.at(k-1))(0);
			lYy(i-1) = mObservationReservoirX[i].pos(1) - (mlC.at(k-1))(1);
			lYz(i-1) = mObservationReservoirX[i].pos(2) - (mlC.at(k-1))(2);

			lCx(i-1,0) = (mlD.at(k-1))(0,0);
			lCx(i-1,1) = (mlD.at(k-1))(0,3);
			lCy(i-1,0) = (mlD.at(k-1))(1,1);
			lCy(i-1,1) = (mlD.at(k-1))(1,4);
			lCz(i-1,0) = (mlD.at(k-1))(2,2);
			lCz(i-1,1) = (mlD.at(k-1))(2,5);

			mReservoirIndexK[i] = k;
		}
		else if(k==0)
		{
			k=1;
			lYx(i-1) = mObservationReservoirX[i].pos(0) - (mlC.at(k-1))(0);
			lYy(i-1) = mObservationReservoirX[i].pos(1) - (mlC.at(k-1))(1);
			lYz(i-1) = mObservationReservoirX[i].pos(2) - (mlC.at(k-1))(2);

			lCx(i-1,0) = (mlD.at(k-1))(0,0);
			lCx(i-1,1) = (mlD.at(k-1))(0,3);
			lCy(i-1,0) = (mlD.at(k-1))(1,1);
			lCy(i-1,1) = (mlD.at(k-1))(1,4);
			lCz(i-1,0) = (mlD.at(k-1))(2,2);
			lCz(i-1,1) = (mlD.at(k-1))(2,5);
			cout << "k error 1" << endl;
		}
		else
		{
			cout << "k error, i: " << i  << "& k : " << k << endl;
			for(int j=0; j<nbObservation; j++)
			{
				cout << mObservationReservoirX[j].time << endl;
			}

		}
	}

	Xx = (lCx.Transpose()*lCx).Inverse()*lCx.Transpose()*lYx;
	Xy = (lCy.Transpose()*lCy).Inverse()*lCy.Transpose()*lYy;
	Xz = (lCz.Transpose()*lCz).Inverse()*lCz.Transpose()*lYz;

	X0(0) = Xx(0);
	X0(1) = Xy(0);
	X0(2) = Xz(0);
	X0(3) = Xx(1);

	X0(4) = Xy(1);
	X0(5) = Xz(1);

	if( k>0 ){

		Xn = mlD.at(k-1)*X0 + mlC.at(k-1);
		/*		cout<<"k "<<k<<endl;
		Xn.Print("Xn");
		X0.Print("X0");*/
		mCurrentPos.Set( Xn(0), Xn(1), Xn(2));
		mCurrentVel.Set( Xn(3), Xn(4), Xn(5));
		/*		mCurrentVel.Set( mObservation.vel(0), mObservation.vel(1), mObservation.vel(2));
		cout<<"Xn "<<Xn(3)<<" "<<Xn(4)<<" "<<Xn(5)<<endl;
		cout<<"mCurrentVel"<<mCurrentVel(0)<<" "<<mCurrentVel(1)<<" "<<mCurrentVel(2)<<endl;*/

		if( DelNoisyData )
		{
			// check nosy data
			for(int j=1; j<nbObservation; j++)
			{
				k = mReservoirIndexK[j];
				Xn = mlD.at(k-1)*X0 + mlC.at(k-1);
				Xpos.Set(Xn(0), Xn(1), Xn(2) );
				error = (mObservationReservoirX[j].pos - Xpos).Norm();
				if( error > maxerror)
				{
					errorIndex = j;
					maxerror = error;
				}
			}
			if( maxerror > MAX_POSITION_NOISE )
			{
				cout << "noisy data found : " << maxerror << endl;
				if( errorIndex == (nbObservation-1) )
				{
					mObservationReservoirX.erase_end(1);
					cout << errorIndex << " is removed!. Success" << endl;
				}
			}
		}
	}
	else{
		cout << "reservoir is zero " << endl;
	}

	if( maxerror > MAX_POSITION_NOISE ) return -1;
	else return 0;
}


void TrajectoryEstimation::PredictNextPos(double dt, MathLib::Vector3& pos)
{
	// position estimation
	mCurrentAccel = mGravity - mCurrentVel*mAirDrag;
	mCurrentPos   += mCurrentVel*dt + mCurrentAccel*dt*dt*0.5;
	mCurrentVel   += mCurrentAccel*dt;
	mCurrentTime  += dt;

	pos.Set(mCurrentPos);
}


int TrajectoryEstimation::PredictNextPos(double duration, MathLib::Matrix& trajectory)
{
	return PredictNextPos(duration, mDt, trajectory);
}

int TrajectoryEstimation::PredictNextPos(double duration, double dt, MathLib::Matrix& trajectory)
{
	unsigned int frame=0;
	MathLib::Vector3 lVel;
	MathLib::Vector3 lAccel;
	MathLib::Vector3 lPos;

	lPos.Set(mCurrentPos);
	lVel.Set(mCurrentVel);	
	while( dt*(frame+1) <= duration ){
		lAccel = mGravity - lVel*mAirDrag;
		lPos   += lVel*dt + lAccel*dt*dt*0.5;
		lVel   += lAccel*dt;

		trajectory.SetRow(lPos, frame++);
	}

	return frame;
}

int TrajectoryEstimation::PredictNextPosVel(double duration, MathLib::Matrix& postrj, MathLib::Matrix& veltrj)
{
	return PredictNextPosVel(duration, mDt, postrj, veltrj)	;
}

int TrajectoryEstimation::PredictNextPosVel(double duration, double dt, MathLib::Matrix& postrj, MathLib::Matrix& veltrj)
{
	unsigned int frame=0;
	MathLib::Vector3 lVel;
	MathLib::Vector3 lAccel;
	MathLib::Vector3 lPos;

	lPos.Set(mCurrentPos);
	lVel.Set(mCurrentVel);
	while( dt*(frame+1) <= duration ){
		lAccel = mGravity - lVel*mAirDrag;
		lPos   += lVel*dt + lAccel*dt*dt*0.5;
		lVel   += lAccel*dt;

		postrj.SetRow(lPos, frame);
		veltrj.SetRow(lVel, frame);
		frame++;
	}

	return frame;
}

int TrajectoryEstimation::PredictNextPosVel(double duration, double dt, Eigen::MatrixXd& postrj, bool constraint, int Direction, double value)
{
	unsigned int frame=0;
	MathLib::Vector3 lVel;
	MathLib::Vector3 lAccel;
	MathLib::Vector3 lPos;

	lPos.Set(mCurrentPos);
	lVel.Set(mCurrentVel);
	postrj(frame,0)=lPos(0);
	postrj(frame,1)=lPos(1);
	postrj(frame,2)=lPos(2);
	frame++;
	postrj.setZero();

	if (lVel.Norm2()<0.0000001)
	{
		return -1;
	}
	if ((mAirDrag==0)&&(mGravity.Norm()==0))
	{
		if (constraint==false)
		{
			while( dt*(frame+1) <= duration ){
				lPos   += lVel*dt;

				postrj(frame,0)=lPos(0);
				postrj(frame,1)=lPos(1);
				postrj(frame,2)=lPos(2);

				frame++;
			}
		}
		else
		{
			while(( dt*(frame+1) <= duration)&&(lPos(Direction)<value) ){
				lPos   += lVel*dt;

				postrj(frame,0)=lPos(0);
				postrj(frame,1)=lPos(1);
				postrj(frame,2)=lPos(2);

				frame++;
		//		cout<<"frame "<<frame<<" "<<lPos(Direction)<<endl;
			}
		}
	}
	else if (mAirDrag==0)
	{
		if (constraint==false)
		{
			while( dt*(frame+1) <= duration ){
				lVel   += mAirDrag*dt;
				lPos   += lVel*dt;

				postrj(frame,0)=lPos(0);
				postrj(frame,1)=lPos(1);
				postrj(frame,2)=lPos(2);

				frame++;
			}
		}
		else
		{
			while(( dt*(frame+1) <= duration)&&(lPos(Direction)>value) ){
				lVel   += mAirDrag*dt;
				lPos   += lVel*dt;

				postrj(frame,0)=lPos(0);
				postrj(frame,1)=lPos(1);
				postrj(frame,2)=lPos(2);

				frame++;
			}
		}

	}
	else
	{
		if (constraint==false)
		{
			while( dt*(frame+1) <= duration ){
				lAccel = mGravity - lVel*mAirDrag;
				lVel   += lAccel*dt;
				lPos   += lVel*dt + lAccel*dt*dt*0.5;

				postrj(frame,0)=lPos(0);
				postrj(frame,1)=lPos(1);
				postrj(frame,2)=lPos(2);

				frame++;
			}
		}
		else
		{
			while(( dt*(frame+1) <= duration)&&(lPos(Direction)>value) ){
				lAccel = mGravity - lVel*mAirDrag;
				lVel   += lAccel*dt;
				lPos   += lVel*dt + lAccel*dt*dt*0.5;

				postrj(frame,0)=lPos(0);
				postrj(frame,1)=lPos(1);
				postrj(frame,2)=lPos(2);

				frame++;
			}
		}
	}



	return frame;
}


void TrajectoryEstimation::PredictEndPos(double duration, MathLib::Vector3& pos)
{
	PredictEndPos(duration, mDt, pos);
}

void TrajectoryEstimation::PredictEndPos(double duration, double dt, MathLib::Vector3& pos )
{
	unsigned int frame=0;
	MathLib::Vector3 lVel;
	MathLib::Vector3 lAccel;
	MathLib::Vector3 lPos;

	lPos.Set(mCurrentPos);
	lVel.Set(mCurrentVel);
	if( duration < dt) duration = dt;

	while( dt*(double)(frame) < duration )
	{
		lAccel = mGravity - lVel*mAirDrag;
		lPos   += lVel*dt + lAccel*dt*dt*0.5;
		lVel   += lAccel*dt;

		frame++;
	}
	pos.Set(lPos);
}
void TrajectoryEstimation::PredictEndPosAndVel(double duration, double dt, MathLib::Vector3& pos, MathLib::Vector3& vel)
{
	unsigned int frame=0;
	MathLib::Vector3 lVel;
	MathLib::Vector3 lAccel;
	MathLib::Vector3 lPos;

	lPos.Set(mCurrentPos);
	lVel.Set(mCurrentVel);
	if( duration < dt) duration = dt;
	//	cout<<" lPos_Before "<<lPos(0);
	while( dt*(double)(frame) <1.0* duration )
	{
		lAccel = mGravity - lVel*mAirDrag;
		lPos   += lVel*dt + lAccel*dt*dt*0.5;
		lVel   += lAccel*dt;

		frame++;
	}
	//	cout<<" lPos_After "<<lPos(0);
	pos.Set(lPos);
	vel.Set(lVel);
}
void TrajectoryEstimation::PredictEndPosAndVel_Predicted(double duration, double dt, MathLib::Vector3& pos, MathLib::Vector3& vel)
{
	unsigned int frame=0;
	MathLib::Vector3 lVel;
	MathLib::Vector3 lAccel;
	MathLib::Vector3 lPos;

	lPos.Set(mCurrentPos);
	lVel.Set(mCurrentVel);
	if( duration < dt) duration = dt;
	//	cout<<" lPos_Before "<<lPos(0);
	while( dt*(double)(frame) < duration )
	{
		lAccel = mGravity - lVel*mAirDrag;
		lPos   += lVel*dt + lAccel*dt*dt*0.5;
		lVel   += lAccel*dt;

		frame++;
	}
	//	cout<<" lPos_After "<<lPos(0);
	pos.Set(lPos);
	vel.Set(lVel);
}


void TrajectoryEstimation::PredictNextOrient(double dt, MathLib::Matrix3& orient)
{
	mRinv = mCurrentR*mSpatialInertia.mInertiaMoment*mCurrentR.Transpose();	
	mAngularVelocity = mRinv.Inverse()*mCurrentAngularMomentum;

	mW.Zero();
	mW(0,1) = -mAngularVelocity(2);
	mW(1,0) =  mAngularVelocity(2);
	mW(0,2) =  mAngularVelocity(1);
	mW(2,0) = -mAngularVelocity(1);
	mW(1,2) = -mAngularVelocity(0);
	mW(2,1) =  mAngularVelocity(0);

	mRdot = mW*mCurrentR;
	mCurrentR = mCurrentR + mRdot*dt;
	mCurrentR.Normalize();

	orient.Set(mCurrentR);
}

int TrajectoryEstimation::PredictNextOrient(double duration, MathLib::Matrix& dirx, MathLib::Matrix& diry)
{
	return PredictNextOrient(duration, mDt, dirx, diry);
}

int TrajectoryEstimation::PredictNextOrient(double duration, double dt, MathLib::Matrix& dirx, MathLib::Matrix& diry)
{
	MathLib::Matrix lRInv(3,3);
	MathLib::Matrix3 lR;

	unsigned int frame=0;

	lR.Set(mCurrentR);
	while( dt*(frame+1) <= duration ){

		lRInv = lR*mSpatialInertia.mInertiaMoment*lR.Transpose();
		mAngularVelocity = lRInv.Inverse()*mCurrentAngularMomentum;

		mW.Zero();
		mW(0,1) = -mAngularVelocity(2);
		mW(1,0) =  mAngularVelocity(2);
		mW(0,2) =  mAngularVelocity(1);
		mW(2,0) = -mAngularVelocity(1);
		mW(1,2) = -mAngularVelocity(0);
		mW(2,1) =  mAngularVelocity(0);

		mRdot = mW*lR;
		lR = lR + mRdot*dt;
		lR.Normalize();

		dirx.SetRow(lR.GetColumn(0), frame);
		diry.SetRow(lR.GetColumn(1), frame);
		frame++;
	}
	return frame;
}

void TrajectoryEstimation::PredictEndOrient(double duration, MathLib::Vector3& dirx, MathLib::Vector3& diry)
{
	PredictEndOrient(duration, mDt, dirx, diry);
}

void TrajectoryEstimation::PredictEndOrient(double duration, double dt, MathLib::Vector3& dirx, MathLib::Vector3& diry)
{
	MathLib::Matrix lRInv(3,3);
	MathLib::Matrix3 lR;

	unsigned int frame=0;

	lR.Set(mCurrentR);
	while( dt*(double)(frame) < duration )
	{
		lRInv = lR*mSpatialInertia.mInertiaMoment*lR.Transpose();
		mAngularVelocity = lRInv.Inverse()*mCurrentAngularMomentum;

		mW.Zero();
		mW(0,1) = -mAngularVelocity(2);
		mW(1,0) =  mAngularVelocity(2);
		mW(0,2) =  mAngularVelocity(1);
		mW(2,0) = -mAngularVelocity(1);
		mW(1,2) = -mAngularVelocity(0);
		mW(2,1) =  mAngularVelocity(0);

		mRdot = mW*lR;
		lR = lR + mRdot*dt;
		lR.Normalize();

		frame++;
	}

	dirx.Set((MathLib::Vector3)lR.GetColumn(0));
	diry.Set((MathLib::Vector3)lR.GetColumn(1));

}


double TrajectoryEstimation::getCurrentTime()
{
	return mCurrentTime;
}
