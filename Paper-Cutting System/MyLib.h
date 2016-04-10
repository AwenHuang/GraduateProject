#ifndef MYLIB_H
#define MYLIB_H

#include "LeapSDK/include/Leap.h"
#include "src/LinearMath/btVector3.h"
#include "gl/glut.h"

using namespace Leap;

class SampleListener : public Listener
{
public:
	virtual void onInit(const Controller &);
	virtual void onConnect(const Controller &);
	virtual void onDisconnect(const Controller &);
	virtual void onExit(const Controller &);
	virtual void onFrame(const Controller &);
	virtual void onFocusGained(const Controller &);
	virtual void onFocusLost(const Controller &);
	virtual void onDeviceChange(const Controller &);
	virtual void onServiceConnect(const Controller &);
	virtual void onServiceDisconnect(const Controller &);

private:
};

extern GLuint texture_paper;
extern GLuint texture_floor;
extern GLuint texture_door;
extern GLuint texture_wall;

// Record the cameraPosition, cameraTargetPosition, and cameraUp
extern float cp[3];
extern float ctp[3];
extern float cu[3];

// Record whether the window has been created by OpenGL
extern bool hasCreatedWindow;

// Store the position of hands provided by Leap Motion
extern Vector fingerBonePosition[5][4][2];

// Record whether hands can be detected by Leap Motion
extern bool hasDisappeared;

// Determine whether the hand detected by Leap Motion is right hand
extern bool isRight;

// Store new position of the distal phalange of index finger
extern btVector3 newPivot;

extern float newFingerBonePositionX;
extern float newFingerBonePositionY;
extern float newIndexPositionX;
extern float newIndexPositionY;

extern int cur_teaching_mode;

extern char outputInfo[1000];

extern int m_azi_now;
extern int m_ele_now;

void DrawJoint( Vector joint );
void DrawBone( Vector prev, Vector next);
void DrawHands( void );

void ApplyLeftFrustum( int type );
void ApplyRightFrustum( int type );

void DrawScene( void );

void DrawInfo( int xStart, int yStart, char *string, std::string color );
void DrawInfo2( int xStart, int yStart, char *string );

void DrawLight( int timeCount, int m_glutScreenHeight );
void DrawCircle( void );
void DrawRoundedRec( void );
void DrawFan( float v0X, float v0Y, float v1X, float v1Y );
btVector3 RotationMatrix( float vec1X, float vec1Y, float vec2X, float vec2Y );

#endif