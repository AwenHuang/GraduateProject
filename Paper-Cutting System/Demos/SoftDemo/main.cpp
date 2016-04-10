/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  http://continuousphysics.com/Bullet/

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#include "SoftDemo.h"
#include "GlutStuff.h"
#include "GLDebugDrawer.h"
#include "btBulletDynamicsCommon.h"

GLDebugDrawer   gDebugDrawer;

SoftDemo *softDemo;

#include <Windows.h>
#include "Leap.h"
#include "../../MyLib.h"
#include "../../command.h"

using namespace Leap;

Leap::Vector fingerBonePosition[5][4][2];

bool hasDisappeared = true;
bool hasPinched     = false;
bool hasCut         = false;
bool hasStarted     = false;
bool isRight        = false;

float newIndexPositionX;
float newIndexPositionY;

int64_t lastFrameId = 0;

int cur_teaching_mode = 0;

char outputInfo[1000];

const std::string fingerNames[] = {"Thumb", "Index", "Middle", "Ring", "Pinky"};
const std::string boneNames[]   = {"Metacarpal", "Proximal", "Middle", "Distal"};
const std::string stateNames[]  = {"STATE_INVALID", "STATE_START", "STATE_UPDATE", "STATE_END"};

bool isScissorhands(FingerList fingers)
{
	Vector thumbPos     = fingers[0].tipPosition();
	Vector indexPos     = fingers[1].tipPosition();
	Vector middlePos    = fingers[2].tipPosition();
	Vector ringPos      = fingers[3].tipPosition();
	Vector pinkyPos     = fingers[4].tipPosition();
	Vector index_pPos   = fingers[1].bone(Bone::TYPE_PROXIMAL).nextJoint();
	Vector middle_pPos  = fingers[2].bone(Bone::TYPE_PROXIMAL).nextJoint();

	if (indexPos.z < middle_pPos.z && middlePos.z < index_pPos.z)
		if ((index_pPos.y > ringPos.y && index_pPos.y > pinkyPos.y) && (middle_pPos.y > ringPos.y && middle_pPos.y > pinkyPos.y))
			if ((index_pPos.z < ringPos.z && index_pPos.z < pinkyPos.z) && (middle_pPos.z < ringPos.z && middle_pPos.z < pinkyPos.z))
				if (thumbPos.x > index_pPos.x && thumbPos.y < index_pPos.y)
					return true;
	return false;
}

bool isPinch(FingerList fingers)
{
	Vector thumbPos = fingers[0].tipPosition();
	Vector indexPos = fingers[1].tipPosition();
	float dist = thumbPos.distanceTo(indexPos);

	if (dist < 23.f)
		return true;
	else
		return false;
}

void SampleListener::onInit(const Controller &controller)
{
	std::cout << "Initialized" << std::endl;
	strcat(outputInfo, "Initialized.\n");
}

void SampleListener::onConnect(const Controller &controller)
{
	std::cout << "Connected" << std::endl;
	controller.enableGesture(Gesture::TYPE_CIRCLE);
	controller.enableGesture(Gesture::TYPE_KEY_TAP);
	controller.enableGesture(Gesture::TYPE_SCREEN_TAP);
	controller.enableGesture(Gesture::TYPE_SWIPE);
	strcat(outputInfo, "Connected.\n");
}

void SampleListener::onDisconnect(const Controller &controller)
{
	// Note: not dispatched when running in a debugger.
	std::cout << "Disconnected" << std::endl;
	strcat(outputInfo, "Disconnected.\n");
}

void SampleListener::onExit(const Controller &controller)
{
	std::cout << "Exited" << std::endl;
	strcat(outputInfo, "Exited.\n");
}

void SampleListener::onFocusGained(const Controller &controller)
{
	std::cout << "Focus Gained" << std::endl;
	strcat(outputInfo, "Focus Gained.\n");
}

void SampleListener::onFocusLost(const Controller &controller)
{
	std::cout << "Focus Lost" << std::endl;
	strcat(outputInfo, "Focus Lost.\n");
}

void SampleListener::onDeviceChange(const Controller &controller)
{
	std::cout << "Device Changed" << std::endl;
	const DeviceList devices = controller.devices();
	strcat(outputInfo, "Device Changed.\n");

	for (int i = 0; i < devices.count(); ++i)
	{
		std::string tmp;

		std::cout << "id: " << devices[i].toString() << std::endl;
		std::cout << "  isStreaming: " << (devices[i].isStreaming() ? "true" : "false") << std::endl;
		tmp = "id: " + devices[i].toString() + "\n" + "  isStreaming: " + (devices[i].isStreaming() ? "true" : "false") + "." + "\n";
		strcat(outputInfo, tmp.c_str());
	}
}

void SampleListener::onServiceConnect(const Controller &controller)
{
	std::cout << "Service Connected" << std::endl;
	strcat(outputInfo, "Initialized.\nService Connected.\n");
}

void SampleListener::onServiceDisconnect(const Controller &controller)
{
	std::cout << "Service Disconnected" << std::endl;
	strcat(outputInfo, "Service Disconnected.\n");
}

void SampleListener::onFrame(const Controller &controller)
{
	// Get the most recent frame and report some basic information
	const Frame frame = controller.frame();

	// Avoid processing the same Frame object twice
	if ( frame.id() == lastFrameId )
		return;
	else
		lastFrameId = frame.id();

	InteractionBox iBox = frame.interactionBox();

	// Get the height and width of screen
	RECT desktop;
	const HWND hDesktop = GetDesktopWindow();
	GetWindowRect(hDesktop, &desktop);
	int screenHeight = desktop.bottom;
	int screenWidth  = desktop.right;

	float indexPositionX = 0;
	float indexPositionY = 0;
	float middlePositionX = 0;
	float middlePositionY = 0;

	HandList hands = frame.hands();
	if ( hands.isEmpty() )
		hasDisappeared = true;

	for (HandList::const_iterator hl = hands.begin(); hl != hands.end(); ++hl)
	{
		hasDisappeared = false;

		// Get the first hand
		const Hand hand = *hl;

		if ( hand.isLeft() )
			isRight = false;
		else
			isRight = true;

        // Only deal with right hand control
		if ( isRight )
		{
			// Get fingers
			const FingerList fingers = hand.fingers();
			std::cout << hand.palmPosition() << std::endl;
			int i = 0;
			for (FingerList::const_iterator fl = fingers.begin(); fl != fingers.end(); ++fl)
			{
				const Finger finger = *fl;

				// Get finger bones
				for (int b = 0; b < 4; ++b)
				{
					Bone::Type boneType = static_cast<Bone::Type>(b);
					Bone bone = finger.bone(boneType);

					// Get the position of distal phalange of index finger
					if ( i == 1 && b == 3 )
					{
						Leap::Vector normalizePoint = iBox.normalizePoint(bone.nextJoint(), true);

						indexPositionX = normalizePoint.x * screenWidth;
						indexPositionY = (1 - normalizePoint.y) * screenHeight;
					}
					// Get the position of distal phalange of middle finger
					else if ( i == 2 && b == 3 )
					{
						Leap::Vector normalizePoint = iBox.normalizePoint(bone.nextJoint(), true);

						middlePositionX = normalizePoint.x * screenWidth;
						middlePositionY = (1 - normalizePoint.y) * screenHeight;
					}

					fingerBonePosition[i][b][0] = bone.prevJoint();
					fingerBonePosition[i][b][1] = bone.nextJoint();

					// Scale and shift every finger bone position
					for ( int c = 0; c < 2; c++ )
					{
						float amt;

						if ( m_azi_now == 0 )
							amt = 0.030f;
						else if ( m_azi_now == 10 || m_azi_now == 350 )
							amt = 0.028f;
						else if ( m_azi_now == 20 || m_azi_now == 340 )
							amt = 0.025f;
						else if ( m_azi_now == 30 || m_azi_now == 330 )
							amt = 0.021f;

						fingerBonePosition[i][b][c].x *= -amt;
						fingerBonePosition[i][b][c].y *=  amt;
						fingerBonePosition[i][b][c].z *= -amt;
					}
				}

				i++;
			} // end for

            // Control the cloth only when window is open, otherwise system may crash
			if ( hasCreatedWindow )
			{
				// Get the position of left and up corner of window
				int screenPositionX = glutGet((GLenum)GLUT_WINDOW_X);
				int screenPositionY = glutGet((GLenum)GLUT_WINDOW_Y);

				indexPositionX  -= screenPositionX;
				indexPositionY  -= screenPositionY;
				middlePositionX -= screenPositionX;
				middlePositionY -= screenPositionY;

				newIndexPositionX = indexPositionX;
				newIndexPositionY = indexPositionY;

				softDemo->mouseMotionFunc2( indexPositionX, indexPositionY );

				// Pinch the cloth
				if ( isPinch( fingers ) && !hasCut && !isScissorhands( fingers ) )
				{
					if ( !hasPinched )
					{
						softDemo->mouseFunc2( 0, indexPositionX, indexPositionY );
						hasPinched = true;
						strcat( outputInfo, "You pinched the paper!\n" );
					}
				}
				else
				{
					softDemo->mouseFunc2( 1, indexPositionX, indexPositionY );
					hasPinched = false;
				}
				
				// Cut the cloth
				if ( isScissorhands( fingers ) && !hasPinched && !isPinch( fingers ) )
				{
					float angle = fingers[1].direction().angleTo( fingers[2].direction() ) * 180 / Leap::PI;
					float dist  = fingers[1].tipPosition().distanceTo( fingers[2].tipPosition() );

					if ( !hasCut && angle > 5.f && dist > 20.f )
					{
						hasCut = true;
					}
					else if ( hasCut && angle < 5.f && dist < 20.f )
					{
						softDemo->mouseFunc2( 0, middlePositionX - 90, middlePositionY );
						softDemo->mouseFunc2( 1, middlePositionX - 90, middlePositionY );
						hasCut = false;
					}
				}
				else
					hasCut = false;
			} // end if ( hasCreatedWindow )
		} // end if ( isRight )
	} // end for

	// Get gestures
	const GestureList gestures = frame.gestures();
	for (int g = 0; g < gestures.count(); ++g)
	{
		Gesture gesture = gestures[g];

		switch (gesture.type())
		{
		// Change teaching mode by swipe
		case Gesture::TYPE_SWIPE:
			{
				if ( teach_mode == true )
				{
					SwipeGesture swipe = gesture;
					Vector direction = swipe.direction();
					Vector startPos  = swipe.startPosition();

					if ( gesture.state() == Gesture::STATE_START )
						hasStarted = true;
					if ( gesture.state() == Gesture::STATE_STOP && hasStarted )
					{
						cur_teaching_mode++;
						strcat( outputInfo, "You swiped!\n" );
						hasStarted = false;
					}
				}

				break;
			} // end case
		default:
			break;
		} // end switch
	} // end for
} // end onFrame

int main(int argc, char **argv)
{

	softDemo = new SoftDemo();

    // Create a sample listener and controller
	SampleListener listener;
	Controller controller;

	// Have the sample listener receive events from the controller
	controller.addListener(listener);

	if (argc > 1 && strcmp(argv[1], "--bg") == 0)
		controller.setPolicyFlags(Leap::Controller::POLICY_BACKGROUND_FRAMES);

	for (int i = 0; i < 500; i++)
		outputInfo[i] = '\0';
	softDemo->initPhysics();
	softDemo->getDynamicsWorld()->setDebugDrawer(&gDebugDrawer);

	glutmain(argc, argv, 1280, 768, "Paper-Cutting Simulation", softDemo);

	delete softDemo;
	// Remove the sample listener when done
	controller.removeListener(listener);
	return 0;

}
