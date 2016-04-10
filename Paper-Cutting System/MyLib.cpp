#include <GL/glut.h>
#include "MyLib.h"

using namespace Leap;

float radius = 0.f;
Vector newFingerBonePosition[5][4][2];
float newFingerBonePositionX = 0;
float newFingerBonePositionY = 0;

void DrawJoint( Vector joint )
{
	GLUquadric *myQuad = gluNewQuadric();

	glPushMatrix();
	glTranslatef( joint.x, joint.y, joint.z );
	gluSphere( myQuad, radius, 100, 100 );
	glPopMatrix();
}

void DrawBone( Vector prev, Vector next )
{
	glBegin( GL_LINES );
	glVertex3f( prev.x, prev.y, prev.z );
	glVertex3f( next.x, next.y, next.z );
	glEnd();

	// Draw joint except tip position of thumb, index and middle finger
	if ( next == newFingerBonePosition[0][3][1] || next == newFingerBonePosition[1][3][1] || next == newFingerBonePosition[2][3][1] )
		return;
	else
		DrawJoint( next );
}

void DrawHands( void )
{
	GLfloat mat_ambient[]   = { 0, 0, 1, 1 };
	GLfloat mat_diffuse[]   = { 0, 0, 1, 1 };
	GLfloat mat_specular[]  = { 1, 1, 1, 1 };
	GLfloat low_shininess[] = { 10 };

	glEnable( GL_LIGHTING );

	glMaterialfv( GL_FRONT,   GL_AMBIENT,   mat_ambient );
	glMaterialfv( GL_FRONT,   GL_DIFFUSE,   mat_diffuse );
	glMaterialfv( GL_FRONT,  GL_SPECULAR,  mat_specular );
	glMaterialfv( GL_FRONT, GL_SHININESS, low_shininess );

	Vector diff[5][4][2];
	float amt = 0.f;

	// Compute the difference of position between distal phalange of index finger and other finger bones
	for ( int j = 0; j < 5; j++ )
		for ( int k = 0; k < 4; k++ )
			for ( int w = 0; w < 2; w++ )
				diff[j][k][w] = fingerBonePosition[j][k][w] - fingerBonePosition[1][3][1];

	// Control size of small hands when changing horizontal viewpoint
	if ( m_azi_now == 0 )
	{
		amt    =  0.5f;
		radius = 0.13f;
	}
	else if ( m_azi_now ==  10 )
	{
		amt    =  1.5f;
		radius = 0.12f;
	}
	else if ( m_azi_now ==  20 )
	{
		amt    =  3.0f;
		radius = 0.11f;
	}
	else if ( m_azi_now ==  30 )
	{
		amt    =  4.5f;
		radius = 0.10f;
	}
	else if ( m_azi_now == 350 )
	{
		amt    =  2.0f;
		radius = 0.12f;
	}
	else if ( m_azi_now == 340 )
	{
		amt    =  3.8f;
		radius = 0.11f;
	}
	else if ( m_azi_now == 330 )
	{
		amt    =  5.6f;
		radius = 0.10f;
	}

	// Set new position of the distal phalange of index finger
	newFingerBonePosition[1][3][1].x = newPivot.getX();
	newFingerBonePosition[1][3][1].y = newPivot.getY();
	newFingerBonePosition[1][3][1].z = newPivot.getZ() - amt;

	newFingerBonePositionX = newFingerBonePosition[1][3][1].x;
	newFingerBonePositionY = newFingerBonePosition[1][3][1].y;

	// Compute new position of finger bones
	for ( int j = 0; j < 5; j++ )
		for ( int k = 0; k < 4; k++ )
			for ( int w = 0; w < 2; w++ )
				newFingerBonePosition[j][k][w] = newFingerBonePosition[1][3][1] + diff[j][k][w];

	// Draw finger bones
	for ( int j = 0; j < 5; j++ )
	{
		for ( int k = 0; k < 4; k++ )
		{
			if ( ( j == 1 || j == 2 || j == 3 ) && k == 0 )
				continue;

			DrawBone( newFingerBonePosition[j][k][0], newFingerBonePosition[j][k][1] );
		}
	}
	for ( int j = 0; j < 4; j++ )
		DrawBone( newFingerBonePosition[j][0][1], newFingerBonePosition[j + 1][0][1] );
	DrawJoint( newFingerBonePosition[4][0][0] );

	glDisable( GL_LIGHTING );


	// Highlight thumb, index and middle finger
	GLfloat mat_ambient2[] = { 1, 0, 0, 1 };
	GLfloat mat_diffuse2[] = { 1, 0, 0, 1 };

	glEnable( GL_LIGHTING );

	glMaterialfv( GL_FRONT,   GL_AMBIENT,  mat_ambient2 );
	glMaterialfv( GL_FRONT,   GL_DIFFUSE,  mat_diffuse2 );
	glMaterialfv( GL_FRONT,  GL_SPECULAR,  mat_specular );
	glMaterialfv( GL_FRONT, GL_SHININESS, low_shininess );

	DrawJoint( newFingerBonePosition[0][3][1] );
	DrawJoint( newFingerBonePosition[1][3][1] );
	DrawJoint( newFingerBonePosition[2][3][1] );

	glDisable( GL_LIGHTING );
}


const float convergence = 2000.0f;
float eyeSeparation;
const float aspectRatio = 1.f;
const float fov = 45.0f * 0.01745329251994329547f;
const float nearClippingDistance = 10.0f;
const float farClippingDistance = 20000.0f;

void ApplyLeftFrustum( int type )
{
	if ( type == 1 )
		eyeSeparation = 0.3f;
	else
		eyeSeparation = 1.0f;

	float top, bottom, left, right;

	top = nearClippingDistance * tan( fov / 2 );
	bottom = -top;

	float a = aspectRatio * tan( fov / 2 ) * convergence;
	float b = a - eyeSeparation / 2;
	float c = a + eyeSeparation / 2;

	left = -b * nearClippingDistance / convergence;
	right = c * nearClippingDistance / convergence;

	glMatrixMode( GL_PROJECTION );
	glPushMatrix();
	glLoadIdentity();
	glFrustum( left, right, bottom, top, nearClippingDistance, farClippingDistance );
	glPopMatrix();

	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	glTranslatef( eyeSeparation / 2, 0.0f, 0.0f );
}

void ApplyRightFrustum( int type )
{
	if ( type == 1 )
		eyeSeparation = 0.3f;
	else
		eyeSeparation = 1.f;

	float top, bottom, left, right;

	top = nearClippingDistance * tan( fov / 2 );
	bottom = -top;

	float a = aspectRatio * tan( fov / 2 ) * convergence;
	float b = a - eyeSeparation / 2;
	float c = a + eyeSeparation / 2;

	left = -c * nearClippingDistance / convergence;
	right = b * nearClippingDistance / convergence;

	glMatrixMode( GL_PROJECTION );
	glPushMatrix();
	glLoadIdentity();
	glFrustum( left, right, bottom, top, nearClippingDistance, farClippingDistance );
	glPopMatrix();

	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	glTranslatef( -eyeSeparation / 2, 0.0f, 0.0f );
}


void DrawScene( void )
{
	glEnable( GL_TEXTURE_2D );
	glBindTexture( GL_TEXTURE_2D, texture_floor );
	glBegin( GL_QUADS );
	glColor4f( 1, 1, 1, 1 );
	glTexCoord2d(  0,  0 ); glVertex3d(  60, -20, -60 );
	glTexCoord2d( 40,  0 ); glVertex3d( -60, -20, -60 );
	glTexCoord2d( 40, 40 ); glVertex3d( -60, -20,  60 );
	glTexCoord2d(  0, 40 ); glVertex3d(  60, -20,  60 );
	glEnd();
	glDisable( GL_TEXTURE_2D );

	glEnable( GL_TEXTURE_2D );
	glBindTexture( GL_TEXTURE_2D, texture_door );
	glBegin( GL_QUADS );
	glTexCoord2d(  0,  0 ); glVertex3d(  60, -20, -60 );
	glTexCoord2d( 40,  0 ); glVertex3d(  60, -20,  60 );
	glTexCoord2d( 40, 40 ); glVertex3d(  60,  80,  60 );
	glTexCoord2d(  0, 40 ); glVertex3d(  60,  80, -60 );
	glEnd();
	glDisable( GL_TEXTURE_2D );

	glEnable( GL_TEXTURE_2D );
	glBindTexture( GL_TEXTURE_2D, texture_wall );
	glBegin( GL_QUADS );
	glColor4f( 1, 1, 1, 1 );
	glTexCoord2d(  0,  0 ); glVertex3d(  60, -20,  60 );
	glTexCoord2d( 40,  0 ); glVertex3d( -60, -20,  60 );
	glTexCoord2d( 40, 40 ); glVertex3d( -60,  80,  60 );
	glTexCoord2d(  0, 40 ); glVertex3d(  60,  80,  60 );
	glEnd();
	glBegin( GL_QUADS );
	glTexCoord2d(  0,  0 ); glVertex3d( -60, -20,  60 );
	glTexCoord2d( 40,  0 ); glVertex3d( -60, -20, -60 );
	glTexCoord2d( 40, 40 ); glVertex3d( -60,  80, -60 );
	glTexCoord2d(  0, 40 ); glVertex3d( -60,  80,  60 );
	glEnd();
	glDisable( GL_TEXTURE_2D );

	glBegin( GL_QUADS );
	glColor4d( 249, 242, 186, 1 );
	glVertex3d(  60, 80,  60 );
	glVertex3d( -60, 80,  60 );
	glVertex3d( -60, 80, -60 );
	glVertex3d(  60, 80, -60 );
	glEnd();
}


void DrawInfo( int xStart, int yStart, char *string, std::string color )
{
	char ch;

	if ( color == "white" )
		glColor3f( 1, 1, 1 );
	else
		glColor3f( 0, 0, 0 );

	glRasterPos3f( xStart, yStart, 0 );
	for ( unsigned int i = 0; i < strlen(string); i++ )
	{
		ch = string[i];
		glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, ch );
	}
}

void DrawInfo2( int xStart, int yStart, char *string )
{
	char ch;

	glColor3f( 0, 0, 0 );

	glRasterPos3f( xStart, yStart, 0 );
	for ( unsigned int i = 0; i < strlen(string); i++ )
	{
		ch = string[i];
		glutBitmapCharacter( GLUT_BITMAP_9_BY_15, ch );
	}
}


void DrawLight( int timeCount, int m_glutScreenHeight )
{
	if ( timeCount >= 15 && timeCount <= 60 )
	{
		glColor4f( 1, 0, 0, 1 );
		glPushMatrix();
		glLoadIdentity();
		glTranslatef( 500, m_glutScreenHeight - 50, 0 );
		glTranslatef( 5, 5, 0 );
		DrawCircle();
		glPopMatrix();
	}
	if ( timeCount >= 30 && timeCount <= 60 )
	{
		glColor4f( 1, 1, 0, 1 );
		glPushMatrix();
		glLoadIdentity();
		glTranslatef( 500, m_glutScreenHeight - 50, 0 );
		glTranslatef( 30, 5, 0 );
		DrawCircle();
		glPopMatrix();
	}
	if ( timeCount >= 45 && timeCount <= 60 )
	{
		glColor4f( 0, 1, 0, 1 );
		glPushMatrix();
		glLoadIdentity();
		glTranslatef( 500, m_glutScreenHeight - 50, 0 );
		glTranslatef( 55, 5, 0 );
		DrawCircle();
		glPopMatrix();
	}
}

void DrawCircle( void )
{
	DrawFan( 10, 10, 10, 0 );
	DrawFan( 10, 10, 20, 10 );
	DrawFan( 10, 10, 10, 20 );
	DrawFan( 10, 10, 0, 10 );
}

void DrawRoundedRec( void )
{
	DrawFan( 20, 20, 0, 20 ); // up-left

	glBegin( GL_QUADS ); // middle-left
	glVertex2f(  0, 20 );
	glVertex2f( 20, 20 );
	glVertex2f( 20, 60 );
	glVertex2f(  0, 60 );
	glEnd();

	DrawFan( 20, 60, 20, 80 ); // down-left

	glBegin( GL_QUADS ); // middle-up
	glVertex2f( 20, 20 );
	glVertex2f( 20,  0 );
	glVertex2f( 70,  0 );
	glVertex2f( 70, 20 );
	glEnd();

	DrawFan( 70, 20, 70, 0 ); // up-right

	glBegin( GL_QUADS ); // middle-right
	glVertex2f( 70, 20 );
	glVertex2f( 90, 20 );
	glVertex2f( 90, 60 );
	glVertex2f( 70, 60 );
	glEnd();

	DrawFan( 70, 60, 90, 60 ); // down-right

	glBegin( GL_QUADS ); // middle-down
	glVertex2f( 20, 60 );
	glVertex2f( 70, 60 );
	glVertex2f( 70, 80 );
	glVertex2f( 20, 80 );
	glEnd();

	glBegin( GL_QUADS ); // middle
	glVertex2f( 20, 20 );
	glVertex2f( 70, 20 );
	glVertex2f( 70, 60 );
	glVertex2f( 20, 60 );
	glEnd();
}

void DrawFan( float v0X, float v0Y, float v1X, float v1Y )
{
	btVector3 trianglePoints[12];

	trianglePoints[0][0] = v0X;
	trianglePoints[0][1] = v0Y;
	trianglePoints[1][0] = v1X;
	trianglePoints[1][1] = v1Y;

	for ( int i = 1; i < 11; i++ )
	{
		btVector3 tmp = RotationMatrix( trianglePoints[i][0], trianglePoints[i][1], trianglePoints[0][0], trianglePoints[0][1] );

		trianglePoints[i + 1][0] = tmp[0];
		trianglePoints[i + 1][1] = tmp[1];
	}

	glBegin( GL_TRIANGLE_FAN );
	for ( int i = 0; i < 12; i++ )
		glVertex2f( trianglePoints[i][0], trianglePoints[i][1] );
	glEnd();
}

btVector3 RotationMatrix( float vec1X, float vec1Y, float vec2X, float vec2Y )
{
	btVector3 before( vec1X - vec2X, vec1Y - vec2Y, 0 );
	btVector3 after( 0, 0, 0 );
	float angle = PI * 9 / 180.0f;

	after[0] = cosf( angle ) * before.getX() - sinf( angle ) * before.getY() + vec2X;
	after[1] = sinf( angle ) * before.getX() + cosf( angle ) * before.getY() + vec2Y;

	return after;
}