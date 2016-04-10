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

#ifndef _WINDOWS

#include "DemoApplication.h"

//glut is C code, this global gDemoApplication links glut to the C++ demo
static DemoApplication *gDemoApplication = 0;


#include "GlutStuff.h"

static  void glutKeyboardCallback(unsigned char key, int x, int y)
{
	gDemoApplication->keyboardCallback(key, x, y);
}

static  void glutKeyboardUpCallback(unsigned char key, int x, int y)
{
	gDemoApplication->keyboardUpCallback(key, x, y);
}

static void glutSpecialKeyboardCallback(int key, int x, int y)
{
	gDemoApplication->specialKeyboard(key, x, y);
}

static void glutSpecialKeyboardUpCallback(int key, int x, int y)
{
	gDemoApplication->specialKeyboardUp(key, x, y);
}


static void glutReshapeCallback(int w, int h)
{
	gDemoApplication->reshape(w, h);
}

static void glutMoveAndDisplayCallback()
{
	gDemoApplication->moveAndDisplay();
}

static void glutMouseFuncCallback(int button, int state, int x, int y)
{
	gDemoApplication->mouseFunc(button, state, x, y);
}


static void glutMotionFuncCallback(int x, int y)
{
	gDemoApplication->mouseMotionFunc(x, y);
}


static void glutDisplayCallback(void)
{
	gDemoApplication->displayCallback();
}

#include "../../MyLib.h"

GLuint texture_paper = 0;
GLuint texture_floor = 0;
GLuint texture_door  = 0;
GLuint texture_wall  = 0;

int texture_type = 0;

bool hasCreatedWindow = false;

unsigned char *LoadBitmapFile( char *fileName, BITMAPINFO *bitmapInfo )
{
	FILE *fp;
	BITMAPFILEHEADER bitmapFileHeader; // Bitmap file header
	unsigned char *bitmapImage; // Bitmap image data
	unsigned int lInfoSize; // Size of information
	unsigned int lBitSize; // Size of bitmap

	unsigned char change;
	int pixel;
	int p = 0;

	fp = fopen( fileName, "rb" );
	fread( &bitmapFileHeader, sizeof( BITMAPFILEHEADER ), 1, fp );

	lInfoSize = bitmapFileHeader.bfOffBits - sizeof( BITMAPFILEHEADER );
	fread( bitmapInfo, lInfoSize, 1, fp );

	lBitSize = bitmapInfo->bmiHeader.biSizeImage;
	bitmapImage = new BYTE[lBitSize];
	fread( bitmapImage, 1, lBitSize, fp );

	fclose( fp );

	pixel = ( bitmapInfo->bmiHeader.biWidth ) * ( bitmapInfo->bmiHeader.biHeight );

	for ( int i = 0; i < pixel; i++, p += 3 )
	{
		change = bitmapImage[p];
		bitmapImage[p] = bitmapImage[p + 2];
		bitmapImage[p + 2]  = change;
	}

	return bitmapImage;
}

void PrepareTexture( void )
{
	int width;
	int height;
	unsigned char *image;
	BITMAPINFO bmpinfo;

	image = LoadBitmapFile( "../../paper.bmp", &bmpinfo );
	width = bmpinfo.bmiHeader.biWidth;
	height = bmpinfo.bmiHeader.biHeight;
	glGenTextures( 1, &texture_paper );
	glBindTexture( GL_TEXTURE_2D, texture_paper );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
	glTexImage2D( GL_TEXTURE_2D, 0, 3, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image );

	image = LoadBitmapFile( "../../floor.bmp", &bmpinfo );
	width = bmpinfo.bmiHeader.biWidth;
	height = bmpinfo.bmiHeader.biHeight;
	glGenTextures( 1, &texture_floor );
	glBindTexture( GL_TEXTURE_2D, texture_floor );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
	glTexImage2D( GL_TEXTURE_2D, 0, 3, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image );

	image = LoadBitmapFile( "../../door.bmp", &bmpinfo );
	width = bmpinfo.bmiHeader.biWidth;
	height = bmpinfo.bmiHeader.biHeight;
	glGenTextures( 1, &texture_door );
	glBindTexture( GL_TEXTURE_2D, texture_door );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
	glTexImage2D( GL_TEXTURE_2D, 0, 3, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image );

	image = LoadBitmapFile( "../../wall.bmp", &bmpinfo );
	width = bmpinfo.bmiHeader.biWidth;
	height = bmpinfo.bmiHeader.biHeight;
	glGenTextures( 1, &texture_wall );
	glBindTexture( GL_TEXTURE_2D, texture_wall );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
	glTexImage2D( GL_TEXTURE_2D, 0, 3, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image );
}

int glutmain(int argc, char **argv, int width, int height, const char *title, DemoApplication *demoApp)
{

	gDemoApplication = demoApp;

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_STENCIL);
	glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH) - width) / 2, (glutGet(GLUT_SCREEN_HEIGHT) - height) / 2);
	glutInitWindowSize(width, height);
	glutCreateWindow(title);

	hasCreatedWindow = true;

#ifdef BT_USE_FREEGLUT
	glutSetOption (GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
#endif
	gDemoApplication->myinit();
	glutKeyboardFunc(glutKeyboardCallback);
	glutKeyboardUpFunc(glutKeyboardUpCallback);
	glutSpecialFunc(glutSpecialKeyboardCallback);
	glutSpecialUpFunc(glutSpecialKeyboardUpCallback);

	glutReshapeFunc(glutReshapeCallback);
	//glutFullScreen(); // FullScreen 
	//createMenu();
	glutIdleFunc(glutMoveAndDisplayCallback);
	glutMouseFunc(glutMouseFuncCallback);
	glutPassiveMotionFunc(glutMotionFuncCallback);
	glutMotionFunc(glutMotionFuncCallback);
	glutDisplayFunc( glutDisplayCallback );

	glutMoveAndDisplayCallback();

	//enable vsync to avoid tearing on Apple (todo: for Windows)

#if defined(__APPLE__) && !defined (VMDMESA)
	int swap_interval = 1;
	CGLContextObj cgl_context = CGLGetCurrentContext();
	CGLSetParameter(cgl_context, kCGLCPSwapInterval, &swap_interval);
#endif

	PrepareTexture();

	glutMainLoop();
	return 0;
}


#endif //_WINDOWS


