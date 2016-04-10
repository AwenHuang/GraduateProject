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

///btSoftBody implementation by Nathanael Presson


#include "btBulletDynamicsCommon.h"
#include "BulletSoftBody/btSoftRigidDynamicsWorld.h"

#include "BulletCollision/CollisionDispatch/btSphereSphereCollisionAlgorithm.h"
#include "BulletCollision/NarrowPhaseCollision/btGjkEpa2.h"
#include "LinearMath/btQuickprof.h"
#include "LinearMath/btIDebugDraw.h"

#include <stdio.h> //printf debugging
#include "LinearMath/btConvexHull.h"
#include "BulletSoftBody/btSoftBodyRigidBodyCollisionConfiguration.h"
#include "BulletSoftBody/btSoftBodyHelpers.h"

#include "SoftDemo.h"
#include "GL_ShapeDrawer.h"
#include "GLDebugFont.h"
#include "GlutStuff.h"
#include "../../command.h"
#include "../../MyLib.h"
#include "time.h"
#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)
//global section
extern btAlignedObjectArray<btVector3> old_position;

extern int cur_step;           //folding step    產生對角線後=1, 造成一摺後=2,二摺後=3,剪取後=4,展開後=5
int record = -1;
int index;
extern int current_mode;
bool teach_mode = false;
bool near_teach_node = false;
btVector3 *tmp;

extern int state;             //紀錄左右摺或上下摺 1:上下摺 2:左右摺 3:右上左下摺 4:左上右下摺
extern btScalar oldmass;
int fold_condition ;
extern int cutend;    //剪斷的狀態


/* ***** ***** ***** ***** ***** added code  end  ***** ***** ***** ***** ***** */
/*********************************一摺***********************************/
void SoftDemo::fold_first()
{
	btSoftBodyArray    &sbs = getSoftDynamicsWorld()->getSoftBodyArray();   //抓目前創造的物件
	btSoftBody *psb = sbs[0];
	int size = 9;   
	for (int i = 0; i < size; i++)
	{
		if (i == (size - 1 - i) || record == 4)break; //中心保持原來位置 對摺後仍一樣
		if ((psb->m_nodes[i].m_x - psb->m_nodes[size - i - 1].m_x).length2() < 0.18)   //當沿著對摺線摺過去時兩點靠近時自動貼齊
		{
			btScalar tmpmass1 = psb->m_nodes[i].m_im, tmpmass2 = psb->m_nodes[size - i - 1].m_im;

			psb->m_nodes[i].m_im = 0;			//m_im代表node的質量 歸零則固定住 這時將貼齊的兩點固定
			psb->m_nodes[size - i - 1].m_im = 0;

			if (tmpmass1 && tmpmass2)
			{
				if ((i >= 0 && i < 2) || (i > 4 && i <= 6)) //將這兩點的質量改為分給摺線的左右兩點
				{
					psb->m_nodes[i + 2].m_im = tmpmass1;
					psb->m_nodes[size - i - 1 - 2].m_im = tmpmass2;
				}
				else
				{
					psb->m_nodes[i - 2].m_im = tmpmass1;
					psb->m_nodes[size - i - 1 + 2].m_im = tmpmass2;
				}

			}

			psb->m_nodes[record].m_x = psb->m_nodes[size - record - 1].m_x; //將要對摺的兩點位置重合 record為當前滑鼠區塊位置
			if (record == 1 || record == 7)
			{
				if (record == 1)
				{
					state = 1; //下往上摺
				}
				else
				{
					state = 7; //上往下摺
				}
				psb->m_nodes[record - 1].m_x = psb->m_nodes[size - record - 2].m_x;  //將同一邊的點一併對摺過去
				psb->m_nodes[record + 1].m_x = psb->m_nodes[size - record].m_x;
			}

			else if (record == 3 || record == 5)
			{
				state = 3; //右往左摺
				if (record == 5)state = 5;
				psb->m_nodes[record - 3].m_x = psb->m_nodes[size - record - 4].m_x;
				psb->m_nodes[record + 3].m_x = psb->m_nodes[size - record + 2].m_x;
			}
			else if (record == 2 || record == 6)
			{

				if (record == 6)
				{
					state = 6; //右上左下摺
					psb->m_nodes[record + 1].m_x = psb->m_nodes[record - 1].m_x;
					psb->m_nodes[size - record ].m_x = psb->m_nodes[size - record - 2].m_x;
				}
				else
				{
					state = 2; //左下右上摺
					psb->m_nodes[record - 1].m_x = psb->m_nodes[record + 1].m_x;
					psb->m_nodes[size - record - 2].m_x = psb->m_nodes[size - record].m_x;
				}
			}
			else if (record == 0 || record == 8)
			{

				if (record == 8)
				{
					state = 8; //右上左下摺
					psb->m_nodes[record - 1].m_x = psb->m_nodes[size - record + 2].m_x;
					psb->m_nodes[record - 3].m_x = psb->m_nodes[size - record].m_x;
				}
				else
				{
					state = 4; //左下右上摺
					psb->m_nodes[record + 1].m_x = psb->m_nodes[size - record - 4].m_x;
					psb->m_nodes[record + 3].m_x = psb->m_nodes[size - record - 2].m_x;
				}
			}
			//printf("state:%d\n", state);
			m_results.fraction = 1.f;      
			m_drag = false;			//對摺後讓滑鼠鬆開
			cur_step = 2;
			fold_condition = 1;  //一摺結束 不管step如何
			//printf("一摺成功\n");
			strcat(outputInfo, "Fold successfully!\n" );

			for (int k = 0; k < psb->m_nodes.size(); k++)
			{
				psb->m_nodes[k].m_v = btVector3(0, 0, 0);  //m_v為node的速度，歸零讓他不亂跑
			}
		}

	}

}
void SoftDemo::fold_second()
{
	btSoftBodyArray    &sbs = getSoftDynamicsWorld()->getSoftBodyArray();
	btSoftBody *psb = sbs[0];

	if (state == 1 || state == 7) //上下摺
	{

		if (record == 2 || record == 5 || record == 8)   //上下摺後左往右摺，若X軸小於-2則直接摺過去
		{
			if(psb->m_nodes[record].m_x.getX() < -2)
			{
				psb->m_nodes[2].m_x = psb->m_nodes[0].m_x;
				psb->m_nodes[5].m_x = psb->m_nodes[3].m_x;
				psb->m_nodes[8].m_x = psb->m_nodes[0].m_x;
				psb->m_nodes[3].m_im = 0; psb->m_nodes[5].m_im = 0;  //將中間點固定住
				m_results.fraction = 1.f;
				m_drag = false;
				cur_step = 3;
				fold_condition = 2;
				printf("二摺成功\n");
				strcat( outputInfo, "Fold successfully!\n" );

			}
			
		}
		else if (record == 0 || record == 3 || record == 6)
		{
			if (psb->m_nodes[record].m_x.getX() > 2)   //右往左摺
			{
				psb->m_nodes[0].m_x = psb->m_nodes[2].m_x;
				psb->m_nodes[3].m_x = psb->m_nodes[5].m_x;
				psb->m_nodes[6].m_x = psb->m_nodes[8].m_x;
				psb->m_nodes[3].m_im = 0; psb->m_nodes[5].m_im = 0;
				m_results.fraction = 1.f;
				m_drag = false;
				cur_step = 3;
				fold_condition = 2;
				printf("二摺成功\n");
				strcat( outputInfo, "Fold successfully!\n" );
			}
			
		}
		
		
		
	}
	else if (state == 3 || state == 5) //左右摺
	{
		if (record <= 2)
		{
			if (psb->m_nodes[record].m_x.getY() > 2)   //下往上摺
			{
				psb->m_nodes[0].m_x = psb->m_nodes[8].m_x;
				psb->m_nodes[1].m_x = psb->m_nodes[7].m_x;
				psb->m_nodes[2].m_x = psb->m_nodes[8].m_x;
				psb->m_nodes[0].m_im = 0; psb->m_nodes[8].m_im = 0;
				psb->m_nodes[1].m_im = 0; psb->m_nodes[7].m_im = 0;
				psb->m_nodes[2].m_im = 0; psb->m_nodes[6].m_im = 0;
				psb->m_nodes[3].m_im = oldmass;
				psb->m_nodes[4].m_im = oldmass;
				psb->m_nodes[5].m_im = oldmass;
				m_results.fraction = 1.f;
				m_drag = false;
				cur_step = 3;
				fold_condition = 2;
				printf("二摺成功\n");
				strcat( outputInfo, "Fold successfully!\n" );
			}
			
		}

		else if (record >= 6)
		{
			if (psb->m_nodes[record].m_x.getY() < -2)   //上往下摺
			{
				psb->m_nodes[6].m_x = psb->m_nodes[0].m_x;
				psb->m_nodes[7].m_x = psb->m_nodes[1].m_x;
				psb->m_nodes[8].m_x = psb->m_nodes[0].m_x;
				m_results.fraction = 1.f;
				m_drag = false;
				cur_step = 3;
				fold_condition = 2;
				printf("二摺成功\n");
				strcat( outputInfo, "Fold successfully!\n" );
			}
			
		}
		
	}

	else if (state == 4 || state == 8) //左上右下摺
	{
			if (psb->m_nodes[6].m_x.getX() > 1.4 && psb->m_nodes[6].m_x.getY() < -1.4) //上往下摺
			{
				psb->m_nodes[3].m_x = psb->m_nodes[1].m_x;
				psb->m_nodes[6].m_x = psb->m_nodes[2].m_x;
				psb->m_nodes[7].m_x = psb->m_nodes[5].m_x;
				m_results.fraction = 1.f;
				m_drag = false;
				cur_step = 3;
				fold_condition = 2;
				printf("二摺成功\n");
				strcat( outputInfo, "Fold successfully!\n" );
			}
			if (psb->m_nodes[2].m_x.getX() < -1.4 && psb->m_nodes[2].m_x.getY() > 1.4) //下往上摺
			{
				psb->m_nodes[1].m_x = psb->m_nodes[3].m_x;
				psb->m_nodes[2].m_x = psb->m_nodes[6].m_x;
				psb->m_nodes[5].m_x = psb->m_nodes[7].m_x;
				m_results.fraction = 1.f;
				m_drag = false;
				cur_step = 3;
				fold_condition = 2;
				psb->m_nodes[0].m_im = 0; psb->m_nodes[2].m_im = 0;
				psb->m_nodes[6].m_im = 0; psb->m_nodes[8].m_im = 0;
				psb->m_nodes[4].m_im = oldmass;
				printf("二摺成功\n");
				strcat( outputInfo, "Fold successfully!\n" );
			}

	}
	else if (state == 2 || state == 6) //右上左下摺
	{
		
			if (psb->m_nodes[8].m_x.getX() < -1.4 && psb->m_nodes[8].m_x.getY() < -1.4) //上往下摺
			{
				psb->m_nodes[5].m_x = psb->m_nodes[1].m_x;
				psb->m_nodes[7].m_x = psb->m_nodes[3].m_x;
				psb->m_nodes[8].m_x = psb->m_nodes[0].m_x;
				psb->m_nodes[0].m_im = 0; psb->m_nodes[8].m_im = 0;
				m_results.fraction = 1.f;
				m_drag = false;
				cur_step = 3;
				fold_condition = 2;
				printf("二摺成功\n");
				strcat( outputInfo, "Fold successfully!\n" );
			}

			if (psb->m_nodes[0].m_x.getX() > 1.4 && psb->m_nodes[0].m_x.getY() > 1.4) //下往上摺
			{
				psb->m_nodes[0].m_x = psb->m_nodes[8].m_x;
				psb->m_nodes[1].m_x = psb->m_nodes[7].m_x;
				psb->m_nodes[3].m_x = psb->m_nodes[5].m_x;
				psb->m_nodes[0].m_im = 0; psb->m_nodes[8].m_im = 0;
				m_results.fraction = 1.f;
				m_drag = false;
				cur_step = 3;
				fold_condition = 2;
				printf("二摺成功\n");
				strcat( outputInfo, "Fold successfully!\n" );
			}
		
	}

	for (int k = 0; k < psb->m_nodes.size(); k++)
		psb->m_nodes[k].m_v = btVector3(0, 0, 0);

}

/*********************************展開***********************************/
void SoftDemo::unfold()
{
	btSoftBodyArray    &sbs = getSoftDynamicsWorld()->getSoftBodyArray();
	btSoftBody *psb = sbs[0];
	int i, size = psb->m_nodes.size();
	teach_mode = false;
	
	if (fold_condition == 1)     //一摺後
	{
		for (i = 0; i < 9; i++)
		{
			psb->m_nodes[i].m_x = old_position[i];
			psb->m_nodes[i].m_im = oldmass;
			psb->m_nodes[i].m_v = btVector3(0, 0, 0);
		}

		if (cur_step == 2) //二摺前
			psb->m_nodes[4].m_im = 0;

		if (cur_step == 4)
		{
			if(psb->m_nodes.size()<=9)  //將教學模式未剪返回原展開狀態
			{
				cur_step=0;
				for (i = 0; i < 9; i++)
				{
					psb->m_nodes[i].m_x = old_position[i];
					psb->m_nodes[i].m_im = oldmass;
					psb->m_nodes[i].m_v = btVector3(0, 0, 0);
				}
				psb->m_nodes[4].m_im =0;
			}
			//if(cur_step!=0)
			else 
			{
			cur_step = 5;
			for (i = 0; i < 9; i++)
			{
				if (i == 4)continue;
				psb->m_nodes[i].m_im = 0;
			}
			for (i = 0; i < psb->m_nodes.size(); ++i)
			{
				btSoftBody::Node &n = psb->m_nodes[i];
				if (n.tag == 1)          //展開前在下方的點
				{
					
					n.m_im = 0;
					if (state == 1)
					{
						n.m_x.setY(-n.m_x.getY());  //若是下往上對摺的情況 展開將下方的點照Y軸鏡像反射
					}
					else if (state == 5)
					{	
						n.m_x.setX(-n.m_x.getX());
					}
				}
				if (n.foldtag == 1)    //展開前在上方的點
				{
					n.m_im = 0;
					if (state == 7)
					{
						n.m_x.setY(-n.m_x.getY());
					}
					else if (state == 3)
					{
						n.m_x.setX(-n.m_x.getX());
					}
				}

			}
		}
		}
		else    cur_step = 0;
		fold_condition = 0;
		state = 0;  //歸零

	}
	else if (fold_condition == 2) //二摺後
	{
		if (state == 1)       //返回第一摺的情形(由下往上摺)
		{
			psb->m_nodes[0].m_x = old_position[6];
			psb->m_nodes[2].m_x = old_position[8];
			psb->m_nodes[3].m_x = old_position[3];
			psb->m_nodes[5].m_x = old_position[5];
			psb->m_nodes[6].m_x = old_position[6];
			psb->m_nodes[8].m_x = old_position[8];
			fold_condition = 1;
			if (cur_step == 3) //剪取前
			{
				for (int i = 0; i < 9; i++)
					psb->m_nodes[i].m_im = oldmass;
				psb->m_nodes[1].m_im = 0; psb->m_nodes[4].m_im = 0; psb->m_nodes[7].m_im = 0;
				cur_step--;
			}
		}
		else if (state == 2)
		{
			psb->m_nodes[0].m_x = old_position[0];
			psb->m_nodes[3].m_x = old_position[3];
			psb->m_nodes[7].m_x = old_position[7];
			psb->m_nodes[8].m_x = old_position[8];

			fold_condition = 1;
			if (cur_step == 3) //剪取前
			{
				cur_step--;
			}

		}
		else if (state == 3)
		{
			psb->m_nodes[6].m_x = old_position[8];
			psb->m_nodes[8].m_x = old_position[8];
			psb->m_nodes[1].m_x = old_position[1];
			psb->m_nodes[7].m_x = old_position[7];
			psb->m_nodes[0].m_x = old_position[2];
			psb->m_nodes[2].m_x = old_position[2];

			fold_condition = 1;
			if (cur_step == 3) //剪取前
			{
				for (int i = 0; i < 9; i++)
					psb->m_nodes[i].m_im = oldmass;
				psb->m_nodes[3].m_im = 0; psb->m_nodes[4].m_im = 0; psb->m_nodes[5].m_im = 0;
				cur_step--;
			}
		}
		else if (state == 4)
		{
			psb->m_nodes[2].m_x = old_position[2];
			psb->m_nodes[5].m_x = old_position[5];
			psb->m_nodes[6].m_x = old_position[6];
			psb->m_nodes[7].m_x = old_position[7];
			fold_condition = 1;
			if (cur_step == 3) //剪取前
			{
				for (int i = 0; i < 9; i++)
					psb->m_nodes[i].m_im = oldmass;
				psb->m_nodes[0].m_im = 0; psb->m_nodes[4].m_im = 0; psb->m_nodes[8].m_im = 0;
				cur_step--;
			}

		}
		else if (state == 5)
		{
			psb->m_nodes[6].m_x = old_position[6];
			psb->m_nodes[8].m_x = old_position[6];
			psb->m_nodes[1].m_x = old_position[1];
			psb->m_nodes[7].m_x = old_position[7];
			psb->m_nodes[0].m_x = old_position[0];
			psb->m_nodes[2].m_x = old_position[0];
			fold_condition = 1;
			if (cur_step == 3) //剪取前
			{
				for (int i = 0; i < 9; i++)
					psb->m_nodes[i].m_im = oldmass;
				psb->m_nodes[3].m_im = 0; psb->m_nodes[4].m_im = 0; psb->m_nodes[5].m_im = 0;
				cur_step--;
			}
		}
		else if (state == 6)
		{
			psb->m_nodes[0].m_x = old_position[0];
			psb->m_nodes[1].m_x = old_position[1];
			psb->m_nodes[5].m_x = old_position[5];
			psb->m_nodes[8].m_x = old_position[8];
			fold_condition = 1;
			if (cur_step == 3) //剪取前
			{
				cur_step--;
			}

		}
		else if (state == 7)
		{
			psb->m_nodes[6].m_x = old_position[0];
			psb->m_nodes[8].m_x = old_position[2];
			psb->m_nodes[3].m_x = old_position[3];
			psb->m_nodes[5].m_x = old_position[5];
			psb->m_nodes[0].m_x = old_position[0];
			psb->m_nodes[2].m_x = old_position[2];
			fold_condition = 1;
			if (cur_step == 3) //剪取前
			{
				for (int i = 0; i < 9; i++)
					psb->m_nodes[i].m_im = oldmass;
				psb->m_nodes[1].m_im = 0; psb->m_nodes[4].m_im = 0; psb->m_nodes[7].m_im = 0;
				cur_step--;
			}
		}
		else if (state == 8)
		{
			psb->m_nodes[2].m_x = old_position[2];
			psb->m_nodes[1].m_x = old_position[1];
			psb->m_nodes[3].m_x = old_position[3];
			psb->m_nodes[6].m_x = old_position[6];
			fold_condition = 1;
			if (cur_step == 3) //剪取前
			{
				for (int i = 0; i < 9; i++)
					psb->m_nodes[i].m_im = oldmass;
				psb->m_nodes[0].m_im = 0; psb->m_nodes[4].m_im = 0; psb->m_nodes[8].m_im = 0;
				cur_step--;
			}

		}

	}

}

extern float eye[3];
extern int glutScreenWidth;
extern int glutScreenHeight;

static bool sDemoMode = false;

const int maxProxies = 32766;
const int maxOverlap = 65535;

static btVector3   *gGroundVertices = 0;
static int *gGroundIndices = 0;
static btBvhTriangleMeshShape *trimeshShape = 0;
static btRigidBody *staticBody = 0;
static float waveheight = 5.f;

const float TRIANGLE_SIZE = 8.f;
int     current_demo = 20;
#define DEMO_MODE_TIMEOUT 15.f //15 seconds for each demo


#ifdef _DEBUG
const int gNumObjects = 1;
#else
const int gNumObjects = 1;//try this in release mode: 3000. never go above 16384, unless you increate maxNumObjects  value in DemoApplication.cp
#endif

const int maxNumObjects = 32760;

#define CUBE_HALF_EXTENTS 1.5
#define EXTRA_HEIGHT -10.f


#ifdef USE_AMD_OPENCL
#include "btOpenCLUtils.h"
#include "BulletMultiThreaded/GpuSoftBodySolvers/OpenCL/btSoftBodySolver_OpenCL.h"
#include "BulletMultiThreaded/GpuSoftBodySolvers/OpenCL/btSoftBodySolver_OpenCLSIMDAware.h"
#include "BulletMultiThreaded/GpuSoftBodySolvers/OpenCL/btSoftBodySolverVertexBuffer_OpenGL.h"

btOpenCLSoftBodySolver *g_openCLSIMDSolver = 0;
btSoftBodySolverOutputCLtoCPU *g_softBodyOutput = 0;

cl_context          g_cxMainContext;
cl_device_id        g_cdDevice;
cl_command_queue    g_cqCommandQue;

void initCL( void *glCtx, void *glDC )
{
	int ciErrNum = 0;

#if defined(CL_PLATFORM_MINI_CL)
	cl_device_type deviceType = CL_DEVICE_TYPE_CPU;//or use CL_DEVICE_TYPE_DEBUG to debug MiniCL
#elif defined(CL_PLATFORM_INTEL)
	cl_device_type deviceType = CL_DEVICE_TYPE_CPU;
#elif defined(CL_PLATFORM_AMD)
	cl_device_type deviceType = CL_DEVICE_TYPE_GPU;
#elif defined(CL_PLATFORM_NVIDIA)
	cl_device_type deviceType = CL_DEVICE_TYPE_GPU;
#else
#ifdef __APPLE__
	cl_device_type deviceType = CL_DEVICE_TYPE_ALL;//GPU;
#else
	cl_device_type deviceType = CL_DEVICE_TYPE_CPU;//CL_DEVICE_TYPE_ALL
#endif//__APPLE__
#endif

	g_cxMainContext = btOpenCLUtils::createContextFromType(deviceType, &ciErrNum, glCtx, glDC);
	oclCHECKERROR(ciErrNum, CL_SUCCESS);


	int numDev = btOpenCLUtils::getNumDevices(g_cxMainContext);
	if (!numDev)
	{
		btAssert(0);
		exit(0);//this is just a demo, exit now
	}

	g_cdDevice = btOpenCLUtils::getDevice(g_cxMainContext, 0);
	oclCHECKERROR(ciErrNum, CL_SUCCESS);

	btOpenCLDeviceInfo clInfo;
	btOpenCLUtils::getDeviceInfo(g_cdDevice, clInfo);
	btOpenCLUtils::printDeviceInfo(g_cdDevice);

	// create a command-queue
	g_cqCommandQue = clCreateCommandQueue(g_cxMainContext, g_cdDevice, 0, &ciErrNum);
	oclCHECKERROR(ciErrNum, CL_SUCCESS);
}

class CachingCLFunctions : public CLFunctions
{
protected:

	cl_device_id        m_device;

	const char *strip(const char *name, const char *pattern);

public:
	CachingCLFunctions(cl_command_queue cqCommandQue, cl_context cxMainContext) :
		CLFunctions(cqCommandQue, cxMainContext)
	{
		size_t actualSize;
		cl_int retval = clGetCommandQueueInfo ( cqCommandQue, CL_QUEUE_DEVICE, sizeof(cl_device_id),
			&m_device, &actualSize);
	}

	/**
	* Compile a compute shader kernel from a string and return the appropriate cl_kernel object.
	*/
	virtual cl_kernel compileCLKernelFromString( const char *kernelSource, const char *kernelName, const char *additionalMacros , const char *orgSrcFileNameForCaching)
	{
		char srcFileNameForCaching[1024];
		sprintf(srcFileNameForCaching, "%s/%s", "../../src/BulletMultiThreaded/GpuSoftBodySolvers/OpenCL", orgSrcFileNameForCaching);

		btAssert(additionalMacros);
		btAssert(srcFileNameForCaching && strlen(srcFileNameForCaching));

		printf("compiling kernelName: %s ", kernelName);
		cl_kernel kernel = 0;
		cl_int ciErrNum;


		size_t program_length = strlen(kernelSource);

		cl_program m_cpProgram = btOpenCLUtils::compileCLProgramFromString(m_cxMainContext, m_device, kernelSource,  &ciErrNum, additionalMacros);


		// Create the kernel
		kernel = clCreateKernel(m_cpProgram, kernelName, &ciErrNum);
		if (ciErrNum != CL_SUCCESS)
		{
			const char *msg = "";
			switch (ciErrNum)
			{
			case CL_INVALID_PROGRAM:
				msg = "Program is not a valid program object.";
				break;
			case CL_INVALID_PROGRAM_EXECUTABLE:
				msg = "There is no successfully built executable for program.";
				break;
			case CL_INVALID_KERNEL_NAME:
				msg = "kernel_name is not found in program.";
				break;
			case CL_INVALID_KERNEL_DEFINITION:
				msg = "the function definition for __kernel function given by kernel_name such as the number of arguments, the argument types are not the same for all devices for which the program executable has been built.";
				break;
			case CL_INVALID_VALUE:
				msg = "kernel_name is NULL.";
				break;
			case CL_OUT_OF_HOST_MEMORY:
				msg = "Failure to allocate resources required by the OpenCL implementation on the host.";
				break;
			default:
				{
				}
			}

			printf("Error in clCreateKernel for kernel '%s', error is \"%s\", Line %u in file %s !!!\n\n", kernelName, msg, __LINE__, __FILE__);

#ifndef BT_SUPPRESS_OPENCL_ASSERTS
			btAssert(0);
#endif //BT_SUPPRESS_OPENCL_ASSERTS
			m_kernelCompilationFailures++;
			return 0;
		}

		printf("ready. \n");
		if (!kernel)
			m_kernelCompilationFailures++;
		return kernel;
	}

};


#endif //USE_AMD_OPENCL

//
void SoftDemo::createStack( btCollisionShape *boxShape, float halfCubeSize, int size, float zPos )
{
	btTransform trans;
	trans.setIdentity();

	for (int i = 0; i < size; i++)
	{
		// This constructs a row, from left to right
		int rowSize = size - i;
		for (int j = 0; j < rowSize; j++)
		{
			btVector3 pos;
			pos.setValue(
				-rowSize * halfCubeSize + halfCubeSize + j * 2.0f * halfCubeSize,
				halfCubeSize + i * halfCubeSize * 2.0f,
				zPos);

			trans.setOrigin(pos);
			btScalar mass = 1.f;

			btRigidBody *body = 0;
			body = localCreateRigidBody(mass, trans, boxShape);

		}
	}
}


////////////////////////////////////

extern int gNumManifold;
extern int gOverlappingPairs;

///for mouse picking
void pickingPreTickCallback (btDynamicsWorld *world, btScalar timeStep)
{
	SoftDemo *softDemo = (SoftDemo *)world->getWorldUserInfo();

	if (softDemo->m_drag)
	{
		const int               x = softDemo->m_lastmousepos[0];
		const int               y = softDemo->m_lastmousepos[1];
		const btVector3         rayFrom = softDemo->getCameraPosition();
		const btVector3         rayTo = softDemo->getRayTo(x, y);
		const btVector3         rayDir = (rayTo - rayFrom).normalized();
		const btVector3         N = (softDemo->getCameraTargetPosition() - softDemo->getCameraPosition()).normalized();
		const btScalar          O = btDot(softDemo->m_impact, N);
		const btScalar          den = btDot(N, rayDir);
		if ((den * den) > 0)
		{
			const btScalar          num = O - btDot(N, rayFrom);
			const btScalar          hit = num / den;
			if ((hit > 0) && (hit < 1500))
			{
				softDemo->m_goal = rayFrom + rayDir * hit;
			}
		}
		btVector3               delta = softDemo->m_goal - softDemo->m_node->m_x;
		static const btScalar   maxdrag = 10;
		if (delta.length2() > (maxdrag * maxdrag))
		{
			delta = delta.normalized() * maxdrag;
		}
		softDemo->m_node->m_v += delta / timeStep;
	}

}



void SoftDemo::displayCallback(void)
{

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	renderme();

	glFlush();
	swapBuffers();
}


//
// ImplicitShape
//

//
struct  ImplicitSphere : btSoftBody::ImplicitFn
{
	btVector3   impact;		 //滑鼠點擊位置
	btScalar    len;		 //長度
	ImplicitSphere() {}
	ImplicitSphere(const btVector3 &c, btScalar r) : impact(c), len(r) {}
	btScalar    Eval(const btVector3 &x)
	{
		if (current_mode == 1)
			return ((x - impact).length2() - len * len);		//以impact為圓心 len為半徑剪圓
		else if (current_mode == 2)
			return impact.getY() - x.getY();					//與impact的Y軸相同則截去
		else if (current_mode == 5)
			return impact.getX() - x.getX();					//與impact的X軸相同則截去
		
	}
};
//
// Random
//

static inline btScalar  UnitRand()
{
	return (rand() / (btScalar)RAND_MAX);
}

static inline btScalar  SignedUnitRand()
{
	return (UnitRand() * 2 - 1);
}

static inline btVector3 Vector3Rand()
{
	const btVector3 p = btVector3(SignedUnitRand(), SignedUnitRand(), SignedUnitRand());
	return (p.normalized());
}

//
// Rb rain
//
static void Ctor_RbUpStack(SoftDemo *pdemo, int count)
{
	float               mass = 10;

	btCompoundShape *cylinderCompound = new btCompoundShape;
	btCollisionShape *cylinderShape = new btCylinderShapeX(btVector3(4, 1, 1));
	btCollisionShape *boxShape = new btBoxShape(btVector3(4, 1, 1));
	btTransform localTransform;
	localTransform.setIdentity();
	cylinderCompound->addChildShape(localTransform, boxShape);
	btQuaternion orn(SIMD_HALF_PI, 0, 0);
	localTransform.setRotation(orn);
	//  localTransform.setOrigin(btVector3(1,1,1));
	cylinderCompound->addChildShape(localTransform, cylinderShape);


	btCollisionShape   *shape[] = {cylinderCompound,
		new btBoxShape(btVector3(1, 1, 1)),
		new btSphereShape(1.5)

	};
	static const int    nshapes = sizeof(shape) / sizeof(shape[0]);
	for (int i = 0; i < count; ++i)
	{
		btTransform startTransform;
		startTransform.setIdentity();
		startTransform.setOrigin(btVector3(0, 2 + 6 * i, 0));
		pdemo->localCreateRigidBody(mass, startTransform, shape[i % nshapes]);
		//pdemo->localCreateRigidBody(mass,startTransform,shape[0]);
	}
}

//
// Big ball
//
static void Ctor_BigBall(SoftDemo *pdemo, btScalar mass = 10)
{
	btTransform startTransform;
	startTransform.setIdentity();
	startTransform.setOrigin(btVector3(0, 13, 0));
	pdemo->localCreateRigidBody(mass, startTransform, new btSphereShape(3));
}

//
// Big plate
//
static btRigidBody *Ctor_BigPlate(SoftDemo *pdemo, btScalar mass = 15, btScalar height = 4)
{
	btTransform startTransform;
	startTransform.setIdentity();
	startTransform.setOrigin(btVector3(0, height, 0.5));
	btRigidBody        *body = pdemo->localCreateRigidBody(mass, startTransform, new btBoxShape(btVector3(5, 1, 5)));
	body->setFriction(1);
	return (body);
}

//
// Linear stair
//
static void Ctor_LinearStair(SoftDemo *pdemo, const btVector3 &org, const btVector3 &sizes, btScalar angle, int count)
{
	btBoxShape *shape = new btBoxShape(sizes);
	for (int i = 0; i < count; ++i)
	{
		btTransform startTransform;
		startTransform.setIdentity();
		startTransform.setOrigin(org + btVector3(sizes.x()*i * 2, sizes.y()*i * 2, 0));
		btRigidBody *body = pdemo->localCreateRigidBody(0, startTransform, shape);
		body->setFriction(1);
	}
}

//
// Softbox
//
static btSoftBody *Ctor_SoftBox(SoftDemo *pdemo, const btVector3 &p, const btVector3 &s)
{
	const btVector3 h = s * 0.5;
	const btVector3 c[] = {   p + h * btVector3(-1, -1, -1),
		p + h * btVector3(+1, -1, -1),
		p + h * btVector3(-1, +1, -1),
		p + h * btVector3(+1, +1, -1),
		p + h * btVector3(-1, -1, +1),
		p + h * btVector3(+1, -1, +1),
		p + h * btVector3(-1, +1, +1),
		p + h * btVector3(+1, +1, +1)
	};
	btSoftBody     *psb = btSoftBodyHelpers::CreateFromConvexHull(pdemo->m_softBodyWorldInfo, c, 8);
	psb->generateBendingConstraints(2);
	pdemo->getSoftDynamicsWorld()->addSoftBody(psb);

	return (psb);

}

//
// SoftBoulder
//
static btSoftBody *Ctor_SoftBoulder(SoftDemo *pdemo, const btVector3 &p, const btVector3 &s, int np, int id)
{
	btAlignedObjectArray<btVector3> pts;
	if (id) srand(id);
	for (int i = 0; i < np; ++i)
	{
		pts.push_back(Vector3Rand()*s + p);
	}
	btSoftBody     *psb = btSoftBodyHelpers::CreateFromConvexHull(pdemo->m_softBodyWorldInfo, &pts[0], pts.size());
	psb->generateBendingConstraints(2);
	pdemo->getSoftDynamicsWorld()->addSoftBody(psb);

	return (psb);
}

//#define TRACEDEMO { pdemo->demoname=__FUNCTION__+5;printf("Launching demo: " __FUNCTION__ "\r\n"); }

//
// Basic ropes
//
static void Init_Ropes(SoftDemo *pdemo)
{
	//TRACEDEMO
	const int n = 15;
	for (int i = 0; i < n; ++i)
	{
		btSoftBody *psb = btSoftBodyHelpers::CreateRope(pdemo->m_softBodyWorldInfo,   btVector3(-10, 0, i * 0.25),
			btVector3(10, 0, i * 0.25),
			16,
			1 + 2);
		psb->m_cfg.piterations      =   4;
		psb->m_materials[0]->m_kLST =   0.1 + (i / (btScalar)(n - 1)) * 0.9;
		psb->setTotalMass(20);
		pdemo->getSoftDynamicsWorld()->addSoftBody(psb);

	}
}

//
// Rope attach
//
static void Init_RopeAttach(SoftDemo *pdemo)
{
	//TRACEDEMO
	pdemo->m_softBodyWorldInfo.m_sparsesdf.RemoveReferences(0);
	struct  Functors
	{
		static btSoftBody *CtorRope(SoftDemo *pdemo, const btVector3 &p)
		{
			btSoftBody *psb = btSoftBodyHelpers::CreateRope(pdemo->m_softBodyWorldInfo, p, p + btVector3(10, 0, 0), 8, 1);
			psb->setTotalMass(50);
			pdemo->getSoftDynamicsWorld()->addSoftBody(psb);
			return (psb);
		}
	};
	btTransform startTransform;
	startTransform.setIdentity();
	startTransform.setOrigin(btVector3(12, 8, 0));
	btRigidBody        *body = pdemo->localCreateRigidBody(50, startTransform, new btBoxShape(btVector3(2, 6, 2)));
	btSoftBody *psb0 = Functors::CtorRope(pdemo, btVector3(0, 8, -1));
	btSoftBody *psb1 = Functors::CtorRope(pdemo, btVector3(0, 8, +1));
	psb0->appendAnchor(psb0->m_nodes.size() - 1, body);
	psb1->appendAnchor(psb1->m_nodes.size() - 1, body);
}

//
// Cloth attach
//
/*
static void Init_ClothAttach(SoftDemo *pdemo)
{
//TRACEDEMO
const btScalar  s = 4;
const btScalar  h = 6;
const int       r = 9;
btSoftBody     *psb = btSoftBodyHelpers::CreatePatch(pdemo->m_softBodyWorldInfo,
btVector3(-s, h - 8, +s), // (-s,h,-s)
btVector3(+s, h - 8, +s), // (+s,h,-s)
btVector3(-s, h, +s),
btVector3(+s, h, +s), r, r, 4 + 8, true); // 1+2+4+8

btSoftBody::Material*   pm=psb->appendMaterial();
pm->m_kLST=1;
pm->m_flags             -=  btSoftBody::fMaterial::DebugDraw;
psb->generateBendingConstraints(2,pm);
psb->m_cfg.kDG          =   0.05; // 3
psb->m_cfg.kLF          =   0.05; // 4
psb->m_cfg.kDF          =   (btScalar)1; // 0.2
psb->m_cfg.piterations  =   25; // 1
psb->m_cfg.diterations  =   16; // 0

pdemo->getSoftDynamicsWorld()->addSoftBody(psb);

btTransform startTransform;
startTransform.setIdentity();
startTransform.setOrigin(btVector3(0, h, -(s + 3.5)));
//btRigidBody*      body=pdemo->localCreateRigidBody(20,startTransform,new btBoxShape(btVector3(s,1,3)));
//psb->appendAnchor(0,body);
//psb->appendAnchor(r-1,body);
pdemo->m_cutting = true;
}
*/

static void Init_ClothAttach_fold(SoftDemo *pdemo)  //創造出一個像布一樣的柔性物體
{
	//TRACEDEMO

	const btScalar  s = 6;			//邊長
	const btScalar  h = 6;			//高度
	const int       r = 3;			//把整個布分成2^r=8塊
	
	//初始一張紙的狀態
	cur_step = 0;
	state = 0;
	current_mode = 0;
	record = -1;
	fold_condition = 0;

	teach_mode = false;

	btSoftBody     *psb = btSoftBodyHelpers::CreatePatch(pdemo->m_softBodyWorldInfo, btVector3(-s, h - 12, +s), //高度從-6~+6
		btVector3(+s, h - 12, +s),
		btVector3(-s, h, +s),
		btVector3(+s, h, +s), r, r, 0, true); //4+8
	psb->m_nodes[4].m_im = 0;
	old_position.resize(200, btVector3(0, 0, 0)); //保存原來點以供展開 最多存到200個點
	for (int j = 0; j < r * r; j++)
		old_position[j] = psb->m_nodes[j].m_x;
	//material changing
	btSoftBody::Material   *pm = psb->appendMaterial();
	pm->m_kLST = 1;								//控制材質硬度，調到1為最接近硬體
	pm->m_flags             -=  btSoftBody::fMaterial::DebugDraw;
	//psb->generateBendingConstraints(2,pm);		//建立端點限制

	psb->m_cfg.aeromodel    =   btSoftBody::eAeroModel::V_TwoSided;
	//psb->m_cfg.kVCF           =   1;//1
	psb->m_cfg.kDG          =   4;//3			//調高會增加無形的拉力阻力
	psb->m_cfg.kLF          =   6;//4			//調高更容易飄動
	psb->m_cfg.kDP          =   0;//0			//空氣阻力
	//psb->m_cfg.kPR            =   0;//0		//空氣壓力
	//psb->m_cfg.kVC            =   0;//0		//當setPose(true,...)被呼叫時，調高後會有浮動在空中的情形
	psb->m_cfg.kDF          =   (btScalar)1;//0.2	//當KDF=0代表會滑動，KDF=1為會固定住
	psb->m_cfg.kMT          =   0;//0			//當setPose(...,true)被呼叫時，定義pose matching的係數
	psb->m_cfg.kCHR         =   (btScalar)1.0;//1	//定義軟體和硬體接觸的情形，0為硬體會溶入軟體，1為硬體完全被彈開
	psb->m_cfg.kKHR         =   (btScalar)0.1;//0.1	//定義軟體和靜止物體接觸的情形
	psb->m_cfg.kSHR         =   (btScalar)1; //1		//定義軟體和其他軟體接觸的情形
	psb->m_cfg.kAHR         =   (btScalar)0.7;//0.7
	psb->m_cfg.kSRHR_CL     =   (btScalar)0.1;//0.1
	psb->m_cfg.kSKHR_CL     =   (btScalar)1;//1
	psb->m_cfg.kSSHR_CL     =   (btScalar)0.5;//0.5
	psb->m_cfg.kSR_SPLT_CL  =   (btScalar)0.5;//0.5
	psb->m_cfg.kSK_SPLT_CL  =   (btScalar)0.5;//0.5
	psb->m_cfg.kSS_SPLT_CL  =   (btScalar)0.5;//0.5
	psb->m_cfg.maxvolume    =   (btScalar)1;//1
	psb->m_cfg.timescale    =   0.70;//1			//將物體速度變慢
	psb->m_cfg.viterations  =   0;//0
	psb->m_cfg.piterations  =   25;//1			//讓紙位置固定
	psb->m_cfg.diterations  =   16;//0
	psb->m_cfg.citerations  =   4;//4
	psb->setTotalMass(5);
	oldmass = psb->m_nodes[0].m_im;
	pdemo->getSoftDynamicsWorld()->addSoftBody(psb);
	btTransform startTransform;
	startTransform.setIdentity();
	startTransform.setOrigin(btVector3(0, h, -(s + 3.5)));
	pdemo->m_cutting = true;


}

static void Init_ClothAttach(SoftDemo *pdemo)
{
	//TRACEDEMO
	const btScalar  s = 6;
	const btScalar  h = 6;
	const int       r = 2;
	btSoftBody     *psb = btSoftBodyHelpers::CreatePatch(pdemo->m_softBodyWorldInfo, btVector3(-s, h - 12, +s),
		btVector3(+s, h - 12, +s),
		btVector3(-s, h, +s),
		btVector3(+s, h, +s), r, r, 4 + 8, true); //4+8
	old_position.resize(r * r, btVector3(0, 0, 0)); //�O�s�����I�H�Ѯi�}
	for (int j = 0; j < r * r; j++)
		old_position[j] = psb->m_nodes[j].m_x;
	//material changing
	btSoftBody::Material   *pm = psb->appendMaterial();
	pm->m_kLST = 1;
	pm->m_flags             -=  btSoftBody::fMaterial::DebugDraw;
	//psb->generateBendingConstraints(2,pm);
	psb->m_cfg.aeromodel    =   btSoftBody::eAeroModel::V_TwoSided;
	//psb->m_cfg.kVCF           =   1;//1
	psb->m_cfg.kDG          =   2;//3
	psb->m_cfg.kLF          =   2;//4
	psb->m_cfg.kDP          =   0;//0
	//psb->m_cfg.kPR            =   0;//0
	//psb->m_cfg.kVC            =   0;//0
	psb->m_cfg.kDF          =   (btScalar)1;//0.2
	psb->m_cfg.kMT          =   0;//0
	psb->m_cfg.kCHR         =   (btScalar)1.0;//1
	psb->m_cfg.kKHR         =   (btScalar)0.1;//0.1
	psb->m_cfg.kSHR         =   (btScalar)1; //1
	psb->m_cfg.kAHR         =   (btScalar)0.7;//0.7
	psb->m_cfg.kSRHR_CL     =   (btScalar)0.1;//0.1
	psb->m_cfg.kSKHR_CL     =   (btScalar)1;//1
	psb->m_cfg.kSSHR_CL     =   (btScalar)0.5;//0.5
	psb->m_cfg.kSR_SPLT_CL  =   (btScalar)0.5;//0.5
	psb->m_cfg.kSK_SPLT_CL  =   (btScalar)0.5;//0.5
	psb->m_cfg.kSS_SPLT_CL  =   (btScalar)0.5;//0.5
	psb->m_cfg.maxvolume    =   (btScalar)1;//1
	psb->m_cfg.timescale    =   0.3;//1
	psb->m_cfg.viterations  =   0;//0
	psb->m_cfg.piterations  =   25;//1
	psb->m_cfg.diterations  =   16;//0
	psb->m_cfg.citerations  =   4;//4
	psb->setTotalMass(5);
	//psb->setPose(true,false);//0 1
	pdemo->getSoftDynamicsWorld()->addSoftBody(psb);
	btTransform startTransform;
	startTransform.setIdentity();
	startTransform.setOrigin(btVector3(0, h, -(s + 3.5)));
	pdemo->m_cutting = true;
}

//
// Impact
//
static void Init_Impact(SoftDemo *pdemo)
{
	//TRACEDEMO
	btSoftBody *psb = btSoftBodyHelpers::CreateRope(pdemo->m_softBodyWorldInfo,   btVector3(0, 0, 0),
		btVector3(0, -1, 0),
		0,
		1);
	pdemo->getSoftDynamicsWorld()->addSoftBody(psb);
	psb->m_cfg.kCHR = 0.5;
	btTransform startTransform;
	startTransform.setIdentity();
	startTransform.setOrigin(btVector3(0, 20, 0));
	pdemo->localCreateRigidBody(10, startTransform, new btBoxShape(btVector3(2, 2, 2)));
}

static void Init_CapsuleCollision(SoftDemo *pdemo)
{
#ifdef USE_AMD_OPENCL
	btAlignedObjectArray<btSoftBody *> emptyArray;
	if (g_openCLSIMDSolver)
		g_openCLSIMDSolver->optimize(emptyArray);
#endif //USE_AMD_OPENCL

	//TRACEDEMO
	const btScalar  s = 4;
	const btScalar  h = 6;
	const int       r = 20;

	btTransform startTransform;
	startTransform.setIdentity();
	startTransform.setOrigin(btVector3(0, h - 2, 0));

	btCollisionShape *capsuleShape = new btCapsuleShapeX(1, 5);
	capsuleShape->setMargin( 0.5 );

	//  capsule->setLocalScaling(btVector3(5,1,1));
	//  btRigidBody*        body=pdemo->localCreateRigidBody(20,startTransform,capsuleShape);
	btRigidBody        *body = pdemo->localCreateRigidBody(0, startTransform, capsuleShape);
	body->setFriction( 0.8f );

	int fixed = 0; //4+8;
	btSoftBody     *psb = btSoftBodyHelpers::CreatePatch(pdemo->m_softBodyWorldInfo, btVector3(-s, h, -s),
		btVector3(+s, h, -s),
		btVector3(-s, h, +s),
		btVector3(+s, h, +s), r, r, fixed, true);
	pdemo->getSoftDynamicsWorld()->addSoftBody(psb);
	psb->setTotalMass(0.1);

	psb->m_cfg.piterations = 10;
	psb->m_cfg.citerations = 10;
	psb->m_cfg.diterations = 10;
	//  psb->m_cfg.viterations = 10;


	//  psb->appendAnchor(0,body);
	//  psb->appendAnchor(r-1,body);
	//  pdemo->m_cutting=true;
}

//
// Collide3
//
static void Init_Collide3(SoftDemo *pdemo)
{
	//TRACEDEMO
	{
		const btScalar  s = 8;
		btSoftBody     *psb = btSoftBodyHelpers::CreatePatch( pdemo->m_softBodyWorldInfo, btVector3(-s, 0, -s),
			btVector3(+s, 0, -s),
			btVector3(-s, 0, +s),
			btVector3(+s, 0, +s),
			15, 15, 1 + 2 + 4 + 8, true);
		psb->m_materials[0]->m_kLST =   0.4;
		psb->m_cfg.collisions       |=  btSoftBody::fCollision::VF_SS;
		psb->setTotalMass(150);
		pdemo->getSoftDynamicsWorld()->addSoftBody(psb);
	}
	{
		const btScalar  s = 4;
		const btVector3 o = btVector3(5, 10, 0);
		btSoftBody     *psb = btSoftBodyHelpers::CreatePatch( pdemo->m_softBodyWorldInfo,
			btVector3(-s, 0, -s) + o,
			btVector3(+s, 0, -s) + o,
			btVector3(-s, 0, +s) + o,
			btVector3(+s, 0, +s) + o,
			7, 7, 0, true);
		btSoftBody::Material   *pm = psb->appendMaterial();
		pm->m_kLST              =   0.1;
		pm->m_flags             -=  btSoftBody::fMaterial::DebugDraw;
		psb->generateBendingConstraints(2, pm);
		psb->m_materials[0]->m_kLST =   0.5;
		psb->m_cfg.collisions       |=  btSoftBody::fCollision::VF_SS;
		psb->setTotalMass(150);
		pdemo->getSoftDynamicsWorld()->addSoftBody(psb);
		pdemo->m_cutting = true;
	}
}

//
// Aerodynamic forces, 50x1g flyers
//
static void Init_Aero(SoftDemo *pdemo)
{
	//TRACEDEMO
	const btScalar  s = 2;
	const btScalar  h = 10;
	const int       segments = 6;
	const int       count = 50;
	for (int i = 0; i < count; ++i)
	{
		btSoftBody     *psb = btSoftBodyHelpers::CreatePatch(pdemo->m_softBodyWorldInfo, btVector3(-s, h, -s),
			btVector3(+s, h, -s),
			btVector3(-s, h, +s),
			btVector3(+s, h, +s),
			segments, segments,
			0, true);
		btSoftBody::Material   *pm = psb->appendMaterial();
		pm->m_flags             -=  btSoftBody::fMaterial::DebugDraw;
		psb->generateBendingConstraints(2, pm);
		psb->m_cfg.kLF          =   0.004;
		psb->m_cfg.kDG          =   0.0003;
		psb->m_cfg.aeromodel    =   btSoftBody::eAeroModel::V_TwoSided;
		btTransform     trs;
		btQuaternion    rot;
		btVector3       ra = Vector3Rand() * 0.1;
		btVector3       rp = Vector3Rand() * 15 + btVector3(0, 20, 80);
		rot.setEuler(SIMD_PI / 8 + ra.x(), -SIMD_PI / 7 + ra.y(), ra.z());
		trs.setIdentity();
		trs.setOrigin(rp);
		trs.setRotation(rot);
		psb->transform(trs);
		psb->setTotalMass(0.1);
		psb->addForce(btVector3(0, 2, 0), 0);
		pdemo->getSoftDynamicsWorld()->addSoftBody(psb);

	}
	pdemo->m_autocam = true;
}

static void Init_Aero2(SoftDemo *pdemo)
{
	//TRACEDEMO
	const btScalar  s = 5;
	const int       segments = 10;
	const int       count = 5;
	btVector3 pos(-s * segments, 0, 0);
	btScalar gap = 0.5;

	for (int i = 0; i < count; ++i)
	{
		btSoftBody     *psb = btSoftBodyHelpers::CreatePatch( pdemo->m_softBodyWorldInfo, btVector3(-s, 0, -s * 3),
			btVector3(+s, 0, -s * 3),
			btVector3(-s, 0, +s),
			btVector3(+s, 0, +s),
			segments, segments * 3,
			1 + 2, true);

		psb->getCollisionShape()->setMargin(0.5);
		btSoftBody::Material *pm = psb->appendMaterial();
		pm->m_kLST      =   0.0004;
		pm->m_flags     -=  btSoftBody::fMaterial::DebugDraw;
		psb->generateBendingConstraints(2, pm);

		psb->m_cfg.kLF          =   0.05;
		psb->m_cfg.kDG          =   0.01;

		//psb->m_cfg.kLF            =   0.004;
		//psb->m_cfg.kDG            =   0.0003;

		psb->m_cfg.piterations = 2;
		psb->m_cfg.aeromodel    =   btSoftBody::eAeroModel::V_TwoSidedLiftDrag;


		psb->setWindVelocity(btVector3(4, -12.0, -25.0));

		btTransform     trs;
		btQuaternion    rot;
		pos += btVector3(s * 2 + gap, 0, 0);
		rot.setRotation(btVector3(1, 0, 0), btScalar(SIMD_PI / 2));
		trs.setIdentity();
		trs.setOrigin(pos);
		trs.setRotation(rot);
		psb->transform(trs);
		psb->setTotalMass(2.0);

		pdemo->getSoftDynamicsWorld()->addSoftBody(psb);
	}

	pdemo->m_autocam = true;
}

//
// Friction
//
static void Init_Friction(SoftDemo *pdemo)
{
	//TRACEDEMO
	const btScalar  bs = 2;
	const btScalar  ts = bs + bs / 4;
	for (int i = 0, ni = 20; i < ni; ++i)
	{
		const btVector3 p(-ni * ts / 2 + i * ts, -10 + bs, 40);
		btSoftBody     *psb = Ctor_SoftBox(pdemo, p, btVector3(bs, bs, bs));
		psb->m_cfg.kDF  =   0.1 * ((i + 1) / (btScalar)ni);
		psb->addVelocity(btVector3(0, 0, -10));
	}
}

//
// Pressure
//
static void Init_Pressure(SoftDemo *pdemo)
{
	//TRACEDEMO
	btSoftBody *psb = btSoftBodyHelpers::CreateEllipsoid(pdemo->m_softBodyWorldInfo, btVector3(35, 25, 0),
		btVector3(1, 1, 1) * 3,
		512);
	psb->m_materials[0]->m_kLST =   0.1;
	psb->m_cfg.kDF              =   1;
	psb->m_cfg.kDP              =   0.001; // fun factor...
	psb->m_cfg.kPR              =   2500;
	psb->setTotalMass(30, true);
	pdemo->getSoftDynamicsWorld()->addSoftBody(psb);

	Ctor_BigPlate(pdemo);
	Ctor_LinearStair(pdemo, btVector3(0, 0, 0), btVector3(2, 1, 5), 0, 10);
	pdemo->m_autocam = true;

}

//
// Volume conservation
//
static void Init_Volume(SoftDemo *pdemo)
{
	//TRACEDEMO
	btSoftBody *psb = btSoftBodyHelpers::CreateEllipsoid(pdemo->m_softBodyWorldInfo, btVector3(35, 25, 0),
		btVector3(1, 1, 1) * 3,
		512);
	psb->m_materials[0]->m_kLST =   0.45;
	psb->m_cfg.kVC              =   20;
	psb->setTotalMass(50, true);
	psb->setPose(true, false);
	pdemo->getSoftDynamicsWorld()->addSoftBody(psb);

	Ctor_BigPlate(pdemo);
	Ctor_LinearStair(pdemo, btVector3(0, 0, 0), btVector3(2, 1, 5), 0, 10);
	pdemo->m_autocam = true;

}

//
// Stick+Bending+Rb's
//
static void Init_Sticks(SoftDemo *pdemo)
{
	//TRACEDEMO
	const int       n = 16;
	const int       sg = 4;
	const btScalar  sz = 5;
	const btScalar  hg = 4;
	const btScalar  in = 1 / (btScalar)(n - 1);
	for (int y = 0; y < n; ++y)
	{
		for (int x = 0; x < n; ++x)
		{
			const btVector3 org(-sz + sz * 2 * x * in,
				-10,
				-sz + sz * 2 * y * in);
			btSoftBody     *psb = btSoftBodyHelpers::CreateRope(  pdemo->m_softBodyWorldInfo, org,
				org + btVector3(hg * 0.001, hg, 0),
				sg,
				1);
			psb->m_cfg.kDP      =   0.005;
			psb->m_cfg.kCHR     =   0.1;
			for (int i = 0; i < 3; ++i)
			{
				psb->generateBendingConstraints(2 + i);
			}
			psb->setMass(1, 0);
			psb->setTotalMass(0.01);
			pdemo->getSoftDynamicsWorld()->addSoftBody(psb);

		}
	}
	Ctor_BigBall(pdemo);
}

//
// Bending
//
static void Init_Bending(SoftDemo *pdemo)
{
	//TRACEDEMO
	const btScalar  s = 4;
	const btVector3 x[] = {   btVector3(-s, 0, -s),
		btVector3(+s, 0, -s),
		btVector3(+s, 0, +s),
		btVector3(-s, 0, +s)
	};
	const btScalar  m[] = {   0, 0, 0, 1};
	btSoftBody     *psb = new btSoftBody(&pdemo->m_softBodyWorldInfo, 4, x, m);
	psb->appendLink(0, 1);
	psb->appendLink(1, 2);
	psb->appendLink(2, 3);
	psb->appendLink(3, 0);
	psb->appendLink(0, 2);

	pdemo->getSoftDynamicsWorld()->addSoftBody(psb);
}

//
// 100kg cloth locked at corners, 10 falling 10kg rb's.
//
static void Init_Cloth(SoftDemo *pdemo)
{
	//TRACEDEMO
	const btScalar  s = 8;
	btSoftBody     *psb = btSoftBodyHelpers::CreatePatch( pdemo->m_softBodyWorldInfo, btVector3(-s, 0, -s),
		btVector3(+s, 0, -s),
		btVector3(-s, 0, +s),
		btVector3(+s, 0, +s),
		31, 31,
		//      31,31,
		1 + 2 + 4 + 8, true);

	psb->getCollisionShape()->setMargin(0.5);
	btSoftBody::Material *pm = psb->appendMaterial();
	pm->m_kLST      =   0.4;
	pm->m_flags     -=  btSoftBody::fMaterial::DebugDraw;
	psb->generateBendingConstraints(2, pm);
	psb->setTotalMass(150);
	pdemo->getSoftDynamicsWorld()->addSoftBody(psb);

	Ctor_RbUpStack(pdemo, 10);
	pdemo->m_cutting = true;
}

//
// Cutting1
//
static void Init_Cutting1(SoftDemo *pdemo)
{
	const btScalar  s = 6;
	const btScalar  h = 2;
	const int       r = 16;
	const btVector3 p[] = {   btVector3(+s, h, -s),
		btVector3(-s, h, -s),
		btVector3(+s, h, +s),
		btVector3(-s, h, +s)
	};
	btSoftBody *psb = btSoftBodyHelpers::CreatePatch(pdemo->m_softBodyWorldInfo, p[0], p[1], p[2], p[3], r, r, 1 + 2 + 4 + 8, true);
	pdemo->getSoftDynamicsWorld()->addSoftBody(psb);
	psb->m_cfg.piterations = 1;
	pdemo->m_cutting = true;
}

//
// Clusters
//

//
static void         Ctor_Gear(SoftDemo *pdemo, const btVector3 &pos, btScalar speed)
{
	btTransform startTransform;
	startTransform.setIdentity();
	startTransform.setOrigin(pos);
	btCompoundShape    *shape = new btCompoundShape();
#if 1
	shape->addChildShape(btTransform(btQuaternion(0, 0, 0)), new btBoxShape(btVector3(5, 1, 6)));
	shape->addChildShape(btTransform(btQuaternion(0, 0, SIMD_HALF_PI)), new btBoxShape(btVector3(5, 1, 6)));
#else
	shape->addChildShape(btTransform(btQuaternion(0, 0, 0)), new btCylinderShapeZ(btVector3(5, 1, 7)));
	shape->addChildShape(btTransform(btQuaternion(0, 0, SIMD_HALF_PI)), new btBoxShape(btVector3(4, 1, 8)));
#endif
	btRigidBody        *body = pdemo->localCreateRigidBody(10, startTransform, shape);
	body->setFriction(1);
	btDynamicsWorld    *world = pdemo->getDynamicsWorld();
	btHingeConstraint  *hinge = new btHingeConstraint(*body, btTransform::getIdentity());
	if (speed != 0) hinge->enableAngularMotor(true, speed, 3);
	world->addConstraint(hinge);
}

//
static struct MotorControl : btSoftBody::AJoint::IControl
{
	MotorControl()
	{
		goal = 0;
		maxtorque = 0;
	}
	btScalar    Speed(btSoftBody::AJoint *, btScalar current)
	{
		return (current + btMin(maxtorque, btMax(-maxtorque, goal - current)));
	}
	btScalar    goal;
	btScalar    maxtorque;
}   motorcontrol;

//
struct SteerControl : btSoftBody::AJoint::IControl
{
	SteerControl(btScalar s)
	{
		angle = 0;
		sign = s;
	}
	void        Prepare(btSoftBody::AJoint *joint)
	{
		joint->m_refs[0][0] = btCos(angle * sign);
		joint->m_refs[0][2] = btSin(angle * sign);
	}
	btScalar    Speed(btSoftBody::AJoint *joint, btScalar current)
	{
		return (motorcontrol.Speed(joint, current));
	}
	btScalar    angle;
	btScalar    sign;
};

static SteerControl steercontrol_f(+1);
static SteerControl steercontrol_r(-1);

static void Init_ClusterCollide1(SoftDemo *pdemo)
{
	const btScalar  s = 8;
	btSoftBody     *psb = btSoftBodyHelpers::CreatePatch( pdemo->m_softBodyWorldInfo, btVector3(-s, 0, -s),
		btVector3(+s, 0, -s),
		btVector3(-s, 0, +s),
		btVector3(+s, 0, +s),
		17, 17, //9,9,//31,31,
		1 + 2 + 4 + 8,
		true);
	btSoftBody::Material *pm = psb->appendMaterial();
	pm->m_kLST      =   0.4;
	pm->m_flags     -=  btSoftBody::fMaterial::DebugDraw;
	psb->m_cfg.kDF  =   1;
	psb->m_cfg.kSRHR_CL     =   1;
	psb->m_cfg.kSR_SPLT_CL  =   0;
	psb->m_cfg.collisions   =   btSoftBody::fCollision::CL_SS +

		btSoftBody::fCollision::CL_RS;
	psb->generateBendingConstraints(2, pm);

	psb->getCollisionShape()->setMargin(0.05);
	psb->setTotalMass(50);

	///pass zero in generateClusters to create  cluster for each tetrahedron or triangle
	psb->generateClusters(0);
	//psb->generateClusters(64);

	pdemo->getSoftDynamicsWorld()->addSoftBody(psb);

	Ctor_RbUpStack(pdemo, 10);
}





/* Init     */
void (*demofncs[])(SoftDemo *) =
{
	Init_ClothAttach_fold,
	Init_ClothAttach

};

void    SoftDemo::clientResetScene()
{
	start_node = -1;
	m_azi = 0;
	m_cameraDistance = 3.f; // 30.f
	m_cameraTargetPosition.setValue(0, 0, 0);

	DemoApplication::clientResetScene();
	/* Clean up */
	for (int i = m_dynamicsWorld->getNumCollisionObjects() - 1; i >= 0; i--)
	{
		btCollisionObject  *obj = m_dynamicsWorld->getCollisionObjectArray()[i];
		btRigidBody        *body = btRigidBody::upcast(obj);
		if (body && body->getMotionState())
		{
			delete body->getMotionState();
		}
		while (m_dynamicsWorld->getNumConstraints())
		{
			btTypedConstraint  *pc = m_dynamicsWorld->getConstraint(0);
			m_dynamicsWorld->removeConstraint(pc);
			delete pc;
		}
		btSoftBody *softBody = btSoftBody::upcast(obj);
		if (softBody)
		{
			getSoftDynamicsWorld()->removeSoftBody(softBody);
		}
		else
		{
			btRigidBody *body = btRigidBody::upcast(obj);
			if (body)
				m_dynamicsWorld->removeRigidBody(body);
			else
				m_dynamicsWorld->removeCollisionObject(obj);
		}
		delete obj;
	}


	//create ground object
	btTransform tr;
	tr.setIdentity();
	tr.setOrigin(btVector3(0, -12, 0));

	btCollisionObject *newOb = new btCollisionObject();
	newOb->setWorldTransform(tr);
	newOb->setInterpolationWorldTransform( tr);
	int lastDemo = (sizeof(demofncs) / sizeof(demofncs[0])) - 1;

	if (current_demo < 0)
		current_demo = lastDemo;
	if (current_demo > lastDemo)
		current_demo = 0;


	if (current_demo > 19)
	{
		newOb->setCollisionShape(m_collisionShapes[0]);
	}
	else
	{
		newOb->setCollisionShape(m_collisionShapes[1]);
	}

	m_dynamicsWorld->addCollisionObject(newOb);

	m_softBodyWorldInfo.m_sparsesdf.Reset();







	motorcontrol.goal = 0;
	motorcontrol.maxtorque = 0;



	m_softBodyWorldInfo.air_density     =   (btScalar)1.2;
	m_softBodyWorldInfo.water_density   =   0;
	m_softBodyWorldInfo.water_offset        =   0;
	m_softBodyWorldInfo.water_normal        =   btVector3(0, 0, 0);
	m_softBodyWorldInfo.m_gravity.setValue(0, -10, 0);


	m_autocam                       =   false;
	m_raycast                       =   false;
	m_cutting                       =   false;
	m_results.fraction              =   1.f;
	demofncs[current_demo](this);
}



void SoftDemo::clientMoveAndDisplay()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);




	float ms = getDeltaTimeMicroseconds();
	float dt = ms / 1000000.f;//1.0/60.;



	if (m_dynamicsWorld)
	{

		if (sDemoMode)
		{
			static float demoCounter = DEMO_MODE_TIMEOUT;
			demoCounter -= dt;
			if (demoCounter < 0)
			{

				demoCounter = DEMO_MODE_TIMEOUT;
				current_demo++;
				current_demo = current_demo % (sizeof(demofncs) / sizeof(demofncs[0]));
				clientResetScene();
			}
		}


		//#define FIXED_STEP
#ifdef FIXED_STEP
		m_dynamicsWorld->stepSimulation(dt = 1.0f / 60.f, 0);

#else
		//during idle mode, just run 1 simulation step maximum, otherwise 4 at max
		//  int maxSimSubSteps = m_idle ? 1 : 4;
		//if (m_idle)
		//  dt = 1.0/420.f;

		int numSimSteps;
		numSimSteps = m_dynamicsWorld->stepSimulation(dt);
		//numSimSteps = m_dynamicsWorld->stepSimulation(dt,10,1./240.f);

#ifdef VERBOSE_TIMESTEPPING_CONSOLEOUTPUT
		if (!numSimSteps)
			printf("Interpolated transforms\n");
		else
		{
			if (numSimSteps > maxSimSubSteps)
			{
				//detect dropping frames
				printf("Dropped (%i) simulation steps out of %i\n", numSimSteps - maxSimSubSteps, numSimSteps);
			}
			else
			{
				printf("Simulated (%i) steps\n", numSimSteps);
			}
		}
#endif //VERBOSE_TIMESTEPPING_CONSOLEOUTPUT

#endif

#ifdef USE_AMD_OPENCL
		if (g_openCLSIMDSolver)
			g_openCLSIMDSolver->copyBackToSoftBodies();
#endif //USE_AMD_OPENCL

		if (m_drag)
		{
			m_node->m_v *= 0;
		}

		m_softBodyWorldInfo.m_sparsesdf.GarbageCollect();

		//optional but useful: debug drawing

	}

#ifdef USE_QUICKPROF
	btProfiler::beginBlock("render");
#endif //USE_QUICKPROF 

	renderme();


	//render the graphics objects, with center of mass shift

	updateCamera();



#ifdef USE_QUICKPROF
	btProfiler::endBlock("render");
#endif
	glFlush();
	//some additional debugging info
#ifdef PRINT_CONTACT_STATISTICS
	printf("num manifolds: %i\n", gNumManifold);
	printf("num gOverlappingPairs: %i\n", gOverlappingPairs);

#endif //PRINT_CONTACT_STATISTICS


	swapBuffers();

}

btVector3 newPivot( 0.f, 0.f, 0.f );
btVector3 *point_set_v[3];
btVector3 *point_set_h[3];
btVector3 star_point[7]		= {btVector3( -0.116006,  5.139614, 6.000000 ),
	btVector3( -1.026965,  2.089578, 5.999999 ),
	btVector3( -4.581097,  1.367989, 6.000000 ),
	btVector3( -1.454528, -0.550541, 6.000000 ),
	btVector3( -3.261538, -4.383197, 6.000000 ),
	btVector3( -1.132061, -3.289138, 6.000000 ),
	btVector3( -0.116006, -2.696469, 6.000000 )
};
btVector3 heart_point[7]    = {btVector3( -0.069894, 3.017673, 5.999999 ),
	btVector3( -1.814291,  3.850095, 5.999999 ),
	btVector3( -4.049088,  3.619138, 5.999999 ),
	btVector3( -4.745137,  1.321365, 6.000000 ),
	btVector3( -2.910381, -0.785488, 6.000000 ),
	btVector3( -1.505490, -1.939827, 5.999999 ),
	btVector3( -0.060064, -3.012280, 6.000001 )
};
btVector3 triangle_point[7] = {btVector3( -0.047926,  3.490256, 5.999999 ),
	btVector3( -2.060899,  3.466443, 5.999966 ),
	btVector3( -4.720665,  3.490234, 5.999941 ),
	btVector3( -3.418616,  1.099191, 5.999981 ),
	btVector3( -2.248087, -0.999268, 6.000000 ),
	btVector3( -0.946964, -3.310328, 6.000000 ),
	btVector3( -0.047159, -4.903529, 6.000000 )
};
btVector3 pacman_point[7]	= {btVector3( -4.404294, 0.176437, 6.000000 ),
	btVector3( -3.849380, 2.321916, 5.999998 ),
	btVector3( -2.419598, 3.734642, 5.999998 ),
	btVector3( -0.092964, 4.265294, 6.000000 ),
	btVector3(  2.419381, 3.780828, 5.999997 ),
	btVector3(  3.759440, 1.857016, 5.999939 ),
	btVector3(  1.874327, 0.170371, 5.999947 )
};
btVector3 cross_point[7]	= {btVector3( -1.921633, 0.100055, 5.999999 ),
	btVector3( -4.162062, 4.034711, 6.000000 ),
	btVector3( -2.274363, 5.001762, 6.000000 ),
	btVector3( -0.233591, 1.670837, 6.000000 ),
	btVector3(  1.995331, 5.139614, 6.000000 ),
	btVector3(  4.069424, 3.988571, 6.000004 ),
	btVector3(  1.522840, 0.120055, 6.000057 )
};
btVector3 square_point[7]	= {btVector3( -4.147284, 0.106176, 6.000268 ),
	btVector3( -4.084215, 2.136092, 6.000130 ),
	btVector3( -4.068505, 4.103907, 6.000000 ),
	btVector3(  0.092977, 4.196140, 6.000000 ),
	btVector3(  3.765425, 4.219193, 6.000000 ),
	btVector3(  3.780316, 2.205812, 6.000159 ),
	btVector3(  3.842906, 0.106176, 6.000381 )
};
int cur_teaching_mode_v = 0;
int cur_teaching_mode_h = 0;

int cur_fold = 1;
int cur_point = 0;
bool isMode0    = false;
bool isMode1    = false;
bool isMode2Or5 = false;
bool isMode3    = false;
bool isMode4    = false;
bool case_s = false;

int timeCount1 = 0;
int timeCount2 = 0;
int timeCount3 = 0;
int timeCount4 = 0;
int timeCount5 = 0;
int timeCount6 = 0;

void    SoftDemo::renderme()
{

	/* ***** ***** ***** ***** ***** Draw scene begin ***** ***** ***** ***** ***** */
	glPushMatrix();

	ApplyLeftFrustum( 2 );
	glColorMask( true, false, false, true );
	gluLookAt( cp[0], cp[1], cp[2], ctp[0], ctp[1], ctp[2], cu[0], cu[1], cu[2] );
	DrawScene();


	glClear( GL_DEPTH_BUFFER_BIT );
	//glEnable( GL_BLEND );
	//glBlendFunc( GL_ONE,GL_ONE );

	ApplyRightFrustum( 2 );
	glColorMask( false, true, true, true );
	gluLookAt( cp[0], cp[1], cp[2], ctp[0], ctp[1], ctp[2], cu[0], cu[1], cu[2] );
	DrawScene();

	//glDisable( GL_BLEND );
	glColorMask( true, true, true, true );

	glPopMatrix();

	/* ***** ***** ***** ***** ***** Draw scene  end  ***** ***** ***** ***** ***** */

	btIDebugDraw   *idraw = m_dynamicsWorld->getDebugDrawer();

	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	m_dynamicsWorld->debugDrawWorld();

	//int debugMode = m_dynamicsWorld->getDebugDrawer()? m_dynamicsWorld->getDebugDrawer()->getDebugMode() : -1;

	btSoftRigidDynamicsWorld *softWorld = (btSoftRigidDynamicsWorld *)m_dynamicsWorld;
	//btIDebugDraw* sdraw = softWorld ->getDebugDrawer();

	for (  int i = 0; i < softWorld->getSoftBodyArray().size(); i++)
	{
		btSoftBody *psb = (btSoftBody *)softWorld->getSoftBodyArray()[i];
		if (softWorld->getDebugDrawer() && !(softWorld->getDebugDrawer()->getDebugMode() & (btIDebugDraw::DBG_DrawWireframe)))
		{
			btSoftBodyHelpers::DrawFrame(psb, softWorld->getDebugDrawer());
			btSoftBodyHelpers::Draw(psb, softWorld->getDebugDrawer(), softWorld->getDrawFlags());
		}
	}

	/* Bodies       */
	btVector3   ps(0, 0, 0);
	int         nps = 0;
	btSoftBodyArray    &sbs = getSoftDynamicsWorld()->getSoftBodyArray();
	for (int ib = 0; ib < sbs.size(); ++ib)
	{
		btSoftBody *psb = sbs[ib];
		nps += psb->m_nodes.size();


		/********************************1摺******************************/

		if (record != -1 && psb->m_nodes[record].m_im == 0) //讓mass=0的點除了剪取之外不能按下
		{
			if (cur_step != 4 && cur_step != 5)
			{
				m_results.fraction = 1.f;
				m_drag = false;
			}
		}

		if (cur_step == 1)
		{
			fold_first();
		}

		if (cutend == 1)		//已剪斷
		{
			m_dynamicsWorld->setGravity(btVector3(0, -10, 0));
			m_softBodyWorldInfo.m_gravity.setValue(0, -10, 0);
			cutend = 0;
		}
		/********************************2摺******************************/
		else if (cur_step == 2)
		{
			fold_second();
		}


		for (int i = 0; i < psb->m_nodes.size(); ++i)
		{
			ps += psb->m_nodes[i].m_x;
		}
	}

	ps /= nps;
	if (m_autocam)
		m_cameraTargetPosition += (ps - m_cameraTargetPosition) * 0.05;
	/* Anm          */
	if (!isIdle())
		m_animtime = m_clock.getTimeMilliseconds() / 1000.f;
	/* Ray cast     */
	if (m_raycast)
	{
		/* Prepare rays */
		const int       res = 64;
		const btScalar  fres = res - 1;
		const btScalar  size = 8;
		const btScalar  dist = 10;
		btTransform     trs;
		trs.setOrigin(ps);
		btScalar rayLength = 1000.f;

		const btScalar  angle = m_animtime * 0.2;
		trs.setRotation(btQuaternion(angle, SIMD_PI / 4, 0));
		btVector3   dir = trs.getBasis() * btVector3(0, -1, 0);
		trs.setOrigin(ps - dir * dist);
		btAlignedObjectArray<btVector3> origins;
		btAlignedObjectArray<btScalar>  fractions;
		origins.resize(res * res);
		fractions.resize(res * res, 1.f);
		for (int y = 0; y < res; ++y)
		{
			for (int x = 0; x < res; ++x)
			{
				const int   idx = y * res + x;
				origins[idx] = trs * btVector3(-size + size * 2 * x / fres, dist, -size + size * 2 * y / fres);
			}
		}
		/* Cast rays    */
		{
			m_clock.reset();
			if (sbs.size())
			{
				btVector3      *org = &origins[0];
				btScalar               *fraction = &fractions[0];
				btSoftBody            **psbs = &sbs[0];
				btSoftBody::sRayCast    results;
				for (int i = 0, ni = origins.size(), nb = sbs.size(); i < ni; ++i)
				{
					for (int ib = 0; ib < nb; ++ib)
					{
						btVector3 rayFrom = *org;
						btVector3 rayTo = rayFrom + dir * rayLength;
						if (psbs[ib]->rayTest(rayFrom, rayTo, results))
						{
							*fraction = results.fraction;
						}
					}
					++org; ++fraction;
				}
				long    ms = btMax<long>(m_clock.getTimeMilliseconds(), 1);
				long    rayperseconds = (1000 * (origins.size() * sbs.size())) / ms;
				printf("%d ms (%d rays/s)\r\n", int(ms), int(rayperseconds));
			}
		}
		/* Draw rays    */
		const btVector3 c[] = {   origins[0],
			origins[res - 1],
			origins[res * (res - 1)],
			origins[res * (res - 1) + res - 1]
		};
		idraw->drawLine(c[0], c[1], btVector3(0, 0, 0));
		idraw->drawLine(c[1], c[3], btVector3(0, 0, 0));
		idraw->drawLine(c[3], c[2], btVector3(0, 0, 0));
		idraw->drawLine(c[2], c[0], btVector3(0, 0, 0));
		for (int i = 0, ni = origins.size(); i < ni; ++i)
		{
			const btScalar      fraction = fractions[i];
			const btVector3    &org = origins[i];
			if (fraction < 1.f)
			{
				idraw->drawLine(org, org + dir * rayLength * fraction, btVector3(1, 0, 0));
			}
			else
			{
				idraw->drawLine(org, org - dir * rayLength * 0.1, btVector3(0, 0, 0));
			}
		}
#undef RES
	}
	/* Water level  */
	static const btVector3  axis[] = {btVector3(1, 0, 0),
		btVector3(0, 1, 0),
		btVector3(0, 0, 1)
	};
	if (m_softBodyWorldInfo.water_density > 0)
	{
		const btVector3 c =  btVector3((btScalar)0.25, (btScalar)0.25, 1);
		const btScalar  a =  (btScalar)0.5;
		const btVector3 n =  m_softBodyWorldInfo.water_normal;
		const btVector3 o =  -n * m_softBodyWorldInfo.water_offset;
		const btVector3 x =  btCross(n, axis[n.minAxis()]).normalized();
		const btVector3 y =  btCross(x, n).normalized();
		const btScalar  s =  25;
		idraw->drawTriangle(o - x * s - y * s, o + x * s - y * s, o + x * s + y * s, c, a);
		idraw->drawTriangle(o - x * s - y * s, o + x * s + y * s, o - x * s + y * s, c, a);
	}
	//

	int lineWidth = 305;
	int xStart = m_glutScreenWidth - lineWidth;
	int yStart = 20;

	
	btVector3 impact_in_bullet;
	impact_in_bullet.setX(impact_x);
	impact_in_bullet.setY(impact_y);
	impact_in_bullet.setZ(impact_z);

	if ((getDebugMode() & btIDebugDraw::DBG_NoHelpText) == 0)
	{
	    /* ***** ***** ***** ***** ***** Draw cutting route begin ***** ***** ***** ***** ***** */
		point_set_v[0] = (btVector3 *)malloc(sizeof(btVector3));
		point_set_v[1] = (btVector3 *)malloc(sizeof(btVector3));
		point_set_v[2] = (btVector3 *)malloc(sizeof(btVector3));
		point_set_h[0] = (btVector3 *)malloc(sizeof(btVector3));
		point_set_h[1] = (btVector3 *)malloc(sizeof(btVector3));
		point_set_h[2] = (btVector3 *)malloc(sizeof(btVector3));
		point_set_v[0] = star_point;
		point_set_v[1] = heart_point;
		point_set_v[2] = triangle_point;
		point_set_h[0] = pacman_point;
		point_set_h[1] = cross_point;
		point_set_h[2] = square_point;

		if ( teach_mode == true )
		{
			m_dynamicsWorld->setGravity(btVector3(0, 0, 0));
			m_softBodyWorldInfo.m_gravity.setValue(0, 0, 0);

			GLUquadric *myQuad = gluNewQuadric();
			
			tmp = (btVector3 *)malloc(sizeof(btVector3) * 7);

			if ( state == 5 )
			{
				cur_teaching_mode_v = cur_teaching_mode % 3;
				tmp = point_set_v[cur_teaching_mode_v];
			}
			else if ( state == 1 )
			{
				cur_teaching_mode_h = cur_teaching_mode % 3;
				tmp = point_set_h[cur_teaching_mode_h];
			}

            // Render cutting points using different colors (white and gray)
			if ( cur_point >= 0 && cur_point <= 6 )
			{
				glPushMatrix();
				glColor4f( 1, 1, 1, 1 );
				glTranslatef( tmp[cur_point].getX(), tmp[cur_point].getY(), tmp[cur_point].getZ() );
				gluSphere( myQuad, 0.1, 100, 100 );
				glPopMatrix();

				for ( int i = cur_point + 1; i < 7; i++ )
				{
					glPushMatrix();
					glColor4f( 0.467, 0.533, 0.6, 1 );
					glTranslatef( tmp[i].getX(), tmp[i].getY(), tmp[i].getZ() );
					gluSphere( myQuad, 0.1, 100, 100 );
					glPopMatrix();
				}
			}
		}
		/* ***** ***** ***** ***** ***** Draw cutting route  end  ***** ***** ***** ***** ***** */

		setOrthographicProjection();

		/* ***** ***** ***** ***** ***** Draw buttons begin ***** ***** ***** ***** ***** */
		std::string color;
		int screenPositionX = glutGet( (GLenum)GLUT_WINDOW_X );
		int screenPositionY = glutGet( (GLenum)GLUT_WINDOW_Y );
		POINT cursorPosition;
		btSoftBody *psb = sbs[0];

		glDisable( GL_LIGHTING );

		GetCursorPos( &cursorPosition );
		cursorPosition.x = cursorPosition.x - screenPositionX;
		cursorPosition.y = cursorPosition.y - screenPositionY;

		glPushMatrix();
		glColor3ub( 0, 0, 0 );
		//////////////////////////////////////// "Punch Hole" button (mode = 1)
		glTranslatef( m_glutScreenWidth - 771, 15, 0 );
		if ( (newFingerBonePositionX >= -3.4 && newFingerBonePositionX <= -1.8 && newFingerBonePositionY >= 6.5 && newFingerBonePositionY <= 7.7) || (cursorPosition.x >= m_glutScreenWidth - 771 && cursorPosition.x <= m_glutScreenWidth - 681 && cursorPosition.y >= 15 && cursorPosition.y <= 95) )
		{ // 23": -4.1 -2.7
			if ( timeCount1 >= 45 )
			{
				glColor3ub( 173, 255, 47 );
				if ( isMode1 == false )
				{
					current_mode = 1;
					cur_step = 4;
					isMode1 = true;

					m_dynamicsWorld->setGravity(btVector3(0, -10, 0));
					m_softBodyWorldInfo.m_gravity.setValue(0, -10, 0);
				}
			}
			if ( current_mode == 1 )
				glColor3ub( 173, 255, 47 );
			else
				glColor3ub( 0, 0, 0 );

			timeCount1++;
		}
		else
		{
			timeCount1 = 0;
			if ( current_mode == 1 )
				glColor3ub( 173, 255, 47 );
			else
				glColor3f( 0, 0, 0 );
			isMode1 = false;
		}
		DrawRoundedRec();

		//////////////////////////////////////// "Slice Paper" button (mode = 2 or mode = 5)
		glTranslatef( 160, 0, 0 );
		if ( (newFingerBonePositionX >= -5.8 && newFingerBonePositionX <= -4.2 && newFingerBonePositionY >= 6.5 && newFingerBonePositionY <= 7.7) || (cursorPosition.x >= m_glutScreenWidth - 611 && cursorPosition.x <= m_glutScreenWidth - 521 && cursorPosition.y >= 15 && cursorPosition.y <= 95) )
		{ // 23": -6.5 -5.1
			if ( timeCount2 >= 45 )
			{
				glColor3ub( 173, 255, 47 );
				if ( isMode2Or5 == false )
				{
					if ( current_mode != 2 && current_mode != 5 )
					{
						current_mode = 2;
						strcat( outputInfo, "You can slice paper horizontally now.\n" );
					}
					else if ( current_mode == 2 )
					{
						current_mode = 5;
						strcat( outputInfo, "You can slice paper vertically now.\n" );
					}
					else
					{
						current_mode = 2;
						strcat( outputInfo, "You can slice paper horizontally now.\n" );
					}
					cur_step = 4;
					isMode2Or5 = true;

					m_dynamicsWorld->setGravity(btVector3(0, -10, 0));
					m_softBodyWorldInfo.m_gravity.setValue(0, -10, 0);
				}
			}
			if ( current_mode == 2 || current_mode == 5 )
				glColor3ub( 173, 255, 47 );
			else
				glColor3ub( 0, 0, 0 );

			timeCount2++;
		}
		else
		{
			timeCount2 = 0;
			if ( current_mode == 2 || current_mode == 5 )
				glColor3ub( 173, 255, 47 );
			else
				glColor3f( 0, 0, 0 );
			isMode2Or5 = false;
		}
		DrawRoundedRec();

		//////////////////////////////////////// "Cut Paper Arbitrarily" button (mode = 3)
		glTranslatef( 160, 0, 0 );
		if ( (newFingerBonePositionX >= -8.2 && newFingerBonePositionX <= -6.6 && newFingerBonePositionY >= 6.5 && newFingerBonePositionY <= 7.7) || (cursorPosition.x >= m_glutScreenWidth - 451 && cursorPosition.x <= m_glutScreenWidth - 361 && cursorPosition.y >= 15 && cursorPosition.y <= 95) )
		{ // 23": -8.9 -7.5
			if ( timeCount3 >= 45 )
			{
				glColor3ub( 173, 255, 47 );
				if ( isMode3 == false && fold_condition != 2)
				{
					current_mode = 3;
					cur_step = 4;
					isMode3 = true;

					m_dynamicsWorld->setGravity(btVector3(0, 0, 0));
					m_softBodyWorldInfo.m_gravity.setValue(0, 0, 0);
					case_s = true;
					last_face[0] = -1;
					last_face[1] = -1;
					start_node = -1;
					psb->m_nodes[6].m_im = 0; psb->m_nodes[8].m_im = 0;
					psb->m_nodes[7].m_im = 0;
					psb->m_nodes[4].m_im = oldmass;
				}
			}
			if ( current_mode == 3 )
				glColor3ub( 173, 255, 47 );
			else
				glColor3ub( 0, 0, 0 );

			timeCount3++;
		}
		else
		{
			timeCount3 = 0;
			if ( current_mode == 3 )
				glColor3ub( 173, 255, 47 );
			else
				glColor3f( 0, 0, 0 );
			isMode3 = false;
		}
		DrawRoundedRec();

		//////////////////////////////////////// "Unfold Paper" button
		glTranslatef( 160, 0, 0 );
		glColor3ub( 0, 0, 0 );
		if ( (newFingerBonePositionX >= -10.6 && newFingerBonePositionX <= -9.0 && newFingerBonePositionY >= 6.5 && newFingerBonePositionY <= 7.7) || (cursorPosition.x >= m_glutScreenWidth - 291 && cursorPosition.x <= m_glutScreenWidth - 201 && cursorPosition.y >= 15 && cursorPosition.y <= 95) )
		{ // 23": -11.3 -9.9
			if ( timeCount4 >= 45 )
			{
				glColor3ub( 173, 255, 47 );
				if ( cur_fold == 1 )
				{
					unfold();
					cur_fold = 0;
				}
				else if ( cur_fold == 2 )
				{
					unfold();
					cur_fold = 0;
				}
			}

			timeCount4++;
		}
		else
		{
			timeCount4 = 0;
			glColor3f( 0, 0, 0 );
			if ( cur_fold == 1 )
				cur_fold = 2;
			else
				cur_fold = 1;
		}
		DrawRoundedRec();

		//////////////////////////////////////// "Teaching Mode" button (mode = 4)
		bool green = false;
		glTranslatef( 160, 0, 0 );
		if ( (newFingerBonePositionX >= -13.0 && newFingerBonePositionX <= -11.4 && newFingerBonePositionY >= 6.5 && newFingerBonePositionY <= 7.7) || (cursorPosition.x >= m_glutScreenWidth - 131 && cursorPosition.x <= m_glutScreenWidth - 41 && cursorPosition.y >= 15 && cursorPosition.y <= 95) )
		{ // 23": -13.7 -12.3
			if ( timeCount5 >= 45 )
			{
				if ( isMode4 == false )
				{
					if (fold_condition == 1)
					{
						if (state == 5 || state == 1)
						{
							teach_mode = true;
							current_mode = 3;
							isMode4 = true;
							cur_step = 4;

							m_dynamicsWorld->setGravity(btVector3(0, 0, 0));
							m_softBodyWorldInfo.m_gravity.setValue(0, 0, 0);
							case_s = true;
							last_face[0] = -1;
							last_face[1] = -1;
							start_node = -1;
							
							psb->m_nodes[6].m_im = 0; psb->m_nodes[8].m_im = 0;
							psb->m_nodes[7].m_im = 0;
							psb->m_nodes[4].m_im = oldmass;
							
						}
					}
				}
			}
			else if ( teach_mode == true )
			{
				glColor3ub( 173, 255, 47 );
				green = true;
			}

			if ( timeCount5 >= 45 && teach_mode == true )
			{
				glColor3ub( 173, 255, 47 );
				green = true;
			}

			timeCount5++;
		}
		else
		{
			timeCount5 = 0;
			if ( teach_mode == true )
			{
				glColor3ub( 173, 255, 47 );
				green = true;

			}
			else
			{
				glColor3f( 0, 0, 0 );
				green = false;
			}
			isMode4 = false;
		}
		DrawRoundedRec();
		glPopMatrix();

		//////////////////////////////////////// "Reset Scene" button (mode = 0)
		glPushMatrix();
		glTranslatef( 15, m_glutScreenHeight - 268, 0 );
		glScalef( 1.35, 1.35, 1.35 );
		glColor3ub( 0, 0, 0 );
		if ( (newIndexPositionX >= 125 && newIndexPositionX <= 224 && newIndexPositionY >= 778 && newIndexPositionY <= 870) || (newFingerBonePositionX >= 12.4 && newFingerBonePositionX <= 14.1 && newFingerBonePositionY >= -5.8 && newFingerBonePositionY <= -4.1) || (cursorPosition.x >= 15 && cursorPosition.x <= 135 && cursorPosition.y >= m_glutScreenHeight - 268 && cursorPosition.y <= m_glutScreenHeight - 160) )
		{
			if ( timeCount6 >= 45 )
			{
				glColor3ub( 173, 255, 47 );
				if ( isMode0 == false )
				{
					isMode0 = true;
					m_dynamicsWorld->setGravity(btVector3(0, -10, 0));
					m_softBodyWorldInfo.m_gravity.setValue(0, -10, 0);
					clientResetScene();
					cur_point = 0;
				}
			}

			timeCount6++;
		}
		else
		{
			timeCount6 = 0;
			glColor3f( 0, 0, 0 );
			isMode0 = false;
		}
		DrawRoundedRec();
		glPopMatrix();

		xStart  = m_glutScreenWidth - 751;
		yStart  = 50;
		if ( current_mode == 1 )
			color = "black";
		else color = "white";
		DrawInfo( xStart, yStart, "Punch", color );
		xStart += 6;
		yStart += 25;
		DrawInfo( xStart, yStart, "Hole", color );

		xStart += 160;
		yStart -= 25;
		if ( current_mode == 2 || current_mode == 5 )
			color = "black";
		else color = "white";
		DrawInfo( xStart, yStart, "Slice", color );
		xStart -= 3;
		yStart += 25;
		DrawInfo( xStart, yStart, "Paper", color );

		xStart += 167;
		yStart -= 35;
		if ( current_mode == 3 )
			color = "black";
		else color = "white";
		DrawInfo( xStart, yStart, "Cut", color );
		xStart -= 6;
		yStart += 20;
		DrawInfo( xStart, yStart, "Paper", color );
		xStart -= 17;
		yStart += 20;
		DrawInfo( xStart, yStart, "Arbitrarily", color );

		xStart += 170;
		yStart  = 50;
		if ( (newFingerBonePositionX >= -11.3 && newFingerBonePositionX <= -9.9 && newFingerBonePositionY >= 6.5 && newFingerBonePositionY <= 7.7) || (cursorPosition.x >= m_glutScreenWidth - 291 && cursorPosition.x <= m_glutScreenWidth - 201 && cursorPosition.y >= 15 && cursorPosition.y <= 95) )
		{
			if ( timeCount4 >= 45 )
				color = "black";
			else
				color = "white";
		}
		else
			color = "white";
		DrawInfo( xStart, yStart, "Unfold", color );
		xStart += 4;
		yStart += 25;
		DrawInfo( xStart, yStart, "Paper", color );

		xStart += 146;
		yStart  = 50;
		if ( green == true )
			color = "black";
		else color = "white";
		DrawInfo( xStart, yStart, "Teaching", color );
		xStart += 14;
		yStart += 25;
		DrawInfo( xStart, yStart, "Mode", color );

        /* ***** ***** ***** ***** ***** Draw shape list begin ***** ***** ***** ***** ***** */
		if ( teach_mode == true )
		{
			xStart = xStart - 50;
			yStart = yStart + 50;
			switch( state )
			{
			case 5:
				DrawInfo2( xStart, yStart      , " --------------" );
				if ( cur_teaching_mode_v == 0 ) DrawInfo2( xStart, yStart +  13, "| V | star     |" );
				else DrawInfo2( xStart, yStart +  13, "|   | star     |" );
				DrawInfo2( xStart, yStart +  26, " --------------" );
				if ( cur_teaching_mode_v == 1 ) DrawInfo2( xStart, yStart +  39, "| V | heart    |" );
				else DrawInfo2( xStart, yStart +  39, "|   | heart    |" );
				DrawInfo2( xStart, yStart +  52, " --------------" );
				if ( cur_teaching_mode_v == 2 ) DrawInfo2( xStart, yStart +  65, "| V | triangle |" );
				else DrawInfo2( xStart, yStart +  65, "|   | triangle |" );
				DrawInfo2( xStart, yStart +  78, " --------------" );
				break;
			case 1:
				DrawInfo2( xStart, yStart      , " --------------" );
				if ( cur_teaching_mode_h == 0 ) DrawInfo2( xStart, yStart +  13, "| V | pacman   |" );
				else DrawInfo2( xStart, yStart +  13, "|   | pacman   |" );
				DrawInfo2( xStart, yStart +  26, " --------------" );
				if ( cur_teaching_mode_h == 1 ) DrawInfo2( xStart, yStart +  39, "| V | cross    |" );
				else DrawInfo2( xStart, yStart +  39, "|   | cross    |" );
				DrawInfo2( xStart, yStart +  52, " --------------" );
				if ( cur_teaching_mode_h == 2 ) DrawInfo2( xStart, yStart +  65, "| V | square   |" );
				else DrawInfo2( xStart, yStart +  65, "|   | square   |" );
				DrawInfo2( xStart, yStart +  78, " --------------" );
				break;
			}
		}
		/* ***** ***** ***** ***** ***** Draw shape list  end  ***** ***** ***** ***** ***** */

		xStart  = 53;
		yStart  = m_glutScreenHeight - 220;
		if ( (newIndexPositionX >= 125 && newIndexPositionX <= 224 && newIndexPositionY >= 778 && newIndexPositionY <= 870) || (newFingerBonePositionX >= 12.4 && newFingerBonePositionX <= 14.1 && newFingerBonePositionY >= -5.8 && newFingerBonePositionY <= -4.1) || (cursorPosition.x >= 15 && cursorPosition.x <= 135 && cursorPosition.y >= m_glutScreenHeight - 268 && cursorPosition.y <= m_glutScreenHeight - 160) )
		{
			if ( timeCount6 >= 45 )
				color = "black";
			else
				color = "white";
		}
		else
			color = "white";
		DrawInfo( xStart, yStart, "Reset", color );
		xStart -= 3;
		yStart += 25;
		DrawInfo( xStart, yStart, "Scene", color );

        /* ***** ***** ***** ***** ***** Draw light begin ***** ***** ***** ***** ***** */
		glColor4f( 0, 0, 0, 1 );
		glPushMatrix();
		glLoadIdentity();
		glTranslatef( 500, m_glutScreenHeight - 50, 0 );
		glBegin( GL_QUADS );
		glVertex2d( 0, 0 );
		glVertex2d( 80, 0 );
		glVertex2d( 80, 30 );
		glVertex2d( 0, 30 );
		glEnd();
		glPopMatrix();

		if ( timeCount1 > 0 )
			DrawLight( timeCount1, m_glutScreenHeight );
		if ( timeCount2 > 0 )
			DrawLight( timeCount2, m_glutScreenHeight );
		if ( timeCount3 > 0 )
			DrawLight( timeCount3, m_glutScreenHeight );
		if ( timeCount4 > 0 )
			DrawLight( timeCount4, m_glutScreenHeight );
		if ( timeCount5 > 0 )
			DrawLight( timeCount5, m_glutScreenHeight );
		if ( timeCount6 > 0 )
			DrawLight( timeCount6, m_glutScreenHeight );

		color = "black";
		/* ***** ***** ***** ***** ***** Draw light  end  ***** ***** ***** ***** ***** */
		/* ***** ***** ***** ***** ***** Draw buttons  end  ***** ***** ***** ***** ***** */

        /* ***** ***** ***** ***** ***** Render output messages begin ***** ***** ***** ***** ***** */
		char tmp1[300];
		char tmp2[300];
		char *ch_ptr;
		int ch_cnt = 0;
		int j = 0;

		ch_ptr = strchr( outputInfo, '\n' );
		while ( ch_ptr )
		{
			ch_cnt++;
			ch_ptr = strchr( ch_ptr + 1, '\n' );
		}
		if ( ch_cnt > 13 ) // Always only render less than or equal to 13 lines
		{
			int diff = ch_cnt - 13;
			ch_ptr = outputInfo;
			while ( diff > 0 )
			{
				ch_ptr = strchr( ch_ptr, '\n' );
				ch_ptr++;
				diff--;
			}
			strcpy( outputInfo, ch_ptr );
		}

		xStart = 10;
		yStart = 20;

		DrawInfo( xStart, yStart, "-------------- Output Messages --------------", color );

		for ( int i = 0; i < 300; i++ )
		{
			tmp1[i] = '\0';
			tmp2[i] = '\0';
		}
		for ( int i = 0; i < 1000; i++ )
		{
			if ( outputInfo[i] == '\n' )
			{
				strcpy( tmp2, tmp1 );
				strcpy( tmp1, ">" );
				strcat( tmp1, tmp2 );
				yStart += 20;
				DrawInfo( xStart, yStart, tmp1, color );

				j = 0;
				for ( int k = 0; k < 300; k++ )
				{
					tmp1[k] = '\0';
					tmp2[k] = '\0';
				}
			}
			else if ( outputInfo[i] == '\0' )
				break;
			else
				tmp1[j++] = outputInfo[i];
		}
		/* ***** ***** ***** ***** ***** Render output messages  end  ***** ***** ***** ***** ***** */

		resetPerspectiveProjection();
		glEnable(GL_LIGHTING);
	}

	/* ***** ***** ***** ***** ***** Draw hands begin ***** ***** ***** ***** ***** */
	if ( !hasDisappeared && isRight )
	{
		glPushMatrix();
		ApplyLeftFrustum( 1 );
		glColorMask( true, false, false, true );
		gluLookAt( cp[0], cp[1], cp[2], ctp[0], ctp[1], ctp[2], cu[0], cu[1], cu[2] );
		DrawHands();

		glClear( GL_DEPTH_BUFFER_BIT );
		//glEnable( GL_BLEND );
		//glBlendFunc( GL_ONE,GL_ONE );

		ApplyRightFrustum( 1 );
		glColorMask( false, true, true, true );
		gluLookAt( cp[0], cp[1], cp[2], ctp[0], ctp[1], ctp[2], cu[0], cu[1], cu[2] );
		DrawHands();

		//glDisable( GL_BLEND );
		glColorMask( true, true, true, true );
		glPopMatrix();
	}
	/* ***** ***** ***** ***** ***** Draw hands  end  ***** ***** ***** ***** ***** */

	DemoApplication::renderme();
}

void    SoftDemo::setDrawClusters(bool drawClusters)
{
	if (drawClusters)
	{
		getSoftDynamicsWorld()->setDrawFlags(getSoftDynamicsWorld()->getDrawFlags() | fDrawFlags::Clusters);
	}
	else
	{
		getSoftDynamicsWorld()->setDrawFlags(getSoftDynamicsWorld()->getDrawFlags() & (~fDrawFlags::Clusters));
	}
}


void    SoftDemo::keyboardCallback(unsigned char key, int x, int y)
{
	/* ***** ***** ***** ***** ***** added code begin ***** ***** ***** ***** ***** */
	btSoftBodyArray    &sbs = getSoftDynamicsWorld()->getSoftBodyArray();
	btSoftBody *psb = sbs[0];
	btVector3 v(0, 0, 0); //速度
	btVector3 gravity = m_dynamicsWorld->getGravity();   //目前的重力值
	switch (key)
	{

	case ' ': case '0':
		strcat( outputInfo, "Enter \"Reset Scene\".\n" );
		clientResetScene();
		m_dynamicsWorld->setGravity(btVector3(0, -10, 0));
		m_softBodyWorldInfo.m_gravity.setValue(0, -10, 0);
		break;
	case    '1':    printf("打洞模式\n");
		m_dynamicsWorld->setGravity(btVector3(0, -10, 0));
		m_softBodyWorldInfo.m_gravity.setValue(0, -10, 0);
		strcat( outputInfo, "Enter \"Punch Hole\".\n" );
		current_mode = 1;
		cur_step = 4;
		break;
	case    '2':    printf("裁切模式\n");
		m_dynamicsWorld->setGravity(btVector3(0, -10, 0));
		m_softBodyWorldInfo.m_gravity.setValue(0, -10, 0);
		strcat( outputInfo, "Enter \"Slice Paper\".\n" );
		if (current_mode == 2)current_mode = 5;
		else current_mode = 2;
		cur_step = 4;
		break;
	case    '3':
		m_dynamicsWorld->setGravity(btVector3(0, 0, 0));
		m_softBodyWorldInfo.m_gravity.setValue(0, 0, 0);
		case_s = true;
		last_face[0] = -1;
		last_face[1] = -1;
		start_node = -1;
		printf("隨意剪取\n");
		strcat( outputInfo, "Enter \"Cut Paper Arbitrarily\".\n" );
		current_mode = 3;
		cur_step = 4;
		psb->m_nodes[6].m_im = 0; psb->m_nodes[8].m_im = 0;
		psb->m_nodes[7].m_im = 0;
		psb->m_nodes[4].m_im = oldmass;
		break;
	case    '4':    printf("切換視角\n");
		strcat( outputInfo, "Enter \"Change Viewpoint\".\n" );
		teach_mode = true;
		break;
	case    '5':
		cur_teaching_mode++;
		break;
	case    '+':    printf("上一步展開\n");
		strcat( outputInfo, "Enter \"Unfold Paper\".\n" );
		unfold();
		break;
	/* ***** ***** ***** ***** ***** added code  end  ***** ***** ***** ***** ***** */
	case    'd':    sDemoMode = !sDemoMode; break;
	case    'n':    motorcontrol.maxtorque = 10; motorcontrol.goal += 1; break;
	case    'm':    motorcontrol.maxtorque = 10; motorcontrol.goal -= 1; break;
	case    'l':    steercontrol_f.angle += 0.1; steercontrol_r.angle += 0.1; break;
	case    'k':    steercontrol_f.angle -= 0.1; steercontrol_r.angle -= 0.1; break;
	case    ']':    ++current_demo; clientResetScene(); break;
	case    '[':    --current_demo; clientResetScene(); break;
	case    ',':    m_raycast = !m_raycast; break;
	case    ';':    m_autocam = !m_autocam; break;
	case    'c':    getSoftDynamicsWorld()->setDrawFlags(getSoftDynamicsWorld()->getDrawFlags()^fDrawFlags::Clusters); break;
	case    '`':
		{
			btSoftBodyArray    &sbs = getSoftDynamicsWorld()->getSoftBodyArray();
			for (int ib = 0; ib < sbs.size(); ++ib)
			{
				btSoftBody *psb = sbs[ib];
				psb->staticSolve(128);
			}
		}
		break;
	default:        DemoApplication::keyboardCallback(key, x, y);
	}
}

//
void    SoftDemo::mouseMotionFunc(int x, int y)
{

	if (m_node && (m_results.fraction < 1.f))
	{
		if (!m_drag)
		{
#define SQ(_x_) (_x_)*(_x_)
			if ((SQ(x - m_lastmousepos[0]) + SQ(y - m_lastmousepos[1])) > 6)
			{
				m_drag = true;
			}
#undef SQ
		}
		if (m_drag)
		{
			m_lastmousepos[0]   =   x;
			m_lastmousepos[1]   =   y;
		}
	}
	else
	{
		DemoApplication::mouseMotionFunc(x, y);
	}
}


void SoftDemo::mouseMotionFunc2(int x, int y)
{
	btVector3 rayTo = getRayTo(x, y);
	btVector3 rayFrom  = m_cameraPosition;
	btVector3 dir = rayTo - rayFrom;
	newPivot = rayFrom + dir * 0.00091;

	if (m_node && (m_results.fraction < 1.f))
	{
		if (!m_drag)
		{
#define SQ(_x_) (_x_)*(_x_)
			if ((SQ(x - m_lastmousepos[0]) + SQ(y - m_lastmousepos[1])) > 6)
			{
				m_drag = true;
			}
#undef SQ
		}
		if (m_drag)
		{
			m_lastmousepos[0]   =   x;
			m_lastmousepos[1]   =   y;
		}
	}
	else
	{
		DemoApplication::mouseMotionFunc(x, y);
	}
}

void SoftDemo::mouseFunc2(int state, int x, int y)
{
	switch (state)
	{
	case 0:
		m_results.fraction = 1.f;
		DemoApplication::mouseFunc2(0, state, x, y);
		if (!m_pickConstraint)
		{
			const btVector3         rayFrom = m_cameraPosition;
			const btVector3         rayTo = getRayTo(x, y);
			const btVector3         rayDir = (rayTo - rayFrom).normalized();
			btSoftBodyArray        &sbs = getSoftDynamicsWorld()->getSoftBodyArray();
			for (int ib = 0; ib < sbs.size(); ++ib)
			{
				btSoftBody             *psb = sbs[ib];



				btSoftBody::sRayCast    res;
				if (psb->rayTest(rayFrom, rayTo, res))
				{
					m_results = res;
				}
			}
			if (m_results.fraction < 1.f)
			{
				m_impact            =   rayFrom + (rayTo - rayFrom) * m_results.fraction;
				btSoftBody             *psb = sbs[0];

				int size = 8;
				int s = 6;
				for (int i = 0; i <= size; i++)
				{
					if ((psb->m_nodes[i].m_x - m_impact).length2() <= s * s / 4) //9=(s/2)^2
					{
						record = i;         //to record the mouse position close to which node
					}
				}
				printf("cur_step:%d\n", cur_step);
				
				if (cur_step == 5) 
				{
					
					int i;
					
					if (record == 4)	//展開後點擊中心則重力恢復且mass恢復不再固定 使邊框掉下
					{
						for (i = 0; i < 9; i++)psb->m_nodes[i].m_im = oldmass;
						m_dynamicsWorld->setGravity(btVector3(0, -10, 0));
						m_softBodyWorldInfo.m_gravity.setValue(0, -10, 0);
					}
					else
					{
						
						m_dynamicsWorld->setGravity(btVector3(0, -10, 0));
						m_softBodyWorldInfo.m_gravity.setValue(0, -10, 0);	//展開後點擊非中心則重力恢復裡面物件掉下
						for (int i = 9; i < psb->m_nodes.size(); ++i)
						{
							btSoftBody::Node &n = psb->m_nodes[i];
							if (n.foldtag == 1 || n.tag == 1)
								n.m_im = oldmass;
						}
					}
				}
				if (cur_step == 0)
				{
					cur_step = 1;
					if (record != -1)
					{
						if ((record >= 0 && record < 2) || (record > 4 && record <= 6)) //產生對角線
						{
							psb->m_nodes[record + 2].m_im = 0;
							psb->m_nodes[size - record - 2].m_im = 0;
						}
						else if (record == 4)
						{
							cur_step = 0;
							m_results.fraction = 1.f;
							m_drag = false;
						}
						else
						{
							psb->m_nodes[record - 2].m_im = 0;
							psb->m_nodes[size - record + 2].m_im = 0;
						}

					}
					else
					{
						cur_step = 0;
						m_results.fraction = 1.f;
						m_drag = false;
					}
				}

				m_drag              =   m_cutting ? false : true;
				m_lastmousepos[0]   =   x;
				m_lastmousepos[1]   =   y;
				m_node              =   0;
				switch (m_results.feature)
				{
				case    btSoftBody::eFeature::Face:
					{
						btSoftBody::Face   &f = m_results.body->m_faces[m_results.index];
						m_node = f.m_n[0];
						for (int i = 1; i < 3; ++i)
						{
							if ( (m_node->m_x - m_impact).length2() >
								(f.m_n[i]->m_x - m_impact).length2())
							{
								m_node = f.m_n[i];
							}
						}
					}
					break;
				}
				if (m_node) m_goal = m_node->m_x;
				return;
			}
		}
		break;
	case 1:
		if ((!m_drag) && m_cutting && (m_results.fraction < 1.f))
		{
			/*ImplicitSphere    isphere(m_impact,1);
			printf("Mass before: %f\r\n",m_results.body->getTotalMass());
			m_results.body->refine(&isphere,0.01,true,m_impact);//0.0001
			printf("Mass after: %f\r\n",m_results.body->getTotalMass());*/



			/************************** 我們自己的剪取function : refine2 *******************************/
			//為設置起始點的狀態
			if (cur_step != 5)
			{
				if (current_mode == 3)
				{

					ImplicitSphere  isphere(m_impact, 0.05);
					if (case_s == true)
					{
						if (teach_mode == true && (tmp[cur_point] - m_impact).length2() < 0.2 && m_results.feature == btSoftBody::eFeature::Face)
						{
								
								m_impact.setX(tmp[cur_point].getX());
								m_impact.setY(tmp[cur_point].getY());
								m_impact.setZ(tmp[cur_point].getZ());
								near_teach_node = true;
						}


						if((teach_mode == true && near_teach_node == true) || (teach_mode == false && near_teach_node == false))
						{
							strcat( outputInfo, "near the teach node\n" );
							m_results.body->refine2(&isphere, 0.8, true, m_impact);
							if(teach_mode == true && near_teach_node == true)
								cur_point++;
						}

						near_teach_node = false;
						
					}
					//沒有設置成功, 請使用者在設置一次 , do nothing
					if (start_node == -1)
					{
						printf("設置失敗!請再近一些\n");
						strcat( outputInfo, "Ouch~ Please close the edge. -.-\n" );
					}
					// 設置成功了, do cutting
					else if (start_node != -1 && case_s == false)
						m_results.body->refine2(&isphere, 0.8, true, m_impact); // set refine function's cut = true

				}
				else
				{
					ImplicitSphere  isphere(m_impact, 0.8);
					m_results.body->refine(&isphere, 0.01, true, m_impact); //0.0001

				}
			}
			else cur_step = 4;
		}
		/************************** **************************** *******************************/
		m_results.fraction = 1.f;
		m_drag = false;
		DemoApplication::mouseFunc2(0, state, x, y);
		break;
	}
}
/* ***** ***** ***** ***** ***** added code  end  ***** ***** ***** ***** ***** */

//
void    SoftDemo::mouseFunc(int button, int state, int x, int y)
{
	if (button == 0)
	{
		switch (state)
		{
		case    0:
			{
				m_results.fraction = 1.f;
				DemoApplication::mouseFunc(button, state, x, y);
				if (!m_pickConstraint)
				{
					const btVector3         rayFrom = m_cameraPosition;
					const btVector3         rayTo = getRayTo(x, y);
					const btVector3         rayDir = (rayTo - rayFrom).normalized();
					btSoftBodyArray        &sbs = getSoftDynamicsWorld()->getSoftBodyArray();
					for (int ib = 0; ib < sbs.size(); ++ib)
					{
						btSoftBody             *psb = sbs[ib];
						btSoftBody::sRayCast    res;
						if (psb->rayTest(rayFrom, rayTo, res))
						{
							m_results = res;
						}
					}
					if (m_results.fraction < 1.f)
					{
						m_impact            =   rayFrom + (rayTo - rayFrom) * m_results.fraction;

						btSoftBody             *psb = sbs[0];

						int size = 8;
						int s = 6;
						for (int i = 0; i <= size; i++)
						{
							if ((psb->m_nodes[i].m_x - m_impact).length2() <= s * s / 4) //9=(s/2)^2
							{
								record = i;         //to record the mouse position close to which node
							}
						}
						printf("cur_step:%d\n", cur_step);
						if (cur_step == 5)
						{

							POINT polygon[1000];
							int k = 0;
							for (int i = 0; i < psb->m_nodes.size(); ++i)
							{
								btSoftBody::Node &n = psb->m_nodes[i];
								if (n.foldtag == 1)
								{
									polygon[k].x = n.m_x.getX();
									polygon[k].y = n.m_x.getY();
									k++;
								}
								if (n.tag == 1)
								{
									polygon[k].x = n.m_x.getX();
									polygon[k].y = n.m_x.getY();
									k++;
								}
								/*
								if(n.route==1)
								{
								polygon[k].x=n.m_x.getX();
								polygon[k].y=n.m_x.getY();
								k++;
								}
								*/
							}
							POINT p;
							p.x = m_impact.getX(); p.y = m_impact.getY();
							int i;
							//if(inside(polygon,k,p)==1)
							if (record == 4)
							{
								//strcat( outputInfo, "inside\n" );
								for (i = 0; i < 9; i++)psb->m_nodes[i].m_im = oldmass;
								m_dynamicsWorld->setGravity(btVector3(0, -10, 0));
								m_softBodyWorldInfo.m_gravity.setValue(0, -10, 0);
							}
							else
							{
								//strcat( outputInfo, "outside\n" );
								m_dynamicsWorld->setGravity(btVector3(0, -10, 0));
								m_softBodyWorldInfo.m_gravity.setValue(0, -10, 0);
								for (int i = 9; i < psb->m_nodes.size(); ++i)
								{
									btSoftBody::Node &n = psb->m_nodes[i];
									if (n.foldtag == 1 || n.tag == 1)
										n.m_im = oldmass;
								}
							}
						}

						if (cur_step == 0)
						{
							cur_step = 1;
							if (record != -1)
							{
								if ((record >= 0 && record < 2) || (record > 4 && record <= 6))
								{
									psb->m_nodes[record + 2].m_im = 0;
									psb->m_nodes[size - record - 2].m_im = 0;
								}
								else if (record == 4)
								{
									cur_step = 0;
									m_results.fraction = 1.f;
									m_drag = false;
								}
								else
								{
									psb->m_nodes[record - 2].m_im = 0;
									psb->m_nodes[size - record + 2].m_im = 0;
								}

							}
							else
							{
								cur_step = 0;
								m_results.fraction = 1.f;
								m_drag = false;
							}
						}

						printf("m_impact x:%f\ty:%f\tz:%f\n", m_impact.getX(), m_impact.getY(), m_impact.getZ());
						m_drag              =   m_cutting ? false : true;
						m_lastmousepos[0]   =   x;
						m_lastmousepos[1]   =   y;
						m_node              =   0;
						switch (m_results.feature)
						{
						case btSoftBody::eFeature::Tetra:
							{
								btSoftBody::Tetra  &tet = m_results.body->m_tetras[m_results.index];
								m_node = tet.m_n[0];
								for (int i = 1; i < 4; ++i)
								{
									if ( (m_node->m_x - m_impact).length2() >
										(tet.m_n[i]->m_x - m_impact).length2())
									{
										m_node = tet.m_n[i];
									}
								}
								break;
							}
						case    btSoftBody::eFeature::Face:
							{
								btSoftBody::Face   &f = m_results.body->m_faces[m_results.index];
								m_node = f.m_n[0];
								for (int i = 1; i < 3; ++i)
								{
									if ( (m_node->m_x - m_impact).length2() >
										(f.m_n[i]->m_x - m_impact).length2())
									{
										m_node = f.m_n[i];
									}
								}
							}
							break;
						}
						if (m_node) m_goal = m_node->m_x;
						return;
					}
				}
			}
			break;
		case    1:
			if ((!m_drag) && m_cutting && (m_results.fraction < 1.f))
			{
				/*ImplicitSphere    isphere(m_impact,1);
				printf("Mass before: %f\r\n",m_results.body->getTotalMass());
				m_results.body->refine(&isphere,0.01,true,m_impact);//0.0001
				printf("Mass after: %f\r\n",m_results.body->getTotalMass());*/



				/************************** 我們自己的剪取function : refine2 *******************************/
				//為設置起始點的狀態
				if (cur_step != 5)
				{
					if (current_mode == 3)
					{
						ImplicitSphere  isphere(m_impact, 0.05);
						if (case_s == true)
						{
							if (teach_mode == true && (tmp[cur_point] - m_impact).length2() < 0.2 && m_results.feature == btSoftBody::eFeature::Face)
							{
								
								m_impact.setX(tmp[cur_point].getX());
								m_impact.setY(tmp[cur_point].getY());
								m_impact.setZ(tmp[cur_point].getZ());
								near_teach_node = true;
							}

							if((teach_mode == true && near_teach_node == true) || (teach_mode == false && near_teach_node == false))
							{
								m_results.body->refine2(&isphere, 0.8, true, m_impact);
								if(teach_mode == true && near_teach_node == true)
									cur_point++;
							}

							near_teach_node = false;
							

						}
						//沒有設置成功, 請使用者在設置一次 , do nothing
						if (start_node == -1)
						{
							printf("設置失敗!請再近一些\n");
							strcat( outputInfo, "Ouch~ Please close the edge. -.-\n" );
						}
						// 設置成功了, do cutting
						else if (start_node != -1 && case_s == false)
							m_results.body->refine2(&isphere, 0.8, true, m_impact); // set refine function's cut = true

					}
					else
					{
						ImplicitSphere  isphere(m_impact, 0.8);
						m_results.body->refine(&isphere, 0.01, true, m_impact); //0.0001

					}
				}
				else cur_step = 4;
			}
			/************************** **************************** *******************************/
			m_results.fraction = 1.f;
			m_drag = false;
			DemoApplication::mouseFunc(button, state, x, y);
			break;
		}
	}
	else
	{
		DemoApplication::mouseFunc(button, state, x, y);
	}
}


void    SoftDemo::initPhysics()
{
	///create concave ground mesh


	m_azi = 0;

	//reset and disable motorcontrol at the start
	motorcontrol.goal = 0;
	motorcontrol.maxtorque = 0;

	btCollisionShape *groundShape = 0;
	{
		int i;
		int j;

		const int NUM_VERTS_X = 30;
		const int NUM_VERTS_Y = 30;
		const int totalVerts = NUM_VERTS_X * NUM_VERTS_Y;
		const int totalTriangles = 2 * (NUM_VERTS_X - 1) * (NUM_VERTS_Y - 1);

		gGroundVertices = new btVector3[totalVerts];
		gGroundIndices = new int[totalTriangles * 3];

		btScalar offset(-50);

		for ( i = 0; i < NUM_VERTS_X; i++)
		{
			for (j = 0; j < NUM_VERTS_Y; j++)
			{
				gGroundVertices[i + j * NUM_VERTS_X].setValue((i - NUM_VERTS_X * 0.5f)*TRIANGLE_SIZE,
					//0.f,
					waveheight * sinf((float)i)*cosf((float)j + offset),
					(j - NUM_VERTS_Y * 0.5f)*TRIANGLE_SIZE);
			}
		}

		int vertStride = sizeof(btVector3);
		int indexStride = 3 * sizeof(int);

		int index = 0;
		for ( i = 0; i < NUM_VERTS_X - 1; i++)
		{
			for (int j = 0; j < NUM_VERTS_Y - 1; j++)
			{
				gGroundIndices[index++] = j * NUM_VERTS_X + i;
				gGroundIndices[index++] = j * NUM_VERTS_X + i + 1;
				gGroundIndices[index++] = (j + 1) * NUM_VERTS_X + i + 1;

				gGroundIndices[index++] = j * NUM_VERTS_X + i;
				gGroundIndices[index++] = (j + 1) * NUM_VERTS_X + i + 1;
				gGroundIndices[index++] = (j + 1) * NUM_VERTS_X + i;
			}
		}

		btTriangleIndexVertexArray *indexVertexArrays = new btTriangleIndexVertexArray(totalTriangles,
			gGroundIndices,
			indexStride,
			totalVerts, (btScalar *) &gGroundVertices[0].x(), vertStride);

		bool useQuantizedAabbCompression = true;

		groundShape = new btBvhTriangleMeshShape(indexVertexArrays, useQuantizedAabbCompression);
		groundShape->setMargin(0.5);
	}

	m_collisionShapes.push_back(groundShape);

	btCollisionShape *groundBox = new btBoxShape (btVector3(100, CUBE_HALF_EXTENTS, 100));
	m_collisionShapes.push_back(groundBox);

	btCompoundShape *cylinderCompound = new btCompoundShape;
	btCollisionShape *cylinderShape = new btCylinderShape (btVector3(CUBE_HALF_EXTENTS, CUBE_HALF_EXTENTS, CUBE_HALF_EXTENTS));
	btTransform localTransform;
	localTransform.setIdentity();
	cylinderCompound->addChildShape(localTransform, cylinderShape);
	btQuaternion orn(btVector3(0, 1, 0), SIMD_PI);
	localTransform.setRotation(orn);
	cylinderCompound->addChildShape(localTransform, cylinderShape);

	m_collisionShapes.push_back(cylinderCompound);


	m_dispatcher = 0;

	///register some softbody collision algorithms on top of the default btDefaultCollisionConfiguration
	m_collisionConfiguration = new btSoftBodyRigidBodyCollisionConfiguration();


	m_dispatcher = new  btCollisionDispatcher(m_collisionConfiguration);
	m_softBodyWorldInfo.m_dispatcher = m_dispatcher;

	////////////////////////////
	///Register softbody versus softbody collision algorithm


	///Register softbody versus rigidbody collision algorithm


	////////////////////////////

	btVector3 worldAabbMin(-1000, -1000, -1000);
	btVector3 worldAabbMax(1000, 1000, 1000);

	m_broadphase = new btAxisSweep3(worldAabbMin, worldAabbMax, maxProxies);

	m_softBodyWorldInfo.m_broadphase = m_broadphase;

	btSequentialImpulseConstraintSolver *solver = new btSequentialImpulseConstraintSolver();

	m_solver = solver;

	btSoftBodySolver *softBodySolver = 0;
#ifdef USE_AMD_OPENCL

	static bool once = true;
	if (once)
	{
		once = false;
		initCL(0, 0);
	}

	if ( g_openCLSIMDSolver  )
		delete g_openCLSIMDSolver;
	if ( g_softBodyOutput )
		delete g_softBodyOutput;

	if (1)
	{
		g_openCLSIMDSolver = new btOpenCLSoftBodySolverSIMDAware( g_cqCommandQue, g_cxMainContext);
		//  g_openCLSIMDSolver = new btOpenCLSoftBodySolver( g_cqCommandQue, g_cxMainContext);
		g_openCLSIMDSolver->setCLFunctions(new CachingCLFunctions(g_cqCommandQue, g_cxMainContext));
	}



	softBodySolver = g_openCLSIMDSolver;
	g_softBodyOutput = new btSoftBodySolverOutputCLtoCPU;
#endif //USE_AMD_OPENCL

	btDiscreteDynamicsWorld *world = new btSoftRigidDynamicsWorld(m_dispatcher, m_broadphase, m_solver, m_collisionConfiguration, softBodySolver);
	m_dynamicsWorld = world;
	m_dynamicsWorld->setInternalTickCallback(pickingPreTickCallback, this, true);


	m_dynamicsWorld->getDispatchInfo().m_enableSPU = true;
	m_dynamicsWorld->setGravity(btVector3(0, -10, 0));
	m_softBodyWorldInfo.m_gravity.setValue(0, -10, 0);

	//  clientResetScene();

	m_softBodyWorldInfo.m_sparsesdf.Initialize();
	clientResetScene();
}






void    SoftDemo::exitPhysics()
{

	//cleanup in the reverse order of creation/initialization

	//remove the rigidbodies from the dynamics world and delete them
	int i;
	for (i = m_dynamicsWorld->getNumCollisionObjects() - 1; i >= 0 ; i--)
	{
		btCollisionObject *obj = m_dynamicsWorld->getCollisionObjectArray()[i];
		btRigidBody *body = btRigidBody::upcast(obj);
		if (body && body->getMotionState())
		{
			delete body->getMotionState();
		}
		m_dynamicsWorld->removeCollisionObject( obj );
		delete obj;
	}

	//delete collision shapes
	for (int j = 0; j < m_collisionShapes.size(); j++)
	{
		btCollisionShape *shape = m_collisionShapes[j];
		m_collisionShapes[j] = 0;
		delete shape;
	}

	//delete dynamics world
	delete m_dynamicsWorld;

	//delete solver
	delete m_solver;

	//delete broadphase
	delete m_broadphase;

	//delete dispatcher
	delete m_dispatcher;



	delete m_collisionConfiguration;


}






