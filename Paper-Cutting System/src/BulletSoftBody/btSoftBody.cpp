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

#include "btSoftBodyInternals.h"
#include "BulletSoftBody/btSoftBodySolvers.h"
#include "btSoftBodyData.h"
#include "LinearMath/btSerializer.h"
#include <stdio.h>
#include <algorithm>
#include <math.h>
#include "../../MyLib.h"
#include "../../command.h"

#define epss 1e-8
#define zero(x) (((x)>0?(x):-(x))<epss)
btAlignedObjectArray<btVector3> old_position;
btScalar oldmass;
int current_mode;
extern int index;      //not in cut
int cur_step;
int state;
int cutend;
//int prev_size;
bool PointinTriangle(btVector3 &A, btVector3 &B, btVector3 &C, btVector3 &P)
{
	btVector3 v0 = C - A ;
	btVector3 v1 = B - A ;
	btVector3 v2 = P - A ;

	float dot00 = btDot(v0, v0) ;
	float dot01 = btDot(v0, v1) ;
	float dot02 = btDot(v0, v2) ;
	float dot11 = btDot(v1, v1) ;
	float dot12 = btDot(v1, v2) ;

	float inverDeno = 1 / (dot00 * dot11 - dot01 * dot01) ;

	float u = (dot11 * dot02 - dot01 * dot12) * inverDeno ;
	if (u < 0 || u > 1) // if u out of range, return directly
	{
		return false ;
	}

	float v = (dot00 * dot12 - dot01 * dot02) * inverDeno ;
	if (v < 0 || v > 1) // if v out of range, return directly
	{
		return false ;
	}

	return u + v <= 1 ;
}
// refine用到的相關自建function
btScalar cross_2(btVector3 &a, btVector3 &b)
{
	return (btFabs(a.getX() * b.getY() - a.getY() * b.getX()));
}
btScalar cross(btVector3 &a, btVector3 &b)
{
	return ((a.getX() * b.getY() - a.getY() * b.getX()));
}
btVector3 getCross(btVector3 &a1, btVector3 &a2, btVector3 &b1, btVector3 &b2)
{
	btVector3 a = a2 - a1, b = b2 - b1, s = b1 - a1;
	btVector3 Cross(0, 0, 0);
	if (cross(a, b) == 0);
	else
		Cross = a1 + a * cross(s, b) / cross(a, b), a1.getZ();

	return Cross;
}
// 參數式
btVector3 getab(btVector3 &line1_a, btVector3 &line1_b)
{

	// y = ax + b
	float a = (line1_a.getY() - line1_b.getY()) / (line1_a.getX() - line1_b.getX());
	btScalar b = line1_a.getY() - a * (line1_a.getX());
	btVector3 ab(a, b, 0);
	return ab;
}
btScalar range(btScalar u, btScalar v, btScalar a, btScalar b, btScalar c, btScalar d)
{
	btScalar t;
	t = (c * u + v - d ) / ((b - d) - (a - c) * u);
	return t;
}

btVector3 getabc(btVector3 &line1_a, btVector3 &line1_b)
{
	// by = ax + c
	btScalar a = (line1_a.getY() - line1_b.getY());
	btScalar c = (line1_a.getY() - a * (line1_a.getX()))*(line1_a.getX() - line1_b.getX());
	btScalar b = (line1_a.getX() - line1_b.getX());
	//printf("%f\n",(line1_a.getX() - line1_b.getX()));
	btVector3 abc(a,b,c);
	return abc;
}

btScalar dis(btVector3& point,btVector3 &line1_a, btVector3 &line1_b)
{
	btScalar a,b,c,d;
	btVector3 abc = getabc(line1_a, line1_b);
	a = abc.getX();
	b = abc.getY();
	c = abc.getZ();	

	d = fabs(b*point.getY()-a*point.getX()-c) / std::sqrt(a*a+b*b);

	return d;
}
//計算交叉乘積(P1-P0)x(P2-P0)
double xmult(btVector3 &p1, btVector3 &p2, btVector3 &p0)
{
	return (p1.getX() - p0.getX()) * (p2.getY() - p0.getY()) - (p2.getX() - p0.getX()) * (p1.getY() - p0.getY());
}

//判斷p點是否在線段上
int dot_online_in(btVector3 &p, btVector3 &l1, btVector3 &l2)
{
	return zero(xmult(p, l1, l2)) && (l1.getX() - p.getX()) * (l2.getX() - p.getX()) < epss && (l1.getY() - p.getY()) * (l2.getY() - p.getY()) < epss;
}

//判斷線段兩端點是否在此線段同側, 點在線段上的話返回0
int same_side(btVector3 &p1, btVector3 &p2, btVector3 &l1, btVector3 &l2)
{
	return xmult(l1, p1, l2) * xmult(l1, p2, l2) > epss;
}

//判斷兩線段是否平行
int parallel(btVector3 &u1, btVector3 &u2, btVector3 &v1, btVector3 &v2)
{
	return zero((u1.getX() - u2.getX()) * (v1.getY() - v2.getY()) - (v1.getX() - v2.getX()) * (u1.getY() - u2.getY()));
}

//判斷三點是否共線
int dots_inline(btVector3 &p1, btVector3 &p2, btVector3 &p3)
{
	return zero(xmult(p1, p2, p3));
}

//判斷兩線段是否相交, 包括端點和部分重合
int intersect_in(btVector3 &u1, btVector3 &u2, btVector3 &v1, btVector3 &v2)
{

	//有交點, 但是為共線
	if (!dots_inline(u1, u2, v1) || !dots_inline(u1, u2, v2))
		return !same_side(u1, u2, v1, v2) && !same_side(v1, v2, u1, u2);
	return dot_online_in(u1, v1, v2) || dot_online_in(u2, v1, v2) || dot_online_in(v1, u1, u2) || dot_online_in(v2, u1, u2); // return 1的話代表有交點
}

//計算兩線段的交點：是否相交 + 是否平行
btVector3 intersection(btVector3 &u1, btVector3 &u2, btVector3 &v1, btVector3 &v2)
{
	btVector3 ret;
	ret = u1;

	btScalar t = ((u1.getX() - v1.getX()) * (v1.getY() - v2.getY()) - (u1.getY() - v1.getY()) * (v1.getX() - v2.getX()))
		/ ((u1.getX() - u2.getX()) * (v1.getY() - v2.getY()) - (u1.getY() - u2.getY()) * (v1.getX() - v2.getX()));

	ret.setX(ret.getX() + (u2.getX() - u1.getX())*t);
	ret.setY(ret.getY() + (u2.getY() - u1.getY())*t);
	return ret;
}

//////////////////////////////////////////////////////////////////////////
//
btSoftBody::btSoftBody(btSoftBodyWorldInfo *worldInfo, int node_count,  const btVector3 *x,  const btScalar *m)
	: m_softBodySolver(0), m_worldInfo(worldInfo)
{
	/* Init     */
	initDefaults();

	/* Default material */
	Material   *pm = appendMaterial();
	pm->m_kLST  =   1;
	pm->m_kAST  =   1;
	pm->m_kVST  =   1;
	pm->m_flags =   fMaterial::Default;

	/* Nodes            */
	const btScalar      margin = getCollisionShape()->getMargin();
	m_nodes.resize(node_count);
	for (int i = 0, ni = node_count; i < ni; ++i)
	{
		Node   &n = m_nodes[i];
		ZeroInitialize(n);
		n.m_x       =   x ? *x++ : btVector3(0, 0, 0);
		n.m_q       =   n.m_x;
		n.m_im      =   m ? *m++ : 1;
		n.m_im      =   n.m_im > 0 ? 1 / n.m_im : 0;
		n.m_leaf    =   m_ndbvt.insert(btDbvtVolume::FromCR(n.m_x, margin), &n);
		n.m_material =   pm;
	}
	updateBounds();

}

btSoftBody::btSoftBody(btSoftBodyWorldInfo *worldInfo)
	: m_worldInfo(worldInfo)
{
	initDefaults();
}


void    btSoftBody::initDefaults()
{
	m_internalType      =   CO_SOFT_BODY;
	m_cfg.aeromodel     =   eAeroModel::V_Point;
	m_cfg.kVCF          =   1;
	m_cfg.kDG           =   0;
	m_cfg.kLF           =   0;
	m_cfg.kDP           =   0;
	m_cfg.kPR           =   0;
	m_cfg.kVC           =   0;
	m_cfg.kDF           =   (btScalar)0.2;
	m_cfg.kMT           =   0;
	m_cfg.kCHR          =   (btScalar)1.0;
	m_cfg.kKHR          =   (btScalar)0.1;
	m_cfg.kSHR          =   (btScalar)1.0;
	m_cfg.kAHR          =   (btScalar)0.7;
	m_cfg.kSRHR_CL      =   (btScalar)0.1;
	m_cfg.kSKHR_CL      =   (btScalar)1;
	m_cfg.kSSHR_CL      =   (btScalar)0.5;
	m_cfg.kSR_SPLT_CL   =   (btScalar)0.5;
	m_cfg.kSK_SPLT_CL   =   (btScalar)0.5;
	m_cfg.kSS_SPLT_CL   =   (btScalar)0.5;
	m_cfg.maxvolume     =   (btScalar)1;
	m_cfg.timescale     =   1;
	m_cfg.viterations   =   0;
	m_cfg.piterations   =   1;
	m_cfg.diterations   =   0;
	m_cfg.citerations   =   4;
	m_cfg.collisions    =   fCollision::Default;
	m_pose.m_bvolume    =   false;
	m_pose.m_bframe     =   false;
	m_pose.m_volume     =   0;
	m_pose.m_com        =   btVector3(0, 0, 0);
	m_pose.m_rot.setIdentity();
	m_pose.m_scl.setIdentity();
	m_tag               =   0;
	m_timeacc           =   0;
	m_bUpdateRtCst      =   true;
	m_bounds[0]         =   btVector3(0, 0, 0);
	m_bounds[1]         =   btVector3(0, 0, 0);
	m_worldTransform.setIdentity();
	setSolver(eSolverPresets::Positions);

	/* Collision shape  */
	///for now, create a collision shape internally
	m_collisionShape = new btSoftBodyCollisionShape(this);
	m_collisionShape->setMargin(0.25f);

	m_initialWorldTransform.setIdentity();

	m_windVelocity = btVector3(0, 0, 0);
	m_restLengthScale = btScalar(1.0);
}

//
btSoftBody::~btSoftBody()
{
	//for now, delete the internal shape
	delete m_collisionShape;
	int i;

	releaseClusters();
	for (i = 0; i < m_materials.size(); ++i)
		btAlignedFree(m_materials[i]);
	for (i = 0; i < m_joints.size(); ++i)
		btAlignedFree(m_joints[i]);
}

//
bool            btSoftBody::checkLink(int node0, int node1) const
{
	return (checkLink(&m_nodes[node0], &m_nodes[node1]));
}

//
bool            btSoftBody::checkLink(const Node *node0, const Node *node1) const
{
	const Node *n[] = {node0, node1};
	for (int i = 0, ni = m_links.size(); i < ni; ++i)
	{
		const Link &l = m_links[i];
		if ( (l.m_n[0] == n[0] && l.m_n[1] == n[1]) ||
			(l.m_n[0] == n[1] && l.m_n[1] == n[0]))
		{
			return (true);
		}
	}
	return (false);
}

//
bool            btSoftBody::checkFace(int node0, int node1, int node2) const
{
	const Node *n[] = {   &m_nodes[node0],
		&m_nodes[node1],
		&m_nodes[node2]
	};
	for (int i = 0, ni = m_faces.size(); i < ni; ++i)
	{
		const Face &f = m_faces[i];
		int         c = 0;
		for (int j = 0; j < 3; ++j)
		{
			if ( (f.m_n[j] == n[0]) ||
				(f.m_n[j] == n[1]) ||
				(f.m_n[j] == n[2])) c |= 1 << j; else break;
		}
		if (c == 7) return (true);
	}
	return (false);
}

//
btSoftBody::Material       *btSoftBody::appendMaterial()
{
	Material   *pm = new(btAlignedAlloc(sizeof(Material), 16)) Material();
	if (m_materials.size() > 0)
		*pm = *m_materials[0];
	else
		ZeroInitialize(*pm);
	m_materials.push_back(pm);
	return (pm);
}

//
void            btSoftBody::appendNote( const char *text,
									   const btVector3 &o,
									   const btVector4 &c,
									   Node *n0,
									   Node *n1,
									   Node *n2,
									   Node *n3)
{
	Note    n;
	ZeroInitialize(n);
	n.m_rank        =   0;
	n.m_text        =   text;
	n.m_offset      =   o;
	n.m_coords[0]   =   c.x();
	n.m_coords[1]   =   c.y();
	n.m_coords[2]   =   c.z();
	n.m_coords[3]   =   c.w();
	n.m_nodes[0]    =   n0; n.m_rank += n0 ? 1 : 0;
	n.m_nodes[1]    =   n1; n.m_rank += n1 ? 1 : 0;
	n.m_nodes[2]    =   n2; n.m_rank += n2 ? 1 : 0;
	n.m_nodes[3]    =   n3; n.m_rank += n3 ? 1 : 0;
	m_notes.push_back(n);
}

//
void            btSoftBody::appendNote( const char *text,
									   const btVector3 &o,
									   Node *feature)
{
	appendNote(text, o, btVector4(1, 0, 0, 0), feature);
}

//
void            btSoftBody::appendNote( const char *text,
									   const btVector3 &o,
									   Link *feature)
{
	static const btScalar   w = 1 / (btScalar)2;
	appendNote(text, o, btVector4(w, w, 0, 0),   feature->m_n[0],
		feature->m_n[1]);
}

//
void            btSoftBody::appendNote( const char *text,
									   const btVector3 &o,
									   Face *feature)
{
	static const btScalar   w = 1 / (btScalar)3;
	appendNote(text, o, btVector4(w, w, w, 0),   feature->m_n[0],
		feature->m_n[1],
		feature->m_n[2]);
}

//
void            btSoftBody::appendNode( const btVector3 &x, btScalar m)
{
	if (m_nodes.capacity() == m_nodes.size()) //capacity為下次剪取保留空間 當capacity等於node數量等於空間不足。
	{
		pointersToIndices();
		m_nodes.reserve(m_nodes.size() * 2 + 1);
		indicesToPointers();
	}
	const btScalar  margin = getCollisionShape()->getMargin();
	m_nodes.push_back(Node());
	Node           &n = m_nodes[m_nodes.size() - 1];
	ZeroInitialize(n);
	n.m_x           =   x;
	n.m_q           =   n.m_x;
	n.m_im          =   m > 0 ? 1 / m : 0;
	n.m_material    =   m_materials[0];
	n.m_leaf        =   m_ndbvt.insert(btDbvtVolume::FromCR(n.m_x, margin), &n);
}

//
void            btSoftBody::appendLink(int model, Material *mat)
{
	Link    l;
	if (model >= 0)
		l = m_links[model];
	else
	{
		ZeroInitialize(l);
		l.m_material = mat ? mat : m_materials[0];
	}
	m_links.push_back(l);
}

//
void            btSoftBody::appendLink( int node0,
									   int node1,
									   Material *mat,
									   bool bcheckexist)
{
	appendLink(&m_nodes[node0], &m_nodes[node1], mat, bcheckexist);
}

//
void            btSoftBody::appendLink( Node *node0,
									   Node *node1,
									   Material *mat,
									   bool bcheckexist)
{
	if ((!bcheckexist) || (!checkLink(node0, node1)))
	{
		appendLink(-1, mat);
		Link   &l = m_links[m_links.size() - 1];
		l.m_n[0]        =   node0;
		l.m_n[1]        =   node1;
		l.m_rl          =   (l.m_n[0]->m_x - l.m_n[1]->m_x).length();
		m_bUpdateRtCst = true;
	}
}

//
void            btSoftBody::appendFace(int model, Material *mat)
{
	Face    f;
	if (model >= 0)
	{
		f = m_faces[model];
	}
	else
	{
		ZeroInitialize(f);
		f.m_material = mat ? mat : m_materials[0];
	}
	m_faces.push_back(f);
}

//
void            btSoftBody::appendFace(int node0, int node1, int node2, Material *mat)
{
	if (node0 == node1)
		return;
	if (node1 == node2)
		return;
	if (node2 == node0)
		return;

	appendFace(-1, mat);
	Face   &f = m_faces[m_faces.size() - 1];
	btAssert(node0 != node1);
	btAssert(node1 != node2);
	btAssert(node2 != node0);
	f.m_n[0]    =   &m_nodes[node0];
	f.m_n[1]    =   &m_nodes[node1];
	f.m_n[2]    =   &m_nodes[node2];
	f.m_ra      =   AreaOf( f.m_n[0]->m_x,
		f.m_n[1]->m_x,
		f.m_n[2]->m_x);
	m_bUpdateRtCst = true;
}

//
void            btSoftBody::appendTetra(int model, Material *mat)
{
	Tetra   t;
	if (model >= 0)
		t = m_tetras[model];
	else
	{
		ZeroInitialize(t);
		t.m_material = mat ? mat : m_materials[0];
	}
	m_tetras.push_back(t);
}

//
void            btSoftBody::appendTetra(int node0,
										int node1,
										int node2,
										int node3,
										Material *mat)
{
	appendTetra(-1, mat);
	Tetra  &t = m_tetras[m_tetras.size() - 1];
	t.m_n[0]    =   &m_nodes[node0];
	t.m_n[1]    =   &m_nodes[node1];
	t.m_n[2]    =   &m_nodes[node2];
	t.m_n[3]    =   &m_nodes[node3];
	t.m_rv      =   VolumeOf(t.m_n[0]->m_x, t.m_n[1]->m_x, t.m_n[2]->m_x, t.m_n[3]->m_x);
	m_bUpdateRtCst = true;
}

//

void            btSoftBody::appendAnchor(int node, btRigidBody *body, bool disableCollisionBetweenLinkedBodies, btScalar influence)
{
	btVector3 local = body->getWorldTransform().inverse() * m_nodes[node].m_x;
	appendAnchor(node, body, local, disableCollisionBetweenLinkedBodies, influence);
}

//
void            btSoftBody::appendAnchor(int node, btRigidBody *body, const btVector3 &localPivot, bool disableCollisionBetweenLinkedBodies, btScalar influence)
{
	if (disableCollisionBetweenLinkedBodies)
	{
		if (m_collisionDisabledObjects.findLinearSearch(body) == m_collisionDisabledObjects.size())
		{
			m_collisionDisabledObjects.push_back(body);
		}
	}

	Anchor  a;
	a.m_node            =   &m_nodes[node];
	a.m_body            =   body;
	a.m_local           =   localPivot;
	a.m_node->m_battach =   1;
	a.m_influence = influence;
	m_anchors.push_back(a);
}

//
void            btSoftBody::appendLinearJoint(const LJoint::Specs &specs, Cluster *body0, Body body1)
{
	LJoint     *pj  =   new(btAlignedAlloc(sizeof(LJoint), 16)) LJoint();
	pj->m_bodies[0] =   body0;
	pj->m_bodies[1] =   body1;
	pj->m_refs[0]   =   pj->m_bodies[0].xform().inverse() * specs.position;
	pj->m_refs[1]   =   pj->m_bodies[1].xform().inverse() * specs.position;
	pj->m_cfm       =   specs.cfm;
	pj->m_erp       =   specs.erp;
	pj->m_split     =   specs.split;
	m_joints.push_back(pj);
}

//
void            btSoftBody::appendLinearJoint(const LJoint::Specs &specs, Body body)
{
	appendLinearJoint(specs, m_clusters[0], body);
}

//
void            btSoftBody::appendLinearJoint(const LJoint::Specs &specs, btSoftBody *body)
{
	appendLinearJoint(specs, m_clusters[0], body->m_clusters[0]);
}

//
void            btSoftBody::appendAngularJoint(const AJoint::Specs &specs, Cluster *body0, Body body1)
{
	AJoint     *pj  =   new(btAlignedAlloc(sizeof(AJoint), 16)) AJoint();
	pj->m_bodies[0] =   body0;
	pj->m_bodies[1] =   body1;
	pj->m_refs[0]   =   pj->m_bodies[0].xform().inverse().getBasis() * specs.axis;
	pj->m_refs[1]   =   pj->m_bodies[1].xform().inverse().getBasis() * specs.axis;
	pj->m_cfm       =   specs.cfm;
	pj->m_erp       =   specs.erp;
	pj->m_split     =   specs.split;
	pj->m_icontrol  =   specs.icontrol;
	m_joints.push_back(pj);
}

//
void            btSoftBody::appendAngularJoint(const AJoint::Specs &specs, Body body)
{
	appendAngularJoint(specs, m_clusters[0], body);
}

//
void            btSoftBody::appendAngularJoint(const AJoint::Specs &specs, btSoftBody *body)
{
	appendAngularJoint(specs, m_clusters[0], body->m_clusters[0]);
}

//
void            btSoftBody::addForce(const btVector3 &force)
{
	for (int i = 0, ni = m_nodes.size(); i < ni; ++i) addForce(force, i);
}

//
void            btSoftBody::addForce(const btVector3 &force, int node)
{
	Node   &n = m_nodes[node];
	if (n.m_im > 0)
	{
		n.m_f   +=  force;
	}
}

void            btSoftBody::addAeroForceToNode(const btVector3 &windVelocity, int nodeIndex)
{
	btAssert(nodeIndex >= 0 && nodeIndex < m_nodes.size());

	const btScalar dt = m_sst.sdt;
	const btScalar kLF = m_cfg.kLF;
	const btScalar kDG = m_cfg.kDG;
	//const btScalar kPR = m_cfg.kPR;
	//const btScalar kVC = m_cfg.kVC;
	const bool as_lift = kLF > 0;
	const bool as_drag = kDG > 0;
	const bool as_aero = as_lift || as_drag;
	const bool as_vaero = as_aero && (m_cfg.aeromodel < btSoftBody::eAeroModel::F_TwoSided);

	Node &n = m_nodes[nodeIndex];

	if ( n.m_im > 0 )
	{
		btSoftBody::sMedium medium;

		EvaluateMedium(m_worldInfo, n.m_x, medium);
		medium.m_velocity = windVelocity;
		medium.m_density = m_worldInfo->air_density;

		/* Aerodynamics         */
		if (as_vaero)
		{
			const btVector3 rel_v = n.m_v - medium.m_velocity;
			const btScalar rel_v_len = rel_v.length();
			const btScalar  rel_v2 = rel_v.length2();

			if (rel_v2 > SIMD_EPSILON)
			{
				const btVector3 rel_v_nrm = rel_v.normalized();
				btVector3   nrm = n.m_n;

				if (m_cfg.aeromodel == btSoftBody::eAeroModel::V_TwoSidedLiftDrag)
				{
					nrm *= (btScalar)( (btDot(nrm, rel_v) < 0) ? -1 : +1);
					btVector3 fDrag(0, 0, 0);
					btVector3 fLift(0, 0, 0);

					btScalar n_dot_v = nrm.dot(rel_v_nrm);
					btScalar tri_area = 0.5f * n.m_area;

					fDrag = 0.5f * kDG * medium.m_density * rel_v2 * tri_area * n_dot_v * (-rel_v_nrm);

					// Check angle of attack
					// cos(10? = 0.98480
					if ( 0 < n_dot_v && n_dot_v < 0.98480f)
						fLift = 0.5f * kLF * medium.m_density * rel_v_len * tri_area * btSqrt(1.0f - n_dot_v * n_dot_v) * (nrm.cross(rel_v_nrm).cross(rel_v_nrm));

					// Check if the velocity change resulted by aero drag force exceeds the current velocity of the node.
					btVector3 del_v_by_fDrag = fDrag * n.m_im * m_sst.sdt;
					btScalar del_v_by_fDrag_len2 = del_v_by_fDrag.length2();
					btScalar v_len2 = n.m_v.length2();

					if (del_v_by_fDrag_len2 >= v_len2 && del_v_by_fDrag_len2 > 0)
					{
						btScalar del_v_by_fDrag_len = del_v_by_fDrag.length();
						btScalar v_len = n.m_v.length();
						fDrag *= btScalar(0.8) * (v_len / del_v_by_fDrag_len);
					}

					n.m_f += fDrag;
					n.m_f += fLift;
				}
				else if (m_cfg.aeromodel == btSoftBody::eAeroModel::V_Point || m_cfg.aeromodel == btSoftBody::eAeroModel::V_OneSided || m_cfg.aeromodel == btSoftBody::eAeroModel::V_TwoSided)
				{
					if (btSoftBody::eAeroModel::V_TwoSided)
						nrm *= (btScalar)( (btDot(nrm, rel_v) < 0) ? -1 : +1);

					const btScalar dvn = btDot(rel_v, nrm);
					/* Compute forces   */
					if (dvn > 0)
					{
						btVector3       force(0, 0, 0);
						const btScalar  c0  =   n.m_area * dvn * rel_v2 / 2;
						const btScalar  c1  =   c0 * medium.m_density;
						force   +=  nrm * (-c1 * kLF);
						force   +=  rel_v.normalized() * (-c1 * kDG);
						ApplyClampedForce(n, force, dt);
					}
				}
			}
		}
	}
}

void            btSoftBody::addAeroForceToFace(const btVector3 &windVelocity, int faceIndex)
{
	const btScalar dt = m_sst.sdt;
	const btScalar kLF = m_cfg.kLF;
	const btScalar kDG = m_cfg.kDG;
	//  const btScalar kPR = m_cfg.kPR;
	//  const btScalar kVC = m_cfg.kVC;
	const bool as_lift = kLF > 0;
	const bool as_drag = kDG > 0;
	const bool as_aero = as_lift || as_drag;
	const bool as_faero = as_aero && (m_cfg.aeromodel >= btSoftBody::eAeroModel::F_TwoSided);

	if (as_faero)
	{
		btSoftBody::Face   &f = m_faces[faceIndex];

		btSoftBody::sMedium medium;

		const btVector3 v = (f.m_n[0]->m_v + f.m_n[1]->m_v + f.m_n[2]->m_v) / 3;
		const btVector3 x = (f.m_n[0]->m_x + f.m_n[1]->m_x + f.m_n[2]->m_x) / 3;
		EvaluateMedium(m_worldInfo, x, medium);
		medium.m_velocity = windVelocity;
		medium.m_density = m_worldInfo->air_density;
		const btVector3 rel_v = v - medium.m_velocity;
		const btScalar rel_v_len = rel_v.length();
		const btScalar  rel_v2 = rel_v.length2();

		if (rel_v2 > SIMD_EPSILON)
		{
			const btVector3 rel_v_nrm = rel_v.normalized();
			btVector3   nrm = f.m_normal;

			if (m_cfg.aeromodel == btSoftBody::eAeroModel::F_TwoSidedLiftDrag)
			{
				nrm *= (btScalar)( (btDot(nrm, rel_v) < 0) ? -1 : +1);

				btVector3 fDrag(0, 0, 0);
				btVector3 fLift(0, 0, 0);

				btScalar n_dot_v = nrm.dot(rel_v_nrm);
				btScalar tri_area = 0.5f * f.m_ra;

				fDrag = 0.5f * kDG * medium.m_density * rel_v2 * tri_area * n_dot_v * (-rel_v_nrm);

				// Check angle of attack
				// cos(10? = 0.98480
				if ( 0 < n_dot_v && n_dot_v < 0.98480f)
					fLift = 0.5f * kLF * medium.m_density * rel_v_len * tri_area * btSqrt(1.0f - n_dot_v * n_dot_v) * (nrm.cross(rel_v_nrm).cross(rel_v_nrm));

				fDrag /= 3;
				fLift /= 3;

				for (int j = 0; j < 3; ++j)
				{
					if (f.m_n[j]->m_im > 0)
					{
						// Check if the velocity change resulted by aero drag force exceeds the current velocity of the node.
						btVector3 del_v_by_fDrag = fDrag * f.m_n[j]->m_im * m_sst.sdt;
						btScalar del_v_by_fDrag_len2 = del_v_by_fDrag.length2();
						btScalar v_len2 = f.m_n[j]->m_v.length2();

						if (del_v_by_fDrag_len2 >= v_len2 && del_v_by_fDrag_len2 > 0)
						{
							btScalar del_v_by_fDrag_len = del_v_by_fDrag.length();
							btScalar v_len = f.m_n[j]->m_v.length();
							fDrag *= btScalar(0.8) * (v_len / del_v_by_fDrag_len);
						}

						f.m_n[j]->m_f += fDrag;
						f.m_n[j]->m_f += fLift;
					}
				}
			}
			else if (m_cfg.aeromodel == btSoftBody::eAeroModel::F_OneSided || m_cfg.aeromodel == btSoftBody::eAeroModel::F_TwoSided)
			{
				if (btSoftBody::eAeroModel::F_TwoSided)
					nrm *= (btScalar)( (btDot(nrm, rel_v) < 0) ? -1 : +1);

				const btScalar  dvn = btDot(rel_v, nrm);
				/* Compute forces   */
				if (dvn > 0)
				{
					btVector3       force(0, 0, 0);
					const btScalar  c0  =   f.m_ra * dvn * rel_v2;
					const btScalar  c1  =   c0 * medium.m_density;
					force   +=  nrm * (-c1 * kLF);
					force   +=  rel_v.normalized() * (-c1 * kDG);
					force   /=  3;
					for (int j = 0; j < 3; ++j) ApplyClampedForce(*f.m_n[j], force, dt);
				}
			}
		}
	}

}

//
void            btSoftBody::addVelocity(const btVector3 &velocity)
{
	for (int i = 0, ni = m_nodes.size(); i < ni; ++i) addVelocity(velocity, i);
}

/* Set velocity for the entire body                                     */
void                btSoftBody::setVelocity(    const btVector3 &velocity)
{
	for (int i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		Node   &n = m_nodes[i];
		if (n.m_im > 0)
		{
			n.m_v   =   velocity;
		}
	}
}


//
void            btSoftBody::addVelocity(const btVector3 &velocity, int node)
{
	Node   &n = m_nodes[node];
	if (n.m_im > 0)
	{
		n.m_v   +=  velocity;
	}
}

//
void            btSoftBody::setMass(int node, btScalar mass)
{
	m_nodes[node].m_im = mass > 0 ? 1 / mass : 0;
	m_bUpdateRtCst = true;
}

//
btScalar        btSoftBody::getMass(int node) const
{
	return (m_nodes[node].m_im > 0 ? 1 / m_nodes[node].m_im : 0);
}

//
btScalar        btSoftBody::getTotalMass() const
{
	btScalar    mass = 0;
	for (int i = 0; i < m_nodes.size(); ++i)
	{
		mass += getMass(i);
	}
	return (mass);
}

//
void            btSoftBody::setTotalMass(btScalar mass, bool fromfaces)
{
	int i;

	if (fromfaces)
	{

		for (i = 0; i < m_nodes.size(); ++i)
		{
			m_nodes[i].m_im = 0;
		}
		for (i = 0; i < m_faces.size(); ++i)
		{
			const Face     &f = m_faces[i];
			const btScalar  twicearea = AreaOf(   f.m_n[0]->m_x,
				f.m_n[1]->m_x,
				f.m_n[2]->m_x);
			for (int j = 0; j < 3; ++j)
			{
				f.m_n[j]->m_im += twicearea;
			}
		}
		for ( i = 0; i < m_nodes.size(); ++i)
		{
			m_nodes[i].m_im = 1 / m_nodes[i].m_im;
		}
	}
	const btScalar  tm = getTotalMass();
	const btScalar  itm = 1 / tm;
	for ( i = 0; i < m_nodes.size(); ++i)
	{
		m_nodes[i].m_im /= itm * mass;
	}
	m_bUpdateRtCst = true;
}

//
void            btSoftBody::setTotalDensity(btScalar density)
{
	setTotalMass(getVolume()*density, true);
}

//
void            btSoftBody::setVolumeMass(btScalar mass)
{
	btAlignedObjectArray<btScalar>  ranks;
	ranks.resize(m_nodes.size(), 0);
	int i;

	for (i = 0; i < m_nodes.size(); ++i)
	{
		m_nodes[i].m_im = 0;
	}
	for (i = 0; i < m_tetras.size(); ++i)
	{
		const Tetra &t = m_tetras[i];
		for (int j = 0; j < 4; ++j)
		{
			t.m_n[j]->m_im += btFabs(t.m_rv);
			ranks[int(t.m_n[j] - &m_nodes[0])] += 1;
		}
	}
	for ( i = 0; i < m_nodes.size(); ++i)
	{
		if (m_nodes[i].m_im > 0)
		{
			m_nodes[i].m_im = ranks[i] / m_nodes[i].m_im;
		}
	}
	setTotalMass(mass, false);
}

//
void            btSoftBody::setVolumeDensity(btScalar density)
{
	btScalar    volume = 0;
	for (int i = 0; i < m_tetras.size(); ++i)
	{
		const Tetra &t = m_tetras[i];
		for (int j = 0; j < 4; ++j)
		{
			volume += btFabs(t.m_rv);
		}
	}
	setVolumeMass(volume * density / 6);
}

//
void            btSoftBody::transform(const btTransform &trs)
{
	const btScalar  margin = getCollisionShape()->getMargin();
	ATTRIBUTE_ALIGNED16(btDbvtVolume)   vol;

	for (int i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		Node   &n = m_nodes[i];
		n.m_x = trs * n.m_x;
		n.m_q = trs * n.m_q;
		n.m_n = trs.getBasis() * n.m_n;
		vol = btDbvtVolume::FromCR(n.m_x, margin);

		m_ndbvt.update(n.m_leaf, vol);
	}
	updateNormals();
	updateBounds();
	updateConstants();
	m_initialWorldTransform = trs;
}

//
void            btSoftBody::translate(const btVector3 &trs)
{
	btTransform t;
	t.setIdentity();
	t.setOrigin(trs);
	transform(t);
}

//
void            btSoftBody::rotate( const btQuaternion &rot)
{
	btTransform t;
	t.setIdentity();
	t.setRotation(rot);
	transform(t);
}

//
void            btSoftBody::scale(const btVector3 &scl)
{

	const btScalar  margin = getCollisionShape()->getMargin();
	ATTRIBUTE_ALIGNED16(btDbvtVolume)   vol;

	for (int i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		Node   &n = m_nodes[i];
		n.m_x *= scl;
		n.m_q *= scl;
		vol = btDbvtVolume::FromCR(n.m_x, margin);
		m_ndbvt.update(n.m_leaf, vol);
	}
	updateNormals();
	updateBounds();
	updateConstants();
}

//
btScalar btSoftBody::getRestLengthScale()
{
	return m_restLengthScale;
}

//
void btSoftBody::setRestLengthScale(btScalar restLengthScale)
{
	for (int i = 0, ni = m_links.size(); i < ni; ++i)
	{
		Link       &l = m_links[i];
		l.m_rl  =   l.m_rl / m_restLengthScale * restLengthScale;
		l.m_c1  =   l.m_rl * l.m_rl;
	}
	m_restLengthScale = restLengthScale;

	if (getActivationState() == ISLAND_SLEEPING)
		activate();
}

//
void            btSoftBody::setPose(bool bvolume, bool bframe)
{
	m_pose.m_bvolume    =   bvolume;
	m_pose.m_bframe     =   bframe;
	int i, ni;

	/* Weights      */
	const btScalar  omass = getTotalMass();
	const btScalar  kmass = omass * m_nodes.size() * 1000;
	btScalar        tmass = omass;
	m_pose.m_wgh.resize(m_nodes.size());
	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		if (m_nodes[i].m_im <= 0) tmass += kmass;
	}
	for ( i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		Node   &n = m_nodes[i];
		m_pose.m_wgh[i] =    n.m_im > 0                    ?
			1 / (m_nodes[i].m_im * tmass)   :
		kmass / tmass;
	}
	/* Pos      */
	const btVector3 com = evaluateCom();
	m_pose.m_pos.resize(m_nodes.size());
	for ( i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		m_pose.m_pos[i] = m_nodes[i].m_x - com;
	}
	m_pose.m_volume =   bvolume ? getVolume() : 0;
	m_pose.m_com    =   com;
	m_pose.m_rot.setIdentity();
	m_pose.m_scl.setIdentity();
	/* Aqq      */
	m_pose.m_aqq[0] =
		m_pose.m_aqq[1] =
		m_pose.m_aqq[2] =   btVector3(0, 0, 0);
	for ( i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		const btVector3    &q = m_pose.m_pos[i];
		const btVector3     mq = m_pose.m_wgh[i] * q;
		m_pose.m_aqq[0] += mq.x() * q;
		m_pose.m_aqq[1] += mq.y() * q;
		m_pose.m_aqq[2] += mq.z() * q;
	}
	m_pose.m_aqq = m_pose.m_aqq.inverse();

	updateConstants();
}

void                btSoftBody::resetLinkRestLengths()
{
	for (int i = 0, ni = m_links.size(); i < ni; ++i)
	{
		Link &l =   m_links[i];
		l.m_rl  =   (l.m_n[0]->m_x - l.m_n[1]->m_x).length();
		l.m_c1  =   l.m_rl * l.m_rl;
	}
}

//
btScalar        btSoftBody::getVolume() const
{
	btScalar    vol = 0;
	if (m_nodes.size() > 0)
	{
		int i, ni;

		const btVector3 org = m_nodes[0].m_x;
		for (i = 0, ni = m_faces.size(); i < ni; ++i)
		{
			const Face &f = m_faces[i];
			vol += btDot(f.m_n[0]->m_x - org, btCross(f.m_n[1]->m_x - org, f.m_n[2]->m_x - org));
		}
		vol /= (btScalar)6;
	}
	return (vol);
}

//
int             btSoftBody::clusterCount() const
{
	return (m_clusters.size());
}

//
btVector3       btSoftBody::clusterCom(const Cluster *cluster)
{
	btVector3       com(0, 0, 0);
	for (int i = 0, ni = cluster->m_nodes.size(); i < ni; ++i)
	{
		com += cluster->m_nodes[i]->m_x * cluster->m_masses[i];
	}
	return (com * cluster->m_imass);
}

//
btVector3       btSoftBody::clusterCom(int cluster) const
{
	return (clusterCom(m_clusters[cluster]));
}

//
btVector3       btSoftBody::clusterVelocity(const Cluster *cluster, const btVector3 &rpos)
{
	return (cluster->m_lv + btCross(cluster->m_av, rpos));
}

//
void            btSoftBody::clusterVImpulse(Cluster *cluster, const btVector3 &rpos, const btVector3 &impulse)
{
	const btVector3 li = cluster->m_imass * impulse;
	const btVector3 ai = cluster->m_invwi * btCross(rpos, impulse);
	cluster->m_vimpulses[0] += li; cluster->m_lv += li;
	cluster->m_vimpulses[1] += ai; cluster->m_av += ai;
	cluster->m_nvimpulses++;
}

//
void            btSoftBody::clusterDImpulse(Cluster *cluster, const btVector3 &rpos, const btVector3 &impulse)
{
	const btVector3 li = cluster->m_imass * impulse;
	const btVector3 ai = cluster->m_invwi * btCross(rpos, impulse);
	cluster->m_dimpulses[0] += li;
	cluster->m_dimpulses[1] += ai;
	cluster->m_ndimpulses++;
}

//
void            btSoftBody::clusterImpulse(Cluster *cluster, const btVector3 &rpos, const Impulse &impulse)
{
	if (impulse.m_asVelocity)    clusterVImpulse(cluster, rpos, impulse.m_velocity);
	if (impulse.m_asDrift)       clusterDImpulse(cluster, rpos, impulse.m_drift);
}

//
void            btSoftBody::clusterVAImpulse(Cluster *cluster, const btVector3 &impulse)
{
	const btVector3 ai = cluster->m_invwi * impulse;
	cluster->m_vimpulses[1] += ai; cluster->m_av += ai;
	cluster->m_nvimpulses++;
}

//
void            btSoftBody::clusterDAImpulse(Cluster *cluster, const btVector3 &impulse)
{
	const btVector3 ai = cluster->m_invwi * impulse;
	cluster->m_dimpulses[1] += ai;
	cluster->m_ndimpulses++;
}

//
void            btSoftBody::clusterAImpulse(Cluster *cluster, const Impulse &impulse)
{
	if (impulse.m_asVelocity)    clusterVAImpulse(cluster, impulse.m_velocity);
	if (impulse.m_asDrift)       clusterDAImpulse(cluster, impulse.m_drift);
}

//
void            btSoftBody::clusterDCImpulse(Cluster *cluster, const btVector3 &impulse)
{
	cluster->m_dimpulses[0] += impulse * cluster->m_imass;
	cluster->m_ndimpulses++;
}

struct NodeLinks
{
	btAlignedObjectArray<int> m_links;
};



//
int             btSoftBody::generateBendingConstraints(int distance, Material *mat)
{
	int i, j;

	if (distance > 1)
	{
		/* Build graph  */
		const int       n = m_nodes.size();
		const unsigned  inf = (~(unsigned)0) >> 1;
		unsigned       *adj = new unsigned[n * n];


#define IDX(_x_,_y_)    ((_y_)*n+(_x_))
		for (j = 0; j < n; ++j)
		{
			for (i = 0; i < n; ++i)
			{
				if (i != j)
				{
					adj[IDX(i, j)] = adj[IDX(j, i)] = inf;
				}
				else
				{
					adj[IDX(i, j)] = adj[IDX(j, i)] = 0;
				}
			}
		}
		for ( i = 0; i < m_links.size(); ++i)
		{
			const int   ia = (int)(m_links[i].m_n[0] - &m_nodes[0]);
			const int   ib = (int)(m_links[i].m_n[1] - &m_nodes[0]);
			adj[IDX(ia, ib)] = 1;
			adj[IDX(ib, ia)] = 1;
		}


		//special optimized case for distance == 2
		if (distance == 2)
		{

			btAlignedObjectArray<NodeLinks> nodeLinks;


			/* Build node links */
			nodeLinks.resize(m_nodes.size());

			for ( i = 0; i < m_links.size(); ++i)
			{
				const int   ia = (int)(m_links[i].m_n[0] - &m_nodes[0]);
				const int   ib = (int)(m_links[i].m_n[1] - &m_nodes[0]);
				if (nodeLinks[ia].m_links.findLinearSearch(ib) == nodeLinks[ia].m_links.size())
					nodeLinks[ia].m_links.push_back(ib);

				if (nodeLinks[ib].m_links.findLinearSearch(ia) == nodeLinks[ib].m_links.size())
					nodeLinks[ib].m_links.push_back(ia);
			}
			for (int ii = 0; ii < nodeLinks.size(); ii++)
			{
				int i = ii;

				for (int jj = 0; jj < nodeLinks[ii].m_links.size(); jj++)
				{
					int k = nodeLinks[ii].m_links[jj];
					for (int kk = 0; kk < nodeLinks[k].m_links.size(); kk++)
					{
						int j = nodeLinks[k].m_links[kk];
						if (i != j)
						{
							const unsigned  sum = adj[IDX(i, k)] + adj[IDX(k, j)];
							btAssert(sum == 2);
							if (adj[IDX(i, j)] > sum)
							{
								adj[IDX(i, j)] = adj[IDX(j, i)] = sum;
							}
						}

					}
				}
			}
		}
		else
		{
			///generic Floyd's algorithm
			for (int k = 0; k < n; ++k)
			{
				for (j = 0; j < n; ++j)
				{
					for (i = j + 1; i < n; ++i)
					{
						const unsigned  sum = adj[IDX(i, k)] + adj[IDX(k, j)];
						if (adj[IDX(i, j)] > sum)
						{
							adj[IDX(i, j)] = adj[IDX(j, i)] = sum;
						}
					}
				}
			}
		}


		/* Build links  */
		int nlinks = 0;
		for (j = 0; j < n; ++j)
		{
			for (i = j + 1; i < n; ++i)
			{
				if (adj[IDX(i, j)] == (unsigned)distance)
				{
					appendLink(i, j, mat);
					m_links[m_links.size() - 1].m_bbending = 1;
					++nlinks;
				}
			}
		}
		delete[] adj;
		return (nlinks);
	}
	return (0);
}

//
void            btSoftBody::randomizeConstraints()
{
	unsigned long   seed = 243703;
#define NEXTRAND (seed=(1664525L*seed+1013904223L)&0xffffffff)
	int i, ni;

	for (i = 0, ni = m_links.size(); i < ni; ++i)
	{
		btSwap(m_links[i], m_links[NEXTRAND % ni]);
	}
	for (i = 0, ni = m_faces.size(); i < ni; ++i)
	{
		btSwap(m_faces[i], m_faces[NEXTRAND % ni]);
	}
#undef NEXTRAND
}

//
void            btSoftBody::releaseCluster(int index)
{
	Cluster    *c = m_clusters[index];
	if (c->m_leaf) m_cdbvt.remove(c->m_leaf);
	c->~Cluster();
	btAlignedFree(c);
	m_clusters.remove(c);
}

//
void            btSoftBody::releaseClusters()
{
	while (m_clusters.size() > 0) releaseCluster(0);
}

//
int             btSoftBody::generateClusters(int k, int maxiterations)
{
	int i;
	releaseClusters();
	m_clusters.resize(btMin(k, m_nodes.size()));
	for (i = 0; i < m_clusters.size(); ++i)
	{
		m_clusters[i]           =   new(btAlignedAlloc(sizeof(Cluster), 16)) Cluster();
		m_clusters[i]->m_collide =   true;
	}
	k = m_clusters.size();
	if (k > 0)
	{
		/* Initialize       */
		btAlignedObjectArray<btVector3> centers;
		btVector3                       cog(0, 0, 0);
		int                             i;
		for (i = 0; i < m_nodes.size(); ++i)
		{
			cog += m_nodes[i].m_x;
			m_clusters[(i * 29873) % m_clusters.size()]->m_nodes.push_back(&m_nodes[i]);
		}
		cog /= (btScalar)m_nodes.size();
		centers.resize(k, cog);
		/* Iterate          */
		const btScalar  slope = 16;
		bool            changed;
		int             iterations = 0;
		do
		{
			const btScalar  w = 2 - btMin<btScalar>(1, iterations / slope);
			changed = false;
			iterations++;
			int i;

			for (i = 0; i < k; ++i)
			{
				btVector3   c(0, 0, 0);
				for (int j = 0; j < m_clusters[i]->m_nodes.size(); ++j)
				{
					c += m_clusters[i]->m_nodes[j]->m_x;
				}
				if (m_clusters[i]->m_nodes.size())
				{
					c           /=  (btScalar)m_clusters[i]->m_nodes.size();
					c           =   centers[i] + (c - centers[i]) * w;
					changed     |=  ((c - centers[i]).length2() > SIMD_EPSILON);
					centers[i]  =   c;
					m_clusters[i]->m_nodes.resize(0);
				}
			}
			for (i = 0; i < m_nodes.size(); ++i)
			{
				const btVector3 nx = m_nodes[i].m_x;
				int             kbest = 0;
				btScalar        kdist = ClusterMetric(centers[0], nx);
				for (int j = 1; j < k; ++j)
				{
					const btScalar  d = ClusterMetric(centers[j], nx);
					if (d < kdist)
					{
						kbest = j;
						kdist = d;
					}
				}
				m_clusters[kbest]->m_nodes.push_back(&m_nodes[i]);
			}
		}
		while (changed && (iterations < maxiterations));
		/* Merge        */
		btAlignedObjectArray<int>   cids;
		cids.resize(m_nodes.size(), -1);
		for (i = 0; i < m_clusters.size(); ++i)
		{
			for (int j = 0; j < m_clusters[i]->m_nodes.size(); ++j)
			{
				cids[int(m_clusters[i]->m_nodes[j] - &m_nodes[0])] = i;
			}
		}
		for (i = 0; i < m_faces.size(); ++i)
		{
			const int idx[] = {   int(m_faces[i].m_n[0] - &m_nodes[0]),
				int(m_faces[i].m_n[1] - &m_nodes[0]),
				int(m_faces[i].m_n[2] - &m_nodes[0])
			};
			for (int j = 0; j < 3; ++j)
			{
				const int cid = cids[idx[j]];
				for (int q = 1; q < 3; ++q)
				{
					const int kid = idx[(j + q) % 3];
					if (cids[kid] != cid)
					{
						if (m_clusters[cid]->m_nodes.findLinearSearch(&m_nodes[kid]) == m_clusters[cid]->m_nodes.size())
						{
							m_clusters[cid]->m_nodes.push_back(&m_nodes[kid]);
						}
					}
				}
			}
		}
		/* Master       */
		if (m_clusters.size() > 1)
		{
			Cluster    *pmaster = new(btAlignedAlloc(sizeof(Cluster), 16)) Cluster();
			pmaster->m_collide  =   false;
			pmaster->m_nodes.reserve(m_nodes.size());
			for (int i = 0; i < m_nodes.size(); ++i) pmaster->m_nodes.push_back(&m_nodes[i]);
			m_clusters.push_back(pmaster);
			btSwap(m_clusters[0], m_clusters[m_clusters.size() - 1]);
		}
		/* Terminate    */
		for (i = 0; i < m_clusters.size(); ++i)
		{
			if (m_clusters[i]->m_nodes.size() == 0)
			{
				releaseCluster(i--);
			}
		}
	}
	else
	{
		//create a cluster for each tetrahedron (if tetrahedra exist) or each face
		if (m_tetras.size())
		{
			m_clusters.resize(m_tetras.size());
			for (i = 0; i < m_clusters.size(); ++i)
			{
				m_clusters[i]           =   new(btAlignedAlloc(sizeof(Cluster), 16)) Cluster();
				m_clusters[i]->m_collide =   true;
			}
			for (i = 0; i < m_tetras.size(); i++)
			{
				for (int j = 0; j < 4; j++)
				{
					m_clusters[i]->m_nodes.push_back(m_tetras[i].m_n[j]);
				}
			}

		}
		else
		{
			m_clusters.resize(m_faces.size());
			for (i = 0; i < m_clusters.size(); ++i)
			{
				m_clusters[i]           =   new(btAlignedAlloc(sizeof(Cluster), 16)) Cluster();
				m_clusters[i]->m_collide =   true;
			}

			for (i = 0; i < m_faces.size(); ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					m_clusters[i]->m_nodes.push_back(m_faces[i].m_n[j]);
				}
			}
		}
	}

	if (m_clusters.size())
	{
		initializeClusters();
		updateClusters();


		//for self-collision
		m_clusterConnectivity.resize(m_clusters.size()*m_clusters.size());
		{
			for (int c0 = 0; c0 < m_clusters.size(); c0++)
			{
				m_clusters[c0]->m_clusterIndex = c0;
				for (int c1 = 0; c1 < m_clusters.size(); c1++)
				{

					bool connected = false;
					Cluster *cla = m_clusters[c0];
					Cluster *clb = m_clusters[c1];
					for (int i = 0; !connected && i < cla->m_nodes.size(); i++)
					{
						for (int j = 0; j < clb->m_nodes.size(); j++)
						{
							if (cla->m_nodes[i] == clb->m_nodes[j])
							{
								connected = true;
								break;
							}
						}
					}
					m_clusterConnectivity[c0 + c1 * m_clusters.size()] = connected;
				}
			}
		}
	}

	return (m_clusters.size());
}

/******************************************** new code ********************************************/
int start_node = -1;  
int start_dir; // 記錄使用者剪取的起始位置 ( 0 = from down , 1 = from right , 2 = from left , 3 = from up)
int last_face[2] = {-1,-1};

btVector3 store_inter,start;

float impact_x=0;
float impact_y=0;
float impact_z=0;

// 結構記錄每次剪取的資訊, 以供下次剪取利用
struct cut_info
{
	int Nimpact;     // 剪取的位置 
	int Nrcut;       // 剪取後設置的right node編號
	int Nlcut;       // 剪取後設置的left node編號 
	int Nright;      // 剪取到的Face(三角形), 其右邊的node編號 
	int Nleft;       // 剪取到的Face(三角形), 其左邊的node編號 
	int Nup;         // 剪取到的Face(三角形), 其上邊的node編號   

};

struct cut_info info[100];

int one,another;

// 仿照原refine內新增node重設質量的code, 為我們新的refine重新定義一個函式出來方便呼叫
btScalar btSoftBody::cal_m(Node &a, Node &b, btScalar t,Node &impact)
{
	btScalar        m = 0;
	if (a.m_im > 0)
	{
		// 如果b點沒有被固定住
		if (b.m_im > 0)
		{
			const btScalar  ma = 1 / a.m_im; // ma = a點的mass
			const btScalar  mb = 1 / b.m_im; // mb = b點的mass
			const btScalar  mc = Lerp(ma, mb, t);
			const btScalar  f = (ma + mb) / (ma + mb + mc);

			a.m_im = 1 / (ma * f);
			b.m_im = 1 / (mb * f);
	
			m = mc * f;
		}
		else
		{
			a.m_im /= 0.5f;
			m = 1 / a.m_im;
		}
	}
	else
	{
		//printf("a.m_im <= 0\n");
		if (b.m_im > 0)
		{
			b.m_im /= 0.5f;
			m = 1 / b.m_im;
		}
		else
		{
			impact.m_im /= 0.5f;
			m = 1 / impact.m_im;
			printf("hey hey hey\n");
		}
	}
	return m;
}

// 我們的cutting method
void			btSoftBody::refine2(ImplicitFn* ifn,btScalar accurary,bool cut,btVector3& impact)
{
	const Node*			nbase = &m_nodes[0];   // 取得m_nodes陣列的位址， tNodeArray m_nodes is btAlignedObjectArray< T > Class Template  
	int					ncount = m_nodes.size();  // size() , return the number of elements in the array, m_nodes陣列的大小
	btSymMatrix<int>	edges(ncount,-2);   
	int					newnodes=0;
	int i,j,k,ni;

	btVector3 V(0,0,0);
	cur_step=4;
	for(i=0;i<9;i++)
		if(m_nodes[i].m_im==0)
			m_nodes[i].m_im=oldmass;

	/* 根據不同對摺狀態, 將該固定的node固定住 */
	if(state==1)  //下往上摺
	{
		m_nodes[0].m_im=0;m_nodes[2].m_im=0;
		m_nodes[6].m_im=0;m_nodes[8].m_im=0;
		m_nodes[1].m_im=0;m_nodes[7].m_im=0;
	}
	else if(state==2)   //左下右上摺
	{
		m_nodes[0].m_im=0;m_nodes[2].m_im=0;
		m_nodes[6].m_im=0;m_nodes[8].m_im=0;
	}
	else if(state==3)  //右往左摺
	{
		m_nodes[3].m_im=0;m_nodes[5].m_im=0;
		m_nodes[1].m_im=0;m_nodes[7].m_im=0;
		//m_nodes[6].m_im=0;m_nodes[8].m_im=0;
	}
	else if(state==4)  //右下往左上摺
	{
		m_nodes[0].m_im=0;m_nodes[2].m_im=0;
		m_nodes[6].m_im=0;m_nodes[8].m_im=0;
	}
	else if(state==5)  //左往右摺
	{
		m_nodes[3].m_im=0;m_nodes[5].m_im=0;
		m_nodes[1].m_im=0;m_nodes[7].m_im=0;
	}
	else if(state==6)   //右上左下摺
	{
		m_nodes[0].m_im=0;m_nodes[2].m_im=0;
		m_nodes[6].m_im=0;m_nodes[8].m_im=0;
	}
	else if(state==7)  //上往下摺
	{
		m_nodes[3].m_im=0;m_nodes[5].m_im=0;
		m_nodes[1].m_im=0;m_nodes[7].m_im=0;
	}
	else if(state==8)  //左上往右下摺
	{
		m_nodes[0].m_im=0;m_nodes[2].m_im=0;
		m_nodes[6].m_im=0;m_nodes[8].m_im=0;
	}
	else
	{
		m_nodes[6].m_im=0;m_nodes[7].m_im=0;m_nodes[8].m_im=0;
	}



	/*************************** 設置剪取起始點 **********************************/

	// 紀錄四個邊界
	btVector3 corner0,corner1,corner2,corner3;

	/* 若是摺一次的情況下做剪取, 必須要知道現在的對摺狀態, 得知邊界位置才能判斷使用者剪取方向 */
	if(state == 0) // 無對折情況
	{
		corner0 = m_nodes[0].m_x;
		corner1 = m_nodes[2].m_x;
		corner2 = m_nodes[6].m_x;
		corner3 = m_nodes[8].m_x;
	}
	if(state == 1) // 由下往上摺
	{
		corner0 = m_nodes[3].m_x;
		corner1 = m_nodes[5].m_x;
		corner2 = m_nodes[6].m_x;

	}
	if(state == 3) // 由右往左摺
	{
		corner0 = m_nodes[1].m_x;
		corner1 = m_nodes[2].m_x;
		corner2 = m_nodes[7].m_x;
	}
	if(state == 5) // 由左往右摺
	{
		corner0 = m_nodes[0].m_x;
		corner1 = m_nodes[1].m_x;
		corner2 = m_nodes[6].m_x;
	}
	if(state == 7) // 由上往下摺
	{
		corner0 = m_nodes[0].m_x;
		corner1 = m_nodes[2].m_x;
		corner2 = m_nodes[3].m_x;

	}

	btScalar allow = 0.7;   // 設置靠近邊界的容忍值

	if(start_node == -1) // 起始點還沒設置成功的狀態
	{
		// 如果剪取的位置離邊界的距離 <= allow值 則判定為靠近邊界

		if(fabs(impact.getY()-corner0.getY()) <= allow)
		{
			start = impact;
			start.setY(corner0.getY());
			start_node = 1;
			start_dir = 0;
			printf("設置成功, 請由你的下方往上剪\n");
			//strcat( outputInfo, "Setting successful! Please cut upwards.\n" );
		}
		// 靠近右邊 (bullet的右邊)
		else if(fabs(impact.getX()-corner1.getX()) <= allow)
		{		
			start = impact;
			start.setX(corner1.getX());
			start_node = 1;
			start_dir = 1;
			printf("設置成功, 請由你的左方往右剪\n");
			//strcat( outputInfo, "Setting successful! Please cut rightwards.\n" );
		}
		// 靠近左邊 (bullet的左邊)
		else if(fabs(impact.getX()-corner0.getX()) <= allow)
		{
			start = impact;
			start.setX(corner0.getX());
			start_node = 1;
			start_dir = 2;
			printf("設置成功, 請由你的右方往左剪\n");
			//strcat( outputInfo, "Setting successful! Please cut leftwards.\n" );
		}
		// 靠近上邊
		else if(fabs(impact.getY()-corner2.getY()) <= allow)
		{
			start = impact;
			start.setY(corner2.getY());
			start_node = 1;
			start_dir = 3;
			printf("設置成功, 請由你的上方往下剪\n");
			//strcat( outputInfo, "Setting successful! Please cut downwards.\n" );
		}
	}
	
	
	// 對摺的情形下, 摺線是否被剪到, 對分割及連接Face有很大的影響

	int cut_foldline = 0;   // 在對摺的情況下, 是否以摺線為剪取起點
	int end_foldline = 0;   // 在對摺的情況下, 是否以摺線為剪取終點 
	

	/***************** 對摺情況下, 如果以摺線當起始點 ******************/

	if((state == 1 && start_dir == 0)
		|| (state == 5 && start_dir == 1)
		|| (state == 3 && start_dir == 2)
		|| (state == 7 && start_dir == 3))
		cut_foldline = 1; // 記錄是否以摺線當起點

	/******************************************************************/

	// 這些是測試及DEBUG用的輸出檔
	FILE *fp,*fr,*frr;
	fp = fopen("EvalData.txt","w");
	fr = fopen("FaceAndLinkData.txt","w");
	frr = fopen("FaceNodeData.txt","w");


	int beding_cnt=0;

	// Filter out & Fill edges 除了做fprintf, 其餘是直接延用refine的東西, code本身重要性不大

	/* Filter out */    
	// 將m_beding的link先pop掉 
	for(i=0;i<m_links.size();++i) 
	{
		Link&	l = m_links[i];
		//   結構中的初始值為 int m_bbending:1 ( Bending link)
		if(l.m_bbending)
		{
			beding_cnt++;

			fprintf(fp,"m_bbending = 1, link %d %d\n",l.m_n[0]-nbase,l.m_n[1]-nbase);
			if(!SameSign(ifn->Eval(l.m_n[0]->m_x), ifn->Eval(l.m_n[1]->m_x)))
			{
				btSwap(m_links[i], m_links[m_links.size()-1]);  // 單純做template型態的swap, 將m_link第i個陣列的東西, 跟最後一個陣列的東西交換
				m_links.pop_back(); // pop_back() 就是將m_size--, 並做destructor
				--i; 
			}
		}	
	}
	fclose(fp);

	/* Fill edges		*/ 
	// 將目前場上有的link都tag上-1
	for(i=0;i<m_links.size();++i)
	{
		Link&	l = m_links[i];
		fprintf(fr,"第%d個link:l.m_n[0]-nbase = %d\tl.m_n[1]-nbase = %d\n",i+1,int(l.m_n[0]-nbase),int(l.m_n[1]-nbase));

		edges(int(l.m_n[0]-nbase),int(l.m_n[1]-nbase)) = -1;
	}
	for(i=0;i<m_faces.size();++i)
	{	
		Face&	f = m_faces[i];

		fprintf(fr,"第%d個face:\n",i+1);
		fprintf(fr,"f.m_n[0]-nbase = %d\tf.m_n[1]-nbase = %d\n",f.m_n[0]-nbase,f.m_n[1]-nbase);
		fprintf(fr,"f.m_n[1]-nbase = %d\tf.m_n[2]-nbase = %d\n",f.m_n[1]-nbase,f.m_n[2]-nbase);
		fprintf(fr,"f.m_n[2]-nbase = %d\tf.m_n[0]-nbase = %d\n",f.m_n[2]-nbase,f.m_n[0]-nbase);

		edges(int(f.m_n[0]-nbase),int(f.m_n[1]-nbase)) = -1;
		edges(int(f.m_n[1]-nbase),int(f.m_n[2]-nbase)) = -1;
		edges(int(f.m_n[2]-nbase),int(f.m_n[0]-nbase)) = -1;
	}

	// 印出每個node中的成員資訊
	for(i=0;i<m_nodes.size();i++)
		fprintf(fr,"%d Node Features:\nm_x:%f %f %f\nm_q:%f %f %f\nm_v:%f %f %f\nm_f:%f %f %f\nm_n:%f %f %f\nm_im:%f\nm_area:%f\nm_leaf:%d\nm_battach:%d\n\n",i,
		m_nodes[i].m_x.getX(),m_nodes[i].m_x.getY(),m_nodes[i].m_x.getZ(),
		m_nodes[i].m_q.getX(),m_nodes[i].m_q.getY(),m_nodes[i].m_q.getZ(),
		m_nodes[i].m_v.getX(),m_nodes[i].m_v.getY(),m_nodes[i].m_v.getZ(),
		m_nodes[i].m_f.getX(),m_nodes[i].m_f.getY(),m_nodes[i].m_f.getZ(),
		m_nodes[i].m_n.getX(),m_nodes[i].m_n.getY(),m_nodes[i].m_n.getZ(),
		m_nodes[i].m_im,
		m_nodes[i].m_area,
		m_nodes[i].m_leaf->dataAsInt,
		m_nodes[i].m_battach);

	fclose(fr);

	
	// impact從SoftDemo傳進來refine2, 是滑鼠點擊的位置座標 => 影響點
	printf("impact: %f %f %f\n",impact.getX(),impact.getY(),impact.getZ());

	// 是否為設置起點狀態, 沒有設置成功不要造成剪取
	if(start_node != -1)
	{
		// 開二維是應付對摺一次的情況, 因為對摺的情況做剪取, 一次會影響到兩塊Face, 也就是說需要跑兩次迴圈去處理形變, 記錄也就需要做兩份
		
		int face[2][3];  // 記錄目前點擊到的face的三個node編號 
		int facenum[2] ; // 記錄目前點擊到的face編號 
		int f=0;		 // 記錄目前點擊到的face有幾個(對摺需要) 

		// 用三角形面積法判定影響點所在的Face是哪幾塊	
		for(i=0;i<m_faces.size();i++)
		{
			if(fabs((cross_2(m_faces[i].m_n[0]->m_x-impact,m_faces[i].m_n[1]->m_x-impact)+cross_2(m_faces[i].m_n[1]->m_x-impact,m_faces[i].m_n[2]->m_x-impact) +cross_2(m_faces[i].m_n[0]->m_x-impact,m_faces[i].m_n[2]->m_x-impact)
				- cross_2(m_faces[i].m_n[1]->m_x-m_faces[i].m_n[0]->m_x,m_faces[i].m_n[2]->m_x-m_faces[i].m_n[0]->m_x))) <= 0.0001)
			{		
				facenum[f] = i; // 點擊到的face的編號

				face[f][0] = m_faces[i].m_n[0]-nbase;   
				face[f][1] = m_faces[i].m_n[1]-nbase;
				face[f][2] = m_faces[i].m_n[2]-nbase;

				f++;
			}		
		}

		// 如果是對摺一次的情況(f=2) 就要記錄正確的last_face(延續上一次影響點在的face, 做連續剪取的形變) 否則會發生剪取錯誤

		// 有幾塊是要被切割的face,就做幾次 (一般情況只跑一次) (對摺=>不只一次)
		for(int fi=0;fi<f;fi++)
		{
			printf("last_face:%d\n",last_face[fi]);
			if(last_face[fi] != -1)
				printf("Face:%d %d %d\n",m_faces[last_face[fi]].m_n[0]-nbase,m_faces[last_face[fi]].m_n[1]-nbase,m_faces[last_face[fi]].m_n[2]-nbase);

			printf("%d*cutface = %d %d %d\n",facenum[fi],face[fi][0],face[fi][1],face[fi][2]);


			// 將cut到的face 三個edge做tag
			edges(face[fi][0],face[fi][1]) = -3;
			edges(face[fi][1],face[fi][2]) = -3;
			edges(face[fi][2],face[fi][0]) = -3;


			int cut_state;						 // 判斷是第一刀 or 非第一刀
			int up_face = -1;					 // 往上剪是剪左半邊還是右半邊的face (0=左 , 1=右)
			int left_node = -1,right_node = -1;	 // 記錄左右兩端點
			int crossface = 0;					 // 記錄是否crossface了
			int cutedge = -1;					 // 記錄是不是剪到底邊, 要進行剪斷的refine face
			btVector3 intersect;                 // 記錄crossface時會有的交點
			int out_cross=0;


			// 求impact(目前點擊的位置)與last_impact(上次點擊的位置),這兩點連線與cutface(目前點擊到的face)的相交點求出來
			
			// 要知道剪取的裂口是該Face上的哪條link，才能決定相對的左右node

			if(last_face[fi] != -1)  // 不為第一刀的狀態
			{
				// 防止bullet容易因為快速點擊, 且點擊的位置相同或是太接近而發生記憶體錯誤
				if(((impact-m_nodes[info[last_face[fi]].Nimpact].m_x).length2()) < 0.001)
					return;

				printf("last_impact = %d\n",info[last_face[fi]].Nimpact);
				printf("last_rcut = %d\n",info[last_face[fi]].Nrcut);
				printf("last_lcut = %d\n",info[last_face[fi]].Nlcut);

				/////////////////////////////////// 避免cross1塊face以上的情形 /////////////////////////////////////
				/* 查看這次點擊的位置與上次點擊的位置, 其連線有沒有橫越兩塊以上的Face,  */
				
				int linkab[50][2];
				int abn=0;

				for(i=0;i<m_faces.size();i++)
				{
					// Face的三個node兩兩配對跑迴圈
					for(j=2,k=0;k<3;j=k++)
					{						
						btVector3 u1,u2,v1,v2;

						//線段1：此face上的某條link兩端點
						u1 = m_faces[i].m_n[j]->m_x;
						u2 = m_faces[i].m_n[k]->m_x;
						
						//線段2：impact與last_impact
						v1 = impact;
						v2 = m_nodes[info[last_face[fi]].Nimpact].m_x;
						
						//此兩線段沒有交點的話
						if (parallel(u1,u2,v1,v2) || !intersect_in(u1,u2,v1,v2)){
							//printf("link %d %d 無交點\n",m_faces[i].m_n[j]-nbase,m_faces[i].m_n[k]-nbase);
						}
						//此兩線段有交點的話
						else
						{
							//printf("%d\n",intersect_in(u1,u2,v1,v2));
							//printf("link %d %d 有交點\n",m_faces[i].m_n[j]-nbase,m_faces[i].m_n[k]-nbase);

							if(f==1)  // 無對摺的情況
							{
								int in;

								/* 這邊做得處理就是把誤判進來的link剔除, 只取被跨的Face上的link做count => count的數量即是跨了幾塊face */ 
								int same=0;
								for(in=0; in<abn; in++)
								{
									if((linkab[in][0] == m_faces[i].m_n[j]-nbase && linkab[in][1] == m_faces[i].m_n[k]-nbase) 
										|| (linkab[in][0] == m_faces[i].m_n[k]-nbase && linkab[in][1] == m_faces[i].m_n[j]-nbase))
										same = 1;
								}

								if(m_faces[i].m_n[j]-nbase == info[last_face[fi]].Nimpact || m_faces[i].m_n[k]-nbase == info[last_face[fi]].Nimpact);				
								else if(same == 1);
								else
								{

									printf("*link %d %d 有交點\n",m_faces[i].m_n[j]-nbase,m_faces[i].m_n[k]-nbase);
									crossface++;
									linkab[abn][0] = m_faces[i].m_n[j]-nbase;
									linkab[abn][1] = m_faces[i].m_n[k]-nbase;
									abn++;
									intersect = intersection(u1,u2,v1,v2);

								}
							}
							else  // 對折的情況
							{
								if(fi==0)
								{
									one = info[last_face[0]].Nimpact;
									another = info[last_face[1]].Nimpact;
								}
								
								int in;
								int same=0;
								for(in=0; in<abn; in++)
								{
									if((linkab[in][0] == m_faces[i].m_n[j]-nbase && linkab[in][1] == m_faces[i].m_n[k]-nbase) 
										|| (linkab[in][0] == m_faces[i].m_n[k]-nbase && linkab[in][1] == m_faces[i].m_n[j]-nbase))
										same = 1;
								}

								if(m_faces[i].m_n[j]-nbase == one || m_faces[i].m_n[k]-nbase == one);
								else if(m_faces[i].m_n[j]-nbase == another || m_faces[i].m_n[k]-nbase == another);
								else if(same == 1);
								else
								{
									if(edges(m_faces[i].m_n[j]-nbase,m_faces[i].m_n[k]-nbase) == -3)
									{
										printf("*link %d %d 有交點\n",m_faces[i].m_n[j]-nbase,m_faces[i].m_n[k]-nbase);									

										crossface++;
										linkab[abn][0] = m_faces[i].m_n[j]-nbase;
										linkab[abn][1] = m_faces[i].m_n[k]-nbase;
										abn++;
										intersect = intersection(u1,u2,v1,v2);
									}
									out_cross++;
								}
							}
						}
					}
				}
				
				printf("crossface = %d\n outcross = %d\n",crossface,out_cross);

				////////////////////////// 判斷是否有跨face，以及tag要產生裂口的link /////////////////////////////////////
				if((out_cross+1)/4 > 1)
				{
					printf("Sorry~請重新再剪一次>口<\n");
					strcat( outputInfo, "Sorry~ Please try again. >O<\n" );
					out_cross = 0;
					//crossface = 0;
					return;
				}
				// 跨1塊face
				else if(crossface == 1)
				{			
					printf("//// crossface ////\n");

					int in;
					for(in=0;in<abn;in++)
					{
						if(edges(int(linkab[in][0]),int(linkab[in][1])) == -3)
						{
							printf("*%d %d\n",linkab[in][0],linkab[in][1]);

							// 把要產生裂口的link做標記, 以便之後能正確分割Face
							edges(int(linkab[in][0]),int(linkab[in][1])) = -4; 
							break;
						}
					}
					// 把impact和last_impact這條線段, 和要產生裂口的link的交點拉出來記錄
					// 多此一舉的原因是不知為何在剪取物體的時後, intersect的值有時會自己跑掉變成亂數, 導致程式判斷就會錯誤
					store_inter.setX(intersect.getX());
					store_inter.setY(intersect.getY());
					store_inter.setZ(intersect.getZ());	
				}
				// 無跨face
				else if(crossface == 0)
				{
					printf("//// not crossface ////\n");

					for(j=2,k=0;k<3;j=k++)
					{

						// last_impact和last_up這條直線為判斷線, 目前點的位置, 是在判斷線的右邊還左邊, 根據此可以知道cutlink到底是哪條
						btVector3 ab = getab(m_nodes[info[last_face[fi]].Nimpact].m_x,m_nodes[info[last_face[fi]].Nup].m_x);
						printf("judgelink m = %f\n impact = %f\n",ab.getX(),(impact.getY()-ab.getX()*impact.getX()-ab.getY()));

						// 要再多加判斷, 此link端點是last_right或是last_left
						if(start_dir == 0 || start_dir == 3)  // cut from down || up
						{
							
							if((m_faces[facenum[fi]].m_n[j]-nbase == info[last_face[fi]].Nimpact && m_faces[facenum[fi]].m_n[k]-nbase == info[last_face[fi]].Nleft) 
								|| (m_faces[facenum[fi]].m_n[j]-nbase == info[last_face[fi]].Nleft && m_faces[facenum[fi]].m_n[k]-nbase == info[last_face[fi]].Nimpact))
								edges(int(m_faces[facenum[fi]].m_n[j]-nbase),int(m_faces[facenum[fi]].m_n[k]-nbase)) = -4;  // 將cutlink標記成-4

							else if((m_faces[facenum[fi]].m_n[j]-nbase == info[last_face[fi]].Nimpact && m_faces[facenum[fi]].m_n[k]-nbase == info[last_face[fi]].Nright) 
								|| (m_faces[facenum[fi]].m_n[j]-nbase == info[last_face[fi]].Nright && m_faces[facenum[fi]].m_n[k]-nbase == info[last_face[fi]].Nimpact))
								edges(int(m_faces[facenum[fi]].m_n[j]-nbase),int(m_faces[facenum[fi]].m_n[k]-nbase)) = -4;  // 將cutlink標記成-4	
							
						}
						else if(start_dir == 1 || start_dir == 2)  // cut from right || left
						{
							if((ab.getX() < 0 && (impact.getY()-ab.getX()*impact.getX()-ab.getY()) < 0)||
								(ab.getX() > 0 && (impact.getY()-ab.getX()*impact.getX()-ab.getY()) < 0) )
							{
								printf("impact on the left face\n");
								if((m_faces[facenum[fi]].m_n[j]-nbase == info[last_face[fi]].Nimpact && m_faces[facenum[fi]].m_n[k]-nbase == info[last_face[fi]].Nleft) 
									|| (m_faces[facenum[fi]].m_n[j]-nbase == info[last_face[fi]].Nleft && m_faces[facenum[fi]].m_n[k]-nbase == info[last_face[fi]].Nimpact))
									edges(int(m_faces[facenum[fi]].m_n[j]-nbase),int(m_faces[facenum[fi]].m_n[k]-nbase)) = -4;  // 將cutlink標記成-4

							}
							else if((ab.getX() < 0 && (impact.getY()-ab.getX()*impact.getX()-ab.getY()) > 0)||
								(ab.getX() > 0 && (impact.getY()-ab.getX()*impact.getX()-ab.getY()) > 0))
							{
								printf("impact on the right face\n");
								if((m_faces[facenum[fi]].m_n[j]-nbase == info[last_face[fi]].Nimpact && m_faces[facenum[fi]].m_n[k]-nbase == info[last_face[fi]].Nright) 
									|| (m_faces[facenum[fi]].m_n[j]-nbase == info[last_face[fi]].Nright && m_faces[facenum[fi]].m_n[k]-nbase == info[last_face[fi]].Nimpact))
									edges(int(m_faces[facenum[fi]].m_n[j]-nbase),int(m_faces[facenum[fi]].m_n[k]-nbase)) = -4;  // 將cutlink標記成-4	
							}
						}

					}

				}


				///////////////////////////////////////////////////////////////////////////////////////////////


				/*******************  不為第一刀的情形下, 接近底邊時要判斷剪斷哪個底邊  *************************/
				int end_dir;

				// 靠近上邊
				if(fabs(impact.getY()-corner2.getY()) <= allow)
				{
					impact.setY(corner2.getY());
					cutedge = 1;
					end_dir = 0;
				}
				// 靠近右邊  (bullet的右邊)
				else if(fabs(impact.getX()-corner1.getX()) <= allow)
				{
					impact.setX(corner1.getX());
					cutedge = 1;
					end_dir = 1;
				}
				// 靠近左邊  (bullet的左邊)
				else if(fabs(impact.getX()-corner0.getX()) <= allow)
				{
					impact.setX(corner0.getX());
					cutedge = 1;
					end_dir = 2;
				}
				// 靠近下邊 
				else if(fabs(impact.getY()-corner0.getY()) <= allow)
				{
					impact.setY(corner0.getY());
					cutedge = 1;
					end_dir = 3;
				}

				/***************** 對摺情況下, 以摺線當終點 ******************/

				if((state == 1 && end_dir == 3)
					||(state == 3 && end_dir == 2)
					||(state == 5 && end_dir == 1)
					||(state == 7 && end_dir == 0))
					end_foldline = 1;
				/**************************************************************/

			}

			if(f > 1 && cutedge == 1 && crossface == 1)
			{
				return;
			}
			if(f == 1 && cutedge == 1)
			{
				cutend=1;
				printf("cutend\n");
			}
			/*************************************************************************************************/

			// 要知道點擊到的這塊face要造成缺口的link是哪條 才能判定right node & left node
			for(j=2,k=0;k<3;j=k++)
			{
				// edges被標記成-3的兩個點, 為cutface上的edge

				int I = m_faces[facenum[fi]].m_n[j]-nbase;
				int J = m_faces[facenum[fi]].m_n[k]-nbase;

				if(edges(I,J) == -3)
				{

					Node&			a = m_nodes[I];
					Node&			b = m_nodes[J];


					// 第一刀 , cutface中的兩個點在impact node底下

					// cut from down
					if(last_face[fi] == -1 && start_dir == 0 && (a.m_x.getY()-impact.getY())<0 && (b.m_x.getY()-impact.getY())<0)
					{

						if(a.m_x.getX()-impact.getX() > 0) // a點在impact右邊
						{

							right_node = I;
							left_node = J;
						}
						else if(a.m_x.getX()-impact.getX() < 0) // idx[0]點在impact點左邊
						{

							right_node = J;
							left_node = I;

						}
						printf("right_node:%d\nleft_node:%d\n",right_node,left_node);					
					}
					// cut from right 
					else if(last_face[fi] == -1 && start_dir == 1 && (a.m_x.getX()-impact.getX())>0 && (b.m_x.getX()-impact.getX())>0)
					{

						if(a.m_x.getY()-impact.getY() > 0) // a點在impact右邊
						{

							right_node = I;
							left_node = J;
						}
						else if(a.m_x.getY()-impact.getY() < 0) // idx[0]點在impact點左邊
						{

							right_node = J;
							left_node = I;

						}
						printf("right_node:%d\nleft_node:%d\n",right_node,left_node);					
					}
					// cut from left
					else if(last_face[fi] == -1 && start_dir == 2 && (a.m_x.getX()-impact.getX())<0 && (b.m_x.getX()-impact.getX())<0)
					{

						if(a.m_x.getY()-impact.getY() > 0) // a點在impact右邊
						{

							right_node = I;
							left_node = J;
						}
						else if(a.m_x.getY()-impact.getY() < 0) // idx[0]點在impact點左邊
						{

							right_node = J;
							left_node = I;

						}
						printf("right_node:%d\nleft_node:%d\n",right_node,left_node);					
					}
					// cut from up
					else if(last_face[fi] == -1 && start_dir == 3 && (a.m_x.getY()-impact.getY())>0 && (b.m_x.getY()-impact.getY())>0)
					{
						if(a.m_x.getX()-impact.getX() > 0) // a點在impact右邊
						{

							right_node = I;
							left_node = J;
						}
						else if(a.m_x.getX()-impact.getX() < 0) // idx[0]點在impact點左邊
						{

							right_node = J;
							left_node = I;

						}
						printf("right_node:%d\nleft_node:%d\n",right_node,left_node);	
					}
				}

				// 不為第一刀
				if(edges(I,J) == -4)
				{
					// edges被標記成-4的兩個點, 為cutedge
					printf("I J %d %d\n",I,J);

					Node&			a = m_nodes[I];
					Node&			b = m_nodes[J];

					// last_impact和impact這條直線為判斷線, 判斷a點和b點在這條線的右邊還左邊
					btVector3 ab = getab(m_nodes[info[last_face[fi]].Nimpact].m_x,impact);
					printf("judgelink m = %f\n",ab.getX());

					printf("a = %f\n",(a.m_x.getY()-ab.getX()*a.m_x.getX()-ab.getY()));
					printf("b = %f\n",(b.m_x.getY()-ab.getX()*b.m_x.getX()-ab.getY()));
					if(I == info[last_face[fi]].Nrcut || I == info[last_face[fi]].Nlcut || J == info[last_face[fi]].Nrcut || J == info[last_face[fi]].Nlcut); // 往回剪的情況
					else
					{
						// 跨face的情況 判斷左右端點
						if(crossface == 1)
						{

							if(start_dir == 0 || start_dir == 3 || start_dir == 1 || start_dir == 2)  // cut from down || up
							{

								if(ab.getX() > 0 && ab.getX() < 1)  
								{
									if(a.m_x.getY()-ab.getX()*a.m_x.getX()-ab.getY() > 0)
									{	
										right_node = I;
										left_node = J;
									}
									else
									{
										right_node = J;
										left_node = I;

									}
								}
								else if((ab.getX() < 0 && (a.m_x.getY()-ab.getX()*a.m_x.getX()-ab.getY()) > 0)||
									(ab.getX() > 0 && (a.m_x.getY()-ab.getX()*a.m_x.getX()-ab.getY()) < 0) )
								{

									right_node = I;
									left_node = J;

								}
								else if((ab.getX() < 0 && (a.m_x.getY()-ab.getX()*a.m_x.getX()-ab.getY()) < 0)||
									(ab.getX() > 0 && (a.m_x.getY()-ab.getX()*a.m_x.getX()-ab.getY()) > 0))
								{										
									right_node = J;
									left_node = I;
								}
								printf("right_node:%d\nleft_node:%d\n",right_node,left_node);
							}
						}
						// 無跨face的情況 判斷左右端點
						else
						{
							if(start_dir == 0 || start_dir == 3) // cut from down
							{
								if((ab.getX() < 0 && (a.m_x.getY()-ab.getX()*a.m_x.getX()-ab.getY()) > 0) ||
									(ab.getX() < 0 && (b.m_x.getY()-ab.getX()*b.m_x.getX()-ab.getY()) < 0))
								{

									right_node = I;
									left_node = J;

								}
								else if((ab.getX() > 0 && (a.m_x.getY()-ab.getX()*a.m_x.getX()-ab.getY()) <= 0) ||
									(ab.getX() > 0 && (b.m_x.getY()-ab.getX()*b.m_x.getX()-ab.getY()) > 0))
								{
									right_node = I;
									left_node = J;
								}

							}
							if(start_dir == 1 || start_dir == 2) // cut from right or left
							{
								if(a.m_x.getY()-b.m_x.getY() > 0) // a點在b右邊
								{					
									right_node = I;
									left_node = J;
								}
								else if(a.m_x.getY()-b.m_x.getY() < 0) // a點在b左邊
								{
									right_node = J;
									left_node = I;
								}
							}
							printf("right_node:%d\nleft_node:%d\n",right_node,left_node);	
						}
					}
				}

			}

			if(right_node == -1 && left_node == -1)
			{
				printf("哀呀~這個不行剪\n");
				strcat( outputInfo, "Oops~ Please try again. ^_^\n" );
				return;

			}

			//////////////////////////////////// 新增impact_node ///////////////////////////////////////////
			int impact_node;
			if(end_foldline == 1 && fi == 1)
			{
				impact_node = info[last_face[fi-1]].Nimpact;
				printf("*impact node %d : %f %f %f\n",impact_node,m_nodes[impact_node].m_x.getX(),m_nodes[impact_node].m_x.getY(),m_nodes[impact_node].m_x.getZ());
			}
			else if(end_foldline == 0 || (end_foldline == 1 && fi==0))
			{
				// 計算新的mass值
				btScalar impact_im, impact_m;

				for(i=0;i<m_nodes.size();i++)
				{

					if(m_nodes[i].m_im != 0)
					{
						impact_im = m_nodes[i].m_im / 0.5f;
						break;
					}
				}
				//impact_im = m_nodes[0].m_im / 0.5f;
				impact_m = 1/impact_im;

				//計算新的m_v
				const btVector3	impact_v = m_nodes[0].m_v;
				// 新增impact_node
				appendNode(impact,impact_m);
				impact_node = m_nodes.size()-1;
				printf("impact node %d : %f %f %f\n",impact_node,m_nodes[impact_node].m_x.getX(),m_nodes[impact_node].m_x.getY(),m_nodes[impact_node].m_x.getZ());
				m_nodes[m_nodes.size()-1].m_v = V;
				++newnodes;
				nbase=&m_nodes[0]; // 每新增一個node 要記得重指定一次nbase的位置

				m_nodes[impact_node].route = 1;
				impact_x = m_nodes[impact_node].m_x.getX();
				impact_y = m_nodes[impact_node].m_x.getY();
				impact_z = m_nodes[impact_node].m_x.getZ();
			}
			///////////////////////////////////////////////////////////////////////////////////////////

			float gap = 0.08;
			int up_node;	
			int rcut_node,lcut_node;
			int edge_rcut,edge_lcut;

			if(crossface == 1)
				printf("store intersect : %f %f %f\n",store_inter.getX(),store_inter.getY(),store_inter.getZ());

		
			btScalar t,m;
			btVector3 x,v; 


			for(i=0;i<3;i++)
			{
				// 找除了左右兩點的那一點 (upnode)
				if((face[fi][i]!= right_node) && (face[fi][i] != left_node))
				{

					up_node = face[fi][i];
					printf("up_node = %d\n",up_node);


					if(end_foldline == 1 && fi==1) // 如果在對摺情況下, 終點在摺線上的話, 第二塊被影響的face, 其邊界點要跟第一塊一樣
					{
						rcut_node = info[last_face[fi-1]].Nrcut;
						lcut_node = info[last_face[fi-1]].Nlcut;

						printf("****end fold line\n");
					}


					if(end_foldline == 0 || (end_foldline == 1 && fi == 0))
					{
						////////////////////////////////////// rcut_node ///////////////////////////////////////////////
						t = Solve(m_nodes[impact_node].m_x, m_nodes[impact_node].m_x ,m_nodes[right_node].m_x, accurary);

						if(t>0)
						{
							x = Lerp(m_nodes[impact_node].m_x,m_nodes[right_node].m_x,t);  // 取得範圍在此link上的交界點
							v = m_nodes[impact_node].m_v;
							m = cal_m(m_nodes[impact_node],m_nodes[right_node],t,m_nodes[impact_node]);

							appendNode(x,m);
							rcut_node = m_nodes.size()-1;
							m_nodes[m_nodes.size()-1].m_v=V;
							++newnodes;
							nbase=&m_nodes[0]; // 每新增一個node 要記得重指定一次nbase的位置
							printf("rcut node %d : %f %f %f\n",rcut_node,m_nodes[rcut_node].m_x.getX(),m_nodes[rcut_node].m_x.getY(),m_nodes[rcut_node].m_x.getZ());

						}


						////////////////////////////////////// lcut_node ///////////////////////////////////////////////

						t = Solve(m_nodes[impact_node].m_x,m_nodes[impact_node].m_x ,m_nodes[left_node].m_x, accurary);
						if(t>0)
						{
							x = Lerp(m_nodes[impact_node].m_x,m_nodes[left_node].m_x,t); 		
							v = m_nodes[rcut_node].m_v;			
							m = cal_m(m_nodes[impact_node],m_nodes[left_node],t,m_nodes[impact_node]);

							appendNode(x,m);
							lcut_node = m_nodes.size()-1;
							m_nodes[m_nodes.size()-1].m_v=V;
							++newnodes;
							nbase=&m_nodes[0]; // 每新增一個node 要記得重指定一次nbase的位置
							printf("lcut node %d : %f %f %f\n",lcut_node,m_nodes[lcut_node].m_x.getX(),m_nodes[lcut_node].m_x.getY(),m_nodes[lcut_node].m_x.getZ());

						}
					}
					// 判斷目前剪取是哪種狀況，根據該狀況新增edge node
					if(last_face[fi] == -1) // 第一刀
					{
						cut_state = 0;			


						////////////////////////////////////// edge_rcut ///////////////////////////////////////////////

						if(cut_foldline == 1 && fi==1) // 如果在對摺情況下, 起始點就在摺線上的話, 第二塊被影響的face, 其邊界點要跟第一塊一樣
						{
							edge_rcut = info[last_face[fi-1]].Nrcut+2;
							edge_lcut = info[last_face[fi-1]].Nlcut+2;
							break;
						}

						t = Solve(start, start,m_nodes[right_node].m_x, accurary);
						if(t>0)
						{
							x = Lerp(start,m_nodes[right_node].m_x,t);  
							v = m_nodes[rcut_node].m_v;
							m = cal_m(m_nodes[impact_node],m_nodes[right_node],t,m_nodes[impact_node]);

							appendNode(x,m);
							edge_rcut = m_nodes.size()-1;
							m_nodes[m_nodes.size()-1].m_v=V;
							++newnodes;
							nbase=&m_nodes[0]; // 每新增一個node 要記得重指定一次nbase的位置
							//printf("edge_rcut %d : %f %f %f\n",edge_rcut,m_nodes[edge_rcut].m_x.getX(),m_nodes[edge_rcut].m_x.getY(),m_nodes[edge_rcut].m_x.getZ());

						}

						////////////////////////////////////// edge_lcut ///////////////////////////////////////////////
						t = Solve(start, start,m_nodes[left_node].m_x, accurary);
						if(t>0)
						{
							x = Lerp(start,m_nodes[left_node].m_x,t); 						
							v = m_nodes[rcut_node].m_v;
							m = cal_m(m_nodes[impact_node],m_nodes[left_node],t,m_nodes[impact_node]);

							appendNode(x,m);
							edge_lcut = m_nodes.size()-1;
							m_nodes[m_nodes.size()-1].m_v=V;
							++newnodes;
							nbase=&m_nodes[0]; // 每新增一個node 要記得重指定一次nbase的位置
							//printf("edge_lcut %d : %f %f %f\n",edge_lcut,m_nodes[edge_lcut].m_x.getX(),m_nodes[edge_lcut].m_x.getY(),m_nodes[edge_lcut].m_x.getZ());

						}

						break;
					}
					// 不為第一刀
					else 
					{

						cut_state = 1;

						// 判斷往上剪的情況, 是往左邊的face, 還是右邊的face (0=左, 1=右)
						if(crossface == 0)
						{
							// 目前點擊的face的right_node為上次點擊的face的impact,代表往左邊的face剪
							if(right_node == info[last_face[fi]].Nimpact)
								up_face = 0;
							else if(left_node == info[last_face[fi]].Nimpact)
								up_face = 1;
						}
						printf("up_face = %d\n",up_face);


						if(crossface == 1)
						{
							///////////////////////////// edge_rcut ////////////////////////////
							t = Solve(store_inter, store_inter,m_nodes[right_node].m_x, accurary);
							printf("crossface t = %f\n",t);
							if(t>0)
							{
								x = Lerp(store_inter,m_nodes[right_node].m_x,t);  
								v = m_nodes[rcut_node].m_v;
								m = cal_m(m_nodes[impact_node],m_nodes[right_node],t,m_nodes[impact_node]);

								appendNode(x,m);
								edge_rcut = m_nodes.size()-1;
								m_nodes[m_nodes.size()-1].m_v=V;
								++newnodes;
								nbase=&m_nodes[0]; // 每新增一個node 要記得重指定一次nbase的位置
								printf("edge_rcut %d : %f %f %f\n",edge_rcut,m_nodes[edge_rcut].m_x.getX(),m_nodes[edge_rcut].m_x.getY(),m_nodes[edge_rcut].m_x.getZ());

							}

							///////////////////////////// edge_lcut ////////////////////////////
							t = Solve(store_inter, store_inter,m_nodes[left_node].m_x, accurary);
							if(t>0)
							{
								x = Lerp(store_inter,m_nodes[left_node].m_x,t);  
								v = m_nodes[rcut_node].m_v;
								m = cal_m(m_nodes[impact_node],m_nodes[left_node],t,m_nodes[impact_node]);

								appendNode(x,m);
								edge_lcut = m_nodes.size()-1;
								m_nodes[m_nodes.size()-1].m_v=V;
								++newnodes;
								nbase=&m_nodes[0]; // 每新增一個node 要記得重指定一次nbase的位置
								printf("edge_lcut %d : %f %f %f\n",edge_lcut,m_nodes[edge_lcut].m_x.getX(),m_nodes[edge_lcut].m_x.getY(),m_nodes[edge_lcut].m_x.getZ());

							}
						}
						break;
					}
				}
			}

		
			if(cut_state == 0)  // 第一刀
			{

				/******************* 分解cutface ********************/
				for(i=0;i<m_faces.size();i++)
				{

					const Face&	feat=m_faces[i];
					const int	idx[]={	int(feat.m_n[0]-nbase),
						int(feat.m_n[1]-nbase),
						int(feat.m_n[2]-nbase)};

					// 找到cutface是陣列中第幾塊face (i) , 才能call appendface()
					if((idx[0] == face[fi][0]) && (idx[1] == face[fi][1]) && (idx[2] == face[fi][2]))
					{

						appendFace(i);

						Face*		pft[]={	&m_faces[i],
							&m_faces[m_faces.size()-1]};

				
						pft[0]->m_n[0]=&m_nodes[left_node];
						pft[0]->m_n[1]=&m_nodes[up_node];
						pft[0]->m_n[2]=&m_nodes[impact_node];
						//printf("new face pft[0] : %d %d %d\n",left_node,up_node,impact_node);
				
						pft[1]->m_n[0]=&m_nodes[impact_node];
						pft[1]->m_n[1]=&m_nodes[up_node];
						pft[1]->m_n[2]=&m_nodes[right_node];
						//printf("new face pft[1] : %d %d %d\n",impact_node,up_node,right_node);


						appendFace(i);

						Face*		pftt[]={	&m_faces[i]};

				
						pftt[0]->m_n[0]=&m_nodes[left_node];
						pftt[0]->m_n[1]=&m_nodes[impact_node];
						pftt[0]->m_n[2]=&m_nodes[right_node];
						//printf("new face pftt[0] : %d %d %d\n",left_node,impact_node,right_node);


						appendLink(impact_node,up_node);  // 4 2
						appendLink(impact_node,left_node);  // 4 0
						appendLink(impact_node,right_node);  // 4 1 



						appendFace(i);

						Face*		pfttt[]={	&m_faces[i],
							&m_faces[m_faces.size()-1]};

					
						pfttt[0]->m_n[0]=&m_nodes[left_node];
						pfttt[0]->m_n[1]=&m_nodes[lcut_node];
						pfttt[0]->m_n[2]=&m_nodes[edge_lcut];
						//printf("new face pfttt[0] : %d %d %d\n",left_node,lcut_node,edge_lcut);
					
						pfttt[1]->m_n[0]=&m_nodes[rcut_node];
						pfttt[1]->m_n[1]=&m_nodes[edge_rcut];
						pfttt[1]->m_n[2]=&m_nodes[right_node];
						//printf("new face pfttt[1] : %d %d %d\n",rcut_node,edge_rcut,right_node);


						appendLink(lcut_node,edge_lcut);  // 6 8
						appendLink(rcut_node,edge_rcut);  // 5 7


						/********* 重新連接被新增node的link  ***********/

						for(j=0,ni=m_links.size();j<ni;++j)
						{
							Link&		feat=m_links[j];
							const int	idx[]={	int(feat.m_n[0]-nbase),
								int(feat.m_n[1]-nbase)};

							if((idx[0] == impact_node && idx[1] == right_node) || (idx[1] == impact_node && idx[0] == right_node))
							{

								appendLink(j);
								// 重新連接link 
								Link*		pft[]={	&m_links[j],
									&m_links[m_links.size()-1]};
							
								pft[0]->m_n[0]=&m_nodes[idx[0]]; 
								pft[0]->m_n[1]=&m_nodes[rcut_node]; 
								//printf("new link pft[0] : %d %d\nlennth:%f\n",idx[0],rcut_node,(impact-m_nodes[rcut_node].m_x).length());

							
								pft[1]->m_n[0]=&m_nodes[rcut_node];   
								pft[1]->m_n[1]=&m_nodes[idx[1]];    
								//printf("new link pft[1] : %d %d\n",rcut_node,idx[1]);


							}
							if((idx[0] == impact_node && idx[1] == left_node) || (idx[1] == impact_node && idx[0] == left_node))
							{

								appendLink(j);
								// 重新連接link 
								Link*		pft[]={	&m_links[j],
									&m_links[m_links.size()-1]};
								
								pft[0]->m_n[0]=&m_nodes[idx[0]]; 
								pft[0]->m_n[1]=&m_nodes[lcut_node]; 
								//printf("new link pft[0] : %d %d\nlength:%f\n",idx[0],lcut_node,(impact-m_nodes[lcut_node].m_x).length());
							
								pft[1]->m_n[0]=&m_nodes[lcut_node];   
								pft[1]->m_n[1]=&m_nodes[idx[1]];    
								//printf("new link pft[1] : %d %d\n",lcut_node,idx[1]);

							}

						}
					}
				}



				// 把底部link pop掉
				for(i=0;i<m_links.size();i++)
				{
					if( ((m_links[i].m_n[0]-nbase) == left_node && (m_links[i].m_n[1]-nbase) == right_node) || 
						((m_links[i].m_n[1]-nbase) == left_node && (m_links[i].m_n[0]-nbase) == right_node) )
					{
						btSwap(m_links[i], m_links[m_links.size()-1]);  // 單純做template型態的swap, 將m_link第i個陣列的東西, 跟最後一個陣列的東西交換
						m_links.pop_back(); // pop_back() 就是將m_size--, 並做destructor
						--i; 
					}

				}


				// 重新連接link
				appendLink(right_node,edge_rcut); // 1 6
				appendLink(edge_lcut,left_node); // 8 0

				//if(f == 1)
				//{
				// 記錄當前impact,rcut,lcut編號

				info[facenum[fi]].Nimpact = impact_node;
				info[facenum[fi]].Nrcut = rcut_node;
				info[facenum[fi]].Nlcut = lcut_node;
				info[facenum[fi]].Nright = right_node;
				info[facenum[fi]].Nleft = left_node;
				info[facenum[fi]].Nup = up_node;

				last_face[fi] = facenum[fi];
				//}


			}
			else if(cut_state == 1)  // 往上剪
			{

				printf("cut_state = 1\n");

				/***************************** 重連被cross過的face ***********************************/
				int poplr = -1; // 記錄到底要pop掉last_rcut還是last_lcut
				if(crossface == 1)
				{
					for(i=0;i<m_faces.size();i++)
					{

						const Face&	feat=m_faces[i];
						const int	idx[]={	int(feat.m_n[0]-nbase),
							int(feat.m_n[1]-nbase),
							int(feat.m_n[2]-nbase)};


						for(j=2,k=0;k<3;j=k++)
						{
							const int	l=(k+1)%3;
							// cross的第一塊face
							if(edges(idx[j],idx[k]) == -4 && (idx[l] == info[last_face[fi]].Nimpact))   
							{
								printf("cross on %d face %d %d %d\n",i,idx[0],idx[1],idx[2]);
								printf("cutdege = %d %d\n",idx[j],idx[k]);


								if(edges(idx[j],info[last_face[fi]].Nlcut) == -1 || edges(idx[k],info[last_face[fi]].Nlcut) == -1)		
								{
									poplr = 0;

									printf("poplr = %d (0:last_lcut,1:last_rcut)\n",poplr);
									appendFace(i);

									Face*		pftttt[]={	&m_faces[i],
										&m_faces[m_faces.size()-1]};

									if(intersect_in(m_nodes[info[last_face[fi]].Nlcut].m_x,m_nodes[edge_lcut].m_x,m_nodes[info[last_face[fi]].Nimpact].m_x,m_nodes[edge_rcut].m_x) == 1)
									{
										pftttt[0]->m_n[0]=&m_nodes[idx[j]];
										pftttt[0]->m_n[1]=&m_nodes[info[last_face[fi]].Nlcut];
										pftttt[0]->m_n[2]=&m_nodes[edge_rcut];
										printf("new face pftttt[0] : %d %d %d\n",idx[j],info[last_face[fi]].Nlcut,edge_lcut);

										pftttt[1]->m_n[0]=&m_nodes[info[last_face[fi]].Nimpact];
										pftttt[1]->m_n[1]=&m_nodes[edge_lcut];
										pftttt[1]->m_n[2]=&m_nodes[idx[k]];
										printf("new face pftttt[1] : %d %d %d\n",info[last_face[fi]].Nimpact,edge_rcut,idx[k]);

										if(fi == 0 )
										{
											m_nodes[info[last_face[fi]].Nlcut].tag = 1;										
											m_nodes[edge_rcut].tag = 1;
										}
										else if(fi == 1)
										{
											m_nodes[info[last_face[fi]].Nlcut].foldtag = 1;										
											m_nodes[edge_rcut].foldtag = 1;
										}

										appendLink(edge_rcut,info[last_face[fi]].Nlcut);  
										appendLink(edge_lcut,info[last_face[fi]].Nimpact);  
									}
									else
									{
										pftttt[0]->m_n[0]=&m_nodes[idx[j]];
										pftttt[0]->m_n[1]=&m_nodes[info[last_face[fi]].Nlcut];
										pftttt[0]->m_n[2]=&m_nodes[edge_lcut];
										printf("new face pftttt[0] : %d %d %d\n",idx[j],info[last_face[fi]].Nlcut,edge_lcut);

										pftttt[1]->m_n[0]=&m_nodes[info[last_face[fi]].Nimpact];
										pftttt[1]->m_n[1]=&m_nodes[edge_rcut];
										pftttt[1]->m_n[2]=&m_nodes[idx[k]];
										printf("new face pftttt[1] : %d %d %d\n",info[last_face[fi]].Nimpact,edge_rcut,idx[k]);

										if(fi == 0)
										{
											m_nodes[info[last_face[fi]].Nlcut].tag = 1;
											m_nodes[edge_lcut].tag = 1;
										}
										else if(fi == 1)
										{
											m_nodes[info[last_face[fi]].Nlcut].foldtag = 1;
											m_nodes[edge_lcut].foldtag = 1;

										}

										appendLink(edge_lcut,info[last_face[fi]].Nlcut);  
										appendLink(edge_rcut,info[last_face[fi]].Nimpact); 
									}
									break;

								}					
								else if(edges(idx[j],info[last_face[fi]].Nrcut) == -1 || edges(idx[k],info[last_face[fi]].Nrcut) == -1)
								{
									poplr = 1;

									printf("poplr = %d (0:last_lcut,1:last_rcut)\n",poplr);
									appendFace(i);

									Face*		pftttt[]={	&m_faces[i],
										&m_faces[m_faces.size()-1]};


									printf("%d\n",intersect_in(m_nodes[info[last_face[fi]].Nimpact].m_x,m_nodes[edge_lcut].m_x,m_nodes[info[last_face[fi]].Nrcut].m_x,m_nodes[edge_rcut].m_x));
									if(intersect_in(m_nodes[info[last_face[fi]].Nimpact].m_x,m_nodes[edge_lcut].m_x,m_nodes[info[last_face[fi]].Nrcut].m_x,m_nodes[edge_rcut].m_x) == 1)
									{
										pftttt[0]->m_n[0]=&m_nodes[idx[j]];
										pftttt[0]->m_n[1]=&m_nodes[info[last_face[fi]].Nimpact];
										pftttt[0]->m_n[2]=&m_nodes[edge_rcut];
										printf("new face pftttt[0] : %d %d %d\n",idx[j],info[last_face[fi]].Nimpact,edge_rcut);

										pftttt[1]->m_n[0]=&m_nodes[info[last_face[fi]].Nrcut];
										pftttt[1]->m_n[1]=&m_nodes[edge_lcut];
										pftttt[1]->m_n[2]=&m_nodes[idx[k]];
										printf("new face pftttt[1] : %d %d %d\n",info[last_face[fi]].Nrcut,edge_lcut,idx[k]);

										if(fi==0)
										{
											m_nodes[info[last_face[fi]].Nrcut].tag = 1;
											m_nodes[edge_lcut].tag = 1;
										}
										else if(fi==1)
										{
											m_nodes[info[last_face[fi]].Nrcut].foldtag = 1;
											m_nodes[edge_lcut].foldtag = 1;

										}

										appendLink(edge_rcut,info[last_face[fi]].Nimpact);  
										appendLink(edge_lcut,info[last_face[fi]].Nrcut);  

									}
									else
									{
										pftttt[0]->m_n[0]=&m_nodes[idx[j]];
										pftttt[0]->m_n[1]=&m_nodes[info[last_face[fi]].Nimpact];
										pftttt[0]->m_n[2]=&m_nodes[edge_lcut];
										printf("new face pftttt[0] : %d %d %d\n",idx[j],info[last_face[fi]].Nimpact,edge_lcut);

										pftttt[1]->m_n[0]=&m_nodes[info[last_face[fi]].Nrcut];
										pftttt[1]->m_n[1]=&m_nodes[edge_rcut];
										pftttt[1]->m_n[2]=&m_nodes[idx[k]];
										printf("new face pftttt[1] : %d %d %d\n",info[last_face[fi]].Nrcut,edge_rcut,idx[k]);

										if(fi==0)
										{
											m_nodes[info[last_face[fi]].Nrcut].tag = 1;
											m_nodes[edge_rcut].tag = 1;
										}
										else if(fi==1)
										{
											m_nodes[info[last_face[fi]].Nrcut].foldtag = 1;
											m_nodes[edge_rcut].foldtag = 1;

										}

										appendLink(edge_lcut,info[last_face[fi]].Nimpact);  
										appendLink(edge_rcut,info[last_face[fi]].Nrcut); 
									}
									break;
								}


							}
						}		
					}
				}
				/***************************** ****************** ***********************************/

				for(i=0;i<m_faces.size();i++)
				{

					const Face&	feat=m_faces[i];
					const int	idx[]={	int(feat.m_n[0]-nbase),
						int(feat.m_n[1]-nbase),
						int(feat.m_n[2]-nbase)};

					// 找到cutface是陣列中第幾塊face (i) , 才能call appendface()
					if((idx[0] == face[fi][0]) && (idx[1] == face[fi][1]) && (idx[2] == face[fi][2]))
					{

						if(cutedge == -1) // 不為底邊
						{
							appendFace(i);

							Face*		pft[]={	&m_faces[i],
								&m_faces[m_faces.size()-1]};

							pft[0]->m_n[0]=&m_nodes[left_node];
							pft[0]->m_n[1]=&m_nodes[up_node];
							pft[0]->m_n[2]=&m_nodes[impact_node];
							printf("new face pft[0] : %d %d %d\n",left_node,up_node,impact_node);

							pft[1]->m_n[0]=&m_nodes[impact_node];
							pft[1]->m_n[1]=&m_nodes[up_node];
							pft[1]->m_n[2]=&m_nodes[right_node];
							printf("new face pft[1] : %d %d %d\n",impact_node,up_node,right_node);


							appendLink(impact_node,up_node);  
							appendLink(impact_node,left_node);  
							appendLink(impact_node,right_node); 
							///////////////////////////////// not crossface ////////////////////////////////////////////
							if(up_face == 0) // last_lcut
							{
								appendFace(i);

								Face*		pftt[]={	&m_faces[i]};

								// 0 11 7
								pftt[0]->m_n[0]=&m_nodes[left_node];
								pftt[0]->m_n[1]=&m_nodes[lcut_node];
								pftt[0]->m_n[2]=&m_nodes[info[last_face[fi]].Nlcut];
								//printf("new face pftt[0] : %d %d %d\n",left_node,lcut_node,info[last_face[fi]].Nlcut);
								if(fi == 0)
									m_nodes[info[last_face[fi]].Nlcut].tag = 1;
								else if(fi == 1)
									m_nodes[info[last_face[fi]].Nlcut].foldtag = 1;

								appendLink(lcut_node,info[last_face[fi]].Nlcut,pftt[0]->m_material);  // 11 7
							}
							else if(up_face == 1) // last_rcut
							{
								appendFace(i);
								Face*		pftt[]={	&m_faces[i]};

								// 1 10 5
								pftt[0]->m_n[0]=&m_nodes[right_node];
								pftt[0]->m_n[1]=&m_nodes[rcut_node];
								pftt[0]->m_n[2]=&m_nodes[info[last_face[fi]].Nrcut];
								//printf("new face pftt[0] : %d %d %d\n",right_node,rcut_node,info[last_face[fi]].Nrcut);
								if(fi == 0)
									m_nodes[info[last_face[fi]].Nrcut].tag = 1;
								else if(fi == 1)
									m_nodes[info[last_face[fi]].Nrcut].foldtag = 1;
								appendLink(rcut_node,info[last_face[fi]].Nrcut);  // 10 5

							}
							/////////////////////////////////// crossface /////////////////////////////////////////
							if(crossface == 1)
							{

								appendFace(i);

								Face*		pftt[]={	&m_faces[i]};


								pftt[0]->m_n[0]=&m_nodes[left_node];
								pftt[0]->m_n[1]=&m_nodes[impact_node];
								pftt[0]->m_n[2]=&m_nodes[right_node];
								printf("new face pftt[0] : %d %d %d\n",left_node,impact_node,right_node);


								appendFace(i);

								Face*		pfttt[]={	&m_faces[i],
									&m_faces[m_faces.size()-1]};


								pfttt[0]->m_n[0]=&m_nodes[left_node];
								pfttt[0]->m_n[1]=&m_nodes[lcut_node];
								pfttt[0]->m_n[2]=&m_nodes[edge_lcut];
								printf("new face pfttt[0] : %d %d %d\n",left_node,lcut_node,edge_lcut);

								pfttt[1]->m_n[0]=&m_nodes[rcut_node];
								pfttt[1]->m_n[1]=&m_nodes[edge_rcut];
								pfttt[1]->m_n[2]=&m_nodes[right_node];
								printf("new face pfttt[1] : %d %d %d\n",rcut_node,edge_rcut,right_node);

								if(fi == 0)
								{
									if(m_nodes[edge_lcut].tag == 1)
										m_nodes[lcut_node].tag = 1;
								}
								else if(fi == 1)
								{
									if(m_nodes[edge_rcut].foldtag == 1)
										m_nodes[rcut_node].foldtag = 1;
								}
								appendLink(lcut_node,edge_lcut);  
								appendLink(rcut_node,edge_rcut);  

							}

							/********* 重新連接被新增node的link  ***********/

							for(j=0,ni=m_links.size();j<ni;++j)
							{
								Link&		feat=m_links[j];
								const int	idx[]={	int(feat.m_n[0]-nbase),
									int(feat.m_n[1]-nbase)};

								if((idx[0] == impact_node && idx[1] == right_node) || (idx[1] == impact_node && idx[0] == right_node))
								{

									appendLink(j);
									// 重新連接link 
									Link*		pft[]={	&m_links[j],
										&m_links[m_links.size()-1]};

									pft[0]->m_n[0]=&m_nodes[idx[0]]; 
									pft[0]->m_n[1]=&m_nodes[rcut_node]; 
									printf("new link pft[0] : %d %d\n",idx[0],rcut_node);

									pft[1]->m_n[0]=&m_nodes[rcut_node];   
									pft[1]->m_n[1]=&m_nodes[idx[1]];    
									printf("new link pft[1] : %d %d\n",rcut_node,idx[1]);

								}

								if((idx[0] == impact_node && idx[1] == left_node) || (idx[1] == impact_node && idx[0] == left_node))
								{
									appendLink(j);
									// 重新連接link 
									Link*		pft[]={	&m_links[j],
										&m_links[m_links.size()-1]};

									pft[0]->m_n[0]=&m_nodes[idx[0]]; 
									pft[0]->m_n[1]=&m_nodes[lcut_node]; 
									printf("new link pft[0] : %d %d\n",idx[0],lcut_node);

									pft[1]->m_n[0]=&m_nodes[lcut_node];   
									pft[1]->m_n[1]=&m_nodes[idx[1]];    
									printf("new link pft[1] : %d %d\n",lcut_node,idx[1]);


								}	
							}
							printf("m_face.size():%d\n",m_faces.size());

							break;
						}
						// 剪到底邊了, 要準備剪斷的動作
						else if(cutedge == 1) 
						{
							if(crossface == 1)
							{
								printf("crossface cut\n");
								if(poplr == 0)
								{
									appendFace(i);

									Face*		pft[]={	&m_faces[i],&m_faces[m_faces.size()-1]};

									pft[0]->m_n[0]=&m_nodes[left_node];
									pft[0]->m_n[1]=&m_nodes[lcut_node];
									pft[0]->m_n[2]=&m_nodes[right_node];

									pft[1]->m_n[0]=&m_nodes[up_node];
									pft[1]->m_n[1]=&m_nodes[lcut_node];
									pft[1]->m_n[2]=&m_nodes[left_node];
									printf("new face pft[0] : %d %d %d\n",left_node,lcut_node,right_node);
									printf("new face pft[1] : %d %d %d\n",up_node,lcut_node,left_node);
									appendLink(lcut_node,left_node);
								}
								else if(poplr == 1)
								{
									appendFace(i);

									Face*		pft[]={	&m_faces[i],&m_faces[m_faces.size()-1]};

									pft[0]->m_n[0]=&m_nodes[left_node];
									pft[0]->m_n[1]=&m_nodes[rcut_node];
									pft[0]->m_n[2]=&m_nodes[right_node];

									pft[1]->m_n[0]=&m_nodes[up_node];
									pft[1]->m_n[1]=&m_nodes[rcut_node];
									pft[1]->m_n[2]=&m_nodes[right_node];
									printf("new face pft[0] : %d %d %d\n",left_node,rcut_node,right_node);
									printf("new face pft[1] : %d %d %d\n",up_node,rcut_node,right_node);
									appendLink(rcut_node,right_node);
								}
								appendFace(i);

								Face*		pftt[]={	&m_faces[i],
									&m_faces[m_faces.size()-1]};

								pftt[0]->m_n[0]=&m_nodes[left_node];
								pftt[0]->m_n[1]=&m_nodes[lcut_node];
								pftt[0]->m_n[2]=&m_nodes[edge_lcut];


								pftt[1]->m_n[0]=&m_nodes[rcut_node];
								pftt[1]->m_n[1]=&m_nodes[edge_rcut];
								pftt[1]->m_n[2]=&m_nodes[right_node];

								printf("new face pftt[0] : %d %d %d\n",left_node,lcut_node,edge_lcut);
								printf("new face pftt[1] : %d %d %d\n",rcut_node,edge_rcut,right_node);
								appendLink(lcut_node,edge_lcut);  
								appendLink(rcut_node,edge_rcut);  

							}
							if(up_face == 0) // last_lcut
							{
								appendFace(i);

								Face*		pft[]={	&m_faces[i],
									&m_faces[m_faces.size()-1]};


								pft[0]->m_n[0]=&m_nodes[left_node];
								pft[0]->m_n[1]=&m_nodes[lcut_node];
								pft[0]->m_n[2]=&m_nodes[info[last_face[fi]].Nlcut];
								//printf("new face pft[0] : %d %d %d\n",left_node,lcut_node,info[last_face[fi]].Nlcut);

								pft[1]->m_n[0]=&m_nodes[impact_node];
								pft[1]->m_n[1]=&m_nodes[up_node];
								pft[1]->m_n[2]=&m_nodes[right_node];
								//printf("new face pft[1] : %d %d %d\n",impact_node,up_node,right_node);

								if(fi == 0)
									m_nodes[lcut_node].tag = 1;
								else if(fi == 1)
									m_nodes[lcut_node].foldtag = 1;

								appendLink(lcut_node,info[last_face[fi]].Nlcut);  
								appendLink(impact_node,info[last_face[fi]].Nimpact);  

							}
							else if(up_face == 1) // last_rcut
							{
								appendFace(i);

								Face*		pft[]={	&m_faces[i],
									&m_faces[m_faces.size()-1]};


								pft[0]->m_n[0]=&m_nodes[right_node];
								pft[0]->m_n[1]=&m_nodes[rcut_node];
								pft[0]->m_n[2]=&m_nodes[info[last_face[fi]].Nrcut];
								//printf("new face pft[0] : %d %d %d\n",right_node,rcut_node,info[last_face[fi]].Nrcut);

								pft[1]->m_n[0]=&m_nodes[impact_node];
								pft[1]->m_n[1]=&m_nodes[up_node];
								pft[1]->m_n[2]=&m_nodes[left_node];
								//printf("new face pft[1] : %d %d %d\n",impact_node,up_node,left_node);

								if(fi == 0)
									m_nodes[rcut_node].tag = 1;
								else if(fi == 1)
									m_nodes[rcut_node].foldtag = 1;

								appendLink(rcut_node,info[last_face[fi]].Nrcut);  
								appendLink(impact_node,info[last_face[fi]].Nimpact);  


							}

							/********* 重新連接被新增node的link  ***********/
							for(j=0,ni=m_links.size();j<ni;++j)
							{
								Link&		feat=m_links[j];
								const int	idx[]={	int(feat.m_n[0]-nbase),
									int(feat.m_n[1]-nbase)};

								if((idx[0] == info[last_face[fi]].Nimpact && idx[1] == left_node) || (idx[1] == info[last_face[fi]].Nimpact && idx[0] == left_node))
								{
									appendLink(j);

									Link*		pft[]={	&m_links[j],
										&m_links[m_links.size()-1]};

									pft[0]->m_n[0]=&m_nodes[idx[0]]; 
									pft[0]->m_n[1]=&m_nodes[info[last_face[fi]].Nlcut]; 
									//printf("new link pft[0] : %d %d\n",idx[0],info[last_face[fi]].Nlcut);

									pft[1]->m_n[0]=&m_nodes[info[last_face[fi]].Nlcut];   
									pft[1]->m_n[1]=&m_nodes[idx[1]];    

									//printf("new link pft[1] : %d %d\n",info[last_face[fi]].Nlcut,idx[1]);
								}
								if((idx[0] == info[last_face[fi]].Nimpact && idx[1] == right_node) || (idx[1] == info[last_face[fi]].Nimpact && idx[0] == right_node))
								{

									appendLink(j);

									Link*		pft[]={	&m_links[j],
										&m_links[m_links.size()-1]};

									pft[0]->m_n[0]=&m_nodes[idx[0]]; 
									pft[0]->m_n[1]=&m_nodes[info[last_face[fi]].Nrcut]; 
									//printf("new link pft[0] : %d %d\n",idx[0],info[last_face[fi]].Nrcut);

									pft[1]->m_n[0]=&m_nodes[info[last_face[fi]].Nrcut];   
									pft[1]->m_n[1]=&m_nodes[idx[1]];    

									//printf("new link pft[1] : %d %d\n",info[last_face[fi]].Nrcut,idx[1]);
								}

							}
						}
					}
				}



				if(up_face == 0)
				{

					// 把底部靠右的link pop掉
					for(i=0;i<m_links.size();i++)
					{
						if( ((m_links[i].m_n[0]-nbase) == right_node && (m_links[i].m_n[1]-nbase) == info[last_face[fi]].Nlcut) || 
							((m_links[i].m_n[1]-nbase) == right_node && (m_links[i].m_n[0]-nbase) == info[last_face[fi]].Nlcut))
						{
							btSwap(m_links[i], m_links[m_links.size()-1]);  // 單純做template型態的swap, 將m_link第i個陣列的東西, 跟最後一個陣列的東西交換
							m_links.pop_back(); // pop_back() 就是將m_size--, 並做destructor
							--i; 
						}

					}

				}
				else if(up_face == 1)
				{

					// 把底部靠左的link pop掉
					for(i=0;i<m_links.size();i++)
					{
						if( ((m_links[i].m_n[0]-nbase) == left_node && (m_links[i].m_n[1]-nbase) == info[last_face[fi]].Nrcut) || 
							((m_links[i].m_n[1]-nbase) == left_node && (m_links[i].m_n[0]-nbase) == info[last_face[fi]].Nrcut))
						{
							btSwap(m_links[i], m_links[m_links.size()-1]);  // 單純做template型態的swap, 將m_link第i個陣列的東西, 跟最後一個陣列的東西交換
							m_links.pop_back(); // pop_back() 就是將m_size--, 並做destructor
							--i; 
						}	
					}

				}

				if(crossface == 1)
				{
					// 把底部link pop掉
					for(i=0;i<m_links.size();i++)
					{
						if( ((m_links[i].m_n[0]-nbase) == left_node && (m_links[i].m_n[1]-nbase) == right_node) || 
							((m_links[i].m_n[1]-nbase) == left_node && (m_links[i].m_n[0]-nbase) == right_node) )
						{
							btSwap(m_links[i], m_links[m_links.size()-1]);  // 單純做template型態的swap, 將m_link第i個陣列的東西, 跟最後一個陣列的東西交換
							m_links.pop_back(); // pop_back() 就是將m_size--, 並做destructor
							--i; 
						}

					}

					// 重新連接link
					appendLink(right_node,edge_rcut); // 1 12
					appendLink(edge_lcut,left_node); // 13 2

					// 把底下link pop掉

					for(i=0;i<m_links.size();i++)
					{
						if(poplr == 1)
						{
							if( ((m_links[i].m_n[0]-nbase) == info[last_face[fi]].Nimpact && (m_links[i].m_n[1]-nbase) == info[last_face[fi]].Nrcut) || 
								((m_links[i].m_n[1]-nbase) == info[last_face[fi]].Nimpact && (m_links[i].m_n[0]-nbase) == info[last_face[fi]].Nrcut) )
							{
								btSwap(m_links[i], m_links[m_links.size()-1]);  // 單純做template型態的swap, 將m_link第i個陣列的東西, 跟最後一個陣列的東西交換
								m_links.pop_back(); // pop_back() 就是將m_size--, 並做destructor
								--i; 
							}

						}
						else if(poplr == 0)
						{
							if( ((m_links[i].m_n[0]-nbase) == info[last_face[fi]].Nimpact && (m_links[i].m_n[1]-nbase) == info[last_face[fi]].Nlcut) || 
								((m_links[i].m_n[1]-nbase) == info[last_face[fi]].Nimpact && (m_links[i].m_n[0]-nbase) == info[last_face[fi]].Nlcut) )
							{
								btSwap(m_links[i], m_links[m_links.size()-1]);  // 單純做template型態的swap, 將m_link第i個陣列的東西, 跟最後一個陣列的東西交換
								m_links.pop_back(); // pop_back() 就是將m_size--, 並做destructor
								--i; 
							}

						}
					}


				}
				if(cutedge == 1) // 剪斷底邊
				{
					// 把底部link pop掉

					if(crossface == 1)
					{
						if(poplr == 0)
						{
							for(i=0;i<m_links.size();i++)
							{
								if( ((m_links[i].m_n[0]-nbase) == up_node && (m_links[i].m_n[1]-nbase) == right_node) || 
									((m_links[i].m_n[1]-nbase) == up_node && (m_links[i].m_n[0]-nbase) == right_node) )
								{
									btSwap(m_links[i], m_links[m_links.size()-1]);  // 單純做template型態的swap, 將m_link第i個陣列的東西, 跟最後一個陣列的東西交換
									m_links.pop_back(); // pop_back() 就是將m_size--, 並做destructor
									--i; 
								}
							}
							// 重新連接link
							appendLink(right_node,rcut_node); 
							appendLink(up_node,lcut_node); 
						}
						else if(poplr == 1)
						{
							for(i=0;i<m_links.size();i++)
							{
								if( ((m_links[i].m_n[0]-nbase) == up_node && (m_links[i].m_n[1]-nbase) == left_node) || 
									((m_links[i].m_n[1]-nbase) == up_node && (m_links[i].m_n[0]-nbase) == left_node) )
								{
									btSwap(m_links[i], m_links[m_links.size()-1]);  // 單純做template型態的swap, 將m_link第i個陣列的東西, 跟最後一個陣列的東西交換
									m_links.pop_back(); // pop_back() 就是將m_size--, 並做destructor
									--i; 
								}
							}
							// 重新連接link
							appendLink(left_node,lcut_node); 
							appendLink(up_node,rcut_node); 

						}
					}
					if(up_face == 0)  // 往左上face
					{
						for(i=0;i<m_links.size();i++)
						{
							if( ((m_links[i].m_n[0]-nbase) == up_node && (m_links[i].m_n[1]-nbase) == left_node) || 
								((m_links[i].m_n[1]-nbase) == up_node && (m_links[i].m_n[0]-nbase) == left_node) )
							{
								btSwap(m_links[i], m_links[m_links.size()-1]);  // 單純做template型態的swap, 將m_link第i個陣列的東西, 跟最後一個陣列的東西交換
								m_links.pop_back(); // pop_back() 就是將m_size--, 並做destructor
								--i; 
							}

						}

						// 重新連接link
						appendLink(left_node,lcut_node); 
						appendLink(up_node,impact_node); 
					}
					else if(up_face == 1)  // 往右上face
					{
						for(i=0;i<m_links.size();i++)
						{
							if( ((m_links[i].m_n[0]-nbase) == up_node && (m_links[i].m_n[1]-nbase) == right_node) || 
								((m_links[i].m_n[1]-nbase) == up_node && (m_links[i].m_n[0]-nbase) == right_node) )
							{
								btSwap(m_links[i], m_links[m_links.size()-1]);  // 單純做template型態的swap, 將m_link第i個陣列的東西, 跟最後一個陣列的東西交換
								m_links.pop_back(); // pop_back() 就是將m_size--, 並做destructor
								--i; 
							}

						}

						// 重新連接link
						appendLink(right_node,rcut_node); 
						appendLink(up_node,impact_node);

					}
				}


				// 記錄當前impact,rcut,lcut編號
				//if(f == 1) // 只有點擊在一塊face上 (非對摺情況的話), 第一刀完要做記錄
				//{

				info[facenum[fi]].Nimpact = impact_node;
				info[facenum[fi]].Nrcut = rcut_node;
				info[facenum[fi]].Nlcut = lcut_node;
				info[facenum[fi]].Nright = right_node;
				info[facenum[fi]].Nleft = left_node;
				info[facenum[fi]].Nup = up_node;

				last_face[fi] = facenum[fi];
				//}

			}
		}


	}
	if(newnodes>0)
	{
		strcat( outputInfo, "You cut the paper!\n" );
	}
	for(i=0;i<m_faces.size();++i)
	{
		fprintf(frr,"%d face : %d %d %d\n",i,m_faces[i].m_n[0]-nbase,m_faces[i].m_n[1]-nbase,m_faces[i].m_n[2]-nbase);
	}

	for(i=0;i<m_links.size();i++)
	{
		fprintf(frr,"%d link : %d %d\n",i,m_links[i].m_n[0]-nbase,m_links[i].m_n[1]-nbase);
	}
	fclose(frr);
	m_bUpdateRtCst=true;

	// 將node的速度歸零
	for(i=0;i<m_nodes.size();i++)
		m_nodes[i].m_v = V;

}


//
void			btSoftBody::refine(ImplicitFn* ifn,btScalar accurary,bool cut,btVector3& impact)
{

	const Node*			nbase = &m_nodes[0];
	int					ncount = m_nodes.size();
	btSymMatrix<int>	edges(ncount,-2);
	int					newnodes=0;
	int i,j,k,ni;

	cur_step=4;
	
	for(i=0;i<9;i++)
		if(m_nodes[i].m_im==0)
			m_nodes[i].m_im=oldmass;


	if(state==1)  //下往上摺
	{
		m_nodes[0].m_im=0;m_nodes[2].m_im=0;
		m_nodes[6].m_im=0;m_nodes[8].m_im=0;
		m_nodes[1].m_im=0;m_nodes[7].m_im=0;
	}
	else if(state==2)   //左下右上摺
	{
		m_nodes[0].m_im=0;m_nodes[2].m_im=0;
		m_nodes[6].m_im=0;m_nodes[8].m_im=0;
	}
	else if(state==3)  //右往左摺
	{
		m_nodes[3].m_im=0;m_nodes[5].m_im=0;
		m_nodes[1].m_im=0;m_nodes[7].m_im=0;
	}
	else if(state==4)  //右下往左上摺
	{
		m_nodes[0].m_im=0;m_nodes[2].m_im=0;
		m_nodes[6].m_im=0;m_nodes[8].m_im=0;
	}
	else if(state==5)  //左往右摺
	{
		m_nodes[3].m_im=0;m_nodes[5].m_im=0;
		m_nodes[1].m_im=0;m_nodes[7].m_im=0;
	}
	else if(state==6)   //右上左下摺
	{
		m_nodes[0].m_im=0;m_nodes[2].m_im=0;
		m_nodes[6].m_im=0;m_nodes[8].m_im=0;
	}
	else if(state==7)  //上往下摺
	{
		m_nodes[3].m_im=0;m_nodes[5].m_im=0;
		m_nodes[1].m_im=0;m_nodes[7].m_im=0;
	}
	else if(state==8)  //左上往右下摺
	{
		m_nodes[0].m_im=0;m_nodes[2].m_im=0;
		m_nodes[6].m_im=0;m_nodes[8].m_im=0;
	}
	else 
	{
		m_nodes[6].m_im=0;
		m_nodes[7].m_im=0;
		m_nodes[8].m_im=0;
	}

	/* Filter out		*/ 
	for(i=0;i<m_links.size();++i)
	{
		Link&	l=m_links[i];

		if(l.m_bbending)
		{
			if(!SameSign(ifn->Eval(l.m_n[0]->m_x),ifn->Eval(l.m_n[1]->m_x)))
			{
				btSwap(m_links[i],m_links[m_links.size()-1]);
				m_links.pop_back();--i;
			}

		}	
	}
	/* Fill edges		*/ 
	for(i=0;i<m_links.size();++i)
	{
		Link&	l=m_links[i];
		edges(int(l.m_n[0]-nbase),int(l.m_n[1]-nbase))=-1;
	}
	for(i=0;i<m_faces.size();++i)
	{	
		Face&	f=m_faces[i];
		edges(int(f.m_n[0]-nbase),int(f.m_n[1]-nbase))=-1;
		edges(int(f.m_n[1]-nbase),int(f.m_n[2]-nbase))=-1;
		edges(int(f.m_n[2]-nbase),int(f.m_n[0]-nbase))=-1;
	}
	/* Intersect		*/ 
	for(i=0;i<ncount;++i)
	{
		for(j=i+1;j<ncount;++j)
		{
			if(edges(i,j)==-1)
			{

				Node&			a=m_nodes[i];
				Node&			b=m_nodes[j];
				const btScalar	t=ImplicitSolve(ifn,a.m_x,b.m_x,accurary);

				if(t>0)
				{
					const btVector3	x=Lerp(a.m_x,b.m_x,t);  //m_x :position
					const btVector3	v=Lerp(a.m_v,b.m_v,t);	//m_v :velocity
					btScalar		m=0;
					if(a.m_im>0)
					{
						if(b.m_im>0)
						{
							const btScalar	ma=1/a.m_im;
							const btScalar	mb=1/b.m_im;
							const btScalar	mc=Lerp(ma,mb,t);
							const btScalar	f=(ma+mb)/(ma+mb+mc);
							a.m_im=1/(ma*f);
							b.m_im=1/(mb*f);
							m=mc*f;
						}
						else
						{ a.m_im/=0.5f;m=1/a.m_im;
						}
					}
					else
					{
						if(b.m_im>0)
						{ b.m_im/=0.5f;m=1/b.m_im; }
						else
							m=oldmass;
					}

					appendNode(x,m);
					edges(i,j)=m_nodes.size()-1;
					m_nodes[edges(i,j)].m_v=v;//v
					++newnodes;

				}
			}
		}
	}
	nbase=&m_nodes[0];

	/* Refine links		*/ 
	//把link剪斷的部分接回來

	for(i=0,ni=m_links.size();i<ni;++i)
	{
		Link&		feat=m_links[i];
		const int	idx[]={	int(feat.m_n[0]-nbase),
			int(feat.m_n[1]-nbase)};

		if((idx[0]<ncount)&&(idx[1]<ncount))
		{
			const int ni=edges(idx[0],idx[1]); //ni為有剪到的邊 若有則不為-1 
			if(ni>0)
			{
				appendLink(i);
				Link*		pft[]={	&m_links[i],
					&m_links[m_links.size()-1]};

				pft[0]->m_n[0]=&m_nodes[idx[0]];
				pft[0]->m_n[1]=&m_nodes[ni];
				pft[1]->m_n[0]=&m_nodes[ni];
				pft[1]->m_n[1]=&m_nodes[idx[1]];
				//pft[0]為還留在布上的link ，pft[1]為掉下布的link。
			}
		}
	}

	/* Refine faces		*/ 
	for(i=0;i<m_faces.size();++i)
	{
		const Face&	feat=m_faces[i];
		const int	idx[]={	int(feat.m_n[0]-nbase),
			int(feat.m_n[1]-nbase),
			int(feat.m_n[2]-nbase)};
		for(j=2,k=0;k<3;j=k++)
		{
			if((idx[j]<ncount)&&(idx[k]<ncount))
			{
				const int ni=edges(idx[j],idx[k]);
				if(ni>0)
				{
					appendFace(i);
					const int	l=(k+1)%3;
					Face*		pft[]={	&m_faces[i],
						&m_faces[m_faces.size()-1]};

					pft[0]->m_n[0]=&m_nodes[idx[l]];
					pft[0]->m_n[1]=&m_nodes[idx[j]];
					pft[0]->m_n[2]=&m_nodes[ni];

					pft[1]->m_n[0]=&m_nodes[ni];
					pft[1]->m_n[1]=&m_nodes[idx[k]];
					pft[1]->m_n[2]=&m_nodes[idx[l]];
					appendLink(ni,idx[l],pft[0]->m_material);
					--i;break;
				}
			}
		}
	}
	/* Cut				*/ 

	if(cut)
	{	
		btAlignedObjectArray<int>	cnodes;  //掉下布的陣列
		const int					pcount=ncount;  //pcount=剪之前的node數量
		int							i;
		ncount=m_nodes.size();
		cnodes.resize(ncount,0);    //重新分配記憶體
		//printf("size:%d\n",cnodes.size());
		/* Nodes		*/ 
		for(i=0;i<ncount;++i)   //ncount=剪完還沒落下時node數量
		{
			const btVector3	x=m_nodes[i].m_x;
			if((i>=pcount)||(btFabs(ifn->Eval(x))<accurary)) //btFabs將浮點數取絕對值 accuracy這邊是抓圓上附近的點而非圓內。
			{
				const btVector3	v=m_nodes[i].m_v;
				btScalar		m=getMass(i);
				if(m>0) { m*=0.5f;
				m_nodes[i].m_im/=0.5f;
				}
				appendNode(x,m);   //將落下的布新增node nodesize會++
				cnodes[i]=m_nodes.size()-1;
				m_nodes[cnodes[i]].m_v=v;
			}
		}
		nbase=&m_nodes[0];
		/* Links		*/ 
		for(i=0,ni=m_links.size();i<ni;++i)   //在切割face時已做完 未處理落下的部分
		{
			const int		id[]={	int(m_links[i].m_n[0]-nbase),
				int(m_links[i].m_n[1]-nbase)};
			int				todetach=0;
			if(cnodes[id[0]]&&cnodes[id[1]])
			{
				// 被剪到的node 彼此之間再新增一條link 
				appendLink(i);
				todetach=m_links.size()-1;

			}
			else
			{
				if((	(ifn->Eval(m_nodes[id[0]].m_x)<accurary)&&
					(ifn->Eval(m_nodes[id[1]].m_x)<accurary)))
					todetach=i;
			}

			if(todetach)    //改動原link的node
			{
				Link&	l=m_links[todetach];

				for(int j=0;j<2;++j)
				{
					int cn=cnodes[int(l.m_n[j]-nbase)];
					if(cn) l.m_n[j]=&m_nodes[cn];
				}
			}
		}
		/* Faces		*/ 
		for(i=0,ni=m_faces.size();i<ni;++i)
		{
			Node**			n=	m_faces[i].m_n;  //剪掉的那塊face
			if(	(ifn->Eval(n[0]->m_x)<accurary)&&  //都在剪取範圍內的三角形
				(ifn->Eval(n[1]->m_x)<accurary)&&
				(ifn->Eval(n[2]->m_x)<accurary))
			{
				for(int j=0;j<3;++j)
				{
					int cn=cnodes[int(n[j]-nbase)];
					if(cn) n[j]=&m_nodes[cn];
				}
			}
		}
		/* Clean orphans	*/ 


		int							nnodes=m_nodes.size();  //分離後的node總數
		btAlignedObjectArray<int>	ranks; //紀錄每個node出現次數 越高代表越多連越多link和face
		btAlignedObjectArray<int>	todelete;
		ranks.resize(nnodes,0);
		for(i=0,ni=m_links.size();i<ni;++i)
		{
			for(int j=0;j<2;++j) 
			{
				ranks[int(m_links[i].m_n[j]-nbase)]++;
			}

		}
		for(i=0,ni=m_faces.size();i<ni;++i)
		{
			for(int j=0;j<3;++j) ranks[int(m_faces[i].m_n[j]-nbase)]++;
		}
		for(i=0;i<m_links.size();++i)
		{
			const int	id[]={	int(m_links[i].m_n[0]-nbase),
				int(m_links[i].m_n[1]-nbase)};
			const bool	sg[]={	ranks[id[0]]==1,
				ranks[id[1]]==1};

			if(sg[0]||sg[1])
			{
				--ranks[id[0]];
				--ranks[id[1]];
				btSwap(m_links[i],m_links[m_links.size()-1]);
				m_links.pop_back();--i;
			}
		}

#if 0	
		for(i=nnodes-1;i>=0;--i)
		{
			if(!ranks[i]) todelete.push_back(i);
		}	
		if(todelete.size())
		{		
			btAlignedObjectArray<int>&	map=ranks;
			for(int i=0;i<nnodes;++i) map[i]=i;
			pointersToIndices();
			for(int i=0,ni=todelete.size();i<ni;++i)
			{
				int		j=todelete[i];
				int&	a=map[j];
				int&	b=map[--nnodes];
				m_ndbvt.remove(m_nodes[a].m_leaf);m_nodes[a].m_leaf=0;
				btSwap(m_nodes[a],m_nodes[b]);
				j=a;a=b;b=j;			
			}
			indicesToPointers(&map[0]);
			//indicesToPointers(this,&map[0]);
			m_nodes.resize(nnodes);
		}
#endif

	}
	if(newnodes>0)
	{
		strcat( outputInfo, "You cut the paper!\n" );
	}
	
	m_bUpdateRtCst=true;
}
//
bool            btSoftBody::cutLink(const Node *node0, const Node *node1, btScalar position)
{
	return (cutLink(int(node0 - &m_nodes[0]), int(node1 - &m_nodes[0]), position));
}

//
bool            btSoftBody::cutLink(int node0, int node1, btScalar position)
{
	bool            done = false;
	int i, ni;
	//  const btVector3 d=m_nodes[node0].m_x-m_nodes[node1].m_x;
	const btVector3 x = Lerp(m_nodes[node0].m_x, m_nodes[node1].m_x, position);
	const btVector3 v = Lerp(m_nodes[node0].m_v, m_nodes[node1].m_v, position);
	const btScalar  m = 1;
	appendNode(x, m);
	appendNode(x, m);
	Node           *pa = &m_nodes[node0];
	Node           *pb = &m_nodes[node1];
	Node           *pn[2] = { &m_nodes[m_nodes.size() - 2],
		&m_nodes[m_nodes.size() - 1]
	};
	pn[0]->m_v = v;
	pn[1]->m_v = v;
	for (i = 0, ni = m_links.size(); i < ni; ++i)
	{
		const int mtch = MatchEdge(m_links[i].m_n[0], m_links[i].m_n[1], pa, pb);
		if (mtch != -1)
		{
			appendLink(i);
			Link   *pft[] = {&m_links[i], &m_links[m_links.size() - 1]};
			pft[0]->m_n[1] = pn[mtch];
			pft[1]->m_n[0] = pn[1 - mtch];
			done = true;
		}
	}
	for (i = 0, ni = m_faces.size(); i < ni; ++i)
	{
		for (int k = 2, l = 0; l < 3; k = l++)
		{
			const int mtch = MatchEdge(m_faces[i].m_n[k], m_faces[i].m_n[l], pa, pb);
			if (mtch != -1)
			{
				appendFace(i);
				Face   *pft[] = {&m_faces[i], &m_faces[m_faces.size() - 1]};
				pft[0]->m_n[l] = pn[mtch];
				pft[1]->m_n[k] = pn[1 - mtch];
				appendLink(pn[0], pft[0]->m_n[(l + 1) % 3], pft[0]->m_material, true);
				appendLink(pn[1], pft[0]->m_n[(l + 1) % 3], pft[0]->m_material, true);
			}
		}
	}
	if (!done)
	{
		m_ndbvt.remove(pn[0]->m_leaf);
		m_ndbvt.remove(pn[1]->m_leaf);
		m_nodes.pop_back();
		m_nodes.pop_back();
	}
	return (done);
}

//
bool            btSoftBody::rayTest(const btVector3 &rayFrom,
									const btVector3 &rayTo,
									sRayCast &results)
{
	//printf("m_facesize:%d\n",m_faces.size());
	if (m_faces.size() && m_fdbvt.empty())
		initializeFaceTree();

	results.body    =   this;
	results.fraction = 1.f;
	results.feature =   eFeature::None;
	results.index   =   -1;

	return (rayTest(rayFrom, rayTo, results.fraction, results.feature, results.index, false) != 0);
}

//
void            btSoftBody::setSolver(eSolverPresets::_ preset)
{
	m_cfg.m_vsequence.clear();
	m_cfg.m_psequence.clear();
	m_cfg.m_dsequence.clear();
	switch (preset)
	{
	case    eSolverPresets::Positions:
		m_cfg.m_psequence.push_back(ePSolver::Anchors);
		m_cfg.m_psequence.push_back(ePSolver::RContacts);
		m_cfg.m_psequence.push_back(ePSolver::SContacts);
		m_cfg.m_psequence.push_back(ePSolver::Linear);
		break;
	case    eSolverPresets::Velocities:
		m_cfg.m_vsequence.push_back(eVSolver::Linear);

		m_cfg.m_psequence.push_back(ePSolver::Anchors);
		m_cfg.m_psequence.push_back(ePSolver::RContacts);
		m_cfg.m_psequence.push_back(ePSolver::SContacts);

		m_cfg.m_dsequence.push_back(ePSolver::Linear);
		break;
	}
}

//
void            btSoftBody::predictMotion(btScalar dt)
{

	int i, ni;

	/* Update               */
	if (m_bUpdateRtCst)
	{
		m_bUpdateRtCst = false;
		updateConstants();
		m_fdbvt.clear();
		if (m_cfg.collisions & fCollision::VF_SS)
		{
			initializeFaceTree();
		}
	}

	/* Prepare              */
	m_sst.sdt       =   dt * m_cfg.timescale;
	m_sst.isdt      =   1 / m_sst.sdt;
	m_sst.velmrg    =   m_sst.sdt * 3;
	m_sst.radmrg    =   getCollisionShape()->getMargin();
	m_sst.updmrg    =   m_sst.radmrg * (btScalar)0.25;
	/* Forces               */
	addVelocity(m_worldInfo->m_gravity * m_sst.sdt);
	applyForces();
	/* Integrate            */
	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		Node   &n = m_nodes[i];
		n.m_q   =   n.m_x;
		btVector3 deltaV = n.m_f * n.m_im * m_sst.sdt;
		{
			btScalar maxDisplacement = m_worldInfo->m_maxDisplacement;
			btScalar clampDeltaV = maxDisplacement / m_sst.sdt;
			for (int c = 0; c < 3; c++)
			{
				if (deltaV[c] > clampDeltaV)
				{
					deltaV[c] = clampDeltaV;
				}
				if (deltaV[c] < -clampDeltaV)
				{
					deltaV[c] = -clampDeltaV;
				}
			}
		}
		n.m_v   +=  deltaV;
		n.m_x   +=  n.m_v * m_sst.sdt;
		n.m_f   =   btVector3(0, 0, 0);
	}
	/* Clusters             */
	updateClusters();
	/* Bounds               */
	updateBounds();
	/* Nodes                */
	ATTRIBUTE_ALIGNED16(btDbvtVolume)   vol;
	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		Node   &n = m_nodes[i];
		vol = btDbvtVolume::FromCR(n.m_x, m_sst.radmrg);
		m_ndbvt.update( n.m_leaf,
			vol,
			n.m_v * m_sst.velmrg,
			m_sst.updmrg);
	}
	/* Faces                */
	if (!m_fdbvt.empty())
	{
		for (int i = 0; i < m_faces.size(); ++i)
		{
			Face           &f = m_faces[i];
			const btVector3 v = ( f.m_n[0]->m_v +
				f.m_n[1]->m_v +
				f.m_n[2]->m_v) / 3;
			vol = VolumeOf(f, m_sst.radmrg);
			m_fdbvt.update( f.m_leaf,
				vol,
				v * m_sst.velmrg,
				m_sst.updmrg);
		}
	}
	/* Pose                 */
	updatePose();
	/* Match                */
	if (m_pose.m_bframe && (m_cfg.kMT > 0))
	{
		const btMatrix3x3   posetrs = m_pose.m_rot;
		for (int i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			Node   &n = m_nodes[i];
			if (n.m_im > 0)
			{
				const btVector3 x = posetrs * m_pose.m_pos[i] + m_pose.m_com;
				n.m_x = Lerp(n.m_x, x, m_cfg.kMT);
			}
		}
	}
	/* Clear contacts       */
	m_rcontacts.resize(0);
	m_scontacts.resize(0);
	/* Optimize dbvt's      */
	m_ndbvt.optimizeIncremental(1);
	m_fdbvt.optimizeIncremental(1);
	m_cdbvt.optimizeIncremental(1);
}

//
void            btSoftBody::solveConstraints()
{

	/* Apply clusters       */
	applyClusters(false);
	/* Prepare links        */

	int i, ni;

	for (i = 0, ni = m_links.size(); i < ni; ++i)
	{
		Link   &l = m_links[i];
		l.m_c3      =   l.m_n[1]->m_q - l.m_n[0]->m_q;
		l.m_c2      =   1 / (l.m_c3.length2() * l.m_c0);
	}
	/* Prepare anchors      */
	for (i = 0, ni = m_anchors.size(); i < ni; ++i)
	{
		Anchor         &a = m_anchors[i];
		const btVector3 ra = a.m_body->getWorldTransform().getBasis() * a.m_local;
		a.m_c0  =   ImpulseMatrix(  m_sst.sdt,
			a.m_node->m_im,
			a.m_body->getInvMass(),
			a.m_body->getInvInertiaTensorWorld(),
			ra);
		a.m_c1  =   ra;
		a.m_c2  =   m_sst.sdt * a.m_node->m_im;
		a.m_body->activate();
	}
	/* Solve velocities     */
	if (m_cfg.viterations > 0)
	{
		/* Solve            */
		for (int isolve = 0; isolve < m_cfg.viterations; ++isolve)
		{
			for (int iseq = 0; iseq < m_cfg.m_vsequence.size(); ++iseq)
			{
				getSolver(m_cfg.m_vsequence[iseq])(this, 1);
			}
		}
		/* Update           */
		for (i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			Node   &n = m_nodes[i];
			n.m_x   =   n.m_q + n.m_v * m_sst.sdt;
		}
	}
	/* Solve positions      */
	if (m_cfg.piterations > 0)
	{
		for (int isolve = 0; isolve < m_cfg.piterations; ++isolve)
		{
			const btScalar ti = isolve / (btScalar)m_cfg.piterations;
			for (int iseq = 0; iseq < m_cfg.m_psequence.size(); ++iseq)
			{
				getSolver(m_cfg.m_psequence[iseq])(this, 1, ti);
			}
		}
		const btScalar  vc = m_sst.isdt * (1 - m_cfg.kDP);
		for (i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			Node   &n = m_nodes[i];
			n.m_v   =   (n.m_x - n.m_q) * vc;
			n.m_f   =   btVector3(0, 0, 0);
		}
	}
	/* Solve drift          */
	if (m_cfg.diterations > 0)
	{
		const btScalar  vcf = m_cfg.kVCF * m_sst.isdt;
		for (i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			Node   &n = m_nodes[i];
			n.m_q   =   n.m_x;
		}
		for (int idrift = 0; idrift < m_cfg.diterations; ++idrift)
		{
			for (int iseq = 0; iseq < m_cfg.m_dsequence.size(); ++iseq)
			{
				getSolver(m_cfg.m_dsequence[iseq])(this, 1, 0);
			}
		}
		for (int i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			Node   &n = m_nodes[i];
			n.m_v   +=  (n.m_x - n.m_q) * vcf;
		}
	}
	/* Apply clusters       */
	dampClusters();
	applyClusters(true);
}

//
void            btSoftBody::staticSolve(int iterations)
{
	for (int isolve = 0; isolve < iterations; ++isolve)
	{
		for (int iseq = 0; iseq < m_cfg.m_psequence.size(); ++iseq)
		{
			getSolver(m_cfg.m_psequence[iseq])(this, 1, 0);
		}
	}
}

//
void            btSoftBody::solveCommonConstraints(btSoftBody ** /*bodies*/, int /*count*/, int /*iterations*/)
{
	/// placeholder
}

//
void            btSoftBody::solveClusters(const btAlignedObjectArray<btSoftBody *> &bodies)
{
	const int   nb = bodies.size();
	int         iterations = 0;
	int i;

	for (i = 0; i < nb; ++i)
	{
		iterations = btMax(iterations, bodies[i]->m_cfg.citerations);
	}
	for (i = 0; i < nb; ++i)
	{
		bodies[i]->prepareClusters(iterations);
	}
	for (i = 0; i < iterations; ++i)
	{
		const btScalar sor = 1;
		for (int j = 0; j < nb; ++j)
		{
			bodies[j]->solveClusters(sor);
		}
	}
	for (i = 0; i < nb; ++i)
	{
		bodies[i]->cleanupClusters();
	}
}

//
void            btSoftBody::integrateMotion()
{
	/* Update           */
	updateNormals();
}

//
btSoftBody::RayFromToCaster::RayFromToCaster(const btVector3 &rayFrom, const btVector3 &rayTo, btScalar mxt)
{
	m_rayFrom = rayFrom;
	m_rayNormalizedDirection = (rayTo - rayFrom);
	m_rayTo = rayTo;
	m_mint  =   mxt;
	m_face  =   0;
	m_tests =   0;
}

//
void                btSoftBody::RayFromToCaster::Process(const btDbvtNode *leaf)
{
	btSoftBody::Face   &f = *(btSoftBody::Face *)leaf->data;
	const btScalar      t = rayFromToTriangle(    m_rayFrom, m_rayTo, m_rayNormalizedDirection,
		f.m_n[0]->m_x,
		f.m_n[1]->m_x,
		f.m_n[2]->m_x,
		m_mint);
	if ((t > 0) && (t < m_mint))
	{
		m_mint = t; m_face = &f;
	}
	++m_tests;
}

//
btScalar            btSoftBody::RayFromToCaster::rayFromToTriangle( const btVector3 &rayFrom,
																   const btVector3 &rayTo,
																   const btVector3 &rayNormalizedDirection,
																   const btVector3 &a,
																   const btVector3 &b,
																   const btVector3 &c,
																   btScalar maxt)
{
	static const btScalar   ceps = -SIMD_EPSILON * 10;
	static const btScalar   teps = SIMD_EPSILON * 10;

	const btVector3         n = btCross(b - a, c - a);
	const btScalar          d = btDot(a, n);
	const btScalar          den = btDot(rayNormalizedDirection, n);
	if (!btFuzzyZero(den))
	{
		const btScalar      num = btDot(rayFrom, n) - d;
		const btScalar      t = -num / den;
		if ((t > teps) && (t < maxt))
		{
			const btVector3 hit = rayFrom + rayNormalizedDirection * t;
			if ( (btDot(n, btCross(a - hit, b - hit)) > ceps)    &&
				(btDot(n, btCross(b - hit, c - hit)) > ceps)    &&
				(btDot(n, btCross(c - hit, a - hit)) > ceps))
			{
				return (t);
			}
		}
	}
	return (-1);
}

//
void                btSoftBody::pointersToIndices()  //將node pointer儲存到face和link
{
#define PTR2IDX(_p_,_b_)    reinterpret_cast<btSoftBody::Node*>((_p_)-(_b_)) //將object轉成node型態
	btSoftBody::Node   *base = m_nodes.size() ? &m_nodes[0] : 0;
	int i, ni;

	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		if (m_nodes[i].m_leaf)
		{
			m_nodes[i].m_leaf->data = *(void **)&i;
		}
	}
	for (i = 0, ni = m_links.size(); i < ni; ++i)
	{
		m_links[i].m_n[0] = PTR2IDX(m_links[i].m_n[0], base);
		m_links[i].m_n[1] = PTR2IDX(m_links[i].m_n[1], base);
	}
	for (i = 0, ni = m_faces.size(); i < ni; ++i)
	{
		m_faces[i].m_n[0] = PTR2IDX(m_faces[i].m_n[0], base);
		m_faces[i].m_n[1] = PTR2IDX(m_faces[i].m_n[1], base);
		m_faces[i].m_n[2] = PTR2IDX(m_faces[i].m_n[2], base);
		if (m_faces[i].m_leaf)
		{
			m_faces[i].m_leaf->data = *(void **)&i;
		}
	}
	for (i = 0, ni = m_anchors.size(); i < ni; ++i)
	{
		m_anchors[i].m_node = PTR2IDX(m_anchors[i].m_node, base);
	}
	for (i = 0, ni = m_notes.size(); i < ni; ++i)
	{
		for (int j = 0; j < m_notes[i].m_rank; ++j)
		{
			m_notes[i].m_nodes[j] = PTR2IDX(m_notes[i].m_nodes[j], base);
		}
	}
#undef  PTR2IDX
}

//
void                btSoftBody::indicesToPointers(const int *map)
{
#define IDX2PTR(_p_,_b_)    map?(&(_b_)[map[(((char*)_p_)-(char*)0)]]): \
	(&(_b_)[(((char*)_p_)-(char*)0)])
	btSoftBody::Node   *base = m_nodes.size() ? &m_nodes[0] : 0;
	int i, ni;

	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		if (m_nodes[i].m_leaf)
		{
			m_nodes[i].m_leaf->data = &m_nodes[i];
		}
	}
	for (i = 0, ni = m_links.size(); i < ni; ++i)
	{
		m_links[i].m_n[0] = IDX2PTR(m_links[i].m_n[0], base);
		m_links[i].m_n[1] = IDX2PTR(m_links[i].m_n[1], base);
	}
	for (i = 0, ni = m_faces.size(); i < ni; ++i)
	{
		m_faces[i].m_n[0] = IDX2PTR(m_faces[i].m_n[0], base);
		m_faces[i].m_n[1] = IDX2PTR(m_faces[i].m_n[1], base);
		m_faces[i].m_n[2] = IDX2PTR(m_faces[i].m_n[2], base);
		if (m_faces[i].m_leaf)
		{
			m_faces[i].m_leaf->data = &m_faces[i];
		}
	}
	for (i = 0, ni = m_anchors.size(); i < ni; ++i)
	{
		m_anchors[i].m_node = IDX2PTR(m_anchors[i].m_node, base);
	}
	for (i = 0, ni = m_notes.size(); i < ni; ++i)
	{
		for (int j = 0; j < m_notes[i].m_rank; ++j)
		{
			m_notes[i].m_nodes[j] = IDX2PTR(m_notes[i].m_nodes[j], base);
		}
	}
#undef  IDX2PTR
}

//
int                 btSoftBody::rayTest(const btVector3 &rayFrom, const btVector3 &rayTo,
										btScalar &mint, eFeature::_& feature, int &index, bool bcountonly) const
{
	int cnt = 0;
	btVector3 dir = rayTo - rayFrom;


	if (bcountonly || m_fdbvt.empty())
	{
		/* Full search */
		//printf("m_face:%d\n",m_faces.size());
		for (int i = 0, ni = m_faces.size(); i < ni; ++i)
		{
			const btSoftBody::Face &f = m_faces[i];

			const btScalar          t = RayFromToCaster::rayFromToTriangle(   rayFrom, rayTo, dir,
				f.m_n[0]->m_x,
				f.m_n[1]->m_x,
				f.m_n[2]->m_x,
				mint);
			if (t > 0)
			{
				++cnt;
				if (!bcountonly)
				{
					feature = btSoftBody::eFeature::Face;
					index = i;
					mint = t;
				}
			}
		}
	}
	else
	{
		/* Use dbvt    */
		RayFromToCaster collider(rayFrom, rayTo, mint);

		btDbvt::rayTest(m_fdbvt.m_root, rayFrom, rayTo, collider);
		//printf("collider.face:%d\n",collider.m_face);
		if (collider.m_face)
		{
			mint = collider.m_mint;
			//printf("mint:%f\n",collider.m_mint);
			feature = btSoftBody::eFeature::Face;
			index = (int)(collider.m_face - &m_faces[0]);
			cnt = 1;

		}

	}
	//m_tetras.size=0;
	for (int i = 0; i < m_tetras.size(); i++)
	{
		const btSoftBody::Tetra &tet = m_tetras[i];
		int tetfaces[4][3] = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};
		for (int f = 0; f < 4; f++)
		{

			int index0 = tetfaces[f][0];
			int index1 = tetfaces[f][1];
			int index2 = tetfaces[f][2];
			btVector3 v0 = tet.m_n[index0]->m_x;
			btVector3 v1 = tet.m_n[index1]->m_x;
			btVector3 v2 = tet.m_n[index2]->m_x;


			const btScalar          t = RayFromToCaster::rayFromToTriangle(   rayFrom, rayTo, dir,
				v0, v1, v2,
				mint);
			if (t > 0)
			{
				++cnt;
				if (!bcountonly)
				{
					feature = btSoftBody::eFeature::Tetra;
					index = i;
					mint = t;
				}
			}
		}
	}

	return (cnt);
}

//
void            btSoftBody::initializeFaceTree()
{
	m_fdbvt.clear();
	for (int i = 0; i < m_faces.size(); ++i)
	{
		Face   &f = m_faces[i];
		f.m_leaf = m_fdbvt.insert(VolumeOf(f, 0), &f);
	}
}

//
btVector3       btSoftBody::evaluateCom() const
{
	btVector3   com(0, 0, 0);
	if (m_pose.m_bframe)
	{
		for (int i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			com += m_nodes[i].m_x * m_pose.m_wgh[i];
		}
	}
	return (com);
}

//
bool                btSoftBody::checkContact(   const btCollisionObjectWrapper *colObjWrap,
											 const btVector3 &x,
											 btScalar margin,
											 btSoftBody::sCti &cti) const
{
	btVector3 nrm;
	const btCollisionShape *shp = colObjWrap->getCollisionShape();
	//  const btRigidBody *tmpRigid = btRigidBody::upcast(colObjWrap->getCollisionObject());
	//const btTransform &wtr = tmpRigid ? tmpRigid->getWorldTransform() : colObjWrap->getWorldTransform();
	const btTransform &wtr = colObjWrap->getWorldTransform();
	//todo: check which transform is needed here

	btScalar dst =
		m_worldInfo->m_sparsesdf.Evaluate(
		wtr.invXform(x),
		shp,
		nrm,
		margin);
	if (dst < 0)
	{
		cti.m_colObj = colObjWrap->getCollisionObject();
		cti.m_normal = wtr.getBasis() * nrm;
		cti.m_offset = -btDot( cti.m_normal, x - cti.m_normal * dst );
		return (true);
	}
	return (false);
}

//
void                    btSoftBody::updateNormals()
{

	const btVector3 zv(0, 0, 0);
	int i, ni;

	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		m_nodes[i].m_n = zv;
	}
	for (i = 0, ni = m_faces.size(); i < ni; ++i)
	{
		btSoftBody::Face   &f = m_faces[i];
		const btVector3     n = btCross(f.m_n[1]->m_x - f.m_n[0]->m_x,
			f.m_n[2]->m_x - f.m_n[0]->m_x);
		f.m_normal = n.normalized();
		f.m_n[0]->m_n += n;
		f.m_n[1]->m_n += n;
		f.m_n[2]->m_n += n;
	}
	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		btScalar len = m_nodes[i].m_n.length();
		if (len > SIMD_EPSILON)
			m_nodes[i].m_n /= len;
	}
}

//
void                    btSoftBody::updateBounds()
{
	/*if( m_acceleratedSoftBody )
	{
	// If we have an accelerated softbody we need to obtain the bounds correctly
	// For now (slightly hackily) just have a very large AABB
	// TODO: Write get bounds kernel
	// If that is updating in place, atomic collisions might be low (when the cloth isn't perfectly aligned to an axis) and we could
	// probably do a test and exchange reasonably efficiently.

	m_bounds[0] = btVector3(-1000, -1000, -1000);
	m_bounds[1] = btVector3(1000, 1000, 1000);

	} else {*/
	if (m_ndbvt.m_root)
	{
		const btVector3    &mins = m_ndbvt.m_root->volume.Mins();
		const btVector3    &maxs = m_ndbvt.m_root->volume.Maxs();
		const btScalar      csm = getCollisionShape()->getMargin();
		const btVector3     mrg = btVector3(  csm,
			csm,
			csm) * 1; // ??? to investigate...
		m_bounds[0] = mins - mrg;
		m_bounds[1] = maxs + mrg;
		if (0 != getBroadphaseHandle())
		{
			m_worldInfo->m_broadphase->setAabb( getBroadphaseHandle(),
				m_bounds[0],
				m_bounds[1],
				m_worldInfo->m_dispatcher);
		}
	}
	else
	{
		m_bounds[0] =
			m_bounds[1] = btVector3(0, 0, 0);
	}
	//}
}


//
void                    btSoftBody::updatePose()
{
	if (m_pose.m_bframe)
	{
		btSoftBody::Pose   &pose = m_pose;
		const btVector3     com = evaluateCom();
		/* Com          */
		pose.m_com  =   com;
		/* Rotation     */
		btMatrix3x3     Apq;
		const btScalar  eps = SIMD_EPSILON;
		Apq[0] = Apq[1] = Apq[2] = btVector3(0, 0, 0);
		Apq[0].setX(eps); Apq[1].setY(eps * 2); Apq[2].setZ(eps * 3);
		for (int i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			const btVector3     a = pose.m_wgh[i] * (m_nodes[i].m_x - com);
			const btVector3    &b = pose.m_pos[i];
			Apq[0] += a.x() * b;
			Apq[1] += a.y() * b;
			Apq[2] += a.z() * b;
		}
		btMatrix3x3     r, s;
		PolarDecompose(Apq, r, s);
		pose.m_rot = r;
		pose.m_scl = pose.m_aqq * r.transpose() * Apq;
		if (m_cfg.maxvolume > 1)
		{
			const btScalar  idet = Clamp<btScalar>(   1 / pose.m_scl.determinant(),
				1, m_cfg.maxvolume);
			pose.m_scl = Mul(pose.m_scl, idet);
		}

	}
}

//
void                btSoftBody::updateArea(bool averageArea)
{
	int i, ni;

	/* Face area        */
	for (i = 0, ni = m_faces.size(); i < ni; ++i)
	{
		Face       &f = m_faces[i];
		f.m_ra  =   AreaOf(f.m_n[0]->m_x, f.m_n[1]->m_x, f.m_n[2]->m_x);
	}

	/* Node area        */

	if (averageArea)
	{
		btAlignedObjectArray<int>   counts;
		counts.resize(m_nodes.size(), 0);
		for (i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			m_nodes[i].m_area   =   0;
		}
		for (i = 0, ni = m_faces.size(); i < ni; ++i)
		{
			btSoftBody::Face   &f = m_faces[i];
			for (int j = 0; j < 3; ++j)
			{
				const int index = (int)(f.m_n[j] - &m_nodes[0]);
				counts[index]++;
				f.m_n[j]->m_area += btFabs(f.m_ra);
			}
		}
		for (i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			if (counts[i] > 0)
				m_nodes[i].m_area /= (btScalar)counts[i];
			else
				m_nodes[i].m_area = 0;
		}
	}
	else
	{
		// initialize node area as zero
		for (i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			m_nodes[i].m_area = 0;
		}

		for (i = 0, ni = m_faces.size(); i < ni; ++i)
		{
			btSoftBody::Face   &f = m_faces[i];

			for (int j = 0; j < 3; ++j)
			{
				f.m_n[j]->m_area += f.m_ra;
			}
		}

		for (i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			m_nodes[i].m_area *= 0.3333333f;
		}
	}
}


void                btSoftBody::updateLinkConstants()
{
	int i, ni;

	/* Links        */
	for (i = 0, ni = m_links.size(); i < ni; ++i)
	{
		Link       &l = m_links[i];
		Material   &m = *l.m_material;
		l.m_c0  =   (l.m_n[0]->m_im + l.m_n[1]->m_im) / m.m_kLST;
	}
}

void                btSoftBody::updateConstants()
{
	resetLinkRestLengths();
	updateLinkConstants();
	updateArea();
}



//
void                    btSoftBody::initializeClusters()
{
	int i;

	for ( i = 0; i < m_clusters.size(); ++i)
	{
		Cluster    &c = *m_clusters[i];
		c.m_imass = 0;
		c.m_masses.resize(c.m_nodes.size());
		for (int j = 0; j < c.m_nodes.size(); ++j)
		{
			if (c.m_nodes[j]->m_im == 0)
			{
				c.m_containsAnchor = true;
				c.m_masses[j]   =   BT_LARGE_FLOAT;
			}
			else
			{
				c.m_masses[j]   =   btScalar(1.) / c.m_nodes[j]->m_im;
			}
			c.m_imass       +=  c.m_masses[j];
		}
		c.m_imass       =   btScalar(1.) / c.m_imass;
		c.m_com         =   btSoftBody::clusterCom(&c);
		c.m_lv          =   btVector3(0, 0, 0);
		c.m_av          =   btVector3(0, 0, 0);
		c.m_leaf        =   0;
		/* Inertia  */
		btMatrix3x3    &ii = c.m_locii;
		ii[0] = ii[1] = ii[2] = btVector3(0, 0, 0);
		{
			int i, ni;

			for (i = 0, ni = c.m_nodes.size(); i < ni; ++i)
			{
				const btVector3 k = c.m_nodes[i]->m_x - c.m_com;
				const btVector3 q = k * k;
				const btScalar  m = c.m_masses[i];
				ii[0][0]    +=  m * (q[1] + q[2]);
				ii[1][1]    +=  m * (q[0] + q[2]);
				ii[2][2]    +=  m * (q[0] + q[1]);
				ii[0][1]    -=  m * k[0] * k[1];
				ii[0][2]    -=  m * k[0] * k[2];
				ii[1][2]    -=  m * k[1] * k[2];
			}
		}
		ii[1][0] = ii[0][1];
		ii[2][0] = ii[0][2];
		ii[2][1] = ii[1][2];

		ii = ii.inverse();

		/* Frame    */
		c.m_framexform.setIdentity();
		c.m_framexform.setOrigin(c.m_com);
		c.m_framerefs.resize(c.m_nodes.size());
		{
			int i;
			for (i = 0; i < c.m_framerefs.size(); ++i)
			{
				c.m_framerefs[i] = c.m_nodes[i]->m_x - c.m_com;
			}
		}
	}
}

//
void                    btSoftBody::updateClusters()
{
	BT_PROFILE("UpdateClusters");
	int i;

	for (i = 0; i < m_clusters.size(); ++i)
	{
		btSoftBody::Cluster    &c = *m_clusters[i];
		const int               n = c.m_nodes.size();
		//const btScalar            invn=1/(btScalar)n;
		if (n)
		{
			/* Frame                */
			const btScalar  eps = btScalar(0.0001);
			btMatrix3x3     m, r, s;
			m[0] = m[1] = m[2] = btVector3(0, 0, 0);
			m[0][0] = eps * 1;
			m[1][1] = eps * 2;
			m[2][2] = eps * 3;
			c.m_com = clusterCom(&c);
			for (int i = 0; i < c.m_nodes.size(); ++i)
			{
				const btVector3     a = c.m_nodes[i]->m_x - c.m_com;
				const btVector3    &b = c.m_framerefs[i];
				m[0] += a[0] * b; m[1] += a[1] * b; m[2] += a[2] * b;
			}
			PolarDecompose(m, r, s);
			c.m_framexform.setOrigin(c.m_com);
			c.m_framexform.setBasis(r);
			/* Inertia          */
#if 1/* Constant    */
			c.m_invwi = c.m_framexform.getBasis() * c.m_locii * c.m_framexform.getBasis().transpose();
#else
#if 0/* Sphere  */
			const btScalar  rk = (2 * c.m_extents.length2()) / (5 * c.m_imass);
			const btVector3 inertia(rk, rk, rk);
			const btVector3 iin(btFabs(inertia[0]) > SIMD_EPSILON ? 1 / inertia[0] : 0,
				btFabs(inertia[1]) > SIMD_EPSILON ? 1 / inertia[1] : 0,
				btFabs(inertia[2]) > SIMD_EPSILON ? 1 / inertia[2] : 0);

			c.m_invwi = c.m_xform.getBasis().scaled(iin) * c.m_xform.getBasis().transpose();
#else/* Actual  */
			c.m_invwi[0] = c.m_invwi[1] = c.m_invwi[2] = btVector3(0, 0, 0);
			for (int i = 0; i < n; ++i)
			{
				const btVector3 k = c.m_nodes[i]->m_x - c.m_com;
				const btVector3     q = k * k;
				const btScalar      m = 1 / c.m_nodes[i]->m_im;
				c.m_invwi[0][0] +=  m * (q[1] + q[2]);
				c.m_invwi[1][1] +=  m * (q[0] + q[2]);
				c.m_invwi[2][2] +=  m * (q[0] + q[1]);
				c.m_invwi[0][1] -=  m * k[0] * k[1];
				c.m_invwi[0][2] -=  m * k[0] * k[2];
				c.m_invwi[1][2] -=  m * k[1] * k[2];
			}
			c.m_invwi[1][0] = c.m_invwi[0][1];
			c.m_invwi[2][0] = c.m_invwi[0][2];
			c.m_invwi[2][1] = c.m_invwi[1][2];
			c.m_invwi = c.m_invwi.inverse();
#endif
#endif
			/* Velocities           */
			c.m_lv = btVector3(0, 0, 0);
			c.m_av = btVector3(0, 0, 0);
			{
				int i;

				for (i = 0; i < n; ++i)
				{
					const btVector3 v = c.m_nodes[i]->m_v * c.m_masses[i];
					c.m_lv  +=  v;
					c.m_av  +=  btCross(c.m_nodes[i]->m_x - c.m_com, v);
				}
			}
			c.m_lv = c.m_imass * c.m_lv * (1 - c.m_ldamping);
			c.m_av = c.m_invwi * c.m_av * (1 - c.m_adamping);
			c.m_vimpulses[0]    =
				c.m_vimpulses[1]    = btVector3(0, 0, 0);
			c.m_dimpulses[0]    =
				c.m_dimpulses[1]    = btVector3(0, 0, 0);
			c.m_nvimpulses      = 0;
			c.m_ndimpulses      = 0;
			/* Matching             */
			if (c.m_matching > 0)
			{
				for (int j = 0; j < c.m_nodes.size(); ++j)
				{
					Node           &n = *c.m_nodes[j];
					const btVector3 x = c.m_framexform * c.m_framerefs[j];
					n.m_x = Lerp(n.m_x, x, c.m_matching);
				}
			}
			/* Dbvt                 */
			if (c.m_collide)
			{
				btVector3   mi = c.m_nodes[0]->m_x;
				btVector3   mx = mi;
				for (int j = 1; j < n; ++j)
				{
					mi.setMin(c.m_nodes[j]->m_x);
					mx.setMax(c.m_nodes[j]->m_x);
				}
				ATTRIBUTE_ALIGNED16(btDbvtVolume)   bounds = btDbvtVolume::FromMM(mi, mx);
				if (c.m_leaf)
					m_cdbvt.update(c.m_leaf, bounds, c.m_lv * m_sst.sdt * 3, m_sst.radmrg);
				else
					c.m_leaf = m_cdbvt.insert(bounds, &c);
			}
		}
	}


}




//
void                    btSoftBody::cleanupClusters()
{
	for (int i = 0; i < m_joints.size(); ++i)
	{
		m_joints[i]->Terminate(m_sst.sdt);
		if (m_joints[i]->m_delete)
		{
			btAlignedFree(m_joints[i]);
			m_joints.remove(m_joints[i--]);
		}
	}
}

//
void                    btSoftBody::prepareClusters(int iterations)
{
	for (int i = 0; i < m_joints.size(); ++i)
	{
		m_joints[i]->Prepare(m_sst.sdt, iterations);
	}
}


//
void                    btSoftBody::solveClusters(btScalar sor)
{
	for (int i = 0, ni = m_joints.size(); i < ni; ++i)
	{
		m_joints[i]->Solve(m_sst.sdt, sor);
	}
}

//
void                    btSoftBody::applyClusters(bool drift)
{
	BT_PROFILE("ApplyClusters");
	//  const btScalar                  f0=m_sst.sdt;
	//const btScalar                    f1=f0/2;
	btAlignedObjectArray<btVector3> deltas;
	btAlignedObjectArray<btScalar> weights;
	deltas.resize(m_nodes.size(), btVector3(0, 0, 0));
	weights.resize(m_nodes.size(), 0);
	int i;

	if (drift)
	{
		for (i = 0; i < m_clusters.size(); ++i)
		{
			Cluster    &c = *m_clusters[i];
			if (c.m_ndimpulses)
			{
				c.m_dimpulses[0] /= (btScalar)c.m_ndimpulses;
				c.m_dimpulses[1] /= (btScalar)c.m_ndimpulses;
			}
		}
	}

	for (i = 0; i < m_clusters.size(); ++i)
	{
		Cluster    &c = *m_clusters[i];
		if (0 < (drift ? c.m_ndimpulses : c.m_nvimpulses))
		{
			const btVector3     v = (drift ? c.m_dimpulses[0] : c.m_vimpulses[0]) * m_sst.sdt;
			const btVector3     w = (drift ? c.m_dimpulses[1] : c.m_vimpulses[1]) * m_sst.sdt;
			for (int j = 0; j < c.m_nodes.size(); ++j)
			{
				const int           idx = int(c.m_nodes[j] - &m_nodes[0]);
				const btVector3    &x = c.m_nodes[j]->m_x;
				const btScalar      q = c.m_masses[j];
				deltas[idx]     +=  (v + btCross(w, x - c.m_com)) * q;
				weights[idx]    +=  q;
			}
		}
	}
	for (i = 0; i < deltas.size(); ++i)
	{
		if (weights[i] > 0)
		{
			m_nodes[i].m_x += deltas[i] / weights[i];
		}
	}
}

//
void                    btSoftBody::dampClusters()
{
	int i;

	for (i = 0; i < m_clusters.size(); ++i)
	{
		Cluster    &c = *m_clusters[i];
		if (c.m_ndamping > 0)
		{
			for (int j = 0; j < c.m_nodes.size(); ++j)
			{
				Node           &n = *c.m_nodes[j];
				if (n.m_im > 0)
				{
					const btVector3 vx = c.m_lv + btCross(c.m_av, c.m_nodes[j]->m_q - c.m_com);
					if (vx.length2() <= n.m_v.length2())
					{
						n.m_v   +=  c.m_ndamping * (vx - n.m_v);
					}
				}
			}
		}
	}
}

//
void                btSoftBody::Joint::Prepare(btScalar dt, int)
{
	m_bodies[0].activate();
	m_bodies[1].activate();
}

//
void                btSoftBody::LJoint::Prepare(btScalar dt, int iterations)
{
	static const btScalar   maxdrift = 4;
	Joint::Prepare(dt, iterations);
	m_rpos[0]       =   m_bodies[0].xform() * m_refs[0];
	m_rpos[1]       =   m_bodies[1].xform() * m_refs[1];
	m_drift         =   Clamp(m_rpos[0] - m_rpos[1], maxdrift) * m_erp / dt;
	m_rpos[0]       -=  m_bodies[0].xform().getOrigin();
	m_rpos[1]       -=  m_bodies[1].xform().getOrigin();
	m_massmatrix    =   ImpulseMatrix(  m_bodies[0].invMass(), m_bodies[0].invWorldInertia(), m_rpos[0],
		m_bodies[1].invMass(), m_bodies[1].invWorldInertia(), m_rpos[1]);
	if (m_split > 0)
	{
		m_sdrift    =   m_massmatrix * (m_drift * m_split);
		m_drift     *=  1 - m_split;
	}
	m_drift /= (btScalar)iterations;
}

//
void                btSoftBody::LJoint::Solve(btScalar dt, btScalar sor)
{
	const btVector3     va = m_bodies[0].velocity(m_rpos[0]);
	const btVector3     vb = m_bodies[1].velocity(m_rpos[1]);
	const btVector3     vr = va - vb;
	btSoftBody::Impulse impulse;
	impulse.m_asVelocity    =   1;
	impulse.m_velocity      =   m_massmatrix * (m_drift + vr * m_cfm) * sor;
	m_bodies[0].applyImpulse(-impulse, m_rpos[0]);
	m_bodies[1].applyImpulse( impulse, m_rpos[1]);
}

//
void                btSoftBody::LJoint::Terminate(btScalar dt)
{
	if (m_split > 0)
	{
		m_bodies[0].applyDImpulse(-m_sdrift, m_rpos[0]);
		m_bodies[1].applyDImpulse( m_sdrift, m_rpos[1]);
	}
}

//
void                btSoftBody::AJoint::Prepare(btScalar dt, int iterations)
{
	static const btScalar   maxdrift = SIMD_PI / 16;
	m_icontrol->Prepare(this);
	Joint::Prepare(dt, iterations);
	m_axis[0]   =   m_bodies[0].xform().getBasis() * m_refs[0];
	m_axis[1]   =   m_bodies[1].xform().getBasis() * m_refs[1];
	m_drift     =   NormalizeAny(btCross(m_axis[1], m_axis[0]));
	m_drift     *=  btMin(maxdrift, btAcos(Clamp<btScalar>(btDot(m_axis[0], m_axis[1]), -1, +1)));
	m_drift     *=  m_erp / dt;
	m_massmatrix =   AngularImpulseMatrix(m_bodies[0].invWorldInertia(), m_bodies[1].invWorldInertia());
	if (m_split > 0)
	{
		m_sdrift    =   m_massmatrix * (m_drift * m_split);
		m_drift     *=  1 - m_split;
	}
	m_drift /= (btScalar)iterations;
}

//
void                btSoftBody::AJoint::Solve(btScalar dt, btScalar sor)
{
	const btVector3     va = m_bodies[0].angularVelocity();
	const btVector3     vb = m_bodies[1].angularVelocity();
	const btVector3     vr = va - vb;
	const btScalar      sp = btDot(vr, m_axis[0]);
	const btVector3     vc = vr - m_axis[0] * m_icontrol->Speed(this, sp);
	btSoftBody::Impulse impulse;
	impulse.m_asVelocity    =   1;
	impulse.m_velocity      =   m_massmatrix * (m_drift + vc * m_cfm) * sor;
	m_bodies[0].applyAImpulse(-impulse);
	m_bodies[1].applyAImpulse( impulse);
}

//
void                btSoftBody::AJoint::Terminate(btScalar dt)
{
	if (m_split > 0)
	{
		m_bodies[0].applyDAImpulse(-m_sdrift);
		m_bodies[1].applyDAImpulse( m_sdrift);
	}
}

//
void                btSoftBody::CJoint::Prepare(btScalar dt, int iterations)
{
	Joint::Prepare(dt, iterations);
	const bool  dodrift = (m_life == 0);
	m_delete = (++m_life) > m_maxlife;
	if (dodrift)
	{
		m_drift = m_drift * m_erp / dt;
		if (m_split > 0)
		{
			m_sdrift    =   m_massmatrix * (m_drift * m_split);
			m_drift     *=  1 - m_split;
		}
		m_drift /= (btScalar)iterations;
	}
	else
	{
		m_drift = m_sdrift = btVector3(0, 0, 0);
	}
}

//
void                btSoftBody::CJoint::Solve(btScalar dt, btScalar sor)
{
	const btVector3     va = m_bodies[0].velocity(m_rpos[0]);
	const btVector3     vb = m_bodies[1].velocity(m_rpos[1]);
	const btVector3     vrel = va - vb;
	const btScalar      rvac = btDot(vrel, m_normal);
	btSoftBody::Impulse impulse;
	impulse.m_asVelocity    =   1;
	impulse.m_velocity      =   m_drift;
	if (rvac < 0)
	{
		const btVector3 iv = m_normal * rvac;
		const btVector3 fv = vrel - iv;
		impulse.m_velocity  +=  iv + fv * m_friction;
	}
	impulse.m_velocity = m_massmatrix * impulse.m_velocity * sor;

	if (m_bodies[0].m_soft == m_bodies[1].m_soft)
	{
		if ((impulse.m_velocity.getX() == impulse.m_velocity.getX()) && (impulse.m_velocity.getY() == impulse.m_velocity.getY()) &&
			(impulse.m_velocity.getZ() == impulse.m_velocity.getZ()))
		{
			if (impulse.m_asVelocity)
			{
				if (impulse.m_velocity.length() < m_bodies[0].m_soft->m_maxSelfCollisionImpulse)
				{

				}
				else
				{
					m_bodies[0].applyImpulse(-impulse * m_bodies[0].m_soft->m_selfCollisionImpulseFactor, m_rpos[0]);
					m_bodies[1].applyImpulse( impulse * m_bodies[0].m_soft->m_selfCollisionImpulseFactor, m_rpos[1]);
				}
			}
		}
	}
	else
	{
		m_bodies[0].applyImpulse(-impulse, m_rpos[0]);
		m_bodies[1].applyImpulse( impulse, m_rpos[1]);
	}
}

//
void                btSoftBody::CJoint::Terminate(btScalar dt)
{
	if (m_split > 0)
	{
		m_bodies[0].applyDImpulse(-m_sdrift, m_rpos[0]);
		m_bodies[1].applyDImpulse( m_sdrift, m_rpos[1]);
	}
}

//
void                btSoftBody::applyForces()
{

	BT_PROFILE("SoftBody applyForces");
	//  const btScalar                  dt =            m_sst.sdt;
	const btScalar                  kLF =           m_cfg.kLF;
	const btScalar                  kDG =           m_cfg.kDG;
	const btScalar                  kPR =           m_cfg.kPR;
	const btScalar                  kVC =           m_cfg.kVC;
	const bool                      as_lift =       kLF > 0;
	const bool                      as_drag =       kDG > 0;
	const bool                      as_pressure =   kPR != 0;
	const bool                      as_volume =     kVC > 0;
	const bool                      as_aero =       as_lift ||
		as_drag     ;
	//const bool                        as_vaero =      as_aero &&
	//                                              (m_cfg.aeromodel < btSoftBody::eAeroModel::F_TwoSided);
	//const bool                        as_faero =      as_aero &&
	//                                              (m_cfg.aeromodel >= btSoftBody::eAeroModel::F_TwoSided);
	const bool                      use_medium =    as_aero;
	const bool                      use_volume =    as_pressure ||
		as_volume   ;
	btScalar                        volume =        0;
	btScalar                        ivolumetp =     0;
	btScalar                        dvolumetv =     0;
	btSoftBody::sMedium medium;
	if (use_volume)
	{
		volume      =   getVolume();
		ivolumetp   =   1 / btFabs(volume) * kPR;
		dvolumetv   =   (m_pose.m_volume - volume) * kVC;
	}
	/* Per vertex forces            */
	int i, ni;

	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		btSoftBody::Node   &n = m_nodes[i];
		if (n.m_im > 0)
		{
			if (use_medium)
			{
				/* Aerodynamics         */
				addAeroForceToNode(m_windVelocity, i);
			}
			/* Pressure             */
			if (as_pressure)
			{
				n.m_f   +=  n.m_n * (n.m_area * ivolumetp);
			}
			/* Volume               */
			if (as_volume)
			{
				n.m_f   +=  n.m_n * (n.m_area * dvolumetv);
			}
		}
	}

	/* Per face forces              */
	for (i = 0, ni = m_faces.size(); i < ni; ++i)
	{
		//  btSoftBody::Face&   f=m_faces[i];

		/* Aerodynamics         */
		addAeroForceToFace(m_windVelocity, i);
	}
}

//
void                btSoftBody::PSolve_Anchors(btSoftBody *psb, btScalar kst, btScalar ti)
{
	const btScalar  kAHR = psb->m_cfg.kAHR * kst;
	const btScalar  dt = psb->m_sst.sdt;
	for (int i = 0, ni = psb->m_anchors.size(); i < ni; ++i)
	{
		const Anchor       &a = psb->m_anchors[i];
		const btTransform  &t = a.m_body->getWorldTransform();
		Node               &n = *a.m_node;
		const btVector3     wa = t * a.m_local;
		const btVector3     va = a.m_body->getVelocityInLocalPoint(a.m_c1) * dt;
		const btVector3     vb = n.m_x - n.m_q;
		const btVector3     vr = (va - vb) + (wa - n.m_x) * kAHR;
		const btVector3     impulse = a.m_c0 * vr * a.m_influence;
		n.m_x += impulse * a.m_c2;
		a.m_body->applyImpulse(-impulse, a.m_c1);
	}
}

//
void btSoftBody::PSolve_RContacts(btSoftBody *psb, btScalar kst, btScalar ti)
{
	const btScalar  dt = psb->m_sst.sdt;
	const btScalar  mrg = psb->getCollisionShape()->getMargin();
	for (int i = 0, ni = psb->m_rcontacts.size(); i < ni; ++i)
	{
		const RContact     &c = psb->m_rcontacts[i];
		const sCti         &cti = c.m_cti;
		btRigidBody *tmpRigid = (btRigidBody *)btRigidBody::upcast(cti.m_colObj);

		const btVector3     va = tmpRigid ? tmpRigid->getVelocityInLocalPoint(c.m_c1) * dt : btVector3(0, 0, 0);
		const btVector3     vb = c.m_node->m_x - c.m_node->m_q;
		const btVector3     vr = vb - va;
		const btScalar      dn = btDot(vr, cti.m_normal);
		if (dn <= SIMD_EPSILON)
		{
			const btScalar      dp = btMin( (btDot(c.m_node->m_x, cti.m_normal) + cti.m_offset), mrg );
			const btVector3     fv = vr - (cti.m_normal * dn);
			// c0 is the impulse matrix, c3 is 1 - the friction coefficient or 0, c4 is the contact hardness coefficient
			const btVector3     impulse = c.m_c0 * ( (vr - (fv * c.m_c3) + (cti.m_normal * (dp * c.m_c4))) * kst );
			c.m_node->m_x -= impulse * c.m_c2;
			if (tmpRigid)
				tmpRigid->applyImpulse(impulse, c.m_c1);
		}
	}
}

//
void                btSoftBody::PSolve_SContacts(btSoftBody *psb, btScalar, btScalar ti)
{
	for (int i = 0, ni = psb->m_scontacts.size(); i < ni; ++i)
	{
		const SContact     &c = psb->m_scontacts[i];
		const btVector3    &nr = c.m_normal;
		Node               &n = *c.m_node;
		Face               &f = *c.m_face;
		const btVector3     p = BaryEval( f.m_n[0]->m_x,
			f.m_n[1]->m_x,
			f.m_n[2]->m_x,
			c.m_weights);
		const btVector3     q = BaryEval( f.m_n[0]->m_q,
			f.m_n[1]->m_q,
			f.m_n[2]->m_q,
			c.m_weights);
		const btVector3     vr = (n.m_x - n.m_q) - (p - q);
		btVector3           corr(0, 0, 0);
		btScalar dot = btDot(vr, nr);
		if (dot < 0)
		{
			const btScalar  j = c.m_margin - (btDot(nr, n.m_x) - btDot(nr, p));
			corr += c.m_normal * j;
		}
		corr            -=  ProjectOnPlane(vr, nr) * c.m_friction;
		n.m_x           +=  corr * c.m_cfm[0];
		f.m_n[0]->m_x   -=  corr * (c.m_cfm[1] * c.m_weights.x());
		f.m_n[1]->m_x   -=  corr * (c.m_cfm[1] * c.m_weights.y());
		f.m_n[2]->m_x   -=  corr * (c.m_cfm[1] * c.m_weights.z());
	}
}

//
void                btSoftBody::PSolve_Links(btSoftBody *psb, btScalar kst, btScalar ti)
{
	for (int i = 0, ni = psb->m_links.size(); i < ni; ++i)
	{
		Link   &l = psb->m_links[i];
		if (l.m_c0 > 0)
		{
			Node           &a = *l.m_n[0];
			Node           &b = *l.m_n[1];
			const btVector3 del = b.m_x - a.m_x;
			const btScalar  len = del.length2();
			if (l.m_c1 + len > SIMD_EPSILON)
			{
				const btScalar  k = ((l.m_c1 - len) / (l.m_c0 * (l.m_c1 + len))) * kst;
				a.m_x -= del * (k * a.m_im);
				b.m_x += del * (k * b.m_im);
			}
		}
	}
}

//
void                btSoftBody::VSolve_Links(btSoftBody *psb, btScalar kst)
{
	for (int i = 0, ni = psb->m_links.size(); i < ni; ++i)
	{
		Link           &l = psb->m_links[i];
		Node          **n = l.m_n;
		const btScalar  j = -btDot(l.m_c3, n[0]->m_v - n[1]->m_v) * l.m_c2 * kst;
		n[0]->m_v += l.m_c3 * (j * n[0]->m_im);
		n[1]->m_v -= l.m_c3 * (j * n[1]->m_im);
	}
}

//
btSoftBody::psolver_t   btSoftBody::getSolver(ePSolver::_ solver)
{
	switch (solver)
	{
	case    ePSolver::Anchors:
		return (&btSoftBody::PSolve_Anchors);
	case    ePSolver::Linear:
		return (&btSoftBody::PSolve_Links);
	case    ePSolver::RContacts:
		return (&btSoftBody::PSolve_RContacts);
	case    ePSolver::SContacts:
		return (&btSoftBody::PSolve_SContacts);
	default:
		{
		}
	}
	return (0);
}

//
btSoftBody::vsolver_t   btSoftBody::getSolver(eVSolver::_ solver)
{
	switch (solver)
	{
	case    eVSolver::Linear:       return (&btSoftBody::VSolve_Links);
	default:
		{
		}
	}
	return (0);
}

//
void            btSoftBody::defaultCollisionHandler(const btCollisionObjectWrapper *pcoWrap)
{

	switch (m_cfg.collisions & fCollision::RVSmask)
	{
	case    fCollision::SDF_RS:
		{
			btSoftColliders::CollideSDF_RS  docollide;
			btRigidBody        *prb1 = (btRigidBody *) btRigidBody::upcast(pcoWrap->getCollisionObject());
			btTransform wtr = pcoWrap->getWorldTransform();

			const btTransform   ctr = pcoWrap->getWorldTransform();
			const btScalar      timemargin = (wtr.getOrigin() - ctr.getOrigin()).length();
			const btScalar      basemargin = getCollisionShape()->getMargin();
			btVector3           mins;
			btVector3           maxs;
			ATTRIBUTE_ALIGNED16(btDbvtVolume)       volume;
			pcoWrap->getCollisionShape()->getAabb(  pcoWrap->getWorldTransform(),
				mins,
				maxs);
			volume = btDbvtVolume::FromMM(mins, maxs);
			volume.Expand(btVector3(basemargin, basemargin, basemargin));
			docollide.psb       =   this;
			docollide.m_colObj1Wrap = pcoWrap;
			docollide.m_rigidBody = prb1;

			docollide.dynmargin =   basemargin + timemargin;
			docollide.stamargin =   basemargin;
			m_ndbvt.collideTV(m_ndbvt.m_root, volume, docollide);
		}
		break;
	case    fCollision::CL_RS:
		{
			btSoftColliders::CollideCL_RS   collider;
			collider.ProcessColObj(this, pcoWrap);
		}
		break;
	}
}

//
void            btSoftBody::defaultCollisionHandler(btSoftBody *psb)
{
	const int cf = m_cfg.collisions & psb->m_cfg.collisions;
	switch (cf & fCollision::SVSmask)
	{
	case    fCollision::CL_SS:
		{

			//support self-collision if CL_SELF flag set
			if (this != psb || psb->m_cfg.collisions & fCollision::CL_SELF)
			{
				btSoftColliders::CollideCL_SS   docollide;
				docollide.ProcessSoftSoft(this, psb);
			}

		}
		break;
	case    fCollision::VF_SS:
		{
			//only self-collision for Cluster, not Vertex-Face yet
			if (this != psb)
			{
				btSoftColliders::CollideVF_SS   docollide;
				/* common                   */
				docollide.mrg =  getCollisionShape()->getMargin() +
					psb->getCollisionShape()->getMargin();
				/* psb0 nodes vs psb1 faces */
				docollide.psb[0] = this;
				docollide.psb[1] = psb;
				docollide.psb[0]->m_ndbvt.collideTT(    docollide.psb[0]->m_ndbvt.m_root,
					docollide.psb[1]->m_fdbvt.m_root,
					docollide);
				/* psb1 nodes vs psb0 faces */
				docollide.psb[0] = psb;
				docollide.psb[1] = this;
				docollide.psb[0]->m_ndbvt.collideTT(    docollide.psb[0]->m_ndbvt.m_root,
					docollide.psb[1]->m_fdbvt.m_root,
					docollide);
			}
		}
		break;
	default:
		{

		}
	}
}



void btSoftBody::setWindVelocity( const btVector3 &velocity )
{
	m_windVelocity = velocity;
}


const btVector3 &btSoftBody::getWindVelocity()
{
	return m_windVelocity;
}



int btSoftBody::calculateSerializeBufferSize()  const
{
	int sz = sizeof(btSoftBodyData);
	return sz;
}

///fills the dataBuffer and returns the struct name (and 0 on failure)
const char *btSoftBody::serialize(void *dataBuffer, class btSerializer *serializer) const
{
	btSoftBodyData *sbd = (btSoftBodyData *) dataBuffer;

	btCollisionObject::serialize(&sbd->m_collisionObjectData, serializer);

	btHashMap<btHashPtr, int>    m_nodeIndexMap;

	sbd->m_numMaterials = m_materials.size();
	sbd->m_materials = sbd->m_numMaterials ? (SoftBodyMaterialData **) serializer->getUniquePointer((void *)&m_materials) : 0;

	if (sbd->m_materials)
	{
		int sz = sizeof(SoftBodyMaterialData *);
		int numElem = sbd->m_numMaterials;
		btChunk *chunk = serializer->allocate(sz, numElem);
		//SoftBodyMaterialData** memPtr = chunk->m_oldPtr;
		SoftBodyMaterialData **memPtr = (SoftBodyMaterialData **)chunk->m_oldPtr;
		for (int i = 0; i < numElem; i++, memPtr++)
		{
			btSoftBody::Material *mat = m_materials[i];
			*memPtr = mat ? (SoftBodyMaterialData *)serializer->getUniquePointer((void *)mat) : 0;
			if (!serializer->findPointer(mat))
			{
				//serialize it here
				btChunk *chunk = serializer->allocate(sizeof(SoftBodyMaterialData), 1);
				SoftBodyMaterialData *memPtr = (SoftBodyMaterialData *)chunk->m_oldPtr;
				memPtr->m_flags = mat->m_flags;
				memPtr->m_angularStiffness = mat->m_kAST;
				memPtr->m_linearStiffness = mat->m_kLST;
				memPtr->m_volumeStiffness = mat->m_kVST;
				serializer->finalizeChunk(chunk, "SoftBodyMaterialData", BT_SBMATERIAL_CODE, mat);
			}
		}
		serializer->finalizeChunk(chunk, "SoftBodyMaterialData", BT_ARRAY_CODE, (void *) &m_materials);
	}




	sbd->m_numNodes = m_nodes.size();
	sbd->m_nodes = sbd->m_numNodes ? (SoftBodyNodeData *)serializer->getUniquePointer((void *)&m_nodes) : 0;
	if (sbd->m_nodes)
	{
		int sz = sizeof(SoftBodyNodeData);
		int numElem = sbd->m_numNodes;
		btChunk *chunk = serializer->allocate(sz, numElem);
		SoftBodyNodeData *memPtr = (SoftBodyNodeData *)chunk->m_oldPtr;
		for (int i = 0; i < numElem; i++, memPtr++)
		{
			m_nodes[i].m_f.serializeFloat( memPtr->m_accumulatedForce);
			memPtr->m_area = m_nodes[i].m_area;
			memPtr->m_attach = m_nodes[i].m_battach;
			memPtr->m_inverseMass = m_nodes[i].m_im;
			memPtr->m_material = m_nodes[i].m_material ? (SoftBodyMaterialData *)serializer->getUniquePointer((void *) m_nodes[i].m_material) : 0;
			m_nodes[i].m_n.serializeFloat(memPtr->m_normal);
			m_nodes[i].m_x.serializeFloat(memPtr->m_position);
			m_nodes[i].m_q.serializeFloat(memPtr->m_previousPosition);
			m_nodes[i].m_v.serializeFloat(memPtr->m_velocity);
			m_nodeIndexMap.insert(&m_nodes[i], i);
		}
		serializer->finalizeChunk(chunk, "SoftBodyNodeData", BT_SBNODE_CODE, (void *) &m_nodes);
	}

	sbd->m_numLinks = m_links.size();
	sbd->m_links = sbd->m_numLinks ? (SoftBodyLinkData *) serializer->getUniquePointer((void *)&m_links[0]) : 0;
	if (sbd->m_links)
	{
		int sz = sizeof(SoftBodyLinkData);
		int numElem = sbd->m_numLinks;
		btChunk *chunk = serializer->allocate(sz, numElem);
		SoftBodyLinkData *memPtr = (SoftBodyLinkData *)chunk->m_oldPtr;
		for (int i = 0; i < numElem; i++, memPtr++)
		{
			memPtr->m_bbending = m_links[i].m_bbending;
			memPtr->m_material = m_links[i].m_material ? (SoftBodyMaterialData *)serializer->getUniquePointer((void *) m_links[i].m_material) : 0;
			memPtr->m_nodeIndices[0] = m_links[i].m_n[0] ? m_links[i].m_n[0] - &m_nodes[0] : -1;
			memPtr->m_nodeIndices[1] = m_links[i].m_n[1] ? m_links[i].m_n[1] - &m_nodes[0] : -1;
			btAssert(memPtr->m_nodeIndices[0] < m_nodes.size());
			btAssert(memPtr->m_nodeIndices[1] < m_nodes.size());
			memPtr->m_restLength = m_links[i].m_rl;
		}
		serializer->finalizeChunk(chunk, "SoftBodyLinkData", BT_ARRAY_CODE, (void *) &m_links[0]);

	}


	sbd->m_numFaces = m_faces.size();
	sbd->m_faces = sbd->m_numFaces ? (SoftBodyFaceData *) serializer->getUniquePointer((void *)&m_faces[0]) : 0;
	if (sbd->m_faces)
	{
		int sz = sizeof(SoftBodyFaceData);
		int numElem = sbd->m_numFaces;
		btChunk *chunk = serializer->allocate(sz, numElem);
		SoftBodyFaceData *memPtr = (SoftBodyFaceData *)chunk->m_oldPtr;
		for (int i = 0; i < numElem; i++, memPtr++)
		{
			memPtr->m_material = m_faces[i].m_material ?  (SoftBodyMaterialData *) serializer->getUniquePointer((void *)m_faces[i].m_material) : 0;
			m_faces[i].m_normal.serializeFloat( memPtr->m_normal);
			for (int j = 0; j < 3; j++)
			{
				memPtr->m_nodeIndices[j] = m_faces[i].m_n[j] ? m_faces[i].m_n[j] - &m_nodes[0] : -1;
			}
			memPtr->m_restArea = m_faces[i].m_ra;
		}
		serializer->finalizeChunk(chunk, "SoftBodyFaceData", BT_ARRAY_CODE, (void *) &m_faces[0]);
	}


	sbd->m_numTetrahedra = m_tetras.size();
	sbd->m_tetrahedra = sbd->m_numTetrahedra ? (SoftBodyTetraData *) serializer->getUniquePointer((void *)&m_tetras[0]) : 0;
	if (sbd->m_tetrahedra)
	{
		int sz = sizeof(SoftBodyTetraData);
		int numElem = sbd->m_numTetrahedra;
		btChunk *chunk = serializer->allocate(sz, numElem);
		SoftBodyTetraData *memPtr = (SoftBodyTetraData *)chunk->m_oldPtr;
		for (int i = 0; i < numElem; i++, memPtr++)
		{
			for (int j = 0; j < 4; j++)
			{
				m_tetras[i].m_c0[j].serializeFloat( memPtr->m_c0[j] );
				memPtr->m_nodeIndices[j] = m_tetras[j].m_n[j] ? m_tetras[j].m_n[j] - &m_nodes[0] : -1;
			}
			memPtr->m_c1 = m_tetras[i].m_c1;
			memPtr->m_c2 = m_tetras[i].m_c2;
			memPtr->m_material = m_tetras[i].m_material ? (SoftBodyMaterialData *)serializer->getUniquePointer((void *) m_tetras[i].m_material) : 0;
			memPtr->m_restVolume = m_tetras[i].m_rv;
		}
		serializer->finalizeChunk(chunk, "SoftBodyTetraData", BT_ARRAY_CODE, (void *) &m_tetras[0]);
	}

	sbd->m_numAnchors = m_anchors.size();
	sbd->m_anchors = sbd->m_numAnchors ? (SoftRigidAnchorData *) serializer->getUniquePointer((void *)&m_anchors[0]) : 0;
	if (sbd->m_anchors)
	{
		int sz = sizeof(SoftRigidAnchorData);
		int numElem = sbd->m_numAnchors;
		btChunk *chunk = serializer->allocate(sz, numElem);
		SoftRigidAnchorData *memPtr = (SoftRigidAnchorData *)chunk->m_oldPtr;
		for (int i = 0; i < numElem; i++, memPtr++)
		{
			m_anchors[i].m_c0.serializeFloat(memPtr->m_c0);
			m_anchors[i].m_c1.serializeFloat(memPtr->m_c1);
			memPtr->m_c2 = m_anchors[i].m_c2;
			m_anchors[i].m_local.serializeFloat(memPtr->m_localFrame);
			memPtr->m_nodeIndex = m_anchors[i].m_node ? m_anchors[i].m_node - &m_nodes[0] : -1;

			memPtr->m_rigidBody = m_anchors[i].m_body ? (btRigidBodyData *)  serializer->getUniquePointer((void *)m_anchors[i].m_body) : 0;
			btAssert(memPtr->m_nodeIndex < m_nodes.size());
		}
		serializer->finalizeChunk(chunk, "SoftRigidAnchorData", BT_ARRAY_CODE, (void *) &m_anchors[0]);
	}


	sbd->m_config.m_dynamicFriction = m_cfg.kDF;
	sbd->m_config.m_baumgarte = m_cfg.kVCF;
	sbd->m_config.m_pressure = m_cfg.kPR;
	sbd->m_config.m_aeroModel = this->m_cfg.aeromodel;
	sbd->m_config.m_lift = m_cfg.kLF;
	sbd->m_config.m_drag = m_cfg.kDG;
	sbd->m_config.m_positionIterations = m_cfg.piterations;
	sbd->m_config.m_driftIterations = m_cfg.diterations;
	sbd->m_config.m_clusterIterations = m_cfg.citerations;
	sbd->m_config.m_velocityIterations = m_cfg.viterations;
	sbd->m_config.m_maxVolume = m_cfg.maxvolume;
	sbd->m_config.m_damping = m_cfg.kDP;
	sbd->m_config.m_poseMatch = m_cfg.kMT;
	sbd->m_config.m_collisionFlags = m_cfg.collisions;
	sbd->m_config.m_volume = m_cfg.kVC;
	sbd->m_config.m_rigidContactHardness = m_cfg.kCHR;
	sbd->m_config.m_kineticContactHardness = m_cfg.kKHR;
	sbd->m_config.m_softContactHardness = m_cfg.kSHR;
	sbd->m_config.m_anchorHardness = m_cfg.kAHR;
	sbd->m_config.m_timeScale = m_cfg.timescale;
	sbd->m_config.m_maxVolume = m_cfg.maxvolume;
	sbd->m_config.m_softRigidClusterHardness = m_cfg.kSRHR_CL;
	sbd->m_config.m_softKineticClusterHardness = m_cfg.kSKHR_CL;
	sbd->m_config.m_softSoftClusterHardness = m_cfg.kSSHR_CL;
	sbd->m_config.m_softRigidClusterImpulseSplit = m_cfg.kSR_SPLT_CL;
	sbd->m_config.m_softKineticClusterImpulseSplit = m_cfg.kSK_SPLT_CL;
	sbd->m_config.m_softSoftClusterImpulseSplit = m_cfg.kSS_SPLT_CL;

	//pose for shape matching
	{
		sbd->m_pose = (SoftBodyPoseData *)serializer->getUniquePointer((void *)&m_pose);

		int sz = sizeof(SoftBodyPoseData);
		btChunk *chunk = serializer->allocate(sz, 1);
		SoftBodyPoseData *memPtr = (SoftBodyPoseData *)chunk->m_oldPtr;

		m_pose.m_aqq.serializeFloat(memPtr->m_aqq);
		memPtr->m_bframe = m_pose.m_bframe;
		memPtr->m_bvolume = m_pose.m_bvolume;
		m_pose.m_com.serializeFloat(memPtr->m_com);

		memPtr->m_numPositions = m_pose.m_pos.size();
		memPtr->m_positions = memPtr->m_numPositions ? (btVector3FloatData *)serializer->getUniquePointer((void *)&m_pose.m_pos[0]) : 0;
		if (memPtr->m_numPositions)
		{
			int numElem = memPtr->m_numPositions;
			int sz = sizeof(btVector3Data);
			btChunk *chunk = serializer->allocate(sz, numElem);
			btVector3FloatData *memPtr = (btVector3FloatData *)chunk->m_oldPtr;
			for (int i = 0; i < numElem; i++, memPtr++)
			{
				m_pose.m_pos[i].serializeFloat(*memPtr);
			}
			serializer->finalizeChunk(chunk, "btVector3FloatData", BT_ARRAY_CODE, (void *)&m_pose.m_pos[0]);
		}
		memPtr->m_restVolume = m_pose.m_volume;
		m_pose.m_rot.serializeFloat(memPtr->m_rot);
		m_pose.m_scl.serializeFloat(memPtr->m_scale);

		memPtr->m_numWeigts = m_pose.m_wgh.size();
		memPtr->m_weights = memPtr->m_numWeigts ? (float *) serializer->getUniquePointer((void *) &m_pose.m_wgh[0]) : 0;
		if (memPtr->m_numWeigts)
		{

			int numElem = memPtr->m_numWeigts;
			int sz = sizeof(float);
			btChunk *chunk = serializer->allocate(sz, numElem);
			float *memPtr = (float *) chunk->m_oldPtr;
			for (int i = 0; i < numElem; i++, memPtr++)
			{
				*memPtr = m_pose.m_wgh[i];
			}
			serializer->finalizeChunk(chunk, "float", BT_ARRAY_CODE, (void *)&m_pose.m_wgh[0]);
		}

		serializer->finalizeChunk(chunk, "SoftBodyPoseData", BT_ARRAY_CODE, (void *)&m_pose);
	}

	//clusters for convex-cluster collision detection

	sbd->m_numClusters = m_clusters.size();
	sbd->m_clusters = sbd->m_numClusters ? (SoftBodyClusterData *) serializer->getUniquePointer((void *)m_clusters[0]) : 0;
	if (sbd->m_numClusters)
	{
		int numElem = sbd->m_numClusters;
		int sz = sizeof(SoftBodyClusterData);
		btChunk *chunk = serializer->allocate(sz, numElem);
		SoftBodyClusterData *memPtr = (SoftBodyClusterData *) chunk->m_oldPtr;
		for (int i = 0; i < numElem; i++, memPtr++)
		{
			memPtr->m_adamping = m_clusters[i]->m_adamping;
			m_clusters[i]->m_av.serializeFloat(memPtr->m_av);
			memPtr->m_clusterIndex = m_clusters[i]->m_clusterIndex;
			memPtr->m_collide = m_clusters[i]->m_collide;
			m_clusters[i]->m_com.serializeFloat(memPtr->m_com);
			memPtr->m_containsAnchor = m_clusters[i]->m_containsAnchor;
			m_clusters[i]->m_dimpulses[0].serializeFloat(memPtr->m_dimpulses[0]);
			m_clusters[i]->m_dimpulses[1].serializeFloat(memPtr->m_dimpulses[1]);
			m_clusters[i]->m_framexform.serializeFloat(memPtr->m_framexform);
			memPtr->m_idmass = m_clusters[i]->m_idmass;
			memPtr->m_imass = m_clusters[i]->m_imass;
			m_clusters[i]->m_invwi.serializeFloat(memPtr->m_invwi);
			memPtr->m_ldamping = m_clusters[i]->m_ldamping;
			m_clusters[i]->m_locii.serializeFloat(memPtr->m_locii);
			m_clusters[i]->m_lv.serializeFloat(memPtr->m_lv);
			memPtr->m_matching = m_clusters[i]->m_matching;
			memPtr->m_maxSelfCollisionImpulse = m_clusters[i]->m_maxSelfCollisionImpulse;
			memPtr->m_ndamping = m_clusters[i]->m_ndamping;
			memPtr->m_ldamping = m_clusters[i]->m_ldamping;
			memPtr->m_adamping = m_clusters[i]->m_adamping;
			memPtr->m_selfCollisionImpulseFactor = m_clusters[i]->m_selfCollisionImpulseFactor;

			memPtr->m_numFrameRefs = m_clusters[i]->m_framerefs.size();
			memPtr->m_numMasses = m_clusters[i]->m_masses.size();
			memPtr->m_numNodes = m_clusters[i]->m_nodes.size();

			memPtr->m_nvimpulses = m_clusters[i]->m_nvimpulses;
			m_clusters[i]->m_vimpulses[0].serializeFloat(memPtr->m_vimpulses[0]);
			m_clusters[i]->m_vimpulses[1].serializeFloat(memPtr->m_vimpulses[1]);
			memPtr->m_ndimpulses = m_clusters[i]->m_ndimpulses;



			memPtr->m_framerefs = memPtr->m_numFrameRefs ? (btVector3FloatData *)serializer->getUniquePointer((void *)&m_clusters[i]->m_framerefs[0]) : 0;
			if (memPtr->m_framerefs)
			{
				int numElem = memPtr->m_numFrameRefs;
				int sz = sizeof(btVector3FloatData);
				btChunk *chunk = serializer->allocate(sz, numElem);
				btVector3FloatData *memPtr = (btVector3FloatData *) chunk->m_oldPtr;
				for (int j = 0; j < numElem; j++, memPtr++)
				{
					m_clusters[i]->m_framerefs[j].serializeFloat(*memPtr);
				}
				serializer->finalizeChunk(chunk, "btVector3FloatData", BT_ARRAY_CODE, (void *)&m_clusters[i]->m_framerefs[0]);
			}

			memPtr->m_masses = memPtr->m_numMasses ? (float *) serializer->getUniquePointer((void *)&m_clusters[i]->m_masses[0]) : 0;
			if (memPtr->m_masses)
			{
				int numElem = memPtr->m_numMasses;
				int sz = sizeof(float);
				btChunk *chunk = serializer->allocate(sz, numElem);
				float *memPtr = (float *) chunk->m_oldPtr;
				for (int j = 0; j < numElem; j++, memPtr++)
				{
					*memPtr = m_clusters[i]->m_masses[j];
				}
				serializer->finalizeChunk(chunk, "float", BT_ARRAY_CODE, (void *)&m_clusters[i]->m_masses[0]);
			}

			memPtr->m_nodeIndices  = memPtr->m_numNodes ? (int *) serializer->getUniquePointer((void *) &m_clusters[i]->m_nodes) : 0;
			if (memPtr->m_nodeIndices )
			{
				int numElem = memPtr->m_numMasses;
				int sz = sizeof(int);
				btChunk *chunk = serializer->allocate(sz, numElem);
				int *memPtr = (int *) chunk->m_oldPtr;
				for (int j = 0; j < numElem; j++, memPtr++)
				{
					int *indexPtr = m_nodeIndexMap.find(m_clusters[i]->m_nodes[j]);
					btAssert(indexPtr);
					*memPtr = *indexPtr;
				}
				serializer->finalizeChunk(chunk, "int", BT_ARRAY_CODE, (void *)&m_clusters[i]->m_nodes);
			}
		}
		serializer->finalizeChunk(chunk, "SoftBodyClusterData", BT_ARRAY_CODE, (void *)m_clusters[0]);

	}



	sbd->m_numJoints = m_joints.size();
	sbd->m_joints = m_joints.size() ? (btSoftBodyJointData *) serializer->getUniquePointer((void *)&m_joints[0]) : 0;

	if (sbd->m_joints)
	{
		int sz = sizeof(btSoftBodyJointData);
		int numElem = m_joints.size();
		btChunk *chunk = serializer->allocate(sz, numElem);
		btSoftBodyJointData *memPtr = (btSoftBodyJointData *)chunk->m_oldPtr;

		for (int i = 0; i < numElem; i++, memPtr++)
		{
			memPtr->m_jointType = (int)m_joints[i]->Type();
			m_joints[i]->m_refs[0].serializeFloat(memPtr->m_refs[0]);
			m_joints[i]->m_refs[1].serializeFloat(memPtr->m_refs[1]);
			memPtr->m_cfm = m_joints[i]->m_cfm;
			memPtr->m_erp = m_joints[i]->m_erp;
			memPtr->m_split = m_joints[i]->m_split;
			memPtr->m_delete = m_joints[i]->m_delete;

			for (int j = 0; j < 4; j++)
			{
				memPtr->m_relPosition[0].m_floats[j] = 0.f;
				memPtr->m_relPosition[1].m_floats[j] = 0.f;
			}
			memPtr->m_bodyA = 0;
			memPtr->m_bodyB = 0;
			if (m_joints[i]->m_bodies[0].m_soft)
			{
				memPtr->m_bodyAtype = BT_JOINT_SOFT_BODY_CLUSTER;
				memPtr->m_bodyA = serializer->getUniquePointer((void *)m_joints[i]->m_bodies[0].m_soft);
			}
			if (m_joints[i]->m_bodies[0].m_collisionObject)
			{
				memPtr->m_bodyAtype = BT_JOINT_COLLISION_OBJECT;
				memPtr->m_bodyA = serializer->getUniquePointer((void *)m_joints[i]->m_bodies[0].m_collisionObject);
			}
			if (m_joints[i]->m_bodies[0].m_rigid)
			{
				memPtr->m_bodyAtype = BT_JOINT_RIGID_BODY;
				memPtr->m_bodyA = serializer->getUniquePointer((void *)m_joints[i]->m_bodies[0].m_rigid);
			}

			if (m_joints[i]->m_bodies[1].m_soft)
			{
				memPtr->m_bodyBtype = BT_JOINT_SOFT_BODY_CLUSTER;
				memPtr->m_bodyB = serializer->getUniquePointer((void *)m_joints[i]->m_bodies[1].m_soft);
			}
			if (m_joints[i]->m_bodies[1].m_collisionObject)
			{
				memPtr->m_bodyBtype = BT_JOINT_COLLISION_OBJECT;
				memPtr->m_bodyB = serializer->getUniquePointer((void *)m_joints[i]->m_bodies[1].m_collisionObject);
			}
			if (m_joints[i]->m_bodies[1].m_rigid)
			{
				memPtr->m_bodyBtype = BT_JOINT_RIGID_BODY;
				memPtr->m_bodyB = serializer->getUniquePointer((void *)m_joints[i]->m_bodies[1].m_rigid);
			}
		}
		serializer->finalizeChunk(chunk, "btSoftBodyJointData", BT_ARRAY_CODE, (void *) &m_joints[0]);
	}


	return btSoftBodyDataName;
}

