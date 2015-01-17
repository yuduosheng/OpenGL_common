#include <png.h>
#include "App.h"
#include "camera.h"
#include <vector>
#include <map>
#include "GameTimer.h"
#include <GL/wglew.h>
#include <GL/freeglut.h>
#include "deformableobject.h"
#include "hash.h"
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
using namespace Eigen;

#pragma comment(lib, "glew32.lib")
#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glfw3.lib")
#pragma comment(lib, "OpenMeshCore.lib")
#pragma comment(lib, "OpenMeshTools.lib")
//#pragma pack(push,1)


extern ofstream debug;
extern ofstream debug2;
float scale = 3.0f;

const int width = 1024, height = 768;

float timeStep = 1 / 60.0f;
float currentTime = 0;
float accumulator = timeStep;

int isMouseButtonDown = 0;


int oldX = 0, oldY = 0;
float rX = 15, rY = 0;
int state = 1;
float dist = -5.0f;
const int GRID_SIZE = 10;


GLint viewport[4];
GLdouble MV[16];
GLdouble P[16];

glm::vec3 Up = glm::vec3(0, 1, 0), Right, viewDir;

LARGE_INTEGER frequency;        // ticks per second
LARGE_INTEGER t1, t2;           // ticks
float frameTimeQP = 0;
float frameTime = 0;

float startTime = 0, fps = 0;
int totalFrames = 0;

char info[MAX_PATH] = { 0 };


void DrawGrid()
{
	glBegin(GL_LINES);
	glColor3f(0.5f, 0.5f, 0.5f);
	for (int i = -GRID_SIZE; i <= GRID_SIZE; i++)
	{
		glVertex3f((float)i, 0, (float)-GRID_SIZE);
		glVertex3f((float)i, 0, (float)GRID_SIZE);

		glVertex3f((float)-GRID_SIZE, 0, (float)i);
		glVertex3f((float)GRID_SIZE, 0, (float)i);
	}
	glEnd();
}
/*
void EigenSolve()
{
	SparseMatrix<double> eA;
	VectorXd             eb, ev;
	eA.resize(3 * total_points, 3 * total_points);
	eA.setZero();
	eb.resize(3 * total_points);
	eb.setZero();
	ev.resize(3 * total_points);

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;

	for (size_t k = 0; k < total_points; k++)
	{
		if (IsFixed[k])
			continue;
		matrix_map tmp = A_row[k];
		matrix_iterator Abegin = tmp.begin();
		matrix_iterator Aend = tmp.end();

		for (matrix_iterator K = Abegin; K != Aend; ++K)
		{
			unsigned int a = K->first;
			glm::mat3 A_ij = K->second;
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					tripletList.push_back(T(k * 3 + i, a * 3 + j, A_ij[i][j]));
				}
			}

		}

	}
	eA.setFromTriplets(tripletList.begin(), tripletList.end());

	for (size_t k = 0; k < total_points; ++k)
	{
		int ind = k * 3;
		eb[ind] = b[k].x;
		eb[ind + 1] = b[k].y;
		eb[ind + 2] = b[k].z;
	}


	ConjugateGradient<SparseMatrix<double> > solver;
	solver.compute(eA);
	if (solver.info() != Success) {
		// decomposition failed
		return;
	}
	ev = solver.solve(eb);
	if (solver.info() != Success) {
		// solving failed
		return;
	}

	for (size_t k = 0; k < total_points; ++k)
	{
		int ind = k * 3;
		V[k].x = ev[ind];
		V[k].y = ev[ind + 1];
		V[k].z = ev[ind + 2];
	}
}*/


__int64 p1 = 73856093;
__int64 p2 = 19349663;
__int64 p3 = 83492791;
float g_l = 0;//global grid size
int n = 199;//hash table size
float a = 9.81;//penalty coefficient
HashMap g_HashTable(n);
vector<DeformableObject> g_models;

class FEMTest : public App
{
public:
	FEMTest();
	~FEMTest();

	bool                    Init();
	void                    UpdateScene();
	void                    Rendering();
	void                    onResize(GLFWwindow* window, int w, int h);

	void                    onMouseWheel(GLFWwindow* window, double x, double y);
	void                    onMouseMove(GLFWwindow* window, double xd, double yd);
	void                    onMouseButton(GLFWwindow* window, int button, int action, int mods);
	void                    onKey(GLFWwindow* window, int key, int scancode, int action, int mods);

private:
	void                    buildGeometryBuffers();
	void                    buildShader();
	void                    setGridCellSize();
	void firstPass(HashMap H, vector<DeformableObject> *g_models);
	void secondPass(HashMap H, vector<DeformableObject> *g_models);
private:

	GameTimer                timer;
	DeformableObject         bunny;
	DeformableObject         bunny2;
	//DeformableObject         block;
	HashMap                  HashTable;
};

int main(void)
{
	FEMTest *theApp = new FEMTest;

	if (!theApp->Init())
		return 0;
	theApp->Run();
	delete theApp;
	return 0;
}

FEMTest::FEMTest() : App(), timer(), HashTable(n), bunny(0.0f, 0.30f, 0.0f, false), bunny2(0.0f, 0.0f, 0.0f, true)//, block(10, 3, 3, 0.1f, 0.1f, 0.1f),
{
	//g_models.push_back(block);
	//g_models.push_back(bunny);
	g_models.push_back(bunny2);
	mWidth = width;
	mHeight = height;
}

FEMTest::~FEMTest()
{
	// Cleanup VBO and shader
}

bool FEMTest::Init()
{
	if (!App::Init())
		return false;

	//startTime = (float)glutGet(GLUT_ELAPSED_TIME);
	//currentTime = startTime;

	// get ticks per second
	//QueryPerformanceFrequency(&frequency);

	// start timer
	//QueryPerformanceCounter(&t1);

	timer.Reset();
	startTime = timer.getCurrenTime();
	currentTime = startTime;

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glPointSize(2); 
	wglSwapIntervalEXT(0);

	onResize(window, width, height);

	setGridCellSize();

	return true;
}

void FEMTest::onResize(GLFWwindow* window, int nw, int nh)
{
	App::onResize(window, nw, nh);

	glViewport(0, 0, nw, nh);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, (GLfloat)nw / (GLfloat)nh, 0.1f, 100.0f);

	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_PROJECTION_MATRIX, P);

	glMatrixMode(GL_MODELVIEW);
}

void FEMTest::UpdateScene()
{
	if (accumulator >= timeStep)
	{
		HashTable.T = timer.TotalTime();
		for (int i = 0; i < g_models.size(); ++i)
		{
			g_models[i].UpdateGradient();
		}
		for (int i = 0; i < g_models.size(); ++i)
		{
			//g_models[i].firstPass(&HashTable, i);
		}
		
		for (int i = 0; i < g_models.size(); ++i)
		{
			//g_models[i].secondPass(&HashTable, i);
		}

		for (int i = 0; i < g_models.size(); ++i)
		{
             g_models[i].StepPhysics(1/1200.0f);
		}
			
		accumulator -= timeStep;
	}
}
void FEMTest::Rendering()
{
	timer.Tick();
	float newTime = timer.getCurrenTime();
	frameTime = newTime - currentTime;
	currentTime = newTime;
	//accumulator += frameTime;

	//Using high res. counter
	//QueryPerformanceCounter(&t2);
	// compute and print the elapsed time in millisec
	//frameTimeQP = (t2.QuadPart - t1.QuadPart) * 1000.0f / frequency.QuadPart;
	frameTimeQP = timer.DeltaTime();
	//t1 = t2;
	accumulator += frameTimeQP;

	++totalFrames;
	if ((newTime - startTime) > 1000)
	{
		float elapsedTime = (newTime - startTime);
		fps = (totalFrames / elapsedTime) * 1000;
		startTime = newTime;
		totalFrames = 0;
	}

	sprintf_s(info, "FPS: %3.2f, Frame time (GLUT): %3.4f msecs, Frame time (QP): %3.3f, Stiffness Warp: %s", fps, frameTime, frameTimeQP, bunny.bUseStiffnessWarping ? "On" : "Off");
	glfwSetWindowTitle(window, info);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	glTranslatef(0, 0, dist);
	glRotatef(rX, 1, 0, 0);
	glRotatef(rY, 0, 1, 0);

	glMatrixMode(GL_MODELVIEW);
	glScaled(scale, scale, scale);


	glGetDoublev(GL_MODELVIEW_MATRIX, MV);
	viewDir.x = (float)-MV[2];
	viewDir.y = (float)-MV[6];
	viewDir.z = (float)-MV[10];
	Right = glm::cross(viewDir, Up);

	//draw grid
	DrawGrid();
	glMatrixMode(GL_MODELVIEW);
	glScaled(2.0f, 2.0f, 2.0f);

	for (int i = 0; i < g_models.size(); ++i)
	     g_models[i].renderModel();
	/*
	//draw grid
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	for (int i = 0; i < g_models[0].gi - 1; ++i)
	//for (int i = 2; i < 3; ++i)
		//for (int j = 7; j < 8; ++j)
		for (int j = 0; j < g_models[0].gj - 1; ++j)
		//for (int k = 1; k < 2; ++k)
			for (int k = 0; k < g_models[0].gk - 1; ++k)
			{
		int n1 = i * g_models[0].gj * g_models[0].gk + j * g_models[0].gk + k;
		int n2 = n1 + 1;
		int n3 = i * g_models[0].gj * g_models[0].gk + (j + 1) * g_models[0].gk + k;
		int n4 = n3 + 1;
		int n5 = (i + 1) * g_models[0].gj * g_models[0].gk + j * g_models[0].gk + k;
		int n6 = n5 + 1;
		int n7 = (i + 1) * g_models[0].gj * g_models[0].gk + (j + 1) * g_models[0].gk + k;
		int n8 = n7 + 1;



		glm::vec3 p1 = g_models[0].grid[n1];
		glm::vec3 p2 = g_models[0].grid[n2];
		glm::vec3 p3 = g_models[0].grid[n3];
		glm::vec3 p4 = g_models[0].grid[n4];
		glm::vec3 p5 = g_models[0].grid[n5];
		glm::vec3 p6 = g_models[0].grid[n6];
		glm::vec3 p7 = g_models[0].grid[n7];
		glm::vec3 p8 = g_models[0].grid[n8];

		glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p2.x, p2.y, p2.z);
		glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p2.x, p2.y, p2.z);
		glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p3.x, p3.y, p3.z);
		glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p3.x, p3.y, p3.z);

		glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p5.x, p5.y, p5.z);
		glVertex3f(p2.x, p2.y, p2.z);		glVertex3f(p6.x, p6.y, p6.z);
		glVertex3f(p7.x, p7.y, p7.z);		glVertex3f(p3.x, p3.y, p3.z);
		glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p8.x, p8.y, p8.z);

		glVertex3f(p5.x, p5.y, p5.z);		glVertex3f(p6.x, p6.y, p6.z);
		glVertex3f(p5.x, p5.y, p5.z);		glVertex3f(p7.x, p7.y, p7.z);
		glVertex3f(p6.x, p6.y, p6.z);		glVertex3f(p8.x, p8.y, p8.z);
		glVertex3f(p7.x, p7.y, p7.z);		glVertex3f(p8.x, p8.y, p8.z);
			}
	glEnd();
	
	
	{
	//draw gradient of grid
	glColor3f(0.0, 0.5, 1.0);
	glBegin(GL_LINES);
	for (int i = 1; i < g_models[0].gi - 1; ++i)
	//for (int i = 2; i < 4; ++i)
		//for (int j = 7; j < 9; ++j)
	for (int j = 1; j < g_models[0].gj - 1; ++j)
	//for (int k = 1; k < 3; ++k)
		for (int k = 1; k < g_models[0].gk - 1; ++k)
		{
		int index = i * g_models[0].gj * g_models[0].gk + j * g_models[0].gk + k;
		glm::vec3 p1, p2;
		p1 = g_models[0].grid[index];
		p2 = g_models[0].gradientFild[index];
		//p2 = glm::normalize(p2);
		p2 = p1 + p2;
		glVertex3f(p1.x, p1.y, p1.z);
		glVertex3f(p2.x, p2.y, p2.z);
		}
	glEnd();
}
	*/
	
	//draw gradient of vertex
	glColor3f(1.0, 1.0, 0.0);
	glBegin(GL_LINES);
	for (int i = 30; i < 31; ++i)
	//for (int i = 0; i < g_models[0].total_points; ++i)
	{
		glm::vec3 p1, p2;
		p1 = g_models[0].X[i];
		p2 = g_models[0].gradientV[i];
		//p2 = glm::normalize(p2);
		p2 = p1 + p2;	
		glVertex3f(p1.x, p1.y, p1.z);
		glVertex3f(p2.x, p2.y, p2.z);
	}
	glEnd();
/*
	
	//is point in tetrahedron

		glColor3f(0.75, 0.75, 0.75);
		{
		glBegin(GL_LINES);

		glm::vec3 p1 = glm::vec3(-0.0158555,   0.149626, - 0.0252782);
		glm::vec3 p2 = glm::vec3(-0.0208207,   0.141843, - 0.0189655);
		glm::vec3 p3 = glm::vec3(-0.0289981,   0.135866, - 0.0183379);
		glm::vec3 p4 = glm::vec3(-0.0131466,   0.140465, - 0.0290179);

			glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p1.x, p1.y, p1.z);
			glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p2.x, p2.y, p2.z);
			glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p3.x, p3.y, p3.z);

			glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p2.x, p2.y, p2.z);
			glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p3.x, p3.y, p3.z);

			glVertex3f(p2.x, p2.y, p2.z);		glVertex3f(p3.x, p3.y, p3.z);
		glEnd();
}
*/
/*
		//draw points	
		glBegin(GL_POINTS);
		    glm::vec3 p = glm::vec3(g_models[0].ABmx, g_models[0].ABmy, g_models[0].ABmz);
			glColor3f(1.0f, 1.0f, 0.0f);
			glVertex3f(p.x, p.y, p.z);

			glm::vec3 p2 = glm::vec3(g_models[0].ABMx, g_models[0].ABMy, g_models[0].ABMz);
			glColor3f(1.0f, 1.0f, 0.0f);
			glVertex3f(p2.x, p2.y, p2.z);
		glEnd();*/

}
void FEMTest::onMouseWheel(GLFWwindow* window, double x, double y)
{
	float mouseWheelScale = 0.3f;
	scale += mouseWheelScale  * (float)y;
	if (scale < 0) scale = 0;
}
void FEMTest::onMouseMove(GLFWwindow* window, double xd, double yd)
{
	double x = xd;
	double y = yd;
	if (isMouseButtonDown)
	{
		if (bunny.GetSelectIndex() == -1) {
			if (state == 0)
				dist *= (1 + (y - oldY) / 60.0f);
			else
			{
				rY += (x - oldX) / 5.0f;
				rX += (y - oldY) / 5.0f;
			}
		}
		else {

			float delta = 1000 / abs(dist);
			float valX = (x - oldX) / delta;
			float valY = (oldY - y) / delta;

			bunny.V[bunny.GetSelectIndex()] = glm::vec3(0);
			bunny.X[bunny.GetSelectIndex()].x += Right[0] * valX;
			float newValue = bunny.X[bunny.GetSelectIndex()].y + Up[1] * valY;
			if (newValue>0)
				bunny.X[bunny.GetSelectIndex()].y = newValue;
			bunny.X[bunny.GetSelectIndex()].z += Right[2] * valX + Up[2] * valY;
		}
		oldX = x;
		oldY = y;
	}

}
void FEMTest::onMouseButton(GLFWwindow* window, int button, int action, int mods)
{
	double xd, yd;

	if ((button == GLFW_MOUSE_BUTTON_1) && (action == GLFW_PRESS))
	{
		glfwGetCursorPos(window, &xd, &yd);
		oldX = xd;
		oldY = yd;
		int window_y = (height - yd);
		float norm_y = float(window_y) / float(height / 2.0);
		int window_x = xd;
		float norm_x = float(window_x) / float(width / 2.0);

		float winZ = 0;
		glReadPixels(xd, height - yd, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
		double objX = 0, objY = 0, objZ = 0;

		gluUnProject(window_x, window_y, winZ, MV, P, viewport, &objX, &objY, &objZ);

		glm::vec3 pt(objX, objY, objZ);
		printf("\nObj [ %3.3f,%3.3f,%3.3f ]", objX, objY, objZ);
		size_t i = 0;
		for (i = 0; i < bunny.total_points; i++) {
			if (glm::distance(bunny.X[i], pt)<0.01) {
				bunny.SetSelectIndex(i);

				printf("Intersected at %d\n", i);
				printf("Pt [ %3.3f,%3.3f,%3.3f ]\n", bunny.X[i].x, bunny.X[i].y, bunny.X[i].z);
				break;
			}
		}
		isMouseButtonDown = 1;
	}
	else if ((button == GLFW_MOUSE_BUTTON_1) && (action == GLFW_RELEASE))
	{
		bunny.SetSelectIndex(-1);
		bunny.UpdateOrientation();
		isMouseButtonDown = 0;
	}

	if ((button == GLFW_MOUSE_BUTTON_2) && (action == GLFW_PRESS))
	{
		state = 0;
	}
	else if ((button == GLFW_MOUSE_BUTTON_2) && (action == GLFW_RELEASE))
	{
		state = 1;
	}

}
void FEMTest::onKey(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	App::onKey(window, key, scancode, action, mods);
	if ((key == GLFW_KEY_SPACE) && (action == GLFW_PRESS))
	{
		//bunny.bUseStiffnessWarping = !bunny.bUseStiffnessWarping;
		//printf("Stiffness Warping %s\n", bunny.bUseStiffnessWarping ? "On" : "Off");
		for (int i = 0; i < g_models.size(); ++i)
			g_models[i].Reset();
	}

}
void FEMTest::buildGeometryBuffers()
{

}
void FEMTest::buildShader()
{

}
void FEMTest::setGridCellSize()
{
	double length = 0;
	double tn = 0;
	vector<DeformableObject>::iterator object_iter;
	for (object_iter = g_models.begin(); object_iter != g_models.end(); object_iter++)
	{
		length += object_iter->totalLength;
		tn += object_iter->total_tetrahedra;
	}
	g_l = length / (6 * tn);
}

void FEMTest::firstPass(HashMap H, vector<DeformableObject> *g_models)
{

	for (int i = 0; i < g_models->size(); ++i)
	{
		DeformableObject *ob = &((*g_models)[i]);
		for (int j = 0; j < (*g_models)[i].total_points; ++j)
		{
			glm::vec3 p = ob->X[j];
			int x = (int)(p.x / g_l);
			int y = (int)(p.y / g_l);
			int z = (int)(p.z / g_l);
			int h = ((x * p1) ^ (y * p2) ^ (z * p3)) % n;
			if (h < 0) h = -h;
			//cout << l << endl;
			//debug << p.x << " " << p.y << " " << p.z << endl;
			//cout << x << " "<< y <<  " "<< z << endl;
			//debug << h << endl;
			if (H.cell[h].T != H.T)
			{
				H.cell[h].nodes.clear();
				H.cell[h].T = H.T;
			}
			Vertex v;
			v.objId = i;
			v.localIndex = j;
			v.p = p;
			H.cell[h].nodes.push_back(v);
		}
	}
}
void FEMTest::secondPass(HashMap H, vector<DeformableObject> *g_models)
{
	for (int mi = 0; mi < g_models->size(); ++mi)
	{
		DeformableObject *ob = &((*g_models)[mi]);
		for (int i = 0; i != ob->total_tetrahedra; ++i)
		{
			float Mx, My, Mz;//max
			float mx, my, mz;//min
			int n1 = ob->tetrahedra[i].indices[0];
			int n2 = ob->tetrahedra[i].indices[1];
			int n3 = ob->tetrahedra[i].indices[2];
			int n4 = ob->tetrahedra[i].indices[3];
			glm::vec3 v0 = ob->X[n1];
			glm::vec3 v1 = ob->X[n2];
			glm::vec3 v2 = ob->X[n3];
			glm::vec3 v3 = ob->X[n4];

			mx = glm::min(v0.x, v1.x);
			mx = glm::min(mx, v2.x);
			mx = glm::min(mx, v3.x);
			my = glm::min(v0.y, v1.y);
			my = glm::min(my, v2.y);
			my = glm::min(my, v3.y);
			mz = glm::min(v0.z, v1.z);
			mz = glm::min(mz, v2.z);
			mz = glm::min(mz, v3.z);
			Mx = glm::max(v0.x, v1.x);
			Mx = glm::max(Mx, v2.x);
			Mx = glm::max(Mx, v3.x);
			My = glm::max(v0.y, v1.y);
			My = glm::max(My, v2.y);
			My = glm::max(My, v3.y);
			Mz = glm::max(v0.z, v1.z);
			Mz = glm::max(Mz, v2.z);
			Mz = glm::max(Mz, v3.z);

			int m_x = int(mx / g_l);
			int m_y = int(my / g_l);
			int m_z = int(mz / g_l);
			int M_x = int(Mx / g_l) + 1;
			int M_y = int(My / g_l) + 1;
			int M_z = int(Mz / g_l) + 1;

			for (int j = m_x; j < M_x; ++j)
				for (int k = m_y; k < M_y; ++k)
					for (int l = m_z; l < M_z; ++l)
					{
				int h = ((j * p1) ^ (k * p2) ^ (l * p3)) % n;
				if (h < 0) h = -h;
				if (H.cell[h].T == H.T)
				{
					list<Vertex> ::iterator node_iter;
					for (node_iter = H.cell[h].nodes.begin(); node_iter != H.cell[h].nodes.end(); node_iter++)
					{
						if (node_iter->objId == mi && ((node_iter->localIndex == n1) || (node_iter->localIndex == n2) || (node_iter->localIndex == n3) || (node_iter->localIndex == n4)))
							continue;
						glm::mat3 A;
						glm::vec3 e1, e2, e3;
						e1 = v1 - v0;
						e2 = v2 - v0;
						e3 = v3 - v0;

						A[0][0] = e1.x;    A[1][0] = e2.x;    A[2][0] = e3.x;//A = glm::inverse(A);
						A[0][1] = e1.y;    A[1][1] = e2.y;    A[2][1] = e3.y;
						A[0][2] = e1.z;    A[1][2] = e2.z;    A[2][2] = e3.z;

						glm::vec3 beta = node_iter->p - v0;
						beta = A * beta;

						if (beta.x >= 0 && beta.y >= 0 && beta.z >= 0 && (beta.x + beta.y + beta.z) <= 1)
						{//contact(n,t)
							//ob->distanceV[n1];
							//ob->distanceV[n2];
							//ob->distanceV[n3];
							//ob->distanceV[n4];
							/*
							glm::vec3 g0(0), g1(0), g2(0), g3(0);

							float d = ob->distanceV[n1] * (1 - beta.x - beta.y - beta.z) + ob->distanceV[n2] * beta.x + ob->distanceV[n3] * beta.y + ob->distanceV[n4] * beta.z;

							for (int numT = 0; numT < ob->tetrahedraOfVertices[n1].size(); ++numT)
							{
								int tetrahedraIndex = ob->tetrahedraOfVertices[n1][numT];
								g0 += ob->tetrahedra[tetrahedraIndex].Re * ob->gradientV[n1];
							}
							for (int numT = 0; numT < ob->tetrahedraOfVertices[n2].size(); ++numT)
							{
								int tetrahedraIndex = ob->tetrahedraOfVertices[n2][numT];
								g1 += ob->tetrahedra[tetrahedraIndex].Re * ob->gradientV[n2];
							}
							for (int numT = 0; numT < ob->tetrahedraOfVertices[n3].size(); ++numT)
							{
								int tetrahedraIndex = ob->tetrahedraOfVertices[n3][numT];
								g2 += ob->tetrahedra[tetrahedraIndex].Re * ob->gradientV[n3];
							}
							for (int numT = 0; numT < ob->tetrahedraOfVertices[n4].size(); ++numT)
							{
								int tetrahedraIndex = ob->tetrahedraOfVertices[n4][numT];
								g3 += ob->tetrahedra[tetrahedraIndex].Re * ob->gradientV[n4];
							}
							g0 = glm::normalize(g0);
							g1 = glm::normalize(g1);
							g2 = glm::normalize(g2);
							g3 = glm::normalize(g3);
							glm::vec3 g = g0 * (1 - beta.x - beta.y - beta.z) + g1 * beta.x + g2 * beta.y + g3 * beta.z;
							glm::vec3 f = a * d * g;
							ob->F[node_iter->localIndex] += f;
							*/
						}

					}
				}
					}
		}
	}
}
