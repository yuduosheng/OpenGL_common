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

extern ofstream debug;

float scale = 5.0f;

const int width = 800, height = 600;

float timeStep = 1 / 60.0f;
float currentTime = 0;
float accumulator = timeStep;

int isMouseButtonDown = 0;


int oldX = 0, oldY = 0;
float rX = 15, rY = 0;
int state = 1;
float dist = -2.5f;
const int GRID_SIZE = 1;


glm::vec3 gravity = glm::vec3(0.0f, 0.0f, 0.0f);


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
float l = 0;
int n = 199;
void firstPass(HashMap H, vector<DeformableObject> objects)
{
	int objId = 0;
	vector<DeformableObject>::iterator object_iter;
	for (object_iter = objects.begin(); object_iter != objects.end(); object_iter++)
	{
		for (int i = 0; i < object_iter->total_points; ++i)
		{
			glm::vec3 p = object_iter->X[i];
			int x = (int)(p.x / l);
			int y = (int)(p.y / l);
			int z = (int)(p.z / l);
			int h = ((x * p1) ^ (y * p2) ^ (z * p3)) % n;
			if (H.cell[h].T != H.T)
			{
				H.cell[h].nodes.clear();
				H.cell[h].T = H.T;
			}
			H.cell[h].nodes.push_back(node(objId, i, object_iter->X[i]));
		}
		++objId;
	}
}
void secondPass(HashMap H, vector<DeformableObject> objects)
{
	int objId = 0;
	vector<DeformableObject>::iterator object_iter;
	for (object_iter = objects.begin(); object_iter != objects.end(); object_iter++)
	{
		for (int i = 0; i < object_iter->total_tetrahedra; ++i)
		{
			float Mx, My, Mz;//max
			float mx, my, mz;//min
			int n1 = object_iter->tetrahedra[i].indices[0];
			int n2 = object_iter->tetrahedra[i].indices[1];
			int n3 = object_iter->tetrahedra[i].indices[2];
			int n4 = object_iter->tetrahedra[i].indices[3];
			glm::vec3 v0 = object_iter->X[n1];
			glm::vec3 v1 = object_iter->X[n2];
			glm::vec3 v2 = object_iter->X[n3];
			glm::vec3 v3 = object_iter->X[n4];

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

			int m_x = int(mx / l);
			int m_y = int(my / l);
			int m_z = int(mz / l);
			int M_x = int(Mx / l) + 1;
			int M_y = int(My / l) + 1;
			int M_z = int(Mz / l) + 1;

			for (int j = m_x; j < M_x; ++j)
				for (int k = m_y; k < M_y; ++k)
					for (int l = m_z; l < M_z; ++l)
					{
				         int h = ((j * p1) ^ (k * p2) ^ (l * p3)) % n;
						 if (H.cell[h].T == H.T)
						 {
							 list<node> ::iterator node_iter;
							 for (node_iter = H.cell[h].nodes.begin(); node_iter != H.cell[h].nodes.end(); node_iter++)
							 {
								 if (node_iter->objID == objId && ((node_iter->index == n1) || (node_iter->index == n2) || (node_iter->index == n3) || (node_iter->index == n4)))
									 break;
								 glm::mat3 A;
								 glm::vec3 e1, e2, e3;
								 e1 = v1 - v0;
								 e2 = v2 - v0;
								 e3 = v3 - v0;
								 A[0][0] = e1.x;    A[1][0] = e2.x;    A[2][0] = e3.x;
								 A[0][1] = e1.y;    A[1][1] = e2.y;    A[2][1] = e3.y;
								 A[0][2] = e1.z;    A[1][2] = e2.z;    A[2][2] = e3.z;

								 A = glm::inverse(A);

								 glm::vec3 beta = node_iter->position - v0;
								 beta = A * beta;

								 if (beta.x >= 0 && beta.y >= 0 && beta.z >= 0 && (beta.x + beta.y + beta.z) <= 1)
									 ;//contact(n,t)
							 }
						 }
					}
		}
	}
}

void spatialHashing(double T, HashMap H, vector<DeformableObject> objects)
{
	H.T = T;
	firstPass(H, objects);
	secondPass(H, objects);
}

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
private:

	GameTimer                timer;
	vector<DeformableObject> models;
	DeformableObject         bunny;
	DeformableObject         bunny2;
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

FEMTest::FEMTest() : App(), timer(), bunny(), bunny2()
{
	models.push_back(bunny);
	models.push_back(bunny2);
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
		bunny.StepPhysics(timeStep);
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
	if ((newTime - startTime)>1000)
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
	bunny.renderModel();

}
void FEMTest::onMouseWheel(GLFWwindow* window, double x, double y)
{
	float mouseWheelScale = 0.3f;
	scale += mouseWheelScale  * (float)y;
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
		bunny.bUseStiffnessWarping = !bunny.bUseStiffnessWarping;
		printf("Stiffness Warping %s\n", bunny.bUseStiffnessWarping ? "On" : "Off");
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
	for (object_iter = models.begin(); object_iter != models.end(); object_iter++)
	{
		length += object_iter->totalLength;
		tn += object_iter->total_tetrahedra;
	}
	l = length / (6 * tn);
}