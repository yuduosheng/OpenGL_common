/// OPENMESH VIEWER source code (using vertex arrays)
/// by Alvaro Cuno
/// bugs? email-me: alvaroecp@gmail.com
/// October, 2008

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <GL/glut.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <queue>
#include <set>
#include <iterator>

#include "selection.h"
#include "arcball.hh"
#include "atmatrix.h"
#include "stroke.h"

#include "GeomTypes.h"
#include "intersections.h"
#include "Cloth.h"

typedef OpenMesh::PolyMesh_ArrayKernelT<>  PolyMesh;

using namespace std;

Cloth *mycloth = 0;

int winw = 700, winh = 700; // window size
int xini, yini; // initial mouse position
double xwini, ywini, zwini; // initial mouse position in world coordinates
int buttonpressed; // 
enum MeshRenderMode{
	POINTS = 0x01, WIREFRAME = 0x02, HIDDENLINES = 0x04, FLATLINES = 0x08,
	FLAT = 0x10, SMOOTH = 0x20, TRANSPARENCY = 0x40
};
MeshRenderMode rendermode; // mesh render mode

AMatrix<GLfloat> sceneT;
AMatrix<GLfloat> sceneIniT;
ArcBall arcball;

Stroke stroke;
Point3 grabStart, grabLast;

/// time info
float      cfps;
time_t     start, end;
unsigned   frame;
/// Info variables
void *fontsmall = (void *)GLUT_BITMAP_HELVETICA_12;
string nvertices("Vertices: ");
string nfaces("Faces: ");
string fps("FPS: ");
string options("Options: p/w/h/l/f/s/+/-");

/// Renders info about the model
void renderInfo();
/// Returns the world coordinates of a point in screen space
void screenToWorld(int x, int y, double &xw, double &yw, double &zw);
/// Implements the "zoomin" operation
void zoomIn();
/// Implements the "zoomout" operation
void zoomOut();

/// OpenGL display function
void display(void) {

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	glMultMatrixf(&sceneT);

	//   glEnable(GL_LIGHTING);
	//   glShadeModel(GL_SMOOTH);
	//   glEnableClientState(GL_VERTEX_ARRAY);
	//   glEnableClientState(GL_NORMAL_ARRAY);
	//   glVertexPointer(3, GL_FLOAT, 0, meshmodel->points());
	//   glNormalPointer(GL_FLOAT, 0, meshmodel->vertex_normals());
	//   glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, &indices[0]);
	//   glDisableClientState(GL_VERTEX_ARRAY);
	//   glDisableClientState(GL_NORMAL_ARRAY);

	/// 
	if (mycloth != 0) {
		glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_LIGHTING);

		GLdouble mMtx[16], pMtx[16];  /// modelview/projection matrix
		GLint viewport[4];         /// the viewport
		glGetIntegerv(GL_VIEWPORT, viewport);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glGetDoublev(GL_MODELVIEW_MATRIX, mMtx);
		glGetDoublev(GL_PROJECTION_MATRIX, pMtx);

		glLineWidth(1.0);
		PolyMesh::ConstFaceIter f_it = mycloth->faces_sbegin();
		while (f_it != mycloth->faces_end()) {
			PolyMesh::ConstFaceVertexIter fv_it(mycloth->cfv_iter(f_it));
			if (mycloth->getType(f_it) == 'B') glColor4f(1.0, 1.0, 0.0, 1.0);
			if (mycloth->getType(f_it) == 'I') glColor4f(1.0, 0.0, 0.0, 1.0);
			if (mycloth->getType(f_it) == 'O') glColor4f(0.0, 1.0, 0.0, 1.0);
			if (mycloth->getType(f_it) == 'X') glColor4f(0.5, 0.5, 0.5, 0.5);
			glBegin(GL_POLYGON);
			do {
				const Point3 &p = mycloth->point(fv_it);
				gluUnProject(p[0], p[1], 1.0, mMtx, pMtx, viewport, &xwini, &ywini, &zwini);
				glVertex3f(xwini, ywini, 0.0);
				//               glVertex3f(p[0], p[1], 0.0);
			} while (++fv_it);
			glEnd();

			fv_it = mycloth->cfv_iter(f_it);
			glColor4f(1.0, 0.5, 0.0, 1.0);
			glBegin(GL_LINE_LOOP);
			do {
				const Point3 &p = mycloth->point(fv_it);
				gluUnProject(p[0], p[1], 1.0, mMtx, pMtx, viewport, &xwini, &ywini, &zwini);
				glVertex3f(xwini, ywini, 0.0);
				//            glVertex3f(p[0], p[1], 0.0);
			} while (++fv_it);
			glEnd();
			++f_it;
		}
		glPopAttrib();
	}

	stroke.draw();

	/// FPS info
	frame++;
	std::time(&end);
	double ddifftime = std::difftime(end, start);
	if (ddifftime > 1) {
		cfps = frame / ddifftime;
		start = end;
		frame = 0;
		ostringstream s1;
		s1 << cfps;
		fps = "FPS: " + s1.str();
	}

	renderInfo();

	glutSwapBuffers();
}

/// Mouse press function
void mouseclick(int button, int state, int x, int y) {

	xini = x; yini = winh - y;
	if (state == GLUT_DOWN) {
		buttonpressed = button;
		if (buttonpressed == GLUT_MIDDLE_BUTTON) {
			arcball.click(xini, yini);
		}
		else if (buttonpressed == GLUT_RIGHT_BUTTON) {
			screenToWorld(x, winh - y, xwini, ywini, zwini);
		}
		else if (buttonpressed == GLUT_LEFT_BUTTON) {
			stroke.clear();
			screenToWorld(x, winh - y, xwini, ywini, zwini);
			grabStart = Point3(xwini, ywini, zwini);
			grabLast = grabStart;
		}
		sceneIniT = sceneT;
	}
	else if (state == GLUT_UP) {
		if (buttonpressed == GLUT_LEFT_BUTTON)
			if (stroke.size() > 3) {
			stroke.subsampleStroke();

			Point2 bbMin = stroke.getBBmin();
			Point2 bbMax = stroke.getBBmax();
			screenToWorld(bbMin[0], bbMin[1], xwini, ywini, zwini);
			Point3 bbmin(xwini, ywini, 0.0);
			screenToWorld(bbMax[0], bbMax[1], xwini, ywini, zwini);
			Point3 bbmax(xwini, ywini, 0.0);

			//if (mycloth != 0) delete mycloth;
			//mycloth = new Cloth("mesh1.off");
			//mycloth = new Cloth;
			//mycloth->intersect(stroke);
			}
	}
	glutPostRedisplay();
}

/// Mouse move function
void mousemove(int x, int y) {

	if (buttonpressed == GLUT_MIDDLE_BUTTON) { // rotations handler
		AMatrix<float> mT;
		arcball.drag(x, winh - y, &mT);
		sceneT = mT*sceneIniT;
	}
	else if (buttonpressed == GLUT_RIGHT_BUTTON) { // translations handler
		double xw, yw, zw;
		screenToWorld(x, winh - y, xw, yw, zw);
		AMatrix<float> mT;
		mT.identity();
		mT.translation(xw - xwini, yw - ywini, zw - zwini);
		sceneT = mT*sceneIniT;
	}
	else if (buttonpressed == GLUT_LEFT_BUTTON) {
		Point2 p(x, winh - y);
		stroke.push_back(p);
	}
	glutPostRedisplay();
}

/// OpenGL key press function
void keypress(unsigned char key, int x, int y) {

	switch (key) {
	case 27: {
		if (mycloth != 0) delete mycloth; exit(1);
	}
	case 'p':
	case 'P': rendermode = POINTS; break;
	case 'w':
	case 'W': rendermode = WIREFRAME; break;
	case 'h':
	case 'H': rendermode = HIDDENLINES; break;
	case 'l':
	case 'L': rendermode = FLATLINES; break;
	case 'f':
	case 'F': rendermode = FLAT; break;
	case 's':
	case 'S': rendermode = SMOOTH; break;

	case '+': { zoomIn(); break; }
	case '-': { zoomOut(); break; }
	}
	glutPostRedisplay();
}

/// Returns the world coordinates of a point in screen space
void screenToWorld(int x, int y, double &xw, double &yw, double &zw) {

	GLdouble mMtx[16], pMtx[16];  /// modelview/projection matrix
	GLint viewport[4];         /// the viewport
	glGetIntegerv(GL_VIEWPORT, viewport);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glGetDoublev(GL_MODELVIEW_MATRIX, mMtx);
	glGetDoublev(GL_PROJECTION_MATRIX, pMtx);

	gluUnProject(x, y, 1.0, mMtx, pMtx, viewport, &xw, &yw, &zw);
}

/// OpenGL reshape funcion
void reshape(int w, int h) {

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1.2, 1.2, -1.2, 1.2, -4, 4);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	winw = w; winh = h;
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glutPostRedisplay();

	arcball.setBounds(winw, winh);
}

/// OpenGL initializations
void init() {

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);

	// lighting setup
	GLfloat mat_A[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat mat_D[] = { 0.5, 0.0, 0.75, 1.0 };
	GLfloat mat_S[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat shine[] = { 128.0 };
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_A);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_D);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_S);
	glMaterialfv(GL_FRONT, GL_SHININESS, shine);

	GLfloat left_light_position[] = { 0.0f, 2.0f, 2.0f, 0.0f };
	GLfloat right_light_position[] = { 0.0f, -2.0f, 2.0f, 0.0f };
	GLfloat left_diffuse_light[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat right_diffuse_light[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	glLightfv(GL_LIGHT0, GL_POSITION, left_light_position);
	glLightfv(GL_LIGHT1, GL_POSITION, right_light_position);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, left_diffuse_light);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, right_diffuse_light);

	glClearColor(1, 1, 1, 1);
	glColor4f(1, 1, 0, 1);
	rendermode = SMOOTH;
	sceneT.identity();

	std::time(&start);
	frame = 0;
}

/// Implements the "zoomin" operation. Narrows the view of the scene.
void zoomIn() {

	AMatrix<float> mT;
	mT.identity();
	mT.scaling(1.2);
	sceneT = mT*sceneT;
}

/// Implements the "zoomout" operation. Enlarges the view of the scene.
void zoomOut() {

	AMatrix<float> mT;
	mT.identity();
	mT.scaling(1.0 / 1.2);
	sceneT = mT*sceneT;
}

/// Renders a string
void renderBitmapString(float x, float y, float z, void *font, const char *string) {
	const char *c;
	glRasterPos3f(x, y, z);
	for (c = string; *c != '\0'; c++) {
		glutBitmapCharacter(font, *c);
	}
}

/// Renders info about the model
void renderInfo() {

	glLoadIdentity();
	glDisable(GL_LIGHTING); glDisable(GL_DEPTH_TEST);
	glColor4f(0, 0, 0, 1);
	renderBitmapString(-0.9, 0.90, 0, (void *)fontsmall, nvertices.c_str());
	renderBitmapString(-0.9, 0.85, 0, (void *)fontsmall, nfaces.c_str());
	renderBitmapString(-0.9, 0.80, 0, (void *)fontsmall, fps.c_str());
	//renderBitmapString(-0.9,0.75,0,(void *)fontsmall, options.c_str());
	glEnable(GL_LIGHTING); glEnable(GL_DEPTH_TEST);
}

/// Process the program's parameters
void input(int argc, char **argv) {

	if (mycloth != 0) delete mycloth;
	mycloth = new Cloth("mesh1.off");

	mycloth->request_face_status();
	mycloth->request_vertex_status();
	mycloth->request_edge_status();
	mycloth->request_halfedge_status();

	cout << "Hello OpenMesh. The mesh has been loaded!" << endl;
	cout << "The mesh has: " << mycloth->n_vertices() << " vertices" << endl;
	cout << "The mesh has: " << mycloth->n_faces() << " faces" << endl;
	cout << "The mesh has: " << mycloth->n_halfedges() << " halfedges" << endl;
	cout << "The mesh has: " << mycloth->n_edges() << " edges" << endl;

	{
		//   Cloth::FaceIter f_it = mycloth->faces_sbegin();
		//   mycloth->split(f_it, Point3(400,400,0));
	}

{
	/*
	Cloth::EdgeIter e_it = mycloth->edges_sbegin();
	Cloth::EdgeHandle e_h0 = e_it.handle();
	++e_it; ++e_it;
	Cloth::EdgeHandle e_h2 = e_it.handle();
	//   mycloth->split_edge(e_h0, Point3(300,100,0));
	//   mycloth->split_edge(e_h2, Point3(300,600,0));
	mycloth->insert_edge(e_h0, e_h2, Point3(350,100,0), Point3(350,600,0));
	*/
}

{
	Cloth::VertexIter v_it = mycloth->vertices_begin();
	while (v_it != mycloth->vertices_end()) {
		cout << " vertex: " << v_it.handle() << endl;
		Cloth::VertexHandle v_h = v_it;
		Cloth::VertexEdgeIter ve_it = mycloth->ve_iter(v_h);
		do {
			cout << "   edge: " << ve_it.handle() << endl;
			//EdgeHandle eh = e_it;
			//Point3 p0 = point(to_vertex_handle(halfedge_handle(eh, 0)));
			//Point3 p1 = point(to_vertex_handle(halfedge_handle(eh, 1)));
		} while (++ve_it);
		++v_it;
	}
}

{
	cout << "inicio ..." << endl;
	Segment<Point2> s(Point2(50, 310), Point2(550, 310));
	Cloth::FaceIter f_it = mycloth->faces_sbegin();
	while (f_it != mycloth->faces_end()) {
		Cloth::FaceHandle fh = f_it;
		cout << "   aaa ..." << endl;
		Cloth::VertexHandle vh = mycloth->intersect(fh, s);
		cout << "   bbb ...: " << mycloth->point(vh) << endl;

		//      Cloth::VertexEdgeIter ve_it = mycloth->ve_iter(vh);
		//      cout<<"   ... :::: "<<ve_it.handle()<<endl;
		//      cout<<"      is_boundary() ...: "<<mycloth->is_boundary(ve_it)<<endl;            

		//mycloth->setSplitFlagOnEdges(vh, true);

		cout << "   ccc ..." << endl;
		++f_it;
		break;
	}
	cout << "fin ..." << endl;
}

	//      Cloth::HalfedgeHandle eh0 = mycloth->halfedge_handle(e_h, 0);
	//      const Point3 &p0 = mycloth->point(mycloth->to_vertex_handle(eh0));
	//      cout<<"p0: "<<p0<<endl;

	/*
	{
	Cloth::VertexIter v_it = mycloth->vertices_begin();
	while (v_it != mycloth->vertices_end()) {
	cout<<" vertex: "<<v_it.handle()<<endl;
	Cloth::VertexHandle v_h = v_it;
	Cloth::VertexEdgeIter ve_it = mycloth->ve_iter(v_h);
	do {
	cout<<"   edge: "<<ve_it.handle()<<endl;
	//EdgeHandle eh = e_it;
	//Point3 p0 = point(to_vertex_handle(halfedge_handle(eh, 0)));
	//Point3 p1 = point(to_vertex_handle(halfedge_handle(eh, 1)));
	} while (++ve_it);
	++v_it;
	}
	}
	*/

{
	//   Cloth::VertexHandle v0, v1, v2, v3;
	//   Cloth::HalfedgeHandle hh = mycloth->vertex_split(v0, v1, v2, v3);
}

	Cloth::FaceIter f_it = mycloth->faces_sbegin();
	while (f_it != mycloth->faces_end()) {
		cout << "new face" << endl;
		Cloth::ConstFaceVertexIter fv_it(mycloth->cfv_iter(f_it));
		do {
			const Point3 &p = mycloth->point(fv_it);
			cout << "   " << p << endl;
		} while (++fv_it);
		++f_it;
	}

	cout << "Hello OpenMesh. The mesh has been loaded!" << endl;
	cout << "The new mesh has: " << mycloth->n_vertices() << " vertices" << endl;
	cout << "The new mesh has: " << mycloth->n_faces() << " faces" << endl;
	cout << "The new mesh has: " << mycloth->n_halfedges() << " halfedges" << endl;
	cout << "The new mesh has: " << mycloth->n_edges() << " edges" << endl;


	/*
	typedef OpenMesh::TriMesh_ArrayKernelT<>   TriMesh;
	TriMesh trimesh;

	if (OpenMesh::IO::read_mesh(trimesh, "meshT1.off")) {
	cout<<"xxx has: "<<trimesh.n_vertices()<<" vertices"<<endl;
	cout<<"xxx has: "<<trimesh.n_faces()<<" faces"<<endl;
	cout<<"xxx has: "<<trimesh.n_halfedges()<<" halfedges"<<endl;
	cout<<"xxx has: "<<trimesh.n_edges()<<" edges"<<endl;

	TriMesh::EdgeIter e_it = trimesh.edges_sbegin();
	TriMesh::EdgeHandle e_h = e_it.handle();
	TriMesh::HalfedgeHandle h0 = trimesh.halfedge_handle(e_h, 0);
	TriMesh::HalfedgeHandle h1 = trimesh.halfedge_handle(e_h, 1);

	if (trimesh.is_boundary(h0)) cout<<"h0 boundary"<<endl;
	if (trimesh.is_boundary(h1)) cout<<"h1 boundary"<<endl;

	const Point3 &p0 = trimesh.point(trimesh.to_vertex_handle(h0));
	const Point3 &p1 = trimesh.point(trimesh.to_vertex_handle(h1));

	cout<<"p0: "<<p0<<endl;
	cout<<"p1: "<<p1<<endl;

	//      TriMesh::EdgeIter e_it = trimesh.edges_sbegin();
	//      TriMesh::EdgeHandle e_h = e_it.handle();
	//      trimesh.split(e_h, Point3(100,400,0));

	cout<<"yyy has: "<<trimesh.n_vertices()<<" vertices"<<endl;
	cout<<"yyy has: "<<trimesh.n_faces()<<" faces"<<endl;
	cout<<"yyy has: "<<trimesh.n_halfedges()<<" halfedges"<<endl;
	cout<<"yyy has: "<<trimesh.n_edges()<<" edges"<<endl;

	//trimesh.split();
	}
	*/
}

// The main function
int main(int argc, char **argv) {

	// Process the input parameters
	input(argc, argv);

	// Glut setup
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(winw, winh);
	glutInitWindowPosition(150, 0);
	glutCreateWindow("Prototype for wearing characters using OpenMesh");

	// OpenGL and data initialization
	init();

	// Callbacks registration
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutMouseFunc(mouseclick);
	glutMotionFunc(mousemove);
	glutKeyboardFunc(keypress);
	glutMainLoop();

	return 0;
}
