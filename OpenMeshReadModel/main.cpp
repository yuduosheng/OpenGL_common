#include <png.h>
#include "App.h"
#include "camera.h"
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

namespace shader
{
	GLuint loadShader(const char * filename, GLenum shader_type, bool check_errors)
	{
		GLuint result = 0;
		string shaderCode = "";
		ifstream codeStream(filename, ios::in);

		if (codeStream.is_open())
		{
			string Line = "";
			while (getline(codeStream, Line))
			{
				shaderCode += Line;
				shaderCode += "\n";
			}

			codeStream.close();
		}

		result = glCreateShader(shader_type);

		// Compile Vertex Shader
		//cout << "Compiling shader : " << filename << endl;
		char const * codePointer = shaderCode.c_str();

		glShaderSource(result, 1, &codePointer, NULL);
		glCompileShader(result);

		if (check_errors)
		{
			GLint status = 0;
			glGetShaderiv(result, GL_COMPILE_STATUS, &status);

			if (!status)
			{
				char buffer[4096];
				glGetShaderInfoLog(result, 4096, NULL, buffer);

				cout << filename << buffer << endl;
			}
		}

		return result;
	}

}

void capture(GLFWwindow *window)
{
	const char filepath[] = "./output.png";
	png_bytep raw1D;
	png_bytepp raw2D;
	int i;
	int width;
	int height;
	glfwGetWindowSize(window, &width, &height);
	// 構造体確保
	FILE *fp = fopen(filepath, "wb");
	png_structp pp = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	png_infop ip = png_create_info_struct(pp);
	// 書き込み準備
	png_init_io(pp, fp);
	png_set_IHDR(pp, ip, width, height,
		8, // 8bit以外にするなら変える
		PNG_COLOR_TYPE_RGBA, // RGBA以外にするなら変える
		PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	// ピクセル領域確保
	raw1D = (png_bytep)malloc(height * png_get_rowbytes(pp, ip));
	raw2D = (png_bytepp)malloc(height * sizeof(png_bytep));
	for (i = 0; i < height; i++)
		raw2D[i] = &raw1D[i*png_get_rowbytes(pp, ip)];
	// 画像のキャプチャ
	glReadBuffer(GL_FRONT);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // 初期値は4
	glReadPixels(0, 0, width, height,
		GL_RGBA, // RGBA以外にするなら変える
		GL_UNSIGNED_BYTE, // 8bit以外にするなら変える
		(void*)raw1D);
	// 上下反転
	for (i = 0; i < height / 2; i++){
		png_bytep swp = raw2D[i];
		raw2D[i] = raw2D[height - i - 1];
		raw2D[height - i - 1] = swp;
	}
	// 書き込み
	png_write_info(pp, ip);
	png_write_image(pp, raw2D);
	png_write_end(pp, ip);
	// 開放
	png_destroy_write_struct(&pp, &ip);
	fclose(fp);
	free(raw1D);
	free(raw2D);

	printf("write out screen capture to '%s'\n", filepath);
}
void printHint()
{
	cout << "Operation:" << endl;
	cout << "1.Click left mouse button and drag for trackball effect." << endl;
	cout << "2.Click right mouse button and drag for translate the object." << endl;
	cout << "3.Press \"1\" to change object." << endl;
	cout << "4.Press \"2\" to change polygon mode." << endl;
	cout << "5.Press \"p\" to get a screenshot." << endl;
}

class glfwTest : public App
{
public:
	glfwTest();
	~glfwTest();

	bool Init();
	void UpdateScene();
	void Rendering();
	void onResize(GLFWwindow* window, int w, int h);

	void onMouseWheel(GLFWwindow* window, double x, double y);
	void onMouseMove(GLFWwindow* window, double xd, double yd);
	void onMouseButton(GLFWwindow* window, int button, int action, int mods);
	void onKey(GLFWwindow* window, int key, int scancode, int action, int mods);

private:
	bool OpenMeshReadFile(const char * filename);
	void PrintMeshStatus()
	{
		cout << "The number of vertices: " << meshVetexNum << endl;
		cout << "The number of faces: " << meshFaceNum << endl;
		cout << "The number of halfedges: " << meshHalfEdgeNum << endl;
		cout << "The number of boundry-edges: " << meshBoundryEdgeNum << endl;
	}
	void buildGeometryBuffers();
	void buildShader();
private:

	GLuint                  vao;
	GLuint                  program;
	GLuint                  mvp_matrix;
	GLuint                  m_matrix;
	GLuint                  v_matrix;
	GLuint                  l_position;

	TrackballCamera         mCamera;
	MyMesh                  mesh;
	GLuint                  meshVBuffer;
	GLuint                  meshVNormal;
	GLuint                  meshFNormal;
	vector<OpenMesh::Vec3f> meshVertexBuffer;//vertex buffer
	vector<OpenMesh::Vec3f> meshVertexNormalBuffer;//vertex normal buffer
	vector<OpenMesh::Vec3f> meshFaceNormalBuffer;//face buffer normal
	GLint                   meshVetexNum = 0;
	GLint                   meshFaceNum = 0;
	GLint                   meshHalfEdgeNum = 0;
	GLint                   meshBoundryEdgeNum = 0;
	enum PolygonMode
	{
		FILL,
		LINE,
		POINT,
		FlatShading,
		SmoothShading
	};
	PolygonMode             mPM;
	PolygonMode             mSM;
};

int main(void)
{
	glfwTest *theApp = new glfwTest;

	if (!theApp->Init())
		return 0;
	theApp->Run();
	delete theApp;
	return 0;
}

glfwTest::glfwTest() : App(), vao(0), program(0), mCamera(800, 600),
mvp_matrix(0),
m_matrix(0),
v_matrix(0),
l_position(0)
{
	mPM = FILL;
	mSM = FlatShading;
}

glfwTest::~glfwTest()
{
	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	// Cleanup VBO and shader
	glDeleteBuffers(1, &meshVBuffer);
	glDeleteBuffers(1, &meshVNormal);
	glDeleteProgram(program);
	glDeleteVertexArrays(1, &vao);
}

bool glfwTest::Init()
{
	if (!App::Init())
		return false;
	printHint();

	OpenMeshReadFile("bun_zipper_res4.ply");
	//OpenMeshReadFile("bun_zipper.ply");

	buildGeometryBuffers();
	buildShader();

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	return true;
}

void glfwTest::onResize(GLFWwindow* window, int w, int h)
{
	App::onResize(window, w, h);
}

void glfwTest::UpdateScene()
{
	static const GLfloat black[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	static const GLfloat one = 1.0f;
	
	glViewport(0, 0, mWidth, mHeight);
	glClearBufferfv(GL_COLOR, 0, black);
	glClearBufferfv(GL_DEPTH, 0, &one);

	glUseProgram(program);
	/*
	glm::vec3 pos = glm::vec3(0.0f, 0.0f, 1.0f);
	glm::vec3 target = glm::vec3(0.0f);
	glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);

	glm::mat4 Projection = glm::perspective(45.0f, (float)mWidth / mHeight, 0.1f, 100.f);
	glm::mat4 View = glm::lookAt(pos, target, up);
	glm::mat4 Model = glm::scale(glm::mat4(1.0f), glm::vec3(2.0f));*/
	//glm::mat4 MVP = Projection * View * Model;
	glm::mat4 MVP = mCamera.getMVP();
	glm::mat4 M = mCamera.getM();
	glm::mat4 V = mCamera.getV();
	glUniformMatrix4fv(mvp_matrix, 1, GL_FALSE, glm::value_ptr(MVP));
	glUniformMatrix4fv(v_matrix, 1, GL_FALSE, glm::value_ptr(V));
	glUniformMatrix4fv(m_matrix, 1, GL_FALSE, glm::value_ptr(M));
	glm::vec3 lightPos = glm::vec3(4, 4, 4);
	glUniform3f(l_position, lightPos.x, lightPos.y, lightPos.z);
}
void glfwTest::Rendering()
{
	if (mSM == SmoothShading)
	{
		glGenBuffers(1, &meshVNormal);
		glBindBuffer(GL_ARRAY_BUFFER, meshVNormal);
		glBufferData(GL_ARRAY_BUFFER, meshVertexNormalBuffer.size() * sizeof(OpenMesh::Vec3f), &meshVertexNormalBuffer[0], GL_STATIC_DRAW);

		// 2rd attribute buffer : normals
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, meshVNormal);
		glVertexAttribPointer(
			1,                                // attribute
			3,                                // size
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			NULL                         // array buffer offset
			);
	}
	if (mSM == FlatShading)
	{ 
	    glGenBuffers(1, &meshFNormal);
	    glBindBuffer(GL_ARRAY_BUFFER, meshFNormal);
	    glBufferData(GL_ARRAY_BUFFER, meshFaceNormalBuffer.size() * sizeof(OpenMesh::Vec3f), &meshFaceNormalBuffer[0], GL_STATIC_DRAW);
	    
	    // 2rd attribute buffer : normals
	    glEnableVertexAttribArray(1);
	    glBindBuffer(GL_ARRAY_BUFFER, meshFNormal);
	    glVertexAttribPointer(
	    	1,                                // attribute
	    	3,                                // size
	    	GL_FLOAT,                         // type
	    	GL_FALSE,                         // normalized?
	    	0,                                // stride
	    	NULL                         // array buffer offset
	    	);
	}
	if (mPM == FILL)
	{
		//glEnable(GL_SMOOTH);
		//glShadeModel(GL_SMOOTH);
		glEnable(GL_FLAT);
		glShadeModel(GL_FLAT);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	else
		if (mPM == LINE)
		{
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}
		else
			if (mPM == POINT)
			{
		glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
		glPointSize(5.0f);
			}
	glDrawArrays(GL_TRIANGLES, 0, meshVertexBuffer.size());
}
void glfwTest::onMouseWheel(GLFWwindow* window, double x, double y)
{
	float scale = 1.0f;
	float mouseWheelScale = 0.1f;
	scale += mouseWheelScale  * (float)y;
	mCamera.setScaleFactor(scale);
	mCamera.setMmworldScle();
}
void glfwTest::onMouseMove(GLFWwindow* window, double xd, double yd)
{
	double x = xd;
	double y = yd;

	if (mCamera.IsMouseLButtonDown())
	{
		glfwGetCursorPos(window, &xd, &yd);
		mCamera.SetCurMousePosition(xd, yd);
		mCamera.computeQuat();
		mCamera.setMmworldQuat();
	}
	if (mCamera.IsMouseRButtonDown())
	{
		glfwGetCursorPos(window, &xd, &yd);
		mCamera.SetCurMousePosition(xd, yd);
		mCamera.computeTran();
		mCamera.setMmworldTran();
	}
}
void glfwTest::onMouseButton(GLFWwindow* window, int button, int action, int mods)
{
	double xd, yd;
	
	if ((button == GLFW_MOUSE_BUTTON_1) && (action == GLFW_PRESS))
	{
		mCamera.SetMouseLButtonStat(true);
		glfwGetCursorPos(window, &xd, &yd);
		mCamera.initMousePosition(xd, yd);

	}
	else if ((button == GLFW_MOUSE_BUTTON_1) && (action == GLFW_RELEASE))
	{
		mCamera.SetMouseLButtonStat(false);
	}

	if ((button == GLFW_MOUSE_BUTTON_2) && (action == GLFW_PRESS))
	{
		mCamera.SetMouseRButtonStat(true);
		glfwGetCursorPos(window, &xd, &yd);
		mCamera.initMousePosition(xd, yd);

	}
	else if ((button == GLFW_MOUSE_BUTTON_2) && (action == GLFW_RELEASE))
	{
		mCamera.SetMouseRButtonStat(false);
	}
}
void glfwTest::onKey(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	App::onKey(window, key, scancode, action, mods);
	//capture
	if ((key == GLFW_KEY_P) && (action == GLFW_PRESS))
	{
		capture(window);
		cout << "capture over." << endl;
	}
	if ((key == GLFW_KEY_S) && (action == GLFW_PRESS))
	{
		PrintMeshStatus();
	}
	if ((key == GLFW_KEY_3) && (action == GLFW_PRESS))
	{
		if (mSM == FlatShading)
			mSM = SmoothShading;
		else
			if (mSM == SmoothShading)
			{
			mSM = FlatShading;
			}
	}
	if ((key == GLFW_KEY_2) && (action == GLFW_PRESS))
	{
		if (mPM == FILL)
			mPM = LINE;
		else
			if (mPM == LINE)
			{
			mPM = POINT;
			}
			else
				if (mPM == POINT)
				{
			mPM = FILL;
				}
	}
	if ((key == GLFW_KEY_1) && (action == GLFW_PRESS))
	{
		//cout << "press 1." << endl;

	}

}
bool glfwTest::OpenMeshReadFile(const char * filename)
{
	
	// request vertex normals, so the mesh reader can use normal information
	// if available
	mesh.request_vertex_normals();
	// request face normals,
	mesh.request_face_normals();
	// assure we have vertex normals
	if (!mesh.has_vertex_normals())
	{
		std::cerr << "ERROR: Standard vertex property 'Normals' not available!\n";
		return 1;
	}

	OpenMesh::IO::Options opt;
	// read mesh from file
	if (!OpenMesh::IO::read_mesh(mesh, filename, opt))
	{
		std::cerr << "Error: Cannot read mesh from " << filename << std::endl;
		return 1;
	}

	// If the file did not provide vertex normals, then calculate them
	if (!opt.check(OpenMesh::IO::Options::VertexNormal))
	{
		// let the mesh update the normals
		mesh.update_normals();
	}

	// Get the face-vertex circulator of face _fh
	// get the vertex nomal
	for (MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		for (MyMesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			meshVertexBuffer.push_back(mesh.point(*fv_it));
			meshVertexNormalBuffer.push_back(mesh.normal(*fv_it));
			meshFaceNormalBuffer.push_back(mesh.normal(*f_it));
		}
	}
	// don't need the normals anymore? Remove them!
	mesh.release_vertex_normals();
	// dispose the face normals, as we don't need them anymore
	mesh.release_face_normals();

	mesh.request_vertex_status();
	meshVetexNum = mesh.n_vertices();
	mesh.request_face_status();
	meshFaceNum = mesh.n_faces();
	mesh.request_halfedge_status();
	meshHalfEdgeNum = mesh.n_halfedges();

	// iterate over all halfedges
	for (MyMesh::HalfedgeIter h_it = mesh.halfedges_begin(); h_it != mesh.halfedges_end(); ++h_it)
	{
		if (!mesh.face_handle(*h_it).is_valid())
		{
			++meshBoundryEdgeNum;
		}
	}

	mesh.release_face_status();
	mesh.release_vertex_status();
	mesh.release_halfedge_status();
}
void glfwTest::buildGeometryBuffers()
{
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	glGenBuffers(1, &meshVBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, meshVBuffer);
	glBufferData(GL_ARRAY_BUFFER,
		meshVertexBuffer.size() * sizeof(OpenMesh::Vec3f),
		&meshVertexBuffer[0],
		GL_STATIC_DRAW);
	// 1rst attribute buffer : vertices
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, meshVBuffer);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);

	/*
	glGenBuffers(1, &meshVNormal);
	glBindBuffer(GL_ARRAY_BUFFER, meshVNormal);
	glBufferData(GL_ARRAY_BUFFER, meshVertexNormalBuffer.size() * sizeof(OpenMesh::Vec3f), &meshVertexNormalBuffer[0], GL_STATIC_DRAW);

	// 2rd attribute buffer : normals
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, meshVNormal);
	glVertexAttribPointer(
		1,                                // attribute
		3,                                // size
		GL_FLOAT,                         // type
		GL_FALSE,                         // normalized?
		0,                                // stride
		NULL                         // array buffer offset
		);
    
	glGenBuffers(1, &meshFNormal);
	glBindBuffer(GL_ARRAY_BUFFER, meshFNormal);
	glBufferData(GL_ARRAY_BUFFER, meshFaceNormalBuffer.size() * sizeof(OpenMesh::Vec3f), &meshFaceNormalBuffer[0], GL_STATIC_DRAW);

	// 2rd attribute buffer : normals
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, meshFNormal);
	glVertexAttribPointer(
		1,                                // attribute
		3,                                // size
		GL_FLOAT,                         // type
		GL_FALSE,                         // normalized?
		0,                                // stride
		NULL                         // array buffer offset
		);*/
}
void glfwTest::buildShader()
{
	GLuint vs;
	GLuint fs;

	vs = shader::loadShader("common.vs.glsl", GL_VERTEX_SHADER, true);
	fs = shader::loadShader("common.fs.glsl", GL_FRAGMENT_SHADER, true);

	if (program)
		glDeleteProgram(program);

	program = glCreateProgram();
	glAttachShader(program, vs);
	glAttachShader(program, fs);

	glLinkProgram(program);

	glDeleteShader(vs);
	glDeleteShader(fs);

	mvp_matrix = glGetUniformLocation(program, "MVP");
	v_matrix = glGetUniformLocation(program, "V");
	m_matrix = glGetUniformLocation(program, "M");
	l_position = glGetUniformLocation(program, "LightPosition_worldspace");
}