#include <png.h>
#include "App.h"
#include "camera.h"
#include "OMmodel.h"

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
		else
		{
			cout << filename <<" can not open! "<< endl;
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
	glFlush();
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

#define	FILL 0x01
#define	LINE 0x02
#define	POINT 0x04

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

	OMmodel                 bunnyMesh;
	//OMmodel                 dragonMesh;
	//OMmodel                 horseMesh;
	OMmodel                 *curModel;
	GLint                   mPM;
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
	curModel = &bunnyMesh;
}

glfwTest::~glfwTest()
{
	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	// Cleanup VBO and shader
	glDeleteProgram(program);
	glDeleteVertexArrays(1, &vao);
}

bool glfwTest::Init()
{
	if (!App::Init())
		return false;
	printHint();

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
	//glm::mat4 M = mCamera.getM();
	//glm::mat4 V = mCamera.getV();
	glUniformMatrix4fv(mvp_matrix, 1, GL_FALSE, glm::value_ptr(MVP));
	//glUniformMatrix4fv(v_matrix, 1, GL_FALSE, glm::value_ptr(V));
	//glUniformMatrix4fv(m_matrix, 1, GL_FALSE, glm::value_ptr(M));
	//glm::vec3 lightPos = glm::vec3(4, 4, 4);
	//glUniform3f(l_position, lightPos.x, lightPos.y, lightPos.z);
}
void glfwTest::Rendering()
{
	if (mPM == FILL)
	{
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
	//curModel->RenderModel();
	curModel->RenderModelWithColor();
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
		curModel->PrintMeshStatus();
	}
	if ((key == GLFW_KEY_3) && (action == GLFW_PRESS))
	{
		curModel->SetSM();
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
		/*
		if (curModel == &bunnyMesh)
			curModel = &dragonMesh;
		else
			if (curModel == &dragonMesh)
				curModel = &horseMesh;
			else
				curModel = &bunnyMesh;*/
	}

}
void glfwTest::buildGeometryBuffers()
{
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);
	bunnyMesh.OpenMeshReadFile("bunny.off");
	//dragonMesh.OpenMeshReadFile("bun_zipper_res4.ply");
	//horseMesh.OpenMeshReadFile("dragon_vrip_res4.ply");
}
void glfwTest::buildShader()
{
	GLuint vs;
	GLuint fs;

	vs = shader::loadShader("ColorVertexShader.glsl", GL_VERTEX_SHADER, true);
	fs = shader::loadShader("ColorFragmentShader.glsl", GL_FRAGMENT_SHADER, true);
	//vs = shader::loadShader("common.vs.glsl", GL_VERTEX_SHADER, true);
	//fs = shader::loadShader("common.fs.glsl", GL_FRAGMENT_SHADER, true);

	if (program)
		glDeleteProgram(program);

	program = glCreateProgram();
	glAttachShader(program, vs);
	glAttachShader(program, fs);

	glLinkProgram(program);

	glDeleteShader(vs);
	glDeleteShader(fs);

	mvp_matrix = glGetUniformLocation(program, "MVP");
	//v_matrix = glGetUniformLocation(program, "V");
	//m_matrix = glGetUniformLocation(program, "M");
	//l_position = glGetUniformLocation(program, "LightPosition_worldspace");
}
