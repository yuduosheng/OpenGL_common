#include "camera.h"
#include "object.h"
#include <png.h>

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
				shaderCode +=  Line;
				shaderCode += "\n";
			}
				
			codeStream.close();
		}

		result = glCreateShader(shader_type);

		// Compile Vertex Shader
		cout << "Compiling shader : " << filename << endl;
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
	GLuint                  mVBuffer;
	GLuint                  mIBuffer;
	GLuint                  vao;
	GLuint                  program;
	GLint                   mvp_matrix;

	bool                    bOne;
	bool                    bTwo;
	bool                    bThree;
	bool                    bFour;
	bool                    bFive;
	object                  one;
	object                  two;
	object                  three;
	object                  four;
	object                  five;
	TrackballCamera         mCamera;
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

glfwTest::glfwTest() : App(), mVBuffer(0),
mIBuffer(0), vao(0), program(0), mCamera(mWidth, mHeight),
bOne(true), bTwo(false), bThree(false), bFour(false), bFive(false)
{
	
}

glfwTest::~glfwTest()
{
}

bool glfwTest::Init()
{
	if (!App::Init())
		return false;
	buildShader();
	buildGeometryBuffers();
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
	
	/*vec3 pos = vec3(1.0f, 1.0f, 1.0f);
	vec3 target = vec3(0.0f);
	vec3 up = vec3(0.0f, 1.0f, 0.0f);

	mat4 Projection = perspective(45.0f, (float)mWidth/mHeight, 0.1f, 100.f);
	mat4 View = lookAt(pos, target, up);
	mat4 Model = scale(mat4(1.0f), vec3(0.5f)); 
	mat4 MVP = Projection * View * Model; */
    mat4 MVP = mCamera.getMVP();
	glUniformMatrix4fv(mvp_matrix, 1, GL_FALSE, glm::value_ptr(MVP));
}

void glfwTest::Rendering()
{
	//glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_SHORT, 0);
	//glDrawElements(GL_TRIANGLES, 24, GL_UNSIGNED_SHORT, 0);
	//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//glPointSize(4);
	
	if (bOne)
	    one.render();
	if (bTwo)
		two.render();
	if(bThree)
		three.render();
	if(bFour)
		four.render();
	if (bFive)
		five.render();
	
}
void glfwTest::onMouseWheel(GLFWwindow* window, double x, double y)
{

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
	if ((key == GLFW_KEY_W) && (action == GLFW_PRESS))
	{

	}
	if ((key == GLFW_KEY_1) && (action == GLFW_PRESS))
	{
		cout << "press 1." << endl;
		if (bOne)
		{
			bOne = false;
			bTwo = true;
			
		}
		else
		if (bTwo)
		{
			bThree = true;
			bTwo = false;
		}
		else
		if (bThree)
		{
			bFour = true;
			bThree = false;
		}
		else
		if (bFour)
		{
			bFive = true;
		    bFour = false;
		}
		else
		if (bFive)
		{
			bOne = true;
			bFive = false;
		}

	}
		
}
void glfwTest::buildGeometryBuffers()
{	
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);
	one.readFile("objects/r01.obj");
	two.readFile("objects/r02.obj");
	three.readFile("objects/r03.obj");
	four.readFile("objects/r04.obj");
	five.readFile("objects/r05.obj");
	
	/*
	// Create vertex buffer
	static const GLushort vertex_indices[] =
	{
		0, 1, 2,
		2, 1, 3,
		2, 3, 4,
		4, 3, 5,
		4, 5, 6,
		6, 5, 7,
		6, 7, 0,
		0, 7, 1,
		6, 0, 2,
		2, 4, 6,
		7, 5, 3,
		7, 3, 1
	};

	static const GLfloat vertex_positions[] =
	{
		-0.25f, -0.25f, -0.25f,
		-0.25f, 0.25f, -0.25f,
		0.25f, -0.25f, -0.25f,
		0.25f, 0.25f, -0.25f,
		0.25f, -0.25f, 0.25f,
		0.25f, 0.25f, 0.25f,
		-0.25f, -0.25f, 0.25f,
		-0.25f, 0.25f, 0.25f,
	};

	glGenBuffers(1, &mVBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, mVBuffer);
	glBufferData(GL_ARRAY_BUFFER,
		sizeof(vertex_positions),
		vertex_positions,
		GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(0);

	glGenBuffers(1, &mIBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mIBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		sizeof(vertex_indices),
		vertex_indices,
		GL_STATIC_DRAW);
		*/
	/*
	static const GLushort vertex_indices[] =
	{
        0, 0, 3,
        0, 4, 3,
        0, 3, 5,
        0, 5, 1,
        1, 2, 4,
        1, 5, 2,
        2, 3, 4,
        2, 5, 3,
	};

	static const GLfloat vertex_positions[] =
	{
	    -0.70710678, -0.70710678,  0.00000000,
	    -0.70710678,  0.70710678,  0.00000000,
	    0.70710678,  0.70710678,  0.00000000,
	    0.70710678, -0.70710678,  0.00000000,
	    0.00000000,  0.00000000, -1.00000000,
	    0.00000000,  0.00000000,  1.00000000,
	};
	glGenBuffers(1, &mVBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, mVBuffer);
	glBufferData(GL_ARRAY_BUFFER,sizeof(vertex_positions),vertex_positions,GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(0);

	glGenBuffers(1, &mIBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mIBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		sizeof(vertex_indices),
		vertex_indices,
		GL_STATIC_DRAW);
	*/
	// glFrontFace(GL_CW);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
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

	mvp_matrix = glGetUniformLocation(program, "mvp");
}