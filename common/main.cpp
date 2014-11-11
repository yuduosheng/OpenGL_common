#include "APP.h"

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
	GLuint buildShader();
private:
	GLuint                  mVBuffer;
	GLuint                  mIBuffer;
	GLuint                  vao;
	GLuint                  program;
	GLint                   mvp_matrix;
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
mIBuffer(0), vao(0), program(0)
{
}

glfwTest::~glfwTest()
{
}

bool glfwTest::Init()
{
	if (!App::Init())
		return false;
	//program = buildShader();
	return true;
}

void glfwTest::onResize(GLFWwindow* window, int w, int h)
{
	App::onResize(window, w, h);
}

void glfwTest::UpdateScene()
{
	float ratio;
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	ratio = width / (float)height;
	glViewport(0, 0, width, height);
	glClear(GL_COLOR_BUFFER_BIT);

	/*
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-ratio, ratio, -1.f, 1.f, 1.f, -1.f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRotatef((float)glfwGetTime() * 50.f, 0.f, 0.f, 1.f);
    */
	glUseProgram(program);

	vmath::vec3 pos = vmath::vec3(5.0f, 5.0f, 5.0f);
	vmath::vec3 target = vmath::vec3(0.0f);
	vmath::vec3 up = vmath::vec3(0.0f, 1.0f, 0.0f);

    vmath::mat4 model = vmath::mat4::identity();
	vmath::mat4 view = vmath::lookat(pos, target, up);
	vmath::mat4 proj = vmath::perspective(60.0f, (float)width / (float)height, 0.1f, 1000.0f);
	vmath::mat4 mvp = model * view * proj;

	glUniformMatrix4fv(mvp_matrix, 1, GL_FALSE, mvp);	
}

void glfwTest::Rendering()
{
	glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_SHORT, 0);
}
void glfwTest::onMouseWheel(GLFWwindow* window, double x, double y)
{

}
void glfwTest::onMouseMove(GLFWwindow* window, double xd, double yd)
{

}
void glfwTest::onMouseButton(GLFWwindow* window, int button, int action, int mods)
{

}
void glfwTest::onKey(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	App::onKey(window, key, scancode, action, mods);
}
void glfwTest::buildGeometryBuffers()
{
	// Create vertex buffer
	static const GLfloat vertices[] =
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
	static const GLushort indices[] =
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

	glGenBuffers(1, &mVBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, mVBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(0);

	glGenBuffers(1, &mIBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mIBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
}
GLuint glfwTest::buildShader()
{
	GLuint vs;
	GLuint fs;

	vs = shader::load("common.vs.glsl", GL_VERTEX_SHADER, true);
	fs = shader::load("conmon.fs.glsl", GL_FRAGMENT_SHADER, true);

	if (program)
		glDeleteProgram(program);

	program = glCreateProgram();
	glAttachShader(program, vs);
	glAttachShader(program, fs);
	glLinkProgram(program);

	mvp_matrix = glGetUniformLocation(program, "mvp");
}