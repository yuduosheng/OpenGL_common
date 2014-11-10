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
};

int main(void)
{
	glfwTest theApp;

	if (!theApp.Init())
		return 0;

	return theApp.Run();
}

glfwTest::glfwTest() : App()
{
}

glfwTest::~glfwTest()
{
}

bool glfwTest::Init()
{
	if (!App::Init())
		return false;

	return true;
}

void glfwTest::onResize(GLFWwindow* window, int w, int h)
{
	App::onResize(window, w, h);
}

void glfwTest::UpdateScene()
{
}

void glfwTest::Rendering()
{
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