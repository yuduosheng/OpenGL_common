#ifndef APP_H
#define APP_H

#if defined(WIN32)
#include <windows.h>
#endif

#include <cstdio>
#include <iostream>
using namespace std;

#include <GLFW/glfw3.h>
#include <string>
using namespace std;

class App
{
public:
	App();
	virtual ~App();

	float AspectRatio()const;

	int Run();
	// Framework methods.  Derived client class overrides these methods to 
	// implement specific application requirements.

	virtual bool Init();
	virtual void UpdateScene();
	virtual void Rendering();
	//callbacks
	virtual void onResize(GLFWwindow* window, int w, int h);
	virtual void onMouseWheel(GLFWwindow* window, double x, double y);
	virtual void onMouseMove(GLFWwindow* window, double xd, double yd);
	virtual void onMouseButton(GLFWwindow* window, int button, int action, int mods);
	virtual void onKey(GLFWwindow* window, int key, int scancode, int action, int mods);

protected:
	bool InitGLFW();
	bool InitGL();

private:
	string mTitle;
	GLFWwindow* window;
	int mWidth;
	int mHeight;
	bool isTransparency;
	bool isLineSmooth;
};

#endif // APP_H