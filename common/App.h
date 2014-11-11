#ifndef APP_H
#define APP_H

#if defined(WIN32)
#include <windows.h>
#endif

#include <cstdio>
#include <iostream>
using namespace std;

#include <GLFW/glfw3.h>
#include "GL\gl3w.h"
#include <string>
using namespace std;

#include <vmath.h>

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

protected:
	string              mTitle;
	GLFWwindow          *window;
	int                 mWidth;
	int                 mHeight;
	bool                isTransparency;
	bool                isLineSmooth;
};

namespace shader
{

	extern
		GLuint load(const char * filename, GLenum shader_type, bool check_errors)
	{
		GLuint result = 0;
		FILE * fp;
		size_t filesize;
		char * data;

		fp = fopen(filename, "rb");

		if (!fp)
			return 0;

		fseek(fp, 0, SEEK_END);
		filesize = ftell(fp);
		fseek(fp, 0, SEEK_SET);

		data = new char[filesize + 1];

		if (!data)
			goto fail_data_alloc;

		fread(data, 1, filesize, fp);
		data[filesize] = 0;
		fclose(fp);

		result = glCreateShader(shader_type);

		if (!result)
			goto fail_shader_alloc;

		glShaderSource(result, 1, &data, NULL);

		delete[] data;

		glCompileShader(result);

		if (check_errors)
		{
			GLint status = 0;
			glGetShaderiv(result, GL_COMPILE_STATUS, &status);

			if (!status)
			{
				char buffer[4096];
				glGetShaderInfoLog(result, 4096, NULL, buffer);
#ifdef _WIN32
				OutputDebugStringA(filename);
				OutputDebugStringA(":");
				OutputDebugStringA(buffer);
				OutputDebugStringA("\n");
#else
				fprintf(stderr, "%s: %s\n", filename, buffer);
#endif
				goto fail_compile_shader;
			}
		}

		return result;

	fail_compile_shader:
		glDeleteShader(result);

	fail_shader_alloc:;
	fail_data_alloc:
		return result;
	}

}

#endif // APP_H