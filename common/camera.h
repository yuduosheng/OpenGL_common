#ifndef CAMERA_H
#define CAMERA_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
using namespace glm;
#include <iostream>
using namespace std;

class TrackballCamera
{
public:
	TrackballCamera(float w, float h);

	void                SetView(vec3 pEye, vec3 pLookat);
	void                SetProj(float fFov, float fAspect, float fNearPlane, float fFarPlane);
	void computeQuat();

	mat4 getMVP()                                         { return mProj * mView * mmWorld; }
	void setMmworld()                                     { mmWorld *= mat4_cast(mQuat); };
	void initMousePosition(float x, float y)              { SetCurMousePosition(x, y); SetPreMousePosition(x, y); }
	void SetMouseLButtonStat(bool stat)                   { mbMouseLButtonDown = stat; }
	void SetMouseLWheelStat(bool stat)                    { mbMouseWheelRoll = stat; }
	void SetMouseRButtonStat(bool stat)                   { mbMouseRButtonDown = stat; }

	void SetCurMousePosition(float x, float y)           { curMousePosition.x = x; curMousePosition.y = y; cout << x << " " << y << endl; }
	void SetPreMousePosition(float x, float y)             { preMousePosition.x = x; preMousePosition.y = y; cout << x << " " << y << endl; }

	void SetWindow(float w, float h)                      { windowWidth = w; windowHeight = h; }

	bool IsMouseLButtonDown() const                       { return mbMouseLButtonDown; }
	bool IsMouseWheelRoll() const                         { return mbMouseRButtonDown; }
	bool IsMouseRButtonDown() const                       { return mbMouseRButtonDown; }

protected:
	vec3 mCenter;
	float mRadius;
	quat mQuat;

	mat4 mView;
	mat4 mProj;
	mat4 mmWorld;

	bool mbMouseLButtonDown;    // True if left button is down 
	bool mbMouseWheelRoll;          // True if middle wheel is roll 
	bool mbMouseRButtonDown;    // True if right button is down 

	float windowWidth;
	float windowHeight;
	vec2 curMousePosition;
	vec2 preMousePosition;
};

#endif 