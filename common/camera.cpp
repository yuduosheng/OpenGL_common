#include "camera.h"

TrackballCamera::TrackballCamera(float w, float h)
{
	SetWindow(w, h);
	vec3 mCenter = vec3(0.0f);
	float mRadius = 1.0f;
	quat mQuat = quat();

	vec3 pos = vec3(.0f, 0.0f, 1.0f);
	vec3 target = vec3(0.0f);

	setMmworld();
	//mmWorld = scale(mat4(1.0f), vec3(0.5f));
	SetView(pos, target);
	SetProj(45.0f, (float)windowWidth / windowHeight, 0.1f, 100.f);

	mbMouseLButtonDown = false;
	mbMouseWheelRoll = false;
	mbMouseRButtonDown = false;
}

void TrackballCamera::SetView(vec3 eye, vec3 target)
{
	vec3 up(0.0f, 1.0f, 0.0f);
	mView = lookAt(eye, target, up);
}

void TrackballCamera::SetProj(float fFov, float fAspect, float fNearPlane, float fFarPlane)
{
	mProj = perspective(fFov, fAspect, fNearPlane, fFarPlane);
}

void TrackballCamera::computeQuat()
{
	float d, ar, as;
	vec3 vMove, axis;
    
	vMove.x = (curMousePosition.x - preMousePosition.x) ;
	vMove.y = (curMousePosition.y - preMousePosition.y) ;
	vMove.z = 0.0f;
	d = vMove.length();
	if (d != 0)
	{
		ar = d * 3.14159;
		as = sin(ar) / d;
		mQuat = angleAxis(ar, vec3(vMove.y * as, vMove.x * as, 0.0f));
	}
	mQuat = angleAxis(10.0f, vec3(0.0f, 0.0f, 1.0f));
	preMousePosition = curMousePosition;
}