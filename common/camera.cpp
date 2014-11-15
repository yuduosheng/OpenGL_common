#include "camera.h"

TrackballCamera::TrackballCamera(float w, float h)
{
	SetWindow(w, h);
	vec3 mCenter = vec3(0.0f);
	float mRadius = 1.0f;
	quat mQuat = quat();

	vec3 pos = vec3(0.0f, 0.0f, 2.0f);
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
	vec3 preVec, curVec, vMove, axis;

	preVec.x = (2.0f * preMousePosition.x / windowWidth - 1.0f) / mProj[0][0];
	preVec.y = (-2.0f * preMousePosition.y / windowHeight + 1.0f) / mProj[1][1];
	preVec.z = 1.0f;
	d = preVec.x * preVec.x + preVec.y * preVec.y;

	curVec.x = (2.0f * curMousePosition.x / windowWidth - 1.0f) / mProj[0][0];
	curVec.y = (-2.0f * curMousePosition.y / windowHeight + 1.0f) / mProj[1][1];
	curVec.z = 1.0f;

	vMove = curVec - preVec;

	axis = cross(curVec, preVec);
	axis = normalize(axis);

	as = (length(vMove) * length(vMove) - length(curVec) * length(curVec) - length(preVec) * length(preVec)) / (2.0f*length(curVec)*length(preVec));

	mQuat = angleAxis(as, axis);

	preMousePosition = curMousePosition;
}