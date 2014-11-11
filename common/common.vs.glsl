#version 430 core

uniform mat4 mvp;
uniform int vid_offset = 0;

void main(void)
{
    const vec4 vertices[] = vec4[](vec4(-0.5, -0.5, 0.0, 1.0),
                                   vec4( 0.5, -0.5, 0.0, 1.0),
                                   vec4( 0.5,  0.5, 0.0, 1.0),
                                   vec4(-0.5,  0.5, 0.0, 1.0));

    gl_Position = mvp * vertices[(gl_VertexID + vid_offset) % 4];
}
