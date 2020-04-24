#version 410

layout (location = 0) in vec3 vtxPosition;
layout (location = 1) in vec3 vtxNormal;
// layout (location = 2) in vec2 vtxUV;

uniform mat4 proj_mat, model_mat, view_mat;

out vec3 normal;
out vec3 frag;

// out vec2 UV;

void main(){
  normal =  mat3(transpose(inverse(model_mat))) * vtxNormal;
  frag = vec3(model_mat * vec4(vtxPosition, 1.0));

  // UV = vtxUV;

  gl_Position = proj_mat * view_mat * model_mat * vec4(vtxPosition, 1.0);

}
