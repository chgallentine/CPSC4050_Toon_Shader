#version 410


// Define INPUTS from fragment shader
uniform mat4 view_mat;

// uniform int DIFFUSE_SHADING, SPECULAR_LIGHTING, DIFFUSE_MAP;
// uniform int SPECULAR_EXPONENT;

// uniform vec3 LIGHT_POSITION;

in vec3 frag;
in vec3 normal;

// These come from the VAO for texture coordinates.
// in vec2 UV;

// And from the uniform outputs for the textures setup in main.cpp.
// uniform sampler2D texture00;
// uniform sampler2D texture01;

out vec4 fragment_color; //RGBA color

void main () {
	vec3 light_color = vec3(0.25,0.25,0.25);

	// Assume following are pointing at origin
	// vec3 light_position = vec3(view_mat[3][0], view_mat[3][1], -view_mat[3][2]);
	vec3 light_position = vec3(0.0,0.0,5.0);
	vec3 view_position = vec3(0.0,0.0,5.0);

	vec3 norm = normalize(normal);

	vec3 light_direction = normalize(light_position - frag);  
	vec3 view_direction = normalize(view_position - frag);
	vec3 half_direction = normalize(light_position + view_direction);

	float diff = max(dot(norm, light_direction), 0.0);
	float lambertian = max(dot(normal, half_direction), 0.0);
	float specular_val = pow(lambertian, 16);

	vec3 diffuse;
	vec3 specular = light_color * specular_val;

	diffuse = diff * light_color;

	fragment_color = vec4(diffuse + specular,1.0);
}


















