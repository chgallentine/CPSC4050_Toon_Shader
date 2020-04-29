#version 410


// Define INPUTS from fragment shader
uniform mat4 view_mat;

// uniform float DIFFUSE_SHADING, SPECULAR_LIGHTING, DIFFUSE_MAP;
// uniform float SPECULAR_EXPONENT;

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
	vec3 light_color = vec3(0.5,0.5,0.5);
//<<<<<<< patch-1
	
	vec3 unlit_color = vec3(0, 76, 153);
	vec3 outline_color = vec3(0,0,0);
	float lit_thickness = 0.06;
	float unlit_thickness = 0.25;
	
	vec3 rim_color = vec3(255,255,255);
	float rim_amount = 0.2;
	float rim_threshold = 0.06;
	
//=======
	vec3 unlit_color = vec3(0, 25, 51);
	vec3 outline_color = vec3(0,0,0);
	float lit_thickness = 0.05;
	float unlit_thickness = 0.2;
//>>>>>>> master
	vec3 ambient = vec3(100.0/255.0,149.0/255.0,237.0/255.0);
	float ambient_gain = 1.0;
	// Assume following are pointing at origin
	// vec3 light_position = vec3(view_mat[3][0], view_mat[3][1], -view_mat[3][2]);
	
	vec3 light_position = vec3(-5.0,-5.0,3.0);
	vec3 view_position = vec3(0.0,0.0,5.0);

	vec3 norm = normalize(normal);

	vec3 light_direction = normalize(light_position - frag);  
	vec3 view_direction = normalize(view_position - frag);
	vec3 half_direction = normalize(light_position + view_direction);

	float diff = max(dot(norm, light_direction), 0.0);
	float lambertian = max(dot(normal, half_direction), 0.0);
	float specular_val = pow(lambertian, 256);

	vec3 diffuse;
	vec3 specular = light_color * specular_val;

	diffuse = diff * light_color;

	ambient = ambient_gain * ambient;
	
	fragment_color = vec4(ambient + specular,1.0);
	
	float rim_dot = 1-dot(view_direction, norm);

	//Draw Outline

	if(dot(view_direction, norm) < mix(unlit_thickness, lit_thickness, max(0.0, dot(norm, light_direction))))
	{
		fragment_color = vec4((light_color * outline_color), 1.0);
	}
	//Draw Rim Light
	if(dot(view_direction, norm) < mix(0, rim_amount, max(0.0, dot(norm, light_direction))))
	{
		float rim_intensity = rim_dot * pow(dot(light_position, norm), rim_threshold);
		rim_intensity = smoothstep(rim_amount - 0.01, rim_amount + 0.01, rim_intensity);
		vec3 rim = rim_intensity * rim_color;
		fragment_color = vec4((ambient + specular + rim)*rim_color, 1.0);
	}
	// fragment_color = vec4(ambient,1.0);
	
	
}
