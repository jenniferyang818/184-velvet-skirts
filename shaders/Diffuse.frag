#version 330

// The camera's position in world-space
uniform vec3 u_cam_pos;

// Color
uniform vec4 u_color;

// Properties of the single point light
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

// We also get the uniform texture we want to use.
uniform sampler2D u_texture_1;

// These are the inputs which are the outputs of the vertex shader.
in vec4 v_position;
in vec4 v_normal;

// This is where the final pixel color is output.
// Here, we are only interested in the first 3 dimensions (xyz).
// The 4th entry in this vector is for "alpha blending" which we
// do not require you to know about. For now, just set the alpha
// to 1.
out vec4 out_color;

void main() {
  // YOUR CODE HERE
  
  // (Placeholder code. You will want to replace it.)
  vec3 diffuse_coeff;
  diffuse_coeff = vec3(1, 1, 1);
  float distance_squared;
  distance_squared = pow(v_position.x - u_light_pos.x, 2) + pow(v_position.y - u_light_pos.y, 2) + pow(v_position.z - u_light_pos.z, 2);
  vec3 illum;
  illum = u_light_intensity / distance_squared;
  vec3 difference = vec3(u_light_pos.x - v_position.x, u_light_pos.y - v_position.y, u_light_pos.z - v_position.z);
  vec3 max_thing;
  max_thing = max(vec3(0.0, 0.0, 0.0), dot(normalize(v_normal), vec4(normalize(difference), 1)));
  out_color = vec4(diffuse_coeff * illum * max_thing, 1);
  out_color.a = 1;
}
