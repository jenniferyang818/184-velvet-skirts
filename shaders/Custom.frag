#version 330

// (Every uniform is available here.)

uniform mat4 u_view_projection;
uniform mat4 u_model;

uniform float u_normal_scaling;
uniform float u_height_scaling;

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

// Feel free to add your own textures. If you need more than 4,
// you will need to modify the skeleton.
uniform sampler2D u_texture_1;
uniform sampler2D u_texture_2;
uniform sampler2D u_texture_3;
uniform sampler2D u_texture_4;

uniform vec2 u_texture_3_size;

// Environment map! Take a look at GLSL documentation to see how to
// sample from this.
uniform samplerCube u_texture_cubemap;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec2 v_uv;

out vec4 out_color;

float h(vec2 uv) {
  // You may want to use this helper function...
  float height;
  height = texture(u_texture_3, uv).r;
  return height;
}

void main() {
  const float backscatter = 0.25;
  const float edginess = 8.0;
  vec3 sheen = vec3(0.75, 0, 0); //sheen is color of highlight
  const float roughness = 0.15;
  
  // include bump mapping for texture
  mat3 tbn;
  vec3 b;
  b = cross(v_normal.xyz, v_tangent.xyz);
  tbn[0] = normalize(v_tangent.xyz);
  tbn[1] = normalize(b);
  tbn[2] = normalize(v_normal.xyz);
  vec2 du_arg = vec2(v_uv[0] + (1 / u_texture_3_size[0]), v_uv[1]);
  float du = (h(du_arg) - h(v_uv)) * u_height_scaling * u_normal_scaling;
  vec2 dv_arg = vec2(v_uv[0], v_uv[1] + (1 / u_texture_3_size[1]));
  float dv = (h(dv_arg) - h(v_uv)) * u_height_scaling * u_normal_scaling;

  vec3 n_o;
  n_o = vec3(-1 * du, -1 * dv, 1);
  vec3 n_d;
  n_d = tbn * n_o;

  //blinn phong shading
  //vec3 normal_normal = normalize(v_normal.xyz);
  vec3 normal_normal = normalize(n_d);
  vec3 l_vector = normalize(vec3(u_light_pos.x - v_position.x, u_light_pos.y - v_position.y, u_light_pos.z - v_position.z));
  vec3 cam_vector = normalize(vec3(u_cam_pos.x - v_position.x, u_cam_pos.y - v_position.y, u_cam_pos.z - v_position.z));
  
  float distance_squared = pow(v_position.x - u_light_pos.x, 2) + pow(v_position.y - u_light_pos.y, 2) + pow(v_position.z - u_light_pos.z, 2);
  vec3 illum = u_light_intensity / distance_squared;
  
  vec3 diffuse = max(dot(l_vector, normal_normal), 0.0) * illum; //diffuse

  float max_thing = max(dot(l_vector, cam_vector), 0.0);
  float p = 1/roughness;
  vec3 shiny = sheen * pow(max_thing, p) * backscatter * illum; //analyze this
  // TODO: why not using the bisector?
  // illumination: sheen * backscatter (backscatter is from other pts)
  // backscatter is over light to cam

  float cosine = max(dot(normal_normal, cam_vector), 0.0);
  float sine = sqrt(1.0 - cosine);
  p = edginess;
  shiny = shiny + sheen * pow(sine, p) * diffuse;
  // p is edginess here, determines how sharp the highlights are going to be
  // here: illumination is sheen * diffuse
  // increasing p makes specular highlight smaller and harsher
  // edginess is only from the visible part
  
  vec3 k_a;
  k_a = vec3(0, 0, 0); // color of ambient (black)
  vec3 k_d;
  k_d = vec3(0.5, 0, 0); // color of diffuse
  float k_s;
  k_s = 0.8;

  out_color = vec4((k_a * illum) + (k_d*diffuse) + (k_s*shiny), 1);
  // Your awesome shader here!
  //out_color = (vec4(1, 1, 1, 0) + v_normal) / 2;
  //out_color.a = 1;
}
