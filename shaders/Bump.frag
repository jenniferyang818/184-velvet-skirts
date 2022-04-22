#version 330

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

uniform vec4 u_color;

uniform sampler2D u_texture_2;
uniform vec2 u_texture_2_size;

uniform float u_normal_scaling;
uniform float u_height_scaling;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec2 v_uv;

out vec4 out_color;

float h(vec2 uv) {
  // You may want to use this helper function...
  float height;
  height = texture(u_texture_2, uv).r;
  return height;
}

void main() {
  // YOUR CODE HERE
  mat3 tbn;
        vec3 b;
        b = cross(v_normal.xyz, v_tangent.xyz);
        tbn[0] = normalize(v_tangent.xyz);
        tbn[1] = normalize(b);
        tbn[2] = normalize(v_normal.xyz);
        vec2 du_arg = vec2(v_uv[0] + (1 / u_texture_2_size[0]), v_uv[1]);
        float du = (h(du_arg) - h(v_uv)) * u_height_scaling * u_normal_scaling;
        vec2 dv_arg = vec2(v_uv[0], v_uv[1] + (1 / u_texture_2_size[1]));
        float dv = (h(dv_arg) - h(v_uv)) *  u_height_scaling * u_normal_scaling;

        vec3 n_o;
        n_o = vec3(-1 * du, -1 * dv, 1);
        vec3 n_d;
        n_d = tbn * n_o;

  // COPIED OVER FROM PHONG
   // AMBIENT
    // we pick this
    float ambient_coeff;
    ambient_coeff = 0.1;
    // we pick this
    vec3 ambient_illum;
    ambient_illum = vec3(1, 1, 1);

    vec4 ambient = vec4(ambient_coeff * ambient_illum, 1);

   // DIFFUSE
   // we pick this
   vec4 diffuse_coeff;
   diffuse_coeff = u_color;

    float distance_squared;
    distance_squared = pow(v_position.x - u_light_pos.x, 2) + pow(v_position.y - u_light_pos.y, 2) + pow(v_position.z - u_light_pos.z, 2);
    vec3 illum;
    illum = u_light_intensity / distance_squared;
    vec3 difference = vec3(u_light_pos.x - v_position.x, u_light_pos.y - v_position.y, u_light_pos.z - v_position.z);
    float max_thing;
    max_thing = max(0.0, dot(normalize(n_d), normalize(difference)));
    vec4 diffuse;
    diffuse = diffuse_coeff * vec4(illum * max_thing, 1);

      // SPECULAR
      // we pick this
      float specular_coeff;
      specular_coeff = 0.5;

      distance_squared = pow(v_position.x - u_light_pos.x, 2) + pow(v_position.y - u_light_pos.y, 2) + pow(v_position.z - u_light_pos.z, 2);
      illum = u_light_intensity / distance_squared;
      vec3 l_vector = vec3(u_light_pos.x - v_position.x, u_light_pos.y - v_position.y, u_light_pos.z - v_position.z);
      vec3 cam_vector = vec3(u_cam_pos.x - v_position.x, u_cam_pos.y - v_position.y, u_cam_pos.z - v_position.z);

      vec3 bisector = (l_vector + cam_vector) / length(l_vector + cam_vector);
      //we pick this
      float p;
      p = 100.0;
      max_thing = pow(max(0.0, dot(normalize(n_d), normalize(bisector))), p);
      vec4 specular;
      specular = vec4(specular_coeff * illum * max_thing, 1);

      out_color = ambient + diffuse + specular;
  out_color.a = 1;
}

