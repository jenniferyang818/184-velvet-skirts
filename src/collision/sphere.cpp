#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../misc/sphere_drawing.h"
#include "sphere.h"

using namespace nanogui;
using namespace CGL;

void Sphere::collide(PointMass &pm) {
  // TODO (Part 3): Handle collisions with spheres.
  // calculate distance from point mass to center of sphere
  Vector3D difference = pm.position - this->origin; // if we hit a bug try changing this to last position
  double distance = difference.norm();
  if (distance <= this->radius) {
      Vector3D tangent_point = radius * difference.unit() + this->origin;
      Vector3D correction_vector = tangent_point - pm.last_position;
      Vector3D new_position = pm.last_position + correction_vector * (1 - friction);
      pm.position = new_position;
  }
}

void Sphere::render(GLShader &shader) {
  // We decrease the radius here so flat triangles don't behave strangely
  // and intersect with the sphere when rendered
  m_sphere_mesh.draw_sphere(shader, origin, radius * 0.92);
}
