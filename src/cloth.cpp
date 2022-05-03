#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
             int num_height_points, float thickness) {
  this->width = width;
  this->height = height;
  this->num_width_points = num_width_points;
  this->num_height_points = num_height_points;
  this->thickness = thickness;

  buildGrid();
  buildClothMesh();
}

Cloth::~Cloth() {
  point_masses.clear();
  springs.clear();

  if (clothMesh) {
    delete clothMesh;
  }
}

void Cloth::buildGrid() {
  // let num_width_points and num_height_points be the max skirt height and width
  // make extra datastructure to build springs
  double height_offset = height / num_height_points;
  double curr_height = height;
  
  double max_radius = width;
  double min_radius = width / 4;
  double curr_radius = min_radius;
  double radius_offset = (max_radius - min_radius) / num_width_points;
  
  double theta_offset = 2 * PI / (num_width_points - 1);
//  double beta_offset = 2 * PI / ((num_width_points - 1) / 2);
  
  // pleat variables
//  double m = 0.1;
//  double beta = 2 * PI;
  
  
  for (int i = 0; i < num_height_points; i++) {
    double theta = 0;
    double gamma = 0;
    for (int j = 0; j < num_width_points; j++) {
      bool is_pinned = false;
      if (i < 2 || j % 7 == 0) {
        is_pinned = true;
      }
      
      double x = curr_radius * cos(theta);
      double y = curr_height;
      double z = curr_radius * sin(theta);

      Vector3D pos = Vector3D(x, y, z);
      PointMass mass = PointMass(pos, is_pinned);
      point_masses.emplace_back(mass);
      
      theta += theta_offset;
    }
    if (i >= 2) {
      curr_radius += radius_offset;
    }
    curr_height -= height_offset;
  }
  
//  bezierCurve(0.25);
  
//  springs.clear();
  for (int i = 0; i < num_height_points; i++) {
    for (int j = 0; j < num_width_points; j++) {
      int index = (i * num_width_points) + j;
      // structural -> left
      if (j > 0) {
        Spring spring = Spring(&point_masses[index - 1], &point_masses[index], STRUCTURAL);
        springs.emplace_back(spring);
      } else {
        int last = (i * num_width_points) + (num_width_points - 1);
        Spring spring = Spring(&point_masses[last], &point_masses[index], STRUCTURAL);
        springs.emplace_back(spring);
      }
      // structural -> above
      if (i > 0) {
        Spring spring = Spring(&point_masses[i - 1], &point_masses[i], STRUCTURAL);
        springs.emplace_back(spring);
      }
      
      // shearing -> diagonal upper right
      if (i > 0 && j < num_width_points - 1) {
        int last = ((i - 1) * num_width_points) + (j + 1);
        Spring spring = Spring(&point_masses[last], &point_masses[index], SHEARING);
        springs.emplace_back(spring);
      } else if (i > 0 && j == num_width_points - 1) {
        int last = ((i - 1) * num_width_points);
        Spring spring = Spring(&point_masses[last], &point_masses[index], SHEARING);
        springs.emplace_back(spring);
      }
      
      // shearing -> diagonal upper left
      if (i > 0 && j > 0) {
        int last = ((i - 1) * num_width_points) + (j - 1);
        Spring spring = Spring(&point_masses[last], &point_masses[index], SHEARING);
        springs.emplace_back(spring);
      } else if (i > 0 && j == 0) {
        int last = ((i - 1) * num_width_points) + (num_width_points - 1);
        Spring spring = Spring(&point_masses[last], &point_masses[index], SHEARING);
        springs.emplace_back(spring);
      }
    }
    
  }
}

//void Cloth::bezierCurve(double t) {
//  std::vector<Vector3D> new_pos = std::vector<Vector3D>();
//  new_pos.reserve(point_masses.size() - 6);
//  for (int i = 0; i < num_width_points; i++) {
//    // start after 2 -> after hemming
//    // start after another 4 bc of control points
//    // j = 6
//    for (int j = 6; j < num_height_points; j++) {
//      int index = (j * num_width_points) + i;
//      int p0 = ((j - 4) * num_width_points) + i;
//      int p1 = ((j - 3) * num_width_points) + i;
//      int p2 = ((j - 2) * num_width_points) + i;
//      int p3 = ((j - 1) * num_width_points) + i;
//      Vector3D pos = pow((1 - t), 3) * point_masses[p0].position;
//      pos += 3 * t * pow(1 - t, 2) * point_masses[p1].position;
//      pos += 3 * pow(t, 2) * (1 - t) * point_masses[p2].position;
//      pos += pow(t, 3) * point_masses[p3].position;
//      new_pos[((j - 6) * num_width_points) + i] = pos;
//    }
//  }
//  for (int i = 0; i < num_width_points; i++) {
//    for (int j = 6; j < num_height_points; j++) {
//      int index = (j * num_width_points) + i;
//      point_masses[index].position = new_pos[((j - 6) * num_width_points) + i];
//    }
//  }
//}

void Cloth::bezierCurve(double t) {
//  std::vector<Vector3D> new_pos = std::vector<Vector3D>();
//  new_pos.reserve(point_masses.size() - 6);
  for (int i = 0; i < num_width_points; i++) {
    // start after 2 -> after hemming
    // start after another 4 bc of control points
    // j = 6
    for (int j = 6; j < num_height_points; j++) {
      int index = (j * num_width_points) + i;
      int p0 = ((j - 4) * num_width_points) + i;
      int p1 = ((j - 3) * num_width_points) + i;
      int p2 = ((j - 2) * num_width_points) + i;
      int p3 = ((j - 1) * num_width_points) + i;
      Vector3D pos = pow((1 - t), 3) * point_masses[p0].position;
      pos += 3 * t * pow(1 - t, 2) * point_masses[p1].position;
      pos += 3 * pow(t, 2) * (1 - t) * point_masses[p2].position;
      pos += pow(t, 3) * point_masses[p3].position;
    }
  }
}


void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
    double mass = width * height * cp->density / num_width_points / num_height_points;
    double delta_t = 1.0f / frames_per_sec / simulation_steps;

    // TODO (Part 2): Compute total force acting on each point mass.
    Vector3D ext_force = Vector3D();
    for (int i = 0; i < external_accelerations.size(); i++) {
      ext_force += (mass * external_accelerations[i]);
    }
    for (int j = 0; j < point_masses.size(); j++) {
      point_masses[j].forces = ext_force;
    }

    for (int i = 0; i < springs.size(); i++) {
      PointMass a = *springs[i].pm_a;
      PointMass b = *springs[i].pm_b;
      if ((cp->enable_structural_constraints && springs[i].spring_type == STRUCTURAL) ||
          (cp->enable_shearing_constraints && springs[i].spring_type == SHEARING) ||
          (cp->enable_bending_constraints && springs[i].spring_type == BENDING)) {
        Vector3D delta_pos = b.position - a.position;
        double mag = delta_pos.norm();
        double spring_const = cp->ks;
        if (springs[i].spring_type == BENDING) {
//          spring_const *= 1 - 0.2;
          spring_const *= 0.2;
        }

        double fs = spring_const * (mag - springs[i].rest_length);
        Vector3D f_ab =  fs * delta_pos.unit();
        springs[i].pm_a->forces += f_ab;
        springs[i].pm_b->forces -= f_ab;
      }
    }

  // TODO (Part 2.2): Use Verlet integration to compute new point mass positions
  for (int i = 0; i < point_masses.size(); i++) {
    if (!point_masses[i].pinned) {
      Vector3D new_pos = point_masses[i].position + (1 - (cp->damping / 100.0)) * (point_masses[i].position - point_masses[i].last_position) + (point_masses[i].forces / mass) * pow(delta_t, 2);
      point_masses[i].last_position = point_masses[i].position;
      point_masses[i].position = new_pos;
    }
  }

  build_spatial_map();
//   TODO (Part 4): Handle self-collisions.
  for (int j = 0; j < point_masses.size(); j++) {
    self_collide(point_masses[j], simulation_steps);
    for (CollisionObject *c : *collision_objects) {
      c->collide(point_masses[j]);
    }
  }

  // TODO (Part 3): Handle collisions with other primitives.

  // TODO (Part 2.3): Constrain the changes to be such that the spring does not change
    for (Spring s : this->springs) {
        PointMass a = *s.pm_a;
        PointMass b = *s.pm_b;
        double length = (b.position - a.position).norm();
        Vector3D direction_vector = (b.position - a.position).unit();

        if (length > s.rest_length * 1.1) {
            double change_distance = (s.rest_length / length * 1.1);
            if (!a.pinned && !b.pinned) {
                a.position -=  direction_vector * change_distance * 0.05;
                b.position += direction_vector * change_distance * 0.05;
            } else if (!a.pinned && b.pinned) {
                a.position -= direction_vector * change_distance;
            } else if (!b.pinned && a.pinned) {
                b.position += direction_vector * change_distance;
            }
        }
    }
}

void Cloth::build_spatial_map() {
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  // TODO (Part 4): Build a spatial map out of all of the point masses.
  for (int i = 0; i < point_masses.size(); i ++) {
    float hash = hash_position(point_masses[i].position);
    if (map.count(hash)) {
      map[hash]->push_back(&point_masses[i]);
    } else {
      vector<PointMass*>* bucket = new vector<PointMass*>({&point_masses[i]});
      map[hash] = bucket;
    }
  }
}


void Cloth::self_collide(PointMass &pm, double simulation_steps) {
    // TODO (Part 4): Handle self-collision for a given point mass.
  float hash = hash_position(pm.position);
  if (map.count(hash)) {
//    cout << ":)";
    vector<PointMass*>* bucket = map[hash];
    // get candidates for collision
    Vector3D correction_vector = Vector3D();
    int count = 0;
//    cout << "no";
    for(PointMass *p : *bucket) {
//      cout << "yes";
      // determine if pm and candidate are w/in 2*thickness apart
      Vector3D dist = pm.position - p->position;
      bool eq = (p == &pm);
      if (!eq && dist.norm() <= 2 * thickness) {
        // if so, compute correction vector
        count += 1;
//        pm.position -= p->position - pm.position;
        correction_vector += dist.unit() * (2 * thickness - dist.norm());
      }
    }
    // final correction vector = avg(correction vectors, scaled down by simu steps)
    if (count > 0) {
      pm.position += correction_vector / count / simulation_steps;
    }
  }
}

float Cloth::hash_position(Vector3D pos) {
  // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.
  float w = 3 * width / num_width_points;
  float h = 3 * height / num_height_points;
  float t = max(w,h);
  int wi = floor(pos.x / w);
  int hi = floor(pos.y / h);
  int ti = floor(pos.z / t);
  
  float unique_pos = float(ti * (w * h)) + ((hi * w) + wi);
  return unique_pos;
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
  PointMass *pm = &point_masses[0];
  for (int i = 0; i < point_masses.size(); i++) {
    pm->position = pm->start_position;
    pm->last_position = pm->start_position;
    pm++;
  }
}

void Cloth::buildClothMesh() {
  if (point_masses.size() == 0) return;

  ClothMesh *clothMesh = new ClothMesh();
  vector<Triangle *> triangles;

  // Create vector of triangles
  for (int y = 0; y < num_height_points - 1; y++) {
    for (int x = 0; x < num_width_points - 1; x++) {
      PointMass *pm = &point_masses[y * num_width_points + x];
      // Get neighboring point masses:
      /*                      *
       * pm_A -------- pm_B   *
       *             /        *
       *  |         /   |     *
       *  |        /    |     *
       *  |       /     |     *
       *  |      /      |     *
       *  |     /       |     *
       *  |    /        |     *
       *      /               *
       * pm_C -------- pm_D   *
       *                      *
       */
      
      float u_min = x;
      u_min /= num_width_points - 1;
      float u_max = x + 1;
      u_max /= num_width_points - 1;
      float v_min = y;
      v_min /= num_height_points - 1;
      float v_max = y + 1;
      v_max /= num_height_points - 1;
      
      PointMass *pm_A = pm                       ;
      PointMass *pm_B = pm                    + 1;
      PointMass *pm_C = pm + num_width_points    ;
      PointMass *pm_D = pm + num_width_points + 1;
      
      Vector3D uv_A = Vector3D(u_min, v_min, 0);
      Vector3D uv_B = Vector3D(u_max, v_min, 0);
      Vector3D uv_C = Vector3D(u_min, v_max, 0);
      Vector3D uv_D = Vector3D(u_max, v_max, 0);
      
      
      // Both triangles defined by vertices in counter-clockwise orientation
      triangles.push_back(new Triangle(pm_A, pm_C, pm_B,
                                       uv_A, uv_C, uv_B));
      triangles.push_back(new Triangle(pm_B, pm_C, pm_D,
                                       uv_B, uv_C, uv_D));
    }
  }

  // For each triangle in row-order, create 3 edges and 3 internal halfedges
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    // Allocate new halfedges on heap
    Halfedge *h1 = new Halfedge();
    Halfedge *h2 = new Halfedge();
    Halfedge *h3 = new Halfedge();

    // Allocate new edges on heap
    Edge *e1 = new Edge();
    Edge *e2 = new Edge();
    Edge *e3 = new Edge();

    // Assign a halfedge pointer to the triangle
    t->halfedge = h1;

    // Assign halfedge pointers to point masses
    t->pm1->halfedge = h1;
    t->pm2->halfedge = h2;
    t->pm3->halfedge = h3;

    // Update all halfedge pointers
    h1->edge = e1;
    h1->next = h2;
    h1->pm = t->pm1;
    h1->triangle = t;

    h2->edge = e2;
    h2->next = h3;
    h2->pm = t->pm2;
    h2->triangle = t;

    h3->edge = e3;
    h3->next = h1;
    h3->pm = t->pm3;
    h3->triangle = t;
  }

  // Go back through the cloth mesh and link triangles together using halfedge
  // twin pointers

  // Convenient variables for math
  int num_height_tris = (num_height_points - 1) * 2;
  int num_width_tris = (num_width_points - 1) * 2;

  bool topLeft = true;
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    if (topLeft) {
      // Get left triangle, if it exists
      if (i % num_width_tris != 0) { // Not a left-most triangle
        Triangle *temp = triangles[i - 1];
        t->pm1->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm1->halfedge->twin = nullptr;
      }

      // Get triangle above, if it exists
      if (i >= num_width_tris) { // Not a top-most triangle
        Triangle *temp = triangles[i - num_width_tris + 1];
        t->pm3->halfedge->twin = temp->pm2->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle to bottom right; guaranteed to exist
      Triangle *temp = triangles[i + 1];
      t->pm2->halfedge->twin = temp->pm1->halfedge;
    } else {
      // Get right triangle, if it exists
      if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
        Triangle *temp = triangles[i + 1];
        t->pm3->halfedge->twin = temp->pm1->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle below, if it exists
      if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
        Triangle *temp = triangles[i + num_width_tris - 1];
        t->pm2->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm2->halfedge->twin = nullptr;
      }

      // Get triangle to top left; guaranteed to exist
      Triangle *temp = triangles[i - 1];
      t->pm1->halfedge->twin = temp->pm2->halfedge;
    }

    topLeft = !topLeft;
  }

  clothMesh->triangles = triangles;
  this->clothMesh = clothMesh;
}
