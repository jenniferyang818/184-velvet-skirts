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

    double theta_offset = 2 * PI / (num_width_points);
    //  double beta_offset = 2 * PI / ((num_width_points - 1) / 2);

      // pleat variables
    //  double m = 0.1;
    //  double beta = 2 * PI;


    for (int i = 0; i < num_height_points; i++) {
        double theta = 0;
        double gamma = 0;
        for (int j = 0; j < num_width_points; j++) {
            bool is_pinned = false;
            if (i < 2 || j % 5 == 0) {
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
            curr_radius += radius_offset * exp(1/i);
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
            }
            else {
                int last = (i * num_width_points) + (num_width_points - 1);
                Spring spring = Spring(&point_masses[last], &point_masses[index], STRUCTURAL);
                springs.emplace_back(spring);
            }
            // structural -> above
            if (i > 0) {
                Spring spring = Spring(&point_masses[index - num_width_points], &point_masses[index], STRUCTURAL);
                springs.emplace_back(spring);
            }

            // shearing -> diagonal upper right
            if (i > 0 && j < num_width_points - 1) {
                int last = ((i - 1) * num_width_points) + (j + 1);
                Spring spring = Spring(&point_masses[last], &point_masses[index], SHEARING);
                springs.emplace_back(spring);
            }
            else if (i > 0 && j == num_width_points - 1) {
                int last = ((i - 1) * num_width_points);
                Spring spring = Spring(&point_masses[last], &point_masses[index], SHEARING);
                springs.emplace_back(spring);
            }

            // shearing -> diagonal upper left
            if (i > 0 && j > 0) {
                int last = ((i - 1) * num_width_points) + (j - 1);
                Spring spring = Spring(&point_masses[last], &point_masses[index], SHEARING);
                springs.emplace_back(spring);
            }
            else if (i > 0 && j == 0) {
                int last = ((i - 1) * num_width_points) + (num_width_points - 1);
                Spring spring = Spring(&point_masses[last], &point_masses[index], SHEARING);
                springs.emplace_back(spring);
            }
        }

    }
}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = width * height * cp->density / num_width_points / num_height_points;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  // TODO (Part 2): Compute total force acting on each point mass.
  Vector3D total_force;
  for (int i = 0; i < external_accelerations.size(); i++) {
      total_force += mass * external_accelerations[i];
  }
  for (int i = 0; i < num_width_points; i++) {
      for (int j = 0; j < num_height_points; j++) {
          PointMass* pm = &point_masses[i * num_height_points + j];
          pm->forces = total_force;
      }
  }
  for (int i = 0; i < springs.size(); i++) {
      Spring* spring = &springs[i];
      if (cp->enable_structural_constraints && spring->spring_type == STRUCTURAL) {
          Vector3D a_b_diff = spring->pm_a->position - spring->pm_b->position;
          double f_s = cp->ks * (a_b_diff.norm() - spring->rest_length);
          Vector3D force_vector = f_s * (spring->pm_b->position - spring->pm_a->position).unit();
          spring->pm_a->forces += force_vector;
          spring->pm_b->forces -= force_vector;
      }
      if (cp->enable_shearing_constraints && spring->spring_type == SHEARING) {
          Vector3D a_b_diff = spring->pm_a->position - spring->pm_b->position;
          double f_s = cp->ks * (a_b_diff.norm() - spring->rest_length);
          Vector3D force_vector = f_s * (spring->pm_b->position - spring->pm_a->position).unit();
          spring->pm_a->forces += force_vector;
          spring->pm_b->forces -= force_vector;
      }
      if (cp->enable_bending_constraints && spring->spring_type == BENDING) {
          Vector3D a_b_diff = spring->pm_a->position - spring->pm_b->position;
          double f_s = cp->ks * 0.2 * (a_b_diff.norm() - spring->rest_length);
          Vector3D force_vector = f_s * (spring->pm_b->position - spring->pm_a->position).unit();
          spring->pm_a->forces += force_vector;
          spring->pm_b->forces -= force_vector;
      }
  }
  
  // TODO (Part 2): Use Verlet integration to compute new point mass positions
  for (int i = 0; i < num_width_points; i++) {
      for (int j = 0; j < num_height_points; j++) {
          PointMass* pm = &point_masses[i * num_height_points + j];
          if (!(pm->pinned)) {
              Vector3D new_pos = pm->position + (1 - cp->damping / 100.0) * (pm->position - pm->last_position) + (pm->forces / mass) * pow(delta_t, 2);
              pm->last_position = pm->position;
              pm->position = new_pos;
          }
      }
  }



  // TODO (Part 4): Handle self-collisions.
  build_spatial_map();
    for (int i = 0; i < num_width_points; i++) {
        for (int j = 0; j < num_height_points; j++) {
            PointMass *curr_pm = &point_masses[i * num_height_points + j];
            this->self_collide(*curr_pm, simulation_steps);
            point_masses[i * num_height_points + j] = *curr_pm;
        }
    }
  // TODO (Part 3): Handle collisions with other primitives.
  for (int i = 0; i < num_width_points; i++) {
      for (int j = 0; j < num_height_points; j++) {
          PointMass *curr_pm = &point_masses[i * num_height_points + j];
          for (CollisionObject *collide : *collision_objects) {
                collide->collide(*curr_pm);
          }
      }
  }

  // TODO (Part 2): Constrain the changes to be such that the spring does not change
  // in length more than 10% per timestep [Provot 1995].
  for (int i = 0; i < springs.size(); i++) {
      Spring *s = &springs[i];
      Vector3D a_b_diff = s->pm_a->position - s->pm_b->position;
      double spring_length = a_b_diff.norm();
      if (spring_length > 1.1 * s->rest_length) {
          double remove_length = spring_length - 1.1 * s->rest_length;
//          double max_length = 1.1 * s->rest_length;
          if (!s->pm_a->pinned && !s->pm_b->pinned) {
              s->pm_a->position -= a_b_diff.unit() * remove_length / 2;
              s->pm_b->position += a_b_diff.unit() * remove_length / 2;
          }
          else if (!s->pm_a->pinned) {
              s->pm_a->position -= a_b_diff.unit() * remove_length;
          }
          else if (!s->pm_b->pinned) {
              s->pm_b->position += a_b_diff.unit() * remove_length;
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
  for (PointMass &pm : this->point_masses) {
      float hash = hash_position(pm.position);
      if (map.find(hash) == map.end()) {
          map[hash] = new vector<PointMass *>;
      }

//          std::vector<PointMass *> curr_list = std::vector<PointMass *>();
//          curr_list.emplace_back(&pm);
//          pair<float, vector<PointMass *> *> insert_pair = {hash, &curr_list};
//          map.insert(insert_pair);
//      } else {
//          std::vector<PointMass *> *curr_list = map[hash];
//          curr_list->emplace_back(&pm);
//          map[hash] = curr_list;
      map[hash]->emplace_back(&pm);
  }
}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
  // TODO (Part 4): Handle self-collision for a given point mass.
  //build_spatial_map();
    float curr_hash = hash_position(pm.position);
    vector<PointMass *> *candidates = map[curr_hash];
    Vector3D sum_correction;
    int count = 0;
    for (int i = 0; i < candidates->size(); i++) {
        PointMass *candidate = (*candidates)[i];
        if (candidate != &pm) {
            if ((candidate->position - pm.position).norm() < 2.0 * thickness) {
                Vector3D correction_vector = ((pm.position - candidate->position)).unit() * (2.0 * thickness - (candidate->position - pm.position).norm());
                sum_correction += correction_vector;
                count++;
            }
        }
    }
    if (count > 0) {
        sum_correction /= (double)count;
        sum_correction /= simulation_steps;
    }
    pm.position = sum_correction + pm.position;
}

float Cloth::hash_position(Vector3D pos) {
  // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.
    double w = 3.0 * width / num_width_points;
    double h = 3.0 * height / num_height_points;
    double t = max(w, h);
    Vector3D hashed_coords;
    if (orientation == HORIZONTAL) {
        int hash_to_z = int(pos.z / w);
        int hash_to_x = int(pos.x / h);
        int hash_to_y = int(pos.y / t);
        hashed_coords = Vector3D(hash_to_x, hash_to_y, hash_to_z);
    }
    else {
        int hash_to_z = int(pos.z / t);
        int hash_to_x = int(pos.x / h);
        int hash_to_y = int(pos.y / w);
        hashed_coords = Vector3D(hash_to_x, hash_to_y, hash_to_z);
    }
//    int hash_to_x = int(pos.x / w);
//    int hash_to_y = int(pos.y / h);
//    int hash_to_z = int(pos.z / t);
//    hashed_coords = Vector3D(hash_to_x, hash_to_y, hash_to_z);
    int star = max(num_width_points, num_height_points);
    int hash = (((hashed_coords.x * star) + hashed_coords.y) * star + hashed_coords.z);
    return float(hash);
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
      
      PointMass *pm_A = pm;
      PointMass* pm_B = pm + 1;
      /* & point_masses[y * num_width_points];
      if (x < num_width_points - 1) {
        pm_B = ;
      }*/
      PointMass *pm_C = pm + num_width_points    ;
      PointMass *pm_D = pm + num_width_points + 1;
      /*&point_masses[(y + 1) * num_width_points];
      if (x < num_width_points - 1) {
          pm_D = 
      }*/
      
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
