#include "kdtree.h"
#include "defines.h"
#include "scene.h"
#include "scene_types.h"
#include <stdio.h>

#include <vector>
#include <stack>

typedef struct s_kdtreeNode KdTreeNode;
Scene *sceneGlobal = nullptr;
int axisGlobal = 0;

struct s_kdtreeNode {
  bool leaf; //! is this node a leaf ?
  int axis;//! axis index of the split, if not leaf
  float split;//!position of the split
  int depth; //!depth in the tree
  std::vector<int> objects;//! index of objects, if leaf
  KdTreeNode* left;//!ptr to left child
  KdTreeNode* right;//! ptr to right child
  vec3 min;//! min pos of node bounding box
  vec3 max;//! max pos of node bounding box
};

KdTreeNode * initNode(bool l, int a, int d) {
    KdTreeNode *ret = new KdTreeNode();
    ret->leaf = l;
    ret->axis = a;
    ret->depth = d;
    ret->left = NULL;
    ret->right = NULL;
    return ret;
}

typedef struct s_stackNode {
    float tmin;
    float tmax;
    KdTreeNode *node;
} StackNode;

struct s_kdtree {
    int depthLimit;
    size_t objLimit;
    KdTreeNode *root;

    std::vector<int> outOfTree;
    std::vector<int> inTree;
};

void subdivide(Scene *scene, KdTree *tree, KdTreeNode *node);

KdTree*  initKdTree(Scene *scene) {
  KdTree* tree = new KdTree();
  tree->root = initNode(false, 0, 0);
  tree->root->objects.clear();
  
  vec3 aabbmin = vec3(0);
  vec3 aabbmax = vec3(0);
  
  for (unsigned int i = 0; i < scene->objects.size(); i++) {
    Object *object = scene->objects.at(i);
    
    tree->root->objects.push_back(i);
    
    Geometry geom = object->geom;
    switch (object->geom.type) {
      case SPHERE:
	float minx = geom.sphere.center.x - geom.sphere.radius;
	float miny = geom.sphere.center.y - geom.sphere.radius;
	float minz = geom.sphere.center.z - geom.sphere.radius;
	
	if (aabbmin.x > minx) aabbmin.x = minx;
	if (aabbmin.y > miny) aabbmin.y = miny;
	if (aabbmin.z > minz) aabbmin.z = minz;
	
	float maxx = geom.sphere.center.x + geom.sphere.radius;
	float maxy = geom.sphere.center.y + geom.sphere.radius;
	float maxz = geom.sphere.center.z + geom.sphere.radius;
	
	if (aabbmax.x < maxx) aabbmax.x = maxx;
	if (aabbmax.y < maxy) aabbmax.y = maxy;
	if (aabbmax.z < maxz) aabbmax.z = maxz;
	break;
    }
  }
  
  tree->root->min = aabbmin;
  tree->root->max = aabbmax;
  tree->depthLimit = ceil(log(scene->objects.size()));
  
  subdivide(scene, tree, tree->root);
  
  return NULL;
}


//from http://blog.nuclex-games.com/tutorials/collision-detection/static-sphere-vs-aabb/
bool intersectSphereAabb(vec3 sphereCenter, float sphereRadius, vec3 aabbMin, vec3 aabbMax) {
    vec3 closestPointInAabb = min(max(sphereCenter, aabbMin), aabbMax);
    vec3 seg = closestPointInAabb -  sphereCenter;
    float distanceSquared = dot(seg, seg);
    // The AABB and the sphere overlap if the closest point within the rectangle is
    // within the sphere's radius
    return distanceSquared < (sphereRadius * sphereRadius);
}


// Generate children, compute split position, move objets to children and subdivide if needed.
void subdivide(Scene *scene, KdTree *tree, KdTreeNode *node) {
  
  if (tree->depthLimit == node->depth) {
    node->leaf = true;
  }
  
  if (node->objects.size() == 1) {
    node->leaf = true;
  }
  
  std::size_t size = scene->objects.size();
  sceneGlobal = scene;
  axisGlobal = node->axis;
  
  std::qsort (node->objects.data(), size, sizeof (int), [](const void* a, const void* b) {
    int arg1 = *static_cast<const int*>(a);
    int arg2 = *static_cast<const int*>(b);
    
    Object* o1 = sceneGlobal->objects.at(arg1);
    Object* o2 = sceneGlobal->objects.at(arg2);
    
    int v1 = 0, v2 = 0;
    
    switch (axisGlobal) {
      case 0: // X axis
	v1 = o1->geom.sphere.center.x;
	v2 = o2->geom.sphere.center.x;
	break;
      case 1: // Y axis
	v1 = o1->geom.sphere.center.y;
	v2 = o2->geom.sphere.center.y;
	break;
      case 2: // Z axis
	v1 = o1->geom.sphere.center.z;
	v2 = o2->geom.sphere.center.z;
	break;
    }
    
    if (v1 < v2) 
      return -1;
    if (v1 > v2)
      return 1;
    return 0;
  });
  
  int med = ceil(size / 2);
  Object *median = scene->objects.at(med);
  
  switch (axisGlobal) {
    case 0:
      node->split = median->geom.sphere.center.x;
      break;
    case 1:
      node->split = median->geom.sphere.center.y;
      break;
    case 2:
      node->split = median->geom.sphere.center.z;
      break;
  }
  
  node->left = initNode(false, ((node->axis + 1) % 3), node->depth+1);
  node->right = initNode(false, ((node->axis + 1) % 3), node->depth+1);
  
  for (size_t i = med; i < node->objects.size(); i++)
    node->right->objects.push_back(node->objects.at(i));
  
  node->left->objects = node->objects;
  node->left->objects.resize(med-1);
  
}

// Traverse kdtree to find intersection
bool traverse(Scene * scene, KdTree * tree, std::stack<StackNode> *stack, StackNode currentNode, Ray * ray, Intersection *intersection) {
  
  
  return false;
}



// from http://www.scratchapixel.com/lessons/3d-basic-lessons/lesson-7-intersecting-simple-shapes/ray-box-intersection/
bool intersectAabb(Ray *theRay,  vec3 min, vec3 max) {
    float tmin, tmax, tymin, tymax, tzmin, tzmax;
    vec3 bounds[2] = {min, max};
    tmin = (bounds[theRay->sign[0]].x - theRay->orig.x) * theRay->invdir.x;
    tmax = (bounds[1-theRay->sign[0]].x - theRay->orig.x) * theRay->invdir.x;
    tymin = (bounds[theRay->sign[1]].y - theRay->orig.y) * theRay->invdir.y;
    tymax = (bounds[1-theRay->sign[1]].y - theRay->orig.y) * theRay->invdir.y;
    if ((tmin > tymax) || (tymin > tmax)) return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;
    tzmin = (bounds[theRay->sign[2]].z - theRay->orig.z) * theRay->invdir.z;
    tzmax = (bounds[1-theRay->sign[2]].z - theRay->orig.z) * theRay->invdir.z;
    if ((tmin > tzmax) || (tzmin > tmax)) return false;
    if (tzmin > tmin) tmin = tzmin;
    if (tzmax < tmax) tmax = tzmax;
    if (tmin > theRay->tmin) theRay->tmin = tmin;
    if (tmax < theRay->tmax) theRay->tmax = tmax;
    return true;
}


bool intersectKdTree(Scene *scene, KdTree *tree, Ray *ray, Intersection *intersection) {
    bool hasIntersection = false;

    //!\todo call vanilla intersection on non kdtree object, then traverse the tree to compute other intersections

    return hasIntersection;
}
