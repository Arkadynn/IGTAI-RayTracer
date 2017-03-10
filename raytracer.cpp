
#include "raytracer.h"
#include "scene_types.h"
#include "ray.h"
#include "image.h"
#include "kdtree.h"
#include <stdio.h>
#include <cmath>

#define MAX_DEPTH 5
#define ANTIALIASING_X 4

/// acne_eps is a small constant used to prevent acne when computing intersection
//  or boucing (add this amount to the position before casting a new ray !
const float acne_eps = 1e-4;
int cpt = 0;

bool intersectTriangle (Ray *ray, Intersection *intersection, Object *triangle) {
  point3 v0 = triangle->geom.triangle.v0;
  point3 v1 = triangle->geom.triangle.v1;
  point3 v2 = triangle->geom.triangle.v2;
  
  vec3 n = cross<float>((v1 - v0), (v2 - v0));
  n = normalize<float>(n);
  float cos_theta = dot<float>(n, ray->dir);
  
  if (cos_theta == 0)	return false;
  
  float D = dot<float>(-n, v0);

  float t = -(dot<float>(n, ray->orig) + D) / cos_theta;
  
  if (t < 0 || t > ray->tmax || t < ray->tmin) return false;
  
  vec3 hitPoint = rayAt(*ray, t);
  
  vec3 edge0 = v1 - v0;
  vec3 edge1 = v2 - v1;
  vec3 edge2 = v0 - v2;
  
  vec3 c0 = hitPoint - v0;
  vec3 c1 = hitPoint - v1;
  vec3 c2 = hitPoint - v2;
  
  
  if (dot<float>(n, cross<float>(edge0, c0)) > 0 &&
      dot<float>(n, cross<float>(edge1, c1)) > 0 &&
      dot<float>(n, cross<float>(edge2, c2)) > 0) {
    intersection->normal = n;
    intersection->position = hitPoint;
    intersection->mat = &triangle->mat;
    ray->tmax = t;
    return true;
  }
  
  return false;
}

bool intersectPlane(Ray *ray, Intersection *intersection, Object *obj) {
  bool hasIntersection = false;
  
  vec3 n = obj->geom.plane.normal;
  vec3 dir = ray->dir;
  
  float temp = dot<float>(n, dir);
  
  
  
  if ((hasIntersection = temp != 0)) {
    float d = obj->geom.plane.dist;
    point3 o = ray->orig;
    float numerator = dot<float>(o, n) + d;
    float denominator = dot<float>(n, dir);
    
    float t = (-numerator / denominator);
    if (t >= ray->tmin && ray->tmax >= t) {
      ray->tmax = t;
      intersection->mat = &obj->mat;
      intersection->normal = n;
      intersection->position = rayAt(*ray, t);
    } else {
      hasIntersection = false;
    }
  }
  
  return hasIntersection;

}

bool intersectEllipsoide(Ray *ray, Intersection *intersection, Object *obj) {
  bool hasIntersection = false;
  
  float t;
  vec3 d = ray->dir;
  point3 o = ray->orig;
  point3 centre_ = obj->geom.ellipsoide.center;
  float Ra = obj->geom.ellipsoide.a;
  float Rb = obj->geom.ellipsoide.b;
  float Rc = obj->geom.ellipsoide.c;
  
  vec3 v = vec3(Rb*Rb*Rc*Rc, Ra*Ra*Rc*Rc, Ra*Ra*Rb*Rb);
  
  float a = dot<float>(dot<float>(v, d), d);
  float b = 2 * dot<float>(dot<float>(v, d), o);
  float c = dot<float>(dot<float>(v, o), o) - a*a*b*b*c*c;
  
  float delta = b * b - 4.0f * a * c;
  
  
  if (delta >= 0) {
    if (delta == 0) {
      t = -b / (2 * a);
      hasIntersection = (t >= ray->tmin && t <= ray->tmax);
    } else {
      float res1 = (-b - sqrt(delta))/(2 * a);
      float res2 = (-b + sqrt(delta))/(2 * a);
      if (res1 >= 0 && res2 >= 0) {
	t = (res1 < res2) ? res1 : res2;
	hasIntersection = (t >= ray->tmin && t <= ray->tmax);
      } else if (res1 < 0 && res2 >= 0) {
	t = res2;
	hasIntersection = (t >= ray->tmin && t <= ray->tmax);
      } else if (res1 >= 0 && res2 < 0) {
	t = res1;
	hasIntersection = (t >= ray->tmin && t <= ray->tmax);
      } else {
	hasIntersection = false;
      }
    }
    
    if (hasIntersection) {
      ray->tmax = t;
      intersection->mat = &obj->mat;
      intersection->position = rayAt(*ray, t);
      vec3 n = v;
      intersection->normal = normalize<float>(n);
    }
  }
}

bool intersectSphere(Ray *ray, Intersection *intersection, Object *obj) {
  bool hasIntersection = false;
  
  // t^2 + 2t (d . (O - C)) - ((O - C) . (O - C) - R^2) = 0
  
  float t;
  vec3 d = ray->dir;
  point3 o = ray->orig;
  point3 centre_ = obj->geom.sphere.center;
  float r = obj->geom.sphere.radius;
  
  float a = 1;
  vec3 tmp = (o - centre_);
  float b = 2 * (dot<float>(d, tmp));
  float c = dot<float>(tmp, tmp) - r * r;
  
  float delta = b * b - 4.0f * a * c;
  
  if (delta >= 0) {
    //! \todo reprendre ici
    if (delta == 0) {
      t = -b / (2 * a);
      hasIntersection = (t >= ray->tmin && t <= ray->tmax);
    } else {
      float res1 = (-b - sqrt(delta))/(2 * a);
      float res2 = (-b + sqrt(delta))/(2 * a);
      if (res1 >= 0 && res2 >= 0) {
	t = (res1 < res2) ? res1 : res2;
	hasIntersection = (t >= ray->tmin && t <= ray->tmax);
      } else if (res1 < 0 && res2 >= 0) {
	t = res2;
	hasIntersection = (t >= ray->tmin && t <= ray->tmax);
      } else if (res1 >= 0 && res2 < 0) {
	t = res1;
	hasIntersection = (t >= ray->tmin && t <= ray->tmax);
      } else {
	hasIntersection = false;
      }
    }
    
    if (hasIntersection) {
      ray->tmax = t;
      intersection->mat = &obj->mat;
      intersection->position = rayAt(*ray, t);
      vec3 n = intersection->position - centre_;
      intersection->normal = normalize<float>(n);
    }
  }
  
  return hasIntersection;
}

bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection) {
  bool hasIntersection = false;

  for (Object *o : scene->objects) {
    switch (o->geom.type) {
      case SPHERE:
	hasIntersection |= intersectSphere(ray, intersection, o);
	break;
      case PLANE:
	hasIntersection |= intersectPlane(ray, intersection, o);
	break;
      case TRIANGLE:
	hasIntersection |= intersectTriangle(ray, intersection, o);
	break;
      default:
	perror("An unhandeld object have been found\n");
    }
  }

  return hasIntersection;
}

/* --------------------------------------------------------------------------- */
/*
 *	The following functions are coded from Cook-Torrance bsdf model description and are suitable only
 *  for rough dielectrics material (RDM. Code has been validated with Mitsuba renderer)
 */

/** Normal Distribution Function : Beckmann
 * NdotH : Norm . Half
 */
float RDM_Beckmann(float NdotH, float alpha) {
    
  float cos2_theta = NdotH * NdotH;
  float cos4_theta = cos2_theta * cos2_theta;
  float tan2_theta = (1-cos2_theta) / cos2_theta;
  
  float alpha2 = alpha * alpha;
  
  float tmp = exp(-tan2_theta/(alpha2));
  
  float ret = tmp / (M_PI * alpha2 * cos4_theta);
  
  return ret;

}

// Fresnel term computation. Implantation of the exact computation. we can use the Schlick approximation
// LdotH : Light . Half
float RDM_Fresnel(float LdotH, float extIOR, float intIOR) {
  float tmp;
  
  float cos_theta_i = LdotH;
  float n1 = extIOR;
  float n2 = intIOR;
  
  float cos2_theta_i = cos_theta_i * cos_theta_i;
  
  tmp = n1 / n2;
  tmp *= tmp;
  float sin2_theta_t = tmp * (1 - cos2_theta_i);
  
  if (sin2_theta_t > 1)
    return 1;
  
  float cos_theta_t = sqrt(1 - sin2_theta_t);
  
  tmp = (n1 * cos_theta_i - n2 * cos_theta_t);
  float rs = tmp * tmp;
  tmp = (n1 * cos_theta_i + n2 * cos_theta_t);
  rs /= (tmp*tmp);
  
  tmp = (n1 * cos_theta_t - n2 * cos_theta_i);
  float rp = tmp * tmp;
  tmp = (n1 * cos_theta_t + n2 * cos_theta_i);
  rp /= (tmp*tmp);
  
  float ret = (rs + rp) / 2;
  
  return ret;
}


// Shadowing and masking function. Linked with the NDF. Here, Smith function, suitable for Beckmann NDF
float RDM_chiplus(float c) {
  return (c > 0.f) ? 1.f : 0.f;
}

// DdotH : Dir . Half
// HdotN : Half . Norm
float RDM_G1(float DdotH, float DdotN, float alpha) {
  float ret = 1;
  
  float cos_theta = DdotN;
  float cos2_theta = cos_theta * cos_theta;
  
  float tmp = sqrt(1 - cos2_theta);
  float tan_theta = tmp / cos_theta;
  
  float b = 1 / (alpha * tan_theta);
  float k = DdotH / DdotN;
  
  
  if (k > 0.0f && b < 1.6f) {
    ret = 3.535f * b + 2.181f * b * b;
    ret /= (1.0f + 2.276f * b + 2.577f * b * b);
  } else {
    ret = 1;
  }
  
  ret *= RDM_chiplus(k);
  
  return ret;
}

// LdotH : Light . Half
// LdotN : Light . Norm
// VdotH : View . Half
// VdotN : View . Norm
float RDM_Smith(float LdotH, float LdotN, float VdotH, float VdotN, float alpha) {
  return RDM_G1(LdotH, LdotN, alpha) * RDM_G1(VdotH, VdotN, alpha);
  
}

// Specular term of the Cook-torrance bsdf
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdotN : View . Norm
color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN, Material *m) {
  float D = RDM_Beckmann(NdotH, m->roughness);
  float F = RDM_Fresnel(LdotH, 1, m->IOR);
  float G = RDM_Smith(LdotH, LdotN, VdotH, VdotN, m->roughness);
  
  color3 tmp = m->specularColor * D * F * G;
  tmp /= (4 * LdotN * VdotN);
  return tmp;

  
}
// diffuse term of the cook torrance bsdf
color3 RDM_bsdf_d(Material *m) {
  return m->diffuseColor / (float)(M_PI);

}

// The full evaluation of bsdf(wi, wo) * cos (thetai)
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdtoN : View . Norm
// compute bsdf * cos(Oi)
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN, Material *m) {
  
  return RDM_bsdf_d(m) + RDM_bsdf_s(LdotH, NdotH, VdotH, LdotN, VdotN, m);
}




/* --------------------------------------------------------------------------- */

color3 shade(vec3 n, vec3 v, vec3 l, color3 lc, Material *mat ){
  color3 ret = color3(0.f);
  
  float cos_theta = dot<float>(n, l);
  if (cos_theta < 0)
    return ret;
  
  vec3 h = (v + l) / (length<float>(v + l));
  h = normalize<float>(h);
  
  float LdotH = dot<float>(l, h);
  float NdotH = dot<float>(n, h);
  float VdotH = dot<float>(v, h);
  float LdotN = dot<float>(l, n);
  float VdotN = dot<float>(v, n);
  
  return RDM_bsdf(LdotH, NdotH, VdotH, LdotN, VdotN, mat) * lc * cos_theta;
	    
}

//! if tree is not null, use intersectKdTree to compute the intersection instead of intersect scene
color3 trace_ray(Scene * scene, Ray *ray, KdTree *tree) {  
  color3 ret = color3(0.f, 0.f, 0.f);
  
  if (ray->depth > MAX_DEPTH) return ret;
  
  Intersection intersection;
  
  if (intersectScene(scene, ray, &intersection)) {
    for (Light *light : scene->lights) {
      vec3 light_dir = light->position - intersection.position;
      vec3 l = normalize<float>(light_dir);
      Ray r;
      rayInit(&r, intersection.position, l, acne_eps, length<float>(light_dir));
      Intersection shadow;
      if (!intersectScene(scene, &r, &shadow)) {
	ret += shade(intersection.normal, -ray->dir, l, light->color, intersection.mat);
      }	
    }
    
    vec3 newDir = normalize<float>(reflect(ray->dir, intersection.normal));
    float LdotH = dot<float>(newDir, normalize<float>(ray->dir + newDir));
    rayInit(ray, intersection.position, newDir, acne_eps, 100000, ray->depth+1);
    ret += RDM_Fresnel(LdotH, 1, intersection.mat->IOR) * trace_ray(scene, ray, tree);
    
  } else {
    ret = scene->skyColor;
  }
  
  return ret;
}

void renderImage(Image *img, Scene *scene) {

  //! This function is already operational, you might modify it for antialiasing and kdtree initializaion
  float aspect = 1.f/scene->cam.aspect;
    
  KdTree *tree =  NULL;

  tree = initKdTree(scene);
  //! \todo initialize KdTree

  float delta_y = 1.f / (img->height * 0.5f); //! one pixel size
  vec3 dy = delta_y * aspect * scene->cam.ydir; //! one pixel step 
  vec3 ray_delta_y = (0.5f - img->height * 0.5f) / (img->height * 0.5f) * aspect * scene->cam.ydir;

  float delta_x = 1.f / (img->width * 0.5f);
  vec3 dx = delta_x * scene->cam.xdir;
  vec3 ray_delta_x = (0.5f - img->width * 0.5f) / (img->width * 0.5f) *scene->cam.xdir;
  
    
  for(size_t j=0; j<img->height; j++) {
    if(j!=0) printf("\033[A\r");
    float progress = (float)j/img->height*100.f;
    printf("progress\t[");
    int cpt = 0;
    for(cpt = 0; cpt<progress; cpt+=5) printf(".");
    for(       ; cpt<100; cpt+=5) printf(" ");
    printf("]\n");
#pragma omp parallel for
    for(size_t i=0; i<img->width; i++) {
      color3 *ptr = getPixelPtr(img, i,j);
      vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y + float(i)*dx + float(j)*dy;

      Ray rx;
      rayInit(&rx, scene->cam.position, normalize(ray_dir));
      *ptr = trace_ray(scene, &rx, tree);

    }
  }
}
