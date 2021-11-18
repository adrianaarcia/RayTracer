#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <cfloat>
#include <math.h>

#include "SETTINGS.h"

using namespace std;

typedef struct _SPHERE
{
  float rad;
  VEC3 center;
  VEC3 color; 
} SPHERE;

typedef struct _LIGHT
{
  VEC3 position;
  VEC3 color; 
} LIGHT;

///////////////////////////////////////////////////////////////////////
// Given Functions
///////////////////////////////////////////////////////////////////////
void readPPM(const string& filename, int& xRes, int& yRes, float*& values);
void writePPM(const string& filename, int& xRes, int& yRes, const float* values);

///////////////////////////////////////////////////////////////////////
// My Functions
///////////////////////////////////////////////////////////////////////
void draw(vector<SPHERE>& spheres, vector<VEC3>& rays, vector<LIGHT>& lights, VEC3 ambient, float phong,
const VEC3 origin, float* values, float res, int mode);
void specular_hl(SPHERE object, VEC3 ray, vector<LIGHT> lights, VEC3 ambient, float phong,
const VEC3 origin, float t, int i, float* values, float res);
void multiple_lights(SPHERE object, VEC3 ray, vector<LIGHT> lights, const VEC3 origin, 
float t, int i, float* values, float res);
void diffuse_single_light(SPHERE object, VEC3 ray, LIGHT& light, const VEC3 origin, 
float t, int i, float* values, float res);
float intersect(const VEC3 ray, const SPHERE object, const VEC3 origin, float eps);
void set_pixel_color(float* values, const VEC3 color, int pixel);
void scale_vecs(vector<VEC3>& vecs, float scalar);
void show_ray_comps(vector<VEC3> rays, float* values, int xRes, int yRes, int size);
void generate_rays(vector<VEC3>& rays, VEC3 eye, VEC3 w, VEC3 u, VEC3 v, float l, float r, float t, float b, float nx, float ny, 
int xRes, int yRes);
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void shadows(SPHERE object, VEC3 ray, vector<LIGHT> lights, VEC3 ambient, float phong,
const VEC3 origin, float t, int i, float* values, float res)
{
  LIGHT light;
  VEC3 color, pt, v, n, l, h; color.setZero();
  pt = origin + t*ray;
  int light_count = lights.size();
  
  for(int j = 0; j < light_count; j++)
  {
    light = lights[j];
    v = origin - pt; v.normalize();
    n = pt - object.center; n.normalize();
    l = light.position - pt; //l.normalize();
    
    if(intersect(l, object, pt, 0.01) != NULL)
    {
      h = v+l; h.normalize();
      l.normalize();
      float d_val = fmax(0, n.dot(l));
      float s_val = pow(fmax(0, n.dot(h)),phong);
    
      for(int k = 0; k < 3; k++)
      {
        color[k] += object.color[k]*(ambient[k] + light.color[k]*d_val + light.color[k]*s_val);
      }
    }
    else
    {
      for(int k = 0; k < 3; k++)
      {
        color[k] += object.color[k]*(ambient[k]);
      }
    }
    
  } 
  set_pixel_color(values, 255*color, i);
}

int main(int argc, char** argv)
{
  /* RESOLUTION */
  float nx = 800.0; float ny = 600.0; int xRes = (int) nx; int yRes = (int) ny;
  int res = xRes*yRes;
  int size = res*3;
  
  /* CAMERA PARAMETERS */
  float near = 1.0;
  float fovy = 65.0; //y direction, top and bottom of the image (in degrees)

  VEC3 eye(0,0,0); VEC3 lookAt(0,0,1); VEC3 up(0,1,0);
  
  float theta = fovy * M_PI/180; //fovy in radians
  float aspect = nx/ny;
  float l,r,t,b; //planes;
  t = near*tanf(theta/2); b = -t;
  r = t*aspect; l = -r;


  //calculate basis vectors
  VEC3 w; w.setZero(); VEC3 u; u.setZero(); VEC3 v; v.setZero();
  VEC3 g; g = eye - lookAt;
  w = g.normalized();
  u = up.cross(w); u.normalize();
  v = w.cross(u);

  /* SCENE GEOMETRY */
  float phong = 10.0;

  /* Spheres (center, color) */
  vector<SPHERE> spheres;
  //sphere 0 -- red, left 
  SPHERE sph0; sph0.rad = 3; 
  sph0.center = VEC3(-3.5,0,10); sph0.color = VEC3(1,0.25,0.25);
  spheres.push_back(sph0);
 
  //sphere 1 -- blue, right
  SPHERE sph1; sph1.rad = 3;
  sph1.center = VEC3(3.5,0,10); sph1.color = VEC3(0.25,0.25,1);
  spheres.push_back(sph1);
  
  //sphere 2 -- gray, floor
  SPHERE sph2; sph2.rad = 997;
  sph2.center = VEC3(0,-1000,10); sph2.color = VEC3(0.5,0.5,0.5);
  spheres.push_back(sph2);
 
  /* Lights (position, color) */
  vector<LIGHT> lights;
  
  LIGHT light0; light0.position = VEC3(10,3,5); light0.color = VEC3(1,1,1); 
  lights.push_back(light0);
  
  LIGHT light1; light1.position = VEC3(-10,3,7.5); light1.color = VEC3(0.5,0,0);
  lights.push_back(light1);

  VEC3 ambient(0,0,0);

  //initalize array of values to zero
  float* values = (float*) calloc(size, sizeof(float));

  /*****          1 Ray Tracing, part I          *****/

  /* 1.1 Eye Ray Generation */
  vector<VEC3> rays;
  generate_rays(rays, eye, w, u, v, l, r, t, b, nx, ny, xRes, yRes);
  show_ray_comps(rays, values, xRes, yRes, size);

  /* 1.2 Next Features */
  //Scene Intersection
  draw(spheres, rays, lights, ambient, phong, eye, values, res, 2); 
  writePPM("2.png", xRes, yRes, values); memset(values,0,sizeof(float)*res*3);
  
  //Shading
  draw(spheres, rays, lights, ambient, phong, eye, values, res, 3); 
  writePPM("3.png", xRes, yRes, values); memset(values,0,sizeof(float)*res*3);

  //Multiple Lights;
  draw(spheres, rays, lights, ambient, phong,  eye, values, res, 4); 
  writePPM("4.png", xRes, yRes, values); memset(values,0,sizeof(float)*res*3);
  
  //Specular Highlights
  draw(spheres, rays, lights, ambient, phong, eye, values, res, 5); 
  writePPM("5.png", xRes, yRes, values); memset(values,0,sizeof(float)*res*3);
  
  //Shadows
  draw(spheres, rays, lights, ambient, phong, eye, values, res, 6); 
  writePPM("6.png", xRes, yRes, values);

  /*****          2 Ray Tracing, part 2          *****/
  light0.position = VEC3(10,10,5); light0.color = VEC3(1,1,1); 
  light1.position = VEC3(-10,10,7.5); light1.color = VEC3(0.5,0.25,0.25);
  
  return 0;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
float intersect(const VEC3 ray, const SPHERE object, const VEC3 origin, float eps)
{
  VEC3 temp(0,0,0);
  float rad, t, t1, t2, A, B, C, discriminant; 
  t = t1 = t2 = A = B = C = discriminant = 0; 
  
  rad = object.rad;
  temp = origin - object.center;

  //use the quadratic formula
  A = ray.dot(ray);
  B = 2*ray.dot(temp);
  C = temp.dot(temp) - (rad*rad);
  discriminant = (B*B)-(4*A*C);
  if(discriminant >= 0)
  {
    discriminant = sqrt(discriminant);
    t1 = fmax(eps, (-B + discriminant)/(2*A));
    t2 = fmax(eps, (-B - discriminant)/(2*A));

    t = fmin(t1,t2);

    return t;
  }
  else return NULL;
}

void draw(vector<SPHERE>& spheres, vector<VEC3>& rays, vector<LIGHT>& lights, VEC3 ambient, float phong, const VEC3 origin, 
float* values, float res, int mode)
{
  SPHERE object, closest_obj; VEC3 ray;
  int sphere_count = spheres.size();
  float t, tmin; t = tmin = INFINITY;
  
  for(int i = 0; i < res; i++) //for each pixel
  {
    ray = rays[i]; 
    tmin = INFINITY;
    for(int j = 0; j < sphere_count; j++) //for each object/sphere
    {
      object = spheres[j];
      if((t = intersect(ray, object, origin, 0)))
      {
        if(t < tmin) 
        {
          closest_obj = object;
          tmin = t;

          //recolor the pixel
          if (mode == 2) set_pixel_color(values, 255*closest_obj.color, i); //simple scene intersection
          if (mode == 3) diffuse_single_light(closest_obj, ray, lights[0], origin, t, i, values, res);
          if (mode == 4) multiple_lights(closest_obj, ray, lights, origin, t, i, values, res);
          if (mode == 5) specular_hl(closest_obj, ray, lights, ambient, phong, origin, t, i, values, res);
          if (mode == 6) shadows(closest_obj, ray, lights, ambient, phong, origin, t, i, values, res);
        } 
      }       
    }
  }
}

void specular_hl(SPHERE object, VEC3 ray, vector<LIGHT> lights, VEC3 ambient, float phong,
const VEC3 origin, float t, int i, float* values, float res)
{
  LIGHT light;
  VEC3 color, pt, v, n, l, h; color.setZero();
  pt = origin + t*ray;
  int light_count = lights.size();
  
  for(int j = 0; j < light_count; j++)
  {
    light = lights[j];
    v = origin - pt; v.normalize();
    n = pt - object.center; n.normalize();
    l = light.position - pt; l.normalize();
    h = v+l; h.normalize();

    float d_val = fmax(0, n.dot(l));
    float s_val = pow(fmax(0, n.dot(h)),phong);
  
    for(int k = 0; k < 3; k++)
    {
      color[k] += object.color[k]*(ambient[k] + light.color[k]*d_val + light.color[k]*s_val);
    }  
  } 
  set_pixel_color(values, 255*color, i);
}

void multiple_lights(SPHERE object, VEC3 ray, vector<LIGHT> lights, const VEC3 origin, 
float t, int i, float* values, float res)
{
  LIGHT light;
  VEC3 color, pt, n, l; color.setZero();
  pt = origin + t*ray;
  int light_count = lights.size();
  
  for(int j = 0; j < light_count; j++)
  {
    light = lights[j];
    n = pt - object.center; n.normalize();
    l = light.position - pt; l.normalize();
  
    float val = fmax(0, n.dot(l));
  
    for(int k = 0; k < 3; k++)
    {
      color[k] += object.color[k]*light.color[k]*val;
    }
  }
  set_pixel_color(values, 255*color, i);
}

void diffuse_single_light(SPHERE object, VEC3 ray, LIGHT& light, const VEC3 origin, 
float t, int i, float* values, float res)
{
  VEC3 color, pt, n, l; color.setZero();
  pt = origin + t*ray;
  n = pt - object.center; n.normalize();
  l = light.position - pt; l.normalize();
  
  float val = fmax(0, n.dot(l));
  
  for(int k = 0; k < 3; k++)
  {
    color[k] = object.color[k]*light.color[k]*val;
  }

  set_pixel_color(values, 255*color, i);
}

void generate_rays(vector<VEC3>&rays, VEC3 eye, VEC3 w, VEC3 u, VEC3 v, float l, float r, float t, float b, float nx, float ny, 
int xRes, int yRes)
{ 
  VEC3 direction; direction.setZero();
  VEC3 s; s.setZero();
  float U, V, d; U = V = d = 0; 

  for(int j = 0; j < yRes; j++) 
  {
    V = b + (t - b)*(j + 0.5)/ny;
    
    for(int i = 0; i < xRes; i++)
    {
      U = l + (r - l)*(i + 0.5)/nx;
      
      s = U*u + V*v - w; s.normalized();
      d = eye.dot(s);
      d = 1.0; // debug: change d so d =1
      
      direction =  d*w + U*u + V*v;  //neg
      rays.push_back(-direction);
    } 
  }

}

void show_xs(vector<VEC3>rays, int xRes, int yRes, int size, float* values, const char* filename)
{
  memset(values,0,sizeof(float)*size);
  int index = 0; float val = 0;
  for(int y = 0; y < yRes; y++)
  {
    for(int x = 0; x < xRes; x++)
    {
      index = y*xRes + x;
      val = rays[index].x(); 
      index *= 3;

      if(val > 0)// > 0)
      {
        values[index] = fmin(255, val);
      }
    }
  }

  writePPM(filename, xRes, yRes, values);
}

void show_xabs(vector<VEC3>rays, int xRes, int yRes, int size, float* values, const char* filename)
{
  memset(values,0,sizeof(float)*size);
  int index = 0; float val = 0;
  for(int y = 0; y < yRes; y++)
  {
    for(int x = 0; x < xRes; x++)
    {
      index = y*xRes + x;
      val = abs(rays[index].x()); 
      index *= 3;

      values[index] = fmin(255, val);
    }
  }

  writePPM(filename, xRes, yRes, values);
}

void show_ys(vector<VEC3>rays, int xRes, int yRes, int size, float* values, const char* filename)
{
  memset(values,0,sizeof(float)*size);
  int index = 0; float val = 0;
  for(int y = 0; y < yRes; y++)
  {
    for(int x = 0; x < xRes; x++)
    {
      index = y*xRes + x;
      val = rays[index].y(); 
      index *= 3;

      if(val > 0)
      {
        values[index+1] = fmin(255, val);
      }
    }
  }

  writePPM(filename, xRes, yRes, values);
}

void show_yabs(vector<VEC3>rays, int xRes, int yRes, int size, float* values, const char* filename)
{
  memset(values,0,sizeof(float)*size);
  int index = 0; float val = 0;
  for(int y = 0; y < yRes; y++)
  {
    for(int x = 0; x < xRes; x++)
    {
      index = y*xRes + x;
      val = abs(rays[index].y()); 
      index *= 3;
      
      values[index+1] = fmin(255, val);
    }
  }

  writePPM(filename, xRes, yRes, values);
}

void show_ray_comps(vector<VEC3>rays, float* values, int xRes, int yRes, int size)
{
  scale_vecs(rays, 255);
  
  show_xs(rays, xRes, yRes, size, values, "1x.png");
  show_xabs(rays, xRes, yRes, size, values, "1xabs.png");
  show_ys(rays, xRes, yRes, size, values, "1y.png");
  show_yabs(rays, xRes, yRes, size, values, "1yabs.png");

  memset(values,0,sizeof(float)*size);
  scale_vecs(rays, 1/255);
}

void set_pixel_color(float* values, const VEC3 color, int pixel)
{
  int index = 3*pixel;
  
  for(int j = 0; j < 3; j++)
  {
    values[index+j] += fmin(255, fmax(0, color[j]));
  }
  
}


/* scale_vecs()
 * Scales all VEC4s in the vector by the scalar
 */
void scale_vecs(vector<VEC3>& vecs, float scalar)
{
  int vec_count = vecs.size();
  
  for(int i = 0; i < vec_count; i++) 
  {
    vecs[i] *= scalar;
  }
}

///////////////////////////////////////////////////////////////////////
// Given
///////////////////////////////////////////////////////////////////////
void readPPM(const string& filename, int& xRes, int& yRes, float*& values)
{
  // try to open the file
  FILE *fp;
  fp = fopen(filename.c_str(), "rb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for reading." << endl;
    cout << " Make sure you're not trying to read from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  // get the dimensions
  fscanf(fp, "P6\n%d %d\n255\n", &xRes, &yRes);
  int totalCells = xRes * yRes;

  // grab the pixel values
  unsigned char* pixels = new unsigned char[3 * totalCells];
  fread(pixels, 1, totalCells * 3, fp);

  // copy to a nicer data type
  values = new float[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = pixels[i];

  // clean up
  delete[] pixels;
  fclose(fp);
  cout << " Read in file " << filename.c_str() << endl;
}

void writePPM(const string& filename, int& xRes, int& yRes, const float* values)
{
  int totalCells = xRes * yRes;
  unsigned char* pixels = new unsigned char[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    pixels[i] = values[i];

  FILE *fp;
  fp = fopen(filename.c_str(), "wb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for writing." << endl;
    cout << " Make sure you're not trying to write from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  fprintf(fp, "P6\n%d %d\n255\n", xRes, yRes);
  fwrite(pixels, 1, totalCells * 3, fp);
  fclose(fp);
  delete[] pixels;
}