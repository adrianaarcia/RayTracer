# ScanlineRenderer
Scanline renderer implemented using C++.
Last Updated: 3/29/2020

## Dependencies
 * Eigen (http://eigen.tuxfamily.org/).

## Description
This is a ray tracer, implemented using C++.

This code includes:

* Eye ray generation
* Scene intersection
* Diffuse shading
* Lighting (multiple lights)
* Specular shading

## Running the Code
You can run the program with the sample scene geometry provided in the code simply by running the executable render.sh. This will produce images demonstrating each of the features mentioned above.

```
chmod +x render.sh
./render.sh
```

## Potential Improvements/Additions
* Shadows
* Reflection
* Refraction
* Fresnel effect
* Triangle Intersection
* Reflections for triangles
* Mandelbrot texturing
