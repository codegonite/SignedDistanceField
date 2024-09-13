# Signed Distance Field (SDF) JavaScript Library

A Signed Distance Field 3D modeler implemented using JavaScript.

## Overview

Signed Distance Fields (SDFs) represent 3D surfaces as implicit functions, which can be used to define and manipulate complex shapes. This library provides tools to perform boolean operations (intersection, union, difference) on SDFs and for convertion to meshes using difference algorithms.

Key features include:

- **Boolean Operations**: Perform intersection, union, and difference on SDFs.
- **Mesh Conversion Using Surface Nets**: Convertion of signed distance fields to meshes using the surface nets algorithm.
- **Mesh Output as STL or OBJ**: Meshes can be exported in stl or obj format.

For more information on implicit functions and SDFs, you can refer to [this Wikipedia article](https://en.wikipedia.org/wiki/Implicit_function).

## License

[GNU GENERAL PUBLIC LICENSE (GNU)](https://github.com/codegonite/SignedDistanceFieldJS/blob/master/LICENSE)