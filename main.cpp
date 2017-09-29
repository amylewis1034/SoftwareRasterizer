#include <iostream>
#include <string>
#include <vector>
#include <memory>

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#include "Image.h"

// This allows you to skip the `std::` in front of C++ standard library
// functions. You can also say `using std::cout` to be more selective.
// You should never do this in a header file.
using namespace std;

int g_width, g_height;

typedef struct _color{
   int r;
   int g;
   int b;
} Color;

typedef struct _vertex {
   float x;
   float y;
   float z;
} Vertex;

typedef struct _bounding_box {
   float min_x;
   float min_y;
   float max_x;
   float max_y;
} BoundingBox;

typedef struct _triangle {
   Vertex v1;
   Vertex v2;
   Vertex v3;
   BoundingBox bb;
} Triangle;

/*
   Helper function you will want all quarter
   Given a vector of shapes which has already been read from an obj file
   resize all vertices to the range [-1, 1]
 */
void resize_obj(std::vector<tinyobj::shape_t> &shapes){
   float minX, minY, minZ;
   float maxX, maxY, maxZ;
   float scaleX, scaleY, scaleZ;
   float shiftX, shiftY, shiftZ;
   float epsilon = 0.001;

   minX = minY = minZ = 1.1754E+38F;
   maxX = maxY = maxZ = -1.1754E+38F;

   //Go through all vertices to determine min and max of each dimension
   for (size_t i = 0; i < shapes.size(); i++) {
      for (size_t v = 0; v < shapes[i].mesh.positions.size() / 3; v++) {
         if(shapes[i].mesh.positions[3*v+0] < minX) minX = shapes[i].mesh.positions[3*v+0];
         if(shapes[i].mesh.positions[3*v+0] > maxX) maxX = shapes[i].mesh.positions[3*v+0];

         if(shapes[i].mesh.positions[3*v+1] < minY) minY = shapes[i].mesh.positions[3*v+1];
         if(shapes[i].mesh.positions[3*v+1] > maxY) maxY = shapes[i].mesh.positions[3*v+1];

         if(shapes[i].mesh.positions[3*v+2] < minZ) minZ = shapes[i].mesh.positions[3*v+2];
         if(shapes[i].mesh.positions[3*v+2] > maxZ) maxZ = shapes[i].mesh.positions[3*v+2];
      }
   }

	//From min and max compute necessary scale and shift for each dimension
   float maxExtent, xExtent, yExtent, zExtent;
   xExtent = maxX-minX;
   yExtent = maxY-minY;
   zExtent = maxZ-minZ;
   if (xExtent >= yExtent && xExtent >= zExtent) {
      maxExtent = xExtent;
   }
   if (yExtent >= xExtent && yExtent >= zExtent) {
      maxExtent = yExtent;
   }
   if (zExtent >= xExtent && zExtent >= yExtent) {
      maxExtent = zExtent;
   }
   scaleX = 2.0 /maxExtent;
   shiftX = minX + (xExtent/ 2.0);
   scaleY = 2.0 / maxExtent;
   shiftY = minY + (yExtent / 2.0);
   scaleZ = 2.0/ maxExtent;
   shiftZ = minZ + (zExtent)/2.0;

   //Go through all verticies shift and scale them
   for (size_t i = 0; i < shapes.size(); i++) {
      for (size_t v = 0; v < shapes[i].mesh.positions.size() / 3; v++) {
         shapes[i].mesh.positions[3*v+0] = (shapes[i].mesh.positions[3*v+0] - shiftX) * scaleX;
         assert(shapes[i].mesh.positions[3*v+0] >= -1.0 - epsilon);
         assert(shapes[i].mesh.positions[3*v+0] <= 1.0 + epsilon);
         shapes[i].mesh.positions[3*v+1] = (shapes[i].mesh.positions[3*v+1] - shiftY) * scaleY;
         assert(shapes[i].mesh.positions[3*v+1] >= -1.0 - epsilon);
         assert(shapes[i].mesh.positions[3*v+1] <= 1.0 + epsilon);
         shapes[i].mesh.positions[3*v+2] = (shapes[i].mesh.positions[3*v+2] - shiftZ) * scaleZ;
         assert(shapes[i].mesh.positions[3*v+2] >= -1.0 - epsilon);
         assert(shapes[i].mesh.positions[3*v+2] <= 1.0 + epsilon);
      }
   }
}

void computeViewandScales(float *scaleX, float *shiftX, float *scaleY, float *shiftY) {
   float left, right, top, bottom;
   
   if (g_width > g_height) {
      left = -1.0 * (g_width / g_height);
      right = g_width / g_height;
      top = 1.0;
      bottom = -1.0;
   } else {
      left = -1.0;
      right = 1.0;
      top = g_height / g_width;
      bottom = -1.0 * (g_height / g_width);
   }
   
   *scaleX = (g_width - 1.0) / (right - left);
   *shiftX = *scaleX * -1.0 * left;
   *scaleY = (g_height - 1.0) / (top - bottom);
   *shiftY = *scaleY * -1.0 * bottom;
}

int worldToPixelX(float inX, float scaleX, float shiftX) {
   return inX * scaleX + shiftX;
}

int worldToPixelY(float inY, float scaleY, float shiftY) {
   return inY * scaleY + shiftY;
}

void fillTrianglesVector(vector<Triangle *> *triangles, vector<float> *posBuf, vector<unsigned int> *triBuf, float *scaleX, float *shiftX, float *scaleY, float *shiftY) {
   for (unsigned int i = 0; i < triBuf->size()/3; i++) {
      Triangle *tri = (Triangle *) malloc(sizeof(Triangle));
      
      float indx0 = triBuf->at(i*3 + 0);
      float indx1 = triBuf->at(i*3 + 1);
      float indx2 = triBuf->at(i*3 + 2);
      
      tri->v1.x = worldToPixelX(posBuf->at(indx0*3 + 0), *scaleX, *shiftX);
      tri->v1.y = worldToPixelY(posBuf->at(indx0*3 + 1), *scaleY, *shiftY);
      tri->v1.z = posBuf->at(indx0*3 + 2);
      
      tri->v2.x = worldToPixelX(posBuf->at(indx1*3 + 0), *scaleX, *shiftX);
      tri->v2.y = worldToPixelY(posBuf->at(indx1*3 + 1), *scaleY, *shiftY);
      tri->v2.z = posBuf->at(indx1*3 + 2);
      
      tri->v3.x = worldToPixelX(posBuf->at(indx2*3 + 0), *scaleX, *shiftX);
      tri->v3.y = worldToPixelY(posBuf->at(indx2*3 + 1), *scaleY, *shiftY);
      tri->v3.z = posBuf->at(indx2*3 + 2);
      
      tri->bb.min_x = min(tri->v1.x, min(tri->v2.x, tri->v3.x));
      tri->bb.min_y = min(tri->v1.y, min(tri->v2.y, tri->v3.y));
      tri->bb.max_x = max(tri->v1.x, max(tri->v2.x, tri->v3.x));
      tri->bb.max_y = max(tri->v1.y, max(tri->v2.y, tri->v3.y)); 

      triangles->push_back(tri);
   }            
}

float computeTriangleArea(vector<Triangle *> *triangles, int i) {
   return .5 * ((triangles->at(i)->v2.x - triangles->at(i)->v1.x) * (triangles->at(i)->v3.y - triangles->at(i)->v1.y) - (triangles->at(i)->v3.x - triangles->at(i)->v1.x) * (triangles->at(i)->v2.y - triangles->at(i)->v1.y));
}

float computeAlpha(vector<Triangle *> *triangles, int x, int y, float area_tri, int i) {
   float area_alpha = .5 * ((x - triangles->at(i)->v3.x) * (triangles->at(i)->v2.y - triangles->at(i)->v3.y) - (triangles->at(i)->v2.x - triangles->at(i)->v3.x) * (y - triangles->at(i)->v3.y));
   return area_alpha / area_tri;
}

float computeBeta(vector<Triangle *> *triangles, int x, int y, float area_tri, int i) {
   float area_beta = .5 * ((triangles->at(i)->v1.x - triangles->at(i)->v3.x) * (y - triangles->at(i)->v3.y) - (x - triangles->at(i)->v3.x) * (triangles->at(i)->v1.y - triangles->at(i)->v3.y));
   return area_beta / area_tri;
}

void computeColoringModeOne(vector<Triangle *> *triangles, float *cur_z, float alpha, float beta, float gamma, Color **cur_color, int i) {
   float z1_color = (triangles->at(i)->v1.z * .5 + .5) * 255;
   float z2_color = (triangles->at(i)->v2.z * .5 + .5) * 255;
   float z3_color = (triangles->at(i)->v3.z * .5 + .5) * 255;
   
   (*cur_color)->r = (alpha * z1_color) + (beta * z2_color) + (gamma * z3_color);
   (*cur_color)->g = 0;
   (*cur_color)->b = 0;
   
   *cur_z = (alpha * triangles->at(i)->v1.z) + (beta * triangles->at(i)->v2.z) + (gamma * triangles->at(i)->v3.z);
}

void computeColoringModeTwo(Color **cur_color, int y) {
   Color *color_1 = (Color *)malloc(sizeof(Color));
   color_1->r = 255;
   color_1->g = 255;
   color_1->b = 0;
   
   Color *color_2 = (Color *)malloc(sizeof(Color));
   color_2->r = 0;
   color_2->g = 0;
   color_2->b = 255;
  
   (*cur_color)->r = color_1->r * y / g_height;
   (*cur_color)->g = color_1->g * y / g_height;
   (*cur_color)->b = 1 - (color_2->b * y / g_height);
   
   free(color_1);
   free(color_2);
}

int main(int argc, char **argv) {
	if(argc < 6) {
      cout << "Usage: Assignment1 meshfile imagefile width height coloring_mode" << endl;
      return 0;
   }
   
	// OBJ filename
	string meshName(argv[1]);
	string imgName(argv[2]);
   
   // Width of image
	g_width = atoi(argv[3]);
	// Height of image
	g_height = atoi(argv[4]);
   // Coloring Mode
   int coloring_mode = atoi(argv[5]);
   
   //create an image
	auto image = make_shared<Image>(g_width, g_height);

	// triangle buffer
	vector<unsigned int> triBuf;
	// position buffer
	vector<float> posBuf;
	// Some obj files contain material information.
	// We'll ignore them for this assignment.
	vector<tinyobj::shape_t> shapes; // geometry
	vector<tinyobj::material_t> objMaterials; // material
	string errStr;
	
   bool rc = tinyobj::LoadObj(shapes, objMaterials, errStr, meshName.c_str());
   
   //keep this code to resize your object to be within -1 -> 1
   resize_obj(shapes); 
   
	/* error checking on read */
	if(!rc) {
		cerr << errStr << endl;
	} else {
		posBuf = shapes[0].mesh.positions;
		triBuf = shapes[0].mesh.indices;
	}
   
   cout << "Number of vertices: " << posBuf.size()/3 << endl;
	cout << "Number of triangles: " << triBuf.size()/3 << endl;
   
	//TODO add code to iterate through each triangle and rasterize it   
   float scaleX, shiftX, scaleY, shiftY;
   float alpha, beta, gamma;
   float epsilon = 0.001;
   float cur_z;
   
   vector<Triangle *> triangles;
   vector<float> z_buffer(g_width * g_height);
   
   computeViewandScales(&scaleX, &shiftX, &scaleY, &shiftY);
   
   fillTrianglesVector(&triangles, &posBuf, &triBuf, &scaleX, &shiftX, &scaleY, &shiftY);
   
   for (unsigned int i = 0; i < z_buffer.size(); i++) {
      z_buffer.at(i) = -1 * INFINITY;
   }
   
   for (unsigned int i = 0; i < triangles.size(); i++) {
      float area_tri = computeTriangleArea(&triangles, i);
      
      for(int y = triangles.at(i)->bb.min_y; y <= triangles.at(i)->bb.max_y; ++y) {
         for(int x = triangles.at(i)->bb.min_x; x <= triangles.at(i)->bb.max_x; ++x) {              
            alpha = computeAlpha(&triangles, x, y, area_tri, i);
            beta = computeBeta(&triangles, x, y, area_tri, i); 
            gamma = 1 - alpha - beta;
            Color *cur_color = (Color *)malloc(sizeof(Color)); 
            
            if (-1.0*epsilon < alpha && alpha < 1+epsilon && -1.0*epsilon < beta && beta < 1+epsilon && -1.0*epsilon < gamma && gamma < 1+epsilon) {
               if (coloring_mode == 1) {
                  computeColoringModeOne(&triangles, &cur_z, alpha, beta, gamma, &cur_color, i);
                  
                  if (z_buffer.at(g_width*y + x) < cur_z) {
                     image->setPixel(x, y, cur_color->r, cur_color->g, cur_color->b);
                     z_buffer.at(g_width*y + x) = cur_z;
                  }
               } else if (coloring_mode == 2) {
                  computeColoringModeTwo(&cur_color, y);
                  image->setPixel(x, y, cur_color->r, cur_color->g, cur_color->b);
               }
            }
            
            free(cur_color);
         }
      }

      free(triangles.at(i));
   }
   	
	//write out the image
   image->writeToFile(imgName);

	return 0;
}
