Source file for mesh class



// Include files

#include "R3Mesh.h"
#include <map>
#include <set>
#include "FileLogger.h"

////////////////////////////////////////////////////////////
// METHODS TO IMPLEMENT -- Twist() (1st Project), Taubin, Loop (2nd Project)
////////////////////////////////////////////////////////////


void R3Mesh::
Twist(double angle)
{
  // Twist mesh by an angle, or other simple mesh warpi25ng of your choice.
  // See Scale() for how to get the vertex positions, and see bbox for the bounding box.
  
  // FILL IN IMPLEMENTATION HERE

  // Update mesh data structures
  Update();
}

void R3Mesh::
Taubin(double lambda, double mu, int iters)
{
  // Apply Taubin smoothing to the mesh for a given number of iterations. I suggest to make a single method that smoothes by a positive or negative amount, and then call Smooth(lambda) and Smooth(mu) in an iterative loop.
  // See Scale() for how to get the vertex positions.
  
  // FILL IN IMPLEMENTATION HERE
	Logger logger;
	logger.set_file("Taubin.log");
	logger.set_level(0);
	logger.debug()<<"Loop through all faces to create the map!"<<std::endl;

	// a map that maps vertex and its neighbors.
	map< R3MeshVertex*,vector<R3MeshVertex*> > taubinMap;
	for (int i = 0; i < NFaces(); i++) {
		R3MeshFace *f = Face(i);
		logger.debug()<<"Face: "<<f->id<<" Num of Vertices: "<<f->vertices.size()<<std::endl;
		for(vector<R3MeshVertex *>::iterator it = f->vertices.begin(); it != f->vertices.end(); ++it ) {
			logger.info()<<"    Vertex "<<(*it)->id<<std::endl;		
			if( it == f->vertices.begin() ) {
				logger.info()<<"        it == begin()"<<endl;
				R3MeshVertex* tmp = f->vertices.back(); 
				if(std::find(taubinMap[(*it)].begin(),taubinMap[(*it)].end(), tmp) == taubinMap[(*it)].end())
					taubinMap[(*it)].push_back((tmp));
				if(std::find(taubinMap[(*it)].begin(),taubinMap[(*it)].end(), *(it+1)) == taubinMap[(*it)].end())
					taubinMap[(*it)].push_back(*(it+1));
			}
			else if( it == f->vertices.end()-1 ) {
				logger.info()<<"        it == end()"<<endl;
				if(std::find(taubinMap[(*it)].begin(),taubinMap[(*it)].end(), *(it-1)) == taubinMap[(*it)].end())
					taubinMap[(*it)].push_back(*(it-1));
				R3MeshVertex* tmp = f->vertices.front(); //  *(f->vertices.begin());
				if(std::find(taubinMap[(*it)].begin(),taubinMap[(*it)].end(), tmp) == taubinMap[(*it)].end())
					taubinMap[(*it)].push_back((tmp));
			}
			else {
				logger.info()<<"        else"<<endl;
				if(std::find(taubinMap[(*it)].begin(),taubinMap[(*it)].end(), *(it-1)) == taubinMap[(*it)].end())
					taubinMap[(*it)].push_back(*(it-1));
				if(std::find(taubinMap[(*it)].begin(),taubinMap[(*it)].end(), *(it+1)) == taubinMap[(*it)].end())
					taubinMap[(*it)].push_back(*(it+1));
			}
			
		}
	}

	map<R3MeshVertex*,R3MeshVertex> vertexMap;  // a copy of vertex that maps the vertex and its copy
	for( int i = 0; i < NVertices(); i++ ) {
		R3MeshVertex * vv = Vertex(i);
		R3MeshVertex tmp;
		tmp.id = vv->id;
		tmp.normal = vv->normal;
		tmp.position = vv->position;
		tmp.texcoords = vv->texcoords;
		vertexMap[vv] = tmp;
	}

	for( int loop = 0; loop < iters; loop ++ ) {
		logger.debug()<<"Loop through all the vertices"<<std::endl;
		for(int i = 0; i < NVertices(); i++ ) {
			R3MeshVertex * vv = Vertex(i);
			logger.info()<<"Vertex: "<<vv->id<<" ("<<vv->position[0]<<","<<vv->position[1]<<","<<vv->position[2]<<")"<<std::endl;
			// delta_Vi = sum(Wij*(Vj-Vi));
			double delta_x = 0.0, delta_y = 0.0, delta_z = 0.0;
			double weight = 1.0/taubinMap[vv].size();
			for (vector<R3MeshVertex *>::iterator it = taubinMap[vv].begin(); it != taubinMap[vv].end(); ++it ) {
				delta_x += weight*(vertexMap[(*it)].position[0]-vv->position[0]);
				delta_y += weight*(vertexMap[(*it)].position[1]-vv->position[1]);
				delta_z += weight*(vertexMap[(*it)].position[2]-vv->position[2]);
			}
			// Vi' = Vi + lamda*delta_Vi;
			double new_x, new_y, new_z;
			new_x = vv->position[0] + lambda*delta_x;
			new_y = vv->position[1] + lambda*delta_y;
			new_z = vv->position[2] + lambda*delta_z;
			
			R3Point newPos(new_x, new_y, new_z);
			vv->position = newPos;
			logger.info()<<"Newvv: "<<vv->id<<" ("<<vv->position[0]<<","<<vv->position[1]<<","<<vv->position[2]<<")"<<std::endl;
		}
		for( int i = 0; i < NVertices(); i++ ) {
			R3MeshVertex * vv = Vertex(i);
			vertexMap[vv].position = vv->position;
		}
		for(int i = 0; i < NVertices(); i++ ) {
			R3MeshVertex * vv = Vertex(i);
			double delta_x = 0.0, delta_y = 0.0, delta_z = 0.0;
			double weight = 1.0/taubinMap[vv].size();
			for (vector<R3MeshVertex *>::iterator it = taubinMap[vv].begin(); it != taubinMap[vv].end(); ++it ) {
				delta_x += weight*(vertexMap[(*it)].position[0]-vv->position[0]);
				delta_y += weight*(vertexMap[(*it)].position[1]-vv->position[1]);
				delta_z += weight*(vertexMap[(*it)].position[2]-vv->position[2]);
			}
			double new_x, new_y, new_z;
			new_x = vv->position[0] + mu*delta_x;
			new_y = vv->position[1] + mu*delta_y;
			new_z = vv->position[2] + mu*delta_z;
			
			R3Point newPos(new_x, new_y, new_z);
			vv->position = newPos;
		}
		for( int i = 0; i < NVertices(); i++ ) {
			R3MeshVertex * vv = Vertex(i);
			vertexMap[vv].position = vv->position;
		}
	}

  // Update mesh data structures
  Update();
}

void R3Mesh::
Loop()
{
  // Perform Loop subdivision on a triangular mesh. Faces that are not triangles can be skipped.
  
  // FILL IN IMPLEMENTATION HERE
	Logger logger;
	logger.set_file("Loop.log");
	logger.set_level(0);
	logger.debug()<<"Loop through all faces to create the map!"<<std::endl;

	// !New mesh faces after subdivision
	vector<R3MeshFace *> facesNew;

	// !record the number of old vertices, id of the vertex before oldVNum is old vertices.
	int oldVNum = NVertices();
	
	// oddMap maps edge( weight of the two point on the edge is 3/8) and two points on the adjacent faces( weight of the two points are 1/8) 
	// use set can avoid duplicate
	map<pair<R3MeshVertex*,R3MeshVertex*>,set<R3MeshVertex*> > oddMap;

	for (int i = 0; i < NFaces(); i++) {
		R3MeshFace *f = Face(i);
		logger.info()<<"Face: "<<f->id<<" has "<<f->vertices.size()<<" vertices."<<std::endl;
		if( f->vertices.size() > 3 ) continue;
		for(vector<R3MeshVertex *>::iterator it = f->vertices.begin(); it != f->vertices.end(); ++it ) {
			R3MeshVertex* tmp;
			if( it == f->vertices.begin() ) tmp = *(f->vertices.end()-1);
			else tmp = *(it-1);
			pair<R3MeshVertex*,R3MeshVertex*> pair_key1 = make_pair(tmp,*it);
			pair<R3MeshVertex*,R3MeshVertex*> pair_key2 = make_pair(*it,tmp);
			pair<R3MeshVertex*,R3MeshVertex*> pair_key = pair_key1;
			if( oddMap.find(pair_key1) == oddMap.end() && !(oddMap.find(pair_key2) == oddMap.end()) )
				pair_key = pair_key2;
			logger.info()<<"    Edge: ("<<pair_key.first->id<<","<<pair_key.second->id<<")"<<std::endl;
			for(vector<R3MeshVertex *>::iterator itt = f->vertices.begin(); itt != f->vertices.end(); ++itt ) {
				if( pair_key.first == *itt || pair_key.second == *itt ) continue;
				oddMap[pair_key].insert(*itt);
			}
		}
	}

	// !vertexMap maps the edge and new subdivision point.
	// check the odd map
	// create the new vertex map ( 1 edge -> 1 vertex ) (pair<R3MeshVertex*,R3MeshVertex*>, R3MeshVertex)
	map<pair<R3MeshVertex*,R3MeshVertex*>,R3MeshVertex*> vertexMap;
	logger.debug()<<"Odd Map"<<std::endl;
	for(map<pair<R3MeshVertex*,R3MeshVertex*>,set<R3MeshVertex*> >::iterator it = oddMap.begin(); it != oddMap.end(); ++it ) {
		logger.info()<<"    Edge: ("<<(*it).first.first->id<<","<<(*it).first.second->id<<")"<<std::endl;
		double weight_n = 3.0/8.0, weight_o = 1.0/8.0;
		if( (*it).second.size() == 1 ) { weight_n = 1.0/2.0; weight_o = 0;} //boundary
		R3Point tmp(0.0,0.0,0.0);
		tmp = tmp + weight_n*(*it).first.first->position + weight_n*(*it).first.second->position;
		for(set<R3MeshVertex*>::iterator itt = (*it).second.begin(); itt != (*it).second.end(); ++itt) {
			logger.info()<<"    "<<(*itt)->id<<" ("<<(*itt)->position[0]<<","<<(*itt)->position[1]<<","<<(*itt)->position[2]<<")"<<std::endl;
			tmp = tmp + weight_o * (*itt)->position;
		}
		R3MeshVertex *vertex = new R3MeshVertex(tmp, R3zero_vector, R2zero_point);
		bbox.Union(tmp);
		vertex->id = vertices.size();
		vertices.push_back(vertex);
		vertexMap[(*it).first] = vertex;
	}

	for (int i = 0; i < NFaces(); i++) {
		// each face can generate 4 sub-faces
		// three of those faces, each face contain one old vertex and two new vertices
		// one face is composed of three new vertices called center_face.
		R3MeshFace *f = Face(i);
		logger.info()<<"Face: "<<f->id<<" has "<<f->vertices.size()<<" vertices."<<std::endl;
		if( f->vertices.size() > 3 ) continue;

		vector<R3MeshVertex *> center_face;

		for(vector<R3MeshVertex *>::iterator it = f->vertices.begin(); it != f->vertices.end(); ++it ) {
			// for each (*it), one new subdivision face can be generated
			vector<R3MeshVertex *> face_vertices; //each face contain 3 vertices.
			face_vertices.push_back((*it));  // add the first vertex, the (*it) itself;
			R3MeshVertex* tmp1;  //the second vertex is the mid point of edge (*it)--(*tmp1);
			R3MeshVertex* tmp2;  //the third vertex is the mid point of edge (*it)--(*tmp2);
			if( it == f->vertices.end()-1 ) tmp1 = *(f->vertices.begin());
			else tmp1 = *(it+1);
			if( it == f->vertices.begin() ) tmp2 = *(f->vertices.end()-1);
			else tmp2 = *(it-1);

			//find the mid-point of edge (*it)--(tmp1)
			pair<R3MeshVertex*,R3MeshVertex*> pair_key1 = make_pair(tmp1,*it);
			pair<R3MeshVertex*,R3MeshVertex*> pair_key2 = make_pair(*it,tmp1);
			pair<R3MeshVertex*,R3MeshVertex*> pair_key = pair_key1;
			if( oddMap.find(pair_key1) == oddMap.end() && !(oddMap.find(pair_key2) == oddMap.end()) )
				pair_key = pair_key2;
			face_vertices.push_back(vertexMap[pair_key]); // add the second  vertex which is the mid-point of vertex (*it) and its succeed neighbor (tmp1)

			center_face.push_back(vertexMap[pair_key]);

			//find the mid-point of edge (*it)--(tmp2)
			pair<R3MeshVertex*,R3MeshVertex*> pair_key3 = make_pair(tmp2,*it);
			pair<R3MeshVertex*,R3MeshVertex*> pair_key4 = make_pair(*it,tmp2);
			pair_key = pair_key3;
			if( oddMap.find(pair_key3) == oddMap.end() && !(oddMap.find(pair_key4) == oddMap.end()) )
				pair_key = pair_key4;
			face_vertices.push_back(vertexMap[pair_key]); // add the third  vertex which is the mid-point of vertex (*it) and its previous neighbor (tmp2)

			//add the face 
			R3MeshFace *face = new R3MeshFace(face_vertices);
			face->id = facesNew.size();
			facesNew.push_back(face);
		}
		//add the center_face 
		R3MeshFace *face = new R3MeshFace(center_face);
		face->id = facesNew.size();
		facesNew.push_back(face);	
	}

	// !Update the old vertices
	// !Even
	map<R3MeshVertex*,set<R3MeshVertex*> > evenMap;
	for( int i = 0; i < facesNew.size(); i++ ) {
		R3MeshFace *f = facesNew[i];
		for(vector<R3MeshVertex *>::iterator it = f->vertices.begin(); it != f->vertices.end(); ++it ) {
			R3MeshVertex *neibor1;
			R3MeshVertex *neibor2;
			if( it == f->vertices.begin() ) { neibor1 = *(f->vertices.end()-1); neibor2 = *(it+1); }
			else if ( it == f->vertices.end()-1 ) { neibor1 = *(it-1); neibor2 = *(f->vertices.begin()); }
			else { neibor1 = *(it-1); neibor2 = *(it+1); }

			evenMap[(*it)].insert(neibor1);
			evenMap[(*it)].insert(neibor2);
		}
	}

	for( int i = 0; i < oldVNum; i++ ) {
		R3MeshVertex* v = Vertex(i);
		// find v's neighbors and update v's position
		double beta; // choose beta: Warren
		if( evenMap[v].size() > 3 ) beta = 3.0/(8.0*evenMap[v].size());
		else beta = 3.0/16.0;
		R3Point pos = (1.0-evenMap[v].size()*beta)*(v->position); //(1-k*beta)
		for(set<R3MeshVertex*>::iterator it = evenMap[v].begin(); it != evenMap[v].end(); ++it ) {
			pos = pos + beta*(*it)->position;
		}
		v->position = pos;
	}
	
	// !delete the old faces
	for( int i = 0; i < NFaces(); i++ )
		DeleteFace(faces[i]);
	faces = facesNew;
  // Update mesh data structures
  Update();
}


////////////////////////////////////////////////////////////
// MESH CONSTRUCTORS/DESTRUCTORS
////////////////////////////////////////////////////////////

R3Mesh::
R3Mesh(void)
  : bbox(R3null_box)
{
}



R3Mesh::
R3Mesh(const R3Mesh& mesh)
  : bbox(R3null_box)
{
  // Create vertices
  for (int i = 0; i < mesh.NVertices(); i++) {
    R3MeshVertex *v = mesh.Vertex(i);
    CreateVertex(v->position, v->normal, v->texcoords);
  }

  // Create faces
  for (int i = 0; i < mesh.NFaces(); i++) {
    R3MeshFace *f = mesh.Face(i);
    vector<R3MeshVertex *> face_vertices;
    for (unsigned int j = 0; j < f->vertices.size(); j++) {
      R3MeshVertex *ov = f->vertices[j];
      R3MeshVertex *nv = Vertex(ov->id);
      face_vertices.push_back(nv);
    }
    CreateFace(face_vertices);
  }
}



R3Mesh::
~R3Mesh(void)
{
  // Delete faces
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *f = Face(i);
    delete f;
  }

  // Delete vertices
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *v = Vertex(i);
    delete v;
  }
}



////////////////////////////////////////////////////////////
// MESH PROPERTY FUNCTIONS
////////////////////////////////////////////////////////////

R3Point R3Mesh::
Center(void) const
{
  // Return center of bounding box
  return bbox.Centroid();
}



double R3Mesh::
Radius(void) const
{
  // Return radius of bounding box
  return bbox.DiagonalRadius();
}



////////////////////////////////////////////////////////////
// MESH PROCESSING FUNCTIONS
////////////////////////////////////////////////////////////

void R3Mesh::
Translate(double dx, double dy, double dz)
{
  // Translate the mesh by adding a 
  // vector (dx,dy,dz) to every vertex

  // This is implemented for you as an example 

  // Create a translation vector
  R3Vector translation(dx, dy, dz);

  // Update vertices
  for (unsigned int i = 0; i < vertices.size(); i++) {
    R3MeshVertex *vertex = vertices[i];
    vertex->position.Translate(translation);
  }

  // Update mesh data structures
  Update();
}




void R3Mesh::
Scale(double sx, double sy, double sz)
{
  // Scale the mesh by increasing the distance 
  // from every vertex to the origin by a factor 
  // given for each dimension (sx, sy, sz)

  // This is implemented for you as an example 

  // Update vertices
  for (unsigned int i = 0; i < vertices.size(); i++) {
    R3MeshVertex *vertex = vertices[i];
    vertex->position[0] *= sx;
    vertex->position[1] *= sy;
    vertex->position[2] *= sz;
  }

  // Update mesh data structures
  Update();
}




void R3Mesh::
Rotate(double angle, const R3Line& axis)
{
  // Rotate the mesh counter-clockwise by an angle 
  // (in radians) around a line axis

  // This is implemented for you as an example 

  // Update vertices
  for (unsigned int i = 0; i < vertices.size(); i++) {
    R3MeshVertex *vertex = vertices[i];
    vertex->position.Rotate(axis, angle);
  }

  // Update mesh data structures
  Update();
}


////////////////////////////////////////////////////////////
// MESH ELEMENT CREATION/DELETION FUNCTIONS
////////////////////////////////////////////////////////////

R3MeshVertex *R3Mesh::
CreateVertex(const R3Point& position, const R3Vector& normal, const R2Point& texcoords)
{
  // Create vertex
  R3MeshVertex *vertex = new R3MeshVertex(position, normal, texcoords);

  // Update bounding box
  bbox.Union(position);

  // Set vertex ID
  vertex->id = vertices.size();

  // Add to list
  vertices.push_back(vertex);

  // Return vertex
  return vertex;
}



R3MeshFace *R3Mesh::
CreateFace(const vector<R3MeshVertex *>& vertices)
{
  // Create face
  R3MeshFace *face = new R3MeshFace(vertices);

  // Set face  ID
  face->id = faces.size();

  // Add to list
  faces.push_back(face);

  // Return face
  return face;
}



void R3Mesh::
DeleteVertex(R3MeshVertex *vertex)
{
  // Remove vertex from list
  for (unsigned int i = 0; i < vertices.size(); i++) {
    if (vertices[i] == vertex) {
      vertices[i] = vertices.back();
      vertices[i]->id = i;
      vertices.pop_back();
      break;
    }
  }

  // Delete vertex
  delete vertex;
}



void R3Mesh::
DeleteFace(R3MeshFace *face)
{
  // Remove face from list
  for (unsigned int i = 0; i < faces.size(); i++) {
    if (faces[i] == face) {
      faces[i] = faces.back();
      faces[i]->id = i;
      faces.pop_back();
      break;
    }
  }

  // Delete face
  delete face;
}



////////////////////////////////////////////////////////////
// UPDATE FUNCTIONS
////////////////////////////////////////////////////////////

void R3Mesh::
Update(void)
{
  // Update everything
  UpdateBBox();
  UpdateFacePlanes();
  UpdateVertexNormals();
  UpdateVertexCurvatures();
}



void R3Mesh::
UpdateBBox(void)
{
  // Update bounding box
  bbox = R3null_box;
  for (unsigned int i = 0; i < vertices.size(); i++) {
    R3MeshVertex *vertex = vertices[i];
    bbox.Union(vertex->position);
  }
}



void R3Mesh::
UpdateVertexNormals(void)
{
  // Update normal for every vertex
  for (unsigned int i = 0; i < vertices.size(); i++) {
    vertices[i]->UpdateNormal();
  }
}




void R3Mesh::
UpdateVertexCurvatures(void)
{
  // Update curvature for every vertex
  for (unsigned int i = 0; i < vertices.size(); i++) {
    vertices[i]->UpdateCurvature();
  }
}




void R3Mesh::
UpdateFacePlanes(void)
{
  // Update plane for all faces
  for (unsigned int i = 0; i < faces.size(); i++) {
    faces[i]->UpdatePlane();
  }
}



////////////////////////////////////////////////////////////////////////
// I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3Mesh::
Read(const char *filename)
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .ply)\n", filename);
    return 0;
  }

  // Read file of appropriate type
  int status = 0;
  if (!strncmp(extension, ".ray", 4)) 
    status = ReadRay(filename);
  else if (!strncmp(extension, ".off", 4)) 
    status = ReadOff(filename);
  else if (!strncmp(extension, ".jpg", 4)) 
    status = ReadImage(filename);
  else if (!strncmp(extension, ".jpeg", 4)) 
    status = ReadImage(filename);
  else if (!strncmp(extension, ".bmp", 4)) 
    status = ReadImage(filename);
  else if (!strncmp(extension, ".ppm", 4)) 
    status = ReadImage(filename);
  else {
    fprintf(stderr, "Unable to read file %s (unrecognized extension: %s)\n", filename, extension);
    return status;
  }

  // Update mesh data structures
  Update();

  // Return success
  return 1;
}



int R3Mesh::
Write(const char *filename)
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .ply)", filename);
    return 0;
  }

  // Write file of appropriate type
  if (!strncmp(extension, ".ray", 4)) 
    return WriteRay(filename);
  else if (!strncmp(extension, ".off", 4)) 
    return WriteOff(filename);
  else {
    fprintf(stderr, "Unable to write file %s (unrecognized extension: %s)", filename, extension);
    return 0;
  }
}



////////////////////////////////////////////////////////////
// IMAGE FILE INPUT/OUTPUT
////////////////////////////////////////////////////////////

int R3Mesh::
ReadImage(const char *filename)
{
  // Create a mesh by reading an image file, 
  // constructing vertices at (x,y,luminance), 
  // and connecting adjacent pixels into faces. 
  // That is, the image is interpretted as a height field, 
  // where the luminance of each pixel provides its z-coordinate.

  // Read image
  R2Image *image = new R2Image();
  if (!image->Read(filename)) return 0;

  // Create vertices and store in arrays
  R3MeshVertex ***vertices = new R3MeshVertex **[image->Width() ];
  for (int i = 0; i < image->Width(); i++) {
    vertices[i] = new R3MeshVertex *[image->Height() ];
    for (int j = 0; j < image->Height(); j++) {
      double luminance = image->Pixel(i, j).Luminance();
      double z = luminance * image->Width();
      R3Point position((double) i, (double) j, z);
      R2Point texcoords((double) i, (double) j);
      vertices[i][j] = CreateVertex(position, R3zero_vector, texcoords);
    }
  }

  // Create faces
  vector<R3MeshVertex *> face_vertices;
  for (int i = 1; i < image->Width(); i++) {
    for (int j = 1; j < image->Height(); j++) {
      face_vertices.clear();
      face_vertices.push_back(vertices[i-1][j-1]);
      face_vertices.push_back(vertices[i][j-1]);
      face_vertices.push_back(vertices[i][j]);
      face_vertices.push_back(vertices[i-1][j]);
      CreateFace(face_vertices);
    }
  }

  // Delete vertex arrays
  for (int i = 0; i < image->Width(); i++) delete [] vertices[i];
  delete [] vertices;

  // Delete image
  delete image;

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////
// OFF FILE INPUT/OUTPUT
////////////////////////////////////////////////////////////

int R3Mesh::
ReadOff(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s\n", filename);
    return 0;
  }

  // Read file
  int nverts = 0;
  int nfaces = 0;
  int nedges = 0;
  int line_count = 0;
  int vertex_count = 0;
  int face_count = 0;
  char buffer[1024];
  char header[64];
  while (fgets(buffer, 1023, fp)) {
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) bufferp++;

    // Skip blank lines and comments
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;

    // Check section
    if (nverts == 0) {
      // Read header keyword
      if (strstr(bufferp, "OFF")) {
        // Check if counts are on first line
        int tmp;
        if (sscanf(bufferp, "%s%d%d%d", header, &tmp, &nfaces, &nedges) == 4) {
          nverts = tmp;
        }
      }
      else {
        // Read counts from second line
        if ((sscanf(bufferp, "%d%d%d", &nverts, &nfaces, &nedges) != 3) || (nverts == 0)) {
          fprintf(stderr, "Syntax error reading header on line %d in file %s\n", line_count, filename);
          fclose(fp);
          return 0;
        }
      }
    }
    else if (vertex_count < nverts) {
      // Read vertex coordinates
      double x, y, z;
      if (sscanf(bufferp, "%lf%lf%lf", &x, &y, &z) != 3) {
        fprintf(stderr, "Syntax error with vertex coordinates on line %d in file %s\n", line_count, filename);
        fclose(fp);
        return 0;
      }

      // Create vertex
      CreateVertex(R3Point(x, y, z), R3zero_vector, R2zero_point);

      // Increment counter
      vertex_count++;
    }
    else if (face_count < nfaces) {
      // Read number of vertices in face 
      int face_nverts = 0;
      bufferp = strtok(bufferp, " \t");
      if (bufferp) face_nverts = atoi(bufferp);
      else {
        fprintf(stderr, "Syntax error with face on line %d in file %s\n", line_count, filename);
        fclose(fp);
        return 0;
      }

      // Read vertex indices for face
      vector<R3MeshVertex *> face_vertices;
      for (int i = 0; i < face_nverts; i++) {
        R3MeshVertex *v = NULL;
        bufferp = strtok(NULL, " \t");
        if (bufferp) v = Vertex(atoi(bufferp));
        else {
          fprintf(stderr, "Syntax error with face on line %d in file %s\n", line_count, filename);
          fclose(fp);
          return 0;
        }

        // Add vertex to vector
        face_vertices.push_back(v);
      }

      // Create face
      CreateFace(face_vertices);

      // Increment counter
      face_count++;
    }
    else {
      // Should never get here
      fprintf(stderr, "Found extra text starting at line %d in file %s\n", line_count, filename);
      break;
    }
  }

  // Check whether read all vertices
  if ((vertex_count != nverts) || (NVertices() < nverts)) {
    fprintf(stderr, "Expected %d vertices, but read %d vertex lines and created %d vertices in file %s\n", 
      nverts, vertex_count, NVertices(), filename);
  }

  // Check whether read all faces
  if ((face_count != nfaces) || (NFaces() < nfaces)) {
    fprintf(stderr, "Expected %d faces, but read %d face lines and created %d faces in file %s\n", 
      nfaces, face_count, NFaces(), filename);
  }

  // Close file
  fclose(fp);

  // Return number of faces read
  return NFaces();
}



int R3Mesh::
WriteOff(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open file %s\n", filename);
    return 0;
  }

  // Write header
  fprintf(fp, "OFF\n");
  fprintf(fp, "%d %d %d\n", NVertices(), NFaces(), 0);

  // Write vertices
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    const R3Point& p = vertex->position;
    fprintf(fp, "%g %g %g\n", p.X(), p.Y(), p.Z());
    vertex->id = i;
  }

  // Write Faces
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *face = Face(i);
    fprintf(fp, "%d", (int) face->vertices.size());
    for (unsigned int j = 0; j < face->vertices.size(); j++) {
      fprintf(fp, " %d", face->vertices[j]->id);
    }
    fprintf(fp, "\n");
  }

  // Close file
  fclose(fp);

  // Return number of faces
  return NFaces();
}



////////////////////////////////////////////////////////////
// RAY FILE INPUT/OUTPUT
////////////////////////////////////////////////////////////

int R3Mesh::
ReadRay(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Read body
  char cmd[128];
  int polygon_count = 0;
  int command_number = 1;
  while (fscanf(fp, "%s", cmd) == 1) {
    if (!strcmp(cmd, "#vertex")) {
      // Read data
      double px, py, pz;
      double nx, ny, nz;
      double ts, tt;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf", &px, &py, &pz, &nx, &ny, &nz, &ts, &tt) != 8) {
        fprintf(stderr, "Unable to read vertex at command %d in file %s", command_number, filename);
        return 0;
      }

      // Create vertex
      R3Point point(px, py, pz);
      R3Vector normal(nx, ny, nz);
      R2Point texcoords(ts, tt);
      CreateVertex(point, normal, texcoords);
    }
    else if (!strcmp(cmd, "#shape_polygon")) {
      // Read data
      int m, nverts;
      if (fscanf(fp, "%d%d", &m, &nverts) != 2) {
        fprintf(stderr, "Unable to read polygon at command %d in file %s", command_number, filename);
        return 0;
      }

      // Get vertices
      vector<R3MeshVertex *> face_vertices;
      for (int i = 0; i < nverts; i++) {
        // Read vertex id
        int vertex_id;
        if (fscanf(fp, "%d", &vertex_id) != 1) {
          fprintf(stderr, "Unable to read polygon at command %d in file %s", command_number, filename);
          return 0;
        }

        // Get vertex
        R3MeshVertex *v = Vertex(vertex_id);
        face_vertices.push_back(v);
      }

      // Create face
      CreateFace(face_vertices);

      // Increment polygon counter
      polygon_count++;
    }
	
    // Increment command number
    command_number++;
  }

  // Close file
  fclose(fp);

  // Return number of faces created
  return polygon_count;
}



int R3Mesh::
WriteRay(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "w"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Write vertices
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    const R3Point& p = vertex->position;
    const R3Vector& n = vertex->normal;
    const R2Point& t = vertex->texcoords;
    fprintf(fp, "#vertex %g %g %g  %g %g %g  %g %g\n", p.X(), p.Y(), p.Z(), 
      n.X(), n.Y(), n.Z(), t.X(), t.Y());
    vertex->id = i;
  }

  // Write faces
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *face = Face(i);
    int nvertices = face->vertices.size();
    fprintf(fp, "#shape_polygon 0 %d ", nvertices);
    for (int j = 0; j < nvertices; j++) {
      R3MeshVertex *v = face->vertices[j];
      fprintf(fp, "%d ", v->id);
    }
    fprintf(fp, "\n");
  }

  // Close file
  fclose(fp);

  // Return number of faces written
  return NFaces();
}



////////////////////////////////////////////////////////////
// MESH VERTEX MEMBER FUNCTIONS
////////////////////////////////////////////////////////////

R3MeshVertex::
R3MeshVertex(void)
  : position(0, 0, 0),
    normal(0, 0, 0),
    texcoords(0, 0),
    curvature(0),
    id(0)
{
}



R3MeshVertex::
R3MeshVertex(const R3MeshVertex& vertex)
  : position(vertex.position),
    normal(vertex.normal),
    texcoords(vertex.texcoords),
    curvature(vertex.curvature),
    id(0)
{
}




R3MeshVertex::
R3MeshVertex(const R3Point& position, const R3Vector& normal, const R2Point& texcoords)
  : position(position),                    
    normal(normal),
    texcoords(texcoords),
    curvature(0),
    id(0)
{
}




double R3MeshVertex::
AverageEdgeLength(void) const
{
  // Return the average length of edges attached to this vertex
  // This feature should be implemented first.  To do it, you must
  // design a data structure that allows O(K) access to edges attached
  // to each vertex, where K is the number of edges attached to the vertex.

  // FILL IN IMPLEMENTATION HERE  (THIS IS REQUIRED)
  // BY REPLACING THIS ARBITRARY RETURN VALUE
  fprintf(stderr, "Average vertex edge length not implemented\n");
  return 0.12345;
}




void R3MeshVertex::
UpdateNormal(void)
{
  // Compute the surface normal at a vertex.  This feature should be implemented
  // second.  To do it, you must design a data structure that allows O(K)
  // access to faces attached to each vertex, where K is the number of faces attached
  // to the vertex.  Then, to compute the normal for a vertex,
  // you should take a weighted average of the normals for the attached faces, 
  // where the weights are determined by the areas of the faces.
  // Store the resulting normal in the "normal"  variable associated with the vertex. 
  // You can display the computed normals by hitting the 'N' key in meshview.

  // FILL IN IMPLEMENTATION HERE (THIS IS REQUIRED)
  // fprintf(stderr, "Update vertex normal not implemented\n");
}




void R3MeshVertex::
UpdateCurvature(void)
{
  // Compute an estimate of the Gauss curvature of the surface 
  // using a method based on the Gauss Bonet Theorem, which is described in 
  // [Akleman, 2006]. Store the result in the "curvature"  variable. 

  // FILL IN IMPLEMENTATION HERE
  // fprintf(stderr, "Update vertex curvature not implemented\n");
}





////////////////////////////////////////////////////////////
// MESH FACE MEMBER FUNCTIONS
////////////////////////////////////////////////////////////

R3MeshFace::
R3MeshFace(void)
  : vertices(),
    plane(0, 0, 0, 0),
    id(0)
{
}



R3MeshFace::
R3MeshFace(const R3MeshFace& face)
  : vertices(face.vertices),
    plane(face.plane),
    id(0)
{
}



R3MeshFace::
R3MeshFace(const vector<R3MeshVertex *>& vertices)
  : vertices(vertices),
    plane(0, 0, 0, 0),
    id(0)
{
  UpdatePlane();
}



double R3MeshFace::
AverageEdgeLength(void) const
{
  // Check number of vertices
  if (vertices.size() < 2) return 0;

  // Compute average edge length
  double sum = 0;
  R3Point *p1 = &(vertices.back()->position);
  for (unsigned int i = 0; i < vertices.size(); i++) {
    R3Point *p2 = &(vertices[i]->position);
    double edge_length = R3Distance(*p1, *p2);
    sum += edge_length;
    p1 = p2;
  }

  // Return the average length of edges attached to this face
  return sum / vertices.size();
}



double R3MeshFace::
Area(void) const
{
  // Check number of vertices
  if (vertices.size() < 3) return 0;

  // Compute area using Newell's method (assumes convex polygon)
  R3Vector sum = R3null_vector;
  const R3Point *p1 = &(vertices.back()->position);
  for (unsigned int i = 0; i < vertices.size(); i++) {
    const R3Point *p2 = &(vertices[i]->position);
    sum += p2->Vector() % p1->Vector();
    p1 = p2;
  }

  // Return area
  return 0.5 * sum.Length();
}



void R3MeshFace::
UpdatePlane(void)
{
  // Check number of vertices
  int nvertices = vertices.size();
  if (nvertices < 3) { 
    plane = R3null_plane; 
    return; 
  }

  // Compute centroid
  R3Point centroid = R3zero_point;
  for (int i = 0; i < nvertices; i++) 
    centroid += vertices[i]->position;
  centroid /= nvertices;
  
  // Compute best normal for counter-clockwise array of vertices using newell's method
  R3Vector normal = R3zero_vector;
  const R3Point *p1 = &(vertices[nvertices-1]->position);
  for (int i = 0; i < nvertices; i++) {
    const R3Point *p2 = &(vertices[i]->position);
    normal[0] += (p1->Y() - p2->Y()) * (p1->Z() + p2->Z());
    normal[1] += (p1->Z() - p2->Z()) * (p1->X() + p2->X());
    normal[2] += (p1->X() - p2->X()) * (p1->Y() + p2->Y());
    p1 = p2;
  }
  
  // Normalize normal vector
  normal.Normalize();
  
  // Update face plane
  plane.Reset(centroid, normal);
}



