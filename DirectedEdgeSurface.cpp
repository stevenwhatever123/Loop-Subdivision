///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  DirectedEdgeSurface.cpp
//  ------------------------
//  
//  Base code for rendering assignments.
//
//  Minimalist (non-optimised) code for reading and 
//  rendering an object file
//  
//  We will make some hard assumptions about input file
//  quality. We will not check for manifoldness or 
//  normal direction, &c.  And if it doesn't work on 
//  all object files, that's fine.
//
//  While I could set it up to use QImage for textures,
//  I want this code to be reusable without Qt, so I 
//  shall make a hard assumption that textures are in 
//  ASCII PPM and use my own code to read them
//  
///////////////////////////////////////////////////

// include the header file
#include "DirectedEdgeSurface.h"

// include the C++ standard libraries we want
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>

// include the Cartesian 3- vector class
#include "Cartesian3.h"
#include "SphereVertices.h"

#define MAXIMUM_LINE_LENGTH 1024
#define PI 3.141592

// constructor will initialise to safe values
DirectedEdgeSurface::DirectedEdgeSurface()
    : centreOfGravity(0.0,0.0,0.0)
    { // DirectedEdgeSurface()
    // force arrays to size 0
    vertices.resize(0);
    normals.resize(0);
	firstDirectedEdge.resize(0);
	faceVertices.resize(0);
	otherHalf.resize(0);
    } // DirectedEdgeSurface()

// read routine returns true on success, failure otherwise
bool DirectedEdgeSurface::ReadObjectStream(std::istream &geometryStream)
    { // ReadObjectStream()
    
    // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];
    
    // the rest of this is a loop reading lines & adding them in appropriate places
    while (true)
        { // not eof
		// token for identifying meaning of line
		std::string token;

        // character to read
        geometryStream >> token;
        
        // check for eof() in case we've run out
        if (geometryStream.eof())
            break;

        // otherwise, switch on the token we read
		if (token == "#")
			{ // comment 
			// read and discard the line
			geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
            } // comment
		else if (token == "Vertex")
			{ // vertex
			// variables for the read
			unsigned int vertexID;
			geometryStream >> vertexID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (vertexID != vertices.size())
				{ // bad vertex ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad vertex ID				
			
			// read in the new vertex position
			Cartesian3 newVertex;
			geometryStream >> newVertex;
			
			// and add it to the vertices
			vertices.push_back(newVertex);
			} // vertex
		else if (token == "Normal")
			{ // normal
			// variables for the read
			unsigned int normalID;
			geometryStream >> normalID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (normalID != normals.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new normal
			Cartesian3 newNormal;
			geometryStream >> newNormal;
			
			// and add it to the vertices
			normals.push_back(newNormal);
			} // normal
		else if (token == "FirstDirectedEdge")
			{ // first directed edge
			// variables for the read
			unsigned int FDEID;
			geometryStream >> FDEID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (FDEID != firstDirectedEdge.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new FDE
			unsigned int newFDE;
			geometryStream >> newFDE;
			
			// and add it to the vertices
			firstDirectedEdge.push_back(newFDE);
			} // first directed edge
		else if (token == "Face")
			{ // face
			// variables for the read
			unsigned int faceID;
			geometryStream >> faceID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (faceID != faceVertices.size()/3)
				{ // bad face ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad face ID				
			
			// read in the new face vertex (3 times)
			unsigned int newFaceVertex;
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			} // face
		else if (token == "OtherHalf")
			{ // other half
			// variables for the read
			unsigned int otherHalfID;
			geometryStream >> otherHalfID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (otherHalfID != otherHalf.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new face vertex (3 times)
			unsigned int newOtherHalf;
			geometryStream >> newOtherHalf;
			otherHalf.push_back(newOtherHalf);
			} // other half
        } // not eof

    // compute centre of gravity
    // note that very large files may have numerical problems with this
    centreOfGravity = Cartesian3(0.0, 0.0, 0.0);

    // if there are any vertices at all
    if (vertices.size() != 0)
        { // non-empty vertex set
        // sum up all of the vertex positions
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            centreOfGravity = centreOfGravity + vertices[vertex];
        
        // and divide through by the number to get the average position
        // also known as the barycentre
        centreOfGravity = centreOfGravity / vertices.size();

        // start with 0 radius
        objectSize = 0.0;

        // now compute the largest distance from the origin to a vertex
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            { // per vertex
            // compute the distance from the barycentre
            float distance = (vertices[vertex] - centreOfGravity).length();         
            
            // now test for maximality
            if (distance > objectSize)
                objectSize = distance;
            } // per vertex
        } // non-empty vertex set

	//sortFaces();
	// Find out all edges
	findEdges();

    // return a success code
    return true;
    } // ReadObjectStream()

// write routine
void DirectedEdgeSurface::WriteObjectStream(std::ostream &geometryStream)
    { // WriteObjectStream()
	geometryStream << "#" << std::endl; 
	geometryStream << "# Created for Leeds COMP 5821M Autumn 2020" << std::endl; 
	geometryStream << "#" << std::endl; 
	geometryStream << "#" << std::endl; 
	geometryStream << "# Surface vertices=" << vertices.size() << " faces=" << faceVertices.size()/3 << std::endl; 
	geometryStream << "#" << std::endl; 

	// output the vertices
    for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        geometryStream << "Vertex " << vertex << " " << std::fixed << vertices[vertex] << std::endl;

    // and the normal vectors
    for (unsigned int normal = 0; normal < normals.size(); normal++)
        geometryStream << "Normal " << normal << " " << std::fixed << normals[normal] << std::endl;

	// and the first directed edges
    for (unsigned int vertex = 0; vertex < firstDirectedEdge.size(); vertex++)
        geometryStream << "FirstDirectedEdge " << vertex<< " " << std::fixed << firstDirectedEdge[vertex] << std::endl;

    // and the faces - increment is taken care of internally
    for (unsigned int face = 0; face < faceVertices.size(); )
        { // per face
        geometryStream << "Face " << face/3 << " ";
        
        // read in three vertices
        geometryStream << faceVertices[face++] << " ";
        geometryStream << faceVertices[face++] << " ";
        geometryStream << faceVertices[face++];
            
        geometryStream << std::endl;
        } // per face

	// and the other halves
	for (unsigned int dirEdge = 0; dirEdge < otherHalf.size(); dirEdge++)
		geometryStream << "OtherHalf " << dirEdge << " " << otherHalf[dirEdge] << std::endl;

	for(unsigned int edge = 0; edge < edges.size()/2; edge++)
		geometryStream << "Edge " << edge << " " << edges[edge*2] << " " << edges[edge*2+1] << std::endl;
    } // WriteObjectStream()

// routine to render
void DirectedEdgeSurface::Render(RenderParameters *renderParameters)
    { // Render()
    // Ideally, we would apply a global transformation to the object, but sadly that breaks down
    // when we want to scale things, as unless we normalise the normal vectors, we end up affecting
    // the illumination.  Known solutions include:
    // 1.   Normalising the normal vectors
    // 2.   Explicitly dividing the normal vectors by the scale to balance
    // 3.   Scaling only the vertex position (slower, but safer)
    // 4.   Not allowing spatial zoom (note: sniper scopes are a modified projection matrix)
    //
    // Inside a game engine, zoom usually doesn't apply. Normalisation of normal vectors is expensive,
    // so we will choose option 2.  

    // Scale defaults to the zoom setting
    float scale = renderParameters->zoomScale;
    scale /= objectSize;
        
    //  now scale everything
    glScalef(scale, scale, scale);

    // apply the translation to the centre of the object if requested
    glTranslatef(-centreOfGravity.x, -centreOfGravity.y, -centreOfGravity.z);

    // start rendering
    glBegin(GL_TRIANGLES);

	// set colour for pick render - ignored for regular render
	glColor3f(1.0, 1.0, 1.0);

    // loop through the faces
	for (unsigned int face = 0; face < faceVertices.size(); face +=3)
		{ // per face
		// if we want flat normals, compute them here
		if (renderParameters->useFlatNormals)
			{ // flat normals
			// find two vectors along edges of the triangle
			Cartesian3 pq = vertices[faceVertices[face+1]] - vertices[faceVertices[face]];
			Cartesian3 pr = vertices[faceVertices[face+2]] - vertices[faceVertices[face]];

			// take their cross product and normalise
			Cartesian3 faceNormal = pq.cross(pr).unit();

			// and use it to set the glNormal
			glNormal3f(faceNormal.x * scale, faceNormal.y * scale, faceNormal.z * scale);
			} // flat normals

		// we have made a HARD assumption that we have enough normals
		for (unsigned int vertex = face; vertex < face+3; vertex++)
			{ // per vertex
		
			// if we are using smooth normals
			if (!renderParameters->useFlatNormals)
				// set the normal vector
				glNormal3f
					(
					normals[faceVertices[vertex]].x * scale,
					normals[faceVertices[vertex]].y * scale,
					normals[faceVertices[vertex]].z * scale
					);
			
			// and set the vertex position
			glVertex3f
				(
				vertices[faceVertices[vertex]].x,
				vertices[faceVertices[vertex]].y,
				vertices[faceVertices[vertex]].z
				);

			} // per vertex

		} // per face

    // close off the triangles
    glEnd();
    
    // now we add a second loop to render the vertices if desired
    if (!renderParameters->showVertices)
    	return;

	glDisable(GL_LIGHTING);

	// loop through the vertices
	for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
		{ // per vertex
		// use modelview matrix (not most efficient solution, but quickest to code)
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glTranslatef(vertices[vertex].x, vertices[vertex].y, vertices[vertex].z);
		glScalef(0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize);
		renderTriangulatedSphere();
		glPopMatrix();
		} // per vertex 
    
    } // Render()

void DirectedEdgeSurface::findEdges()
{
	edges.clear();
	for(int i = 0; i < faceVertices.size(); i+=3)
	{
		for(int j = 0; j < 3; j++)
		{
			int index1 = i + ((j + 2) % 3);
			int index2 = i + ((j + 3) % 3);

			//std::cout << index1 << "\n";
			//std::cout << index2 << "\n\n";

			// Check if edge exist
			bool edgeExist = false;
			for(int k = 0; k < edges.size(); k+=2)
			{
				if(edges[k] == faceVertices[index1] && edges[k+1] == faceVertices[index2])
				{
					edgeExist = true;
					//std::cout << edgeExist << "\n";
					break;
				}
			}

			if(!edgeExist)
			{	
				edges.push_back(faceVertices[index1]);
				edges.push_back(faceVertices[index2]);
			}
		}
	}
}

void DirectedEdgeSurface::updateFirstDirectedEdge()
{
	firstDirectedEdge.clear();
	for(int i = 0; i < vertices.size(); i++)
	{
		unsigned int min = -1;
		for(int j = 0; j < edges.size(); j+=2)
		{
			if(i == edges[j])
			{
				firstDirectedEdge.push_back(j/2);
				break;
			}
		}
	}
}

void DirectedEdgeSurface::updateOtherHalf()
{
	// Initialise all to -1;
	otherHalf.resize(edges.size()/2, -1);

	for(int i = 0; i < otherHalf.size(); i++)
	{
		bool firstEdgeFound = false;
		for(int j = 0; j < edges.size()/2; j++)
		{
			if(edges[i*2] == edges[j*2+1] && edges[i*2+1] == edges[j*2]
			  && !firstEdgeFound)
			{
				firstEdgeFound = true;
				otherHalf[i] = j;
			}
		}
	}
}

float DirectedEdgeSurface::getAlpha(int n)
{
	float alpha;
	float center = 0.375f + (0.25f * cos(2 * PI / (float) n));
	alpha = (1 / (float) n) * (0.625f - (center * center));
	return alpha;
}

void DirectedEdgeSurface::subdivision()
{
 	verticesA = vertices;
 	normalsA = normals;
	faceVerticesA.clear();
	faceVerticesA.shrink_to_fit();
	//edgesA.clear();
	edgesA.resize(edges.size() * 4, -1);

	// Subdivision
	for(int face = 0; face < faceVertices.size()/3; face++)
	{
		int currentVertex0 = faceVertices[face*3];
		int currentVertex1 = faceVertices[face*3 + 1];
		int currentVertex2 = faceVertices[face*3 + 2];
		std::vector<unsigned int> newVertices;

		for(int i = face * 6; i < face * 6 + 6; i+=2)
		{
			Cartesian3 *vertex0 = new Cartesian3();
			Cartesian3 *vertex1 = new Cartesian3();
			*vertex0 = vertices[edges[i]];
			*vertex1 = vertices[edges[i+1]];

			Cartesian3 *normal0 = new Cartesian3();
			Cartesian3 *normal1 = new Cartesian3();
			*normal0 = normals[edges[i]];
			*normal1 = normals[edges[i+1]];

			// Get the other vertex from the current face
			//int currentFaceIndex = std::round(i/3);
			int currentFaceIndex = face;
			int thirdVertex = -1;
			for(int k = 0; k < 3; k++)
			{
				// std::cout << face << "\n";
				// std::cout << i << "\n";
				// std::cout << currentFaceIndex << "\n\n";
				if(faceVertices[currentFaceIndex * 3 + k] != edges[i] 
					&& faceVertices[currentFaceIndex * 3 + k] != edges[i+1])
				{
					thirdVertex = faceVertices[currentFaceIndex * 3 + k];
					//std::cout << thirdVertex << "\n";
					break;
				}
			}
			
			//std::cout << "Third: " << thirdVertex << "\n";

			Cartesian3 *vertex2 = new Cartesian3();
			Cartesian3 *normal2 = new Cartesian3();
			*vertex2 = vertices[thirdVertex];
 			*normal2 = normals[thirdVertex];

			// Find the other vertex that is adjacent to the current face
			int otherV = -1;
			int currentEdgeIndex = i/2;
			int otherHalfEdgeIndex = otherHalf[currentEdgeIndex];

			// Select the end vertex from the next edge
			if(otherHalfEdgeIndex % 3 == 2)
				otherV = edges[((otherHalfEdgeIndex - 2) * 2) + 1];
			else
				otherV = edges[((otherHalfEdgeIndex + 1) * 2) + 1];

			Cartesian3 *vertex3 = new Cartesian3();
 			Cartesian3 *normal3 = new Cartesian3();
			*vertex3 = vertices[otherV];
 			*normal3 = normals[otherV];

			// Create a new vertex in between the edges
			// using edge vertex computation
			// Cartesian3 newVertex = 0.375f * vertex0
			// 				 + 0.375f * vertex1
			// 				 + 0.125f * vertex2
			// 				 + 0.125f * vertex3;

			// Cartesian3 newNormal = 0.375f * normal0
			// 				 + 0.375f * normal1
			// 				 + 0.125f * normal2
			// 				 + 0.125f * normal3;

			Cartesian3 *newVertex = new Cartesian3();
			Cartesian3 *newNormal = new Cartesian3();

			*newVertex = 0.375f * (*vertex0)
							 + 0.375f * (*vertex1)
							 + 0.125f * (*vertex2)
							 + 0.125f * (*vertex3);

			*newNormal = 0.375f * (*normal0)
							+ 0.375f * (*normal1)
							+ 0.125f * (*normal2)
							+ 0.125f * (*normal3);

			// Check if this new vertex exist
			bool vertexExist = false;
			int newVertexId = -1;
			for(int k = 0; k < verticesA.size(); k++)
			{
				if(*newVertex == verticesA[k])
				{
					vertexExist = true;
					newVertexId = k;
					break;
				}
			}

			if(!vertexExist)
			{
				newVertexId = verticesA.size();

				verticesA.push_back(*newVertex);

 				normalsA.push_back(*newNormal);
			}

			newVertices.push_back(newVertexId);

			// Dellocate memory
			delete vertex0;
			delete normal0;
			delete vertex1;
			delete normal1;
			delete vertex2;
			delete normal2;
			delete vertex3;
			delete normal3;
			delete newVertex;
			delete newNormal;
		}

		faceVerticesA.insert(faceVerticesA.end(), {currentVertex0, newVertices[1], newVertices[0]});
 		faceVerticesA.insert(faceVerticesA.end(), {currentVertex1, newVertices[2], newVertices[1]});
		faceVerticesA.insert(faceVerticesA.end(), {currentVertex2, newVertices[0], newVertices[2]});
 		faceVerticesA.insert(faceVerticesA.end(), {newVertices[0], newVertices[1], newVertices[2]});

		newVertices.clear();
		newVertices.shrink_to_fit();
	}

	// std::vector<bool> edgeVisited;
	// edgeVisited.resize(otherHalf.size(), false);

	// for(int face = 0; face < faceVertices.size()/3; face++)
	// {
	// 	int currentVertex0 = faceVertices[face*3];
	// 	int currentVertex1 = faceVertices[face*3 + 1];
	// 	int currentVertex2 = faceVertices[face*3 + 2];
	// 	std::vector<unsigned int> newVertices;

	// 	// Loop through every edge of the face
	// 	for(int i = face * 6; i < face * 6 + 6; i+=2)
	// 	{
	// 		Cartesian3 vertex0 = vertices[edges[i]];
	// 		Cartesian3 vertex1 = vertices[edges[i+1]];

	// 		Cartesian3 normal0 = normals[edges[i]];
	// 		Cartesian3 normal1 = normals[edges[i+1]];

	// 		// Get the other vertex from the current face
	// 		//int currentFaceIndex = std::round(i/3);
	// 		int currentFaceIndex = face;
	// 		int thirdVertex = -1;
	// 		for(int k = 0; k < 3; k++)
	// 		{
	// 			// std::cout << face << "\n";
	// 			// std::cout << i << "\n";
	// 			// std::cout << currentFaceIndex << "\n\n";
	// 			if(faceVertices[currentFaceIndex * 3 + k] != edges[i] 
	// 				&& faceVertices[currentFaceIndex * 3 + k] != edges[i+1])
	// 			{
	// 				thirdVertex = faceVertices[currentFaceIndex * 3 + k];
	// 				//std::cout << thirdVertex << "\n";
	// 				break;
	// 			}
	// 		}
			
	// 		//std::cout << "Third: " << thirdVertex << "\n";

	// 		Cartesian3 vertex2 = vertices[thirdVertex];
 	// 		Cartesian3 normal2 = normals[thirdVertex];

	// 		// Find the other vertex that is adjacent to the current face
	// 		int otherV = -1;
	// 		int currentEdgeIndex = i/2;
	// 		int otherHalfEdgeIndex = otherHalf[currentEdgeIndex];

	// 		// Select the end vertex from the next edge
	// 		if(otherHalfEdgeIndex % 3 == 2)
	// 			otherV = edges[((otherHalfEdgeIndex - 2) * 2) + 1];
	// 		else
	// 			otherV = edges[((otherHalfEdgeIndex + 1) * 2) + 1];

	// 		Cartesian3 vertex3 = vertices[otherV];
 	// 		Cartesian3 normal3 = normals[otherV];

	// 		// Create a new vertex in between the edges
	// 		// using edge vertex computation
	// 		Cartesian3 newVertex = 0.375f * vertex0
	// 						 + 0.375f * vertex1
	// 						 + 0.125f * vertex2
	// 						 + 0.125f * vertex3;

	// 		Cartesian3 newNormal = 0.375f * normal0
	// 						 + 0.375f * normal1
	// 						 + 0.125f * normal2
	// 						 + 0.125f * normal3;

	// 		// Cartesian3 newVertex = 0.5f * vertex0
	// 		// 				 + 0.5f * vertex1;

	// 		// Cartesian3 newNormal = 0.5f * normal0
	// 		// 				 + 0.5f * normal1;

	// 		// Check if this new vertex exist
	// 		bool vertexExist = false;
	// 		int newVertexId = -1;
	// 		for(int k = 0; k < verticesA.size(); k++)
	// 		{
	// 			if(newVertex == verticesA[k])
	// 			{
	// 				vertexExist = true;
	// 				newVertexId = k;
	// 				break;
	// 			}
	// 		}

	// 		if(!vertexExist)
	// 		{
	// 			newVertexId = verticesA.size();

	// 			verticesA.push_back(newVertex);

 	// 			normalsA.push_back(newNormal);
	// 		}

	// 		newVertices.push_back(newVertexId);

	// 		if(!edgeVisited[i/2])
	// 		{
	// 			int firstHalfEdgeIndex = i / 2;
	// 			int secondHalfEdgeIndex = edgesA.size()/2 - (firstHalfEdgeIndex + 1);

	// 			int otherHalfIndex = otherHalf[i/2];
	// 			int firstOtherHalfEdgeIndex = (otherHalfIndex + 1) * 4 - (otherHalfIndex + 1);
	// 			int secondOtherHalfEdgeIndex = (otherHalfIndex + 1) * 4;

	// 			// First Half
	// 			edgesA[firstHalfEdgeIndex * 2] = edges[i];
	// 			edgesA[firstHalfEdgeIndex * 2 + 1] = newVertexId;

	// 			edgesA[(secondHalfEdgeIndex - 1) * 2] = newVertexId;
	// 			edgesA[(secondHalfEdgeIndex - 1) * 2 + 1] = edges[i+1];

	// 			// Other Half
	// 			edgesA[(firstOtherHalfEdgeIndex - 1) * 2] = newVertexId;
	// 			edgesA[(firstOtherHalfEdgeIndex - 1) * 2 + 1] = edges[i];

	// 			edgesA[(secondOtherHalfEdgeIndex - 1) * 2] = edges[i+1];
	// 			edgesA[(secondOtherHalfEdgeIndex - 1) * 2 + 1] = newVertexId;

	// 			edgeVisited[i/2] = true;
	// 			edgeVisited[otherHalf[i/2]] = true;
	// 		}
	// 	}

	// 	// //Find three empty space for edge
	// 	std::vector<unsigned int> space;
	// 	int co = 0;
	// 	while(space.size() < 3)
	// 	{
	// 		if(edgesA[co] == -1 && edgesA[co + 1] == -1)
	// 		{
	// 			space.push_back(co);
	// 		}
	// 		co+=2;
	// 	}

	// 	edgesA[space[0]] = newVertices[0];
	// 	edgesA[space[0] + 1] = newVertices[1];

	// 	edgesA[space[1]] = newVertices[1];
	// 	edgesA[space[1] + 1] = newVertices[2];

	// 	edgesA[space[2]] = newVertices[2];
	// 	edgesA[space[2] + 1] = newVertices[0];

	// 	// Complete face
	// 	for(int i = 0; i < edges.size()/2; i++)
	// 	{
	// 		if(edges[i*2 + 1])
	// 		{

	// 		}
	// 		for(int j = 0; j < edges.size()/2; j++)
	// 		{
	// 			if(edges[i*2 + 1])
	// 		}
	// 	}

	// 	newVertices.clear();
	// }

	// Update position of the original vertices
	//std::vector<std::vector<Cartesian3>> neighbourList;

	for(int vertex = 0; vertex < vertices.size(); vertex++)
	{
		// Find all neighbour vertices
		std::vector<Cartesian3> list;
		int halfEdgeIndex = firstDirectedEdge[vertex];
		//std::cout << halfEdgeIndex << "\n";
		//std::cout << otherHalf[halfEdgeIndex] << "\n";

		int n = 0;
		do
		{
			//list.push_back(vertices[edges[halfEdgeIndex * 2 + 1]]);
			list.push_back(vertices[edges[halfEdgeIndex * 2 + 1]]);

			if(otherHalf[halfEdgeIndex] % 3 == 2)
			{
				halfEdgeIndex = otherHalf[halfEdgeIndex] - 2;
			}
			else
			{
				halfEdgeIndex = otherHalf[halfEdgeIndex] + 1;
			}
			n++;

		} while (halfEdgeIndex != firstDirectedEdge[vertex]);

		float alpha = getAlpha(n);

		// Cartesian3 neighbourVerticesSum;
		// Cartesian3 neighbourNormalsSum;
		Cartesian3 *neighbourVerticesSum = new Cartesian3();
		Cartesian3 *neighbourNormalsSum = new Cartesian3();
		for(int i = 0; i < n; i++)
		{
			//std::cout << list[i] << "\n";
			*neighbourVerticesSum = *neighbourVerticesSum + list[i];
			*neighbourNormalsSum = *neighbourNormalsSum + list[i];
		}
		list.clear();
		list.shrink_to_fit();
		//std::cout << "\n";
		
		verticesA[vertex] = (vertices[vertex] * (1.0f - n * alpha)) + (alpha * *neighbourVerticesSum);
		// std::cout << n << "\n";
		// std::cout << alpha << "\n";
		// std::cout << neighbourVerticesSum << "\n";
		// std::cout << vertices[vertex] << "\n";
		// std::cout << verticesA[vertex] << "\n\n";
		normalsA[vertex] = (normalsA[vertex] * (1.0f - n * alpha)) + (alpha * *neighbourNormalsSum);

		delete neighbourVerticesSum;
		delete neighbourNormalsSum;
	}
	
	vertices = verticesA;
 	normals = normalsA;
	// Create new faces
 	faceVertices = faceVerticesA;

	verticesA.clear();
	normalsA.clear();
	faceVerticesA.clear();
}