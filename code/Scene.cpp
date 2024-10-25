#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"


#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
*/
void Scene::convertPPMToPNG(string ppmFileName)
{
	string command;

	// Change implementation if necessary.
	command = "./magick convert " + ppmFileName + " " + ppmFileName + ".png";
	int removeWarning = system(command.c_str());
}




































/*
	Transformations, clipping, culling, rasterization are done here.
*/



//helper classes
class Vec3i
{
	public:
		int x, y;
		double z, r, g, b;

	Vec3i()
	{
		this->x = 0;
		this->y = 0;
		this->z = 0.0;
		this->r = 0.0;
		this->g = 0.0;
		this->b = 0.0;
	}

	Vec3i(int x, int y, double z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
		this->r = 0.0;
		this->g = 0.0;
		this->b = 0.0;
	}

	Vec3i(int x, int y, double z, double r, double g, double b)
	{
		this->x = x;
		this->y = y;
		this->z = z;
		this->r = r;
		this->g = g;
		this->b = b;
	}

	Vec3i(const Vec3i &other)
	{
		this->x = other.x;
		this->y = other.y;
		this->z = other.z;
		this->r = other.r;
		this->g = other.g;
		this->b = other.b;
	}
	
	friend std::ostream &operator<<(std::ostream &os, const Vec3i &v)
	{
		os << std::fixed << std::setprecision(6) << "[" << v.x << ", " << v.y << ", " << v.z << "]";
		return os;
	}
	
};


class myVec4
{

	public:
    double x, y, z, t;
    double r, g, b;

	myVec4()
	{
		this->x = 0.0;
		this->y = 0.0;
		this->z = 0.0;
		this->t = 0.0;
		this->r = 0.0;
		this->g = 0.0;
		this->b = 0.0;
	}

	myVec4(double x, double y, double z, double t)
	{
		this->x = x;
		this->y = y;
		this->z = z;
		this->t = t;
		this->r = 0.0;
		this->g = 0.0;
		this->b = 0.0;
	}

	myVec4(double x, double y, double z, double t, double r, double g, double b)
	{
		this->x = x;
		this->y = y;
		this->z = z;
		this->t = t;
		this->r = r;
		this->g = g;
		this->b = b;
	}

	myVec4(const myVec4 &other)
	{
		this->x = other.x;
		this->y = other.y;
		this->z = other.z;
		this->t = other.t;
		this->r = other.r;
		this->g = other.g;
		this->b = other.b;
	}
	
	double getNthComponent(int n)
	{
		switch (n)
		{
		case 0:
		    return this->x;

		case 1:
		    return this->y;

		case 2:
		    return this->z;

		case 3:
		default:
		    return this->t;
		}
	}

	friend std::ostream &operator<<(std::ostream &os, const myVec4 &v)
	{
		os << std::fixed << std::setprecision(6) << "[" << v.x << ", " << v.y << ", " << v.z << ", " << v.t << "]";
		return os;
	}
	
};










// class helpers


Vec3 Vec4ToVec3(Vec4 v4) {
	Vec3 v3 = {v4.x, v4.y, v4.z, v4.colorId};
	return v3;
}


Vec3i Vec4ToVec3i(Vec4 v4, Scene &scene) {
	Color c = *scene.colorsOfVertices[v4.colorId - 1];
	Vec3i v3 = {(int)v4.x, (int)v4.y, v4.z, c.r, c.g, c.b};
	return v3;
}



myVec4 Vec4ToMyVec4(Vec4 v4, Scene &scene) {
	Color c = *scene.colorsOfVertices[v4.colorId - 1];
	myVec4 v = {v4.x, v4.y, v4.z, v4.t, c.r, c.g, c.b};
	return v;
		
}


Vec3i myVec4ToVec3i(myVec4 v4) {
	Vec3i v3 = {(int)v4.x, (int)v4.y, v4.z, v4.r, v4.g, v4.b};
	return v3;
		
}




























//helpers


//2 functions to copy values from Vec3 to Matrix4
void RowToMatrix(Vec3 &v, Matrix4 &M, int i) {
	M.values[i][0] = v.x;
    M.values[i][1] = v.y;
    M.values[i][2] = v.z;
}

void ColToMatrix(Vec3 &v, Matrix4 &M, int j) {
	M.values[0][j] = v.x;
    M.values[1][j] = v.y;
    M.values[2][j] = v.z;
}




// Multiply 4x4 matrix4 with myVec4
myVec4 MTimesMyVec4(Matrix4 m, myVec4 v)
{
    double values[4];
    double total;

    for (int i = 0; i < 4; i++)
    {
        total = 0;
        for (int j = 0; j < 4; j++)
        {
            total += m.values[i][j] * v.getNthComponent(j);
        }
        values[i] = total;
    }

    return myVec4(values[0], values[1], values[2], values[3], v.r, v.g, v.b);
}






//calculates the camera and projection transformations
Matrix4 calcProjT(Camera &cam) {
	//camera transform
    Matrix4 projT = getIdentityMatrix();
    auto temp = inverseVec3(cam.position);
    ColToMatrix(temp, projT, 3);
	
	Matrix4 tempT = getIdentityMatrix();
	RowToMatrix(cam.u, tempT, 0);
	RowToMatrix(cam.v, tempT, 1);
	RowToMatrix(cam.w, tempT, 2);
	
	projT = multiplyMatrixWithMatrix(tempT, projT);
	
	
	//projection transformation
	
	double rl = cam.right - cam.left;
	double tb = cam.top - cam.bottom;
	double fn = cam.far - cam.near;
	
	if (cam.projectionType) { // if perspective projection enabled
		
		tempT = getIdentityMatrix();
		
		tempT.values[0][0] = (2*cam.near)/rl;
		tempT.values[0][2] = (cam.right + cam.left)/rl;
		tempT.values[1][1] = (2*cam.near)/tb;
		tempT.values[1][2] = (cam.top + cam.bottom)/tb;
		tempT.values[2][2] = -(cam.far + cam.near)/fn;
		tempT.values[2][3] = ((-2)*cam.far*cam.near)/fn;
		tempT.values[3][2] = -1;
		tempT.values[3][3] = 0;
	}
	
	else { // only orthographic transform
		tempT = getIdentityMatrix();
		tempT.values[0][0] = 2/rl;
		tempT.values[0][3] = -(cam.right + cam.left)/rl;
		tempT.values[1][1] = 2/tb;
		tempT.values[1][3] = -(cam.top + cam.bottom)/tb;
		tempT.values[2][2] = (-2)/fn;
		tempT.values[2][3] = -(cam.far + cam.near)/fn;
	}
	
	projT = multiplyMatrixWithMatrix(tempT, projT);
	return projT;

}


//calculates viewport transformation
Matrix4 calcViewT(Camera &cam) {

	Matrix4 viewT = getIdentityMatrix();
	viewT.values[0][0] = ((double)cam.horRes)/2;
	viewT.values[0][3] = ((double)(cam.horRes - 1))/2;
	viewT.values[1][1] = ((double)cam.verRes)/2;
	viewT.values[1][3] = ((double)(cam.verRes - 1))/2;
	viewT.values[2][2] = 0.5;
	viewT.values[2][3] = 0.5;
	
	return viewT;
}

Matrix4 calcModelT(Mesh &mesh, Scene &scene) {

    Matrix4 modelT = getIdentityMatrix();
    
    for (int i = 0; i < mesh.numberOfTransformations; i++) {
            Matrix4 temp = getIdentityMatrix();
        
        if (mesh.transformationTypes[i] == 't') { // translation
            Translation tr = *(scene.translations[mesh.transformationIds[i]-1]);
            temp.values[0][3] = tr.tx;
            temp.values[1][3] = tr.ty;
            temp.values[2][3] = tr.tz;
            modelT = multiplyMatrixWithMatrix(temp, modelT);
        }
        
        
        
        else if (mesh.transformationTypes[i] == 'r') { // rotation

			Rotation rt = *(scene.rotations[mesh.transformationIds[i]-1]);
			double alpha = (rt.angle*M_PI)/180;
			
			
			Vec3 u = {rt.ux, rt.uy, rt.uz};
			Vec3 v;
			if (rt.uz <= rt.ux && rt.uz <= rt.uy) // z is smallest
				v = {-rt.uy, rt.ux, 0};
			else if (rt.uy <= rt.ux && rt.uy <= rt.uz) // y is smallest
				v = {-rt.uz, 0, rt.ux};
			else // x is smallest
				v = {0, -rt.uz, rt.uy};
			Vec3 w = crossProductVec3(u, v);
			
			normalizeVec3(u);
			normalizeVec3(v);
			normalizeVec3(w);

			RowToMatrix(u, temp, 0);
			RowToMatrix(v, temp, 1);
			RowToMatrix(w, temp, 2);
			
			modelT = multiplyMatrixWithMatrix(temp, modelT); //rotating basis
			
			temp = getIdentityMatrix();
			double c = cos(alpha);
			double s = sin(alpha);
			temp.values[1][1] = c;
			temp.values[1][2] = -s;
			temp.values[2][1] = s;
			temp.values[2][2] = c;
			
			modelT = multiplyMatrixWithMatrix(temp, modelT); // performing rotation
			
			temp = getIdentityMatrix();
			ColToMatrix(u, temp, 0);
			ColToMatrix(v, temp, 1);
			ColToMatrix(w, temp, 2);
			
			modelT = multiplyMatrixWithMatrix(temp, modelT); //inverting the basis
        }
        
        
        else { // scaling
            Scaling sc = *(scene.scalings[mesh.transformationIds[i]-1]);
            temp.values[0][0] = sc.sx;
            temp.values[1][1] = sc.sy;
            temp.values[2][2] = sc.sz;
            modelT = multiplyMatrixWithMatrix(temp, modelT);
        }
        
    }
    
    return modelT;
}



void persDiv(Vec4 &curr) {
	curr.x /= curr.t;
	curr.y /= curr.t;
	curr.z /= curr.t;
	curr.t = 1;
}

void myPersDiv(myVec4 &curr) {
	curr.x /= curr.t;
	curr.y /= curr.t;
	curr.z /= curr.t;
	curr.t = 1;
}



bool visible(double d, double num, double &te, double &tl) {

	if (d > 0) { //potentially entering
		double t = num/d;
		if (t > tl)
			return false;
		if (t > te)
			te = t;
	}
	
	else if (d < 0) { //potentially leaving
		double t = num/d;
		if (t < te)
			return false;
		if (t < tl)
			tl = t;
	}
	else if (num > 0) //line parallel to edge
		return false;
	
	return true;
	
}


bool clipHelper(Vec4 v0, Vec4 v1, vector<myVec4> &curr, Camera &cam, Scene &scene) {
	

	double te = 0, tl = 1;
	double dx = v1.x - v0.x, dy = v1.y - v0.y, dz = v1.z - v0.z;
	if (!visible(dx, -1 - v0.x, te, tl)) return false;
	if (!visible(-dx, v0.x - 1, te, tl)) return false;
	if (!visible(dy, -1 - v0.y, te, tl)) return false;
	if (!visible(-dy, v0.y - 1, te, tl)) return false;
	if (!visible(dz, -1 - v0.z, te, tl)) return false;
	if (!visible(-dz, v0.z - 1, te, tl)) return false;
	
	curr.clear();
	
	Color c0 = *scene.colorsOfVertices[v0.colorId - 1];
	Color c1 = *scene.colorsOfVertices[v1.colorId - 1];
	myVec4 myv0;
	myVec4 myv1;
	
	
	if (tl < 1) {
		myv1.x = v0.x + dx*tl; myv1.y = v0.y + dy*tl; myv1.z = v0.z + dz*tl;
		myv1.t = v1.t;
		
		myv1.r = (1 - tl)*c0.r + tl*c1.r;
		myv1.g = (1 - tl)*c0.g + tl*c1.g;
		myv1.b = (1 - tl)*c0.b + tl*c1.b;
	}
	
	else myv1 = Vec4ToMyVec4(v1, scene);
	
	if (te > 0) {
		myv0.x = v0.x + dx*te; myv0.y = v0.y + dy*te; myv0.z = v0.z + dz*te;
		myv0.t = v0.t;
		
		myv0.r = (1 - te)*c0.r + te*c1.r;
		myv0.g = (1 - te)*c0.g + te*c1.g;
		myv0.b = (1 - te)*c0.b + te*c1.b;
	}
	
	else myv0 = Vec4ToMyVec4(v0, scene);
	curr.push_back(myv0);
	curr.push_back(myv1);
	return true;
}



vector<vector<myVec4>> clipping(vector<Vec4> &face, Camera &cam, Scene &scene) {
	
	vector<myVec4> curr;
	vector<vector<myVec4>> result;
	if (clipHelper(face[0], face[1], curr, cam, scene)) result.push_back(curr);
	if (clipHelper(face[1], face[2], curr, cam, scene)) result.push_back(curr);
	if (clipHelper(face[2], face[0], curr, cam, scene)) result.push_back(curr);
	
	return result;
}





bool backface(vector<Vec4> &face) {
	
	Vec3 v0 = Vec4ToVec3(face[0]);
	Vec3 v1 = Vec4ToVec3(face[1]);
	Vec3 v2 = Vec4ToVec3(face[2]);
	Vec3 triN = subtractVec3(v1, v0);
	Vec3 tempV = subtractVec3(v2, v0);
	triN = crossProductVec3(triN, tempV); // normal of triangle
	
	Vec3 triMid = multiplyVec3WithScalar(v0, 1/3);
	tempV = multiplyVec3WithScalar(v1, 1/3);
	Vec3 tempV2 = multiplyVec3WithScalar(v2, 1/3);
	triMid = addVec3(addVec3(triMid, tempV), tempV2); //midpoint of triangle
	
	if (dotProductVec3({0, 0, 1}, triN) < 0) return true;
	else return false;
	
}











void drawVert(Vec3i &v0, Vec3i &v1, Scene &scene) {

	int x = v0.x;
	int dy = v1.y - v0.y;
	
	double z = v0.z;
	double dz = (v1.z - v0.z) / (double)dy;
	
	Color c = {v0.r, v0.g, v0.b};
	double dcr = (v1.r - v0.r) / (double)dy;
	double dcg = (v1.g - v0.g) / (double)dy;
	double dcb = (v1.b - v0.b) / (double)dy;
	
	for (int y = v0.y; y <= v1.y; y++) {
	
		if (z < scene.depth[x][y]) {
			scene.assignColorToPixel(x, y, c);
			scene.depth[x][y] = z;
		}
		z += dz; // increment depth
		c.r += dcr; c.g += dcg; c.b += dcb; // increment color
	}
	
}


void drawLower(Vec3i &v0, Vec3i &v1, Scene &scene) {

	int dx = v1.x - v0.x;
	
	int y = v0.y;
	int dy = v0.y - v1.y;
	int yi = 1;
	
	if (dy > 0) { // y0 > y1
		yi = -1;
		dy = -dy;
	}
	
	int d = 2*dy + dx;
	
	double z = v0.z;
	double dz = (v1.z - v0.z) / (double)dx;
	
	Color c = {v0.r, v0.g, v0.b};
	double dcr = (v1.r - v0.r) / (double)dx;
	double dcg = (v1.g - v0.g) / (double)dx;
	double dcb = (v1.b - v0.b) / (double)dx;
	
	for (int x = v0.x; x <= v1.x; x++) {
	
		if (z < scene.depth[x][y]) {
			scene.assignColorToPixel(x, y, c);
			scene.depth[x][y] = z;
		}
		
		if (d < 0) { //choose NE
			y += yi;
			d += 2*(dy + dx);
		}
		
		else //choose E
			d+= 2*dy;
		
		z += dz; // increment depth
		c.r += dcr; c.g += dcg; c.b += dcb; // increment color
	}
	
	
}


void drawUpper(Vec3i &v0, Vec3i &v1, Scene &scene) {

	int dy = v1.y - v0.y;
	
	int x = v0.x;
	int dx = v0.x - v1.x;
	int xi = 1;
	
	if (dx > 0) { // x0 > x1
		xi = -1;
		dx = -dx;
	}
	
	int d = 2*dx + dy;
	
	double z = v0.z;
	double dz = (v1.z - v0.z) / (double)dy;
	
	Color c = {v0.r, v0.g, v0.b};
	double dcr = (v1.r - v0.r) / (double)dy;
	double dcg = (v1.g - v0.g) / (double)dy;
	double dcb = (v1.b - v0.b) / (double)dy;
	
	for (int y = v0.y; y <= v1.y; y++) {
	
		if (z < scene.depth[x][y]) {
			scene.assignColorToPixel(x, y, c);
			scene.depth[x][y] = z;
		}
		
		if (d < 0) { //choose NE
			x += xi;
			d += 2*(dx + dy);
		}
		
		else //choose E
			d+= 2*dx;
		z += dz; // increment depth
		c.r += dcr; c.g += dcg; c.b += dcb; // increment color
	}
	
}




void rasterWire(vector<Vec3i> &line, Camera &cam, Scene &scene) {

	double dx = line[1].x - line[0].x;
	double dy = line[1].y - line[0].y;
	
	if (dx == 0) {// vertical line
		if (dy < 0) drawVert(line[1], line[0], scene);
		else drawVert(line[0], line[1], scene);
	}
	
	
	else if (abs(dy) < abs(dx)) { // gradient between abs 0 and 1
		if (dx < 0) drawLower(line[1], line[0], scene);
		else drawLower(line[0], line[1], scene);
	}
	
	else { // gradient between abs 1 and infinity
		if (dy < 0) drawUpper(line[1], line[0], scene);
		else drawUpper(line[0], line[1], scene);
	}
		
}




















void rasterSolid(vector<Vec3i> &face, Camera &cam, Scene &scene) {

	Vec3i v0 = face[0], v1 = face[1], v2 = face[2];
	
	int xmax = max({v0.x, v1.x, v2.x});
	int ymax = max({v0.y, v1.y, v2.y});
	
	double f01 = v2.x*(v0.y - v1.y) + v2.y*(v1.x - v0.x) + v0.x*v1.y - v0.y*v1.x;
	double f12 = v0.x*(v1.y - v2.y) + v0.y*(v2.x - v1.x) + v1.x*v2.y - v1.y*v2.x;
	double f20 = v1.x*(v2.y - v0.y) + v1.y*(v0.x - v2.x) + v2.x*v0.y - v2.y*v0.x;
	
	for (int y = min({v0.y, v1.y, v2.y}); y <= ymax; y++) {
	for (int x = min({v0.x, v1.x, v2.x}); x <= xmax; x++) {
	
		double alpha = ((double)(x*(v1.y - v2.y) + y*(v2.x - v1.x) + v1.x*v2.y - v1.y*v2.x)) / f12;
		double beta = ((double)(x*(v2.y - v0.y) + y*(v0.x - v2.x) + v2.x*v0.y - v2.y*v0.x)) / f20;
		double gamma = ((double)(x*(v0.y - v1.y) + y*(v1.x - v0.x) + v0.x*v1.y - v0.y*v1.x)) / f01;
		
		if (alpha >= 0 && beta >= 0 && gamma >= 0) {
			
			double z = alpha*v0.z + beta*v1.z + gamma*v2.z;
			
			if (x >= 0 && x < cam.horRes && y >= 0 && y < cam.verRes && z < scene.depth[x][y]) {
				Color c;
				c.r = alpha*v0.r + beta*v1.r + gamma*v2.r;
				c.g = alpha*v0.g + beta*v1.g + gamma*v2.g;
				c.b = alpha*v0.b + beta*v1.b + gamma*v2.b;
				scene.assignColorToPixel(x, y, c);
				scene.depth[x][y] = z;
			}
		}
		
		
	}
	}

}











void mainHelper(Mesh &mesh, Camera &cam, Scene &scene, Matrix4 &projT, Matrix4 &viewT) {
	
	for (auto tri : mesh.triangles) { // loop over all triangles
	
		Vec3 v[3];
		v[0] = *scene.vertices[tri.vertexIds[0]-1];
		v[1] = *scene.vertices[tri.vertexIds[1]-1];
		v[2] = *scene.vertices[tri.vertexIds[2]-1];
		
		
		vector<Vec4> face;
		for (int i = 0; i < 3; i++){ //iterating over three vertices of triangle
		
			Vec4 curr = {v[i].x, v[i].y, v[i].z, 1, v[i].colorId};
			curr = multiplyMatrixWithVec4(projT, curr); //Vec4 after model, cam, and proj transforms
			if (cam.projectionType) //if perspective projection, do perspective divide
					persDiv(curr);
			face.push_back(curr); //face is triangle after perspective divide, ready for clipping
			
		}
		
		
		if (scene.cullingEnabled && backface(face))//check for backface
			continue;

		
		if (mesh.type) { // IF MESH IS SOLID
			vector<Vec3i> triangle;
			for (int i = 0; i < 3; i++){
				Vec3i currI = Vec4ToVec3i(multiplyMatrixWithVec4(viewT, face[i]), scene);
					triangle.push_back(currI);
				}
			rasterSolid(triangle, cam, scene);
		}
		
		
		else { // MESH IS WIREFRAME
			vector<vector<myVec4>> lines = clipping(face, cam, scene);
			
			for (int i = 0; i < lines.size(); i++) { //for each line (0 to 3 lines)
			vector<Vec3i> line;
			
				for (int j = 0; j < 2; j++) { //for each of the two vertexes in line
					Vec3i currI = myVec4ToVec3i(MTimesMyVec4(viewT, lines[i][j]));
					line.push_back(currI);
				}
				
			rasterWire(line, cam, scene);
			
			}
			
		}
	
	} //end of main for loop
	
}

















//main HW function
void Scene::forwardRenderingPipeline(Camera *camera)
{

    //calculating camera and projection transformations
    Matrix4 projT = calcProjT(*camera);
    
    //calculating viewport transform
    Matrix4 viewT = calcViewT(*camera);
    
    
    for (auto mesh: meshes) { //iterating over meshes, 1st loop
    
	    Matrix4 modelT = calcModelT(*mesh, *this);
	    modelT = multiplyMatrixWithMatrix(projT, modelT);
	    
	    mainHelper(*mesh, *camera, *this, modelT, viewT);
	    
	    //SHOULD BE DONE AT THIS POINT
    }
    
}



