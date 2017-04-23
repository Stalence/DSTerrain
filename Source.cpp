//2184 Fwteinh Bekiaroudh, 2167 Nikolaos Karalias
//Procedural terrain

#include <stdlib.h>
#include <glut.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <random>
#include <time.h> 

const int bounds = 513;

float left = -(bounds - 1) / 2.0 + 1000;
float right = (bounds - 1) / 2.0 - 1000;
float bottom = -1 * (bounds - 1) - 800;
float top = 1 * (bounds - 1) + 700;
float near = (bounds - 1)+1500 ;
float far = -(bounds - 1)-5500 ;

GLint window; 

float recupower = 0;

int mode = 0;
int style = 0;

float terrain[bounds][bounds] ;
float tempo[bounds - 1][bounds - 1];
float height =1;

int smoothing = 1;
int elevationobs = 1;
int topocheck = 1;
int diamondsquare = 1;

float toporadius = 150;


std::default_random_engine generator;
std::mt19937 seeder;


std::uniform_real_distribution<float> seeds(0.0,0.45);
float genseed=0;


int iterationcount = 0;
float obsx = 0;
float obsy = -700.5;
float obsz = 0;
float pi = 3.141592;

float viewangle =45;


float Roughness = 10;
float Roughness2 = 1;

int limit = 1;

//int readbmp()
//{
//		// Data read from the header of the BMP file
//	unsigned char header[54]; // Each BMP file begins by a 54-bytes header
//	unsigned int dataPos;     // Position in the file where the actual data begins
//	unsigned int width, height;
//	unsigned int imageSize;   // = width*height*3
//	// Actual RGB data
//	unsigned char * data;
//
//	//FILE * file = fopen(imagepath, "rb");
//	if (!file)         
//	{ 
//		printf("Image could not be opened\n"); 
//		return 0; 
//	}
//
//	if (fread(header, 1, 54, file) != 54){ // If not 54 bytes read : problem
//		printf("Not a correct BMP file\n");
//		return false;
//	}
//
//
//		// Read ints from the byte array
//	dataPos = *(int*)&(header[0x0A]);
//	imageSize = *(int*)&(header[0x22]);
//	width = *(int*)&(header[0x12]);
//	height = *(int*)&(header[0x16]);
//
//	if (imageSize == 0)    imageSize = width*height * 3; // 3 : one byte for each Red, Green and Blue component
//	if (dataPos == 0)      dataPos = 54;
//
//}


void observer(int key, int x, int y) //move around on scene with keyboard 
{
	if (key == GLUT_KEY_RIGHT)
	{
		viewangle += 1.5;
	}
	if (key == GLUT_KEY_LEFT)
	{
		viewangle -= 1.5;
	}
	if (key == GLUT_KEY_UP)
	{
		if (elevationobs)
		{ 
		height += 30.1;
		}
		else
		{
			obsy -= 30;
		}
	}
	if (key == GLUT_KEY_DOWN)
	{
		if (elevationobs)
		{
			height -= 30.1;
		}
		else
		{
			if (obsy<-100)
			{
			obsy += 30;
			}
		}
	}
	if (key == GLUT_KEY_PAGE_DOWN)
	{
		//if (Roughness > 0.012) 
		{
			Roughness2 -= 0.06;
		}
		//else
		/*{
			Roughness -= 0.0001;
		}*/
		std::cout << "Roughness : " << Roughness2 << "\n";

	}
	if (key == GLUT_KEY_PAGE_UP)
	{
		
		/*if(Roughness<1)*/ Roughness2 += 0.06;
		//height += 5;
		std::cout << "Roughness : " << Roughness2 << "\n";


	}


	glutPostRedisplay();
}

void vecdivide(GLfloat vec1[3], float divisor)
{
	vec1[0] = vec1[0] / divisor;
	vec1[1] = vec1[1] / divisor;
	vec1[2] = vec1[2] / divisor;
}

float vecdot(GLfloat vec1[3], GLfloat vec2[3])
{
	return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

void vecscale(GLfloat vec1[3], float scalefactor)
{

	vec1[0] = vec1[0] * scalefactor;
	vec1[1] = vec1[1] * scalefactor;
	vec1[2] = vec1[2] * scalefactor;
}

void vecsubtract(GLfloat vec1[3], GLfloat vec2[3], GLfloat vec3[3])
{
	vec3[0] = vec2[0] - vec1[0];
	vec3[1] = vec2[1] - vec1[1];
	vec3[2] = vec2[2] - vec1[2];
}

void normalize(float vector[3]) //function for vector normalization
{
	double norm = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
	vector[0] /= norm;
	vector[1] /= norm;
	vector[2] /= norm;

}

void flipnormal(float vector[3])  //flips a vector/normal on the opposite direction
{
	vector[0] = -vector[0];
	vector[1] = -vector[1];
	vector[2] = -vector[2];


}

void  vecavg(GLfloat vec1[3], GLfloat vec2[3], GLfloat vec3[3], GLfloat normed[3], int normflag)
{
	normed[0] = vec1[0] / 2.0 + vec2[0] / 2.0 + vec3[0] / 2.0;
	normed[1] = vec1[1] / 2.0 + vec2[1] / 2.0 + vec3[1] / 2.0;
	normed[2] = vec1[2] / 2.0 + vec2[2] / 2.0 + vec3[2] / 2.0;


	if (normflag)normalize(normed);
}

void crossprod(GLfloat vect1[3], GLfloat vect2[3], GLfloat crossvect[3])  //cross product
{

	crossvect[0] = vect1[1] * vect2[2] - vect1[2] * vect2[1];
	crossvect[1] = vect1[2] * vect2[0] - vect1[0] * vect2[2];
	crossvect[2] = vect1[0] * vect2[1] - vect1[1] * vect2[0];

	normalize(crossvect);


}

void init()
{
	/* attributes */
	glClearColor(1.0, 1.0, 1.0, 0.0); /* white background */
	/* orthographic projection */

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(left, right, bottom, top, near, far);
	glMatrixMode(GL_MODELVIEW);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_FLAT);
	glEnable(GL_NORMALIZE);
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	glEnable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	

}

void movingavg()
{


	for (int i = 1; i < bounds-1; ++i)
	{
		for (int j = 1; j < bounds-1; ++j)
		{
			
			tempo[i-1][j-1] = (1.0 / 10.0)*(terrain[i][j] + terrain[i - 1][j] + terrain[i + 1][j] + terrain[i][j - 1] + terrain[i][j + 1] + terrain[i - 1][j - 1] + terrain[i + 1][j + 1] +terrain[i-1][j+1]+terrain[i+1][j-1]);
			
		}
	}

	for (int i = 1; i < bounds - 1; ++i)
	{
		for (int j = 1; j < bounds - 1; ++j)
		{
			terrain[i][j] = tempo[i - 1][j - 1];
		}
	}
	

}

void median()
{

	for (int i = 1; i < bounds - 1; ++i)
	{
		for (int j = 1; j < bounds - 1; ++j)
		{
			std::vector<float> mask;

			for (int k = -1; k < 2; ++k)
			{
				for (int l = -1; l < 2; ++l)
				{
					mask.push_back(terrain[i + k][j + l]);
				}
			}
			std::sort(mask.begin(), mask.end());
			terrain[i][j] = mask[3];

		}
	}

}



void clampintensity()
{
  float 	max = -999999999;
  float     min = 999999999;
  float range = 0;
	for (int k = 0; k < bounds - 1; ++k)
	{
		for (int j = 0; j < bounds - 1;++j)
		{
			if (terrain[k][j]>max)
			{
				max = terrain[k][j];
			}
			if (terrain[k][j] < min)
			{
				min = terrain[k][j];
			}
		}
	}

	range = max - min;
	for (int k = 0; k < bounds - 1; ++k)
	{
		for (int j = 0; j < bounds - 1; ++j)
		{
			terrain[k][j] = abs(terrain[k][j]);
			terrain[k][j] = terrain[k][j]/(range);
			/*if (terrain[k][j] < 0.1*(max/range))
			{
				terrain[k][j] =0;
			}*/
		}
	}
}

void generatev2()
{
   int	iterations = log2(bounds - 1);
   int squaresize = bounds - 1;
 // pow(2, recupower);

   for (int q = 0; q< limit;++q)
   { 

   while (squaresize >limit)
   {
	   
	   for (int i = 0; i <bounds-1; i=i+squaresize)
	   {
		   for (int j = 0; j < bounds - 1; j = j + squaresize)
		   {
			   int maxx = i + squaresize;
			   int maxy = j + squaresize;
			   int minx = i;
			   int miny = j;

			   int centerx = (i + i + squaresize) / 2;
			   int centery = (j + j + squaresize) / 2;

			   int rangex = squaresize;
			   int rangey = squaresize;
			   float randmean = (float)(rangex) / (float)(bounds - 1);
			   float reclevel = log2(bounds - 1) - log2(rangex) + 1;
			   float fracdim = log(4.0) / log(1.0 / 4.0);
			   float scalef = 0.1*(float)(rangex);

			   float kay = 1;
			   float mean = 0.00;
			   float expo = sqrt(pow(squaresize, 2 * (3 - (Roughness2))));
			   float std = kay*expo;// roughness*expo;//pow(rangex, 6-2*fracdim);
			   float unibound = 1.0;

			   if (topocheck)
			   {
				   float  topocircdist = sqrt((centerx - (bounds - 1) / 2)*(centerx - (bounds - 1) / 2) + (centery - (bounds - 1) / 2)*(centery - (bounds - 1) / 2));
				   if (topocircdist <= toporadius)
				   {
					   mean = 0;
					   std = kay*expo;
				   }
				   else
				   {
					   mean = 0;
					   std = kay*expo*0.75;

				   }
			   }

			   std::normal_distribution<float> distribution(mean, std);


			   float topocircfact = 1;
			   if (topocheck)
			   {
				   float  topocircdist = sqrt((centerx - (bounds - 1) / 2)*(centerx - (bounds - 1) / 2) + (centery - (bounds - 1) / 2)*(centery - (bounds - 1) / 2));
				   if (topocircdist <= toporadius)
				   {
					   topocircfact = topocircdist*0.1;
				   }
			   }

			   float randfact1 = distribution(generator)*Roughness*topocircfact;

			   terrain[centerx][centery] =(((terrain[maxx][maxy] + terrain[maxx][miny] + terrain[minx][maxy] + terrain[minx][miny]) / 4.0) + randfact1 / 1);

			 
		   }
	   }

	   for (int i = 0; i <bounds - 1; i = i + squaresize)
	   {
		   for (int j = 0; j < bounds - 1; j = j + squaresize)
		   {
			   int maxx = i + squaresize;
			   int maxy = j + squaresize;
			   int minx = i;
			   int miny = j;

			   int centerx = (i + i + squaresize) / 2;
			   int centery = (j + j + squaresize) / 2;

			   int rangex = squaresize;
			   int rangey = squaresize;
			   float randmean = (float)(rangex) / (float)(bounds - 1);
			   float reclevel = log2(bounds - 1) - log2(rangex) + 1;
			   float fracdim = log(4.0) / log(1.0 / 4.0);
			   float scalef = 0.1*(float)(rangex);
			   float smoothness =0.1;
			   float roughness = 1.0 / smoothness;
			   float kay =1;
			   float mean = 0.00;
			   float expo = sqrt(pow(squaresize, 2 * (3 - (Roughness2))));
			   float std = kay*expo;//pow(rangex, 6-2*fracdim);
			   float unibound = 10.0;
			   //std::normal_distribution<float> distribution(mean, std);
			   //distribution.reset();



			   float topocircfact =1;

			   if (topocheck)
			   {
				   float  topocircdist = sqrt((centerx - (bounds - 1) / 2)*(centerx - (bounds - 1) / 2) + (centery - (bounds - 1) / 2)*(centery - (bounds - 1) / 2));
				   if (topocircdist <= toporadius)
				   {
					   mean = 0;
					   std = kay*expo;
				   }
				   else
				   {
					   mean = 0;
					   std = kay*expo*0.75;

				   }
			   }

			   std::normal_distribution<float> distribution(mean, std);




			   if (topocheck)
			   {

				   float  topocircdist = sqrt((centerx - (bounds - 1) / 2)*(centerx - (bounds - 1) / 2) + (centery - (bounds - 1) / 2)*(centery - (bounds - 1) / 2));
				   if (topocircdist <= toporadius)
				   {
					   topocircfact = topocircdist*0.1;
				   }
			   }


			   float randfact1 = distribution(generator)*Roughness*topocircfact;
			   float randfact2 = distribution(generator)*Roughness*topocircfact;
			   float randfact3 = distribution(generator)*Roughness*topocircfact;
			   float randfact4 = distribution(generator)*Roughness*topocircfact;
			   float randfact5 = distribution(generator)*Roughness*topocircfact;


			   float diamond1 = 0;
			   float diamond2 = 0;
			   float diamond3 = 0;
			   float diamond4 = 0;
			   float divdiam = 2;
			   float dsq = 0;

			   if (diamondsquare)
			   {
				   divdiam = 3;
				   dsq = 1;
			   }
			  

			   if ((centerx + rangex < bounds) && (centery + rangey < bounds) && (centerx - rangex >= 0) && (centery - rangey >= 0))
			   {
				   if (diamondsquare)
				   {
					   diamond4 =  terrain[centerx + rangex][centery] + 50 * distribution(generator);
					   diamond3 =  terrain[centerx - rangex][centery] + 50 * distribution(generator);
					   diamond2 =  terrain[centerx][centery + rangey] + 50 * distribution(generator);
					   diamond1 =  terrain[centerx][centery - rangey] + 50 * distribution(generator);
					   divdiam = 4;
				   }
				   else
				   {

				   
				   diamond4 = 0;
				   diamond3 = 0;
				   diamond2 = 0;
				   diamond1 = 0;
				   divdiam = 2;

				   }
			   }


			   {
				   //  std::cout << "Centerx : " << centerx << " Centery : " << centery << " Miny : " << miny << " Maxy : " << maxy << " Minx : " << minx << " Maxx : " << maxx << "\n";
			   }

			   terrain[centerx][miny] = (((terrain[minx][miny] * 1 + terrain[maxx][miny] * 1 + terrain[centerx][centery] * dsq + diamond1 * 1) / divdiam)*1.0 + randfact2 * 1.0 / 1.0);
			   terrain[centerx][maxy] = (((terrain[minx][maxy] * 1 + terrain[maxx][maxy] * 1 + terrain[centerx][centery] * dsq + diamond2 * 1) / divdiam)*1.0 + randfact3 * 1.0 / 1.0);
			   terrain[minx][centery] = (((terrain[minx][miny] * 1 + terrain[minx][maxy] * 1 + terrain[centerx][centery] * dsq + diamond3 * 1) / divdiam)*1.0 + randfact4 * 1.0 / 1.0);
			   terrain[maxx][centery] = (((terrain[maxx][miny] * 1 + terrain[maxx][maxy] * 1 + terrain[centerx][centery] * dsq + diamond4 * 1) / divdiam)*1.0 + randfact5 * 1.0 / 1.0);



		   }
	   }



	   squaresize = squaresize / 2;
   }

   }

}

void generate(int maxx, int minx, int maxy, int miny)
{
	iterationcount += 1;
	//std::cout << " Iteration count : " << iterationcount << "\n";

	int rangex = maxx - minx;
	int rangey = maxy - miny;
	float randmean = (float)(rangex) / (float)(bounds - 1);
	float reclevel = log2(bounds-1)- log2(rangex)+1;
	float fracdim =log(4.0)/log(1.0/4.0);
	float scalef = 0.1*(float)(rangex);
	float roughness = 1;
	float kay = 0.04;
	float mean = 0.00;
	float expo = pow(2,reclevel*fracdim);
	float scaled = 1*((float)rangex); /*pow(2, -roughness*reclevel)*(float)(rangex)*0.1;*/
	float std = 1.0;//pow(rangex, 6-2*fracdim);
	float unibound =10.0;

	if (rangex == 0)
	{
		std = 0.1;
	}
	//std::cout << "STD : " << std << "\n";


	std::normal_distribution<float> distribution(mean,std);
	//std::uniform_real_distribution<float> distribution(-unibound, unibound);
	std::uniform_int_distribution<int> distribution2(-rangex, rangex);
	std::uniform_int_distribution<int> shuff(1,5);

	float randfact1 = distribution(generator)*roughness*rangex;
	float randfact2 = distribution(generator)*roughness*rangex;
	float randfact3 = distribution(generator)*roughness*rangex;
	float randfact4 = distribution(generator)*roughness*rangex;
	float randfact5 = distribution(generator)*roughness*rangex;




	int centerx = floor((float)(minx+maxx) / 2.0);
	int centery = floor((float)(miny+maxy) / 2.0);
	



	/*std::cout << " Rangex : " << rangex << "\n";
	std::cout << " Rangey : " << rangey << "\n";
	std::cout << " Centerx : " << centerx << "\n";
	std::cout << " Centery : " << centery << "\n";
	std::cout << " Maxx : " << maxx << "\n";
	std::cout << " Maxy : " << maxy << "\n";
	std::cout << " Minx : " << minx << "\n";
	std::cout << " Miny : " << miny << "\n";


*/
	


	if ((rangex > 1) && (rangey > 1))
	{
		

		/*	if ((maxx < 0 || maxx>128) || (maxy < 0 || maxy>128) || (minx < 0 || minx>128) || (miny < 0 || miny>128) || (centerx < 0 || centerx>128) || (centery < 0 || centery>128))
			{
			std::cout << " Out of bounds. " << "\n";
			}*/


		float diamond1 = 0;
		float diamond2 = 0;
		float diamond3 = 0;
		float diamond4 = 0;

		if ((centerx + rangex < bounds) && (centery + rangey < bounds)&&(centerx-rangex>=0)&&(centery-rangey>=0))
		{
			diamond4 = terrain[centerx + rangex][centery];
			diamond3 = terrain[centerx - rangex][centery];
			diamond2 = terrain[centerx][centery + rangey];
			diamond1 = terrain[centerx][centery - rangey];
		}

	
		//+ terrain[centerx][miny] + terrain[centerx][maxy] + terrain[minx][centery] + terrain[maxx][centery]
	
		terrain[centerx][centery] = terrain[centerx][centery] + (((terrain[maxx][maxy] + terrain[maxx][miny] + terrain[minx][maxy] + terrain[minx][miny]) / 4.0) + randfact1 / 1);


		terrain[centerx][miny] = (((terrain[minx][miny] * 1 + terrain[maxx][miny] * 1 /*+ terrain[centerx][centery] * 1 + diamond1*1*/) / 2)*1.0 +randfact2 * 1 / 1);
		terrain[centerx][maxy] = (((terrain[minx][maxy] * 1 + terrain[maxx][maxy] * 1 /*+ terrain[centerx][centery] * 1 + diamond2*1*/) / 2)*1.0 + randfact3 * 1 / 1);
		terrain[minx][centery] = (((terrain[minx][miny] * 1 + terrain[minx][maxy] * 1/* + terrain[centerx][centery] * 1 + diamond3*1*/) / 2)*1.0 + randfact4 * 1 / 1);
		terrain[maxx][centery] = (((terrain[maxx][miny] * 1 + terrain[maxx][maxy] * 1/* + terrain[centerx][centery] * 1 + diamond4*1*/) / 2)*1.0 + randfact5 * 1 / 1);

	


//		std::cout << "Terrain center miny: " << terrain[centerx][centery] << "\n";
	/*	std::cout << "Terrain centerx miny: " << terrain[centerx][miny] << "\n";
		std::cout << "Terrain centerx maxy: " << terrain[centerx][maxy] << "\n";
		std::cout << "Terrain centery maxy: " << terrain[centery][maxy] << "\n";
		std::cout << "Terrain centery miny: " << terrain[centery][miny] << "\n";*/
		style = shuff(generator);

			if (style==1)
			{ 
				generate(centerx, minx, maxy, centery);
				generate(centerx, minx, centery, miny);
				generate(maxx, centerx, maxy, centery);
				generate(maxx, centerx, centery, miny);
			}
			else if (style==2)
			{ 
				generate(centerx, minx, centery, miny);
				generate(centerx, minx, maxy, centery);
				generate(maxx, centerx, maxy, centery);
				generate(maxx, centerx, centery, miny);
			}
			else if (style == 3)
			{
				generate(centerx, minx, centery, miny);
				generate(centerx, minx, maxy, centery);
				generate(maxx, centerx, centery, miny);
				generate(maxx, centerx, maxy, centery);
			}
			else if (style == 4)
			{
				generate(maxx, centerx, maxy, centery);
				generate(maxx, centerx, centery, miny);
				generate(centerx, minx, maxy, centery);
				generate(centerx, minx, centery, miny);
			}
			else if (style == 5)
			{
				generate(maxx, centerx, maxy, centery);
				generate(centerx, minx, centery, miny);
				generate(maxx, centerx, centery, miny);
				generate(centerx, minx, maxy, centery);
			}
		
	}
	else
	{

			
		return;
	}



}

void stochastic()
{
	float mean = 0;
	float roughness = 0.05;
	float fracdim = log(4.0) / log(1.0 / 4.0);


	std::uniform_int_distribution<int> blocksize(1, log2(bounds - 1));
	std::uniform_int_distribution<int> points(1, bounds - 1);




	for (int i = 0; i < 50000; ++i)
	{
		int x = points(generator);
		int y = points(generator);

		float std = roughness*pow(2, -fracdim*i*0.000005);
	    std::normal_distribution<float> distribution(mean, std);

		
		for (int j = 0; j < log2(bounds - 1); ++j)
		{

			int sqsize = blocksize(generator);
			if ((x + sqsize < bounds - 1) && (x - sqsize>0) && (y - sqsize > 0) && (y + sqsize < bounds - 1))
			{
				terrain[x][y] = (terrain[x + sqsize][y + sqsize] + terrain[x - sqsize][y - sqsize] + terrain[x + sqsize][y - sqsize] + terrain[x - sqsize][y + sqsize]) / 4.0 + distribution(generator);

				int halfblock = floor((float)(sqsize / 2));
				int r = x + halfblock;
				int l = x - halfblock;
				int u = y - halfblock;
				int d = y + halfblock;


				terrain[r][y] = (terrain[r + halfblock][y + halfblock] + terrain[r + halfblock][y - halfblock] + terrain[r - halfblock][y + halfblock] + terrain[r - halfblock][y - halfblock])/4.0+distribution(generator) ;
				terrain[l][y] = (terrain[l + halfblock][y + halfblock] + terrain[l - halfblock][y + halfblock] + terrain[l + halfblock][y - halfblock] + terrain[l - halfblock][y - halfblock])/4.0+distribution(generator);
				terrain[x][u] = (terrain[x+halfblock][u+halfblock] + terrain[x-halfblock][u+halfblock] + terrain[x+halfblock][u-halfblock] + terrain[x-halfblock][u-halfblock])/4.0 +distribution(generator);
				terrain[x][d] = terrain[x + halfblock][d+halfblock] + terrain[x - halfblock][d+halfblock] + terrain[x + halfblock][d-halfblock] + terrain[x - halfblock][d-halfblock]/4.0+distribution(generator);
			}

		}

	}

}



void terraingen()
{	
	terrain[0][0] = 0.15;
	terrain[bounds-1][0] = 0.15;
	terrain[0][bounds-1] = 0.15;
	terrain[bounds-1][bounds-1] = 0.15;
	terrain[(bounds - 1) / 2][(bounds - 1) / 2] = 0.15;

	//stochastic();
	//generate(bounds-1,0,(bounds-1),0);
	generatev2();

	clampintensity();
	//median();
	if (smoothing)
	{
		movingavg();
	}
	/*movingavg();
	movingavg();*/
   
	GLfloat dif1[4] = { 0.3, 0.45, 0.0, 1.0 };
	
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, dif1);
	glPushMatrix();
	glRotatef(0, 1, 0, 0);

	for (int k = 3; k < bounds-3; ++k)
	{
		glBegin(GL_TRIANGLE_STRIP);

		for (int j = 3; j < bounds-3; ++j)
		{	
			

			float intensity = terrain[k][j];
			if (intensity > 1)
			{
			//	std::cout << " Coordx : " << terrain[k][j] << "\n";
			}

			glColor3f( terrain[k][j], terrain[k][j] , 0.0);
			GLfloat norm[3] = {};
			GLfloat norm1[3] = {};
			GLfloat norm2[3] = {};
			GLfloat norm3[3] = {};

			GLfloat finorm[3] = {};
			GLfloat finorm1[3] = {};
			GLfloat finorm2[3] = {};
			GLfloat finorm3[3] = {};

			GLfloat temp1[3] = {};
			GLfloat temp2[3] = {};


			GLfloat point[3] = {};
			point[0] = 1 * (k - (bounds - 1) / 2);
			point[1] = terrain[k][j] * height;
			point[2] = 1 * (j - (bounds - 1) / 2);
			//std::cout << "Point x: " << point[0] << " Point y: " << point[1] << " Point z: " << point[2] << "\n";
			
	     /*	vecdivide(point, sqrt(vecdot(point, point)));
			vecscale(point, terrain[k][j] * height);*/


			GLfloat point3[3] = {};
			point3[0] = (float)k - 1 - ((bounds - 1) / 2);
			point3[1] = terrain[k - 1][j] * height;
			point3[2] = (float)j - ((bounds - 1) / 2);
			if ((k == 200) && (j == 1))	std::cout << "Point3 x: " << point3[0] << " Point y: " << point3[1] << " Point z: " << point3[2] << " obsx : " << obsx << " obsy : " << obsy << " obsz " << obsz << "\n";

			//vecdivide(point3, sqrt(vecdot(point3, point3)));
			//vecscale(point3, terrain[k - 1][j] * height);


			GLfloat point2[3] = {};
			point2[0] = (float)k - ((bounds - 1) / 2);
			point2[1] = terrain[k][j - 1] * height;
			point2[2] = (float)j - 1 - ((bounds - 1) / 2);
			//std::cout << "Point2 x: " << point2[0] << " Point y: " << point2[1] << " Point z: " << point2[2] << "\n";

			/*vecdivide(point2, sqrt(vecdot(point2, point2)));
			vecscale(point2, terrain[k][j - 1] * height);*/



			GLfloat point1[3] = {};
			point1[0] = (float)k - 1 - ((bounds - 1) / 2);
			point1[1] = terrain[k - 1][j - 1] * height;
			point1[2] = (float)j - 1 - ((bounds - 1) / 2);
			//std::cout << "Point1 x: " << point1[0] << " Point y: " << point1[1] << " Point z: " << point1[2] << "\n";

		/*	vecdivide(point1, sqrt(vecdot(point1, point1)));
			vecscale(point1, terrain[k - 1][j - 1] * height);
*/



			////////////////////////////////////////////////////////


		//	float pi = 3.14159265f;
		//	float thetaDelta = ((2 * pi) / (bounds - 1));
		//	float phiDelta = ((pi) / (bounds - 5.0));
		//	float RadiusBase = ((bounds) / pi / 2.0);
		//	float theta = (point[2] * thetaDelta);
		//	float phi = (point[0] * phiDelta);

		//	//Limit the map to half a sphere
		//	if (theta > pi) { theta = theta - (pi); }

		//	if (theta < 0.0)    { theta = theta + (pi); }

		//	if (phi > 2 * pi) { phi = phi - (2 * pi); }

		//	if (phi < 0.0) { phi = phi + (2 * pi); }

		//	//GLfloat  coords2[3] = {};
		//	point[0] = (float)(((RadiusBase)*sin(theta) * cos(phi)) + bounds / 2.0);
		//	point[1] = (float)((RadiusBase)* sin(theta) * sin(phi));
		//	point[2] = (float)(((RadiusBase)* cos(theta)) + bounds / 2.0);



		////	float pi = 3.14159265f;
		//	 thetaDelta = ((2 * pi) / (bounds - 1));
		//	 phiDelta = ((pi) / (bounds - 5.0));
		//	 RadiusBase = ((bounds) / pi / 2.0);
		//	 theta = (point1[2] * thetaDelta);
		//	 phi = (point1[0] * phiDelta);

		//	//Limit the map to half a sphere
		//	if (theta > pi) { theta = theta - (pi); }

		//	if (theta < 0.0)    { theta = theta + (pi); }

		//	if (phi > 2 * pi) { phi = phi - (2 * pi); }

		//	if (phi < 0.0) { phi = phi + (2 * pi); }

		//	//GLfloat  coords21[3] = {};
		//	point1[0] = (float)(((RadiusBase)*sin(theta) * cos(phi)) + bounds / 2.0);
		//	point1[1] = (float)((RadiusBase)* sin(theta) * sin(phi));
		//	point1[2] = (float)(((RadiusBase)* cos(theta)) + bounds / 2.0);




		////	float pi = 3.14159265f;
		//	 thetaDelta = ((2 * pi) / (bounds - 1));
		//	 phiDelta = ((pi) / (bounds - 5.0));
		//	 RadiusBase = ((bounds) / pi / 2.0);
		//	 theta = (point2[2] * thetaDelta);
		//	 phi = (point2[0] * phiDelta);

		//	//Limit the map to half a sphere
		//	if (theta > pi) { theta = theta - (pi); }

		//	if (theta < 0.0)    { theta = theta + (pi); }

		//	if (phi > 2 * pi) { phi = phi - (2 * pi); }

		//	if (phi < 0.0) { phi = phi + (2 * pi); }

		//	//GLfloat  coords22[3] = {};
		//	point2[0] = (float)(((RadiusBase)*sin(theta) * cos(phi)) + bounds / 2.0);
		//	point2[1] = (float)((RadiusBase)* sin(theta) * sin(phi));
		//	point2[2] = (float)(((RadiusBase)* cos(theta)) + bounds / 2.0);




		////	float pi = 3.14159265f;
		//	 thetaDelta = ((2 * pi) / (bounds - 1));
		//	 phiDelta = ((pi) / (bounds - 5.0));
		//	 RadiusBase = ((bounds) / pi / 2.0);
		//	 theta = (point3[2] * thetaDelta);
		//	 phi = (point3[0] * phiDelta);

		//	//Limit the map to half a sphere
		//	if (theta > pi) { theta = theta - (pi); }

		//	if (theta < 0.0)    { theta = theta + (pi); }

		//	if (phi > 2 * pi) { phi = phi - (2 * pi); }

		//	if (phi < 0.0) { phi = phi + (2 * pi); }

		//	//GLfloat  coords23[3] = {};
		//	point3[0] = (float)(((RadiusBase)*sin(theta) * cos(phi)) + bounds / 2.0);
		//	point3[1] = (float)((RadiusBase)* sin(theta) * sin(phi));
		//	point3[2] = (float)(((RadiusBase)* cos(theta)) + bounds / 2.0);



			///////////////////////////////////////////////////////





			vecsubtract(point, point3, temp1);
			vecsubtract(point, point2, temp2);
			crossprod(temp1, temp2, norm);
			flipnormal(norm);


			vecsubtract(point1, point3, temp1);
			vecsubtract(point1, point2, temp2);
			crossprod(temp1, temp2, norm1);
			flipnormal(norm);



			vecsubtract(point2, point, temp1);
			vecsubtract(point2, point3, temp2);
			crossprod(temp1, temp2, norm2);
			flipnormal(norm);



			vecsubtract(point3, point, temp1);
			vecsubtract(point3, point2, temp2);
			crossprod(temp1, temp2, norm3);
			//flipnormal(norm);



			vecavg(norm, norm1, norm2, finorm, 1);
			vecavg(norm1, norm, norm2, finorm1, 1);
			vecavg(norm2, norm, norm1, finorm2, 1);
			vecavg(norm3, norm1, norm2, finorm3, 1);


			glVertex3fv(point);
			glNormal3fv(finorm);
			glVertex3fv(point3);
			glNormal3fv(finorm3);
			glVertex3fv(point2);
			glNormal3fv(finorm2);
			glVertex3fv(point1);
			glNormal3fv(finorm1);


		}
		glEnd();
	}


	glPopMatrix();

}

void Jerry(int state, int button,int x, int y)
{
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		memset(tempo, 0.00, sizeof(tempo[0][0]) * bounds-1 *bounds-1);
		memset(terrain, 0.00, sizeof(terrain[0][0]) * bounds *bounds);
		genseed = seeds(seeder);
		std::cout << genseed << "\n";
		glutPostRedisplay();

	}
}


void planetes()
{
	terraingen();
	/*glPushMatrix();
	glTranslatef(0.0, 60.0, 60.0);
	glRotatef(90.0, 1.0, 0.0, 0.0);
	terraingen();
	glPopMatrix();

	glPushMatrix();

	glRotatef(90.0, 0.0, 0.0, 1.0);
	glTranslatef(60.0, 60.0, 0.0);
	terraingen();
	glPopMatrix();


	glPushMatrix();

	glRotatef(90.0, 0.0, 0.0, 1.0);
	glTranslatef(60.0, -60.0, 0.0);
	terraingen();
	glPopMatrix();


	glPushMatrix();

	glRotatef(-90.0, 1.0, 0.0, 0.0);
	glTranslatef(0.0, 60.0, 60.0);
	terraingen();
	glPopMatrix();

	glPushMatrix();
	glTranslatef(0.0, 120.0, 0.0);
	glRotatef(90.0, 0.0, 1.0, 0.0);
	terraingen();
	glPopMatrix();
*/
}

void display()
{

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/*for (int i = 0; i < 2;++i)
	{ 
		std::uniform_real_distribution<float> distribution(0.0, i);
	std::cout << "Gener1 : " << distribution(generator)<<" " ; 
	}
	*/
	generator.seed(genseed*100);



	glPushMatrix();
	float camerax = 650 * cos(viewangle*pi / 180.0);
	float cameraz = 650 * sin(viewangle*pi / 180.0);
	gluLookAt(camerax, obsy, cameraz, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);        //set up viewing

	GLfloat pos[] = { 270, 2500, 270, 1 };//position vector, 4th slot = 1--> positional light, 4th slot = 0 -->directional light
	GLfloat dif[] = { 0.3, 0.85, 0.0, 1.0 }; //diffuse vector
	GLfloat spec[] = { 0.00, 0.00, 0.0, 1.0 }; // specular vector
	GLfloat amb[] = { 0.1, 0.1, 0.1, 1.0 };// //ambient vector
	glLightfv(GL_LIGHT0, GL_POSITION, pos);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
	glLightfv(GL_LIGHT0, GL_SPECULAR, spec);
	glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0.0);
	glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.000);
	glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.0000005);
	glEnable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, amb);

	glScalef(3.5, 0.5, 3.5);

	planetes();
	
	glPopMatrix();
	glutSwapBuffers();
}





/* main menu */
void menu(int Selection)
{
	/* check if the user chose exit */
	if (Selection == 0)
	{
		glutDestroyWindow(window);
		exit(0);
	}
	/* else choose the right curve */
	else if (Selection == 1)
	{
		smoothing = (smoothing + 1) % 2;
	}
	else if (Selection == 2)
	{
		topocheck = (topocheck + 1) % 2;
	}
	else if (Selection == 3)
	{
		elevationobs = (elevationobs + 1) % 2;
	}
	else if (Selection == 4)
	{
		diamondsquare = (diamondsquare + 1) % 2;
	}
	glutPostRedisplay(); /* the window needs redrawing */
}

/* menus */
void create_menu(void)
{
	int app_menu;
	app_menu = glutCreateMenu(menu);
	glutAddMenuEntry("Smoothing", 1);
	glutAddMenuEntry("Circle topography", 2);
	glutAddMenuEntry("Elevation/Observer", 3);
	glutAddMenuEntry("Diamond/Square", 4);
	glutAddMenuEntry("Quit", 0);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
}


int main(int argc, char** argv)
{

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(bounds+800, bounds+800);
	glutInitWindowPosition(100, 100);
	window = glutCreateWindow("Fractal Terrain - Project 2");


	glutMouseFunc(Jerry);
	glutSpecialFunc(observer);


	create_menu();
	glutDisplayFunc(display);

	generator.seed(time(NULL));
	genseed = seeds(generator);

	init();
	glutMainLoop();             /* loop forever & start displaying */
	return 1;


}