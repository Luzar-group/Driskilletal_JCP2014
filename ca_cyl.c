// Calculate the contact angle of a cylindrical drop oriented along
// the x axis by using the circular fit by Taubin.
//
// Ver. 0.4 - 20/09/2013


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define X 10
#define Y 11
#define Z 12

typedef struct 
{
  double coord[3];
  int molecule;
  int atom_type;
  int index;
} DATA;

typedef struct
{
  int type;
  int atom[3];
} LIST;


void circlefit(double pts[][2], int npoints, double *radius, double *center, double *detcoeff)
{

  /*
    Find the best circle that fit a given set of points describing an arc 
    using Taubin's method.

    G. Taubin, IEEE Trans. PAMI, Vol.13, pp.1115-1138, 1991 

    Translation from Matlab code written by Nikolai Chernov
    http://www.math.uab.edu/~chernov/cl/
   */
  
  int point, iter;
  double centroid[2];
  double Mxx, Myy, Mzz, Mxy, Mxz, Myz;
  double Mz, Cov_xy;
  double A0, A1, A2, A3, A22, A33;
  double norm_pt[3];
  double new_pt[2], old_pt[2], epsilon, maxIter, Dy;
  double Det;
  double diff[2];

  // Calculate data set centroid

  for (point = 0; point < npoints; point++) 
    {
      centroid[0] += pts[point][0];
      centroid[1] += pts[point][1];
    }

  centroid[0] /= (double)npoints;
  centroid[1] /= (double)npoints;

  // Calculate normalized moments

  Mxx = 0.0;
  Myy = 0.0;
  Mzz = 0.0;
  Mxy = 0.0;
  Mxz = 0.0;
  Myz = 0.0;

  for (point = 0; point < npoints; point++) 
    {
      norm_pt[0] = pts[point][0] - centroid[0];
      norm_pt[1] = pts[point][1] - centroid[1];
      norm_pt[2] = norm_pt[0] * norm_pt[0] + norm_pt[1] * norm_pt[1];
      
      Mxx += norm_pt[0] * norm_pt[0];
      Myy += norm_pt[1] * norm_pt[1];
      Mzz += norm_pt[2] * norm_pt[2];
      Mxy += norm_pt[0] * norm_pt[1];
      Mxz += norm_pt[0] * norm_pt[2];
      Myz += norm_pt[1] * norm_pt[2];
    }

  Mxx /= (double)npoints;
  Myy /= (double)npoints;
  Mzz /= (double)npoints;
  Mxy /= (double)npoints;
  Mxz /= (double)npoints;
  Myz /= (double)npoints;

  // Calculate polynomial coefficients

  Mz = Mxx + Myy;
  Cov_xy = Mxx * Myy - Mxy * Mxy;

  A0 = Mxz * Mxz * Myy + Myz * Myz * Mxx - Mzz * Cov_xy -
       2.0 * Mxz * Myz * Mxy + Mz * Mz * Cov_xy;
  A1 = Mzz * Mz + 4 * Cov_xy * Mz - Mxz * Mxz - Myz * Myz - Mz * Mz * Mz;
  A2 = -3 * Mz * Mz - Mzz;
  A3 = 4 * Mz;
  
  A22 = 2.0 * A2;
  A33 = 3.0 * A3;

  // Find zeroes using Newton-Rhapson method

  new_pt[0] = 0.0;
  new_pt[1] = 1e+20;
  
  epsilon = 1e-12;
  maxIter = 20;
  iter = 0;

  do
    {
      old_pt[1] = new_pt[1];
      new_pt[1] = A0 + new_pt[0] * (A1 + new_pt[0] * (A2 + new_pt[0] * A3));

      if (fabs(new_pt[1]) > fabs(old_pt[1]))
	{
	  new_pt[0] = 0.0;
	  printf("WARNING: Newton-Taubin goes in wrong direction.\n");
	  break;
	}

      Dy = A1 + new_pt[0] * (A22 + new_pt[0] * A33);
      
      old_pt[0] = new_pt[0];
      new_pt[0] = old_pt[0] - new_pt[1] / Dy;

      if ( new_pt[0] < 0.0) 
	{
	  new_pt[0] = 0.0;
	  //printf("WARNING: Newton-Taubin negative root.\n");
	}
    } 
  while ((fabs((new_pt[0] - old_pt[0]) / new_pt[0]) > epsilon) && (++iter < maxIter)); 
  
  if (iter == maxIter)
    {
      new_pt[0] = 0.0;
      printf("Newton-Taubin don't converge.\n");
    }

  // Compute circle center and radius

  Det = new_pt[0] * new_pt[0] - new_pt[0] * Mz + Cov_xy;

  center[0] = (Mxz * (Myy - new_pt[0]) - Myz * Mxy) / (2.0 * Det);
  center[1] = (Myz * (Mxx - new_pt[0]) - Mxz * Mxy) / (2.0 * Det);

  *radius = sqrt(((center[0] * center[0]) + (center [1] * center[1])) + Mz);

  center[0] += centroid[0];
  center[1] += centroid[1];

  // Calculate coefficient of determination

  *detcoeff = 0.0;

  for (point = 0; point < npoints; point++) 
    {

      diff[0] = pts[point][0] - center[0];
      diff[1] = pts[point][1] - center[1];
      *detcoeff += fabs((sqrt(diff[0] * diff[0] + diff[1] * diff[1])) - *radius);
    }

  *detcoeff /= npoints;

}


double com_pbc_Y(DATA *molset, int natoms, int selspe, double box_size[3])
{
  // Calculate center of mass of periodic system by using circular 
  // projection of cartesian coordinates.

  int i;
  int counter;
  double com_y;
  double alpha, beta, gamma;
  double twopi;
  const double pi=3.1415926535897932384626433832795;

  twopi = 2.0 * pi;

  beta = 0.0;
  gamma = 0.0;

  counter = 0;

  for (i = 0; i < natoms; i++)
    {
      if (molset[i].atom_type == selspe)
	{
	  alpha = (molset[i].coord[1] / box_size[1]) * twopi;
	  
	  beta += (box_size[1] / twopi) * cos(alpha);
	  gamma += (box_size[1] / twopi) * sin(alpha);

	  counter++;
	}
    }

  beta /= counter;
  beta *= -1.0;

  gamma /= counter;
  gamma *= -1.0;

  alpha = atan2(gamma, beta) + pi;

  com_y = box_size[1] * (alpha / twopi);

  return com_y;
}

double contact_angle(double center[], double radius, double line[])
{
  double dist;
  double a, b, c, delta;
  double secant[2][2], vec1[2], vec2[2];
  double dotprod, modprod;
  double angle;
  const double pi=3.1415926535897932384626433832795;

  // line[0] = Line slope
  // line[1] = Line intercept

  // Check that line intersect the circle

  dist = fabs(center[1] - (line[0] * center[0]) - line[1]);
  dist /= sqrt(1.0 + (line[0] * line[0]));
  
  if (dist < radius) 
    {
      
      a = (1.0 + (line[0] * line[0]));
      
      b = 2.0 * ((line[0] * line[1]) - center[0] + (line[0] * center[1]));
      
      c = (center[0] * center[0]) + (center[1] * center[1])
	+ (line[1] * line[1]) - (radius * radius) - (2.0 * line[1] * center[1]);
      
      delta = (b * b) - (4.0 * a * c);
      
      if (delta > 0.0)
	{
	  secant[0][0] = ((-b) + sqrt(delta)) / (2.0 * a);
	  secant[1][0] = ((-b) - sqrt(delta)) / (2.0 * a);  // useless
	  
	  secant[0][1] = (line[0] * secant[0][0]) + line[1];
	  secant[1][1] = (line[0] * secant[1][0]) + line[1];  // useless
	  
	  vec1[0] = center[0] - secant[0][0];
	  vec1[1] = center[1] - secant[0][1];
	  
	  vec2[0] = center[0];
	  vec2[1] = (line[0] * vec2[0]) + line[1];
	  
	  vec2[0] -= secant[0][0];
	  vec2[1] -= secant[0][1];
	  
	  dotprod = (vec1[0] * vec2[0]) + (vec1[1] * vec2[1]);
	  modprod = (sqrt((vec1[0] * vec1[0]) + (vec1[1] * vec1[1]))) *
	    (sqrt((vec2[0] * vec2[0]) + (vec2[1] * vec2[1])));
	  
	  angle = (acos(dotprod / modprod) * (180.0 / pi));
	  
	  if (center[1] > line[1]) // Circle center over or below the surface
	    angle += 90.0;	
	  else
	    angle = 90.0 - angle;
	  
	  //printf("\n%f\t%f\n%f\t%f\n", secant[0][0], secant[0][1], secant[1][0], secant[1][1]);
	  //printf("\n%f\t%f\n%f\t%f\n", vec1[0], vec1[1], vec2[0], vec2[1]);
	  //printf("dotprod= %f\tmodprod= %f\n", dotprod, modprod);
	  
	}
      else
	{
	  printf("\nWARNING: Negative delta parameter!\n");
	  printf("\na= %f\tb= %f\tc =%f\tdelta= %f\n ", a, b, c, delta); 
	}

      return angle;
    }
  else
    return -1; // Line don't intersect circle
}


void dens_matrix(DATA *set, int nmol, int **densmat, int resolution, int lowest)
{
  int i;
  int row, col;

  for (i = 0; i < nmol; i++)
    {
      row = (int)((set[i].coord[2] - set[lowest].coord[2]) / resolution);
      
      col = (int)(set[i].coord[1] / resolution);

      densmat[row][col]++;
    }
}


int lowest_atom(DATA *set, int natom, int axis)
{
  int i;
  int lowest;

  lowest = 0;

  for(i = 1; i < natom; i++)
    {
      lowest = (set[i].coord[axis] < set[lowest].coord[axis]) ? i : lowest;
    } 

  return lowest;
}


int main(int argc, char *argv[])
{
  int i, j, row, col;
  int timestep;
  int direction;
  int selspe;
  int tstp, type;
  int found;
  int numgas;
  int natoms, totatoms;
  int neighbours;
  int natoms_slice[11];
  int numslices;
  int slc;
  int **densmat;
  int nrow, ncol;
  int dens_ave[11], cells;
  int delta_first, delta_next;
  int num_profile, old_num_profile;
  int lowest;
  int index;
  int print;
  double resolution, half_res, surface_ref[2];
  double box[3][2], box_size[3];
  double coord[3];
  double radius;
  double slice_thickness;
  double profile[500][2], profile_old[500][2];
  double r_fit, center[2], detcoeff[11];
  double ca[11], ca_ave, ca_std, diff;
  double ref_z;
  double cutoff;
  double ycom, delta_y;
  char dumpname[50];
  char check[50], useless[101];
  DATA *molset;
  DATA *slice_1, *slice_2, *slice_3, *slice_4, *slice_5, *slice_6, *slice_7, *slice_8, *slice_9, *slice_10;
  FILE *dumpfile, *ca_file;

  //Read and check command line arguments

  if (argc != 9)
    {
      printf("\nUsage: ./ca_cyl [dump file] [from timestep] [specie] [slices] [resolution] [fit cutoff (-1 = none)] [surface z] [print fit?]\n\n");
      exit(1);
    }
  
  strcpy(dumpname, *(++argv));

  if (isdigit(**(++argv)) && atoi(*argv) > 0)
     timestep = atoi(*argv);
  else
    {
      printf("\nWrong timestep value.\n\n");
      exit(1);
    }

  /*
  if (strcmp(*(++argv), "x") == 0) 
    {
      direction = X;
    }
  else if (strcmp(*argv, "y") == 0) 
    {
      direction = Y;
    }
  else if (strcmp(*argv, "z") == 0) 
    {
      direction = Z;
    }
  else 
    {
      printf("\nInvalid axis.\n\n");
      exit(1);
    }
  */

  if (isdigit(**(++argv)) && atoi(*argv) > 0)
    selspe = atoi(*argv);
  else
    {
      printf("\nWrong specie value.\n\n");
      exit(1);
    }
  
  if (isdigit(**(++argv)) && atof(*argv) >= 1)
    numslices = atof(*argv);
  else
    {
      printf("\nWrong number of slices value.\n\n");
      exit(1);
    }
  
  if (isdigit(**(++argv)) && atof(*argv) > 0)
    {
      resolution = atof(*argv);
      half_res = resolution / 2.0;
    }
  else
    {
      printf("\nWrong resolution value.\n\n");
      exit(1);
    }
  
  cutoff = atof(*(++argv));  // Fit cutoff

  surface_ref[0] = 0.0;          // Reference line slope
  surface_ref[1] = atof(*(++argv));  // Reference line intercept

  if (isdigit(**(++argv)) && (atoi(*argv) == 0 || atoi(*argv) == 1))
    {
      print = atoi(*argv);
    }
  else
    {
      printf("\nWrong print switch value.\n\n");
      exit(1);
    }
  
  // Open dump and results file

  if ((dumpfile = fopen(dumpname, "r")) == NULL)
    {
      printf("Dump file \"%s\" not found.", dumpname);
      exit(1);
    }

  ca_file = fopen("ca_results.dat", "w");

  // Skip until selected starting timestep is found

  fscanf(dumpfile, "%s\n", check);

  found = 0;

  do
    {
       if (strcmp(check, "TIMESTEP") == 0)
	{
	  fscanf(dumpfile, "%i\n", &tstp);
	  
	  printf("Timestep %i: ", tstp);

	  if (tstp >= timestep)
	    {
	      printf("In process\n");
	      
	      found = 1;
	      
	      fscanf(dumpfile, "%*s%*s%*s%*s\n");
	      fscanf(dumpfile, "%i\n", &totatoms);
	      fscanf(dumpfile, "%*[^\n]");
	      
	      molset = malloc(totatoms * sizeof(DATA));
	      slice_1 = malloc(totatoms * sizeof(DATA));
	      slice_2 = malloc(totatoms * sizeof(DATA));
	      slice_3 = malloc(totatoms * sizeof(DATA));
	      slice_4 = malloc(totatoms * sizeof(DATA));
	      slice_5 = malloc(totatoms * sizeof(DATA));
	      slice_6 = malloc(totatoms * sizeof(DATA));
	      slice_7 = malloc(totatoms * sizeof(DATA));
	      slice_8 = malloc(totatoms * sizeof(DATA));
	      slice_9 = malloc(totatoms * sizeof(DATA));
	      slice_10 = malloc(totatoms * sizeof(DATA));

	      for (i = 0; i < 3; i++)
		{
		  fscanf(dumpfile, "%lf%lf\n", &box[i][0], &box[i][1]);
		  //printf("%f\t%f\n",box[i][0], box[i][1]);  //CTRL
		}
	      
	      box_size[0] = box[0][1] - box[0][0];
	      box_size[1] = box[1][1] - box[1][0];
	      box_size[2] = box[2][1] - box[2][0];
	      
	      fscanf(dumpfile, "%*s%*s%*s%*s%*s%*s%*s\n");
	      
	      natoms = 0;

	      for (i = 0; i < totatoms; i++)
		{
		  fscanf(dumpfile, "%i%i%lf%lf%lf\n", 
			 &index, &type, &coord[0], &coord[1], &coord[2]);

		  // Store only selected specie atoms, 
		  // with z coordinate higher than the surface 
		 		  
		  if (type == selspe && (coord[2] * box_size[2]) > surface_ref[1])
		    {
		      molset[natoms].atom_type = type;
		      molset[natoms].coord[0] = (coord[0] * box_size[0]) + box[0][0];
		      molset[natoms].coord[1] = (coord[1] * box_size[1]) + box[1][0];
		      molset[natoms].coord[2] = (coord[2] * box_size[2]) + box[2][0];
		      
		      molset[natoms].index = index;

		      natoms++;
		    }
		}

	      if (natoms == 0)
		{
		   printf("\nSelected atoms type not found.\n\n");
		   exit(1);
		}
	      
	      // Remove gas molecules (set to 0 the atom type)

	      radius = 3;  // It uses a cubic volume to reduce computational time
	      
	      numgas = 0;

	      for (i = 0; i < natoms; i++)
		{
		  neighbours = 0;

		  for (j = 0; j < natoms; j++)
		    {
		      if (fabs(molset[i].coord[0] - molset[j].coord[0]) <= radius &&
			  fabs(molset[i].coord[1] - molset[j].coord[1]) <= radius &&
			  fabs(molset[i].coord[2] - molset[j].coord[2]) <= radius &&
			  i != j)
			{
			  neighbours++;
			}
		    }

		  if (neighbours == 0) 
		    {
		      molset[i].atom_type = -1;
		      numgas++;
		    }
		}

	      printf("Number of gas molecules: %i\n\n", numgas);

	      // Translate the drop so that the y component of its approximate center of mass lies at Ly/2
	      
	      ycom = com_pbc_Y(molset, natoms, selspe, box_size);
	      delta_y = (box_size[1] / 2.0) - ycom;

	      for (i = 0; i < natoms; i++)
		{
		  molset[i].coord[1] += delta_y;
		  molset[i].coord[1] -= box_size[1] * floor(molset[i].coord[1] / box_size[1]);
		}
	      
	      // Slice the sample

	      for (i = 0; i <= numslices; i++)
		{
		  natoms_slice[i] = 0;
		}

	      slice_thickness = box_size[0] / numslices;

	      //printf("\nSlice thickness: %f\n", slice_thickness);
	      
	      for (i = 0; i < natoms; i++)
		{
		  if (molset[i].atom_type != -1)   // Skip gas molecules
		    {
		      slc = (int)((molset[i].coord[0] - box[0][0]) / slice_thickness) + 1;

		      //printf("%i\t%f\t%i\n", i, molset[i].coord[0], slc);
		      
		      if (slc > numslices)
			{
			  printf("\nWARNING: Wrong slice found.\n\n");
			  //printf("%i\t%f\t%i\n", molset[i].index, molset[i].coord[0], slc);
			  continue; // Jump to next iteration
			}
		      
		      switch (slc)
			{
			case (1):
			  {
			    slice_1[natoms_slice[1]].coord[0] = molset[i].coord[0];
			    slice_1[natoms_slice[1]].coord[1] = molset[i].coord[1];
			    slice_1[natoms_slice[1]].coord[2] = molset[i].coord[2];
			    
			    natoms_slice[1]++;
			    
			    break;
			  }
			case (2):
			  {
			    slice_2[natoms_slice[2]].coord[0] = molset[i].coord[0];
			    slice_2[natoms_slice[2]].coord[1] = molset[i].coord[1];
			    slice_2[natoms_slice[2]].coord[2] = molset[i].coord[2];
			    
			    natoms_slice[2]++;
			    
			    break;
			  }
			case (3):
			  {
			    slice_3[natoms_slice[3]].coord[0] = molset[i].coord[0];
			    slice_3[natoms_slice[3]].coord[1] = molset[i].coord[1];
			    slice_3[natoms_slice[3]].coord[2] = molset[i].coord[2];
			    
			    natoms_slice[3]++;
			    
			    break;
			  }
			case (4):
			  {
			    slice_4[natoms_slice[4]].coord[0] = molset[i].coord[0];
			    slice_4[natoms_slice[4]].coord[1] = molset[i].coord[1];
			    slice_4[natoms_slice[4]].coord[2] = molset[i].coord[2];
			    
			    natoms_slice[4]++;
			    
			    break;
			  }
			case (5):
			  {
			    slice_5[natoms_slice[5]].coord[0] = molset[i].coord[0];
			    slice_5[natoms_slice[5]].coord[1] = molset[i].coord[1];
			    slice_5[natoms_slice[5]].coord[2] = molset[i].coord[2];
			    
			    natoms_slice[5]++;
			    
			    break;
			  }
			case (6):
			  {
			    slice_6[natoms_slice[6]].coord[0] = molset[i].coord[0];
			    slice_6[natoms_slice[6]].coord[1] = molset[i].coord[1];
			    slice_6[natoms_slice[6]].coord[2] = molset[i].coord[2];
			    
			    natoms_slice[6]++;
			    
			    break;
			  }
			case (7):
			  {
			    slice_7[natoms_slice[7]].coord[0] = molset[i].coord[0];
			    slice_7[natoms_slice[7]].coord[1] = molset[i].coord[1];
			    slice_7[natoms_slice[7]].coord[2] = molset[i].coord[2];
			    
			    natoms_slice[7]++;
			    
			    break;
			  }
			case (8):
			  {
			    slice_8[natoms_slice[8]].coord[0] = molset[i].coord[0];
			    slice_8[natoms_slice[8]].coord[1] = molset[i].coord[1];
			    slice_8[natoms_slice[8]].coord[2] = molset[i].coord[2];
			    
			    natoms_slice[8]++;
			    
			    break;
			  }
			case (9):
			  {
			    slice_9[natoms_slice[9]].coord[0] = molset[i].coord[0];
			    slice_9[natoms_slice[9]].coord[1] = molset[i].coord[1];
			    slice_9[natoms_slice[9]].coord[2] = molset[i].coord[2];
			    
			    natoms_slice[9]++;
			    
			    break;
			  }
			case (10):
			  {
			    slice_10[natoms_slice[10]].coord[0] = molset[i].coord[0];
			    slice_10[natoms_slice[10]].coord[1] = molset[i].coord[1];
			    slice_10[natoms_slice[10]].coord[2] = molset[i].coord[2];
			    
			    natoms_slice[10]++;
			    
			    break;
			  }
			}
		    }
		    
		}
	      /*
	      for (i = 1; i <= 10; i++)
		printf("%i\n", natoms_slice[i]);
	      */

	      // Calculate density matrix and contact angle
	      // Use the lowest atom as reference for the bottom of the density matrix

	      for (slc = 1; slc <= numslices; slc++)
		{
		  switch (slc)
		    {
		    case(1):
		      {
			lowest = lowest_atom(slice_1, natoms_slice[1], 2);
			ref_z = slice_1[lowest].coord[2];

			nrow = (int)((box_size[2] - slice_1[lowest].coord[2]) / resolution) + 1;
			ncol = (int)(box_size[1] / resolution) + 1;

			// Allocate density matrix

			densmat = (int**)calloc(nrow, sizeof(int*));   // Rows
			for (j = 0; j < nrow; j++)
			  densmat[j] = calloc(ncol, sizeof(int));      // Coloumns

			dens_matrix(slice_1, natoms_slice[1], densmat, resolution, lowest);

			break;
		      }
		    case(2):
		      {
			lowest = lowest_atom(slice_2, natoms_slice[2], 2);
			ref_z = slice_2[lowest].coord[2];

			nrow = (int)((box_size[2] - slice_2[lowest].coord[2]) / resolution) + 1;
			ncol = (int)(box_size[1] / resolution) + 1;

			// Allocate density matrix

			densmat = (int**)calloc(nrow, sizeof(int*));   // Rows
			for (j = 0; j < nrow; j++)
			  densmat[j] = calloc(ncol, sizeof(int));      // Coloumns

			dens_matrix(slice_2, natoms_slice[2], densmat, resolution, lowest);
			break;
		      }
		    case(3):
		      {
			lowest = lowest_atom(slice_3, natoms_slice[3], 2);
			ref_z = slice_3[lowest].coord[2];

			nrow = (int)((box_size[2] - slice_3[lowest].coord[2]) / resolution) + 1;
			ncol = (int)(box_size[1] / resolution) + 1;

			// Allocate density matrix

			densmat = (int**)calloc(nrow, sizeof(int*));   // Rows
			for (j = 0; j < nrow; j++)
			  densmat[j] = calloc(ncol, sizeof(int));      // Coloumns

			dens_matrix(slice_3, natoms_slice[3], densmat, resolution, lowest);
			break;
		      }
		    case(4):
		      {
			lowest = lowest_atom(slice_4, natoms_slice[4], 2);
			ref_z = slice_4[lowest].coord[2];

			nrow = (int)((box_size[2] - slice_4[lowest].coord[2]) / resolution) + 1;
			ncol = (int)(box_size[1] / resolution) + 1;

			// Allocate density matrix

			densmat = (int**)calloc(nrow, sizeof(int*));   // Rows
			for (j = 0; j < nrow; j++)
			  densmat[j] = calloc(ncol, sizeof(int));      // Coloumns

			dens_matrix(slice_4, natoms_slice[4], densmat, resolution, lowest);
			break;
		      }
		    case(5):
		      {
			lowest = lowest_atom(slice_5, natoms_slice[5], 2);
			ref_z = slice_5[lowest].coord[2];

			nrow = (int)((box_size[2] - slice_5[lowest].coord[2]) / resolution) + 1;
			ncol = (int)(box_size[1] / resolution) + 1;

			// Allocate density matrix

			densmat = (int**)calloc(nrow, sizeof(int*));   // Rows
			for (j = 0; j < nrow; j++)
			  densmat[j] = calloc(ncol, sizeof(int));      // Coloumns

			dens_matrix(slice_5, natoms_slice[5], densmat, resolution, lowest);
			break;
		      }
		    case(6):
		      {
			lowest = lowest_atom(slice_6, natoms_slice[6], 2);
			ref_z = slice_6[lowest].coord[2];

			nrow = (int)((box_size[2] - slice_6[lowest].coord[2]) / resolution) + 1;
			ncol = (int)(box_size[1] / resolution) + 1;

			// Allocate density matrix

			densmat = (int**)calloc(nrow, sizeof(int*));   // Rows
			for (j = 0; j < nrow; j++)
			  densmat[j] = calloc(ncol, sizeof(int));      // Coloumns

			dens_matrix(slice_6, natoms_slice[6], densmat, resolution, lowest);
			break;
		      }
		    case(7):
		      {
			lowest = lowest_atom(slice_7, natoms_slice[7], 2);
			ref_z = slice_7[lowest].coord[2];

			nrow = (int)((box_size[2] - slice_7[lowest].coord[2]) / resolution) + 1;
			ncol = (int)(box_size[1] / resolution) + 1;

			// Allocate density matrix

			densmat = (int**)calloc(nrow, sizeof(int*));   // Rows
			for (j = 0; j < nrow; j++)
			  densmat[j] = calloc(ncol, sizeof(int));      // Coloumns

			dens_matrix(slice_7, natoms_slice[7], densmat, resolution, lowest);
			break;
		      }
		    case(8):
		      {
			lowest = lowest_atom(slice_8, natoms_slice[8], 2);
			ref_z = slice_8[lowest].coord[2];

			nrow = (int)((box_size[2] - slice_8[lowest].coord[2]) / resolution) + 1;
			ncol = (int)(box_size[1] / resolution) + 1;

			// Allocate density matrix

			densmat = (int**)calloc(nrow, sizeof(int*));   // Rows
			for (j = 0; j < nrow; j++)
			  densmat[j] = calloc(ncol, sizeof(int));      // Coloumns

			dens_matrix(slice_8, natoms_slice[8], densmat, resolution, lowest);
			break;
		      }
		    case(9):
		      {
			lowest = lowest_atom(slice_9, natoms_slice[9], 2);
			ref_z = slice_9[lowest].coord[2];

			nrow = (int)((box_size[2] - slice_9[lowest].coord[2]) / resolution) + 1;
			ncol = (int)(box_size[1] / resolution) + 1;

			// Allocate density matrix

			densmat = (int**)calloc(nrow, sizeof(int*));   // Rows
			for (j = 0; j < nrow; j++)
			  densmat[j] = calloc(ncol, sizeof(int));      // Coloumns

			dens_matrix(slice_9, natoms_slice[9], densmat, resolution, lowest);
			break;
		      }
		    case(10):
		      {
			lowest = lowest_atom(slice_10, natoms_slice[10], 2);
			ref_z = slice_10[lowest].coord[2];

			nrow = (int)((box_size[2] - slice_10[lowest].coord[2]) / resolution) + 1;
			ncol = (int)(box_size[1] / resolution) + 1;

			// Allocate density matrix

			densmat = (int**)calloc(nrow, sizeof(int*));   // Rows
			for (j = 0; j < nrow; j++)
			  densmat[j] = calloc(ncol, sizeof(int));      // Coloumns

			dens_matrix(slice_10, natoms_slice[10], densmat, resolution, lowest);
			break;
		      }
		    }
		  
		  // Calculate average density and warn if too low

		  dens_ave[slc] = 0;
		  cells = 0;

		  for (row = 0; row < nrow; row++)
		    {
		      for (col = 0; col < ncol; col++)
			{
			  if (densmat[row][col] != 0)
			    {
			      dens_ave[slc] += densmat[row][col];
			      cells++;
			    }
			}
		    }
		
		  dens_ave[slc] /= cells;

		  if (dens_ave[slc] < 4)
		    {
		       printf("\nWARNING: Average density too low: %i\n\n", dens_ave[slc]);
		    }
		  

		  /*
		  // Print density matrix on stdout 
		  
		  int sum;

		  sum = 0;

		  for (row = 0; row < nrow; row++)
		    {
		      for (col = 0; col < ncol; col++)
			{
			  printf("%i  ", densmat[row][col]);
			  sum += densmat[row][col];
			}
		      printf("\n");
		    }
		  
		  printf("\nSum: %i\n\n", sum);
		  */


		  // Find drop profile points

		  num_profile = 0;
	
		  for (row = 0; row < nrow; row++) 
		    {
		      for(col = 0; col < ncol; col++)        // Find sx point
			{
			  if (densmat[row][col + 1] != 0)    // Skip zeroes
			    {
			      delta_first = abs(densmat[row][col] - densmat[row][col + 1]);
			      delta_next = abs(densmat[row][col + 1] - densmat[row][col + 2]);
			      
			      if((delta_first > delta_next && delta_first > 2) || densmat[row][col] > 2)
				{
				  profile[num_profile][0] = ((col + 1) * resolution);
				  profile[num_profile][1] = ((row + 1) * resolution) - half_res + ref_z;

				  num_profile++;
				  
				  break;  // Point found; exit cycle
				}
			    }
			} 

		      for (col = (ncol - 1); col > 0; col--)     // Find dx point
			{
			  if (densmat[row][col - 1] != 0)    // Skip zeroes
			    {
			      delta_first = abs(densmat[row][col] - densmat[row][col - 1]);
			      delta_next = abs(densmat[row][col - 1] - densmat[row][col - 2]);
			      
			      if((delta_first > delta_next && delta_first > 2) || densmat[row][col] > 2)
				{
				  profile[num_profile][0] = (col * resolution);
				  profile[num_profile][1] = ((row + 1) * resolution) - half_res + ref_z;

				  num_profile++;
				  
				  break;  // Point found; exit cycle
				}
			    }
			}
		    }

		  // Remove points below cutoff before fitting

		  if (cutoff > 0)
		    {
		      old_num_profile = num_profile;
		      num_profile = 0;

		      // Copy profile points

		      for (row = 0; row < old_num_profile; row++)
			{
			  profile_old[row][0] = profile[row][0];
			  profile_old[row][1] = profile[row][1];
			}		

		      // Check

		      for (row = 0; row < old_num_profile; row++)
			{
			  
			  if (profile_old[row][1] >= cutoff)
			    {
			      profile[num_profile][0] = profile_old[row][0];
			      profile[num_profile][1] = profile_old[row][1];
			      
			      num_profile++;
			    }
			}
		    }

		  /*
		  for (row = 0; row < num_profile; row++) 
		    {
		      printf("%i\t%f\t%f\n", row, profile[row][0], profile[row][1]);
		    }
		  */

		  // Calculate best fitting circle

		  circlefit(profile, num_profile, &r_fit, center, &detcoeff[slc]);
		  
		  //printf("\nRadius: %f\nCenter: %f\t%f\nError: %f %%\n", r_fit, center[0], center[1], detcoeff[slc]);

		  // Calculate contact angle

		  ca[slc] = contact_angle(center, r_fit, surface_ref);

		  //printf("Contact angle: %f\n\n", ca[slc]);
		  
		  
		  // Print fit if requested

		  if (print == 1 && slc == 1)
		    {
		      FILE *contour, *fit, *atoms, *surface;
		      double x, y, rad;
		      
		      contour = fopen("contour.dat", "w");
		      fit = fopen("fit.dat", "w");
		      atoms = fopen("atoms.dat", "w");
		      surface = fopen("surface.dat", "w");
		      
		      for (i = 0; i < num_profile; i++)
			{
			  fprintf(contour, "%f\t%f\n", profile[i][0], profile[i][1]);
			}
		      
		      for (i = 0; i <= 1000; i++)
			{
			  rad = ((double)i / (double)1000) * 3.14;
			  
			  x = (r_fit * cos(rad)) + center[0];
			  y = (r_fit * sin(rad)) + center[1];
			  fprintf(fit, "%f\t%f\n",x ,y);
			}
		      
		      for (i = 1000; i >= 0; i--)
			{
			  rad = ((double)i / (double)1000) * 3.14;
			  
			  x = (r_fit * cos(rad)) + center[0];
			  y = (r_fit * sin(rad) * -1.0) + center[1];
			  fprintf(fit, "%f\t%f\n",x , y);
			}
		      
		      for (i = 0; i <= 1000; i++)
			{
			  fprintf(surface, "%f\t%f\n", ((box_size[1]/1000) * (double)i), surface_ref[1]);
			}
		      
		      for (i = 0; i < natoms_slice[1]; i++)
			{
			  fprintf(atoms, "%f\t%f\n", slice_1[i].coord[1], slice_1[i].coord[2]);
			}
		      
		      fclose(contour);
		      fclose(fit);
		      fclose(atoms);
		      fclose(surface);
		    
		      return 0;
		    }
		    

		  // Deallocate density matrix
		  
		  for (j = 0; j < nrow; j++)
		    free(densmat[j]);
		  
		  free(densmat);
		}
     
	      // Calculate average contact angle and standard deviation

	      ca_ave = 0.0;

	      for (slc = 1; slc <= numslices; slc++)
		{
		  ca_ave += ca[slc];
		}

	      ca_ave /= (double)numslices;

	      ca_std = 0.0;
	      
	       for (slc = 1; slc <= numslices; slc++)
		{
		  diff = ca[slc] - ca_ave;
		  ca_std += diff * diff;
		}
	      
	       ca_std /= (double)numslices;

	       // Print results

	       if (tstp == timestep)
		 {
		   fprintf(ca_file, "# Timestep\t");
		   
		   for (slc = 1; slc <= numslices; slc++)
		     {
		       fprintf(ca_file, "CA_%i\tFit_err_%i\t", slc, slc);
		     }
		   
		   fprintf(ca_file, "CA_ave\tCA_std\n");
		 }
	       
	       fprintf(ca_file, "%7i", tstp);
	      
	        for (slc = 1; slc <= numslices; slc++)
		{
		  fprintf(ca_file, "%10.3f%10.3f", ca[slc], detcoeff[slc]);
		}

		fprintf(ca_file, "%10.3f%10.3f\n", ca_ave, ca_std);

	      // Deallocate slices arrays
	      
	      free(molset);
	      free(slice_1);
	      free(slice_2);
	      free(slice_3);
	      free(slice_4);
	      free(slice_5);
	      free(slice_6);
	      free(slice_7);
	      free(slice_8);
	      free(slice_9);
	      free(slice_10);


	    }
	  else
	    printf("Skipped\n");

	}

 
    }
  while (fscanf(dumpfile, "%s\n", check) != EOF);

  if (found == 0)
    {
      printf("ERROR: Selected timestep not found in file \"%s\"\n\n", dumpname);
      exit(1);
    }

  fclose(dumpfile);
  fclose(ca_file);

  return 0;
}


