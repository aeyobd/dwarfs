#include "allvars.h"
#include "proto.h"
#include <math.h>
#include <gsl/gsl_sf_gamma.h> //Added by JL 07/26 for incomplete gamma function

#ifdef EXTERNAL_POTENTIAL
#ifndef PERT_PROFILE
#error Specify perturbation mass profile.
#endif

#if PERT_PROFILE == POINTMASS
FLOAT f_ext(int i, int j)
{
  return ext_amplitude() * (-All.G) * All.pert_M * P[i].Pos[j] / pow(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2],1.5); //no softening!
}

//--------------------------------------------

#elif PERT_PROFILE == PLUMMER
FLOAT f_ext(int i, int j)
{
  FLOAT r = pow(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2],0.5);
  return ext_amplitude() * (-All.G) * All.pert_M * (pow(r,3)/pow(pow(r,2) + pow(a,2), 1.5)) * P[i].Pos[j] / pow(pow(r,2) + pow(All.SofteningHalo,2), 1.5);
}

//--------------------------------------------

#elif PERT_PROFILE == HERNQUIST
FLOAT f_ext(int i, int j)
{
  FLOAT r = pow(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2],0.5);
  return ext_amplitude() * (-All.G) * All.pert_M * P[i].Pos[j] / ( r * pow(r+a, 2) ) ;				//James Lane edits removal of softening
}

//--------------------------------------------

#elif PERT_PROFILE == SPHERICALEXPDISC
FLOAT f_ext(int i, int j)
{
  FLOAT r = pow(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2],0.5);
  return ext_amplitude() * (-All.G) * All.pert_M * (1-exp(-r/All.pert_a)*(1+r/All.pert_a)) * P[i].Pos[j] / pow(pow(r,2) + pow(All.SofteningHalo,2), 1.5);
}

//--------------------------------------------

#elif PERT_PROFILE == ISOTHERMAL //truncated isothermal sphere
FLOAT f_ext(int i, int j)
{
  FLOAT r = pow(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2],0.5);
  return ext_amplitude() * (-All.G) * All.pert_M  * ((r < All.pert_a) ? (r / All.pert_a) : (1.0)) *  P[i].Pos[j] / pow(pow(r, 2) + pow(All.SofteningHalo, 2), 1.5);
}

//--------------------------------------------

#elif PERT_PROFILE == MIYAMOTONAGAI //Miyamoto-Nagai disk potential, begin James Lane edit 04/28
FLOAT f_ext(int i, int j)
{
  FLOAT r = pow( P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] , 0.5); //Cylindrical coordinate r
  FLOAT z_scale = pow( P[i].Pos[2] * P[i].Pos[2] + All.pert_b * All.pert_b , 0.5);
  FLOAT z_comp = 1.0;
  if ( j == 2 ) //When coordinate is z there is an extra term
  {
    z_comp += ( All.pert_a)/(z_scale);
  }
  return ext_amplitude() * (-All.G) * All.pert_M * ( ( P[i].Pos[j] ) / pow( pow( All.pert_a + z_scale , 2) + pow(r,2) , 1.5) ) * z_comp;
}

//-------------------------------------------- end James Lane edit 04/28

#elif PERT_PROFILE == MILKYWAY //Milky-Way potential, begin James Lane edit 05/19
FLOAT f_ext(int i, int j)
{
  // Disk: cylindrical radial scale = a, height scale = b, mass = Md
  // Bulge: scale length = c, mass = Mb
  // Halo: scale length = d, mass = Mh
  FLOAT r_cyl = pow( P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] , 0.5); //Cylindrical coordinate r
  FLOAT r_sph = pow(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2],0.5); //Spherical coordinate r
  FLOAT z_scale = pow( P[i].Pos[2] * P[i].Pos[2] + All.pert_b * All.pert_b , 0.5);
  FLOAT z_comp = 1.0;
  if ( j == 2 ) //When coordinate is z there is an extra term
  {
    z_comp += (All.pert_a)/(z_scale);
  }
  FLOAT disk_comp = All.pert_Md * ( ( P[i].Pos[j] ) / pow( pow( All.pert_a + z_scale , 2) + pow(r_cyl,2) , 1.5) ) * z_comp; //Miyamoto-Nagai disk
  FLOAT bulge_comp = All.pert_Mb * P[i].Pos[j] / ( r_sph * pow(r_sph + All.pert_c, 2) ); //Hernquist bulge
  FLOAT halo_comp = All.pert_Mh * P[i].Pos[j] / ( r_sph * pow(r_sph + All.pert_d, 2) ); //Hernquist halo
  return ext_amplitude() * (-All.G) * ( disk_comp + bulge_comp + halo_comp );
}

//-------------------------------------------- end James Lane edit 05/19

#elif PERT_PROFILE == MWPOTENTIAL2014 //Bovy+ 2015 Milky Way potential, Begin James Lane edit 07/26
FLOAT f_ext(int i, int j)
{
  // Disk, Miyamoto-Nagai: cylindrical radial scale = a, height scale = b, mass = Md
  // Bulge, Power-law w/ Exponential cutoff: power-law index = alpha, exponential scale = c, Mass = Mb
  // Halo, NFW: concentration = conc, scale length = d, Virial mass = Mhvir
  FLOAT r_cyl = pow( P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] , 0.5); //Cylindrical coordinate r
  FLOAT r_sph = pow(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2],0.5); //Spherical coordinate r
  FLOAT z_scale = pow( P[i].Pos[2] * P[i].Pos[2] + All.pert_b * All.pert_b , 0.5);
  FLOAT z_comp = 1.0;
  FLOAT g_conc = pow( log( 1 + All.pert_conc ) - ( All.pert_conc / ( All.pert_conc + 1 ) )  , -1 );
  //FLOAT gamma_comp = gsl_sf_gamma_inc( 1.5-(All.pert_alpha/2) , pow(r_sph/All.pert_c,2) ) / gsl_sf_gamma( 1.5-(All.pert_alpha/2) );
  if ( j == 2 )
  {
    z_comp += (All.pert_a)/(z_scale);
  }
  FLOAT disk_comp = All.pert_Md * ( ( P[i].Pos[j] ) / pow( pow( All.pert_a + z_scale , 2) + pow(r_cyl,2) , 1.5) ) * z_comp; //Miyamoto-Nagai disk
  FLOAT bulge_comp = All.pert_Mb * P[i].Pos[j] * ( gsl_sf_gamma_inc_P( 1.5-(All.pert_alpha/2) , pow(r_sph/All.pert_c,2)) ) / pow(r_sph,3);
  FLOAT halo_comp = All.pert_Mhvir * g_conc * P[i].Pos[j] * ( ( log(1+r_sph/All.pert_d) / pow(r_sph,3) )  - 1 / ( pow(r_sph,2) * (r_sph+All.pert_d ) ) );
  return ext_amplitude() * (-All.G) * ( disk_comp + bulge_comp + halo_comp );
}

//-------------------------------------------- end James Lane edit 07/26


//-------------------------------------------- begin RE edit 01/20
#elif PERT_PROFILE == RAPHAMW // as in Errani+20

FLOAT f_ext(int i, int j)
{
  // particle ID i ; coordinate j
  // Disk Thin   Miyamoto-Nagai: cylindrical radial scale = a_thin, height scale = b_thin, mass = Md_thin
  // Disk Thick  Miyamoto-Nagai: cylindrical radial scale = a_thick, height scale = b_thick, mass = Md_thick
  
  // Bulge, Hernquist Mass = Mb  scale abulge
  // Halo, NFW: concentration = conc, scale length = d, Virial mass = Mhvir
  FLOAT r_cyl = pow( P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] , 0.5); //Cylindrical coordinate r
  FLOAT r_sph = pow(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2],0.5); //Spherical coordinate r
  FLOAT z_scale_thin  = pow( P[i].Pos[2] * P[i].Pos[2] + All.pert_b_thin  * All.pert_b_thin , 0.5);
  FLOAT z_scale_thick = pow( P[i].Pos[2] * P[i].Pos[2] + All.pert_b_thick * All.pert_b_thick , 0.5);
  FLOAT z_comp_thin = 1.0;
  FLOAT z_comp_thick = 1.0;
  FLOAT g_conc = pow( log( 1 + All.pert_conc ) - ( All.pert_conc / ( All.pert_conc + 1 ) )  , -1 );
  //FLOAT gamma_comp = gsl_sf_gamma_inc( 1.5-(All.pert_alpha/2) , pow(r_sph/All.pert_c,2) ) / gsl_sf_gamma( 1.5-(All.pert_alpha/2) );
  if ( j == 2 )
  {
    z_comp_thin  += (All.pert_a_thin)/(z_scale_thin);
    z_comp_thick += (All.pert_a_thick)/(z_scale_thick);
  }
  FLOAT thin_disk_comp  = All.pert_Md_thin  * ( ( P[i].Pos[j] ) / pow( pow( All.pert_a_thin  + z_scale_thin  , 2) + pow(r_cyl,2) , 1.5) ) * z_comp_thin; 
  FLOAT thick_disk_comp = All.pert_Md_thick * ( ( P[i].Pos[j] ) / pow( pow( All.pert_a_thick + z_scale_thick , 2) + pow(r_cyl,2) , 1.5) ) * z_comp_thick; 
  // Hernquist   
  FLOAT bulge_comp = All.pert_Mb * P[i].Pos[j] / ( r_sph * pow(r_sph+All.pert_abulge, 2) ) ;
  FLOAT halo_comp = All.pert_Mhvir * g_conc * P[i].Pos[j] * ( ( log(1+r_sph/All.pert_d) / pow(r_sph,3) )  - 1 / ( pow(r_sph,2) * (r_sph+All.pert_d ) ) );
  return ext_amplitude() * (-All.G) * ( thick_disk_comp + thin_disk_comp + bulge_comp + halo_comp );
}

//-------------------------------------------- end RE edit 01/20

#elif PERT_PROFILE == CUSTOMTABLE
FLOAT f_ext(int i, int j)
{
  FLOAT r = pow(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2],0.5);
  FLOAT M;
  int r_index;
  if(r < All.pert_interp_r[0])
    {
      //do linear interpolation in time
      M = All.pert_interp_M[All.pert_t_index][0] + (All.pert_interp_M[All.pert_t_index+1][0] - All.pert_interp_M[All.pert_t_index][0]) * (All.Time - All.pert_interp_t[All.pert_t_index]) / (All.pert_interp_t[All.pert_t_index+1] - All.pert_interp_t[All.pert_t_index]);
      return ext_amplitude() * (-All.G) * M * P[i].Pos[j] / pow(pow(r,2) + pow(All.SofteningHalo,2), 1.5);
    }
  for(r_index = 0; (r > All.pert_interp_r[r_index+1]) && (All.pert_interp_r[r_index+1] > 0); r_index++);
  if(All.pert_interp_r[r_index+1] < 0)
    {
      //do linear interpolation in time
      M = All.pert_interp_M[All.pert_t_index][r_index] + (All.pert_interp_M[All.pert_t_index+1][r_index] - All.pert_interp_M[All.pert_t_index][r_index]) * (All.Time - All.pert_interp_t[All.pert_t_index]) / (All.pert_interp_t[All.pert_t_index+1] - All.pert_interp_t[All.pert_t_index]);
      return ext_amplitude() * (-All.G) * M * P[i].Pos[j] / pow(pow(r,2) + pow(All.SofteningHalo,2), 1.5);
    }
  //do bilinear interpolation
  M = (All.pert_interp_M[All.pert_t_index][r_index] * (All.pert_interp_t[All.pert_t_index+1] - All.Time) * (All.pert_interp_r[r_index+1] - r) + All.pert_interp_M[All.pert_t_index+1][r_index] * (All.Time - All.pert_interp_t[All.pert_t_index]) * (All.pert_interp_r[r_index+1] - r) + All.pert_interp_M[All.pert_t_index][r_index+1] * (All.pert_interp_t[All.pert_t_index+1] - All.Time) * (r - All.pert_interp_r[r_index]) + All.pert_interp_M[All.pert_t_index+1][r_index+1] * (All.Time - All.pert_interp_t[All.pert_t_index]) * (r - All.pert_interp_r[r_index])) / ((All.pert_interp_t[All.pert_t_index+1] - All.pert_interp_t[All.pert_t_index]) * (All.pert_interp_r[r_index+1] - All.pert_interp_r[r_index]));
  return ext_amplitude() * (-All.G) * M * P[i].Pos[j] / pow(pow(r,2) + pow(All.SofteningHalo,2), 1.5);
}

//--------------------------------------------

#else
#error Unknown perturbation mass profile.
#endif

//--------------------------------------------


#ifndef PERT_TIMEDEP
#error Specify perturbation time dependence.
#endif

#if PERT_TIMEDEP == LINEAR
FLOAT ext_amplitude()
{
  unsigned int loop;
  for(loop = 1; loop < All.n_anchors; loop++)
    {
      if(All.Time > All.t_anchors[loop])
	continue;
      return All.amp_anchors[loop-1] + (All.Time - All.t_anchors[loop-1]) * (All.amp_anchors[loop] - All.amp_anchors[loop-1]) / (All.t_anchors[loop] - All.t_anchors[loop-1]);
    }
  printf("Error in ext_amplitude!"); //should abort here perhaps
  return -1.0; //should never be reached
}

#else
#error Unknown perturbation time dependence.
#endif

#endif


