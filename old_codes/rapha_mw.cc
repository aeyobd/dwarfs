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


