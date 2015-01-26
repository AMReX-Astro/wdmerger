#include "mangle_names.h"
#include <hdf5.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

#include <stdio.h>
#include <unistd.h>

/* Parallel IO read doesn't work on Intrepid for some reason */
#ifndef IBM
#define TURB_PARALLELIO
#endif

#define H5_FAIL -1

int Driver_abortFlashC(char* message);

/* globals shared between functions */

hid_t turb_vel_file_id, turb_vel_dset_id, turb_vel_dspace_id, dxfer_template;
int USE_COLLECTIVE = 0;
#ifdef TURB_PARALLELIO
static MPI_Info info;
#endif

/* this bounding box specifies the region which the data in the
   turbulence field file covers.  Data is assumed to be
   a cell-centered grid of values within this box.
 */
double tf_bbox[3][2];

/* size of turbulence field grid */
int tf_size[3];

/* some convenient metadata 
   corner is the x,y,z coordinates of the center of the corner cell
          in the turbulence field
   delta  is the x,y,z grid spacing in the turbulence field
*/
double tf_corner[3], tf_delta[3];

/* number of cells corresponding to change in refinement level
 * to smooth turbulence field */
int s;

/* 
  Setup function opens file and dataset, gets a dataspace and calculates
  metadata based on size of turbulence field grid.
*/
void FTOC(sim_turb_field_setup)(double* xmin, double* xmax, double* ymin, double* ymax,
                                double* zmin, double* zmax, int* smooth,
                                char filename[4096], int* useCollectiveHDF5) {

  int i,d;
  hsize_t dset_size[4];
  hid_t file_access_id;
  herr_t ret;

  char c_filename[4096];

  tf_bbox[0][0] = *xmin; tf_bbox[0][1] = *xmax;
  tf_bbox[1][0] = *ymin; tf_bbox[1][1] = *ymax;
  tf_bbox[2][0] = *zmin; tf_bbox[2][1] = *zmax;

  /* make a null terminated filename string */
  i=0;
  while ( filename[i] != ' ' && i<4095 ) {
    c_filename[i] = filename[i];
    i++;
  }
  c_filename[i] = '\0';

  file_access_id = H5Pcreate(H5P_FILE_ACCESS);
  if ( file_access_id == H5_FAIL )
     Driver_abortFlashC("H5Pcreate: file access failed\n");
  dxfer_template = H5Pcreate(H5P_DATASET_XFER);
  if ( dxfer_template == H5_FAIL )
     Driver_abortFlashC("H5Pcreate: dataset xfer failed\n");

#ifdef IBM
  /* Set sieve size for Intrepid
   * https://wiki.alcf.anl.gov/index.php/I_O_Tuning
   */
  ret = H5Pset_sieve_buf_size(file_access_id, 4*1024*1024);
  if ( ret == H5_FAIL )
     Driver_abortFlashC("H5Pset_sieve_buf_size: failed\n");
#endif

#ifdef TURB_PARALLELIO
  MPI_Info_create(&info);
  /* no particular keys to set here, but could do so in info*/
/*acc  ret = H5Pset_fapl_mpio(file_access_id, MPI_COMM_WORLD, info);*/
  if ( ret == H5_FAIL )
     Driver_abortFlashC("H5Pset_fapl_mpio: failed\n");

  USE_COLLECTIVE = *useCollectiveHDF5;
  if (USE_COLLECTIVE) {
/*acc     ret = H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_COLLECTIVE);*/
     if ( ret == H5_FAIL )
        Driver_abortFlashC("H5Pset_dxpl_mpio: failed\n");
  }
#endif
  /* tell calling subroutine if collective io is enabled */
  *useCollectiveHDF5 = USE_COLLECTIVE;

  /* open file, will be parallel if parallel is enabled above, 
     otherwise will be serial*/
  turb_vel_file_id = H5Fopen(c_filename, H5F_ACC_RDONLY, file_access_id);
  if ( turb_vel_file_id == H5_FAIL ) 
     Driver_abortFlashC("H5Fopen: file open failed\n");

  ret = H5Pclose(file_access_id);
  if ( ret == H5_FAIL )
     Driver_abortFlashC("H5Pclose: file access close failed\n");

  /* open the dataset and keep it open */
  turb_vel_dset_id = H5Dopen(turb_vel_file_id, "/velocity");
  if ( turb_vel_dset_id == H5_FAIL )
     Driver_abortFlashC("H5Dopen: dataset open failed\n");

  /* open a dataspace for it that we will use for reads */
  turb_vel_dspace_id = H5Dget_space(turb_vel_dset_id);
  if ( turb_vel_dset_id == H5_FAIL )
     Driver_abortFlashC("H5Dget_space: dataset space open failed\n");

  /* get and store some metadata about the turbulence field layout 
   * for convenience.
   * we want size, delta, and corner to correspond to target smoothed field
   * where the interpolation or averaging will be done. */
  s = (int) pow(2,*smooth);
  ret = H5Sget_simple_extent_dims(turb_vel_dspace_id, dset_size, NULL);
  if ( ret == H5_FAIL )
     Driver_abortFlashC("H5Sget_simple_extent_dims: failed\n");

  for (d=0;d<3;d++) {
     /* dimension order in storage is k,j,i */
     i = 2-d;
     tf_size[d]=(int) dset_size[i];
     if ( tf_size[d] % s != 0 )
        Driver_abortFlashC("error smoothing turbulence field\n");
     tf_size[d] /= s;
     tf_delta[d] = (tf_bbox[d][1]-tf_bbox[d][0])/tf_size[d];
     tf_corner[d] = tf_bbox[d][0] + 0.5*tf_delta[d];
  }


  return;
}


/* need to perform 8 reads for collective mode to work right */
void sim_turb_field_null_read_c() {

  herr_t ret;
  int i;
  ret = H5Sselect_none( turb_vel_dspace_id );
  if ( ret == H5_FAIL )
     Driver_abortFlashC("H5Sselect_none: failed\n");
  for (i = 0; i < 8; i++) {
      MPI_Barrier( MPI_COMM_WORLD );
      ret = H5Dread( turb_vel_dset_id, H5T_NATIVE_DOUBLE,
                     H5S_ALL, turb_vel_dspace_id, 
                     dxfer_template, NULL );
      if ( ret == H5_FAIL )
         Driver_abortFlashC("H5Dread: failed\n");
  }
  return;
}
/* this is the fortran interface */
void FTOC(sim_turb_field_null_read)(void) {
  sim_turb_field_null_read_c();
}


/*
  get_vel function retrieves data for this piece of turbulence field
 
  For efficiency only the necessary piece of the file is read.  This is then
  trilinearly interpolated to the positions provided.
*/

void FTOC(sim_turb_field_get_vel)(int* isize, int* jsize, int* ksize,
                                  double* velx, double* vely, double* velz, 
                                  double* iCoords, double* jCoords, double* kCoords) {

  int min_i_t, max_i_t, min_j_t, max_j_t, min_k_t, max_k_t, test;
  unsigned int size_i_t, size_j_t, size_k_t;
  unsigned int size_is_t, size_js_t, size_ks_t;
  unsigned int on_edge_i, on_edge_j, on_edge_k;
  unsigned int use_cell_center_x, use_cell_center_y, use_cell_center_z;
  int i_t, j_t, k_t, i_t_r, i_t_l, j_t_l, j_t_r, k_t_l, k_t_r;
  int i,j,k;

  hsize_t start[4], start_buf[4], count[4], count_buf[4];
  hid_t vel_buf_space_id;
  herr_t status;

  double* vel_buf;
  double* vel_buf_read;
  int ic, jc, kc, a;
  int ic_t, jc_t, kc_t;
  unsigned int index, index_buf;
  double x,y,z,zl,zr,yl,yr,xl,xr;
  double dx,dy,dz;
  double norm, tmp;

  char err_mesg[4096];

  /* first get deltas of simulation grid
   * Coords contain cell-centered coords */
  if (*isize > 1) 
     dx = (iCoords[*isize-1] - iCoords[0]) / (*isize-1);
  else
     dx = 0.0;

  if (*jsize > 1)
     dy = (jCoords[*jsize-1] - jCoords[0]) / (*jsize-1);
  else
     dy = 0.0;

  if (*ksize > 1)
     dz = (kCoords[*ksize-1] - kCoords[0]) / (*ksize-1);
  else
     dz = 0.0;

  /* turbulence should be resolved with refinement level of max-1 */
  if ( ( (int) (dx / tf_delta[0]) > 2 ) ||
       ( (int) (dy / tf_delta[1]) > 2 ) ||
       ( (int) (dz / tf_delta[2]) > 2 ) ) {

     if (USE_COLLECTIVE)
        sim_turb_field_null_read_c();

     for (index = 0; index < (*isize)*(*jsize)*(*ksize); index++) {
        velx[index] = 0.0;
        vely[index] = 0.0;
        velz[index] = 0.0;
     }

     return;
  }

  /* find indices of turbulence grid to read in */
  /* first determine whether to use cell centered coords */
  /* left edge */
  min_i_t = (int) floor( ( iCoords[0] - dx/2.0 -
                           ( tf_corner[0] - tf_delta[0]/2.0 ) ) /
                         tf_delta[0] );
  /* right edge */
  max_i_t = (int) ceil( ( iCoords[0] + dx/2.0 -
                           ( tf_corner[0] + tf_delta[0]/2.0 ) ) /
                         tf_delta[0] );
  size_is_t = max_i_t - min_i_t + 1;
  /* if less than 2, then move to cell-centered indices */
  use_cell_center_x = size_is_t < 2;
  if ( use_cell_center_x ) {
     min_i_t = (int) floor( (iCoords[0] - tf_corner[0]) / tf_delta[0] );
     max_i_t = (int) floor( (iCoords[*isize-1] - tf_corner[0]) / tf_delta[0] )
             + 1;
  } else {
     max_i_t = (int) ceil( ( iCoords[*isize-1] + dx/2.0 -
                             ( tf_corner[0] + tf_delta[0]/2.0 ) ) /
                           tf_delta[0] );
  }
  size_is_t = max_i_t - min_i_t + 1;
  /* handle periodicity */
  while (max_i_t < 0) max_i_t += tf_size[0];
  while (min_i_t < 0) min_i_t += tf_size[0];
  max_i_t %= tf_size[0];
  min_i_t %= tf_size[0];
  on_edge_i = min_i_t > max_i_t;
  /* convert to unsmoothed grid to be read in */
  min_i_t *= s;
  max_i_t *= s;
  max_i_t += s - 1;
  size_i_t = s * size_is_t;

  /* left edge */
  min_j_t = (int) floor( ( jCoords[0] - dy/2.0 -
                           ( tf_corner[1] - tf_delta[1]/2.0 ) ) /
                         tf_delta[1] );
  /* right edge */
  max_j_t = (int) ceil( ( jCoords[0] + dy/2.0 -
                          ( tf_corner[1] + tf_delta[1]/2.0 ) ) /
                        tf_delta[1] );
  size_js_t = max_j_t - min_j_t + 1;
  /* if less than 2, then move to cell-centered indices */
  use_cell_center_y = size_js_t < 2;
  if ( use_cell_center_y ) {
     min_j_t = (int) floor( (jCoords[0] - tf_corner[1]) / tf_delta[1] );
     max_j_t = (int) floor( (jCoords[*jsize-1] - tf_corner[1]) / tf_delta[1] )
             + 1;
  } else {
     max_j_t = (int) ceil( ( jCoords[*jsize-1] + dy/2.0 -
                             ( tf_corner[1] + tf_delta[1]/2.0 ) ) /
                           tf_delta[1] );
  }
  size_js_t = max_j_t - min_j_t + 1;
  /* handle periodicity */
  while (max_j_t < 0) max_j_t += tf_size[1];
  while (min_j_t < 0) min_j_t += tf_size[1];
  max_j_t %= tf_size[1];
  min_j_t %= tf_size[1];
  on_edge_j = min_j_t > max_j_t;
  /* convert to unsmoothed grid to be read in */
  min_j_t *= s;
  max_j_t *= s;
  max_j_t += s - 1;
  size_j_t = s * size_js_t;

  /* left edge */
  min_k_t = (int) floor( ( kCoords[0] - dz/2.0 -
                           ( tf_corner[2] - tf_delta[2]/2.0 ) ) /
                         tf_delta[2] );
  /* right edge */
  max_k_t = (int) ceil( ( kCoords[0] + dz/2.0 -
                          ( tf_corner[2] + tf_delta[2]/2.0 ) ) /
                        tf_delta[2] );
  size_ks_t = max_k_t - min_k_t + 1;
  /* if less than 2, then move to cell-centered indices */
  use_cell_center_z = size_ks_t < 2;
  if ( use_cell_center_z ) {
     min_k_t = (int) floor( (kCoords[0] - tf_corner[2]) / tf_delta[2] );
     max_k_t = (int) floor( (kCoords[*ksize-1] - tf_corner[2]) / tf_delta[2] )
             + 1;
  } else {
     max_k_t = (int) ceil( ( kCoords[*ksize-1] + dz/2.0 -
                             ( tf_corner[2] + tf_delta[2]/2.0 ) ) /
                           tf_delta[2] );
  }
  size_ks_t = max_k_t - min_k_t + 1;
  /* handle periodicity */
  while (max_k_t < 0) max_k_t += tf_size[2];
  while (min_k_t < 0) min_k_t += tf_size[2];
  max_k_t %= tf_size[2];
  min_k_t %= tf_size[2];
  on_edge_k = min_k_t > max_k_t;
  /* convert to unsmoothed grid to be read in */
  min_k_t *= s;
  max_k_t *= s;
  max_k_t += s - 1;
  size_k_t = s * size_ks_t;

  /* read data */

  /* the data is ordered with the velocity component iterating the
   * fastest (so that values for a single cell are contiguous) and 
   * then with the x/i coordinate index incrementing next fastest, 
   * followed by y/j and finally z/k.
   * Note HDF5 convection is that the higher index dimensions
   * increment faster. */

  count_buf[0] = count[0] = (hsize_t) size_k_t;
  count_buf[1] = count[1] = (hsize_t) size_j_t;
  count_buf[2] = count[2] = (hsize_t) size_i_t;
  count_buf[3] = count[3] = (hsize_t) 3;

  /* our buffer is only the read block, 
   * so it needs its own dataspace descriptor */
  vel_buf_space_id = H5Screate_simple(4, count_buf, NULL);

  vel_buf = malloc(sizeof(double)*count_buf[0]*
                                  count_buf[1]*
                                  count_buf[2]*
                                  count_buf[3]);
  if (vel_buf == NULL)
    Driver_abortFlashC("error allocating memory for turbulence field\n");

  start[0] = (hsize_t) min_k_t;
  start[1] = (hsize_t) min_j_t;
  start[2] = (hsize_t) min_i_t;
  start[3] = (hsize_t) 0;

  start_buf[0] = (hsize_t) 0;
  start_buf[1] = (hsize_t) 0;
  start_buf[2] = (hsize_t) 0;
  start_buf[3] = (hsize_t) 0;

  /* if region does not overlap an edge in the turbulence data,
   * just read it all at once.
   * if it does overlap an edge, read is done in multiple stages
   * starting with the portion on the "low" side of the edge. */
  if ( on_edge_k )
     count_buf[0] = count[0] = (hsize_t) (s*tf_size[2]-min_k_t);

  if ( on_edge_j )
     count_buf[1] = count[1] = (hsize_t) (s*tf_size[1]-min_j_t);

  if ( on_edge_i )
     count_buf[2] = count[2] = (hsize_t) (s*tf_size[0]-min_i_t);

  /* select what we want to read */
  status = H5Sselect_hyperslab( turb_vel_dspace_id, H5S_SELECT_SET,
                                start, NULL, count, NULL);
  if ( status == H5_FAIL )
     Driver_abortFlashC("H5Sselect_hyperslab: file\n");

  /* select where we want to write */
  status = H5Sselect_hyperslab( vel_buf_space_id, H5S_SELECT_SET,
                                start_buf, NULL, count_buf, NULL); 
  if ( status == H5_FAIL )
     Driver_abortFlashC("H5Sselect_hyperslab: buffer\n");

  /* perform read */
  if ( USE_COLLECTIVE ) MPI_Barrier( MPI_COMM_WORLD );
  status = H5Dread( turb_vel_dset_id, H5T_NATIVE_DOUBLE,
                    vel_buf_space_id, turb_vel_dspace_id, 
                    dxfer_template, vel_buf);
  if ( status == H5_FAIL )
     Driver_abortFlashC("H5Dread: error reading from turbulence field\n");

  if ( on_edge_k ) {
     start[0] = (hsize_t) 0;      // on an edge
     start_buf[0] = count_buf[0]; // start where we left off 
     count_buf[0] = count[0] = (hsize_t) (max_k_t+1);
     status = H5Sselect_hyperslab( turb_vel_dspace_id, 
                                   H5S_SELECT_SET, start, 
                                   NULL, count, NULL);
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_hyperslab: file\n");
     status = H5Sselect_hyperslab( vel_buf_space_id, H5S_SELECT_SET,
                                   start_buf, NULL, count_buf, NULL); 
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_hyperslab: buffer\n");
     vel_buf_read = vel_buf;

     /* reset to original values for possible additional slices */
     start[0] = (hsize_t) min_k_t;
     start_buf[0] = (hsize_t) 0;
     count[0] = count_buf[0] = (hsize_t) (s*tf_size[2]-min_k_t);
  } else {
     status = H5Sselect_none( turb_vel_dspace_id );
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_none: failed\n");
     status = H5Sselect_none( vel_buf_space_id );
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_none: failed\n");
     vel_buf_read = NULL;
  }
  if ( USE_COLLECTIVE ) MPI_Barrier( MPI_COMM_WORLD );
  status = H5Dread( turb_vel_dset_id, H5T_NATIVE_DOUBLE,
                    vel_buf_space_id, turb_vel_dspace_id, 
                    dxfer_template, vel_buf_read);
  if ( status == H5_FAIL )
     Driver_abortFlashC("H5Dread: error reading from turbulence field\n");

  if ( on_edge_j ) {
     start[1] = (hsize_t) 0;
     start_buf[1] = count_buf[1];
     count_buf[1] = count[1] = (hsize_t) (max_j_t+1);
     status = H5Sselect_hyperslab( turb_vel_dspace_id, H5S_SELECT_SET,
                                   start, NULL, count, NULL);
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_hyperslab: file\n");
     status = H5Sselect_hyperslab( vel_buf_space_id, H5S_SELECT_SET,
                                   start_buf, NULL, count_buf, NULL); 
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_hyperslab: buffer\n");
     vel_buf_read = vel_buf;

     /* reset to original values for possible additional slices */
     start[1] = (hsize_t) min_j_t;
     start_buf[1] = (hsize_t) 0;
     count[1] = count_buf[1] = (hsize_t) (s*tf_size[1]-min_j_t);
  } else {
     status = H5Sselect_none( turb_vel_dspace_id );
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_none: failed\n");
     status = H5Sselect_none( vel_buf_space_id );
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_none: failed\n");
     vel_buf_read = NULL;
  }
  if ( USE_COLLECTIVE ) MPI_Barrier( MPI_COMM_WORLD );
  status = H5Dread( turb_vel_dset_id, H5T_NATIVE_DOUBLE,
                    vel_buf_space_id, turb_vel_dspace_id, 
                    dxfer_template, vel_buf_read);
  if ( status == H5_FAIL )
     Driver_abortFlashC("H5Dread: error reading from turbulence field\n");

  if ( on_edge_i ) {
     start[2] = (hsize_t) 0;
     start_buf[2] = count_buf[2];
     count_buf[2] = count[2] = (hsize_t) (max_i_t+1);
     status = H5Sselect_hyperslab( turb_vel_dspace_id, H5S_SELECT_SET,
                                   start, NULL, count, NULL);
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_hyperslab: file\n");
     status = H5Sselect_hyperslab( vel_buf_space_id, H5S_SELECT_SET,
                                   start_buf, NULL, count_buf, NULL); 
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_hyperslab: buffer\n");
     vel_buf_read = vel_buf;

     /* reset to original values for possible additional slices */
     start[2] = (hsize_t) min_i_t;
     start_buf[2] = (hsize_t) 0;
     count[2] = count_buf[2] = (hsize_t) (s*tf_size[0]-min_i_t);
  } else {
     status = H5Sselect_none( turb_vel_dspace_id );
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_none: failed\n");
     status = H5Sselect_none( vel_buf_space_id );
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_none: failed\n");
     vel_buf_read = NULL;
  }
  if ( USE_COLLECTIVE ) MPI_Barrier( MPI_COMM_WORLD );
  status = H5Dread( turb_vel_dset_id, H5T_NATIVE_DOUBLE,
                    vel_buf_space_id, turb_vel_dspace_id,
                    dxfer_template, vel_buf_read);
  if ( status == H5_FAIL )
     Driver_abortFlashC("H5Dread: error reading from turbulence field\n");

  if ( on_edge_k && on_edge_j ) {
     start[0] = (hsize_t) 0;
     start_buf[0] = count_buf[0];
     count_buf[0] = count[0] = (hsize_t) (max_k_t+1);
     start[1] = (hsize_t) 0;
     start_buf[1] = count_buf[1];
     count_buf[1] = count[1] = (hsize_t) (max_j_t+1);
     status = H5Sselect_hyperslab( turb_vel_dspace_id, H5S_SELECT_SET,
                                   start, NULL, count, NULL);
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_hyperslab: file\n");
     status = H5Sselect_hyperslab( vel_buf_space_id, H5S_SELECT_SET,
                                   start_buf, NULL, count_buf, NULL); 
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_hyperslab: buffer\n");
     vel_buf_read = vel_buf;

     /* reset to original values for possible additional slices */
     start[0] = (hsize_t) min_k_t;
     start_buf[0] = (hsize_t) 0;
     count[0] = count_buf[0] = (hsize_t) (s*tf_size[2]-min_k_t);
     start[1] = (hsize_t) min_j_t;
     start_buf[1] = (hsize_t) 0;
     count[1] = count_buf[1] = (hsize_t) (s*tf_size[1]-min_j_t);
  } else {
     status = H5Sselect_none( turb_vel_dspace_id );
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_none: failed\n");
     status = H5Sselect_none( vel_buf_space_id );
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_none: failed\n");
     vel_buf_read = NULL;
  }
  if ( USE_COLLECTIVE ) MPI_Barrier( MPI_COMM_WORLD );
  status = H5Dread( turb_vel_dset_id, H5T_NATIVE_DOUBLE,
                    vel_buf_space_id, turb_vel_dspace_id,
                    dxfer_template, vel_buf_read);
  if ( status == H5_FAIL )
     Driver_abortFlashC("H5Dread: error reading from turbulence field\n");

  if ( on_edge_k && on_edge_i ) {
     start[0] = (hsize_t) 0;
     start_buf[0] = count_buf[0];
     count_buf[0] = count[0] = (hsize_t) (max_k_t+1);
     start[2] = (hsize_t) 0;
     start_buf[2] = count_buf[2];
     count_buf[2] = count[2] = (hsize_t) (max_i_t+1);
     status = H5Sselect_hyperslab( turb_vel_dspace_id, H5S_SELECT_SET,
                                   start, NULL, count, NULL);
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_hyperslab: file\n");
     status = H5Sselect_hyperslab( vel_buf_space_id, H5S_SELECT_SET,
                                   start_buf, NULL, count_buf, NULL); 
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_hyperslab: buffer\n");
     vel_buf_read = vel_buf;

     /* reset to original values for possible additional slices */
     start[0] = (hsize_t) min_k_t;
     start_buf[0] = (hsize_t) 0;
     count[0] = count_buf[0] = (hsize_t) (s*tf_size[2]-min_k_t);
     start[2] = (hsize_t) min_i_t;
     start_buf[2] = (hsize_t) 0;
     count[2] = count_buf[2] = (hsize_t) (s*tf_size[0]-min_i_t);
  } else {
     status = H5Sselect_none( turb_vel_dspace_id );
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_none: failed\n");
     status = H5Sselect_none( vel_buf_space_id );
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_none: failed\n");
     vel_buf_read = NULL;
  }
  if ( USE_COLLECTIVE ) MPI_Barrier( MPI_COMM_WORLD );
  status = H5Dread( turb_vel_dset_id, H5T_NATIVE_DOUBLE,
                    vel_buf_space_id, turb_vel_dspace_id, 
                    dxfer_template, vel_buf_read);
  if ( status == H5_FAIL )
     Driver_abortFlashC("H5Dread: error reading from turbulence field\n");

  if ( on_edge_j && on_edge_i ) {
     start[1] = (hsize_t) 0;
     start_buf[1] = count_buf[1];
     count_buf[1] = count[1] = (hsize_t) (max_j_t+1);
     start[2] = (hsize_t) 0;
     start_buf[2] = count_buf[2];
     count_buf[2] = count[2] = (hsize_t) (max_i_t+1);
     status = H5Sselect_hyperslab( turb_vel_dspace_id, H5S_SELECT_SET,
                                   start, NULL, count, NULL);
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_hyperslab: file\n");
     status = H5Sselect_hyperslab( vel_buf_space_id, H5S_SELECT_SET,
                                   start_buf, NULL, count_buf, NULL); 
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_hyperslab: buffer\n");
     vel_buf_read = vel_buf;

     /* reset to original values for possible additional slices */
     start[1] = (hsize_t) min_j_t;
     start_buf[1] = (hsize_t) 0;
     count[1] = count_buf[1] = (hsize_t) (s*tf_size[1]-min_j_t);
     start[2] = (hsize_t) min_i_t;
     start_buf[2] = (hsize_t) 0;
     count[2] = count_buf[2] = (hsize_t) (s*tf_size[0]-min_i_t);
  } else {
     status = H5Sselect_none( turb_vel_dspace_id );
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_none: failed\n");
     status = H5Sselect_none( vel_buf_space_id );
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_none: failed\n");
     vel_buf_read = NULL;
  }
  if ( USE_COLLECTIVE ) MPI_Barrier( MPI_COMM_WORLD );
  status = H5Dread( turb_vel_dset_id, H5T_NATIVE_DOUBLE,
                    vel_buf_space_id, turb_vel_dspace_id,
                    dxfer_template, vel_buf_read);
  if ( status == H5_FAIL )
     Driver_abortFlashC("H5Dread: error reading from turbulence field\n");

  if ( on_edge_k && on_edge_j && on_edge_i ) {
     start[0] = (hsize_t) 0;
     start_buf[0] = count_buf[0];
     count_buf[0] = count[0] = (hsize_t) (max_k_t+1);
     start[1] = (hsize_t) 0;
     start_buf[1] = count_buf[1];
     count_buf[1] = count[1] = (hsize_t) (max_j_t+1);
     start[2] = (hsize_t) 0;
     start_buf[2] = count_buf[2];
     count_buf[2] = count[2] = (hsize_t) (max_i_t+1);
     status = H5Sselect_hyperslab( turb_vel_dspace_id, H5S_SELECT_SET,
                                   start, NULL, count, NULL);
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_hyperslab: file\n");
     status = H5Sselect_hyperslab( vel_buf_space_id, H5S_SELECT_SET,
                                  start_buf, NULL, count_buf, NULL); 
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_hyperslab: buffer\n");
     vel_buf_read = vel_buf;
  } else {
     status = H5Sselect_none( turb_vel_dspace_id );
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_none: failed\n");
     status = H5Sselect_none( vel_buf_space_id );
     if ( status == H5_FAIL )
        Driver_abortFlashC("H5Sselect_none: failed\n");
     vel_buf_read = NULL;
  }
  if ( USE_COLLECTIVE ) MPI_Barrier( MPI_COMM_WORLD );
  status = H5Dread( turb_vel_dset_id, H5T_NATIVE_DOUBLE,
                    vel_buf_space_id, turb_vel_dspace_id,
                    dxfer_template, vel_buf_read);
  if ( status == H5_FAIL )
     Driver_abortFlashC("H5Dread: error reading from turbulence field\n");

  status = H5Sclose(vel_buf_space_id);
  if ( status == H5_FAIL )
     Driver_abortFlashC("H5Sclose: buffer space close failed\n");
  /* turb_vel_dspace_id will be reused */

  /* now vel_buf contains unsmoothed velocity field
   * smooth it. */
  for (k_t = 0, kc_t = 0; kc_t < size_k_t; kc_t+=s, k_t++) {
     for (j_t = 0, jc_t = 0; jc_t < size_j_t; jc_t+=s, j_t++) {
        for (i_t = 0, ic_t = 0; ic_t < size_i_t; ic_t+=s, i_t++) {
           for (a = 0; a < 3; a++) {
              tmp = 0.0;
              for (kc = kc_t; kc < kc_t+s; kc++) {
                 for (jc = jc_t; jc < jc_t+s; jc++) {
                    for (ic = ic_t; ic < ic_t+s; ic++) {
                       index_buf = kc*size_j_t*size_i_t*3 +
                                   jc*size_i_t*3 + ic*3 + a;
                       tmp += vel_buf[index_buf];
                    }
                 }
              }
              index_buf = k_t*size_j_t*size_i_t*3 +
                          j_t*size_i_t*3 + i_t*3 + a;
              vel_buf[index_buf] = tmp / pow(s,3);
           }
        }
     }
  }

  /* convert min indices back to smoothed turb field */
  min_i_t /= s;
  min_j_t /= s;
  min_k_t /= s;
  max_i_t = size_is_t + min_i_t - 1;
  max_j_t = size_js_t + min_j_t - 1;
  max_k_t = size_ks_t + min_k_t - 1;

  /* find index of grid interval in turbulence field grid
   * containting this point. note that the grid interval indices 
   * are different from the cell indices */
  /* line up cell edges of smoothed turbulence field */
  for (k=0; k<*ksize; k++) {
     if ( use_cell_center_z ) {
        k_t_l = (int) floor( (kCoords[k] - tf_corner[2]) / tf_delta[2] );
        k_t_r = k_t_l + 1;
        zl = (kCoords[k] - tf_corner[2]) / tf_delta[2] - k_t_l; 
        zr = 1.0 - zl;
     } else {
        // find index of turbulence field for left edge
        k_t_l = (int) floor( ( kCoords[k] - dz/2.0 - 
                               ( tf_corner[2] - tf_delta[2]/2.0 ) ) / 
                             tf_delta[2] );
        // find index for right edge
        k_t_r = (int)  ceil( ( kCoords[k] + dz/2.0 - 
                               ( tf_corner[2] + tf_delta[2]/2.0 ) ) /
                             tf_delta[2] );
        zl = ( kCoords[k] - dz/2.0 - ( tf_corner[2] - tf_delta[2]/2.0 ) ) / 
             tf_delta[2] - k_t_l;
        zr = k_t_r -
             ( kCoords[k] + dz/2.0 - ( tf_corner[2] + tf_delta[2]/2.0 ) ) / 
             tf_delta[2];
     }
     size_ks_t = k_t_r - k_t_l + 1;

     // convert to read in indices
     k_t_l -= min_k_t;
     while (k_t_l < 0) k_t_l += tf_size[2];
     k_t_l %= tf_size[2];
     k_t_r = size_ks_t + k_t_l - 1;
     if ( ( size_ks_t < 2 )   ||
          ( k_t_l < 0 ) ||
          ( k_t_r >= size_k_t/s ) ) {
        sprintf(err_mesg,"ktl=%d, ktr=%d, sizek=%d\n",
                k_t_l, k_t_r, size_k_t/s);
        Driver_abortFlashC(err_mesg);
     }
     for (j=0; j<*jsize; j++) {
        if ( use_cell_center_y ) {
           j_t_l = (int) floor( ( jCoords[j] - tf_corner[1] ) / tf_delta[1] );
           j_t_r = j_t_l + 1;
           yl = ( jCoords[j] - tf_corner[1] ) / tf_delta[1] - j_t_l;
           yr = 1.0 - yl;
        } else {
           j_t_l = (int) floor( ( jCoords[j] - dy/2.0 - 
                                  ( tf_corner[1] - tf_delta[1]/2.0 ) ) / 
                                tf_delta[1] );
           j_t_r = (int)  ceil( ( jCoords[j] + dy/2.0 - 
                                  ( tf_corner[1] + tf_delta[1]/2.0 ) ) /
                                tf_delta[1] );
           yl = ( jCoords[j] - dy/2.0 - ( tf_corner[1] - tf_delta[1]/2.0 ) ) /
                tf_delta[1] - j_t_l;
           yr = j_t_r -
                ( jCoords[j] + dy/2.0 - ( tf_corner[1] + tf_delta[1]/2.0 ) ) /
                tf_delta[1];
        }
        size_js_t = j_t_r - j_t_l + 1;

        j_t_l -= min_j_t;
        while (j_t_l < 0) j_t_l += tf_size[1];
        j_t_l %= tf_size[1];
        j_t_r = size_js_t + j_t_l - 1;
        if ( ( size_js_t < 2 )   ||
             ( j_t_l < 0 ) ||
             ( j_t_r >= size_j_t/s ) ) {
           sprintf(err_mesg,"jtl=%d, jtr=%d, sizej=%d\n",
                   j_t_l, j_t_r, size_j_t/s );
           Driver_abortFlashC(err_mesg);
        }
        for (i=0; i<*isize; i++) {
           if ( use_cell_center_x ) {
              i_t_l = (int) floor( ( iCoords[i] - tf_corner[0] ) / tf_delta[0] );
              i_t_r = i_t_l + 1;
              xl = ( iCoords[i] - tf_corner[0] ) / tf_delta[0] - i_t_l;
              xr = 1.0 - xl;
           } else {
              i_t_l = (int) floor( ( iCoords[i] - dx/2.0 -
                                     ( tf_corner[0] - tf_delta[0]/2.0 ) ) /
                                   tf_delta[0] );
              i_t_r = (int)  ceil( ( iCoords[i] + dx/2.0 -
                                     ( tf_corner[0] + tf_delta[0]/2.0 ) ) /
                                   tf_delta[0] );
              xl = ( iCoords[i] - dx/2.0 - ( tf_corner[0] - tf_delta[0]/2.0 ) ) /
                   tf_delta[0] - i_t_l;
              xr = i_t_r -
                   ( iCoords[i] + dx/2.0 - ( tf_corner[0] + tf_delta[0]/2.0 ) ) /
                   tf_delta[0];
           }
           size_is_t = i_t_r - i_t_l + 1;

           i_t_l -= min_i_t;
           while (i_t_l < 0) i_t_l += tf_size[0];
           i_t_l %= tf_size[0];
           i_t_r = size_is_t + i_t_l - 1;
           if ( ( size_is_t < 2 )   ||
                ( i_t_l < 0 ) ||
                ( i_t_r >= size_i_t/s ) ) {
              sprintf(err_mesg,"itl=%d, itr=%d, sizei=%d\n",
                      i_t_l, i_t_r, size_i_t/s );
              Driver_abortFlashC(err_mesg);
           }

           /* index in output arrays */
           index = k*(*isize)*(*jsize) + j*(*isize) + i;

           /* initialize values */
           velx[index] = 0.0;
           vely[index] = 0.0;
           velz[index] = 0.0;

           /* calculate value for return */
           for (k_t = k_t_l; k_t <= k_t_r; k_t++) {
              if (k_t == k_t_l) z = 1.0-zl;
              else if (k_t == k_t_r) z = 1.0-zr;
              else z = 1.0; 
              for (j_t = j_t_l; j_t <= j_t_r; j_t++) {
                 if (j_t == j_t_l) y = 1.0-yl;
                 else if (j_t == j_t_r) y = 1.0-yr;
                 else y = 1.0;
                 for (i_t = i_t_l; i_t <= i_t_r; i_t++) {
                     if (i_t == i_t_l) x = 1.0-xl;
                     else if (i_t == i_t_r) x = 1.0-xr;
                     else x = 1.0;

                     // indexing still in read-in space
                     index_buf = k_t*size_j_t*size_i_t*3
                               + j_t*size_i_t*3
                               + i_t*3;
                     velx[index] += z * y * x * vel_buf[index_buf];
                     vely[index] += z * y * x * vel_buf[index_buf+1];
                     velz[index] += z * y * x * vel_buf[index_buf+2];


                 }
              }
           }
           /* norm = ( 1 - xr + 1 - xl + size_i_t/s - 2 ) * ( j ) * ( k ) 
            * product of the sum of the weights in each direction */
           norm = ( (double) size_is_t - xl - xr )
                * ( (double) size_js_t - yl - yr )
                * ( (double) size_ks_t - zl - zr );
           velx[index] /= norm;
           vely[index] /= norm;
           velz[index] /= norm;

        }
     }
  }

  if (vel_buf!=NULL) free(vel_buf);

  return;
}

/* teardown function closes the data file */

void FTOC(sim_turb_field_teardown)(void) {
  herr_t ret;

#ifdef TURB_PARALLELIO
  MPI_Info_free(&info);
#endif
  ret = H5Pclose(dxfer_template);
  if ( ret == H5_FAIL )
     Driver_abortFlashC("H5Pclose: dataset transfer close failed\n");
  ret = H5Sclose(turb_vel_dspace_id);
  if ( ret == H5_FAIL )
     Driver_abortFlashC("H5Sclose: dataspace close failed\n");
  ret = H5Dclose(turb_vel_dset_id);
  if ( ret == H5_FAIL )
     Driver_abortFlashC("H5Dclose: dataset close failed\n");
  ret = H5Fclose(turb_vel_file_id);
  if ( ret == H5_FAIL )
     Driver_abortFlashC("H5Fclose: file close failed\n");

}
