//******************************************************************************
//Name:  field.hpp                                                             *
//                                                                             *
//Purpose:  header file for defining the field class                           *
//                                                                             *
//Modification History:  10/24/98 - added modification history field           *
//                    :  10/24/98 - changed SetTensor to use int indices       *
//                    :             rather than a double position array        *
//                    :  11/12/98 - promoted temporary variables in            *
//                    :             GetInterpTensor to protected member data   *
//                    :             in interest of performance                 *
//                    :  12/10/98 - implemented the user-defined grid smoothing*
//                                  radii
// *****************************************************************************
#ifndef _FIELDHPP

#include "tensor.hpp"

#define X_MIN 0
#define X_MAX 1
#define Y_MIN 2
#define Y_MAX 3
#define Z_MIN 4
#define Z_MAX 5

#define X_ID 0
#define Y_ID 1
#define Z_ID 2

#define SPATIAL_TOL  1e-10

#define PARTICLE_SMOOTHING 0
#define GRID_SMOOTHING     1

#define CWM_SMOOTHER 0
#define DEFAULT_SMOOTHING 9999
#define SMOOTHING_TOL 1e-7


class field : public tensor
{
  protected:

    //parameters related to the time of the current spatial field slice
    double  current_time;
    int	   is_time_set;

    //parameters related to the spatial grid structure of the field
    double  spatial_extent[6];
    int     is_pos_set;
    double  Dist(double *position);
    double  Deltas[3];
    double  Delta_mag;

    //parameters related to the smoothing performed on the interpolated tensors
    //associated with a position
    int     particle_smoothing_method;
    int     grid_smoothing_method;
    int     out_of_range;
    int     is_periodic;
    int     grid_radius;
    int     grid_lower_bound;
    int     grid_upper_bound;
    double  smoothing_length;
    double  Smoother(int type);
    tensor  Smoother_Deriv(int type);
    double  reference_pos[3];
    double  current_pos[3];
    double  pos_diff[3];
    double  y;
    double  normalization;
    struct  series_parms {int order; field **list;};
    struct  series_parms series;
    int     is_series_set;
    int     is_deriv_field;
    int     is_dep_field;
	//tensors defined in for use in GetInterpTensor herein declared in
	//the interest of speed
    tensor  *tmp;
    tensor  Wm, a2, sum_Wm, sum_a2, sum_da1, d_sum_a2;
    tensor d, ret_dummy, temp0, temp1, temp2;



  public:

    field(void);
    field(int num_indices, int range1, ...);

    void    SetSpatialExtent(int id, double value);
    void    SetSpatialExtent(double *value);
    void    ResizeField(int *num_grid_pts);
    void    GetIndices(double *position, int *indices);
    void    GetIndices(double *position, double *indices);
    void    GetPosition(int *indices, double *position);
    tensor  GetPosition(int *indices);
    void    SetSmoother(int id, int smoothing_method);
    void    SetGridRadius(int grid_rad);
    void    SetPeriodic(int choice);
    tensor  GetTensor(double *position);
    tensor  GetTensor(const tensor &pos);
    tensor  GetTensor(int *indices);
    void    SetTensor(double *position, const tensor &rval);
    void    SetTensor(const tensor &pos, const tensor &rval);
    void    SetTensor(int *indices, const tensor &rval);
    void    SetCurrentTime(double time);
    double  GetCurrentTime(void);
    void    SetSeriesOrder(int order);
    void    SetSeriesList(field *list[]);
    void    SetDerivField(int type);
    void    SetFieldDep(int type);
    tensor  GetInterpTensor(tensor &pos);
    tensor  GetInterpTensor(double *position);
	 void    SetInterpTensor(double *position, const tensor &rval);
    int     GetNumIndices(int id);
};

#define _FIELDHPP
#endif



