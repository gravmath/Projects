//**************************************************************
//crank nicholson driver
//
//
//**************************************************************
#include "crnk_nch.h"

void main(void);

tensor anal_schw(double time, const tensor &state, struct sim_parms *parms);
tensor grid_schw(double time, const tensor &state, struct sim_parms *parms);

void   std_schw(struct sim_parms *parms, const tensor &pos,
   			    tensor &lapse,           tensor &shift,     tensor &h,
				    tensor &d_lapse,         tensor &d_shift,   tensor &d_h);

void   schw_higher_order(struct sim_parms *parms, const tensor &pos,
                         tensor &dd_lapse,        tensor &dd_shift,
                         tensor &dd_h);

void load_field(struct sim_parms *parms);

void inspect_smoothing(struct sim_parms *parms);

void initialize(struct sim_parms *parms);
void output(void);
void GetLine(FILE *file, char *string);

#define STATE_SIZE 6
#define ANAL_SCHW  0
#define GRID_SCHW  1
#define INSPECT    2

//global variables to keep track of whether the dd_lapse, etc. terms are needed
int max_basic_order = 0;
int max_deriv_order = 0;
int max_order = 0;



void main(void)
{
   tensor (*rhs)(double time, const tensor &state, struct sim_parms *parms);
   int method;
   struct sim_parms parms;

   parms.M = 0.0;
   parms.m = 0.0;
   parms.H = 0.0;
   parms.part_sm_meth = 1;
   parms.grid_sm_meth = 1;


   cout << "Welcome to the geodesic playroom" <<endl;
   cout << endl;
   cout << "Please enter:  0 to use analytic Schwarzschild in cartesian coords"
        << endl;
   cout << "            :  1 to use gridded Schwarzschild in cartesian coords"
        << endl;
   cout << "            :  2 to inspect the Schwarzschild fields             "
        << endl;
   cin >> method;
   switch(method)
   {
     case ANAL_SCHW: rhs = anal_schw;                               break;
     case GRID_SCHW: rhs = grid_schw;    initialize(&parms);        break;
     case INSPECT:   initialize(&parms); inspect_smoothing(&parms); break;
   }


   char out_file[80];
   FILE *out;
   cout << "Enter the output filename" << endl;
   cin >> out_file;
   out = fopen(out_file,"w");


   double cur_time = 0.0;
   double final_time = 10.0;
   double time_step = 0.1;
   cout << "Please enter the central body mass and the particle mass" << endl;
   cin >> parms.M >> parms.m;
   cout << "Enter the initial time, final time,and the time step" << endl;
   cin >> cur_time >> final_time >> time_step;
   int num_steps;
   num_steps = (final_time - cur_time)/time_step;

   int frequency;
   int counter = 0;
   cout <<"Enter the output frequency" << endl;
   cin >> frequency;

   double state[STATE_SIZE];
   double new_state[STATE_SIZE];
   int i;
   cout << "Enter the state:  position first followed by the momenta" << endl;
   for(i = 0; i < STATE_SIZE; i++)  cin >> state[i];

   fprintf(out,"%19.16g ", cur_time);
   for(i = 0; i < STATE_SIZE; i++) fprintf(out,"%19.16g ", state[i]);
   fprintf(out,"%19.16g ", parms.H);

   while( cur_time <= final_time )
    {
     crank_nicholson(cur_time, time_step, STATE_SIZE, state, new_state, &parms, rhs);
	  counter++;
	  for(i = 0; i < STATE_SIZE; i++) state[i] = new_state[i];
	  cur_time = cur_time + time_step;
	  if( !(counter%frequency) )
	  {
        fprintf(out,"\n%6.4g ", cur_time);
        for(i = 0; i < STATE_SIZE; i++) fprintf(out,"%19.16g ", state[i]);
        fprintf(out,"%19.16g", parms.H);
	    printf(".");
	  }
    }

    fclose(out);

}

/******************************************************************************/
/*     Function to calculate the time derivative of the state given the       */
/*     standard Schwarzschild metric                                          */
/******************************************************************************/
tensor anal_schw(double time, const tensor &state, struct sim_parms *parms)
{

  tensor lapse(1,1),   shift(1,3),     h(2,3,3);
  tensor d_lapse(1,3), d_shift(2,3,3), d_h(3,3,3,3);
  tensor temp;
  int i;

  lapse.SetName("lapse");
  shift.SetName("shift");
  h.SetName("h");
  d_lapse.SetName("d_lapse");
  d_shift.SetName("d_shift");
  d_h.SetName("d_h");
  temp.SetName("temp");

  //Pack the position
  tensor pos;
  pos.SetName("pos");
  for( i = 0; i < 3; i++) pos.Set(state.Val(i),i);

  //Pack the momentum
  tensor momentum;
  momentum.SetName("momentum");
  for( i = 3; i < 6; i++) momentum.Set(state.Val(i),i-3);

  //Get the metric and its first derivatives
  std_schw(parms, pos, lapse, shift, h, d_lapse, d_shift, d_h);

  //Calculate the forces
  tensor vel;
  double mu;

  //form mu
  temp <= (h.Contract(momentum,1,0)).Contract(momentum,0,0);
  mu = sqrt(parms->m*parms->m + temp.Val(0));

  //Calculate the energy
  parms->H = -(shift.Contract(momentum,0,0)).Val(0) + lapse.Val(0)*mu;

  //calculate the position updates
  vel <= -1.0*shift + (1.0/mu)*lapse.Val(0)*h.Contract(momentum,1,0);

  //calculate the momentum updates
  tensor f;
  f <= d_shift.Contract(momentum,0,0) - mu*d_lapse -
       1.0/(2.0*mu)*lapse.Val(0)*
       (d_h.Contract(momentum,1,0)).Contract(momentum,0,0);

  //repack temp
  temp <= state;
  for( i = 0; i < 3; i++)
   {
     temp.Set(vel.Val(i),i);
     temp.Set(f.Val(i),i+3);
   }

  return temp;
}

/******************************************************************************/
/*     Function to calculate the time derivative of the state given the       */
/*     standard Schwarzschild metric on a grid                                */
/******************************************************************************/
tensor grid_schw(double time, const tensor &state, struct sim_parms *parms)
{
  tensor lapse(1,1),   shift(1,3),     h(2,3,3);
  tensor d_lapse(1,3), d_shift(2,3,3), d_h(3,3,3,3);
  tensor temp;
  int i;
  static first_time = TRUE;

  lapse.SetName("lapse");
  shift.SetName("shift");
  h.SetName("h");
  d_lapse.SetName("d_lapse");
  d_shift.SetName("d_shift");
  d_h.SetName("d_h");
  temp.SetName("temp");


  if( first_time == TRUE )
   {
      //load the values of the various fields against their known
      //Schwarzshild values
      load_field(parms);
      first_time = FALSE;
   }

  //Pack the position
  tensor pos;
  pos.SetName("pos");
  for( i = 0; i < 3; i++) pos.Set(state.Val(i),i);

  //Pack the momentum
  tensor momentum;
  momentum.SetName("momentum");
  for( i = 3; i < 6; i++) momentum.Set(state.Val(i),i-3);

  //Get the tensors at the position
  lapse   <= parms->  lapse.GetInterpTensor(pos);
  shift   <= parms->  shift.GetInterpTensor(pos);
  h       <= parms->      h.GetInterpTensor(pos);
  d_lapse <= parms->d_lapse.GetInterpTensor(pos);
  d_shift <= parms->d_shift.GetInterpTensor(pos);
  d_h     <= parms->    d_h.GetInterpTensor(pos);

  //Calculate the forces
  tensor vel;
  double mu;

  //form mu
  temp <= (h.Contract(momentum,1,0)).Contract(momentum,0,0);
  mu = sqrt(parms->m*parms->m + temp.Val(0));

  //Calculate the energy
  parms->H = (shift.Contract(momentum,0,0)).Val(0) + lapse.Val(0)*mu;

  //calculate the position updates
  vel <= -1.0*shift + (1.0/mu)*lapse.Val(0)*h.Contract(momentum,1,0);

  //calculate the momentum updates
  tensor f;
  f <= d_shift.Contract(momentum,0,0) - mu*d_lapse -
       1.0/(2.0*mu)*lapse.Val(0)*
       (d_h.Contract(momentum,1,0)).Contract(momentum,0,0);

  //repack temp
  temp <= state;
  for( i = 0; i < 3; i++)
   {
     temp.Set(vel.Val(i),i);
     temp.Set(f.Val(i),i+3);
   }

  return temp;
}

/******************************************************************************/
/*     Function to calculate the metric and its first derivatives for the     */
/*     standard Schwarzschild metric                                          */
/******************************************************************************/
void   std_schw(struct sim_parms *parms, const tensor &pos,
   			    tensor &lapse,           tensor &shift,     tensor &h,
				    tensor &d_lapse,         tensor &d_shift,   tensor &d_h)
{
  int i, j, k;
  double M = parms->M;
  double m = parms->m;
  double sniggle = 1e-4;
  tensor temp;
  int *range;

  double r = 0.0;
  for( i = 0; i < 3; i++)r += pos.Val(i)*pos.Val(i);
  r = sqrt(r);

  //If within a sniggle of the event horizon zero-out the fields
  if( r <= 2*M + sniggle )
  {
     lapse   <= 0.0*lapse;
     shift   <= 0.0*shift;
     h       <= 0.0*h;
     d_lapse <= 0.0*d_lapse;
     d_shift <= 0.0*d_shift;
     d_h     <= 0.0*d_h;
	  return;
  }

  //Calculate the lapse
  lapse.Set(sqrt(1 - 2*M/r),1);

  //Calculate the shift
     //do nothing since zero

  //Calculate the spatial metric
  tensor delta(2,3,3);
  delta.SetName("delta");
  for(i = 0; i < 3; i++) delta.Set(1.0,i,i);
  h <= (-2.0*M/pow(r,3))*pos*pos + delta;

  //Calculate the gradient of the lapse
  d_lapse <= (M/( pow(r,3.0)*sqrt(1-2*M/r) ) ) * pos;

  //Calculate the gradient of the shift
    //do nothing since zero

  //Calculate the gradient of the netric
  d_h <= ( 3.0/pow(r,3.0) )*pos*pos*pos + (-1.0/r)*pos*delta;
  range = d_h.Ranges();
  temp.Resize(d_h.NumIndices(),range);
  delete [] range;

  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      for(k = 0; k < 3; k++)
         temp.Set(-1.0/r*delta.Val(i,k)*pos.Val(j),i,j,k);
  d_h <= (2.0*M/pow(r,2.0))* (d_h + temp);
}

/******************************************************************************/
/*     Function to calculate the second derivatives of the metric for the     */
/*     standard Schwarzschild metric                                          */
/******************************************************************************/
void   schw_higher_order(struct sim_parms *parms, const tensor &pos,
                         tensor &dd_lapse,        tensor &dd_shift,
                         tensor &dd_h)
{
  int i, j, k, l;
  double M = parms->M;
  double m = parms->m;
  double sniggle = 1e-4;
  tensor temp;

  double r = 0.0;
  for( i = 0; i < 3; i++)r += pos.Val(i)*pos.Val(i);
  r = sqrt(r);

  if( r <= 2*M + sniggle )
  {
     dd_lapse = 0.0*dd_lapse;
     dd_shift = 0.0*dd_shift;
     dd_h     = 0.0*dd_h;
	 return;
  }

  //Calculate the second derivative of the lapse
  double alpha;
  double prefac;
  tensor delta(2,3,3);
  delta.SetName("delta");
  for(i = 0; i < 3; i++) delta.Set(1.0,i,i);

  alpha = sqrt( 1.0 - 2.0*M/r );
  prefac = M / ( pow(alpha,3.0) * pow(r,6.0) );

  dd_lapse <= prefac*(  pow(r,3.0)*alpha*alpha*delta 
                        - pos*pos*(3.0*r*alpha*alpha + M));

  //Calculate the second derivative of the shift
     //do nothing since zero


  //Calculate the second derivative of the spatial metric
  //create similar sized tensors to dd_h for convenience
  tensor part1, part2;
  part1.Resize(3,3,3,3);
  part2 <= dd_h;

  //Calculate part1 according to 10/25/98 notes
  prefac = -10.0*M/pow(r,7.0);
  for(i = 0; i < 3; i++)
   for(j = 0; j < 3; j++)
    for(k = 0; k < 3; k++)
      part1.Set(
                  3*pos.Val(i)*pos.Val(j)*pos.Val(k)
                - r*r*pos.Val(j)*delta.Val(i,k)
                - r*r*pos.Val(i)*delta.Val(j,k),
                  i,j,k
               );

   part1 <= prefac*part1*pos;

  //Calculate part2 according to 10/25/98 notes
  prefac = 2.0*M/pow(r,5.0);
  for(i = 0; i < 3; i++)
   for(j = 0; j < 3; j++)
    for(k = 0; k < 3; k++)
     for(l = 0; l < 3; l++)
       part2.Set(
                   3.0*delta.Val(i,l)*pos.Val(j)*pos.Val(k)
                 + 3.0*delta.Val(j,l)*pos.Val(i)*pos.Val(k)
                 + 3.0*delta.Val(k,l)*pos.Val(i)*pos.Val(j)
                 - 2.0*pos.Val(l)*pos.Val(j)*delta.Val(i,k)
                 - 2.0*pos.Val(l)*pos.Val(i)*delta.Val(j,k)
                 - r*r*delta.Val(j,l)*delta.Val(i,k)
                 - r*r*delta.Val(i,l)*delta.Val(j,k),
                   i,j,k,l
                );

  part2 <=prefac*part2;

  //Combine part1 and part2 to form dd_h
  dd_h <= part1 + part2;
}


/******************************************************************************/
/*     Function to calculate the metric and its first derivatives for the     */
/*     standard Schwarzschild metric                                          */
/******************************************************************************/
void load_field(struct sim_parms *parms)
{

  tensor lapse(1,1),      shift(1,3),        h(2,3,3);
  tensor d_lapse(1,3),    d_shift(2,3,3),    d_h(3,3,3,3);
  tensor dd_lapse(2,3,3), dd_shift(3,3,3,3), dd_h(4,3,3,3,3);

  tensor temp;
  tensor pos;
  int i, j, k, m, indices[3];


  //load the tensor portion of the field
  for(i = 0; i < parms->lapse.GetNumIndices(X_ID); i++)
   {
    for(j = 0; j < parms->lapse.GetNumIndices(Y_ID); j++)
     {
      for(k = 0; k < parms->lapse.GetNumIndices(Z_ID); k++)
       {
         indices[0] = i;
         indices[1] = j;
         indices[2] = k;
         //Pack the position
         pos <= parms->lapse.GetPosition(indices);

         //Get the local value of the fields
         std_schw(parms, pos, lapse, shift, h, d_lapse, d_shift, d_h);

         //Set the local value of the fields
         parms->   lapse.SetTensor(indices, lapse);
         parms->   shift.SetTensor(indices, shift);
         parms->       h.SetTensor(indices, h);
         parms-> d_lapse.SetTensor(indices, d_lapse);
         parms-> d_shift.SetTensor(indices, d_shift);
         parms->     d_h.SetTensor(indices, d_h);
         if(max_order == 2)
          {
            schw_higher_order(parms, pos, dd_lapse, dd_shift, dd_h);
            parms->dd_lapse.SetTensor(indices, dd_lapse);
            parms->dd_shift.SetTensor(indices, dd_shift);
            parms->    dd_h.SetTensor(indices, dd_h);
          }
       }
     }
   }
  printf("Finished loading fields\n");
}

/******************************************************************************/
/*     Function to initialize the fields required for the grid method ---     */
/*     (not called for the analytic method)                                   */
/******************************************************************************/
void initialize(struct sim_parms *parms)
{
  FILE *input_file;
  char dummy[80];
  double spatial_extent[6];
  int num_grid_pts[3];
  int part_smth;
  int grid_smth;
  int basic_order;
  int d_order;
  int dep_field;
  int grid_smth_rad;

  input_file = fopen("schw.dat","r");

  //Resize the fields due to struct problem
  parms->   lapse.Resize(4,2,2,2,1);
  parms->   shift.Resize(4,2,2,2,3);
  parms->       h.Resize(5,2,2,2,3,3);
  parms-> d_lapse.Resize(4,2,2,2,3);
  parms-> d_shift.Resize(5,2,2,2,3,3);
  parms->     d_h.Resize(6,2,2,2,3,3,3);
  parms->dd_lapse.Resize(5,2,2,2,3,3);
  parms->dd_shift.Resize(6,2,2,2,3,3,3);
  parms->    dd_h.Resize(7,2,2,2,3,3,3,3);

  //Set the names
  parms->   lapse.SetName("lapse");
  parms->   shift.SetName("shift");
  parms->       h.SetName("h");
  parms-> d_lapse.SetName("d_lapse");
  parms-> d_shift.SetName("d_shift");
  parms->     d_h.SetName("d_h");
  parms->dd_lapse.SetName("dd_lapse");
  parms->dd_shift.SetName("dd_shift");
  parms->    dd_h.SetName("dd_h");

  //read the banner and blank line
  GetLine(input_file, dummy);
  GetLine(input_file, dummy);

  //read in the number of grid points in x direction
  GetLine(input_file, dummy);
  GetLine(input_file, dummy);
  num_grid_pts[0] = atoi(dummy);

  //read in the number of grid points in y direction
  GetLine(input_file, dummy);
  GetLine(input_file, dummy);
  num_grid_pts[1] = atoi(dummy);

  //read in the number of grid points in z direction
  GetLine(input_file, dummy);
  GetLine(input_file, dummy);
  num_grid_pts[2] = atoi(dummy);

  //read in the x spatial extent
  GetLine(input_file, dummy);
  GetLine(input_file, dummy);
  spatial_extent[0] =  atof(dummy);
  GetLine(input_file, dummy);
  spatial_extent[1] =  atof(dummy);

  //read in the y spatial extent
  GetLine(input_file, dummy);
  GetLine(input_file, dummy);
  spatial_extent[2] =  atof(dummy);
  GetLine(input_file, dummy);
  spatial_extent[3] =  atof(dummy);

  //read in the z spatial extent
  GetLine(input_file, dummy);
  GetLine(input_file, dummy);
  spatial_extent[4] =  atof(dummy);
  GetLine(input_file, dummy);
  spatial_extent[5] =  atof(dummy);

  //read in the particle smoothing method
  GetLine(input_file, dummy);
  GetLine(input_file, dummy);
  part_smth = atoi(dummy);
  parms->part_sm_meth = part_smth;
  parms->   lapse.SetSmoother(PARTICLE_SMOOTHING, part_smth);
  parms->   shift.SetSmoother(PARTICLE_SMOOTHING, part_smth);
  parms->       h.SetSmoother(PARTICLE_SMOOTHING, part_smth);
  parms-> d_lapse.SetSmoother(PARTICLE_SMOOTHING, part_smth);
  parms-> d_shift.SetSmoother(PARTICLE_SMOOTHING, part_smth);
  parms->     d_h.SetSmoother(PARTICLE_SMOOTHING, part_smth);
  parms->dd_lapse.SetSmoother(PARTICLE_SMOOTHING, part_smth);
  parms->dd_shift.SetSmoother(PARTICLE_SMOOTHING, part_smth);
  parms->    dd_h.SetSmoother(PARTICLE_SMOOTHING, part_smth);

  //read in the grid smoothing method
  GetLine(input_file, dummy);
  GetLine(input_file, dummy);
  grid_smth = atoi(dummy);
  parms->part_sm_meth = grid_smth;
  parms->   lapse.SetSmoother(GRID_SMOOTHING, grid_smth);
  parms->   shift.SetSmoother(GRID_SMOOTHING, grid_smth);
  parms->       h.SetSmoother(GRID_SMOOTHING, grid_smth);
  parms-> d_lapse.SetSmoother(GRID_SMOOTHING, grid_smth);
  parms-> d_shift.SetSmoother(GRID_SMOOTHING, grid_smth);
  parms->     d_h.SetSmoother(GRID_SMOOTHING, grid_smth);
  parms->dd_lapse.SetSmoother(GRID_SMOOTHING, grid_smth);
  parms->dd_shift.SetSmoother(GRID_SMOOTHING, grid_smth);
  parms->    dd_h.SetSmoother(GRID_SMOOTHING, grid_smth);

  //read in the grid smoothing radius
  GetLine(input_file, dummy);
  GetLine(input_file, dummy);
  grid_smth_rad = atoi(dummy);
  parms->   lapse.SetGridRadius(grid_smth_rad);
  parms->   shift.SetGridRadius(grid_smth_rad);
  parms->       h.SetGridRadius(grid_smth_rad);
  parms-> d_lapse.SetGridRadius(grid_smth_rad);
  parms-> d_shift.SetGridRadius(grid_smth_rad);
  parms->     d_h.SetGridRadius(grid_smth_rad);
  parms->dd_lapse.SetGridRadius(grid_smth_rad);
  parms->dd_shift.SetGridRadius(grid_smth_rad);
  parms->    dd_h.SetGridRadius(grid_smth_rad);

  //create some useful variables on the fly
  field *list1[3];
  int i;
  for( i = 0; i < 3; i++) list1[i] = NULL;

  //read in the dependency flag for the derivatives of the basic fields
  GetLine(input_file, dummy);
  GetLine(input_file, dummy);
  dep_field = atoi(dummy);

  //read in the order for the basic fields
  GetLine(input_file, dummy);
  GetLine(input_file, dummy);
  basic_order = atoi(dummy);
  if( basic_order > max_basic_order ) max_basic_order = basic_order;

  //read in the order for the first derivatives of the basic fields
  GetLine(input_file, dummy);
  GetLine(input_file, dummy);
  d_order = atoi(dummy);
  if( d_order > max_deriv_order ) max_deriv_order = d_order;

  if( max_basic_order > 1 || max_deriv_order > 0 ) max_order = 2;

  //Set the fields taylor's series expansion order
  parms->   lapse.SetSeriesOrder(basic_order);
  parms->   shift.SetSeriesOrder(basic_order);
  parms->       h.SetSeriesOrder(basic_order);
  parms-> d_lapse.SetSeriesOrder(d_order);
  parms-> d_shift.SetSeriesOrder(d_order);
  parms->     d_h.SetSeriesOrder(d_order);
  parms->dd_lapse.SetSeriesOrder(0);
  parms->dd_shift.SetSeriesOrder(0);
  parms->    dd_h.SetSeriesOrder(0);

  //Set the field as basic or as a derivative
  parms->   lapse.SetFieldDep(FALSE);
  parms->   shift.SetFieldDep(FALSE);
  parms->       h.SetFieldDep(FALSE);
  parms-> d_lapse.SetFieldDep(dep_field);
  parms-> d_shift.SetFieldDep(dep_field);
  parms->     d_h.SetFieldDep(dep_field);
  parms->dd_lapse.SetFieldDep(FALSE);
  parms->dd_shift.SetFieldDep(FALSE);
  parms->    dd_h.SetFieldDep(FALSE);

  //Set the field dependency
  parms->   lapse.SetDerivField(FALSE);
  parms->   shift.SetDerivField(FALSE);
  parms->       h.SetDerivField(FALSE);
  parms-> d_lapse.SetDerivField(TRUE);
  parms-> d_shift.SetDerivField(TRUE);
  parms->     d_h.SetDerivField(TRUE);
  parms->dd_lapse.SetDerivField(FALSE);
  parms->dd_shift.SetDerivField(FALSE);
  parms->    dd_h.SetDerivField(FALSE);

  list1[0] = &(parms->lapse);
  list1[1] = &(parms->d_lapse);
  list1[2] = &(parms->dd_lapse);
  parms->lapse.SetSeriesList(list1);
  parms->dd_lapse.SetSeriesList(list1);
  if( d_order < 2 && dep_field == FALSE )
  {
	  list1[0] = &(parms->d_lapse);
	  list1[1] = &(parms->dd_lapse);
  }
  if( d_order == 2 && dep_field == FALSE )
  {
	  cout << "Second order on derived fields - setting field as dependent" <<endl;
	  parms-> d_lapse.SetSeriesOrder(2);
     parms-> d_lapse.SetFieldDep(TRUE);
  }

  parms->d_lapse.SetSeriesList(list1);


  list1[0] = &(parms->shift);
  list1[1] = &(parms->d_shift);
  list1[2] = &(parms->dd_shift);
  parms->shift.SetSeriesList(list1);
  parms->dd_shift.SetSeriesList(list1);
  if( d_order < 2 && dep_field == FALSE )
  {
	  list1[0] = &(parms->d_shift);
	  list1[1] = &(parms->dd_shift);
  }
  if( d_order == 2 && dep_field == FALSE )
  {
	  cout << "Second order on derived fields - setting field as dependent" <<endl;
	  parms-> d_shift.SetSeriesOrder(2);
     parms-> d_shift.SetFieldDep(TRUE);
  }
  parms->d_shift.SetSeriesList(list1);

  list1[0] = &(parms->h);
  list1[1] = &(parms->d_h);
  list1[2] = &(parms->dd_h);
  parms->h.SetSeriesList(list1);
  parms->dd_h.SetSeriesList(list1);
  if( d_order <2 && dep_field == FALSE )
  {
	  list1[0] = &(parms->d_h);
	  list1[1] = &(parms->dd_h);
  }
  if( d_order == 2 && dep_field == FALSE )
  {
	  cout << "Second order on derived fields - setting field as dependent" <<endl;
	  parms-> d_h.SetSeriesOrder(2);
     parms-> d_h.SetFieldDep(TRUE);
  }
  parms->d_h.SetSeriesList(list1);

  //Resize the fields
  parms->   lapse.ResizeField(num_grid_pts);
  parms->   shift.ResizeField(num_grid_pts);
  parms->       h.ResizeField(num_grid_pts);
  parms-> d_lapse.ResizeField(num_grid_pts);
  parms-> d_shift.ResizeField(num_grid_pts);
  parms->     d_h.ResizeField(num_grid_pts);
  if( max_order > 1 )
   {
    parms->dd_shift.ResizeField(num_grid_pts);
    parms->dd_lapse.ResizeField(num_grid_pts);
    parms->    dd_h.ResizeField(num_grid_pts);
   }

  //Now fill in the spatial extents of the fields
  parms->   lapse.SetSpatialExtent(spatial_extent);
  parms->   shift.SetSpatialExtent(spatial_extent);
  parms->       h.SetSpatialExtent(spatial_extent);
  parms-> d_lapse.SetSpatialExtent(spatial_extent);
  parms-> d_shift.SetSpatialExtent(spatial_extent);
  parms->     d_h.SetSpatialExtent(spatial_extent);
  parms->dd_lapse.SetSpatialExtent(spatial_extent);
  parms->dd_shift.SetSpatialExtent(spatial_extent);
  parms->    dd_h.SetSpatialExtent(spatial_extent);

  //close the input file
  fclose(input_file);
}

/******************************************************************************/
/*     Function to inspect the smoothing method                               */
/******************************************************************************/
void inspect_smoothing(struct sim_parms *parms)
{
  tensor temp, zero_pos;
  char *dump;
  int i;

  dump = NULL;

  char out_file[80];
  FILE *out;
  FILE *fields;
  cout << "\nEnter the output filename" << endl;
  cin >> out_file;
  out = fopen(out_file,"w");

  cout << "Please enter the central body mass and the particle mass" << endl;
  cin >> parms->M >> parms->m;

  load_field(parms);

  int choice;
  cout << "Print the fields to a file? (0 = no; 1 = yes)" << endl;
  cin >> choice;

  if( choice == 1 )
   {
    fields = fopen("fields.prn","w");
    parms->lapse.print("lapse",&dump);
    fprintf(fields,"%s\n\n",dump);
//    parms->shift.print("shift",&dump);
//    fprintf(fields,"%s\n\n",dump);
//    parms->h.print("h",&dump);
//    fprintf(fields,"%s\n\n",dump);
//    parms->d_lapse.print("d_lapse",&dump);
//    fprintf(fields,"%s\n\n",dump);
//    parms->d_shift.print("d_shift",&dump);
//    fprintf(fields,"%s\n\n",dump);
//    parms->d_h.print("d_h",&dump);
//    fprintf(fields,"%s\n\n",dump);
   }

  double position[3];
  for(i = 0; i < 3; i++) position[i] = 0.0;

  while( position[0] != -9999)
   {
     cout << "\nInput the position vector (enter -9999 to end)" << endl;
     cin >> position[0];
     if(position[0] == -9999) break;
     cin >> position[1];
     cin >> position[2];

     temp <= zero_pos;
     for( i = 0; i < 3; i++) temp.Set(position[i],i);
     temp.print("position",&dump);
     fprintf(out,"%s",dump);

     temp <= parms->lapse.GetInterpTensor(position);
     temp.print("lapse",&dump);
     fprintf(out,"%s",dump);

     temp <= parms->h.GetInterpTensor(position);
     temp.print("h",&dump);
     fprintf(out,"%s",dump);

     temp <= parms->d_lapse.GetInterpTensor(position);
     temp.print("d_lapse",&dump);
     fprintf(out,"%s",dump);

     temp <= parms->d_h.GetInterpTensor(position);
     temp.print("d_h",&dump);
     fprintf(out,"%s",dump);

   }
   delete [] dump;
   exit(1);
}

/******************************************************************************/
/*     Utility function to read in a dummy string from an input file          */
/******************************************************************************/
void GetLine(FILE *file, char *string)
{
   int i = 0;
   char c;
	while( (c = fgetc(file)) != '\n' && i < 80 )
     {
       string[i] = c;
       i++;
     }
   string[i] = '\0';

}


