//**************************************************************
//crank nicholson algorithm for ODEs
//
//
//**************************************************************
#include "crnk_nch.h"


void crank_nicholson(double current_time, double time_step,
                     int    state_size,   double *state,
                     double *new_state,   struct sim_parms *parms,
                     tensor (*rhs)(double time, const tensor &state,
                                   struct sim_parms *parms))
{
    tensor prev_state(1,state_size);
    tensor last_guess(1,state_size);
    tensor current_guess(1,state_size);
    tensor temp(1,state_size);

    double delta;
    double time_step_2 = time_step/2.0;
    int i, counter;

    for( i = 0; i < state_size; i++ ) prev_state.Set(state[i],i);

    last_guess <= prev_state;
    current_guess <= prev_state;
    delta = 1e200;
    counter = 0;

    while( delta > CN_TOL && counter <= CN_MAX_ITERATES)
     {
       temp <= time_step_2*( rhs(current_time, prev_state, parms) +
                             rhs(current_time + time_step, last_guess, parms));
       current_guess <= prev_state + temp;
       temp <= current_guess - last_guess;
       delta = temp.max();
       last_guess <= current_guess;
       counter++;
     }

    for( i = 0; i < state_size; i++ )
      new_state[i] =  current_guess.Val(i);

}