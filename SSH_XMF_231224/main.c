/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief PLUTO main function.
 
  The file main.c contains the PLUTO main function and several other 
  top-level routines.
  main() provides basic code initialization, handles the the principal 
  integration loop and calls the output driver write_data.c.
  Other useful functions contained in this file are Integrate() which does
  the actual integration, GetNextTimeStep() responsible for computing the
  next time step based on the information available at the last time
  level.
  
  We use two slightly different integration loops depending on whether
  asnchrounous I/O has to be performed (macro USE_ASYNC_IO).
  If the macro USE_ASYNC_IO is not defined, the standard 
  integration loop consists of the following steps:
 
   - Check for last step & adjust dt if necessary
   - Dump log information, n, t(n), dt(n), MAX_MACH(n-1), etc..
   - Check output/analysis:  t(n) < tout < t(n)+dt(n)
   - write to disk/call analysis using {U(n), t(n), dt(n)}
   - Advance solution using dt(n): U(n) --> U(n+1)
   - Increment t(n+1) = t(n) + dt(n)
   - [MPI] Show dominant time step (n)
   - [MPI] Get next time step dt(n+1)
   - [MPI] reduction operations (n)
   - Increment n --> n+1
 
  Otherwise, using Asynchrounous I/O:
 
   - Check for last step & adjust dt
   - check for output/analysis:   t(n) < tout < t(n+1)
     - Write data/call analysis using {U(n), t(n), dt(n)}
   - [MPI] Show dominant time step (n-1)
   - [MPI] reduction operations (n-1)
   - [AIO]: finish writing
   - Dump log information, n, t(n), dt(n), MAX_MACH(n-1), etc..
   - Advance solution using dt(n), U(n) --> U(n+1)
   - Increment t(n+1) = t(n) + dt(n)
   - [MPI] Get next time step dt(n)
   - Increment n --> n+1

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "globals.h"

#include "mmd_sm_ports.h"

#ifndef SHOW_TIME_STEPS
  #define SHOW_TIME_STEPS  NO   /* -- show time steps due to advection,
                                       diffusion and cooling */
#endif

static char *TotalExecutionTime (double);
static void CheckForOutput (Data *, Runtime *, Grid *);
static void CheckForAnalysis (Data *, Runtime *, Grid *);

/* ********************************************************************* */
int main (int argc, char *argv[])
/*!
 * Start PLUTO, initialize functions, define data structures and 
 * handle the main integration loop.
 *
 * \param [in] argc Argument counts.
 * \param [in] argv Array of pointers to the strings.
 * \return This function return 0 on normal exit.
 *
 *********************************************************************** */
{
  int    n, nv, idim, err, f;
  char   first_step=1, last_step = 0;
  double scrh;

  Cmd_Line cmd_line;
  Grid     grd[3];

  Runtime  ini;  

  time_t tbeg, tend;  

  Model model;

  #ifdef PARALLEL
   AL_Init (&argc, &argv);
   MPI_Comm_rank (MPI_COMM_WORLD, &prank);
  #endif

//public Initialize
  mmd_initialize_public(argc, argv, &ini, grd, &cmd_line);

//static machine Initialize
  mmd_sm_initialize(&ini);

//model Initialize
  mmd_model_initialize(&model, &ini);

  mmd_model_IC(&model, grd, &ini);
  
  mmd_model_BC(&model, grd);
  
  time (&tbeg);
  g_stepNumber = 0;

//restart deleted

  print1 ("> Starting computation... \n\n");


/*
   FLD_LOOP(f)CheckForOutput  (&(fluid[f].d), ini+f, grd);
   CheckForAnalysis(&(fluid[f].d), ini, grd);
*/

/* =====================================================================
          M A I N      L O O P      S T A R T S      H E R E
   ===================================================================== */

  while (!last_step){//

  /* ------------------------------------------------------
      Check if this is the last integration step:
      - final tstop has been reached: adjust time step 
      - or max number of steps has been reached
     ------------------------------------------------------ */

    if ((g_time + g_dt) >= ini.tstop*(1.0 - 1.e-8)) {
      g_dt   = (ini.tstop - g_time);
      last_step = 1;
    }
    if (g_stepNumber == cmd_line.maxsteps && cmd_line.maxsteps > 0) {
      last_step = 1;
    }


  /* ------------------------------------------------------
                Dump log information
     ------------------------------------------------------ */
    mmd_model_output(&model, grd, &ini, g_inputParam[OUTPUT_INTERVAL], ini.output_dir);
 
    if (g_stepNumber%ini.log_freq == 0) {
      print1 ("step:%d ; t = %10.4e ; dt = %10.4e ; %d %% ; [%f, %d",
               g_stepNumber, g_time, g_dt, (int)(100.0*g_time/ini.tstop), 
               g_maxMach, g_maxRiemannIter);
      print1 ("]\n");
    }

    g_time += g_dt;
    scrh = g_dt;//save previous time step

    mmd_model_BC(&model, grd);

    mmd_AdvanceStep (&model, grd);

 
  /* ------------------------------------------------------
          Global MPI reduction operations
     ------------------------------------------------------ */
  
    #ifdef PARALLEL
     MPI_Allreduce (&g_maxMach, &scrh, 1, 
                    MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
     g_maxMach = scrh;

     MPI_Allreduce (&g_maxRiemannIter, &nv, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
     g_maxRiemannIter = nv;
    #endif

    g_stepNumber++;
    
    first_step = 0;
  }

/* =====================================================================
          M A I N       L O O P      E N D S       H E R E 
   ===================================================================== */
/*
  if (cmd_line.write){
    FLD_LOOP(f)CheckForOutput (&(fluid[f].d), ini+f, grd);
    CheckForAnalysis (&(fluid[f].d), ini, grd);
  }
*/
  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);
   print1  ("\n> Total allocated memory  %6.2f Mb (proc #%d)\n",
             (float)g_usedMemory/1.e6,prank);
   MPI_Barrier (MPI_COMM_WORLD);
  #else
   print1  ("\n> Total allocated memory  %6.2f Mb\n",(float)g_usedMemory/1.e6);
  #endif

  time(&tend);
  g_dt = difftime(tend, tbeg);
  print1("> Elapsed time             %s\n", TotalExecutionTime(g_dt));
  print1("> Average time/step        %10.2e  (sec)  \n", 
          difftime(tend,tbeg)/(double)g_stepNumber);
  print1("> Local time               %s",asctime(localtime(&tend)));
  print1("> Done\n");

 
  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);
   AL_Finalize ();
  #endif

  return (0);
}
#undef SHOW_TIME_STEPS

/* ********************************************************************* */
char *TotalExecutionTime (double dt)
/*!
 *
 *   convert a floating-point variable (dt, in seconds) to a string 
 *   displaying days:hours:minutes:seconds
 *   
 *********************************************************************** */
{
  static char c[128];
  int days, hours, mins, secs;

  days  = (int) (dt/86400.0);
  hours = (int) ((dt - 86400.0*days)/3600.0);
  mins  = (int) ((dt - 86400.0*days - 3600.0*hours)/60.);
  secs  = (int) (dt - 86400.0*days - 3600.0*hours - 60.0*mins);

  sprintf (c, " %dd:%dh:%dm:%ds", days,hours, mins, secs);
  return (c);
}


/* ********************************************************************* */
void CheckForOutput (Data *d, Runtime *ini, Grid *grid)
/*!
 *  Check if file output has to be performed.
 *  
 *********************************************************************** */
{
  static int first_call = 1;
  int  n, check_dt, check_dn, check_dclock;
  int  restart_update, last_step;
  double t, tnext;
  Output *output;
  static time_t clock_beg[MAX_OUTPUT_TYPES], clock_end;
  static double tbeg[MAX_OUTPUT_TYPES], tend;
  double dclock;

  restart_update = 0;
  t     = g_time;
  tnext = t + g_dt;
  
  last_step = (fabs(t-ini->tstop) < 1.e-12 ? 1:0);

/* -- on first execution initialize
      current beginning time for all output types -- */

  if (first_call){
    #ifdef PARALLEL
     if (prank == 0){
       double tstart;
       tstart = MPI_Wtime();
       for (n = 0; n < MAX_OUTPUT_TYPES; n++) tbeg[n] = tstart;
     }
     MPI_Bcast(tbeg, MAX_OUTPUT_TYPES, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #else
     for (n = 0; n < MAX_OUTPUT_TYPES; n++) time(clock_beg + n);
    #endif
  } 

/* -- get current time -- */

  #ifdef PARALLEL  
   if (prank == 0) tend = MPI_Wtime();
   MPI_Bcast(&tend, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #else
   time(&clock_end);
  #endif
  
/* -------------------------------------------------------
          start main loop on outputs
   ------------------------------------------------------- */
   
  for (n = 0; n < MAX_OUTPUT_TYPES; n++){
    output = ini->output + n;    
    check_dt = check_dn = check_dclock = 0; 

  /* -- check time interval in code units (dt) -- */

    if (output->dt > 0.0){
      check_dt = (int) (tnext/output->dt) - (int)(t/output->dt);
      check_dt = check_dt || g_stepNumber == 0 || last_step;
    }   

  /* -- check time interval in number of steps (dn) -- */

    if (output->dn > 0){
      check_dn = (g_stepNumber%output->dn) == 0;
      check_dn = check_dn || g_stepNumber == 0 || last_step;
    }

  /* -- check time interval in clock time (dclock) -- */

    if (output->dclock > 0.0){
      #ifdef PARALLEL
       dclock = tend - tbeg[n];
      #else
       dclock = difftime(clock_end, clock_beg[n]);
      #endif
      if (dclock >= output->dclock) {
        check_dclock = 1;
        #ifdef PARALLEL
         tbeg[n] = tend;
        #else
         time(clock_beg + n);
        #endif
      }else{ 
        check_dclock = 0;
      }
      check_dclock = check_dclock || g_stepNumber == 0 || last_step;
    }

  /* -- if any of the previous is true dump data to disk -- */

    if (check_dt || check_dn || check_dclock) {  
       /* --------------------------------------------------------
            Get user var if necessary 
        -------------------------------------------------------- */
         
       if(d->Vuser!=0) ComputeUserVar(d, grid);
  
       WriteData(d, output, grid);

    /* ----------------------------------------------------------
        save the file number of the dbl and dbl.h5 output format
        for writing restart.out once we exit the loop.
       ---------------------------------------------------------- */

      if ((output->type == DBL_OUTPUT) ||
          (output->type == DBL_H5_OUTPUT)) restart_update = 1;
    }
  }

  first_call = 0;
}

/* ******************************************************************** */
void CheckForAnalysis (Data *d, Runtime *ini, Grid *grid)
/*
 *
 * PURPOSE 
 *
 *   Check if Analysis needs to be called
 *
 ********************************************************************** */
{
  int check_dt, check_dn;
  double t, tnext;

  t     = g_time;
  tnext = t + g_dt;
  check_dt = (int) (tnext/ini->anl_dt) - (int)(t/ini->anl_dt);
  check_dt = check_dt || g_stepNumber == 0 || fabs(t - ini->tstop) < 1.e-9; 
  check_dt = check_dt && (ini->anl_dt > 0.0);

  check_dn = (g_stepNumber%ini->anl_dn) == 0;
  check_dn = check_dn && (ini->anl_dn > 0);

  if (check_dt || check_dn) Analysis (d, grid);

  return;
}


//trash can
/*
  print1 ("> Basic data type:\n");
  print1 ("  sizeof (char)     = %d\n", sizeof(char));
  print1 ("  sizeof (uchar)    = %d\n", sizeof(unsigned char));
  print1 ("  sizeof (short)    = %d\n", sizeof(short));
  print1 ("  sizeof (ushort)   = %d\n", sizeof(unsigned short));
  print1 ("  sizeof (int)      = %d\n", sizeof(int));
  print1 ("  sizeof (*int)     = %d\n", sizeof(*int));
  print1 ("  sizeof (float)    = %d\n", sizeof(float));
  print1 ("  sizeof (double)   = %d\n", sizeof(double));
  print1 ("  sizeof (*double)  = %d\n", sizeof(*double));

  print1 ("\n> Structure data type:\n");
  print1 ("  sizeof (CMD_LINE)   = %d\n", sizeof(Cmd_Line));
  print1 ("  sizeof (DATA)       = %d\n", sizeof(Data));
  print1 ("  sizeof (STATE_1D)   = %d\n", sizeof(State_1D));
  print1 ("  sizeof (GRID)       = %d\n", sizeof(Grid));
  print1 ("  sizeof (TIME_STEP)  = %d\n", sizeof(Time_Step));
  print1 ("  sizeof (OUTPUT)     = %d\n", sizeof(Output));
  print1 ("  sizeof (RUNTIME)    = %d\n", sizeof(Runtime));
  print1 ("  sizeof (RESTART)    = %d\n", sizeof(Restart));
  print1 ("  sizeof (RGB)        = %d\n", sizeof(RGB));
  print1 ("  sizeof (IMAGE)      = %d\n", sizeof(Image));
  print1 ("  sizeof (FLOAT_VECT) = %d\n", sizeof(Float_Vect));
  print1 ("  sizeof (INDEX)      = %d\n", sizeof(Index));
  print1 ("  sizeof (RBOX)       = %d\n", sizeof(RBox));
*/



