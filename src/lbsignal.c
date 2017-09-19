/*
 * This file is part of SunlightLB - a 3D Lattice Boltzmann code
 * Copyright (C) 2005 Unilever UK Central Resources Ltd.
 *  
 * SunlightLB is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2 
 * of the License, or (at your option) any later version. 
 *
 * SunlightLB is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details. 
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. 
 */

#define LBSIGNAL_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <unistd.h>
#include "sunlightlb.h"

/* Signalling routines can be switched out by compiling with -D_NO_SIGNALLING */

#ifdef SIGNAL

/* This routine is called at a convenient point to handle a
 * signal received by the program.  Its action depends on the
 * value of slb_signal_flag which indicates what kind of interrupt
 * was received.  This flag is reset to zero at the end.
 */

void sighandle(LBSys lb) {
  char *dmsfile;
  char *objfile;
  if (slb_signal_flag == 0) return;
  switch (slb_signal_flag) {
  case 1:
    printf("Acting on -USR1 signal\n");
    if (lb->name != NULL) wrmon(lb, "running");
    break;
  case 2:
    printf("Acting on -USR2 signal\n");
    if (lb->name != NULL) {
      dmsfile = mycat(lb->name, ".dms");
      objfile = mycat(lb->name, ".obj");
      lbsave(lb, dmsfile);
      if (lb->obj_head != NULL) lbsaveobj(lb, objfile);
      free(dmsfile); free(objfile);
    }
    break;
  }
  slb_signal_flag = 0;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : sighandle done\n", lb->step);
}


/* This routine is called when the process receives a signal.
 * It sets the value of slb_signal_flag (but only if it has been reset to zero
 * after a previous interrupt) so that the main loop can deal with it
 * at the most convenient time. Note that after a signal has been
 * received, it must be re-registered.
 */

void sighook(int sig) {
  int flag = 0;
  switch (sig) {
    case SIGUSR1: flag = 1; signal(SIGUSR1, sighook); break;
    case SIGUSR2: flag = 2; signal(SIGUSR2, sighook); break;
  }
  if (slb_signal_flag == 0) slb_signal_flag = flag;
}


/* This routine initialises the signal handling. */

void siginit() {
  slb_signal_flag = 0;
  signal(SIGUSR1, sighook);
  signal(SIGUSR2, sighook);
}

#endif /* SIGNAL */

/* End of lbsignal.c */
