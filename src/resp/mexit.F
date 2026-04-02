C
C************************************************************************
C                              AMBER                                   **
C                                                                      **
C                  Copyright (c) 1986, 1991, 1995                      **
C             Regents of the University of California                  **
C                       All Rights Reserved.                           ** 
C                                                                      **
C  This software provided pursuant to a license agreement containing   **
C  restrictions on its disclosure, duplication, and use. This software **
C  contains confidential and proprietary information, and may not be   **
C  extracted or distributed, in whole or in part, for any purpose      **
C  whatsoever, without the express written permission of the authors.  **
C  This notice, and the associated author list, must be attached to    **
C  all copies, or extracts, of this software. Any additional           **
C  restrictions set forth in the license agreement also apply to this  **
C  software.                                                           **
C************************************************************************
C
C  mexit() - machine-dependent exit() procedure, designed to return an 
C            appropriate (success/failure) value to the operating system.
C
      subroutine mexit(ifil, i)
#ifdef MPI
# include "parallel.h"
# include "mpif.h"
#endif
#ifdef IBM_VM_CMS
        if (ifil.ne.0 .and. ifil.ne.6) close(unit=ifil)
#else
        if (ifil.ne.0) close(unit=ifil)
#endif
#ifdef MPI
c
c       ...i .gt. 0 implies and error condition, therefore we want to
c       kill all the nodes.  This is accomplished with mpi_abort.  If
c       it is not an error, exit gracefully with the mpi_finalize.
c       NOTE: no mpi functions may be called after a call to 
c       mpi_finalize.
c
        if (i.ne.0) then
          call MPI_ABORT(MPI_COMM_WORLD, i, ier)
        else
          call MPI_FINALIZE(ier)
        endif
#endif
#ifdef VMS
c       if anyone knows how to get a reasonable %DCL msg w/ exit
c       status, please send it in. this will at least cause DCL
c       "on error then goto BOMB" to trap.
c
        if (i.ne.0) then
c         exit with an _even_ value to signal error
          call exit(2)
        else
          call exit(1)
        endif
        
#else 
#ifdef F2C
        if (i.ne.0) then
          call exit(1)
        else
          stop 0
        endif
#else
#if rs6000 || IBM3090
        if (i.ne.0) then
          stop 1
        else
          stop 0
        endif
#else
#ifdef HP
        if (i.ne.0) then
          call exit(1)
        else
          stop 'ok'
        endif
#else
#ifdef UXP
        call setrcd(i)
#else
        call exit(i)
#endif
#endif
#endif
#endif
#endif
      end



