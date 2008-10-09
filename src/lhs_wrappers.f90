!   _______________________________________________________________________
!
!   DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
!   Copyright (c) 2006, Sandia National Laboratories.
!   This software is distributed under the GNU General Public License.
!   For more information, see the README file in the top Dakota directory.
!   _______________________________________________________________________
!
C These Fortran wrappers circumvent problems with implicit string sizes
C in f90.

C -----------------------------
C Wrapper for LHS's lhs_options
C -----------------------------
      subroutine lhs_options2( lhsreps, lhspval, lhsopts, ierror )

C Fix the string size and always call lhs_options2 from C++ with strings of
C length 32
      character*32 lhsopts
      integer      lhsreps, lhspval, ierror

C Since calling from F90 now, the implicit string size passing should work
      call lhs_options( lhsreps, lhspval, lhsopts, ierror )

      end

C --------------------------
C Wrapper for LHS's lhs_dist
C --------------------------
      subroutine lhs_dist2( namvar, iptflag, ptval, distype, aprams,
     1                      numprms, ierror, idistno, ipvno )

C Fix the string size and always call lhs_dist2 from C++ with strings of
C length 32
      character*16     namvar
      character*32     distype
      integer          iptflag, numprms, ierror, idistno, ipvno
      double precision ptval, aprams(1) 

C Since calling from F90 now, the implicit string size passing should work
      call lhs_dist( namvar, iptflag, ptval, distype, aprams,
     1               numprms, ierror, idistno, ipvno )

      end

C ---------------------------
C Wrapper for LHS's lhs_udist
C ---------------------------
      subroutine lhs_udist2( namvar, iptflag, ptval, distype, numpts,
     1                       xval, yval, ierror, idistno, ipvno )

C Fix the string size and always call lhs_udist2 from C++ with strings of
C length 32
      character*16     namvar
      character*32     distype
      integer          iptflag, numpts, ierror, idistno, ipvno
      double precision ptval, xval(1), yval(1)

C Since calling from F90 now, the implicit string size passing should work
      call lhs_udist( namvar, iptflag, ptval, distype, numpts,
     1                xval, yval, ierror, idistno, ipvno )

      end

C --------------------------
C Wrapper for LHS's lhs_corr
C --------------------------
      subroutine lhs_corr2( nam1, nam2, corrval, ierror )

C Fix the string size and always call lhs_corr2 from C++ with strings of
C length 32
      character*16     nam1, nam2
      integer          ierror
      double precision corrval

C Since calling from F90 now, the implicit string size passing should work
      call lhs_corr( nam1, nam2, corrval, ierror )

      end

C --------------------------
C Wrapper for LHS's lhs_run
C --------------------------
      subroutine lhs_run2( max_var, max_obs, max_names, ierror, 
     1                     dist_names, name_order, pt_vals, num_names,
     2                     sample_matrix, num_vars, rank_matrix, rflag )

      integer          max_var, max_obs, max_names, num_names, num_vars
      integer          rflag, ierror, name_order(1) 
      character*16     dist_names(1) 
      double precision pt_vals(1), sample_matrix(1),rank_matrix(1)

      call lhs_run( max_var, max_obs, max_names, ierror, 
     1              dist_names, name_order, pt_vals, num_names,
     2              sample_matrix, num_vars, rank_matrix, rflag )

      end

C ---------------------------
C Wrapper for LHS's lhs_files
C ---------------------------
      subroutine lhs_files2( lhsout, lhsmsg, lhstitl, lhsopts, ierror )

C Fix the string size and always call lhs_files from C++ with strings of
C length 32
      character*32 lhsout, lhsmsg, lhstitl, lhsopts
      integer      ierror

C Since calling from F90 now, the implicit string size passing should work
      call lhs_files( lhsout, lhsmsg, lhstitl, lhsopts, ierror )
      
      end
