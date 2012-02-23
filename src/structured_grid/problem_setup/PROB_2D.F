#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "BC_TYPES_PM.H"
#include "PROB_AMR_F.H"
#include "PROB_PM_F.H"
#include "ArrayLim.H"

#define SDIM 2

c ::: -----------------------------------------------------------
c ::: This routine is called at problem initialization time
c ::: and when restarting from a checkpoint file.
c ::: The purpose is (1) to specify the initial time value
c ::: (not all problems start at time=0.0) and (2) to read
c ::: problem specific data from a namelist or other input
c ::: files and possibly store them or derived information
c ::: in FORTRAN common blocks for later use.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: init      => TRUE if called at start of problem run
c :::              FALSE if called from restart
c ::: name      => name of "probin" file
c ::: namlen    => length of name
c ::: strttime <=  start problem with this time variable
c ::: 
c ::: -----------------------------------------------------------
      subroutine FORT_PROBINIT (init,name,namlen,problo,probhi)

      implicit none

      integer init,namlen
      integer name(namlen)
      integer untin, i
      REAL_T  problo(SDIM), probhi(SDIM)

#include "probdata.H"

c
c Dimensions of the Inflow file.
c

      integer maxlen, isioproc
      parameter (maxlen=256)

      character probin*(maxlen)

      call bl_pd_is_ioproc(isioproc)

      if (namlen .gt. maxlen) call bl_error('probin file name too long')

      do i = 1, namlen
         probin(i:i) = char(name(i))
      end do

#if 0
c     REMOVE NEED FOR PROBIN...WILL SET THESE IN c++ EVENTUALLY ANYWAY
      untin = 9
      if (namlen .eq. 0) then
         open(untin,file='probin',form='formatted',status='old')
      else
         open(untin,file=probin(1:namlen),form='formatted',status='old')
      end if

      read(untin,fortin)
      if (isioproc .eq. 1) write(6,fortin)
      close(unit=untin)
#endif

      probtype = 1
      denerr = 1.d0
      adverr = 10.d0

      end

c ::: -----------------------------------------------------------
c ::: This user-defined routine is called at setup and 
c ::: initializes components.
c ::: 
c ::: NOTE:  all arrays have one cell of ghost zones surrounding
c :::        the grid interior.  Values in these cells need not
c :::        be set here.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: level     => amr level of grid
c ::: time      => time at which to init data          
c ::: nscal     => number of scalar quantities. 
c ::: scal     <=  Scalar array
c ::: dx        => cell size
c ::: -----------------------------------------------------------
      subroutine FORT_INITDATA(level,time,
     &     scal,DIMS(state),rhoinit,nscal,dx)
      
      integer level, nscal
      integer DIMDEC(state)
      REAL_T  time, dx(SDIM)
      REAL_T  scal(DIMV(state),1:nscal)
      REAL_T  rhoinit(1:nscal)
      
      end

c ::: -----------------------------------------------------------
c ::: This routine is called at setup and initializes tracers
c ::: 
c ::: NOTE:  all arrays have one cell of ghost zones surrounding
c :::        the grid interior.  Values in these cells need not
c :::        be set here.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: level     => amr level of grid
c ::: time      => time at which to init data  
c ::: nscal     => number of scalar quantities. 
c ::: ntracer   => number of tracers
c ::: scal     <=  Scalar array
c ::: dx        => cell size
c ::: -----------------------------------------------------------
      subroutine FORT_INIT_TRACER(level,time,
     &     scal,DIMS(state),tinit,nscal,ntracer,dx)
      
      integer level, nscal, ntracer
      integer DIMDEC(state)
      REAL_T  time, dx(SDIM)
      REAL_T  scal(DIMV(state),1:nscal+ntracer)
      REAL_T  tinit(1:ntracer)
      
      end

c ::: -----------------------------------------------------------
c ::: This routine will tag high error cells based on the 
c ::: density gradient
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: tag      <=  integer tag array
c ::: DIMS(tag) => index extent of tag array
c ::: set       => integer value to tag cell for refinement
c ::: clear     => integer value to untag cell
c ::: adv       => scalar array
c ::: DIMS(adv) => index extent of scalar array
c ::: lo,hi     => index extent of grid
c ::: nvar      => number of components in rho array (should be 1)
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: -----------------------------------------------------------
      subroutine FORT_ADVERROR (tag,DIMS(tag),set,clear,
     &                          adv,DIMS(adv),lo,hi,nvar,
     &                          domlo,domhi,dx,xlo,
     &			        problo,time,level)
#include "probdata.H"
      
      integer   DIMDEC(tag)
      integer   DIMDEC(adv)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      integer   lo(SDIM), hi(SDIM)
      integer   tag(DIMV(tag))
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      REAL_T    adv(DIMV(adv),1)

      integer   i, j

      if (model .ge. 0) then
         do i = lo(1), hi(1)
            do j = lo(2), hi(2)
              tag(i,j) = merge(set,tag(i,j),adv(i,j,1) .gt. adverr)
           end do
        end do
      else
        call bl_abort('model not defined IN FORT_ADVERROR ')
      end if

c     the tagging of inflow and outflow should be handled through 
c     manual_tags_placement
      
      end

c ::: -----------------------------------------------------------
c ::: This routine will tag high error cells based on the 
c ::: value provided
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: tag      <=  integer tag array
c ::: DIMS(tag) => index extent of tag array
c ::: set       => integer value to tag cell for refinement
c ::: clear     => integer value to untag cell
c ::: adv       => scalar array
c ::: DIMS(adv) => index extent of scalar array
c ::: lo,hi     => index extent of grid
c ::: nvar      => number of components in rho array (should be 1)
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: value     => threshold value
c ::: min_time  => minimum time to apply this condition
c ::: max_time  => minimum time to apply this condition
c ::: max_level => maximum level to be generated by this condition
c ::: -----------------------------------------------------------
      subroutine FORT_VALGTERROR (tag,DIMS(tag),set,clear,
     &                            adv,DIMS(adv),mask,DIMS(mask),
     &                            lo,hi,nvar,
     &                            domlo,domhi,dx,xlo,
     &			          problo,time,level,value)

      implicit none

#include "probdata.H"
      
      integer   DIMDEC(tag)
      integer   DIMDEC(adv)
      integer   DIMDEC(mask)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      integer   lo(SDIM), hi(SDIM)
      integer   tag(DIMV(tag))
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      REAL_T    adv(DIMV(adv),1)
      REAL_T    mask(DIMV(mask))
      REAL_T    value

      integer   i, j

      do i = lo(1), hi(1)
         do j = lo(2), hi(2)
            tag(i,j) = merge(set,tag(i,j), mask(i,j).eq.1 .and. adv(i,j,1).gt.value)
         end do
      end do

      end

c ::: -----------------------------------------------------------
c ::: This routine will tag high error cells based on the 
c ::: value provided
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: tag      <=  integer tag array
c ::: DIMS(tag) => index extent of tag array
c ::: set       => integer value to tag cell for refinement
c ::: clear     => integer value to untag cell
c ::: adv       => scalar array
c ::: DIMS(adv) => index extent of scalar array
c ::: lo,hi     => index extent of grid
c ::: nvar      => number of components in rho array (should be 1)
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: value     => threshold value
c ::: min_time  => minimum time to apply this condition
c ::: max_time  => minimum time to apply this condition
c ::: max_level => maximum level to be generated by this condition
c ::: -----------------------------------------------------------
      subroutine FORT_VALLTERROR (tag,DIMS(tag),set,clear,
     &                            adv,DIMS(adv),mask,DIMS(mask),
     &                            lo,hi,nvar,
     &                            domlo,domhi,dx,xlo,
     &			          problo,time,level,value)

      implicit none

#include "probdata.H"
      
      integer   DIMDEC(tag)
      integer   DIMDEC(adv)
      integer   DIMDEC(mask)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      integer   lo(SDIM), hi(SDIM)
      integer   tag(DIMV(tag))
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      REAL_T    adv(DIMV(adv),1)
      REAL_T    mask(DIMV(mask))
      REAL_T    value

      integer   i, j

      do i = lo(1), hi(1)
         do j = lo(2), hi(2)
            tag(i,j) = merge(set,tag(i,j), mask(i,j).eq.1 .and. adv(i,j,1).lt.value)
         end do
      end do

      end

c ::: -----------------------------------------------------------
c ::: This routine will tag high error cells based on the 
c ::: value provided
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: tag      <=  integer tag array
c ::: DIMS(tag) => index extent of tag array
c ::: set       => integer value to tag cell for refinement
c ::: clear     => integer value to untag cell
c ::: adv       => scalar array
c ::: DIMS(adv) => index extent of scalar array
c ::: lo,hi     => index extent of grid
c ::: nvar      => number of components in rho array (should be 1)
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: value     => threshold value
c ::: min_time  => minimum time to apply this condition
c ::: max_time  => minimum time to apply this condition
c ::: max_level => maximum level to be generated by this condition
c ::: -----------------------------------------------------------
      subroutine FORT_DIFFGTERROR (tag,DIMS(tag),set,clear,
     &                             adv,DIMS(adv),mask,DIMS(mask),
     &                             lo,hi,nvar,
     &                             domlo,domhi,dx,xlo,
     &			           problo,time,level,value)

      implicit none

#include "probdata.H"
      
      integer   DIMDEC(tag)
      integer   DIMDEC(adv)
      integer   DIMDEC(mask)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      integer   lo(SDIM), hi(SDIM)
      integer   tag(DIMV(tag))
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      REAL_T    adv(DIMV(adv),1)
      REAL_T    mask(DIMV(mask))
      REAL_T    value
      REAL_T    ax, ay, aerr

      integer   i, j

      do i = lo(1), hi(1)
         do j = lo(2), hi(2)
            ax = ABS(adv(i+1,j,1) - adv(i,j,1))
            ay = ABS(adv(i,j+1,1) - adv(i,j,1))
            ax = MAX(ax,ABS(adv(i,j,1) - adv(i-1,j,1)))
            ay = MAX(ay,ABS(adv(i,j,1) - adv(i,j-1,1)))
            aerr = MAX(ax,ay)
            tag(i,j) = merge(set,tag(i,j), mask(i,j).eq.1 .and. aerr.gt.value)
         end do
      end do
      
      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: E
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: rho      <=  density array
c ::: DIMS(rho) => index extent of rho array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_ONE_N_FILL (rho,DIMS(rho),domlo,domhi,dx,xlo,time,bc)

#include "probdata.H"

      integer    DIMDEC(rho)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     rho(DIMV(rho))
      integer    bc(SDIM,2)

      call filcc(rho(ARG_L1(rho),ARG_L2(rho)),
     &            DIMS(rho),domlo,domhi,dx,xlo,bc)   
  
      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: E
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: rho      <=  density array
c ::: DIMS(rho) => index extent of rho array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_ALL_N_FILL (rho,DIMS(rho),domlo,domhi,dx,xlo,time,bc)

      implicit none

#include "probdata.H"

      integer  DIMDEC(rho)
      integer  domlo(SDIM), domhi(SDIM)
      integer  bc(SDIM,2)
      REAL_T   dx(SDIM), xlo(SDIM), time
      REAL_T   rho(DIMV(rho),ncomps)

      integer  i,j,n
      integer  bctmp(SDIM,2)

      !
      ! I'll call filcc() with a slightly different bc.
      ! filcc() doesn't know about SEEPAGE.
      !
      do j=1,2
         do i=1,SDIM
            bctmp(i,j) = bc(i,j)
            if (bc(i,j) .eq. SEEPAGE .or. 
     &           bc(i,j) .eq. REFLECT_ODD) then
               bctmp(i,j) = FOEXTRAP
            end if
         end do
      end do
      

      do n = 1,ncomps
        call filcc(rho(ARG_L1(rho),ARG_L2(rho),n),
     &             DIMS(rho),domlo,domhi,dx,xlo,bctmp)
      end do

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation of Tracers 
c ::: when the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: rho      <=  density array
c ::: DIMS(rho) => index extent of rho array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_ALL_T_FILL (rho,DIMS(rho),domlo,domhi,dx,xlo,time,bc)

#include "probdata.H"

      integer  DIMDEC(rho)
      integer  domlo(SDIM), domhi(SDIM)
      integer  bc(SDIM,2)
      REAL_T   dx(SDIM), xlo(SDIM), time
      REAL_T   rho(DIMV(rho),ntracers)
      
      integer  i,j,n
      integer  bctmp(SDIM,2)

      !
      ! I'll call filcc() with a slightly different bc.
      ! filcc() doesn't know about SEEPAGE.
      !
      do j=1,2
         do i=1,SDIM
            bctmp(i,j) = bc(i,j)
            if (bc(i,j) .eq. SEEPAGE) bctmp(i,j) = FOEXTRAP
         end do
      end do
      
      do n = 1,ntracers
        call filcc(rho(ARG_L1(rho),ARG_L2(rho),n),
     &             DIMS(rho),domlo,domhi,dx,xlo,bctmp)
      end do

      end

c ::: -----------------------------------------------------------
c ::: Boundary function for COREREACT
c ::: -----------------------------------------------------------
      subroutine FORT_REC_N_FILL (rho,DIMS(rho),domlo,domhi,dx,xlo,time,bc)

#include "probdata.H"

      integer    DIMDEC(rho)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     rho(DIMV(rho),nreact)
      integer    bc(SDIM,2)
      integer    n

      do n = 1,nreact
        call filcc(rho(ARG_L1(rho),ARG_L2(rho),n),
     &             DIMS(rho),domlo,domhi,dx,xlo,bc)
      end do
     
      end

      subroutine FORT_XVELFILL (vel,DIMS(vel),domlo,domhi,dx,
     &                          xlo,time,bc )

      integer    DIMDEC(vel)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     vel(DIMV(vel))
      integer    bc(SDIM,2)

      integer    i, j

#include "probdata.H"

      call filcc(vel,DIMS(vel),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(vel).lt.domlo(1)) then
         do i = ARG_L1(vel), domlo(1)-1
            do j = ARG_L2(vel), ARG_H2(vel)
               vel(i,j) = zero
	    end do
	 end do
      end if            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(vel).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(vel)
            do j = ARG_L2(vel), ARG_H2(vel)
	       vel(i,j) = zero
	    end do
	 end do
      end if            


      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(vel).lt.domlo(2)) then
           do j = ARG_L2(vel), domlo(2)-1
              do i = ARG_L1(vel), ARG_H1(vel)
	        vel(i,j) = zero
	      end do
	   end do
      end if            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(vel).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(vel)
            do i = ARG_L1(vel), ARG_H1(vel)
	       vel(i,j) = zero
	    end do
	 end do
      end if            

      end

      subroutine FORT_YVELFILL (vel,DIMS(vel),domlo,domhi,dx,
     &                          xlo,time,bc )

      integer    DIMDEC(vel)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     vel(DIMV(vel))
      integer    bc(SDIM,2)

      integer    i, j

#include "probdata.H"

      call filcc(vel,DIMS(vel),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(vel).lt.domlo(1)) then
         do i = ARG_L1(vel), domlo(1)-1
            do j = ARG_L2(vel), ARG_H2(vel)
               vel(i,j) = zero
	    end do
	 end do
      end if            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(vel).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(vel)
            do j = ARG_L2(vel), ARG_H2(vel)
	       vel(i,j) = zero
	    end do
	 end do
      end if            


      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(vel).lt.domlo(2)) then
           do j = ARG_L2(vel), domlo(2)-1
              do i = ARG_L1(vel), ARG_H1(vel)
	        vel(i,j) = zero
	      end do
	   end do
      end if            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(vel).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(vel)
            do i = ARG_L1(vel), ARG_H1(vel)
	       vel(i,j) = zero
	    end do
	 end do
      end if            

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: adv      <=  advected quantity array
c ::: DIMS(adv) => index extent of adv array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of adv array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_ADVFILL (adv,DIMS(adv),domlo,domhi,dx,xlo,time,bc)

      integer    DIMDEC(adv)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     adv(DIMV(adv))
      integer    bc(SDIM,2)

      integer    i, j

#include "probdata.H"

      call filcc(adv,DIMS(adv),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(adv).lt.domlo(1)) then
         do i = ARG_L1(adv), domlo(1)-1
            do j = ARG_L2(adv), ARG_H2(adv)
	       adv(i,j) = zero
	    end do
	 end do
      end if            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(adv).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(adv)
            do j = ARG_L2(adv), ARG_H2(adv)
	       adv(i,j) = zero
	    end do
	 end do
      end if            

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(adv).lt.domlo(2)) then

         do j = ARG_L2(adv), domlo(2)-1
            do i = ARG_L1(adv), ARG_H1(adv)
	       adv(i,j) = zero
	    end do
	 end do

      end if            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(adv).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(adv)
            do i = ARG_L1(adv), ARG_H1(adv)
	       adv(i,j) = zero
	    end do
	 end do
      end if            

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.
c :::
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: temperature <=  temperature array
c ::: DIMS(temp)   => index extent of adv array
c ::: domlo,hi     => index extent of problem domain
c ::: dx           => cell spacing
c ::: xlo          => physical location of lower left hand
c :::                 corner of temperature array
c ::: time         => problem evolution time
c ::: bc           => array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_ENTHFILL (enthalpy,DIMS(temp),domlo,domhi,dx,
     &                          xlo,time,bc )

      integer    DIMDEC(temp)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     enthalpy(DIMV(temp))
      integer    bc(SDIM,2)

      integer    i, j

#include "probdata.H"

      call filcc(enthalpy,DIMS(temp),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(temp).lt.domlo(1)) then
         do i = ARG_L1(temp), domlo(1)-1
           do j = ARG_L2(temp), ARG_H2(temp)
               enthalpy(i,j) = one
           end do
         end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(temp).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(temp)
           do j = ARG_L2(temp), ARG_H2(temp)
               enthalpy(i,j) = one
           end do
         end do
      end if    

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(temp).lt.domlo(2)) then
         do j = ARG_L2(temp), domlo(2)-1
           do i = ARG_L1(temp), ARG_H1(temp)
               enthalpy(i,j) = one
          end do
       end do
      end if    

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(temp).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(temp)
           do i = ARG_L1(temp), ARG_H1(temp)
               enthalpy(i,j) = one
           end do
         end do
      end if    

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: p        <=  pressure array
c ::: DIMS(p)   => index extent of p array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi) 
c ::: -----------------------------------------------------------

      subroutine FORT_PRESFILL (p,DIMS(p),domlo,domhi,dx,xlo,time,bc)

      integer    DIMDEC(p)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     p(DIMV(p))
      integer    bc(SDIM,2)

      integer    i, j

      call filcc(p,DIMS(p),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(p).lt.domlo(1)) then
         do i = ARG_L1(p), domlo(1)-1
            do j = ARG_L2(p), ARG_H2(p)
               p(i,j) = p(i+1,j)
	    end do
	 end do
      end if            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(p).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(p)
            do j = ARG_L2(p), ARG_H2(p)
	       p(i,j) = p(i-1,j)
	    end do
	 end do
      end if            


      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(p).lt.domlo(2)) then
           do j = ARG_L2(p), domlo(2)-1
              do i = ARG_L1(p), ARG_H1(p)
	        p(i,j) = p(i,j+1)
	      end do
	   end do
      end if            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(p).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(p)
            do i = ARG_L1(p), ARG_H1(p)
	       p(i,j) = p(i,j-1)
	    end do
	 end do
      end if            

      end

c     NOTE: we hope that the ghost cells of pressure are correctly filled
c           by the solver, since we don't have good information here to fill them

c ::: -----------------------------------------------------------
c ::: This routine initializes the velocity at the beginning of 
c ::: each pressure solve.  It also construct RhoG for subsequent 
c ::: computation of \nabla RhoG.
c ::: -----------------------------------------------------------
      subroutine FORT_INIT_UMAC(umac,DIMS(umac),vmac,DIMS(vmac),
     &                          pc,DIMS(pc), lbd, DIMS(lbd),
     &                          ldx,DIMS(ldx),ldy,DIMS(ldy),
     &                          kpx,DIMS(kpx),kpy,DIMS(kpy),
     &                          RhoG,DIMS(rg),
     &                          RhoDX,DIMS(rdx),RhoDY,DIMS(rdy),
     &                          n,DIMS(n),
     &                          phi,DIMS(phi),
     &                          const_diff_coef,sbc,pbc,
     &                          dom_lo,dom_hi,dx,lo,hi,
     &                          wt_lo, wt_hi,
     &                          inflow_bc_lo,inflow_bc_hi,
     &                          inflow_vel_lo,inflow_vel_hi)

      implicit none

#include "probdata.H"

      integer DIMDEC(umac), DIMDEC(vmac)
      integer DIMDEC(ldx), DIMDEC(ldy)
      integer DIMDEC(kpx), DIMDEC(kpy)
      integer DIMDEC(pc)
      integer DIMDEC(lbd)
      integer DIMDEC(rg)
      integer DIMDEC(rdx)
      integer DIMDEC(rdy)
      integer DIMDEC(n)
      integer DIMDEC(phi)
      integer lo(SDIM),hi(SDIM)
      integer dom_lo(SDIM),dom_hi(SDIM)
      integer sbc(SDIM,2)
      integer pbc(SDIM,2)
      integer inflow_bc_lo(SDIM),inflow_bc_hi(SDIM)
      REAL_T dx(SDIM)
      REAL_T const_diff_coef(ncomps)
      REAL_T umac(DIMV(umac))
      REAL_T vmac(DIMV(vmac))
      REAL_T ldx(DIMV(ldx)),ldy(DIMV(ldy))
      REAL_T kpx(DIMV(kpx)),kpy(DIMV(kpy))  
      REAL_T pc(DIMV(pc))
      REAL_T lbd(DIMV(lbd),ncomps)
      REAL_T RhoG(DIMV(rg))
      REAL_T RhoDX(DIMV(rdx))
      REAL_T RhoDY(DIMV(rdy))
      REAL_T n(DIMV(n),ncomps)
      REAL_T phi(DIMV(phi))
      REAL_T wt_lo, wt_hi
      REAL_T inflow_vel_lo(SDIM),inflow_vel_hi(SDIM)

      integer i,j,k,l,nc
      integer compcount(nphases)
      integer compidx(nphases,ncomps)
      REAL_T  mu(nphases)
      REAL_T  y
      REAL_T  mucomp(ncomps),rhocomp(ncomps)
      REAL_T  dpc, stmp(ncomps),stmq(ncomps)
      REAL_T  ld(ncomps),kedge,ledge,tmp(2)
      REAL_T  rho, rho1, rho2
      REAL_T  s1,s2
      logical hv_diffuse


c     initialize compcount and compidx (F90 syntax)
      call get_compidx(compcount,compidx,phaseidx,nphases,ncomps)
      do i = 1,ncomps
         mucomp(i)  = muval(i)
         rhocomp(i) = rhoval(i)
      end do

c ::: If RHS = -\nabla \cdot U, the following constructs U
      do i = lo(1),hi(1)+1
         do j = lo(2),hi(2)
            if (model .eq. 2) then
               dpc = (pc(i,j) - pc(i-1,j))/dx(1)
               if (dpc .ne. zero) then
                  tmp(1) = lbd(i,j,1)
                  tmp(2) = lbd(i-1,j,1)
                  ld(1)  = (tmp(1)+tmp(2))*half

                  if (i .eq. dom_lo(1) .and. sbc(1,1) .ne. INT_DIR) then
                     dpc = dpc*two
                     ld  = lbd(i-1,j,:)
                  else if (i .eq. dom_hi(1)+1 .and. sbc(1,2) .ne. INT_DIR) then
                     dpc = dpc*2d0
                     ld = lbd(i,j,:)
                  end if
                  kedge = kpx(i,j)*ld(1) 
                  umac(i,j) = umac(i,j) + kedge*dpc
               end if
            end if
         end do
      end do

      do j = lo(2),hi(2)+1
         do i = lo(1),hi(1)

            stmp = n(i,j,:)
            stmq = n(i,j-1,:)

            if (model .eq. 0 .or. model .eq. 1) then
               s1 = 0
               s2 = 0
               do k = 1,ncomps
                  s1 = s1 + stmp(k)/rhocomp(k)
                  s2 = s2 + stmq(k)/rhocomp(k)
               end do
               rho = 0
               do k = 1,ncomps
                  rho = rho + half*(stmp(k)/s1+stmq(k)/s2)
               end do
               vmac(i,j) = vmac(i,j) - kpy(i,j)*rho*abs(gravity)
               
            elseif (model .eq. 2) then  
               dpc = (pc(i,j) - pc(i,j-1))/dx(2)
               ld  = lbd(i,j,:)
               tmp(1) = gravity*(ld(1)*rhocomp(1)+ld(2)*rhocomp(2))-ld(1)*dpc
               ld  = lbd(i,j-1,:)
               tmp(2) = gravity*(ld(1)*rhocomp(1)+ld(2)*rhocomp(2))-ld(1)*dpc
               rho = (tmp(1)+tmp(2))*half

               if (j .eq. dom_lo(2) .and. sbc(2,1) .ne. INT_DIR) then 
                  dpc = dpc*two
                  ld  = lbd(i,j-1,:)
                  rho = gravity*(ld(1)*rhocomp(1)+ld(2)*rhocomp(2))-ld(1)*dpc
               elseif (j .eq. dom_hi(2)+1 .and. sbc(2,2) .ne. INT_DIR) then
                  dpc = dpc*two
                  ld = lbd(i,j,:)
                  rho = gravity*(ld(1)*rhocomp(1)+ld(2)*rhocomp(2))-ld(1)*dpc
               end if

               kedge = kpy(i,j)*rho
               vmac(i,j) = vmac(i,j) - kedge

            elseif (model .eq. 3) then
c           NEED REVIEW
               call get_mu(mu,stmp,mucomp,compcount,compidx,nphases,ncomps)
               do l = 0,1
                  do nc = 1,nphases
                     if (compcount(nc) .eq. 2) then
                        s1 = n(i,j-l,compidx(nc,1))/rhocomp(compidx(nc,1)) + 
     &                       n(i,j-l,compidx(nc,2))/rhocomp(compidx(nc,2))
                        rho = half*s1**2*rhocomp(compidx(nc,1))/ mu(nc)
                     else if (compcount(nc) .eq. 1) then
                        s1 = n(i,j-l,compidx(nc,1))/rhocomp(compidx(nc,1))
                        rho = rho + half*s1**2*rhocomp(compidx(nc,1))/mu(nc)
                     end if
                  end do
               end do

               kedge = kpy(i,j)
               vmac(i,j) = vmac(i,j) - kedge*rho*abs(gravity)

            end if         
         end do
      end do

c ::: Adding inhomogenous term due to density variation
      hv_diffuse= .false.
      do k = 1,ncomps
         if (const_diff_coef(k) .gt. zero) then
            hv_diffuse = .true.
         end if
      end do
         
      if (hv_diffuse) then
c     ::: In y direction
         do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)
               rho  = zero
               rho1 = zero
               rho2 = zero
               do k = 1,ncomps
                  rho  = rho  + half*(n(i,j,k)+n(i,j-1,k))
                  rho1 = rho1 + n(i,j,k)
                  rho2 = rho2 + n(i,j-1,k)
               end do

               do k = 1,ncomps
                  RhoDY(i,j) = RhoDY(i,j) +
     &                 0.5*(phi(i,j)+phi(i,j-1))*const_diff_coef(k)*rho*
     &                 (n(i,j,k)/rho1 - n(i,j-1,k)/rho2)/(dx(2)*rhocomp(k))   
               end do
               
            end do
         end do
         
c     ::: In x direction         
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)+1
               rho  = zero
               rho1 = zero
               rho2 = zero
               do k = 1,ncomps
                  rho  = rho  + half*(n(i,j,k)+n(i-1,j,k))
                  rho1 = rho1 + n(i,j,k)
                  rho2 = rho2 + n(i-1,j,k)
               end do
               
               do k = 1,ncomps
                  RhoDX(i,j) = RhoDX(i,j) + 
     &                 0.5*(phi(i,j)+phi(i-1,j))*const_diff_coef(k)*rho*
     &                 (n(i,j,k)/rho1 - n(i-1,j,k)/rho2)/(dx(1)*rhocomp(k))
               end do
            end do
         end do
      end if

c     Dirichlet BC for the pressure equation
      if (pbc(1,1).eq. EXT_DIR .and.lo(1) .eq. dom_lo(1)) then
         i = lo(1)-1
         if (sbc(1,1) .eq. SEEPAGE .or. sbc(1,1) .eq. HOEXTRAP) then
            do j = lo(2),hi(2)+1
               stmp = zero
               y = (dble(j)+half)*dx(2)
               if (y .lt. wt_lo) then
                  stmp(1) = rhocomp(1)
               else
                  stmp(1) = rhocomp(2)
               end if
               if (y - dx(2) .lt. wt_lo) then
                  stmp(2) = rhocomp(1)
               else
                  stmp(2) = rhocomp(2)
               end if
               RhoG(lo(1),j)= -half*(stmp(1)+stmp(2))*gravity
            end do
         else
            do j = lo(2),hi(2)+1
               dpc = (pc(i,j) - pc(i,j-1))/dx(2)
               if (j .eq. dom_lo(2)) then
c                 based on dom_lo(2) data
                  ld = lbd(i,j,:)
                  kedge = gravity*(ld(1)*rhocomp(1) + ld(2)*rhocomp(2)) - ld(1)*dpc*2d0
                  ledge = ld(1) + ld(2)
               elseif (j .eq. dom_hi(2)+1) then
c                 basd on dom_hi(2) data
                  ld = lbd(i,j-1,:)
                  kedge = gravity*(ld(1)*rhocomp(1) + ld(2)*rhocomp(2)) - ld(1)*dpc*2d0 
                  ledge = ld(1) + ld(2)
               else
                  ld = lbd(i,j,:)
                  tmp(1) = gravity*(ld(1)*rhocomp(1) + ld(2)*rhocomp(2)) - ld(1)*dpc
                  ledge = ld(1) + ld(2)
                  ld = lbd(i,j-1,:)
                  tmp(2) = gravity*(ld(1)*rhocomp(1) + ld(2)*rhocomp(2)) - ld(1)*dpc
                  kedge = (tmp(1)+tmp(2))*half 
                  ledge = (ld(1) + ld(2) + ledge)*half
               end if
               RhoG(lo(1),j)= RhoG(lo(1),j)-kedge/ledge
            end do
         end if
      end if  

      if (pbc(1,2).eq. EXT_DIR .and. hi(1) .eq. dom_hi(1)) then
         i = hi(1)+1
         if (sbc(1,2) .eq. SEEPAGE .or. sbc(1,2) .eq. HOEXTRAP) then
            do j = lo(2),hi(2)+1
               stmp = zero
               y = (dble(j)+half)*dx(2)
               if (y .lt. wt_hi) then
                  stmp(1) = rhocomp(1)
               else
                  stmp(1) = rhocomp(2)
               end if
               if (y - dx(2) .lt. wt_hi ) then
                  stmp(2) = rhocomp(1)
               else
                  stmp(2) = rhocomp(2)
               end if
               RhoG(hi(1),j)= -half*(stmp(1)+stmp(2))*gravity
            end do
         else
            do j = lo(2),hi(2)+1
               dpc = (pc(i,j) - pc(i,j-1))/dx(2)
               if (j .eq. dom_lo(2)) then
                  ld = lbd(i,j,:)
                  kedge = gravity*(ld(1)*rhocomp(1) + ld(2)*rhocomp(2)) - ld(1)*dpc*2d0
                  ledge = ld(1) + ld(2)
               elseif (j .eq. dom_hi(2)+1) then
                  ld = lbd(i,j-1,:)
                  kedge = gravity*(ld(1)*rhocomp(1) + ld(2)*rhocomp(2)) - ld(1)*dpc*2d0 
                  ledge = ld(1) + ld(2)
               else
                  ld = lbd(i,j,:)
                  tmp(1) = gravity*(ld(1)*rhocomp(1) + ld(2)*rhocomp(2)) - ld(1)*dpc
                  ledge = ld(1) + ld(2)
                  ld = lbd(i,j-1,:)
                  tmp(2) = gravity*(ld(1)*rhocomp(1) + ld(2)*rhocomp(2)) - ld(1)*dpc
                  kedge = (tmp(1)+tmp(2))*half 
                  ledge = (ld(1) + ld(2) + ledge)*half
               end if
               RhoG(hi(1),j)= RhoG(hi(1),j)-kedge/ledge
            end do
         end if
      end if

c ::: Set no slip b.c / inflow velocity b.c.
      if (pbc(1,1).eq. FOEXTRAP .and. lo(1) .eq. dom_lo(1)) then
           do j = lo(2),hi(2)
              if (inflow_bc_lo(1) .eq. 1) then
                 umac(lo(1),j) = inflow_vel_lo(1)
              else
                 umac(lo(1),j) = zero
              end if
           end do
      end if

      if (pbc(1,2).eq. FOEXTRAP .and. hi(1) .eq. dom_hi(1)) then
           do j = lo(2),hi(2)
              if (inflow_bc_hi(1) .eq. 1) then
                 umac(hi(1)+1,j) = inflow_vel_hi(1)
              else
                 umac(hi(1)+1,j) = zero
              end if
           end do
      end if

      if (pbc(2,1).eq. FOEXTRAP .and. lo(2) .eq. dom_lo(2)) then
           do i = lo(1),hi(1)
              if (inflow_bc_lo(2) .eq. 1) then
                 vmac(i,lo(2)) = inflow_vel_lo(2)
              else
                 vmac(i,lo(2)) = zero
              end if
           end do
      end if

      if (pbc(2,2).eq. FOEXTRAP .and. hi(2) .eq. dom_hi(2)) then
           do i = lo(1),hi(1)
              if (inflow_bc_hi(2) .eq. 1) then
                 vmac(i,hi(2)+1) = inflow_vel_hi(2)
              else
                 vmac(i,hi(2)+1) = zero
              end if
           end do
      end if

      end

c ::: -----------------------------------------------------------
c ::: User-specified kappa distirbution
c ::: -----------------------------------------------------------
      subroutine FORT_USER_KAPPA(kappa,DIMS(kappa),
     &     kpval,nkpval,dkp,dx,problo,probhi)

      integer DIMDEC(kappa)
      integer nkpval
      REAL_T  dkp
      REAL_T  kappa(DIMV(kappa))
      REAL_T  kpval(nkpval)
      REAL_T  dx(SDIM)
      REAL_T  problo(SDIM),probhi(SDIM)

      integer i,j

      do j=ARG_L2(kappa),ARG_H2(kappa)
         do i=ARG_L1(kappa),ARG_H1(kappa)
               kappa(i,j) = kpval(1)
         end do
      end do

      end

c ::: -----------------------------------------------------------
c ::: This routine tags cells for refinement based on the
c ::: kappa values.
c ::: -----------------------------------------------------------
      subroutine FORT_KPERROR (tag,DIMS(tag),
     &     kappa,DIMS(kappa),set,val,dx, problo,probhi,
     &     lo,hi,domlo,domhi,level)
      
      integer   DIMDEC(tag)
      integer   DIMDEC(kappa)
      integer   set
      integer   tag(DIMV(tag))
      integer   lo(SDIM),hi(SDIM)
      integer   domlo(SDIM),domhi(SDIM)
      integer   level
      REAL_T    kappa(DIMV(kappa))
      REAL_T    val
      REAL_T    dx(SDIM)
      REAL_T    problo(SDIM),probhi(SDIM)

      integer   i, j
      REAL_T    grad

      do i = lo(1), hi(1)
         do j = lo(2), hi(2)
           grad = max(dabs(kappa(i+1,j)-kappa(i-1,j))/(2*dx(1)),
     &            dabs(kappa(i,j+1)-kappa(i,j-1))/(2*dx(2))) 
           tag(i,j) = merge(set,tag(i,j),grad .gt. val)
        end do
      end do
      
      end


c ::: --------------------------------------------------------------------
c ::: This routine add the user-defined forcing terms to the flow equation
c ::: --------------------------------------------------------------------
      subroutine FORT_MK_FORCE(time,
     &     force,DIMS(force),
     &     dx,problo,probhi,scomp,ncomp,
     &     do_rho_scale)

      implicit none

#include "probdata.H"

      integer    do_rho_scale
      integer    DIMDEC(force)
      integer    scomp, ncomp
      REAL_T     time, dx(SDIM)
      REAL_T     problo(SDIM), probhi(SDIM)
      REAL_T     force(DIMV(force),scomp:scomp+ncomp-1)


c     local variables
      integer i, j
      REAL_T  x, y

      do i = ARG_L1(force),ARG_H1(force)
         x = (dble(i)+half)*dx(1)
         do j = ARG_L2(force),ARG_H2(force)
            y = (dble(j)+half)*dx(2)

            force(i,j,:) = zero

c           infiltration
            if (y .gt. probhi(2)-dx(2)) then
               force(i,j,scomp) = force(i,j,scomp) + 1.18398d-5/dx(2)
            end if

c           discharge from basin
            if (y .gt. probhi(2)-dx(2) .and. 
     &          x .ge. 550d0 .and. x .lt. 650d0) then
               force(i,j,scomp) = force(i,j,scomp) + 6.02486d-4/dx(2)
            end if

            if (do_rho_scale .eq. 1) then
               force(i,j,scomp) = force(i,j,scomp)/rhoval(1)
            end if

         end do
      end do

      end

c ::: ------------------------------------------------------------------
c ::: This routine add user-defined forcing terms to the tracer equation
c ::: ------------------------------------------------------------------
      subroutine FORT_MK_TFORCE(time,
     &     force,DIMS(force),
     &     dx,problo,probhi,scomp,ncomp,
     &     do_rho_scale)
      
      implicit none

#include "probdata.H"

      integer    do_rho_scale
      integer    DIMDEC(force)
      integer    scomp, ncomp
      REAL_T     time, dx(SDIM)
      REAL_T     problo(SDIM), probhi(SDIM)
      REAL_T     force(DIMV(force),scomp:scomp+ncomp-1)


c     local variables
      integer i, j, nc
      REAL_T  x, y, rhoinflow

      rhoinflow = zero
 
      do i = ARG_L1(force),ARG_H1(force)
         x = (dble(i)+half)*dx(1)
         do j = ARG_L2(force),ARG_H2(force)
            y = (dble(j)+half)*dx(2)

           force(i,j,:) = zero

c           discharge from basin
            if (y .gt. probhi(2)-dx(2) .and. 
     &          x .ge. 550d0 .and. x .lt. 650d0) then
               do nc = scomp,scomp+ncomp-1
                  force(i,j,nc) = force(i,j,nc) + 
     &                 6.02486d-7*rhoinflow/dx(2)
               end do
            end if

            force(i,j,:) = zero

         end do
      end do

      end
