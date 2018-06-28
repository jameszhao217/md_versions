	subroutine getEnergiesAndForcesNoCl (id, atomCoord, ix, f, 
     &		atomDispl, aveDispl, atomForce, systemEnergy, 
     &		MoveAtoms, MoveDisl, FullField, solveFEM, strainE0)

	use mod_global
        use mod_timming
	implicit none

	double precision atomDispl(ndf,*), atomCoord(ndf,*),f(ndf,*),
     &		atomForce(ndf,*), systemEnergy, aveDispl(ndf,*)
	integer  id(ndf,*), ix(nen1,*)
	logical MoveAtoms,MoveDisl, FullField, solveFEM
	integer iAtom
        double precision strainE0

	if (solveFEM .eq. .true.) then

!!	   Solve FEM
           call vafuncMD(id, atomCoord, ix, f, aveDispl, 
     &	    atomForce, systemEnergy, MoveAtoms, MoveDisl, 
     &      FullField, strainE0)

!!	   Get Displacements on CONTINUUM and PAD atoms:
	   do iAtom = 1, numnp

   	     if (isRelaxed(iAtom) .eq. indexContinuum
     &	     .or. isRelaxed(iAtom) .eq. indexPad) then

 	       atomDispl(1:ndf, iAtom) = aveDispl(1:ndf, iAtom)

	     endif

	   enddo

	endif

	return
	end

 !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc

	subroutine getEnergiesAndForces (id, atomCoord, ix, f, 
     &		atomDispl, aveDispl, atomForce, systemEnergy, 
     &		MoveAtoms, MoveDisl, FullField, solveFEM, strainE0)

	use mod_global
        use mod_timming
	implicit none

	double precision atomDispl(ndf,*), atomCoord(ndf,*),f(ndf,*),
     &		atomForce(ndf,*), systemEnergy, aveDispl(ndf,*)
	integer  id(ndf,*), ix(nen1,*)
	logical MoveAtoms,MoveDisl, FullField, solveFEM
	integer iAtom, i, j, wcnt, wcnt_max, clsr_cnt,ctod_pair_i
        double precision strainE0, dErr, dErrTol, systemEnergy_last
        double precision :: clMg,partDist,fmag,f_finite,dx,d_use
        double precision :: oxyThick_f,ctod_pair_x
	   
        integer, allocatable :: clsr_n(:)
	   ! clsr_n ! Stores closure node #

        double precision, allocatable :: clsr_x(:)
        double precision, allocatable :: clsr_d(:)
        double precision, allocatable :: clsr_k(:)
        double precision, allocatable :: clsr_l(:)
	   ! clsr_x ! Stores original position
	   ! clsr_d ! Stores closure "overlap"
	   ! clsr_k ! Stores closure "stiffness"
	   ! clsr_l ! Stores last "force"

	f_finite = closureMag
	dErrTol = 0.5d0
	wcnt = 0
	wcnt_max = 25

	if (solveFEM .eq. .true.) then

	 allocate(clsr_l(numnp)) ! Stores last "force"
	 clsr_l(1:numnp) = 0.d0

         ! Zeroth:
         if(firstAt)then
          do i = 1, numnp
	   if(id_closure(i).eq.1)then
	    f(2,i) = 0.d0 ! Make sure there is no "closure" force set
	   endif
          enddo
	  firstAt = .false.
	 endif

	 dErr = dErrTol+1.d0
	 do while ((dErr.gt.dErrTol).and.(wcnt.lt.wcnt_max))   

	  ! First: Run FEM
          call vafuncMD(id, atomCoord, ix, f, aveDispl, 
     &	    atomForce, systemEnergy, MoveAtoms, MoveDisl, 
     &      FullField, strainE0)
	  do iAtom = 1, numnp
   	    if (isRelaxed(iAtom) .eq. indexContinuum
     &	     .or. isRelaxed(iAtom) .eq. indexPad) then
 	       atomDispl(1:ndf, iAtom) = aveDispl(1:ndf, iAtom)
	    endif
	  enddo

	  ! Second: Search for "closure" nodes that are too close
	  clsr_cnt = 0
	  dErr = 0.d0
	  ctod_pair_i = 0
	  ctod_pair_x = -1000000.d0
	  do i = 1, numnp
	    if((id_closure(i).eq.1).and.
     &			(atomCoord(1,i).le.-closureDist))then
	     if(nodePartner(i).ne.0)then
	       partDist = (atomCoord(2,i)+atomDispl(2,i))
     &	 	-(atomCoord(2,nodePartner(i))+
     &			atomDispl(2,nodePartner(i)))

	       if(atomCoord(2,i).gt.0.d0)then
	        if(partDist.lt.(2.d0*oxyThick))then
		 dErr = dErr + ((2.d0*oxyThick-partDist)*0.5)**2
		 clsr_cnt = clsr_cnt+1
	         clsr_l(i) = f(2,i)
	        endif
	       else
	        if(partDist.gt.(-2.d0*oxyThick))then
		 dErr = dErr + ((partDist+2.d0*oxyThick)*0.5)**2
		 clsr_cnt = clsr_cnt+1
	         clsr_l(i) = f(2,i)
	        endif
	       endif

	       if(atomCoord(2,i).gt.0.d0)then
	         if(atomCoord(1,i).gt.ctod_pair_x)then
	           ctod_pair_i = i
	           ctod_pair_x = atomCoord(1,i)
	         endif
	       endif

	     endif
	    endif
	  enddo
	  dErr = sqrt(dErr)
	  !if(rank.eq.0) write(*,*) 'dErr = ',dErr 

	  if(ctod_pair_i.gt.0)then
	   ctod = (atomCoord(2,ctod_pair_i)+atomDispl(2,ctod_pair_i))
     &	 	-(atomCoord(2,nodePartner(ctod_pair_i))+
     &			atomDispl(2,nodePartner(ctod_pair_i)))
	  endif
	
	  ! If there are "contacting" closure nodes...
	  ! 		Need to find Approx. stiffness of each node.
	  if(clsr_cnt.gt.0)then

!	   if(rank.eq.0) write(*,*) 'Closure!!.. 
!     &				Need to rerun FEM a few times..'
!	   if(rank.eq.0) write(*,*) 'dErr = ',dErr 
!	   if(rank.eq.0) write(*,*) 'wcnt = ',wcnt

	   allocate(clsr_x(clsr_cnt)) ! Stores closure node #
	   allocate(clsr_n(clsr_cnt)) ! Stores closure node #
	   allocate(clsr_d(clsr_cnt)) ! Stores closure "overlap"
	   allocate(clsr_k(clsr_cnt)) ! Stores closure "stiffness"
!	   if(rank.eq.0) write(*,*) 'clsr_cnt = ',clsr_cnt

	   clsr_cnt = 0
	   do i = 1, numnp
	    if((id_closure(i).eq.1).and.
     &			(atomCoord(1,i).le.-closureDist))then
	     if(nodePartner(i).ne.0)then
	       partDist = (atomCoord(2,i)+atomDispl(2,i))
     &	 	-(atomCoord(2,nodePartner(i))+
     &			atomDispl(2,nodePartner(i)))
               if((atomCoord(2,i).gt.0.d0).and.
     &		(partDist.lt.(2.d0*oxyThick)))then
		 clsr_cnt = clsr_cnt+1
		 clsr_n(clsr_cnt) = i
		 clsr_d(clsr_cnt) = (2.d0*oxyThick-partDist)*0.5
		 clsr_x(clsr_cnt) = atomCoord(2,i)+atomDispl(2,i)
	       elseif((atomCoord(2,i).lt.0.d0).and.
     &		(partDist.gt.(-2.d0*oxyThick)))then
		 clsr_cnt = clsr_cnt+1
		 clsr_n(clsr_cnt) = i
		 clsr_d(clsr_cnt) = -(partDist+2.d0*oxyThick)*0.5
		 clsr_x(clsr_cnt) = atomCoord(2,i)+atomDispl(2,i)
	       endif
	     endif
	    endif
	   enddo
!	   if(rank.eq.0) write 'clsr_cnt = ',clsr_cnt

	   ! Rerun FEM to get "Stiffnesses"
	   do i = 1,clsr_cnt
	     if(atomCoord(2,clsr_n(i)).gt.0.d0)then
	       f(2,clsr_n(i))=f_finite
     &		*abs(nodeArea(clsr_n(i)))/time+clsr_l(clsr_n(i))
	     else
	       f(2,clsr_n(i))=-f_finite
     &		*abs(nodeArea(clsr_n(i)))/time+clsr_l(clsr_n(i))
	     endif
	   enddo

           call vafuncMD(id, atomCoord, ix, f, aveDispl, 
     &	     atomForce, systemEnergy, MoveAtoms, MoveDisl, 
     &       FullField, strainE0)
	   do iAtom = 1, numnp
   	     if (isRelaxed(iAtom) .eq. indexContinuum
     &	      .or. isRelaxed(iAtom) .eq. indexPad) then
 	       atomDispl(1:ndf, iAtom) = aveDispl(1:ndf, iAtom)
	     endif
	   enddo

	   do i = 1, clsr_cnt
	     dx = (atomCoord(2,clsr_n(i))
     &			+atomDispl(2,clsr_n(i)))-clsr_x(i)
	     if(abs(dx).gt.(0.001d0))then
	       if(atomCoord(2,clsr_n(i)).gt.0.d0)then
		  clsr_k(i)=(f_finite*nodeArea(clsr_n(i))/time)/dx
	       else
		  clsr_k(i)=(-f_finite*nodeArea(clsr_n(i))/time)/dx
	       endif
	     else
	        clsr_k(i) = 0.d0
	     endif
	   enddo

	   ! Set appropriate forces:
	   do i = 1,clsr_cnt
	     f(2,clsr_n(i)) = clsr_l(clsr_n(i))
	     f(2,clsr_n(i)) = f(2,clsr_n(i)) + clsr_k(i)*clsr_d(i)
	   enddo

	   ! deallocate temp arrays:
	   deallocate(clsr_x) 
	   deallocate(clsr_n) ! Stores closure node #
	   deallocate(clsr_d) ! Stores closure "overlap"
	   deallocate(clsr_k) ! Stores closure "stiffness"

	   ! Go to start of while loop to test out new forces...

	  else

	    dErr=dErrTol-1.d0

	  endif ! IF clsr_cnt .gt. 0

	  wcnt=wcnt+1
	 enddo !WHILE

	 deallocate(clsr_l) ! Stores last "force"

	endif ! IF FEM

	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
	
	double precision function rotate3D(a,b,c,x,y,z,outdir)

	implicit none
        double precision newx, a, b, c, x, y, z
	integer outdir

        newx=0.0
        if (outdir.eq.0) then
          rotate3D = x*(cos(b)*cos(c)) + y*(cos(b)*sin(c)) + 
     &		     z*(-sin(b))
        else if (outdir.eq.1) then
          rotate3D = x*(sin(a)*sin(b)*cos(c)-cos(a)*sin(c)) + 
     &		     y*(sin(a)*sin(b)*sin(c)+cos(a)*cos(c)) + 
     &		     z*(sin(a)*cos(b))
        else if (outdir.eq.2) then
          rotate3D = x*(cos(a)*sin(b)*cos(c)+sin(a)*sin(c)) + 
     &		     y*(cos(a)*sin(b)*sin(c)-sin(a)*cos(c)) + 
     &		     z*(cos(a)*cos(b))
        else
            write(*,*)'ERROR... rotate3D takes 0,1 or 2 ONLY!'
        endif 

	return
	end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc

	subroutine dosteps(n,atomDispl,atomForce,db,fl,epps,
     &	  iprint,dsmax,dfn,atomID,atomCoord,ix,f,itx,
     &	  CheckSlip,AddedSlip,LostSlip,MoveDisl,MoveAtoms) 

	use mod_global
	use mod_crack
	use mod_poten
	use mod_dynamo
	use mod_output
        use mod_material
        use mod_timming
        use lammps
        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int

        implicit none

        include 'environment.h' ! <- interfaces for passing allocatable arrays
	include 'common.f'

	! Might remove later...
	integer :: maxAtoms
	parameter (maxAtoms = 180000)
	double precision :: systemEnergy

!--	Passed variables
	integer :: n, ix(*), atomID(ndf,*), iprint, itx(*)
	double precision :: atomDispl(ndf, numnp),db(n)
	double precision :: atomForce(ndf,*),fl,epps,dsmax,dfn
	double precision :: atomCoord(ndf,*),f(*)
	logical :: AddedSlip, CheckSlip, MoveDisl,LostSlip,MoveAtoms

!--	Local variables
        !type (C_ptr) :: lmp

        integer, allocatable :: type_lmp(:)
        integer, allocatable :: type_Imp(:)
        integer, allocatable :: type_Imp_temp(:)
        integer, allocatable :: CV_cell_free(:,:)
        integer, allocatable :: CV_list(:,:)

!ER- I made these to split Al into top and bottom regions
        integer, allocatable :: CV_cell_top(:,:)
        integer, allocatable :: CV_cell_bot(:,:)
        integer, allocatable :: Atom_top(:)
        integer, allocatable :: Atom_bot(:)

        double precision, allocatable :: x_lmp(:)
        double precision, allocatable :: v_lmp(:)
        double precision, allocatable :: f_lmp(:)
        double precision, allocatable :: q_lmp(:)
        double precision, allocatable :: x_Imp(:)
        double precision, allocatable :: v_Imp(:)
        double precision, allocatable :: x_Imp_temp(:)
        double precision, allocatable :: v_Imp_temp(:)
        double precision, allocatable :: maxforcelist(:)
        double precision, allocatable :: dis_12(:)

        integer :: jj,iii,ii,i_cntr
        integer :: loopsPerQeq
        integer :: irep, i, j, iAtom, Nsteps, FEMSteps, FEMStepMin
        integer :: padIDmin, padIDmax, n_remove, n_temp, zr, rz
        integer :: n_O, n_H, n_Imp, dumpcg, n_h2o_need, n_o2_need
        integer :: n_atoms, n_atoms_lmp, stepsPerFE, velFile, n_CV
	integer :: ndis_checked, MDSteps, nnsteps, FEMStepRange
	integer :: iStep, FEMStepCounter,n_new_Imp,itype,n_Al

	double precision :: cx,cy,cz,o2_conc,h2o_conc,n_surf
	double precision :: o2_conc_max, charge_max, charge_min
	double precision :: fx,fy,fz,dfx,dfy,dfz,fmag,kB
        double precision :: monothick_max,deltK,cycleStart
	double precision :: fx1_i,fy1_i,fz1_i,fx1_f,fy1_f,fz1_f
        double precision :: interval_x,interval_y,interval_z,vstart
	double precision :: timestep,mdTemp,ForceMax,ForceMaxP,fgrad
	double precision :: x_cracktip, y_cracktip,adiff_T,adiff_n
	double precision :: aveDispl(ndf,maxAtoms),strainE0,Ftol
	double precision :: cx1,cx2,cy1,cy2,cz1,cz2,dx,dy,dz,dmag 
	double precision :: dmagC,dxC,dyC,randomNum,K_max,conc_scale

        logical :: NeedList,EAMdebug,switchok
	logical :: FullField,solveFEM,reStartAveraging
        logical :: dchk,redistribute,keeprunning,newfix
        logical :: quickminOn,hftnOn,writefiles

	character *1000 :: str001,str002,str003

!--	Functions
	integer dislcheck
	double precision rotate3D

	data MDSteps /0/
	data nnsteps /0/
        data FEMStepCounter /0/
        data ndis_checked /0/

        if (maxAtoms.lt.numnp) then
	  print*, 'numnp: ', numnp
          print*,'increase size of maxAtoms in md.f'
          stop
        endif
	
	reStartAveraging = .false.
	firstAt = .true.
	
	call lostslipinit(LostSlip)
	
	! Some Assignements..
	FEMSteps = 1
	FullField = .true.
	FEMStepMin = 5
	FEMStepRange = 0
	newfix = .false.
	switchok = .false.
        quickminOn = .false.
	hftnOn = .false.
	charge_max = 0.8870
	charge_min = -0.8138 

	! WRITE FILE TO CHECK:
!	if(rank.eq.0) then
!	  open(unit=2117, file='dump_CHECK_4', status='unknown')
!	  do i = 1, numnp
!	    if(id_closure(i).eq.1)then
!	      write(2117,'(F16.8,1X,F16.8,1X,F16.8,1X,I5,1X,I5)') 
!     &		atomCoord(1,i)+atomDispl(1,i),
!     &		atomCoord(2,i)+atomDispl(2,i),nodeArea(i),
!     &		id_closure(i),nodePartner(i)
!	    else
!	      write(2117,'(F16.8,1X,F16.8,1X,F16.8,1X,I5,1X,I5)') 
!     &		atomCoord(1,i)+atomDispl(1,i),
!     &		atomCoord(2,i)+atomDispl(2,i),nodeArea(i),
!     &		id_closure(i),nodePartner(i)
!	    endif
!          enddo
!	  close(2117)
!	endif
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	! Read in 'md.inp' input file:
!	if(rank.eq.0) write(*,*)"call ReadMDinp()"
        call ReadMDinp(timestep,FEMStepMin,Nsteps,
     &	o2_conc,h2o_conc,stepsPerFE,dumpcg,EAMdebug,
     &	mdTemp,Ftol,K_max)

	writefiles = .false.
	if(MDSteps.eq.0)then
          loadChange = 1.d0 ! Assumes load increment is positive in beginning
          loadLast = time
          nFCGcycle = 0.5
          writefiles = .true.
	  if(rank.eq.0)then
	   open(unit=1115,file='da_ctod.out',status='replace')
	   write(1115,*) 'MDSteps,x_cracktip,ctod -> '
	   close(1115)
	  endif
	else
	  if(time.ne.loadLast)then
	     deltK = time - loadLast
	     if(deltK*loadChange.lt.(0.d0))then
		nFCGcycle = nFCGcycle + 0.5
	        loadChange = -loadChange
		writefiles = .true.
	     endif
	     loadLast = time
	  endif
	endif
	if(rank.eq.0) write(*,*)"Doing FCG cycle #",nFCGcycle

	if(oxyThick_max.gt.0.d0)then
	  cycleStart = 3.5
	else
	  cycleStart = 0.5
	endif

	! Adjust o2_conc to increase linearly with load (if EAM):
	if(startIMP)then

	 o2_conc_max = o2_conc
	 monoThick_max = monoThick

         if(Monolayer)then
	  if(time.le.(K_max*fstart))then
            conc_scale = 0.d0
	  else
            conc_scale = (time-(K_max*fstart)) / (K_max*(1.d0-fstart))
	  endif
          if(conc_scale.gt.(1.d0)) conc_scale = 1.d0
          if(conc_scale.lt.(0.d0)) conc_scale = 0.d0
          o2_conc = o2_conc_max * conc_scale
	  monoThick = monoThick_max * conc_scale + 2.0*grid_dx
         endif

	 ! Adjust Crack-Closure Params...
         conc_scale = time / K_max
	 if(conc_scale.lt.0.d0) conc_scale = 0.d0
	 if(oxythick.lt.(oxyThick_max * conc_scale))then
	  oxyThick = oxyThick_max * conc_scale
	 endif
	 oxyThick = max(min(r_crack,2.d0),oxyThick)

	else

          o2_conc = 0.d0
	  if(nFCGcycle.ge.(cycleStart))then ! <- Change cycle to start adding environment HERE.
	    startIMP = .true.
	  endif
	  monoThick = 2.0*grid_dx 
	  oxyThick = min(r_crack,2.d0)

	endif

	if(rank.eq.0) write(*,*) 'o2_conc currently: ',o2_conc
	if(rank.eq.0) write(*,*) 'oxyThick currently: ',oxyThick
	if(rank.eq.0) write(*,*) 'monoThick currently: ',monoThick
	if(rank.eq.0) write(*,*) 'startIMP currently: ',startIMP
	if(rank.eq.0) print*, 'MDSteps: ', MDSteps

        ! Load Existing Impurities:
!	if(rank.eq.0) write(*,*)"call ReadIMP()"
	call ReadIMP(velFile,n_O,n_H,n_Imp,x_Imp,v_Imp,type_Imp,
     &	redistribute,adiff_T,adiff_n)

	! Set up control volume (CV):
!	if(rank.eq.0) write(*,*)"call CV_setup()"
        call CV_setup(atomCoord,atomDispl,n_CV,
     &		CV_cell_free,CV_list,dis_12,CV_cell_top,CV_cell_bot,
     &		Atom_top,Atom_bot)
!	if(rank.eq.0) write(*,*)"CV_setup() Done."
!        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!        if(rank.eq.0)then
!          open(unit=6116, file='dump_Surface', status='replace')
!	  i=0
!          do irep = 1, numnp
!	   if(isRelaxed(irep).ne.0)then
!	    i=i+1
!            write(6116,'(F16.8,1X,F16.8,1X,F16.8,
!     &      1X,F16.8,1X)')
!     &      atomCoord(1,irep)+atomDispl(1,irep),
!     &      atomCoord(2,irep)+atomDispl(2,irep),
!     &      atomCoord(3,irep)+atomDispl(3,irep),
!     &      dis_12(i)
!	   endif
!          enddo
!          close(6116)
!        endif

	! Find Crack-tip location...
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	if(rank.eq.0) write(*,*)"call find_crack_tip()"
	call find_crack_tip(n_CV,CV_list,x_cracktip,y_cracktip,
     &	EAMdebug,nnsteps,atomCoord,atomDispl,dis_12)

	! Plan number of impurity molecules to add:
	n_surf = 0.d0
!	if(rank.eq.0) write(*,*)"call plan_IMP()"
        call plan_IMP(n_remove,n_Imp,
     &	CV_cell_free,x_Imp,type_Imp,
     &	n_o2_need,n_h2o_need,o2_conc,h2o_conc,n_CV,EAMdebug,
     &  MDSteps,atomCoord,atomDispl,n_surf,x_cracktip,dis_12)

	! Define Interval where impurities will go:
    	interval_x = grid_dx 
    	interval_y = grid_dx 
    	interval_z = z_length

	! Allocate temp arrays:
	! n_temp = 2*n_o2_need+3*n_h2o_need+n_Imp-n_remove
	if(Monolayer)then
	  n_h2o_need = 0 ! and n_o2_need is number of "single" o atoms to add (not O2)
	  n_temp = n_o2_need+n_Imp
	else
	  n_temp = 2*n_o2_need+3*n_h2o_need+n_Imp
	endif
	if(allocated(x_Imp_temp))deallocate(x_Imp_temp)
	if(allocated(v_Imp_temp))deallocate(v_Imp_temp)
	if(allocated(type_Imp_temp))deallocate(type_Imp_temp)
	allocate(x_Imp_temp(n_temp*3))
	allocate(v_Imp_temp(n_temp*3))
	allocate(type_Imp_temp(n_temp))

        ! Figure Out Reasonable Velocity to Send Impurities Toward Crack-Tip
	kB = 8.314462175 ! Boltzmann Const. [J/mol/K]
	vstart = sqrt( 3.0 * kB * mdTemp / 15.9994 ) ! [m/s]
	vstart = vstart * 0.00001 ! [Ang/fs]
	vstart = 0.d0 ! <- Input atoms at rest for now.

	i_cntr=1 ! <- START COUNTER

	 ! Generate O2 molecules:
!	 if(rank.eq.0) write(*,*)"call add_O2()"
	 call add_O2(i_cntr,vstart,CV_list,n_CV,interval_x,
     &		interval_y,interval_z, x_Imp,x_Imp_temp,
     &		v_Imp_temp,type_Imp_temp,n_o2_need,n_Imp,
     &		EAMdebug,n_surf,atomCoord,atomDispl,
     &		x_cracktip,MDSteps,n_temp,dis_12)

	 ! Generate H2O molecules:
!	 if(rank.eq.0) write(*,*)"call add_H2O()"
	 call add_H2O(i_cntr,vstart,CV_list,n_CV,interval_x,
     &		interval_y,interval_z, x_Imp,x_Imp_temp,
     &		v_Imp_temp,type_Imp_temp,n_h2o_need,n_Imp)

	 ! Copy over temp arrays:
!	 if(rank.eq.0) write(*,*)"cp_temp_arrays()"
	 call cp_temp_arrays(i_cntr,n_Imp,x_Imp,v_Imp,type_Imp,
     &		x_Imp_temp,v_Imp_temp,type_Imp_temp,n_H,n_O)

	! Write Lammps Geometry Input File (called 'POSCAR' for now):
!	if(rank.eq.0) write(*,*)"write_geometry_input()"
        call write_geometry_input(n_atoms,n_Imp,padIDmin,
     &		padIDmax,velFile,EAMdebug,n_atoms_lmp,n_O,n_H,
     &		atomCoord,atomDispl,x_Imp,v_Imp,type_Imp)

	! Setup up Lammps:
        if(MDSteps.eq.0)then

	  ! If this is the first time we will be using Lammps...
	  ! Write Lammps "Instruction" Input File (in.lmp):
!	  if(rank.eq.0) write(*,*)"write_lammps_input()"
          call write_lammps_input(EAMdebug,padIDmin,
     &		padIDmax,mdTemp,timestep,dumpcg,MDSteps)
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	  if(EAMdebug.and.explicitOxy)then

	    ! Read in.lmp into 'lmp'
            call lammps_file (qeq, 'in.qeq')

	    str001=''
	    write(str001,*) 'run 0'
	    CALL lammps_command(qeq,str001)

	  endif

	  ! Read in.lmp into 'lmp'
          call lammps_file (lmp, 'in.lmp')
	  str001=''
	  write(str001,*) 'run 0'
	  CALL lammps_command(lmp,str001)

	  if(EAMdebug.and.explicitOxy.and.(o2_conc.gt.0.d0))then

	    ! Transfer Coords from 'lmp' to qeq'
            call lammps_gather_atoms (lmp, 'x', 3, x_lmp)
            call lammps_scatter_atoms (qeq, 'x', x_lmp)

	    ! run qeq to get chrages
	    str001=''
	    write(str001,*) 'run 1'
	    CALL lammps_command(qeq,str001)

	    ! Transfer charges from 'qeq' to 'lmp'
            call lammps_gather_atoms (qeq, 'q', 1, q_lmp)
	    do i = 1, n_atoms_lmp
	      if(q_lmp(i).gt.charge_max)then
		q_lmp(i) = charge_max
	      elseif(q_lmp(i).lt.charge_min)then
		q_lmp(i) = charge_min
	      endif
	    enddo
	    q_lmp(:) = q_lmp(:) * scaleChg
            call lammps_scatter_atoms (lmp, 'q', q_lmp)

	  endif

	  ! Run minimization in 'lmp'
	  str001=''
	  write(str001,*) 'minimize 0.0 ',Ftol,' 100 1000'
	  CALL lammps_command(lmp,str001)

          n_atoms_lmp = lammps_get_natoms(lmp)
          if(rank.eq.0)write(*,*)'#Atoms in Lammps: ',n_atoms_lmp

	  n_Al = n_atoms_lmp - n_Imp !

	else

	  ! If Lammps has already read in an input file (in.lmp)...
	  ! Adjust atom info, and create new impurities (if necessary):

	  n_Al = n_atoms_lmp - n_Imp ! <- Use 'n_atoms_lmp' from POSCAR write!

	  ! Now get the actual number of atoms that Lammps knows about:
          n_atoms_lmp = lammps_get_natoms(lmp)
!          if(rank.eq.0)write(*,*)'#Atoms in Lammps: ',n_atoms_lmp

	  ! Determine # of new impurity atoms to add to Lammps:
	  n_new_Imp = (n_Imp+n_Al)-n_atoms_lmp
          if(rank.eq.0) write(*,*) 'n_new_Imp=',n_new_Imp
          if(rank.eq.0) write(*,*) 'n_Imp=',n_Imp
          if(rank.eq.0) write(*,*) 'n_Al=',n_Al

	  ! Add new impurity atoms to Lammps:
	  if((n_new_Imp.gt.0).and.explicitOxy)then
	   do i = 1, n_new_Imp

	    irep = n_Imp - n_new_Imp + i

!            if(rank.eq.0) write(*,*) 'creating new atom...'

	    itype=type_Imp(irep)
	    cx = x_Imp((irep-1)*3+1)
	    cy = x_Imp((irep-1)*3+2)
	    cz = x_Imp((irep-1)*3+3)
	    if(cz.lt.(0.d0)) x_Imp((irep-1)*3+3) = cz + z_length
	    if(cz.gt.z_length) x_Imp((irep-1)*3+3) = cz - z_length
	    cz = x_Imp((irep-1)*3+3)

	    if(itype.eq.2)then
             write(str001,'(A21,1X,F16.8,1X,F16.8,1X,F16.8)') 
     &	'create_atoms 2 single',cx,cy,cz
	    else
             write(str001,'(A21,1X,F16.8,1X,F16.8,1X,F16.8)') 
     &	'create_atoms 3 single',cx,cy,cz
	    endif
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	    call lammps_command(lmp,str001)
	    if(EAMdebug) call lammps_command(qeq,str001)

	   enddo
	  endif

	  if(EAMdebug.and.explicitOxy.and.(n_Imp.gt.0))then
	    call lammps_gather_atoms (lmp, 'x', 3, x_lmp)
            call lammps_scatter_atoms (qeq, 'x', x_lmp)
	    call lammps_command(qeq,'run 1')
            call lammps_gather_atoms (qeq, 'q', 1, q_lmp)
            do i = 1, n_atoms_lmp
              if(q_lmp(i).gt.charge_max)then
                q_lmp(i) = charge_max
              elseif(q_lmp(i).lt.charge_min)then
                q_lmp(i) = charge_min
              endif
            enddo
            q_lmp(:) = q_lmp(:) * scaleChg
            call lammps_scatter_atoms (lmp, 'q', q_lmp)
	  endif
	  call lammps_command(lmp,'run 0')

	  ! How many atoms does lammps know about:
          n_atoms_lmp = lammps_get_natoms(lmp)
!          if(rank.eq.0)write(*,*)'New #Atoms in Lammps: ',n_atoms_lmp

	  ! Gather these atoms from Lammps:
          call lammps_gather_atoms (lmp, 'type', 1, type_lmp)
          call lammps_gather_atoms (lmp, 'x', 3, x_lmp)
          call lammps_gather_atoms (lmp, 'v', 3, v_lmp)

	  ! Adjust pad atoms for new load:
	  i=0
	  do irep = 1, numnp
            NeedList=IsRelaxed(irep).ne.0
            if(NeedList)then
	     i=i+1
             if(IsRelaxed(irep).eq.-1)then
	      x_lmp((i-1)*3+1) = atomCoord(1,irep)+atomDispl(1,irep)
	      x_lmp((i-1)*3+2) = atomCoord(2,irep)+atomDispl(2,irep)
	      x_lmp((i-1)*3+3) = atomCoord(3,irep)+atomDispl(3,irep)
	     endif
	    endif
	  enddo

	  ! Set all info for impurity atoms:
	  do irep = 1, n_Imp
	      i=irep+n_Al
	      type_lmp(i) = type_Imp(irep)
	      x_lmp((i-1)*3+1) = x_Imp((irep-1)*3+1)
	      x_lmp((i-1)*3+2) = x_Imp((irep-1)*3+2)
	      x_lmp((i-1)*3+3) = x_Imp((irep-1)*3+3)
	      v_lmp((i-1)*3+1) = v_Imp((irep-1)*3+1)
	      v_lmp((i-1)*3+2) = v_Imp((irep-1)*3+2)
	      v_lmp((i-1)*3+3) = v_Imp((irep-1)*3+3)
	  enddo

	  ! Scatter atoms back to Lammps:
	  if(EAMdebug.and.(n_Imp.gt.0).and.explicitOxy)then
            call lammps_scatter_atoms (qeq, 'type', type_lmp)
            call lammps_scatter_atoms (qeq, 'x', x_lmp)
	  endif
          call lammps_scatter_atoms (lmp, 'type', type_lmp)
          call lammps_scatter_atoms (lmp, 'x', x_lmp)
          call lammps_scatter_atoms (lmp, 'v', v_lmp)
	  if(EAMdebug.and.(n_Imp.gt.0).and.explicitOxy) 
     &			deallocate(q_lmp)
	  !deallocate(type_lmp)
	  !deallocate(x_lmp)
	  !deallocate(v_lmp)

	  if(EAMdebug.and.(n_Imp.gt.0).and.explicitOxy)then
	    call lammps_command(qeq,'run 1')
            call lammps_gather_atoms (qeq, 'q', 1, q_lmp)
            do i = 1, n_atoms_lmp
              if(q_lmp(i).gt.charge_max)then
                q_lmp(i) = charge_max
              elseif(q_lmp(i).lt.charge_min)then
                q_lmp(i) = charge_min
              endif
            enddo
	    q_lmp(:) = q_lmp(:) * scaleChg
	    if(n_Imp.gt.0)then
              call lammps_scatter_atoms (lmp, 'q', q_lmp)
	    endif
	  endif

	  call lammps_command(lmp,'run 0')

          n_atoms_lmp = lammps_get_natoms(lmp)
!          if(rank.eq.0)write(*,*)'checking n_atoms: ',n_atoms_lmp

        endif


!!====================!! ENTER MAIN LOOP !!====================!!
        FEMSteps=FEMStepMin
	iStep = 0
	keeprunning = .true. 

	! Start with Steepest Decent Scheme:
	str001=''
	write(str001,*)'min_style sd'
	call lammps_command(lmp,str001)

        allocate(maxforcelist(Nsteps+1))
	maxforcelist(1:Nsteps+1)=0.d0

	do while (keeprunning)  !do iStep = 0, Nsteps

	  ! At the very first step, aveDispl = atomDispl for all nodes
          if (nnsteps .eq. 0) then
	    aveDispl(1:ndf, 1:numnp) = atomDispl(1:ndf, 1:numnp)
	  endif

	  ! Find Average positions of free and interface atoms
	  if (mod(FEMStepCounter, FEMSteps) .eq. 0) then	
		solveFEM = .true.
		restartAveraging = .true.
	  else
		solveFEM = .false.
	  endif

          !! Get Energies and Forces on MD atoms
	  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	  if(rank.eq.0) write(*,*)'Solve FEM System..'

          call getEnergiesAndForces (atomID, atomCoord, ix, f,
     &       atomDispl, aveDispl, atomForce,  
     &       systemEnergy, MoveAtoms, MoveDisl, FullField, 
     &       solveFEM, strainE0)

	  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	  if(rank.eq.0) write(*,*)'R0 - Done w/ FEM System. ',iStep
!	  if(rank.eq.1) write(*,*)'R1 - Done w/ FEM System. ',iStep
          call lammps_gather_atoms (lmp, 'x', 3, x_lmp)
          call lammps_gather_atoms (lmp, 'type', 1, type_lmp)
	  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	  if(rank.eq.0) write(*,*)'Adjusting displacements..'

          !! Adjust Atom Types and Atom Displacements for Lammps:
	  call adjust_displ(n_atoms_lmp,
     &		atomCoord,atomDispl,type_lmp,x_lmp,
     &		n_Imp,x_Imp,o2_conc,iStep,dis_12,Atom_top,Atom_bot)

	  call lammps_scatter_atoms(lmp, 'x', x_lmp)
          call lammps_scatter_atoms(lmp, 'type', type_lmp)
	  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	  if(rank.eq.0) write(*,*)'Done adjusting displacements..'
	  deallocate(x_lmp)
	  deallocate(type_lmp)

            if(EAMdebug.and.(n_Imp.gt.0).and.explicitOxy)then
	     if(quickminOn)then
	       loopsPerQeq = max(2/stepsPerFE,1)
	     elseif(hftnOn)then
	       loopsPerQeq = max(2/stepsPerFE,1)
	     else
	       loopsPerQeq = max(2/stepsPerFE,1)
	     endif
	     if(mdTemp.gt.(2.d0))then
	       loopsPerQeq = max(4/stepsPerFE,1)
	     endif
	     if(MOD(iStep,loopsPerQeq).eq.0)then
	     if(scaleChg.gt.(0.001))then
              call lammps_gather_atoms (lmp, 'x', 3, x_lmp)
              call lammps_scatter_atoms (qeq, 'x', x_lmp)
              call lammps_command(qeq,'run 1')
              call lammps_gather_atoms (qeq, 'q', 1, q_lmp)
              do i = 1, n_atoms_lmp
                if(q_lmp(i).gt.charge_max)then
                  q_lmp(i) = charge_max
                elseif(q_lmp(i).lt.charge_min)then
                  q_lmp(i) = charge_min
                endif
              enddo
              q_lmp(:) = q_lmp(:) * scaleChg
              call lammps_scatter_atoms (lmp, 'q', q_lmp)
	     endif
	     endif
            endif

	    ! RUN LAMMPS:
            if((mdTemp.le.2.d0).and.(EAMdebug))then

	        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	        if(rank.eq.0) write(*,*)'Run Lammps..'

	        if(quickminOn)then
	          j = stepsPerFE * 10
	        elseif(hftnOn)then
	          j = stepsPerFE * 2
	        else
	          j = stepsPerFE
	        endif
	        str001=''
	        write(str001,*)'minimize 0.0 0.0 ',j,' 1000'
	        CALL lammps_command(lmp,str001)

	        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	        if(rank.eq.0) write(*,*)'Done running Lammps..'

	    else

	        str001=''
	        write(str001,*) 'run ',stepsPerFE
	        CALL lammps_command(lmp,str001)

            endif

	  ! Gather Lammps Results:
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          call lammps_gather_atoms (lmp, 'type', 1, type_lmp)
          call lammps_gather_atoms (lmp, 'x', 3, x_lmp)
          call lammps_gather_atoms (lmp, 'v', 3, v_lmp)
          call lammps_gather_atoms (lmp, 'f', 3, f_lmp)
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	  ! Update CADD arrays after LAMMPS run:
	  call UpdateCadd(atomDispl,x_lmp,atomForce,
     &	    f_lmp,atomCoord,ForceMax,ForceMaxP,EAMdebug,type_lmp)

	  if(rank.eq.0)then
	   !write(*,*) ''
	   write(*,*) '* * * * iStep            = ',iStep
	   write(*,*) '* * * * Maximum Force    = ',ForceMax
	   !write(*,*) '* * * * Maximum Force Al = ',ForceMaxP
	   !write(*,*) ''
	  endif

	  maxforcelist(iStep+1)=max(ForceMax,ForceMaxP)
	  if(iStep.gt.1)then
	    if(iStep.eq.2) switchok = .true.
	    fgrad = abs(maxforcelist(iStep+1)
     &		-maxforcelist(iStep-1))
	  endif

	  ! Update Impurity Arrays:
	  j=0
	  do i=1,n_atoms_lmp
	    if((type_lmp(i).ne.1).and.(type_lmp(i).lt.4))then
	      j=j+1
	      x_Imp((j-1)*3+1) = x_lmp((i-1)*3+1)
	      x_Imp((j-1)*3+2) = x_lmp((i-1)*3+2)
	      x_Imp((j-1)*3+3) = x_lmp((i-1)*3+3)
	      v_Imp((j-1)*3+1) = v_lmp((i-1)*3+1)
	      v_Imp((j-1)*3+2) = v_lmp((i-1)*3+2)
	      v_Imp((j-1)*3+3) = v_lmp((i-1)*3+3)
	      type_Imp(j) = type_lmp(i)
	      if(x_Imp((j-1)*3+1).lt.grid_xmin)then
		v_Imp((j-1)*3+1)=vstart
		v_lmp((i-1)*3+1)=vstart
	      endif
	    endif
	  enddo

	  ! Scatter Results:
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          call lammps_scatter_atoms (lmp, 'type', type_lmp)
          call lammps_scatter_atoms (lmp, 'x', x_lmp)
          call lammps_scatter_atoms (lmp, 'v', v_lmp)
          call lammps_scatter_atoms (lmp, 'f', f_lmp)
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	  if ( ( (fgrad.lt.0.000001).and.switchok )
     &         .or.( (iStep.gt.25).and.switchok ) ) then

	    switchok = .false.
	    quickminOn = .true.
	    str001=''
	    write(str001,*)'min_style quickmin'
	    call lammps_command(lmp,str001)
	    str001=''
	    write(str001,*)'min_modify dmax 0.1'
	    call lammps_command(lmp,str001)

	  elseif ((fgrad.lt.0.000001).and.quickminOn) then

	    quickminOn = .false.
	    hftnOn = .true.
	    str001=''
	    write(str001,*)'min_style hftn'
	    call lammps_command(lmp,str001)

	  endif

	  if (ForceMaxP.le.Ftol) then

	    !if (rank.eq.0) write(*,*) ' Forces are converged!'
	    keeprunning = .false.

	  endif

	  ! Restart Averaging:
 	  if (reStartAveraging .eq. .true.) then
	      do iAtom = 1, numnp	
		if ((isRelaxed(iAtom).eq.indexAtom)
     &          .or.(isRelaxed(iAtom).eq.indexInterface)) then
		  aveDispl(1:ndf, iAtom) = 0.d0
		endif	
	      enddo
	      FEMSteps = FEMStepMin
	      if (FEMSteps.gt.Nsteps)then
	        FEMSteps = Nsteps+1
	      endif
	      FEMStepCounter = 0
	      reStartAveraging = .false.
          endif

	  ! Calculate average displacements for interface and free atoms
     	  do 30 iAtom = 1, numnp
	     if (isRelaxed(iAtom) .eq. indexAtom
     &	      .or. isRelaxed(iAtom) .eq. indexInterface) then	
     	         aveDispl(1:ndf, iAtom) =  
     &		 aveDispl(1:ndf, iAtom) + 
     &		 atomDispl(1:ndf, iAtom)/FEMSteps
	     endif	
30	  continue

          FEMStepCounter = FEMStepCounter + 1

          ! Check for emitted Dislocations:
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          !if(rank.eq.0) write(*,*)'R0-Start dislcheck stuff. ',iStep
          !if(rank.eq.1) write(*,*)'R1-Start dislcheck stuff. ',iStep
	  dchk=dislcheck(CheckSlip,LostSlip,AddedSlip,MoveDisl,ix,
     &      atomCoord,atomDispl,itx,IsRelaxed,numnp,ndf,
     &      nxdm,numel,nen1,newmesh,nnsteps)
          !if(rank.eq.0) write(*,*)'R0-Done dislcheck stuff. ',iStep
          !if(rank.eq.1) write(*,*)'R1-Done dislcheck stuff. ',iStep
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          if (dchk) then
             ndis_checked=ndis_checked+1
!             if(rank.eq.0) write(*,*) nnsteps,'ps rank=',rank,
!     &                   ' disnumchecked ',ndis_checked
          endif
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          !if(rank.eq.0) write(*,*)'R0-Done w/ dowhile loop. ',iStep
          !if(rank.eq.1) write(*,*)'R1-Done w/ dowhile loop. ',iStep

          iStep = iStep + 1
	  nnsteps=nnsteps+1
	  if (iStep.gt.Nsteps) keeprunning = .false.

        enddo
!!====================!! END MAIN LOOP !!====================!!

        !! Build nearest neighbor lists...
        !! I have removed this because I am using the CV grid 
        !! to detect "surface" atoms instead.
        !! call build_nnlists(n_Al,x_lmp) 

	! Write IMP File:
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(rank.eq.0)then
	open(unit=215, file='IMP', status='replace')  
	  write(215,*) n_O
	  write(215,*) n_H
	  if(n_Imp.gt.0)then
            do i = 1, n_Imp
	      write(215,'(I6,1X,F16.8,1X,F16.8,1X,F16.8,
     &	    1X,F16.8,1X,F16.8,1X,F16.8)') type_Imp(i), 
     &	    x_Imp((i-1)*3+1),x_Imp((i-1)*3+2),x_Imp((i-1)*3+3),
     &	    v_Imp((i-1)*3+1),v_Imp((i-1)*3+2),v_Imp((i-1)*3+3)
	    enddo
	  endif
	close(215)
	endif

	! Free Impurity Arrays:
        if(ALLOCATED(x_Imp)) deallocate(x_Imp)
        if(ALLOCATED(v_Imp)) deallocate(v_Imp)
        if(ALLOCATED(type_Imp)) deallocate(type_Imp)

	! Write VEL File:
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(rank.eq.0)then
	!write(*,*) 'Writing VEL file..'
	open(unit=214, file='VEL', status='replace')  
         i=0
         do irep = 1, numnp
            NeedList=IsRelaxed(irep).ne.0
            if(NeedList)then
                i=i+1      
                write(214,'(F16.8,1X,F16.8,1X,F16.8)')
     &		v_lmp((i-1)*3+1),
     &		v_lmp((i-1)*3+2), v_lmp((i-1)*3+3)
            endif
	 enddo
	close(214)
	!write(*,*) 'VEL file Written.'
	endif
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	! Free Remaining Allocated Arrays:
        if(ALLOCATED(dis_12)) deallocate(dis_12)
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'deallocating dis_12'
        if(ALLOCATED(x_lmp)) deallocate(x_lmp)
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'de-aloc 1'
        if(ALLOCATED(q_lmp)) deallocate(q_lmp)
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'de-aloc 2'
        if(ALLOCATED(v_lmp)) deallocate(v_lmp)
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'de-aloc 3'
        if(ALLOCATED(type_lmp)) deallocate(type_lmp)
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'de-aloc 4'
        if(ALLOCATED(f_lmp)) deallocate(f_lmp)
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'de-aloc 5'
        if(ALLOCATED(CV_cell_free)) deallocate(CV_cell_free)
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'de-aloc 6'
        if(ALLOCATED(CV_list)) deallocate(CV_list)
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'de-aloc 7'
	if(ALLOCATED(maxforcelist)) deallocate(maxforcelist)

        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'Remaining Arrays Deallocated.'

	! CHECK DEALLOCATION.. Just in case...
        if(ALLOCATED(type_Imp)) deallocate(type_Imp)
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'de-aloc 8'
        if(ALLOCATED(type_Imp_temp)) deallocate(type_Imp_temp)
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'de-aloc 9'
        if(ALLOCATED(x_Imp)) deallocate(x_Imp)
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'de-aloc 10'
        if(ALLOCATED(v_Imp)) deallocate(v_Imp)
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'de-aloc 11'
        if(ALLOCATED(x_Imp_temp)) deallocate(x_Imp_temp)
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'de-aloc 12'
        if(ALLOCATED(v_Imp_temp)) deallocate(v_Imp_temp)

        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'Other arrays checked.'

	! change name of dump_CV file
	if(writefiles)then
	  if(MDSteps.lt.10)then
	    str001='mv dump_CV dump_CV.000'
	  elseif(MDSteps.lt.100)then
	    str001='mv dump_CV dump_CV.00'
	  elseif(MDSteps.lt.1000)then
	    str001='mv dump_CV dump_CV.0'
	  else
	    str001='mv dump_CV dump_CV.'
	  endif
	  str002=''
	  write(str002,*) MDSteps
	  str003=''
 	  str003=trim(adjustl(str001))//trim(adjustl(str002))
	  if(rank.eq.0) call system(str003)

	  if(MDSteps.lt.10)then
	    str001='mv dump_type dump_type.000'
	  elseif(MDSteps.lt.100)then
	    str001='mv dump_type dump_type.00'
	  elseif(MDSteps.lt.1000)then
	    str001='mv dump_type dump_type.0'
	  else
	    str001='mv dump_type dump_type.'
	  endif
	  str002=''
	  write(str002,*) MDSteps
	  str003=''
 	  str003=trim(adjustl(str001))//trim(adjustl(str002))
	  if(rank.eq.0) call system(str003)
	endif

	if(rank.eq.0)then
	  write(*,*) 'MDSteps, x0crack --> ',
     &		MDSteps, x_cracktip
	  write(*,*) 'o2_conc_max: ',o2_conc_max
	  open(unit=1115,file='da_ctod.out',
     &			access='append',status='unknown')
	  write(1115,*) MDSteps,x_cracktip,ctod
	  close(1115)
	endif
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	MDSteps = MDSteps + 1
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	return
	end

