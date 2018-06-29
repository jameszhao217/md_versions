!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
!FFFF!                                                          !FFFF!
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
	subroutine UpdateCadd(atomDispl,x_lmp,atomForce,
     &	 f_lmp,atomCoord,ForceMax,ForceMaxP,EAMdebug,type_lmp)

	use mod_global
	use mod_grain
	implicit none
        double precision, allocatable, INTENT(INOUT) :: x_lmp(:)
        double precision, allocatable, INTENT(INOUT) :: f_lmp(:)
        integer, allocatable, INTENT(INOUT) :: type_lmp(:)
	integer iStep,maxtype,maxtypeP
	double precision atomDispl(ndf,*), atomForce(ndf,*), 
     &	atomCoord(ndf,*)
	double precision displ,force,maxDispl,ForceMax,ForceMaxP
	logical debug,EAMdebug
	common/debugger/debug
	
C--	Local variables	
	integer iAtom, j, nAtom
	double precision  rmsForce

	nAtom = 0
	ForceMax = 0.d0
	ForceMaxP = 0.d0

	do iAtom = 1, numnp
		
!-	  Skip near continuum nodes and pad atoms
	  if (isRelaxed(iAtom) .ne. indexContinuum) then

           nAtom = nAtom + 1
	   if (isRelaxed(iAtom) .ne. indexPad) then

	    rmsForce=0.d0
!-	    Repeat for each dimension
	    do j =1, ndf
	      if(EAMdebug)then
	        force = f_lmp((nAtom-1)*3+j)
	      else
	        force = f_lmp((nAtom-1)*3+j)*0.04336
	      endif
              displ = x_lmp((nAtom-1)*3+j)-
     &		atomCoord(j,iAtom)
	      if(j.lt.3) then
		atomDispl(j,iAtom) = displ
	      else
		do while(displ.gt.(0.75*z_length))
		  displ=displ-z_length
		enddo
		do while(displ.lt.(-0.75*z_length))
		  displ=displ+z_length
		enddo
		atomDispl(j,iAtom) = displ
	      endif

	      atomForce(j,iAtom) = force
	      rmsForce = rmsForce + atomForce(j, iAtom)**2
	    enddo
	    rmsForce=sqrt(rmsForce)

	    if(rmsForce.gt.ForceMax) then 
	      ForceMax=rmsForce
	      maxtype = type_lmp(nAtom)
	    endif

	    if((rmsForce.gt.ForceMaxP)
     &	    .and.(type_lmp(nAtom).ne.4)) then 
	      ForceMaxP=rmsForce
	      maxtypeP = type_lmp(nAtom)
	    endif

	  endif

	 endif

	enddo

	!write(*,*) ' ForceMax  is type: ',maxtype
	!write(*,*) ' ForceMaxP is type: ',maxtypeP

        return
	end

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
!FFFF!                                                          !FFFF!
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
	subroutine ReadMDinp(timestep,FEMStepMin,Nsteps,
     &	o2_conc,h2o_conc,stepsPerFE,dumpcg,EAMdebug,
     &	mdTemp,Ftol,K_max)
	use mod_global
	implicit none
	
	integer FEMStepMin,Nsteps,stepsPerFE,dumpcg
	logical EAMdebug
	double precision timestep,Ftol,
     &	o2_conc,h2o_conc,mdTemp,K_max

!	Read Data	
	if(rank.eq.0) print *, 'reading imp data'
	if(rank.eq.0) then		
	 open(unit=200, file='md.inp', status='old')
	 read(200,*) mdTemp
         read(200,*) timestep
	 read(200,*) FEMStepMin
	 read(200,*) Nsteps 
	 read(200,*) o2_conc
	 read(200,*) h2o_conc
	 read(200,*) stepsPerFE
	 read(200,*) dumpcg 
	 if(rank.eq.0) print *, 'read up to Adhes imp data'
	 read(200,*) EAMdebug, Monolayer, monoThick, Adhes
	 if(rank.eq.0) print *, 'read Adhes '
	 read(200,*) Ftol ! Should be in [eV/Ang] units
	 read(200,*) Sbarrier
	 read(200,*) scaleChg
	 read(200,*) K_max, fstart

	 close(200)

	 if(rank.eq.0) print *, 'read imp data'
	 if(EAMdebug) timestep=timestep*0.001

	 if(EAMdebug.and.(h2o_conc.gt.(0.000001)))then
	  write(*,*) 'Cannot use EAMdebug with H2O for now...'
	  h2o_conc = 0.d0
	 endif

	 !write(*,*) 'o2_conc ',o2_conc
	 !write(*,*) 'h2o_conc ',h2o_conc
	 !write(*,*) 'stepsPerFE ',stepsPerFE
	 !write(*,*) 'mdTemp ',mdTemp
	 !write(*,*) 'dumpcg ',dumpcg
	 !write(*,*) 'Debugging with EAM? -> ',EAMdebug
         !write(*,*) 'Monolayer: ', Monolayer
         !write(*,*) 'MonoThick: ', monoThick
         !write(*,*) 'Adhes: ', Adhes
         !write(*,*) 'timestep: ', timestep
	 !write(*,*) 'FEMStepMin: ', FEMStepMin
	 !write(*,*) 'Nsteps: ', Nsteps
	 !write(*,*) 'Ftol: ', Ftol 
	 !write(*,*) 'Sbarrier: ', Sbarrier
	 !write(*,*) 'scaleChg: ', scaleChg
	 !write(*,*) 'K_max: ', K_max
	 !write(*,*) 'fstart: ', fstart

	endif
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mdTemp,1,MPI_DOUBLE,0,MPI_COMM_WORLD,
     & ierr)
      CALL MPI_BCAST(timestep,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(FEMStepMin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Nsteps,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(o2_conc,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(h2o_conc,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(stepsPerFE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mdTemp,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dumpcg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(EAMdebug,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Monolayer,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Adhes,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(monoThick,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Ftol,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Sbarrier,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(scaleChg,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(K_max,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fstart,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

        return
	end

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
!FFFF!                                                          !FFFF!
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
	subroutine ReadIMP(velFile,n_O,n_H,n_Imp,
     &	x_Imp,v_Imp,type_Imp,redistribute,adiff_T,adiff_n)
	use mod_global
	implicit none
	
	logical redistribute
	integer i,velFile,n_O,n_H,n_Imp,adiff_n
	double precision adiff_T
        double precision, allocatable, INTENT(INOUT) :: x_Imp(:)
        double precision, allocatable, INTENT(INOUT) :: v_Imp(:)
        integer, allocatable, INTENT(INOUT) :: type_Imp(:)

        velFile=0
	if(rank.eq.0) write(*,*) 'firstIMP: ',firstIMP
	if(firstIMP.eq.0)then
	  !redistribute = .false.
          velFile=1
	  if(rank.eq.0) open(unit=215, file='IMP', status='unknown')
	  if(rank.eq.0) read(215,*) n_O
	  if(rank.eq.0) read(215,*) n_H
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(n_O,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(n_H,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	  n_Imp = n_O + n_H
	  if(.not.ALLOCATED(x_Imp)) allocate(x_Imp(n_Imp*3))
	  if(.not.ALLOCATED(v_Imp)) allocate(v_Imp(n_Imp*3))
	  if(.not.ALLOCATED(type_Imp)) allocate(type_Imp(n_Imp))
	  if(n_Imp.gt.0)then
            do i = 1, n_Imp
	     if(rank.eq.0)read(215,*) type_Imp(i), x_Imp((i-1)*3+1),
     &	     x_Imp((i-1)*3+2), x_Imp((i-1)*3+3),
     &	     v_Imp((i-1)*3+1),v_Imp((i-1)*3+2),v_Imp((i-1)*3+3)
             CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(type_Imp(i),1,MPI_INTEGER,0,
     &		MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(x_Imp((i-1)*3+1),1,MPI_DOUBLE,0,
     &		MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(x_Imp((i-1)*3+2),1,MPI_DOUBLE,0,
     &		MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(x_Imp((i-1)*3+3),1,MPI_DOUBLE,0,
     &		MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(v_Imp((i-1)*3+1),1,MPI_DOUBLE,0,
     &		MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(v_Imp((i-1)*3+2),1,MPI_DOUBLE,0,
     &		MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(v_Imp((i-1)*3+3),1,MPI_DOUBLE,0,
     &		MPI_COMM_WORLD,ierr)
	    enddo
	  endif
	  if(rank.eq.0) close(215)
	else
	  if(rank.eq.0) open(unit=218,file='impurities.init',
     &	  status='unknown')
	  if(rank.eq.0) read(218,*) n_O
	  if(rank.eq.0) read(218,*) n_H
	  if(rank.eq.0) then
	    read(218,*) redistribute,adiff_T,adiff_n
	    if(redistribute)then
	      write(*,*)''
	      write(*,*)'Redistributing initial impurities...'
	      write(*,*)' This is only being done for H atoms!! '
	      write(*,*)''
	    endif
	  endif
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(n_O,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(n_H,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(adiff_T,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(adiff_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(redistribute,1,MPI_LOGICAL,0,
     &		MPI_COMM_WORLD,ierr)
	  n_Imp = n_O + n_H
	!if(rank.eq.0) write(*,*) 'Start Allocating... '
	  if(.not.ALLOCATED(x_Imp)) allocate(x_Imp(n_Imp*3))
	  if(.not.ALLOCATED(v_Imp)) allocate(v_Imp(n_Imp*3))
	  if(.not.ALLOCATED(type_Imp)) allocate(type_Imp(n_Imp))
	!if(rank.eq.0) write(*,*) 'Finish Allocating... '
	  if(n_Imp.gt.0)then
            do i = 1, n_Imp
	     if((i.gt.n_O).and.(redistribute))then
	      type_Imp(i)=3
	      x_Imp((i-1)*3+1)=0.d0
	      x_Imp((i-1)*3+2)=0.d0
	      x_Imp((i-1)*3+3)=0.d0
	      v_Imp((i-1)*3+1)=0.d0
	      v_Imp((i-1)*3+2)=0.d0
	      v_Imp((i-1)*3+3)=0.d0
	     else
	      if(rank.eq.0)read(218,*) type_Imp(i), x_Imp((i-1)*3+1),
     &	      x_Imp((i-1)*3+2), x_Imp((i-1)*3+3),
     &	      v_Imp((i-1)*3+1),v_Imp((i-1)*3+2),v_Imp((i-1)*3+3)
              CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
              CALL MPI_BCAST(type_Imp(i),1,MPI_INTEGER,0,
     &		MPI_COMM_WORLD,ierr)
              CALL MPI_BCAST(x_Imp((i-1)*3+1),1,MPI_DOUBLE,0,
     &		MPI_COMM_WORLD,ierr)
              CALL MPI_BCAST(x_Imp((i-1)*3+2),1,MPI_DOUBLE,0,
     &		MPI_COMM_WORLD,ierr)
              CALL MPI_BCAST(x_Imp((i-1)*3+3),1,MPI_DOUBLE,0,
     &		MPI_COMM_WORLD,ierr)
              CALL MPI_BCAST(v_Imp((i-1)*3+1),1,MPI_DOUBLE,0,
     &		MPI_COMM_WORLD,ierr)
              CALL MPI_BCAST(v_Imp((i-1)*3+2),1,MPI_DOUBLE,0,
     &		MPI_COMM_WORLD,ierr)
              CALL MPI_BCAST(v_Imp((i-1)*3+3),1,MPI_DOUBLE,0,
     &		MPI_COMM_WORLD,ierr)
	     endif
	    enddo
	  endif
	  if(rank.eq.0) close(218)
	  firstIMP=0
	endif
	if(rank.eq.0) print *, 'finishing ReadIMP'
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

        return
	end

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
!FFFF!                                                          !FFFF!
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
	subroutine CV_setup(atomCoord,atomDispl,n_CV,
     &			CV_cell_free,CV_list,dis_12,CV_cell_top,
     &                  CV_cell_bot,Atom_top,Atom_bot)
	use mod_global
	implicit none
	logical Needlist,notfreenomo
	integer irep,i,j,ii,jj,xsite_i,ysite_i,iin,jin,ifn,jfn,n_CV,
     &          dir,n_edge,n_found,flip_top,flip_bot,iii,jjj
        integer nlook,n_Al,MDSteps
        double precision cx,cy,atomCoord(ndf,*),atomDispl(ndf,*),
     &		xsite_d,ysite_d,cxg,cyg,nearsurf,xTip

        integer, allocatable :: sumit(:,:)
        integer, allocatable :: proxim(:,:)
        integer, allocatable, INTENT(INOUT) :: CV_cell_free(:,:)
	integer, allocatable, INTENT(INOUT) :: CV_list(:,:)

!ER- I made these to try to distinguish top and bottom layers of Al
        integer, allocatable :: sumit_surf(:,:)
        integer, allocatable, INTENT(INOUT) :: CV_cell_top(:,:)
        integer, allocatable, INTENT(INOUT) :: CV_cell_bot(:,:)
!        integer, allocatable :: CV_cell_edge(:,:)

        double precision, allocatable, INTENT(INOUT) :: dis_12(:)

!ER- I made these to try to distinguish top and bottom layers of Al        
        double precision, allocatable :: surf_x(:,:)
        integer, allocatable, INTENT(INOUT) :: Atom_top(:)
        integer, allocatable, INTENT(INOUT) :: Atom_bot(:)
	
	! Search all CV grid cells for Al atoms...
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'aloc A'
        if(.not.allocated(CV_list)) allocate(CV_list(n_CV_cells,2))
        
	!CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'aloc B'
        if(.not.allocated(CV_cell_free)) 
     &		allocate(CV_cell_free(nx_grid,ny_grid))

!ER- I made these to try to distinguish top and bottom layers of Al
	if(.not.allocated(CV_cell_top)) 
     &          allocate(CV_cell_top(nx_grid,ny_grid))
        if(.not.allocated(CV_cell_bot)) 
     &          allocate(CV_cell_bot(nx_grid,ny_grid))
	
	print *, 'declared and allocated' 
        
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!if(rank.eq.0)write(*,*) 'aloc C'
	CV_cell_free(:,:) = 1

!ER- I made these to try to distinguish top and bottom layers of Al
           CV_cell_top(:,:)=0
           CV_cell_top(:,:)=0
!	print *, 'In CV_setup, here are all the input variables:'


!MZ- modified code to produce less contents in file 'output'
	n_Al = 0
	
	do irep = 1, numnp    

	    NeedList=IsRelaxed(irep).ne.0

	    if(NeedList) then
	        
		n_Al = n_Al + 1  ! count number of Al atoms

	     	cx = atomCoord(1,irep) + atomDispl(1,irep) ! find locations
	     	cy = atomCoord(2,irep) + atomDispl(2,irep)

		!Check if Al atom is in a CV_cell at all...
	     	if(cx.le.grid_xmin) goto 2171
		if(cx.gt.grid_xmax-grid_dx) goto 2171
		if(cy.le.grid_ymin) goto 2171
		if(cy.gt.grid_ymax-grid_dx) goto 2171
	
		! Determine which cell we are in:
	     	xsite_d = (cx-grid_xmin)/(nx_grid*grid_dx)
	     	xsite_i = CEILING(xsite_d*nx_grid)
	     	ysite_d = (cy-grid_ymin)/(ny_grid*grid_dx)
	     	ysite_i = CEILING(ysite_d*ny_grid)

! 	     	if(rank.eq.0)write(*,*)'xsite_i,ysite_i=',xsite_i,ysite_i

	     	if(xsite_i.lt.1)then
!		    if(rank.eq.0)write(*,*)'xsite_i,ysite_i=',xsite_i,ysite_i
	            go to 2171
	     	elseif(ysite_i.lt.1)then
!                   if(rank.eq.0)write(*,*)'xsite_i,ysite_i=',xsite_i,ysite_i
	            go to 2171
             	elseif(xsite_i.gt.nx_grid)then
!                   if(rank.eq.0)write(*,*)'xsite_i,ysite_i=',xsite_i,ysite_i
	       	    go to 2171
             	elseif(ysite_i.gt.ny_grid)then
!                   if(rank.eq.0)write(*,*)'xsite_i,ysite_i=',xsite_i,ysite_i
	            go to 2171
             	endif

	     ! Take occupied CV_cell out of play:
	     	CV_cell_free(xsite_i,ysite_i) = 0
!	     	print *, 'cell ',xsite_i,' ',ysite_i,' is occupied.'	
            
	    endif

2171	  continue
        
	enddo


        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	if(rank.eq.0)write(*,*) 'CV_cell_free Made.'
	if(rank.eq.0)write(*,*) 'CV_cell_free Made.'

	! Check Round 2:
!        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	if(rank.eq.0)write(*,*) 'aloc D'
        if(.not.allocated(sumit)) allocate(sumit(nx_grid,ny_grid))
!        if(.not.allocated(sumit)) allocate(sumit_surf(nx_grid,ny_grid))
!        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	if(rank.eq.0)write(*,*) 'aloc E'
        if(.not.allocated(proxim)) allocate(proxim(nx_grid,ny_grid))
!        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	if(rank.eq.0)write(*,*) 'aloc F'
	sumit(1:nx_grid,1:ny_grid)=0
!	sumit_surf(1:nx_grid,1:ny_grid)=0
	!proxim(1:nx_grid,1:ny_grid)=0
	proxim(1:nx_grid,1:ny_grid)=1
	nearsurf = 2.3
	nlook = int( monoThick / grid_dx )
        
	if(rank.eq.0)print *, 'nlook is ',nlook	
	if(rank.eq.0)print *, 'checking if cells are free'	
	do i = 1, nx_grid
	  do j = 1, ny_grid
	    if(CV_cell_free(i,j).eq.1)then
	      ! check that all surrounding cells are also free
	      do ii = (i-1), (i+1)
		do jj = (j-1), (j+1)
                  if((ii.ge.1).and.(ii.le.nx_grid))then
                    if((jj.ge.1).and.(jj.le.ny_grid))then
		      sumit(i,j)=sumit(i,j)+CV_cell_free(ii,jj)
                    endif
                  endif
	        enddo
	      enddo
!	    elseif(CV_cell_free(i,j).eq.0)then
!              ! check that all surrounding cells are also free
!	      do ii = (i-1), (i+1)
!		do jj = (j-1), (j+1)
!                  if((ii.ge.1).and.(ii.le.nx_grid))then
!                    if((jj.ge.1).and.(jj.le.ny_grid))then
!		    sumit_surf(i,j)=sumit_surf(i,j)+CV_cell_free(ii,jj)
!                    endif
!                  endif
!	        enddo
!	      enddo

	    endif
	  enddo
	enddo

        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	if(rank.eq.0)write(*,*) 'sumit Made.'

	! Make CV_list:
        CV_list(:,:) = 0
	n_CV = 0
	do i = 1, nx_grid
	  do j = 1, ny_grid
	    if(CV_cell_free(i,j).eq.1)then
	     if((sumit(i,j).gt.sumReq).and.(proxim(i,j).eq.1))then
	       n_CV = n_CV + 1
	       CV_list(n_CV,1) = i
	       CV_list(n_CV,2) = j
	     else
	       CV_cell_free(i,j)=0
	     endif
	    endif
	  enddo
	enddo
	
!	! Make surface list using my sumit surf:
!        CV_cell_edge(:,:) = 0
!	xTip=0
!	n_edge=0
!	do i = 1, nx_grid
!	  do j = 1, ny_grid
!	    if(CV_cell_free(i,j).eq.0)then
!		!set the vector to 1 if its on the edge
!	     if((sumit_surf(i,j).lt.sumReq).and.
!     &		(sumit_surf(i,j).gt.0))then
!	       CV_cell_edge(i,j)=1
!	       n_edge=n_edge+1
!		!now see if its at the tip
!		if(i.gt.xTip)then
!			xTip=i
!		endif
!	     endif
!	    endif
!	  enddo
!	enddo

	! find the start of the top crack face 
!	n_found=0
!	do i = 1, nx_grid
!	  do j = 1, ny_grid
!	    if(CV_cell_free(i,j).eq.0)then
!		!set the vector to 1 if its on the edge
!	     if((sumit(i,j).lt.sumReq).and.(sumit(i,j).gt.0))then
!	       CV_cell_edge(i,j)=1
!	     endif
!	    endif
!	  enddo
!	enddo


        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	if(rank.eq.0)write(*,*) 'CV_cell_free Modified #1.'

        allocate(dis_12(n_Al))
	if(rank.eq.0) write(*,*) 'Allocating dis_12...'
	dis_12(1:n_Al)=0

	if(rank.eq.0) write(*,*) 'Allocated dis_12...'
!ER- I made these to try to distinguish top and bottom layers of Al
        if(.not.allocated(Atom_top)) allocate(Atom_top(n_Al))
        if(.not.allocated(Atom_bot)) allocate(Atom_bot(n_Al))
	Atom_top(1:n_Al)=0
	Atom_bot(1:n_Al)=0

	if(rank.eq.0) write(*,*) 'Allocated Atom vecs...'

	do irep = 1, n_CV
	   i = CV_list(irep,1)
	   j = CV_list(irep,2)
	   
   
	   do ii = i-nlook,i+nlook
             do jj = j-nlook,j+nlook
               if((ii.ge.1).and.(ii.le.nx_grid))then
                 if((jj.ge.1).and.(jj.le.ny_grid))then
                   if(CV_cell_free(ii,jj).eq.0)then
		     CV_cell_free(ii,jj) = -1 

!                     if(rank.eq.0) print*, ii,' ',jj,'was named edge'

!ER- I made these to try to distinguish top and bottom layers of Al
!			!if the Al you saw was above you-	     
!			if(rank.eq.0) print*, 'jj is: ',jj,', j: ',j
!			if(jj.gt.j) then
!				CV_cell_top(i,j)=1
!			else if(jj.lt.j) then
!				CV_cell_bot(i,j)=1
!			endif


		   endif 
                 endif
               endif
             enddo
           enddo
	enddo

!ER- for now, a very rough way of assigning CVcell top & bot
	do iii=1, nx_grid
	  flip_top=0
	  flip_bot=0

	  do jjj=2,(ny_grid-1)
	  if(flip_bot.eq.0)then
	    if((CV_cell_free(iii,jjj).eq.-1).and.
     &          (CV_cell_free(iii,jjj-1).eq.0))then
		flip_bot=1
	    endif
	  endif
	 
          if(flip_top.eq.0)then
	   if((CV_cell_free(iii,jjj).eq.-1).and.
     &          (CV_cell_free(iii,jjj-1).eq.1))then
		flip_top=1
	    endif
	  endif
	
          if((flip_bot.eq.1).and.(CV_cell_free(iii,jjj).eq.-1).and.
     &          (flip_top.eq.0))then
!		if(rank.eq.0)print *, iii,' ',jjj,' is bot'
		CV_cell_bot(iii,jjj)=1
	  elseif((flip_top.eq.1).and.(CV_cell_free(iii,jjj).eq.-1))then
!		if(rank.eq.0)print *, iii,' ',jjj,' is top'
		CV_cell_top(iii,jjj)=1
	  endif
	  enddo

	!if we make it through a whole column and never found the top
	!then go through and un-flip everythink we said was bot
	if((flip_bot.eq.1).and.(flip_top.eq.0))then
	  do jjj=2,(ny_grid-1)
		if(CV_cell_bot(iii,jjj).eq.1)then
			CV_cell_bot(iii,jjj)=0
		endif
	  enddo
	endif

	enddo
	

!this is going to store the horizontal coord of cv on crack surfaceOB
	allocate(surf_x(2,n_Al))
	surf_x(:,:)=0

        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	if(rank.eq.0)write(*,*) 'CV_cell_free Modified #2.'

	i = 0
        do irep = 1, numnp
           NeedList=IsRelaxed(irep).ne.0
           if(NeedList)then
	     i = i+1

             cx = atomCoord(1,irep) + atomDispl(1,irep)
             cy = atomCoord(2,irep) + atomDispl(2,irep)

             ! Check if Al atom is in a CV_cell at all...
             if(cx.le.grid_xmin) goto 2172
             if(cx.gt.grid_xmax-grid_dx) goto 2172
             if(cy.le.grid_ymin) goto 2172
             if(cy.gt.grid_ymax-grid_dx) goto 2172

             ! Determine which cell we are in:
             xsite_d = (cx-grid_xmin)/(nx_grid*grid_dx)
             xsite_i = CEILING(xsite_d*nx_grid)
             ysite_d = (cy-grid_ymin)/(ny_grid*grid_dx)
             ysite_i = CEILING(ysite_d*ny_grid)

	     if(xsite_i.lt.1)then
!		if(rank.eq.0)write(*,*)'xsite_i,ysite_i='
!     &				,xsite_i,ysite_i
	        go to 2172
	     elseif(ysite_i.lt.1)then
!                if(rank.eq.0)write(*,*)'xsite_i,ysite_i='
!     &                          ,xsite_i,ysite_i
	        go to 2172
             elseif(xsite_i.gt.nx_grid)then
!                if(rank.eq.0)write(*,*)'xsite_i,ysite_i='
!     &                          ,xsite_i,ysite_i
	        go to 2172
             elseif(ysite_i.gt.ny_grid)then
!                if(rank.eq.0)write(*,*)'xsite_i,ysite_i='
!     &                          ,xsite_i,ysite_i
	        go to 2172
             endif

             ! Set dis_12 to well below Sbarrier if not near surface:
             if(CV_cell_free(xsite_i,ysite_i).eq.-1)then
		dis_12(i) = 2.0*Sbarrier 
		
!Here we set surf_x
		surf_x(1,i)=1
		surf_x(2,i)=xsite_d
		
		if(CV_cell_top(xsite_i,ysite_i).eq.1)then
!                if(rank.eq.0)write(*,*) i,' is top'
		   Atom_top(i)=1
		endif

		if(CV_cell_bot(xsite_i,ysite_i).eq.1)then
!               if(rank.eq.0)write(*,*) i,' is bot'
		   Atom_bot(i)=1
		endif


	     else
	        dis_12(i) = 0.d0
	     endif

           endif
2172      continue
        enddo

!Now we're going to use surf_x to find the tip of the crack
	!initialize xTip to a very negative # (this could mess up)
!	do i=1,n_Al
!		if(surf_x(1,i).eq.1)then !it is on surface
!			if(surf_x(2,i).gt.xTip)then
!				xTip=surf_x(2,i)
!			endif
!		endif 
!	enddo

!Now name the surface atoms as top or bot
!	do i=1,n_Al
!		if(surf_x(1,i).eq.1)then !it is on surface
!			if(surf_x(2,i).gt.xTip)then
!				xTip=surf_x(2,i)
!			endif
!		endif 
!	enddo


        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	if(rank.eq.0)write(*,*) 'dis_12 Modified.'

!        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	if(rank.eq.0)write(*,*) 'de-aloc G'
        if(allocated(sumit)) deallocate(sumit)
!        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	if(rank.eq.0)write(*,*) 'de-aloc H'
        if(allocated(proxim)) deallocate(proxim)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	if(rank.eq.0)write(*,*) 'de-aloc I'

	! WRITE FILE TO CHECK THAT CV GRID IS SET UP CORRECTLY:
	if(rank.eq.0) then
	  open(unit=2121, file='dump_CV', status='replace')
	  do i = 1, nx_grid
		    do j = 1, ny_grid
	      cx = grid_xmin + (i-1) * grid_dx + grid_dx/2.0
	      cy = grid_ymin + (j-1) * grid_dx + grid_dx/2.0
              write(2121,*) cx, cy, CV_cell_free(i,j)
	    enddo
          enddo
	  write(2121,*) 'ZONE'
	  do irep = 1, numnp
            NeedList=IsRelaxed(irep).ne.0
            if(NeedList)then
	     cx = atomCoord(1,irep) + atomDispl(1,irep)
	     cy = atomCoord(2,irep) + atomDispl(2,irep)
              write(2121,*) cx, cy, 0.0
            endif
          enddo
	  close(2121)
	endif
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

        return
	end

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
!FFFF!                                                          !FFFF!
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
	subroutine plan_IMP(n_remove,n_Imp,
     &		CV_cell_free,x_Imp,type_Imp,
     &		n_o2_need,n_h2o_need,o2_conc,h2o_conc,n_CV,
     &		EAMdebug,nnsteps,atomCoord,atomDispl,n_surf,
     &		x_cracktip,dis_12)
	use mod_global
	implicit none
	
	integer n_remove,n_Imp,xsite_i,ysite_i,
     &		n_o2_need,n_h2o_need,
     &		i,n_CV,nnsteps,j,irep,n_surface
        double precision n_O_exp,n_H_exp,n_H_in,n_O_in,
     &		n_o2_exp,n_h2o_exp,n_O_need,n_H_need,n_surf
        double precision xsite_d,ysite_d,c_vol,x_cracktip,
     &		o2_conc,h2o_conc,randomNum1,randomNum2
        double precision cx,cy,cz,cxO,cyO,czO,dmag,imp_min
        double precision atomCoord(ndf,*),atomDispl(ndf,*)
	logical EAMdebug,NeedList
        integer, allocatable, intent(inout) :: CV_cell_free(:,:)
        double precision, allocatable, intent(inout) :: x_Imp(:)
        integer, allocatable, intent(inout) :: type_Imp(:)
        double precision, allocatable, INTENT(INOUT) :: dis_12(:)

	! Count number of O and H in control vol.
	n_remove=0
	n_H_in=0.0
	n_O_in=0.0
	
	if(nnsteps.eq.1)then
	  imp_min = (grid_xmin)
	else
	  imp_min = (x_cracktip-20.d0)
	endif

	print *, 'in planImp'	
	if(Monolayer)then
	
	print *, 'Monolayer is true'	
	 if((nnsteps.gt.0).and.(explicitOxy))then
	
	print *, 'nnsteps > 0 and explicitOxy is True'	
	print *, 'or is it?'	
	print *, 'nnsteps.gt.0: ', nnsteps.gt.0	
	print *, 'explicityOxy: ', explicitOxy	

	print *, '(nnsteps.gt.0).and.(explicitOxy) is true if were in here'
	  ! Using EAM+LJ... So assum o2_conc is approx monolayer coverage needed
	  ! ML = 1.0 means all surface atoms > grid_xmin are within delta_ml (2.3) of an O atom
	  i=0
	  n_surf = 0.0
          do irep = 1, numnp
            NeedList=IsRelaxed(irep).ne.0
            if(NeedList)then
	 print *, 'NeedList is true'	
	      i=i+1
	      if(dis_12(i).ge.(Sbarrier))then !IS ON SURFACE
	 print *, 'found a surface atom at line 648'	
	        cx = atomDispl(1,irep)+atomCoord(1,irep)
	        cy = atomDispl(2,irep)+atomCoord(2,irep)
	        cz = atomDispl(3,irep)+atomCoord(3,irep)
	        if((cx.gt.imp_min).and.(cx.lt.grid_xmax).and.
     &		   (cy.gt.grid_ymin).and.(cy.lt.grid_ymax))then
	        ! Know we are looking at 'in-play' surface atom now

		  n_surf = n_surf + 1.0

	        endif !((cx.gt.grid_xmin).and.(cx.lt.grid_xmax))
	      endif ! dis_12(i).ge.(0.25)
	    endif ! Needlist
	  enddo !irep = 1, numnp
	
	print *, 'cycling through impurities, n_Imp=',n_Imp	
	  do j = 1, n_Imp
	    if(type_Imp(j).eq.2)then
              cxO = x_Imp((j-1)*3+1)
	      if(cxO.gt.(imp_min))then
		n_O_in = n_O_in + 1.0
	      endif
	    endif
	  enddo

	  endif !nnsteps
	  else !Monolayer

	  print *, 'Monolayer is false'	
          do i = 1, n_Imp
            cx = x_Imp((i-1)*3+1)
	    cy = x_Imp((i-1)*3+2)

	    if(cx.gt.grid_xmax-grid_dx) goto 217
	    if(cy.gt.grid_ymax-grid_dx) goto 217
	    if(cy.lt.grid_ymin) goto 217
	    if(cx.lt.grid_xmin) then
	      !n_remove=n_remove+1
	      goto 217
	    endif

	    ! Determine which cell we are in:
	    xsite_d = (cx-grid_xmin)/(nx_grid*grid_dx)
	    xsite_i = CEILING(xsite_d*nx_grid)
	    ysite_d = (cy-grid_ymin)/(ny_grid*grid_dx)
	    ysite_i = CEILING(ysite_d*ny_grid)
	    if(CV_cell_free(xsite_i,ysite_i).le.0) goto 217

	    if(type_Imp(i).eq.2) n_O_in = n_O_in + 1.0
	    if(type_Imp(i).eq.3) n_H_in = n_H_in + 1.0
217	    continue
	  enddo

	endif !Monolayer
	!endif !nnstep.gt.0

	if(rank.eq.0)then
	 write(*,*) ''
	 write(*,*) n_O_in,' O atoms already in CV'
	 write(*,*) n_H_in,' H atoms already in CV'
	 write(*,*) ''
	 write(*,*) n_remove,' Imp atoms will be REMOVED'
	 write(*,*) ''
	endif

	n_o2_exp = 0.0
	n_O_exp = 0.0
	n_h2o_exp = 0.0
	n_H_exp = 0.0
	!if(nnsteps.gt.0)then
	 !if(EAMdebug)then
	if(Monolayer)then
	 if((nnsteps.gt.0).and.(explicitOxy))then

	   n_O_exp = n_surf
	   n_O_exp = n_O_exp * o2_conc
	   n_o2_exp = n_O_exp * 2.0

	 endif !nnsteps
	else !Monolayer

	   ! Calculate Control Volume:
	   c_vol = REAL(n_CV) * grid_dx * grid_dx
	   c_vol = c_vol * z_length
	   ! Calculate Required # H2O and O2 in CV:
	   n_o2_exp = o2_conc * c_vol
	   n_h2o_exp = h2o_conc * c_vol
	   ! Calculate Required # O and H in CV:
	   n_O_exp = 2*n_o2_exp + 1*n_h2o_exp
	   n_H_exp = 0*n_o2_exp + 2*n_h2o_exp

	endif !Monolayer

	if(rank.eq.0)then
	 write(*,*) ''
	 write(*,*) n_o2_exp,' is n_o2_exp'
	 write(*,*) o2_conc,' is o2_conc'
	 write(*,*) c_vol,' is c_vol'
	 write(*,*) n_O_exp,' O atoms EXPECTED in CV'
	 write(*,*) n_H_exp,' H atoms EXPECTED in CV'
	 write(*,*) ''
	endif

	! Calculate # O and H atoms we need to add:
	n_O_need = n_O_exp - n_O_in
	if(n_O_need.lt.0) n_O_need=0.0
	n_H_need = n_H_exp - n_H_in
	if(n_H_need.lt.0) n_H_need=0.0
	if(rank.eq.0)then
	 write(*,*) ''
	 write(*,*) n_O_exp,' n_O_exp'
	 write(*,*) n_O_in,' n_O_in'
	 write(*,*) n_O_need,' O MORE atoms NEEDED in CV'
	 write(*,*) n_H_need,' H MORE atoms NEEDED in CV'
	 write(*,*) ''
        endif
	! Determine # of H2O molecules need to add:
	n_h2o_need = int(ceiling(n_H_need/2.0))
	if(n_h2o_need.lt.0) n_h2o_need=0

	if(Monolayer)then
	  n_o2_need = int(n_O_need) ! For EAM+LJ System n_o2_need is actually n_O_need
	else
	  ! Determine # of O2 molecules need to add:
	  n_o2_need = int(n_O_need-n_h2o_need)
	  n_o2_need = n_o2_need/2
	endif

	if(n_o2_need.lt.0) n_o2_need=0
!	if(rank.eq.0)then
!	 write(*,*) ''
!	 write(*,*) n_o2_need,' O2 MORE molecules NEEDED in CV'
!	 write(*,*) n_h2o_need,' H2O MORE molecules NEEDED in CV'
!	 write(*,*) ''
!	endif

        return
	end

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
!FFFF!                                                          !FFFF!
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
	subroutine add_O2(i_cntr,vstart,CV_list,n_CV,interval_x,
     &		interval_y,interval_z,x_Imp,x_Imp_temp,
     &		v_Imp_temp,type_Imp_temp,n_o2_need,n_Imp,
     &		EAMdebug,n_surf,atomCoord,atomDispl,x_cracktip,
     &		nnsteps,n_temp,dis_12)
	use mod_global
	implicit none
	
	logical EAMdebug
	integer i_cntr,j,irep,allgood,n_CV,n_o2_need,n_Imp,n_Al,ii
	integer nnsteps,n_temp,iii,refnum
        double precision atomCoord(ndf,*),atomDispl(ndf,*)
        double precision O_BOND,rand_x,rand_y,rand_z,CV_choice,
     &		CV_xadd,CV_yadd,vstart,n_mag,dx1,dy1,dz1,d1,dx2,
     &		dy2,dz2,d2,px,py,pz,interval_x,interval_y,interval_z,
     &		TOO_CLOSE,n_surf,cx,cy,cz,x_cracktip,imp_min
        integer, allocatable, intent(inout) :: CV_list(:,:)
        double precision, allocatable, intent(inout) :: x_Imp(:)
        double precision, allocatable, intent(inout) :: x_Imp_temp(:)
        double precision, allocatable, intent(inout) :: v_Imp_temp(:)
        double precision, allocatable, INTENT(INOUT) :: dis_12(:)
        integer, allocatable, intent(inout) :: type_Imp_temp(:)
	integer, allocatable :: addlist(:)
	logical :: NeedList,skipnode

	if(nnsteps.eq.1)then
	  imp_min = (grid_xmin)
	else
	  imp_min = (x_cracktip-20.d0)
	endif

	! Generate Random O2 Molecules:
!	if(rank.eq.0) write(*,*)' Start Adding O2...'

        if(Monolayer.and.(n_o2_need.gt.0))then

	  if(rank.eq.0) call random_seed

	  allocate(addlist(int(n_surf)))
	  addlist(:) = 0

	  ! Make list of all nodes (in 1-numnp list) that are available for Oxygen
	  j=0
          iii=0
	  do irep = 1, numnp
            NeedList=IsRelaxed(irep).ne.0
            if(NeedList)then
	      iii = iii + 1
	      skipnode = .false.
	      if(dis_12(iii).ge.(Sbarrier))then !IS ON SURFACE
	        cx = atomDispl(1,irep)+atomCoord(1,irep)
	        cy = atomDispl(2,irep)+atomCoord(2,irep)
	        cz = atomDispl(3,irep)+atomCoord(3,irep)
	        if((cx.gt.imp_min).and.(cx.lt.grid_xmax).and.
     &		   (cy.gt.grid_ymin).and.(cy.lt.grid_ymax))then
	          do ii = 1,n_Imp
                    dx1 = x_Imp((ii-1)*3+1)-cx
                    dy1 = x_Imp((ii-1)*3+2)-cy
                    dz1 = x_Imp((ii-1)*3+3)-cz
	            d1 = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1 )
	            if(d1.lt.(3.0))then
		      skipnode = .true.
		      exit
	            endif
	          enddo
	          if(.not.skipnode)then
		    j=j+1
		    addlist(j) = irep
	            if(j.eq.int(n_surf)) exit
	          endif
	        endif
	      endif
	    endif
	  enddo
	  if(j<n_surf) n_surf = j

!	 if(rank.eq.0) write(*,*)' List Made...'

	  ! Sort list RANDOMLY
	  do ii = 1, int(n_surf)-1
	     if(rank.eq.0) call RANDOM_NUMBER(rand_x)
             CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(rand_x,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
	     rand_x = rand_x * (n_surf-ii+1) ! Get random number ( 0 , (nsurf-ii+1) ]
	     refnum = int(floor(rand_x)+1) + ii - 1 ! gives some Random integer [ ii , n_surf ]
	     j = addlist(ii)
	     addlist(ii) = addlist(refnum)
	     addlist(refnum) = j
	  enddo

!	  if(rank.eq.0) write(*,*)' List Sorted...'

	  if(n_o2_need.gt.n_surf) then 
!            if(rank.eq.0) write(*,*)' Decreasing n_temp!! '
	    n_temp = n_temp - (n_o2_need-n_surf)
	    n_o2_need = n_surf
	  endif

	  do irep = 1, n_o2_need

	   if(rank.eq.0) write(*,*)' Adding O atom #',i_cntr

	   allgood=0

	   do while (allgood.eq.0)

	    allgood = 1

            !! Define Random Normal for O atom away from Al surface atom.
	     if(rank.eq.0)then
	      call RANDOM_NUMBER(rand_x)
	      call RANDOM_NUMBER(rand_y)
	      call RANDOM_NUMBER(rand_z)
	     endif
             CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(rand_x,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(rand_y,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(rand_z,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
             rand_x = rand_x * 2.0 - 1.0
             rand_y = rand_y * 2.0 - 1.0
             rand_z = rand_z * 2.0 - 1.0
             n_mag = sqrt(rand_x*rand_x+rand_y*rand_y+rand_z*rand_z)
             px = rand_x / n_mag
             py = rand_y / n_mag
             pz = rand_z / n_mag


            x_Imp_temp((i_cntr-1)*3+1)=(atomDispl(1,addlist(irep))
     &		+atomCoord(1,addlist(irep)))+px*(1.75)
            x_Imp_temp((i_cntr-1)*3+2)=(atomDispl(2,addlist(irep))
     &		+atomCoord(2,addlist(irep)))+py*(1.75)
            x_Imp_temp((i_cntr-1)*3+3)=(atomDispl(3,addlist(irep))
     &		+atomCoord(3,addlist(irep)))+pz*(1.75)

	    v_Imp_temp((i_cntr-1)*3+1) = 0.0
	    v_Imp_temp((i_cntr-1)*3+2) = 0.0
	    v_Imp_temp((i_cntr-1)*3+3) = 0.0
            type_Imp_temp(i_cntr) = 2

	    ! Check that new Imp Atoms are not too close to existing atoms:
	    do j=1,n_Imp
              dx1 = x_Imp_temp((i_cntr-1)*3+1) - x_Imp((j-1)*3+1)
              dy1 = x_Imp_temp((i_cntr-1)*3+2) - x_Imp((j-1)*3+2)
              dz1 = x_Imp_temp((i_cntr-1)*3+3) - x_Imp((j-1)*3+3)
	      d1 = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1 )
	      if(d1.le.(2.5))then
	        allgood = 0
	      endif
	    enddo
	    do j=1,i_cntr-1
              dx1 = x_Imp_temp((i_cntr-1)*3+1) - x_Imp_temp((j-1)*3+1)
              dy1 = x_Imp_temp((i_cntr-1)*3+2) - x_Imp_temp((j-1)*3+2)
              dz1 = x_Imp_temp((i_cntr-1)*3+3) - x_Imp_temp((j-1)*3+3)
	      d1 = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1 )
	      if(d1.le.(2.5))then
	        allgood = 0
	      endif
	    enddo
	    n_Al = 0
	    do j = 1, numnp
              NeedList=IsRelaxed(irep).ne.0
              if(NeedList)then
	        n_Al=n_Al+1
		if(j.ne.addlist(i_cntr))then
	          dx1 = x_Imp_temp((i_cntr-1)*3+1) - 
     &			(atomDispl(1,n_Al)+atomCoord(1,n_Al))
                  dy1 = x_Imp_temp((i_cntr-1)*3+2) - 
     &			(atomDispl(2,n_Al)+atomCoord(2,n_Al))
                  dz1 = x_Imp_temp((i_cntr-1)*3+3) - 
     &			(atomDispl(3,n_Al)+atomCoord(3,n_Al))
	          d1 = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1 )
	          if(d1.le.(2.5))then
	            allgood = 0
	          endif
		endif
	      endif
	    enddo

	  enddo ! do while

          i_cntr=i_cntr+1

	 enddo

	 deallocate(addlist)

	elseif(n_o2_need.gt.0)then ! Monolayer

	O_BOND = 1.21
	if(rank.eq.0) call random_seed

	do irep = 1,n_o2_need
	 allgood=0
	 do while (allgood.eq.0)
	  allgood=1

	  ! Choose Random CV_cell for O2 Molecule
	  if(rank.eq.0)then
	   call RANDOM_NUMBER(rand_x)
	  endif
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(rand_x,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
	  CV_choice = CEILING(rand_x * n_CV)
	  CV_xadd = grid_xmin + (CV_list(CV_choice,1)-1) * grid_dx
	  CV_yadd = grid_ymin + (CV_list(CV_choice,2)-1) * grid_dx

          ! Define Random Coordinate for 1st O atom
	  if(rank.eq.0)then
	   call RANDOM_NUMBER(rand_x)
	   call RANDOM_NUMBER(rand_y)
	   call RANDOM_NUMBER(rand_z)
	  endif
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(rand_x,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(rand_y,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(rand_z,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

	  x_Imp_temp((i_cntr-1)*3+1) = rand_x*interval_x + CV_xadd
	  x_Imp_temp((i_cntr-1)*3+2) = rand_y*interval_y + CV_yadd
	  x_Imp_temp((i_cntr-1)*3+3) = rand_z*interval_z + 0.0
	  v_Imp_temp((i_cntr-1)*3+1) = vstart 
	  v_Imp_temp((i_cntr-1)*3+2) = 0.0
	  v_Imp_temp((i_cntr-1)*3+3) = 0.0
          type_Imp_temp(i_cntr) = 2

          ! Define Random Normal for 2nd O atom
	  if(.not.Monolayer)then
	    if(rank.eq.0)then
	     call RANDOM_NUMBER(rand_x)
	     call RANDOM_NUMBER(rand_y)
	     call RANDOM_NUMBER(rand_z)
	    endif
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(rand_x,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(rand_y,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(rand_z,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
            rand_x = rand_x * 2.0 - 1.0
            rand_y = rand_y * 2.0 - 1.0
            rand_z = rand_z * 2.0 - 1.0
            n_mag = sqrt(rand_x*rand_x+rand_y*rand_y+rand_z*rand_z)
            px = rand_x / n_mag
            py = rand_y / n_mag
            pz = rand_z / n_mag
            x_Imp_temp(i_cntr*3+1) = x_Imp_temp((i_cntr-1)*3+1) 
     &		+ px * O_BOND
            x_Imp_temp(i_cntr*3+2) = x_Imp_temp((i_cntr-1)*3+2) 
     &		+ py * O_BOND
            x_Imp_temp(i_cntr*3+3) = x_Imp_temp((i_cntr-1)*3+3) 
     &		+ pz * O_BOND
	    v_Imp_temp(i_cntr*3+1) = vstart
	    v_Imp_temp(i_cntr*3+2) = 0.0
	    v_Imp_temp(i_cntr*3+3) = 0.0
            type_Imp_temp(i_cntr+1) = 2
	  endif

	  ! Check that new Imp Atoms are not too close to existing atoms:
	  TOO_CLOSE = 2.0 
	  do j=1,n_Imp
            dx1 = x_Imp_temp((i_cntr-1)*3+1) - x_Imp((j-1)*3+1)
            dy1 = x_Imp_temp((i_cntr-1)*3+2) - x_Imp((j-1)*3+2)
            dz1 = x_Imp_temp((i_cntr-1)*3+3) - x_Imp((j-1)*3+3)
	    d1 = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1 )
	    if(Monolayer)then
	      d2 = 2.0 * TOO_CLOSE
	    else
              dx2 = x_Imp_temp(i_cntr*3+1) - x_Imp((j-1)*3+1)
              dy2 = x_Imp_temp(i_cntr*3+2) - x_Imp((j-1)*3+2)
              dz2 = x_Imp_temp(i_cntr*3+3) - x_Imp((j-1)*3+3)
	      d2 = sqrt( dx2*dx2 + dy2*dy2 + dz2*dz2 )
	    endif
	    if((d1.le.TOO_CLOSE).or.(d2.le.TOO_CLOSE))then
	      allgood = 0
	    endif
	  enddo
	  do j=1,i_cntr-1
            dx1 = x_Imp_temp((i_cntr-1)*3+1) - x_Imp_temp((j-1)*3+1)
            dy1 = x_Imp_temp((i_cntr-1)*3+2) - x_Imp_temp((j-1)*3+2)
            dz1 = x_Imp_temp((i_cntr-1)*3+3) - x_Imp_temp((j-1)*3+3)
	    d1 = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1 )
	    if(Monolayer)then
	      d2 = 2.0 * TOO_CLOSE
	    else
              dx2 = x_Imp_temp(i_cntr*3+1) - x_Imp_temp((j-1)*3+1)
              dy2 = x_Imp_temp(i_cntr*3+2) - x_Imp_temp((j-1)*3+2)
              dz2 = x_Imp_temp(i_cntr*3+3) - x_Imp_temp((j-1)*3+3)
	      d2 = sqrt( dx2*dx2 + dy2*dy2 + dz2*dz2 )
	    endif
	    if((d1.le.TOO_CLOSE).or.(d2.le.TOO_CLOSE))then
	      allgood = 0
	    endif
	  enddo

	 enddo

         i_cntr=i_cntr+2

	enddo

	endif ! Monolayer

!	if(rank.eq.0) write(*,*)' Done Adding O2...'

        return
	end

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
!FFFF!                                                          !FFFF!
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
	subroutine add_H2O(i_cntr,vstart,CV_list,n_CV,interval_x,
     &		interval_y,interval_z, x_Imp,x_Imp_temp,
     &		v_Imp_temp,type_Imp_temp,n_h2o_need,n_Imp)
	use mod_global
	implicit none
	
	integer i_cntr,j,irep,allgood,n_CV,n_h2o_need,n_Imp
        double precision OH_BOND,rand_x,rand_y,rand_z,CV_choice,
     &		CV_xadd,CV_yadd,vstart,n_mag,dx1,dy1,dz1,d1,dx2,
     &		dy2,dz2,d2,px,py,pz,interval_x,interval_y,interval_z,
     &		TOO_CLOSE,PI,H2O_ANGLE,dx3,dy3,dz3,d3,
     &		tx1,ty1,tz1,tx2,ty2,tz2
        integer, allocatable, intent(inout) :: CV_list(:,:)
        double precision, allocatable, intent(inout) :: x_Imp(:)
        double precision, allocatable, intent(inout) :: x_Imp_temp(:)
        double precision, allocatable, intent(inout) :: v_Imp_temp(:)
        integer, allocatable, intent(inout) :: type_Imp_temp(:)

!--	Functions
	double precision rotate3D 

	! Generate Random H2O Molecules:
!	if(rank.eq.0) write(*,*)' Start Adding H2O...'
	OH_BOND = 0.96            ! Ang
	H2O_ANGLE = 1.823         ! Radians
	PI = 3.1415926535897932   ! Radians
	if(n_h2o_need.gt.0)then
	do irep = 1,n_h2o_need
	 allgood=0
	 do while (allgood.eq.0)
	  allgood=1

	  ! Choose Random CV_cell for H2O Molecule
	  if(rank.eq.0)then
	   call RANDOM_NUMBER(rand_x)
	  endif
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(rand_x,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
	  CV_choice = CEILING(rand_x * n_CV)
	  CV_xadd = grid_xmin + (CV_list(CV_choice,1)-1) * grid_dx
	  CV_yadd = grid_ymin + (CV_list(CV_choice,2)-1) * grid_dx

          ! Define Random Coordinate for O atom
	  if(rank.eq.0)then
	   call RANDOM_NUMBER(rand_x)
	   call RANDOM_NUMBER(rand_y)
	   call RANDOM_NUMBER(rand_z)
	  endif
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(rand_x,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(rand_y,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(rand_z,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
	  x_Imp_temp((i_cntr-1)*3+1) = rand_x*interval_x + CV_xadd
	  x_Imp_temp((i_cntr-1)*3+2) = rand_y*interval_y + CV_yadd
	  x_Imp_temp((i_cntr-1)*3+3) = rand_z*interval_z + 0.0
	  v_Imp_temp((i_cntr-1)*3+1) = vstart 
	  v_Imp_temp((i_cntr-1)*3+2) = 0.0
	  v_Imp_temp((i_cntr-1)*3+3) = 0.0
          type_Imp_temp(i_cntr) = 2

          ! Insert Hydrogen atoms in at UNROTATED Position */
          tx1 = OH_BOND 
          ty1 = 0.0 
          tz1 = 0.0 
          type_Imp_temp(i_cntr+1) = 3
 
          tx2 = OH_BOND * cos(H2O_ANGLE) 
          ty2 = OH_BOND * sin(H2O_ANGLE) 
          tz2 = 0.0 
          type_Imp_temp(i_cntr+2) = 3 

          ! Define 3 Random Rotation Angles*/
	  if(rank.eq.0)then
	   call RANDOM_NUMBER(rand_x)
	   call RANDOM_NUMBER(rand_y)
	   call RANDOM_NUMBER(rand_z)
	  endif
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(rand_x,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(rand_y,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(rand_z,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
          rand_x = rand_x * PI 
          rand_y = rand_y * PI 
          rand_z = rand_z * PI 

          ! Create ROTATED H Atoms */
          x_Imp_temp((i_cntr+0)*3+1) = rotate3D(rand_x,rand_y,rand_z,
     &			tx1,ty1,tz1,0) + x_Imp_temp((i_cntr-1)*3+1)
          x_Imp_temp((i_cntr+0)*3+2) = rotate3D(rand_x,rand_y,rand_z,
     &			tx1,ty1,tz1,1) + x_Imp_temp((i_cntr-1)*3+2)
          x_Imp_temp((i_cntr+0)*3+3) = rotate3D(rand_x,rand_y,rand_z,
     &			tx1,ty1,tz1,2) + x_Imp_temp((i_cntr-1)*3+3)
	  v_Imp_temp((i_cntr)*3+1) = vstart 
	  v_Imp_temp((i_cntr)*3+2) = 0.0
	  v_Imp_temp((i_cntr)*3+3) = 0.0

          x_Imp_temp((i_cntr+1)*3+1) = rotate3D(rand_x,rand_y,rand_z,
     &			tx2,ty2,tz2,0) + x_Imp_temp((i_cntr-1)*3+1)
          x_Imp_temp((i_cntr+1)*3+2) = rotate3D(rand_x,rand_y,rand_z,
     &			tx2,ty2,tz2,1) + x_Imp_temp((i_cntr-1)*3+2)
          x_Imp_temp((i_cntr+1)*3+3) = rotate3D(rand_x,rand_y,rand_z,
     &			tx2,ty2,tz2,2) + x_Imp_temp((i_cntr-1)*3+3)
	  v_Imp_temp((i_cntr+1)*3+1) = vstart
	  v_Imp_temp((i_cntr+1)*3+2) = 0.0
	  v_Imp_temp((i_cntr+1)*3+3) = 0.0

	  ! Check that new Imp Atoms are not too close to existing atoms:
	  TOO_CLOSE = 1.75
	  do j=1,n_Imp
            dx1 = x_Imp_temp((i_cntr-1)*3+1) - x_Imp((j-1)*3+1)
            dy1 = x_Imp_temp((i_cntr-1)*3+2) - x_Imp((j-1)*3+2)
            dz1 = x_Imp_temp((i_cntr-1)*3+3) - x_Imp((j-1)*3+3)
	    d1 = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1 )
            dx2 = x_Imp_temp((i_cntr)*3+1) - x_Imp((j-1)*3+1)
            dy2 = x_Imp_temp((i_cntr)*3+2) - x_Imp((j-1)*3+2)
            dz2 = x_Imp_temp((i_cntr)*3+3) - x_Imp((j-1)*3+3)
	    d2 = sqrt( dx2*dx2 + dy2*dy2 + dz2*dz2 )
            dx3 = x_Imp_temp((i_cntr+1)*3+1) - x_Imp((j-1)*3+1)
            dy3 = x_Imp_temp((i_cntr+1)*3+2) - x_Imp((j-1)*3+2)
            dz3 = x_Imp_temp((i_cntr+1)*3+3) - x_Imp((j-1)*3+3)
	    d3 = sqrt( dx3*dx3 + dy3*dy3 + dz3*dz3 )
	    if((d1.le.TOO_CLOSE).or.(d2.le.TOO_CLOSE)
     &			.or.(d3.le.TOO_CLOSE))then
	      allgood = 0
	    endif
	  enddo
	  do j=1,i_cntr-1
            dx1 = x_Imp_temp((i_cntr-1)*3+1) - x_Imp_temp((j-1)*3+1)
            dy1 = x_Imp_temp((i_cntr-1)*3+2) - x_Imp_temp((j-1)*3+2)
            dz1 = x_Imp_temp((i_cntr-1)*3+3) - x_Imp_temp((j-1)*3+3)
	    d1 = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1 )
            dx2 = x_Imp_temp((i_cntr)*3+1) - x_Imp_temp((j-1)*3+1)
            dy2 = x_Imp_temp((i_cntr)*3+2) - x_Imp_temp((j-1)*3+2)
            dz2 = x_Imp_temp((i_cntr)*3+3) - x_Imp_temp((j-1)*3+3)
	    d2 = sqrt( dx2*dx2 + dy2*dy2 + dz2*dz2 )
            dx3 = x_Imp_temp((i_cntr+1)*3+1) - x_Imp_temp((j-1)*3+1)
            dy3 = x_Imp_temp((i_cntr+1)*3+2) - x_Imp_temp((j-1)*3+2)
            dz3 = x_Imp_temp((i_cntr+1)*3+3) - x_Imp_temp((j-1)*3+3)
	    d3 = sqrt( dx3*dx3 + dy3*dy3 + dz3*dz3 )
	    if((d1.le.TOO_CLOSE).or.(d2.le.TOO_CLOSE)
     &			.or.(d3.le.TOO_CLOSE))then
	      allgood = 0
	    endif
	  enddo

	 enddo
         i_cntr=i_cntr+3
	enddo
	endif
!	if(rank.eq.0) write(*,*)' Done Adding H2O...'

        return
	end

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
!FFFF!                                                          !FFFF!
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
	subroutine cp_temp_arrays(i_cntr,n_Imp,x_Imp,v_Imp,type_Imp,
     &		x_Imp_temp,v_Imp_temp,type_Imp_temp,n_H,n_O)
	use mod_global
	implicit none
	
	integer i_cntr,j,n_Imp,n_new,n_H,n_O
        double precision, allocatable, intent(inout) :: x_Imp_temp(:)
        double precision, allocatable, intent(inout) :: v_Imp_temp(:)
        integer, allocatable, intent(inout) :: type_Imp_temp(:)
        double precision, allocatable, intent(inout) :: x_Imp(:)
        double precision, allocatable, intent(inout) :: v_Imp(:)
        integer, allocatable, intent(inout) :: type_Imp(:)

	! Copy temp arrays over to perminant arrays:
!	if(rank.eq.0) write(*,*)' Copy Over Arrays...'
	if(n_Imp.gt.0)then
         do j = 1, n_Imp
	  x_Imp_temp((i_cntr-1)*3+1) = x_Imp((j-1)*3+1)
	  x_Imp_temp((i_cntr-1)*3+2) = x_Imp((j-1)*3+2)
	  x_Imp_temp((i_cntr-1)*3+3) = x_Imp((j-1)*3+3)
	  v_Imp_temp((i_cntr-1)*3+1) = v_Imp((j-1)*3+1)
	  v_Imp_temp((i_cntr-1)*3+2) = v_Imp((j-1)*3+2)
	  v_Imp_temp((i_cntr-1)*3+3) = v_Imp((j-1)*3+3)
	  type_Imp_temp(i_cntr) = type_Imp(j)
	  i_cntr=i_cntr+1
215	 continue
	enddo
	endif
	n_new=i_cntr-1 ! New total # of Imp atoms
        if(allocated(x_Imp)) deallocate(x_Imp)
        if(allocated(v_Imp)) deallocate(v_Imp)
        if(allocated(type_Imp)) deallocate(type_Imp)
	n_Imp = n_new
        if(.not.allocated(x_Imp)) allocate(x_Imp(n_Imp*3))
        if(.not.allocated(v_Imp)) allocate(v_Imp(n_Imp*3))
        if(.not.allocated(type_Imp)) allocate(type_Imp(n_Imp))
	n_H=0
	n_O=0
	if(n_Imp.gt.0)then
	 do j=1,n_Imp
	  x_Imp((j-1)*3+1) = x_Imp_temp((j-1)*3+1)
	  x_Imp((j-1)*3+2) = x_Imp_temp((j-1)*3+2)
	  x_Imp((j-1)*3+3) = x_Imp_temp((j-1)*3+3)
	  v_Imp((j-1)*3+1) = v_Imp_temp((j-1)*3+1)
	  v_Imp((j-1)*3+2) = v_Imp_temp((j-1)*3+2)
	  v_Imp((j-1)*3+3) = v_Imp_temp((j-1)*3+3)
	  type_Imp(j) = type_Imp_temp(j)
	  if(type_Imp(j).eq.2) n_O=n_O+1
	  if(type_Imp(j).eq.3) n_H=n_H+1
	 enddo
	endif
        if(allocated(x_Imp_temp)) deallocate(x_Imp_temp)
        if(allocated(v_Imp_temp)) deallocate(v_Imp_temp)
        if(allocated(type_Imp_temp)) deallocate(type_Imp_temp)
        if(rank.eq.0) then
!	 write(*,*) ''
!	 write(*,*) '  New n_Imp: ',n_Imp
!	 write(*,*) '  New n_O: ',n_O
!	 write(*,*) '  New n_H: ',n_H
!	 write(*,*) ''
!	 write(*,*)' Done Copying Over Arrays...'
	endif
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

        return
	end

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
!FFFF!                                                          !FFFF!
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
	subroutine write_geometry_input(n_atoms,n_Imp,padIDmin,
     &		padIDmax,velFile,EAMdebug,n_atoms_lmp,n_O,n_H,
     &		atomCoord,atomDispl,x_Imp,v_Imp,type_Imp)
	use mod_global
	implicit none
	
	logical Needlist,EAMdebug
	integer i,irep,n_Imp,n_atoms,velFile,itype,n_atoms_lmp,
     &		padIDmin,padIDmax,n_O,n_H
        double precision xmaxx,xminx,ymaxy,yminy,atomCoord(ndf,*),
     &		atomDispl(ndf,*),cx,cy,cz,chg,vx,vy,vz
        double precision, allocatable, intent(inout) :: x_Imp(:)
        double precision, allocatable, intent(inout) :: v_Imp(:)
        integer, allocatable, intent(inout) :: type_Imp(:)

	n_atoms=0
	padIDmin=numnp
	padIDmax=0
	xmaxx=-1000.d0
	xminx= 1000.d0
	ymaxy=-1000.d0
	yminy= 1000.d0
        do irep = 1, numnp
           NeedList=IsRelaxed(irep).ne.0
           if(NeedList)then
	     n_atoms=n_atoms+1
	     if(IsRelaxed(irep).eq.-1)then
               if(n_atoms.gt.padIDmax)padIDmax=n_atoms
               if(n_atoms.lt.padIDmin)padIDmin=n_atoms
	     endif
	     if(atomCoord(1,irep).gt.xmaxx) xmaxx=atomCoord(1,irep)
	     if(atomCoord(1,irep).lt.xminx) xminx=atomCoord(1,irep)
	     if(atomCoord(2,irep).gt.ymaxy) ymaxy=atomCoord(2,irep)
	     if(atomCoord(2,irep).lt.yminy) yminy=atomCoord(2,irep)
           endif
        enddo
	print *, 'xmaxx=',xmaxx,' xminx=',xminx,' ymaxy=',ymaxy,
     &  'yminy=',yminy      
        if(rank.eq.0) then
	 write(*,*) ' Writing POSCAR'
	 write(*,*) ' Starting with ',n_atoms,' Al atoms'
	 write(*,*) ' Starting with ',n_Imp,' Impurity atoms'
	 write(*,*) '               ',n_O,'   O atoms'
	 write(*,*) '               ',n_H,'   H atoms'
	endif
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	n_atoms=n_atoms+n_Imp

	! Write POSCAR
	if(rank.eq.0)then
	 open(unit=217, file='POSCAR', status='replace')  
         if(velFile.eq.1)then
	  open(unit=214, file='VEL', status='unknown')  
         endif
         write(217,*) 'Lammps-CADD_EAM+LJ_version'
	 write(217,*) ''
	 write(217,*) n_atoms,' atoms'
	 write(217,*) '14 atom types'
	 write(217,*) ''
	 write(217,*) xminx-30.0,xmaxx+30.0,' xlo xhi'
	 write(217,*) yminy-30.0,ymaxy+30.0,' ylo yhi'
	 write(217,*) 0.0,z_length,' zlo zhi'
!	 write(217,*) -z_length,z_length,' zlo zhi'
	 write(217,*) ''
	 write(217,*) 'Masses'
	 write(217,*) ''
	 write(217,*) '1 26.9820'
	 write(217,*) '2 15.9994'
	 write(217,*) '3 1.0030'
	 write(217,*) '4 26.9820'
	 write(217,*) '5 26.9820'
	 write(217,*) '6 26.9820'
	 write(217,*) '7 26.9820'
	 write(217,*) '8 26.9820'
	 write(217,*) '9 26.9820'
	 write(217,*) '10 26.9820'
	 write(217,*) '11 26.9820'
	 write(217,*) '12 26.9820'
	 write(217,*) '13 26.9820'
	 write(217,*) '14 26.9820'
	 write(217,*) ''
	 write(217,*) 'Atoms'
	 write(217,*) ''
	endif

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	! print Al coordinates
	i=0
	do irep = 1, numnp
          NeedList=IsRelaxed(irep).ne.0
          if(NeedList)then
	    i=i+1
	    itype=1 !Al_type(i)
	    cx = atomCoord(1,irep) + atomDispl(1,irep)
	    cy = atomCoord(2,irep) + atomDispl(2,irep)
	    cz = atomCoord(3,irep) + atomDispl(3,irep)
	    chg=0.d0
	    if(rank.eq.0) write(217,'(I6,1X,I6,1X,F16.8,1X,F16.8,
     &	    1X,F16.8,1X,F16.8)')
     &      i,itype,chg,cx,cy,cz
	  endif
	enddo

	! Now, add impurity corrdinates...
	do irep = 1, n_Imp
	    i=i+1
	    itype=type_Imp(irep)
	    cx = x_Imp((irep-1)*3+1)
	    cy = x_Imp((irep-1)*3+2)
	    cz = x_Imp((irep-1)*3+3)
	    chg=0.d0
	    if(rank.eq.0) write(217,'(I6,1X,I6,1X,F16.8,1X,F16.8,
     &	    1X,F16.8,1X,F16.8)')
     &      i,itype,chg,cx,cy,cz
	enddo

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	if(rank.eq.0) write(217,*) ''
	if(rank.eq.0) write(217,*) 'Velocities'
	if(rank.eq.0) write(217,*) ''
	! print Al Velocities
	i=0
	do irep = 1, numnp
          NeedList=IsRelaxed(irep).ne.0
          if(NeedList)then
	    i=i+1
            if(velFile.eq.1)then
              if(rank.eq.0) read(214,*) vx, vy, vz
            else
              vx=0.0
              vy=0.0
              vz=0.0
            endif 
	    if(rank.eq.0)write(217,'(I6,1X,F16.8,1X,F16.8,1X,F16.8)')
     &			i,vx,vy,vz
	  endif
	enddo
	! Now, add impurity velocities...
	do irep = 1, n_Imp
	    i=i+1
	    vx = v_Imp((irep-1)*3+1)
	    vy = v_Imp((irep-1)*3+2)
	    vz = v_Imp((irep-1)*3+3)
	    if(rank.eq.0)write(217,'(I6,1X,F16.8,1X,F16.8,1X,F16.8)')
     &			i,vx,vy,vz
	enddo

        if(velFile.eq.1)then
	 if(rank.eq.0) close(214)  
        endif
	if(rank.eq.0) close(217)
	n_atoms_lmp = i

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        return
	end

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
!FFFF!                                                          !FFFF!
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
	subroutine write_lammps_input(EAMdebug,padIDmin,
     &		padIDmax,mdTemp,timestep,dumpcg,mdstep)
	use mod_global
	use lammps
        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
	implicit none
	
	logical EAMdebug, ompON
	integer padIDmin,padIDmax,dumpcg,mdstep,qeqi
        double precision mdTemp,timestep,gammao,gamsca
        character *80 :: str001
        double precision kappa,cutoff1,cutoff2
        double precision OOepsilon,OOsigma,OOcutoff1,OOcutoff2
        double precision OxOxepsilon,OxOxsigma,OxOxcutoff1,OxOxcutoff2
        double precision AlOepsilon,AlOsigma,AlOcutoff1,AlOcutoff2

	ompON = .false.

	gammao = 1000.0
	if(EAMdebug.or.(qeqi.eq.1))then
	  gamsca = 0.001
	else
	  gamsca = 1.0
	endif

	if(Adhes)then

	if (EAMdebug.and.(rank.eq.0)) then
	 open(unit=711, file='lj.params', status='unknown')
	 read(711,*)kappa,cutoff1,cutoff2
	 read(711,*)OxOxepsilon,OxOxsigma,OxOxcutoff1,OxOxcutoff2
         read(711,*)AlOepsilon,AlOsigma,AlOcutoff1,AlOcutoff2
	 close(711)
	endif

	else !not adhes

	if (EAMdebug.and.(rank.eq.0)) then
	 open(unit=711, file='lj.params', status='unknown')
	 read(711,*)kappa,cutoff1,cutoff2
	 read(711,*)OOepsilon,OOsigma,OOcutoff1,OOcutoff2
         read(711,*)AlOepsilon,AlOsigma,AlOcutoff1,AlOcutoff2
	 close(711)
	endif

	endif !Adhes


	if (EAMdebug) then
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(kappa,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(cutoff1,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(cutoff2,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(OOepsilon,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(OOsigma,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(OOcutoff1,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(OOcutoff2,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

	CALL MPI_BCAST(OxOxepsilon,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(OxOxsigma,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(OxOxcutoff1,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(OxOxcutoff2,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

        CALL MPI_BCAST(AlOepsilon,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(AlOsigma,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(AlOcutoff1,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(AlOcutoff2,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
	endif

	do qeqi = 0, 1

	! Write in.lmp file
	if(rank.eq.0) then
	if(qeqi.eq.0)then
	  open(unit=217, file='in.lmp', status='replace')  
	else
	  open(unit=217, file='in.qeq', status='replace') 
	endif
        write(217,*) '# CADD + LAMMPS'
        write(217,*) '# .....'
	write(217,*) ''
	!if((qeqi.eq.0).and.ompON) then
	if((qeqi.eq.0)) then
	  write(217,*) 'package omp *'
	endif
	write(217,*) ''
	if((.not.EAMdebug).or.(qeqi.eq.1))then
	  write(217,*) 'units		real'
	else
	  write(217,*) 'units		metal'
	endif
	write(217,*) 'atom_modify map array'
	write(217,*) 'dimension	        3'
	write(217,*) 'boundary	        p p p'
	write(217,*) 'atom_style	charge'
	write(217,*) 'read_data	        POSCAR'
	write(217,*) ''
	if((.not.EAMdebug).or.(qeqi.eq.1))then
	 write(217,*) 'pair_style
     &	reax/c NULL safezone 3.2 mincap 200'
	 write(217,*) 'pair_coeff * * 
     &	materials/ffield.reax.AlO &'
        write(217,*) 'Al O H Al Al Al Al Al Al Al Al Al Al Al'
!        write(217,*) 'Al O H Al Al Al Al Al Al Al Ox Ox Al Al'

	else !meaning we are writing in.lmpOB

         if(explicitOxy)then
	   write(217,*) 'pair_style hybrid &'
	 else if(Adhes)then
	   write(217,*) 'pair_style hybrid &'

	 else
	   write(217,*) 'pair_style hybrid/overlay &'
	 endif
	 if(ompON)then
	   write(217,*) 'lj/cut/coul/debye/omp &'
	 else
	   write(217,*) 'lj/cut/coul/debye &'
	 endif
	 write(217,*) kappa,' &'
	 write(217,*) cutoff1,' &'
	 write(217,*) cutoff2,' &'
	 if(ompON)then
	   write(217,*) ' eam/alloy/omp'
	 else
	   write(217,*) ' eam/alloy'
	 endif

	 !! LJ Law for UNWANTED interactions:
	 !if(ompON)then
	 !  write(217,*) 'pair_coeff * * lj/cut/coul/debye/omp &'
	 !else
	 !  write(217,*) 'pair_coeff * * lj/cut/coul/debye &'
	 !endif
	 !write(217,*) OOepsilon*0.01d0,' &'
	 !write(217,*) OOsigma*0.01d0,' &'
	 !write(217,*) 0.1d0,' &'
	 !write(217,*) 0.1d0 

	 ! LJ Law for Unwanted Interactions...
	 if(.not.explicitOxy)then
	  if(ompON)then
	     write(217,*) 'pair_coeff * * lj/cut/coul/debye/omp &'
	  else
	     write(217,*) 'pair_coeff * * lj/cut/coul/debye &'
	  endif
	  write(217,*) 0.000,' &'
	  write(217,*) 0.001,' &'
	  write(217,*) 0.001,' &'
	  write(217,*) 0.001
	 endif

	 ! LJ Law for Al-O interactions:
	 if(Adhes)then
	 	 if(ompON)then
		   write(217,*) 'pair_coeff * * eam/alloy/omp &'
		 else
		   write(217,*) 'pair_coeff * * eam/alloy &'
		 endif
		 write(217,*) './materials/Al_EA.eam.alloy &' 
		 write(217,*) 'Al NULL Al Al Al Al Al Al Al Al Al Al Al Al'
!		 write(217,*) 'Al NULL Al Al Al Al Al Al Al Al Ox Ox Al Al'

		   if(ompON)then
		     write(217,*) 'pair_coeff 11 14 lj/cut/coul/debye/omp &'
			 write(217,*) AlOepsilon,' &'
			 write(217,*) AlOsigma,' &'
			 write(217,*) AlOcutoff1,' &'
			 write(217,*) AlOcutoff2


		     write(217,*) 'pair_coeff 12 13 lj/cut/coul/debye/omp &'
			 write(217,*) AlOepsilon,' &'
			 write(217,*) AlOsigma,' &'
			 write(217,*) AlOcutoff1,' &'
			 write(217,*) AlOcutoff2


		   else
		     write(217,*) 'pair_coeff 11 14 lj/cut/coul/debye &'
			 write(217,*) AlOepsilon,' &'
			 write(217,*) AlOsigma,' &'
			 write(217,*) AlOcutoff1,' &'
			 write(217,*) AlOcutoff2


		     write(217,*) 'pair_coeff 12 13 lj/cut/coul/debye &'
			 write(217,*) AlOepsilon,' &'
			 write(217,*) AlOsigma,' &'
			 write(217,*) AlOcutoff1,' &'
			 write(217,*) AlOcutoff2


		   endif
		
		 ! LJ Law for O-O interactions (in "atomistic" stadium region only):
!		    if(ompON)then
!		     write(217,*) 'pair_coeff 11 11 lj/cut/coul/debye/omp &'
!			 write(217,*) OOepsilon,' &'
!			 write(217,*) OOsigma,' &'
!			 write(217,*) OOcutoff1,' &'
!			 write(217,*) OOcutoff2
!
!
!		     write(217,*) 'pair_coeff 12 12 lj/cut/coul/debye/omp &'
!			 write(217,*) OOepsilon,' &'
!			 write(217,*) OOsigma,' &'
!			 write(217,*) OOcutoff1,' &'
!			 write(217,*) OOcutoff2
!
!
!		   else
!		     write(217,*) 'pair_coeff 11 11 lj/cut/coul/debye &'
!			 write(217,*) OOepsilon,' &'
!			 write(217,*) OOsigma,' &'
!			 write(217,*) OOcutoff1,' &'
!			 write(217,*) OOcutoff2
!
!
!		     write(217,*) 'pair_coeff 12 12 lj/cut/coul/debye &'
!			 write(217,*) OOepsilon,' &'
!			 write(217,*) OOsigma,' &'
!			 write(217,*) OOcutoff1,' &'
!			 write(217,*) OOcutoff2
!
!
!		   endif

		 ! LJ Law for cross-crack ox interactions 
		    if(ompON)then
		     write(217,*) 'pair_coeff 11 12 lj/cut/coul/debye/omp &'
			 write(217,*) OxOxepsilon,' &'
			 write(217,*) OxOxsigma,' &'
			 write(217,*) OxOxcutoff1,' &'
			 write(217,*) OxOxcutoff2


!		     write(217,*) 'pair_coeff 12 11 lj/cut/coul/debye/omp &'
!			 write(217,*) OOepsilon,' &'
!			 write(217,*) OOsigma,' &'
!			 write(217,*) OOcutoff1,' &'
!			 write(217,*) OOcutoff2


		   else
		     write(217,*) 'pair_coeff 11 12 lj/cut/coul/debye &'
			 write(217,*) OxOxepsilon,' &'
			 write(217,*) OxOxsigma,' &'
			 write(217,*) OxOxcutoff1,' &'
			 write(217,*) OxOxcutoff2


!		     !write(217,*) 'pair_coeff 12 11 lj/cut/coul/debye &'
!		     write(217,*) 'pair_coeff 12 12 lj/cut/coul/debye &'
!			 write(217,*) OOepsilon,' &'
!			 write(217,*) OOsigma,' &'
!			 write(217,*) OOcutoff1,' &'
!			 write(217,*) OOcutoff2


		   endif
		 
	 else
		 if(explicitOxy)then
		   if(ompON)then
		     write(217,*) 'pair_coeff 2 * lj/cut/coul/debye/omp &'
		   else
		     write(217,*) 'pair_coeff 2 * lj/cut/coul/debye &'
		   endif
		 else
		   if(ompON)then
		     write(217,*) 'pair_coeff 4 * lj/cut/coul/debye/omp &'
		   else
		     write(217,*) 'pair_coeff 4 * lj/cut/coul/debye &'
		   endif
		 endif
		 write(217,*) AlOepsilon,' &'
		 write(217,*) AlOsigma,' &'
		 write(217,*) AlOcutoff1,' &'
		 write(217,*) AlOcutoff2

		 ! LJ Law for O-O interactions (in "atomistic" stadium region only):
		 if(explicitOxy)then
		   if(ompON)then
		     write(217,*) 'pair_coeff 2 2 lj/cut/coul/debye/omp &'
		   else
		     write(217,*) 'pair_coeff 2 2 lj/cut/coul/debye &'
		   endif
		 else
		   if(ompON)then
		     write(217,*) 'pair_coeff 4 4 lj/cut/coul/debye/omp &'
		   else
		     write(217,*) 'pair_coeff 4 4 lj/cut/coul/debye &'
		   endif
		 endif
		 write(217,*) OOepsilon,' &'
		 write(217,*) OOsigma,' &'
		 write(217,*) OOcutoff1,' &'
		 write(217,*) OOcutoff2

		 if(ompON)then
		   write(217,*) 'pair_coeff * * eam/alloy/omp &'
		 else
		   write(217,*) 'pair_coeff * * eam/alloy &'
		 endif
		 write(217,*) './materials/Al_EA.eam.alloy &' 
		 write(217,*) 'Al NULL Al Al Al Al Al Al Al Al Al Al Al Al'
!		 write(217,*) 'Al NULL Al Al Al Al Al Al Al Al Ox Ox Al Al'

		endif
	endif !Adhes
        write(217,*) 'compute 1 all stress/atom NULL'

        if(qeqi.eq.0)then
	  write(217,*) ''
	  write(217,*) 'neighbor	        2.0 bin'
	  write(217,*) 'neigh_modify      every 1 delay 0 check no'
	  write(217,*) ''
	  write(217,*) 'group	gpad id >= ',padIDmin
	  write(217,*) 'group	lpad id <= ',padIDmax
	  write(217,*) 'group	pad intersect gpad lpad'
	  write(217,*) ''
          write(217,*) 'fix 1 all nve'
	endif

	if((.not.EAMdebug).or.(qeqi.eq.1))
     &	write(217,*) 'fix 2 all qeq/reax 1 0.0 10.0 1e-6 reax/c'

	if(qeqi.eq.0)then
	  write(217,*) 'fix 3 all langevin &'
	  write(217,*) mdTemp, mdTemp,' &'
	  write(217,*) gamsca*gammao,' 48279 &'
	  write(217,*) ' scale 1 ',1.0,' &' ! Normal Al atoms, dont want to dampen much
	  if(EAMdebug)then
	   write(217,*) ' scale 2 ',(4.0/6.0),' &' ! O atoms (need to dampen a lot)
	   write(217,*) ' scale 3 ',(4.0/6.0),' &' ! H atoms (need to dampen a lot)
	   write(217,*) ' scale 4 ',(4.0/6.0),' &' ! Al atoms near an Impurity (need to dampen)
          endif
	  ! Start linear gradient of "stadium" damping:
	  write(217,*) ' scale 5 ',(5.0/6.0),' &'
	  write(217,*) ' scale 6 ',(4.0/6.0),' &'
	  write(217,*) ' scale 7 ',(3.0/6.0),' &'
	  write(217,*) ' scale 8 ',(2.0/6.0),' &'
	  write(217,*) ' scale 9 ',(1.0/6.0)
	  !write(217,*) ' scale 10 ',(1.0/(gamsca*gammao)),' &' ! Al pad atoms (need to dampen a lot)
	  write(217,*) ''
	  write(217,*) 'fix 4 pad setforce 0.0 0.0 0.0'
	  write(217,*) ''
	  write(217,*) 'min_style hftn' 
	  write(217,*) 'min_modify dmax 0.25' 
	  write(217,*) 'timestep ',timestep
	  write(217,*) ''
	  write(217,*) 'dump 1 all cfg ',dumpcg,
     &    ' dump_lmp.*.cfg mass type xs ys zs vx vy vz'
	  write(217,*) 'dump_modify 1 element &'
	  write(217,*) ' Al O H Al Al Al Al Al Al Al Al Al Al Al'
!	  write(217,*) ' Al O H Al Al Al Al Al Al Al Ox Ox Al Al'
	  write(217,*) 'dump_modify 1 sort id'
	  write(217,*) ''
	endif
!	write(217,*) 'dump 2 all local ',dumpcg,'&'
!	write(217,*) '  dump_pair index c_1[1] c_1[2] c_1[3] c_1[4]'
!	write(217,*) ''
!	write(217,*) 'dump 2 all custom ',dumpcg,' &'
!	write(217,*) ' dump_STRESS.* x y z &'
!	write(217,*) ' c_1[1] c_1[2] c_1[3] c_1[4] c_1[5] c_1[6]'
!	write(217,*) 'dump_modify 2 sort id'
!	write(217,*) ''
!	if((.not.EAMdebug).or.(qeqi.eq.1))then
!	 write(217,*) 'dump 3 all custom ',dumpcg,' dump_CHARGE.* x y z type q'
!	 write(217,*) 'dump_modify 3 sort id'
!	 write(217,*) ''
!	endif
	close(217)
	endif
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	enddo

        return
	end

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
!FFFF!                                                          !FFFF!
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
	subroutine adjust_displ(n_atoms_lmp,
     &		atomCoord,atomDispl,type_lmp,x_lmp,
     &		n_Imp,x_Imp,o2_conc,iStep,dis_12,
     &          Atom_top,Atom_bot)
	use mod_global
	implicit none

	logical NeedList,keepTrying
	integer n_atoms_lmp,i,irep,n_Imp,j,s_cnt,o_cnt,o_need
	integer iStep, rand_i, n_Al
        double precision cx,cy,cz,atomCoord(ndf,*),
     &	 atomDispl(ndf,*),axmax,axmin,aymax,aymin,stad_dx,
     &	 dx,dy,dz,dmag,impcut,o2_conc,rand_d
        integer, allocatable, intent(inout) :: type_lmp(:)
        double precision, allocatable, intent(inout) :: x_lmp(:)
        double precision, allocatable, intent(inout) :: x_Imp(:)
        double precision, allocatable, INTENT(INOUT) :: dis_12(:)
        integer, allocatable :: surf_list(:)

!ER- I added this
        integer, allocatable, INTENT(INOUT) :: Atom_top(:)
        integer, allocatable, INTENT(INOUT) :: Atom_bot(:)

          !! Setup "stadium" damping (should only need to do this once):
	  stad_dx = 2.5
	  axmax=-10000.0
	  axmin= 10000.0
	  aymax=-10000.0
	  aymin= 10000.0

	  n_Al=0
	  do irep = 1, numnp
            NeedList=IsRelaxed(irep).ne.0
            if(NeedList)then
	      n_Al=n_Al+1
	      if(IsRelaxed(irep).eq.-1)then
		type_lmp(n_Al) = 10 ! pad atom
	      else
	        cx = atomCoord(1,irep)
	        cy = atomCoord(2,irep)
		if(cx.gt.axmax)axmax=cx
		if(cx.lt.axmin)axmin=cx
		if(cy.gt.aymax)aymax=cy
		if(cy.lt.aymin)aymin=cy
	      endif
	    endif
	  enddo
	  i=0
	  do irep = 1, numnp
            NeedList=IsRelaxed(irep).ne.0
            if(NeedList)then
	      i=i+1
	      if(i.gt.n_atoms_lmp) write(*,*) ' ERROR.. TOO MANY Al!'
	      if(IsRelaxed(irep).ne.-1)then
	        cx = atomCoord(1,irep)
	        cy = atomCoord(2,irep)

		if( (cx.gt.(axmax-1.0*stad_dx))
     &		.or.(cx.lt.(axmin+1.0*stad_dx))
     &		.or.(cy.gt.(aymax-1.0*stad_dx))
     &		.or.(cy.lt.(aymin+1.0*stad_dx)) )then

		  type_lmp(i) = 9

	        elseif((cx.gt.(axmax-2.0*stad_dx))
     &		.or.(cx.lt.(axmin+2.0*stad_dx))
     &		.or.(cy.gt.(aymax-2.0*stad_dx))
     &		.or.(cy.lt.(aymin+2.0*stad_dx)))then

		  type_lmp(i) = 8 

	        elseif((cx.gt.(axmax-3.0*stad_dx))
     &		.or.(cx.lt.(axmin+3.0*stad_dx))
     &		.or.(cy.gt.(aymax-3.0*stad_dx))
     &		.or.(cy.lt.(aymin+3.0*stad_dx)))then

		  type_lmp(i) = 7 

	        elseif((cx.gt.(axmax-4.0*stad_dx))
     &		.or.(cx.lt.(axmin+4.0*stad_dx))
     &		.or.(cy.gt.(aymax-4.0*stad_dx))
     &		.or.(cy.lt.(aymin+4.0*stad_dx)))then

		  type_lmp(i) = 6 

	        elseif((cx.gt.(axmax-5.0*stad_dx))
     &		.or.(cx.lt.(axmin+5.0*stad_dx))
     &		.or.(cy.gt.(aymax-5.0*stad_dx))
     &		.or.(cy.lt.(aymin+5.0*stad_dx)))then

		  type_lmp(i) = 5 

	 	endif

	      endif
	    endif
	  enddo

	  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	  if(rank.eq.0) write(*,*) 'were in adjust displacement'
	  if(rank.eq.0) write(*,*) 'Adjust Step#1 Done.'
	  if(rank.eq.0) write(*,*) 'Currently ',n_Imp,' impurity atoms'

	  ! Now, must find all Al atoms near (closer than "impcut") O or H and set to type "4"
	  impcut = 5.0
	  if(n_Imp.gt.0)then
	    do i = 1, (n_atoms_lmp - n_Imp)
	      do j = (n_atoms_lmp-n_Imp)+1,n_atoms_lmp
		  if((type_lmp(j).eq.2).or.(type_lmp(j).eq.3))then
	          dx = x_lmp((i-1)*3+1) - x_lmp((j-1)*3+1)
	          dy = x_lmp((i-1)*3+2) - x_lmp((j-1)*3+2)
	          dz = x_lmp((i-1)*3+3) - x_lmp((j-1)*3+3)
	          dmag = sqrt(dx*dx+dy*dy+dz*dz)
		  if(dmag.le.impcut)then
		    type_lmp(i) = 4
		    exit
		  endif
	        endif
	      enddo
	    enddo
	  endif

	  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	  if(rank.eq.0) write(*,*) 'Adjust Step#2 Done.'

          if(Monolayer.and.(.not.explicitOxy)
     &				.and.(iStep.eq.0))then

	   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	   if(rank.eq.0) write(*,*) 'doing Monolayer Stuff...'

	   allocate(surf_list(numnp))
	   surf_list(1:numnp) = 0

           i = 0
	   j = 0
	   s_cnt = 0
	   o_cnt = 0
	if(Adhes)then
           do irep = 1, numnp
            if(IsRelaxed(irep).ne.0)then
              i = i + 1
              if(IsRelaxed(irep).eq.1)then
                if( (dis_12(i).gt.Sbarrier)
     &		.and.((type_lmp(i).eq.1)
     &		.or.(type_lmp(i).eq.11).or.(type_lmp(i).eq.12)) )then
	          s_cnt = s_cnt + 1
	          if((type_lmp(i).eq.11).or.(type_lmp(i).eq.12))then
	            o_cnt = o_cnt + 1
		  else
		    j = j + 1
		    surf_list(j) = i
                  endif
                endif
              endif
            endif
           enddo
	else
	    do irep = 1, numnp
            if(IsRelaxed(irep).ne.0)then
              i = i + 1
              if(IsRelaxed(irep).eq.1)then
                if( (dis_12(i).gt.Sbarrier)
     &		.and.((type_lmp(i).eq.1)
     &		.or.(type_lmp(i).eq.4) ) )then
	          s_cnt = s_cnt + 1
	          if(type_lmp(i).eq.4)then
	            o_cnt = o_cnt + 1
		  else
		    j = j + 1
		    surf_list(j) = i
                  endif
                endif
              endif
            endif
           enddo
	endif !adhes

	  o_need = INT(o2_conc*s_cnt)
	  o_need = max(0,o_need-o_cnt)
	  if(rank.eq.0) write(*,*) 's_cnt = ',s_cnt
	  if(rank.eq.0) write(*,*) 'o_cnt = ',o_cnt
	  if(rank.eq.0) write(*,*) 'o_need = ',o_need

	  if(o_need.gt.0)then

	   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	   if(rank.eq.0) write(*,*) 'adding some Implicit Oxy...'

	   if(rank.eq.0) call RANDOM_SEED
	   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

           do i= 1, o_need
	    keepTrying = .true.
	    o_cnt = 0
	    ! Look for a random surface atom of type-1 to "stiffin"...
		if(rank.eq.0) write(*,*)' lookin 4 atom to stiffen !'
	    do while (keepTrying)
	      if(rank.eq.0) call RANDOM_NUMBER(rand_d)
	      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	      CALL MPI_BCAST(rand_d,1,MPI_DOUBLE,0,
     &					MPI_COMM_WORLD,ierr)
	      rand_i = int(floor(rand_d*j)+1)
	      if((surf_list(rand_i).lt.1).or.
     &			(surf_list(rand_i).gt.n_Al))then
		if(rank.eq.0) write(*,*)' Error in adjust_displ !'
	      endif
              if(dis_12(surf_list(rand_i)).gt.Sbarrier)then
	        if((type_lmp(surf_list(rand_i)).eq.1).or.
     &			(type_lmp(surf_list(rand_i)).eq.13).or.
     &			(type_lmp(surf_list(rand_i)).eq.14))then
		if(rank.eq.0) write(*,*)'we found',surf_list(rand_i), '!'
!	          type_lmp(surf_list(rand_i)) = 4  
		  if(Atom_top(surf_list(rand_i)).eq.1)then	
		    if(rank.eq.0)then
			 write(*,*)'it equals ',
     &                      Atom_top(surf_list(rand_i))
                    endif
		    if(rank.eq.0) write(*,*)'and its on top'
	            type_lmp(surf_list(rand_i)) = 12 
		  elseif(Atom_bot(surf_list(rand_i)).eq.1)then
		    if(rank.eq.0)then
			 write(*,*)'it equals ',
     &                      Atom_bot(surf_list(rand_i))
		    endif	

		    if(rank.eq.0) write(*,*)'and its on bot'
	            type_lmp(surf_list(rand_i)) = 11
		  endif 
		  keepTrying = .false.
                endif
              endif
	      o_cnt = o_cnt + 1
	      if(o_cnt.gt.1000)then
		keepTrying = .false.
		if(rank.eq.0) write(*,*) 'Failed finding surface atom...'
              endif
            enddo
	   enddo


	   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	   if(rank.eq.0) write(*,*) 'Implicit Atoms Set.'

	  endif
	
	!ER- Now label the top and bottom surface atoms
	print *, 'about to set surf Al!'
	if(Adhes)then
	print *, 'setting surf Al!'
		 do irep = 1, n_Al
!		 if(rank.eq.0)then
!			print *,'atom number ',irep
!			print *,'has type ',type_lmp(irep)
!			print *,'topflag ',Atom_top(irep)
!			print *,'botflag ',Atom_bot(irep)
!		 endif
			if(type_lmp(irep).eq.1)then ! it's not a pad or ox
			  if(Atom_top(irep).eq.1)then
				type_lmp(irep)=14
			  elseif(Atom_bot(irep).eq.1)then
				type_lmp(irep)=13
			  endif ! top or bot

			endif !if at first type 1
		  enddo
	endif !Adhes


	  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	  deallocate(surf_list)

	  endif !Monolayer

	  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	  if(rank.eq.0) write(*,*) 'Adjust Step#3 done.'

	  ! Correct for pad atom motion:
	  i=0
	  do irep = 1, numnp
            NeedList=IsRelaxed(irep).ne.0
            if(NeedList)then
	      i=i+1
	      if(IsRelaxed(irep).eq.-1)then
	      x_lmp((i-1)*3+1) = atomDispl(1,irep)+atomCoord(1,irep)
	      x_lmp((i-1)*3+2) = atomDispl(2,irep)+atomCoord(2,irep)
	      x_lmp((i-1)*3+3) = atomDispl(3,irep)+atomCoord(3,irep)
	      endif
	    endif
	  enddo

	  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	  if(rank.eq.0) write(*,*) 'Adjust Step#4 done.'

	  if(iStep.eq.0)then
	  if(rank.eq.0) 
     &	  open(unit=216, file='dump_type', status='replace') 
	  do i = 1, n_atoms_lmp
           if(rank.eq.0) write(216,'(F16.8,1X,F16.8,1X,F16.8,1X,I6)')
     &  x_lmp((i-1)*3+1),x_lmp((i-1)*3+2),x_lmp((i-1)*3+3),type_lmp(i)
	  enddo
	  if(rank.eq.0) close(216)
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	  if(rank.eq.0) write(*,*) 'dump_type written...'

	  endif

        return
	end

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
!FFFF!                                                          !FFFF!
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!
	subroutine find_crack_tip(n_CV,CV_list,x_cracktip,
     &	y_cracktip,EAMdebug,nnsteps,atomCoord,atomDispl,dis_12)
	use mod_global
	implicit none

        double precision, allocatable, INTENT(INOUT) :: dis_12(:)
	double precision cx, cy, x_cracktip, y_cracktip
	double precision x_cracktip_cv, y_cracktip_cv
	double precision xcrack0, ycrack0
        double precision atomCoord(ndf,*),atomDispl(ndf,*)
	integer i,n_CV,nnsteps,irep
        integer, allocatable, intent(inout) :: CV_list(:,:)
	logical EAMdebug

	x_cracktip = 0.01
	y_cracktip = 0.01

          do i = 1, n_CV
	    cx = grid_xmin + (CV_list(i,1)-1) * grid_dx
	    cy = grid_ymin + (CV_list(i,2)-1) * grid_dx
	    if(cx.gt.x_cracktip) then
	      x_cracktip_cv = cx
	      y_cracktip_cv = cy
	    endif
	  enddo

!MZ- I am printing x_crack_cv to file 'output'	
	if(rank.eq.0) then
	    print *,'MINGJIE - x_cracktip_cv = ',x_cracktip_cv
	endif	

	  i = 0
	  do irep = 1, numnp
	    if(IsRelaxed(irep).ne.0)then
	      i = i + 1
              if(IsRelaxed(irep).eq.1)then
	      if(dis_12(i).gt.Sbarrier)then
	        cx = atomDispl(1,irep)+atomCoord(1,irep)
	        cy = atomDispl(2,irep)+atomCoord(2,irep)
	        if(cx-monoThick.gt.x_cracktip) then
	          x_cracktip = cx !- monoThick
	          y_cracktip = cy
	        endif
	      endif
	      endif
	    endif
	  enddo

!MZ- I am printing x_cracktip to file 'output'	
	if(rank.eq.0) then
	    print *,'MINGJIE - x_cracktip = ',x_cracktip
	endif	


	if((x_cracktip-x_cracktip_cv).gt.(4.d0))then
	  x_cracktip = x_cracktip_cv
	  y_cracktip = y_cracktip_cv
	endif

        ! READ "cracktip" file: 
        if(rank.eq.0) then
          open(unit=2117, file='cracktip', status='unknown')
	  read(2117,*) xcrack0, ycrack0
	  close(2117)
	endif
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	!if(x_cracktip.gt.xcrack0)then
	if(x_cracktip.gt.(0.d0))then
	  xcrack0=x_cracktip
	  ycrack0=y_cracktip
	endif

        ! WRITE "cracktip" file: 
        if(rank.eq.0) then
          open(unit=2117, file='cracktip', status='replace')
	  write(2117,*) xcrack0, ycrack0
	  close(2117)
!	  write(*,*)''
!	  write(*,*)'New Cracktip x-coord: ',xcrack0
!	  write(*,*)''
	endif
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	x_cracktip = xcrack0
	y_cracktip = ycrack0

	closureDist = closureDist_o + x_cracktip

        return
	end
