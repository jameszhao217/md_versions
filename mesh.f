C Remove 3 layers of atoms around crack faces, for nickel study. x=[-2 -1 1], y=[1 -1 1]
C Define detection band rings around crack tip
C mesh generator for atomically sharp crack tip
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
c build a blunt *center, crack
      subroutine mp01(id,x,ix,f,b,itx)
      use mod_grain
      use mod_global
      use mod_file
      use mod_boundary
      use mod_crack
      implicit none

!     Input variables
      integer id,ix,itx(3,*)
      double precision x,f,b
      dimension id(ndf,1),x(nxdm,1),ix(nen1,1),f(ndf,1),b(ndf,1)
      integer nodalmap(1,maxnp)

!     common declarations
      common/debugger/debug
      logical debug
      double precision tol, stopfine

      double precision, allocatable :: IsRelaxed_TMP(:)
      double precision, allocatable :: x_TMP(:,:)
      double precision, allocatable :: f_TMP(:,:)

!     local variables
      logical n1,n2,n3,m1,m2,m3, useDetectionBand
      double precision, pointer:: nodeAngle(:)
      integer nXRegions,nYRegions, icell,ndxdy,numNodes,numx,numy,
     $     iSpace,ixrem,iyrem,iGrain,coincidentNodes,nodeStart,inode,
     $     nsort,np1,numnp0, countLast, i,j, k, node1, node2, node3
      integer nr1,nr2,nr3,nr4
      double precision XMax(0:20),YMax(0:20),nodeSite(3),
     $     dx,dy,dxdy,xxmin,xxmax,yymin,yymax,yyMinOrig,
     $     delx,dely,xx,yy, minDb, rcutmesh,Dist2, large
      double precision minDbx,minDby
      double precision ro_x,ro_y,rdist,minfoundx
      data tol /1.d-6/
      logical placeNode,top,bot,left,right,mirror
      
      type Region
      double precision xmin, xmax, ymin, ymax
      end type Region
      TYPE(Region) atomRegion, detectionBand, innerRegion, 
     $	mirrorAtomRegion, simulationCell
      logical insideRegion,inside,crack111   
      
C! VBS added this to read in tolerances
      common/elemconv/numelold
      integer lcrack,numelold,numnpc
      integer dwnumx,dwnumy,tst,minfoundi,bndry
      double precision dwfactorx, dwfactory

      dwfactorx=2.0
      dwfactory=1.0


      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

      if(rank.eq.0)then
      write(*,*)
      write(*,*) 'Generating mesh containing an embedded crack'
      write(*,*)        
      endif

	firstIMP=1

      x0crack = 0.01 
      y0crack = 0.01
!
      XMax(0)=-1.
      YMax(0)=-1.
      
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
! Read mesh data            
      if(rank.eq.0)then
       read(5,*) nXRegions,nYRegions
       read(5,*) (XMax(i),i=1,nXRegions)
       read(5,*) (YMax(i),i=1,nYRegions)
       read(5,*) minDb
       read(5,*) rcutmesh
       read(5,*) x0crack, y0crack, r_crack
       read(5,*) closureMag, oxyThick, closureDist
       read(5,*) dwfactorx, dwfactory
       read(5,*) grid_dx
       read(5,*) sumReq
       read(5,*) grid_xmin,grid_xmax
       read(5,*) grid_ymin,grid_ymax
	crack111 = .true.

	minDbx = XMax(1)-rcutmesh
	minDby = YMax(1)-rcutmesh

      endif

      !if(rank.eq.0) write(*,*) 'Starting mpi broadcast...' 
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nXRegions,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nYRegions,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(minDb,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(minDbx,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(minDby,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rcutmesh,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(x0crack,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(y0crack,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(r_crack,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(closureMag,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(oxyThick,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      oxyThick_max = oxyThick
      oxyThick = 0.d0
      CALL MPI_BCAST(closureDist,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dwfactorx,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dwfactory,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

	closureDist_o = closureDist

      CALL MPI_BCAST(grid_dx,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sumReq,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(grid_xmin,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(grid_xmax,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(grid_ymin,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(grid_ymax,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(crack111,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

      do k=1,nXRegions
      CALL MPI_BCAST(XMax(k),1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      enddo
      do k=1,nYRegions
      CALL MPI_BCAST(YMax(k),1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      enddo
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      !if(rank.eq.0) write(*,*) 'mesh mpi broadcast done...' 

      mirror=.false.      

      do k=min(nXRegions,nYRegions)+1,max(nXRegions,nYRegions)
         if (nXRegions.lt.nYRegions) then
            XMax(k)=XMax(nXRegions)
         else
            YMax(k)=YMax(nYRegions)
         endif
      enddo

        ! write "cracktip" file: 
        if(rank.eq.0) then
          open(unit=2117, file='cracktip', status='replace')
	  write(2117,*) x0crack, y0crack
	  close(2117)
	endif
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

c**---------------------------------------------------------------------
!RZ!  MAKE Adaptive CV Grid:
c**---------------------------------------------------------------------

      nx_grid = INT((grid_xmax-grid_xmin)/grid_dx)
      ny_grid = INT((grid_ymax-grid_ymin)/grid_dx)
      n_CV_cells = nx_grid*ny_grid
      if(rank.eq.0)write(*,*)'nx_grid=',nx_grid
      if(rank.eq.0)write(*,*)'ny_grid=',ny_grid

c**---------------------------------------------------------------------
c** Define the cracktip origin and R-value to use...
c**---------------------------------------------------------------------
      ro_x = x0crack+0.01
      ro_y = y0crack
      ro_x = ro_x - r_crack

!     Extract lattice data
!     Normally, the second atom is dx away from the first due to the
!     way cell is sorted.
!     if there is only 1 atom on the lowest y-plane of cell, then dx is the
!     cell width in the x dirn.
      dx=grains(1)%cell(1,2)
      if (grains(1)%cell(2,2).gt.tol) dx=grains(1)%dcell(1)

      do icell=1,grains(1)%ncell
         if (grains(1)%cell(2,icell).gt.tol) then
            dy=grains(1)%cell(2,icell)
            dxdy=grains(1)%cell(1,icell)
            goto 9
         endif
      enddo
 9    continue
 
      if (dxdy.ne.0.) then
         ndxdy=nint(dx/dxdy)
      else
         ndxdy=1
      endif
      
      if(rank.eq.0) print*,'Generating nodes'
! Generate the nodes in each box, the mesh is coarsened with increasing
! distance from the center of the box.
      numx=int(XMax(nXRegions)/dx)+1
      numy=int(YMax(nYRegions)/dy)+1
      numNodes=0
      do 10 i=-numx,numx
      
         xx=i*dx
         if (abs(xx).gt.XMax(nXRegions)) go to 10
	 
         do 20 j=-numy,numy
	 
            yy=j*dy
            if (abs(yy).gt.YMax(nYRegions)) go to 20
	    
!! Determine the region in which the node is being placed	    
            do k=max(nXRegions,nYRegions)-1,1,-1
!               if ( (abs(xx).gt.XMax(k))
!     $              .or.(abs(yy).gt.YMax(k)) ) go to 21

		if (k.eq.1)then
		  stopfine = -XMax(k)/dwfactorx
		else
		  stopfine = -XMax(k)
		endif

                if (     (xx.gt.XMax(k))
     $              .or. (yy.gt.YMax(k))
     $              .or. (xx.lt.(stopfine))
     $              .or. (yy.lt.(-YMax(k)/dwfactory))
     $              ) go to 21

!                if ( (xx.gt.XMax(k))
!     $              .or.(yy.gt.YMax(k)) ) go to 21   !
!                if ( (xx.lt.-XMax(k)/dwfactor)
!     $              .or.(yy.gt.YMax(k)) ) go to 21   !
!                if ( (xx.gt.XMax(k))
!     $              .or.(yy.lt.-YMax(k)/dwfactor) ) go to 21  !
!               if ( (xx.lt.-XMax(k)/dwfactor)  
!     $              .or.(yy.lt.-YMax(k)/dwfactor) ) go to 21  !

            enddo
 21         iSpace=2**k
            if(iSpace.eq.0) iSpace=1

!! Decide if a node should be placed in this region
            ixrem=mod(abs(i),iSpace)
            iyrem=mod(abs(j),iSpace)
            placeNode=(ixrem+iyrem).eq.0
            if (.not.placeNode) go to 20

! Assign the node a position
            numNodes=numNodes+1
            x(1,numNodes)=xx+mod(j,ndxdy)*dxdy
            x(2,numNodes)=yy


c**---------------------------------------------------------------------
c** Remove node if it is in crack-region...
c**---------------------------------------------------------------------
c! Qu modification to remove extra layers of atoms around crack faces
c            if((j.eq.0.and.x(1,numNodes).lt.-tol).or.
c     &      (j.eq.1.and.x(1,numNodes).lt.-0.5*dx*2).or.
c  !kris added *2 to match current work, otherwise original may crack!
c     &      (j.eq.-1.and.x(1,numNodes).lt.-dx) )then
c              numNodes=numNodes-1
c            endif

	    rdist = sqrt( (x(1,numNodes)-ro_x)**2 + 
     &			(x(2,numNodes)-ro_y)**2 )

            if(    (rdist.lt.r_crack)
     &	       .or.( (x(1,numNodes).lt.ro_x)
     &	      .and.( abs( x(2,numNodes)-ro_y ).lt.r_crack ) ) )then
              numNodes=numNodes-1
            endif

! FOR NOW ONLY, assign nodes to the crack faces. Leo mentions
! this is a hack. But this appears legitimate.
            !if(j.eq.0.and.x(1,numNodes+1).lt.ro_x) then
            if((j.eq.0).and.
     &		(x(1,numNodes+1).lt.-XMax(1)/dwfactorx)) then

	      ! What is the number of steps to go in each dir?
	      do tst = 1,100
		if ((yy+tst*dy).ge.(r_crack-ro_y))then
		  exit
		endif
	      enddo

              numNodes = numNodes+2
              x(1,numNodes-1) = xx+mod(j,ndxdy)*dxdy
              x(2,numNodes-1) = yy+tst*dy
              !x(2,numNodes-1) = r_crack - ro_y 

              x(1,numNodes) = xx+mod(j,ndxdy)*dxdy
              x(2,numNodes) = yy-tst*dy 
              !x(2,numNodes) = ro_y - r_crack

!              print *, i, x(1, numNodes)
              !if(i.eq.-320) then
                countLast = numNodes-1
              !endif
            endif
C Qu modification ends

 20      continue
 10   continue
 
 
	if(rank.eq.0) print*, 'Moving nodes to the nearest atomic sites'
! Move nodes to the nearest atomic sites
        large = 1.e30
	xxmax=-large
	xxmin=large
	yymax=-large
	yymin=large	
	do i=1,numNodes
	
		nodeSite(1)=x(1,i)
		nodeSite(2)=x(2,i)
		if(i.ne.countLast) then
		call NearestBsite(nodeSite,1,.false.,x(1,i),iGrain)
		endif

!! find xxmax, xxmin, etc.
		xxmax = max(x(1,i), xxmax)
		xxmin = min(x(1,i), xxmin)
		yymax = max(x(2,i), yymax)
		yymin = min(x(2,i), yymin)
	enddo

      simulationCell%xmin = xxmin
      simulationCell%xmax = xxmax
      simulationCell%ymin = yymin
      simulationCell%ymax = yymax
            
      if(mirror) then
         numnp=numNodes
         do i=1,numnp
	 
            if(x(2,i).lt.yymin+tol) then
               nodeSite(1)=x(1,i)
               nodeSite(2)=-136.85
               call NearestBsite(nodeSite,1,.false.,x(1,i),iGrain)
            endif
	    
            numNodes=numNodes+1
            nodeSite(1)=x(1,i)
            nodeSite(2)=2*yymin-x(2,i)
            call NearestBsite(nodeSite,1,.false.,x(1,numNodes),iGrain)
	    
         enddo
	 
         yyMinOrig=yymin
         yymin=2*yymin-yymax
	 
	 simulationCell%ymin = yymin
      endif

! Remove coincident nodes
      if(rank.eq.0) print*, 'Removing coincident nodes'
      coincidentNodes=0
      numnp=numNodes
      do 30 i=numnp,2,-1
         if (i.gt.numNodes) goto 30
         do j=1,i-1
            delx=abs(x(1,i)-x(1,j))
            dely=abs(x(2,i)-x(2,j))
            if (delx+dely.lt.2.*tol) then
               x(1,j)=x(1,numNodes)
               x(2,j)=x(2,numNodes)
               numNodes=numNodes-1
               coincidentNodes=coincidentNodes+1
            endif
         enddo
 30   enddo
 
      if (coincidentNodes.ne.0) then
        if(rank.eq.0)then
	write(6,*) coincidentNodes,' coincident nodes removed'
        endif
      endif
      
      numnp=numNodes
      if (numnp.gt.maxnp) then
         write(6,*)'***ERROR: Insufficient storage for nodes'
         write(6,*) '         numnp = ',numnp
         stop
      endif
      if(rank.eq.0)then
       write(6,*) 'Total nodes: numnp = ',numnp
      endif

       !write(*,*) ' Processor ',rank,' is done..'
       !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
       !call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
       !call MPI_FINALIZE(ierr)
       !stop

!     Apply boundary conditions and define the boundary for the
!     triangulator

      nce=0
      nodeStart=0

      if(rank.eq.0) print*, 'Detecting the lower crack ledge'
!     FIND THE LOWER CRACK LEDGE
      do i =1, numnp
	!if(dabs(x(2,i)+tst*dy).lt.tol.and.x(1,i).lt.-dx) then
	if(dabs(x(2,i)+tst*dy).lt.tol.and.x(1,i).le.ro_x+1.1*dx) then
        nce=nce+1
        if (nce.gt.NCEMAX) then
           if(nce.gt.NCEMAX) call IncreaseElist(100)
        endif
        elist(1,nce)=i
        endif
      enddo

! sort the lower ledge so that is goes CW (from right to left).
      allocate(nodeAngle(numnp))
      nodeAngle=0.
      do i=1,nce
         inode=elist(1,i)
         nodeAngle(inode)= datan2(dy,x(1,inode))
      enddo
      nsort=nce
      call qsortr(nsort,elist(1,1),nodeAngle,1,1,1.d0)
      deallocate(nodeAngle)
      
! remove the last node from the list
Cc!!! Bill's changes!!!!
      nce=nce-1
      nodeStart = nce
      
!     find all external boundary nodes
      simulationCell%xmin = simulationCell%xmin + 10.d0
      simulationCell%xmax = simulationCell%xmax - 10.d0
      simulationCell%ymin = simulationCell%ymin + 10.d0
      simulationCell%ymax = simulationCell%ymax - 10.d0

      if(rank.eq.0) print*, 'Detecting the outer cell boundary'
      do 77 i=1,numnp

         top=(x(2,i) .gt. simulationCell%ymax)
         bot=(x(2,i) .lt. simulationCell%ymin)
         left=(x(1,i) .lt. simulationCell%xmin)
         right=(x(1,i) .gt. simulationCell%xmax)

!     store all boundary points, but put crack faces at the beginning of
!     elist.  While you are at it, apply the b.c.s

         if (top.or.bot.or.right.or.left) then
            nce=nce+1
            if (nce.gt.NCEMAX) then
               if(nce.gt.NCEMAX) call IncreaseElist(100)
            endif
            elist(1,nce)=i

! apply the b.c's
            if (top.or.bot.or.right.or.left) then
               id(1,i)=1
               id(2,i)=1
		if(abs(x(2,i)).lt.2000.d0)then
		  id(1,i)=0
                  id(2,i)=0
		endif
            endif
	    	    
         endif
 77   continue

! sort the boundary so that is goes CW.
      allocate(nodeAngle(numnp))
      nodeAngle=0.
      do i=nodeStart+1,nce
         inode=elist(1,i)
! YET ANOTHER HACK
         nodeAngle(inode)=datan2(x(2,inode)-tol,x(1,inode))
      enddo
      nsort=nce-nodeStart
      call qsortr(nsort,elist(1,nodeStart+1),nodeAngle,1,1,1.d0)
      deallocate(nodeAngle)

! remove the last node from the list
      nce=nce-1
      nodeStart = nce


!     FIND THE UPPER CRACK LEDGE
      if(rank.eq.0) print*, 'Detecting the upper crack ledge'
      do i =1, numnp
! Qu modification to remove extra layers of atoms around crack faces
        !if(dabs(x(2,i)-tst*dy).lt.tol.and.x(1,i).lt.tol) then
        if(dabs(x(2,i)-tst*dy).lt.tol.and.x(1,i).le.ro_x+1.1*dx) then
! Qu modification ends
!        if(dabs(x(2,i)-dy).lt.tol.and.x(1,i).lt.XMax(1)-4*rcutmesh) then
!        if(dabs(x(2,i)-dy).lt.tol.and.x(1,i).lt.-XMax(1)+2*dx) then
        nce=nce+1
        if (nce.gt.NCEMAX) then
           if(nce.gt.NCEMAX) call IncreaseElist(100)
        endif
        elist(1,nce)=i
        endif
      enddo

! sort the upper ledge so that is goes CW (from left to right).
      allocate(nodeAngle(numnp))
      nodeAngle=0.
      do i=nodeStart+1,nce
         inode=elist(1,i)
         nodeAngle(inode)= -datan2(x(2,inode),x(1,inode))
      enddo
      nsort=nce-nodeStart
      call qsortr(nsort,elist(1,nodeStart+1),nodeAngle,1,1,1.d0)
      deallocate(nodeAngle)

c**---------------------------------------------------------------------
c** Finish defining the boundary
c**---------------------------------------------------------------------
      do i=1,nce-1
         elist(2,i)=elist(1,i+1)
      enddo
      elist(2,nce)=elist(1,1)
      ncb=nce

      if(rank.eq.0) write(*,*) 'tst = ',tst

      do j = (tst-1),-(tst-1),-1

	minfoundx = 10000.d0
	minfoundi = 1
        do i=1,numnp
          if( dabs(x(2,i)-j*dy).lt.tol )then
	     if(x(1,i).lt.minfoundx)then
		minfoundx = x(1,i)
		minfoundi = i
	     endif
          endif
        enddo

        nce=nce+1
        if(nce.gt.NCEMAX) then
           if(nce.gt.NCEMAX) call IncreaseElist(100)
        endif

        elist(1,nce)=minfoundi
        ncb=nce
        elist(2,nce)=elist(1,1)
        elist(2,nce-1)=elist(1,nce)

      enddo

c     Triangulate, sets all elements to material 1 for this mesh
      if(rank.eq.0) print*, 'Triangulating'
      numnpc = numnp
      call delaunay(id,x,ix,f,b,itx)
      if(rank.eq.0) print*, 'Done'

      if(rank.eq.0)then
       write(*,*) 'BEFORE adding overlap'
       write(6,*) 'Number of nodes: numnp = ',numnp
       write(6,*) 'Number of elements: numel = ',numel
      endif
      if(numel.gt.maxel) stop 'too many elements'
      if(numnp.gt.maxnp) stop 'too many nodes'

!     Find max/min of atomistic region
	innerRegion%xmin = -Xmax(1)/dwfactorx
	innerRegion%xmax = Xmax(1)
	innerRegion%ymin = -Ymax(1)/dwfactory   ! 
	innerRegion%ymax = Ymax(1)   !
	call findAtomRegionSize(x,dx,dy,tol,innerRegion,atomRegion)

       do 75 i = 1,numel
         ix(nen1,i) = 1
 75   continue

!     Find continuum region elements
      if (mirror) then
		mirrorAtomRegion%xmin = atomRegion%xmin
		mirrorAtomRegion%xmax = atomRegion%xmax
		mirrorAtomRegion%ymin = 2*yyMinOrig-atomRegion%ymin
		mirrorAtomRegion%ymax = 2*yyMinOrig-atomRegion%ymax
      endif
	
      
!	Check if each node of any continuum element is in the
!	atomistic region      
      do i=1,numel

      	node1 = ix(1,i)
	node2 = ix(2,i)
	node3 = ix(3,i)
	
!! Determine if any node is in the atomistic region      
         n1=.not.insideRegion(x(1:2,node1),atomRegion)
         n2=.not.insideRegion(x(1:2,node2),atomRegion)
         n3=.not.insideRegion(x(1:2,node3),atomRegion)

!! Determine if any node is in the mirrored atomistic region
         if(mirror) then	 	     
             m1=.not.insideRegion(x(1:2,node1),mirrorAtomRegion)
             m2=.not.insideRegion(x(1:2,node2),mirrorAtomRegion)
             m3=.not.insideRegion(x(1:2,node3),mirrorAtomRegion)	     
         else
            m1=.true.
            m2=.true.
            m3=.true.
         endif
	 
         if((n1.and.n2.and.n3).and.(m1.and.m2.and.m3)) ix(nen1,i)=0
      enddo
           
!     Add pad atoms in the interface region. Note that elements are 
!     not needed in this region, so just add atoms.
!      go to 1233
!      go to 1234
      if(rank.eq.0) print*, 'Adding pad atoms'
      numNodes=numnp
!      XMax(2)=XMax(1)+2.d0*rcutmesh
!      YMax(2)=YMax(1)+2.d0*rcutmesh
!rz!      XMax(2)=XMax(1)+6.d0
!rz!      YMax(2)=YMax(1)+6.d0
!rz!      numx=int(XMax(2)/dx)+1
!rz!      numy=int(YMax(2)/dy)+1
!rz!      dwnumx=int((XMax(1)+6.d0)/dwfactor/dx)+1
!rz!      dwnumy=int((YMax(1)+6.d0)/dwfactor/dy)+1

      XMax(2)=XMax(1)+10.d0
      YMax(2)=YMax(1)+10.d0
      numx=int(XMax(2)/dx)+1
      numy=int(YMax(2)/dy)+1
      dwnumx=int(((XMax(1)/dwfactorx)+10.d0)/dx)+1
      dwnumy=int(((YMax(1)/dwfactory)+10.d0)/dy)+1

      do i=-dwnumx,numx
         xx=i*dx
         do j=-dwnumy,numy !was numy,dwnumy
            yy=j*dy
	    
            numNodes=numNodes+1
	    
            nodeSite(1)=xx+mod(j,ndxdy)*dxdy
            nodeSite(2)=yy
	    
            call NearestBsite(nodeSite,1,.false.,x(1,numNodes),iGrain)

!!          Skip this node if it is the atomistic region
            if(abs(j).lt.(tst).and.i.lt.0)then
              numNodes=numNodes-1
            elseif(insideRegion(x(1:2,numNodes),atomRegion)) then
               numNodes=numNodes-1
               endif
	    
            if (mirror) then
               nodeSite(1)=x(1,numNodes)
               nodeSite(2)=2*yyMinOrig-x(2,numNodes)
               numNodes=numNodes+1
               call NearestBsite(nodeSite,1,.false.,x(1,numNodes),
     $				iGrain)
            endif	    
         enddo
      enddo
      if(rank.eq.0) print*, 'Done adding pad atoms'      
            
! 1233 if(rank.eq.0) write(*,*) 'Forget the pad atoms...'
      np1=numnp+1
      numnp0=numnp
      numnp=numNodes
      if(numel.gt.maxel) stop 'Too many elements'
      if(numnp.gt.maxnp) stop 'Too many nodes' 
      
!     Determine nodal character -- continuum, interface, atomistic etc.
      num2Dnode=numnp
      call StatusCalc(x,ix,.true.)

!     remove the coincident nodes on the interface
      if(rank.eq.0) then
	print*, 'Removing coincident nodes near the interface'
      endif
      do i=numnp,np1,-1
         do j=1,numnp0
            if(IsRelaxed(j).ne.0) then
               if(Dist2(x(1:2,i),x(1:2,j),2).lt.1.e-6) then
                  f(1:ndf,i)=f(1:ndf,numnp)
                  x(1:nxdm,i)=x(1:nxdm,numnp)
                  numnp=numnp-1
                  go to 300
               endif
            endif
         enddo
 300     continue
      enddo
      num2Dnode=numnp
      if(numperiodz.gt.1)then
         call IncAtoms(x)
      endif

      if(numnp.gt.maxnp) stop 'Too many nodes'
      allocate(atomSpecie(numnp))
      do i=1,numnp
         atomSpecie(i)=1
      enddo
      if(rank.eq.0) print*, 'done'

 1234 if(rank.eq.0) write(*,*) 'Final mesh size'
      if(rank.eq.0) write(6,*) 'Total nodes: numnp = ',numnp
      if(rank.eq.0) write(6,*) 'Total elements: numel = ',numel
      if(rank.eq.0) write(6,*) 'rcutmesh = ', rcutmesh
      numnpp1=-1

      allocate(id_closure(numnp))
      id_closure(1:numnp) = 0
      allocate(nodeArea(numnp))
      nodeArea(1:numnp) = 1.d0
      do i=1,nce
	if((elist(1,i).gt.numnp).or.(elist(1,i).lt.1))then
	  if(rank.eq.0) write(*,*) ' ERROR IN MESH.F (A) !'
	endif
	if((isrelaxed(elist(1,i)).eq.0)
     &	.and.(abs(x(2,elist(1,i))).lt.(2.d0*r_crack))
     &	.and.(x(1,elist(1,i)).lt.0.d0) !-closureDist)
     &	)then
	  id(1,elist(1,i)) = 0
	  id(2,elist(1,i)) = 0
	  id_closure(elist(1,i)) = 1
	  do j = 1,nce
	    if(elist(2,j).eq.elist(1,i))then
	     nodeArea(elist(1,i)) = 
     &		abs(x(1,elist(1,i))-x(1,elist(2,i)))/2.d0
     &		+ abs(x(1,elist(1,i))-x(1,elist(1,j)))/2.d0
	     exit
	    endif
	  enddo
	endif
      enddo

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!     Create a detection band
!     Identify a path, defined by a closed polygon, ccw around the
!     vertices, along which the detection band elements will be placed.
!     If the atomistic continuum interface is not a closed path, define
!     this polygon to extend outside the mesh and surround the atomistic
!     region.

	useDetectionBand = .true.
	!useDetectionBand = .false.
	if (.not. useDetectionBand) then
		ndbpoly = 0
		return
	endif

	if(rank.eq.0) print*, 'Using a detection band'
!     There is one detection band with 4 points
! Qu modified detection band rings starts
	detectionBand%xmin = atomRegion%xmin + (rcutmesh)
	detectionBand%xmax = atomRegion%xmax - (rcutmesh)
	detectionBand%ymin = atomRegion%ymin + (rcutmesh)
	detectionBand%ymax = atomRegion%ymax - (rcutmesh)
        minDby=int(minDby/dy)*dy
        minDbx=int(minDbx/dx)*dx
	nr1=int((abs(detectionBand%xmin)-minDbx)/dx)+1
	nr2=int((abs(detectionBand%xmAX)-minDbx)/dx)+1
	nr3=int((abs(detectionBand%ymin)-minDby)/dy)+1
	nr4=int((abs(detectionBand%xmax)-minDby)/dy)+1
	ndbpoly=max(nr1,nr2,nr3,nr4)+1
	allocate(ndbvtx(ndbpoly)) 
	allocate(dbpoly(2,4,ndbpoly))
	ndbvtx(1:ndbpoly)=4

      do i=1,ndbpoly
         detectionBand%xmin = atomRegion%xmin + rcutmesh
         detectionBand%xmax = min( minDbx+(i-1)*dx
     $                            ,atomRegion%xmax - rcutmesh)
         detectionBand%ymin = max(-minDby-(i-1)*dy
     $                            ,atomRegion%ymin + rcutmesh)
         detectionBand%ymax = min( minDby+(i-1)*dy
     $                            ,atomRegion%ymax - rcutmesh)
      
   
         dbpoly(1,1,i)=detectionBand%xmin
         dbpoly(2,1,i)=detectionBand%ymin
      
         dbpoly(1,2,i)=detectionBand%xmax
         dbpoly(2,2,i)=detectionBand%ymin
      
         dbpoly(1,3,i)=detectionBand%xmax
         dbpoly(2,3,i)=detectionBand%ymax
      
         dbpoly(1,4,i)=detectionBand%xmin
         dbpoly(2,4,i)=detectionBand%ymax

	if(rank.eq.0)then
	  write(6,*) ''
	  write(6,*) 'dbpoly(1,1,',i,') = ',dbpoly(1,1,i)
	  write(6,*) 'dbpoly(2,1,',i,') = ',dbpoly(2,1,i)
	  write(6,*) ''
	  write(6,*) 'dbpoly(1,2,',i,') = ',dbpoly(1,2,i)
	  write(6,*) 'dbpoly(2,2,',i,') = ',dbpoly(2,2,i)
	  write(6,*) ''
	  write(6,*) 'dbpoly(1,3,',i,') = ',dbpoly(1,3,i)
	  write(6,*) 'dbpoly(2,3,',i,') = ',dbpoly(2,3,i)
	  write(6,*) ''
	  write(6,*) 'dbpoly(1,4,',i,') = ',dbpoly(1,4,i)
	  write(6,*) 'dbpoly(2,4,',i,') = ',dbpoly(2,4,i)
	  write(6,*) ''
	endif

      end do
! Qu modified detection band rings ends
      
      if(rank.eq.0) write(6,*) 
     $ 'width of the detection band: rcutmesh = ', rcutmesh

      return
      end
      
      

	 
************************************************************************       
!	Checks to see if a 2D point x is located inside thisRegion.    
      logical function insideRegion(x, thisRegion)
      implicit none

      type Region
      	double precision xmin, xmax, ymin, ymax
      end type Region
      TYPE(Region) thisRegion
      double precision x(2)
      
      insideRegion = (x(1).gt. thisRegion%xmin .and.
     $			x(1) .lt. thisRegion%xmax .and.
     $			x(2) .gt. thisRegion%ymin .and.
     $			x(2) .lt. thisRegion%ymax)
      
      return
      end
      
************************************************************************
      subroutine qsortr(n,list,xkey,nxdm,ind,sign)
      implicit double precision (a-h,o-z)
c
      dimension xkey(nxdm,1)
      integer list(2,*),n,ll,lr,lm,nl,nr,ltemp,stktop,maxstk
c
      parameter (maxstk=32)
c
      integer lstack(maxstk), rstack(maxstk)
c
      ll = 1
      lr = n
      stktop = 0
   10 if (ll.lt.lr) then
         nl = ll
         nr = lr
         lm = (ll+lr) / 2
         guess = sign*xkey(ind,list(1,lm))
c
c     Find xkeys for exchange
c
   20    if (sign*xkey(ind,list(1,nl)).lt.guess) then
            nl = nl + 1
            go to 20
         endif
   30    if (guess.lt.sign*xkey(ind,list(1,nr))) then
            nr = nr - 1
            go to 30
         endif
         if (nl.lt.(nr-1)) then
            ltemp    = list(1,nl)
            list(1,nl) = list(1,nr)
            list(1,nr) = ltemp
            nl = nl + 1
            nr = nr - 1
            go to 20
         endif
c
c     Deal with crossing of pointers
c
         if (nl.le.nr) then
            if (nl.lt.nr) then
               ltemp    = list(1,nl)
               list(1,nl) = list(1,nr)
               list(1,nr) = ltemp
            endif
            nl = nl + 1
            nr = nr - 1
         endif
c
c     Select sub-list to be processed next
c
         stktop = stktop + 1
         if (nr.lt.lm) then
            lstack(stktop) = nl
            rstack(stktop) = lr
            lr = nr
         else
            lstack(stktop) = ll
            rstack(stktop) = nr
            ll = nl
         endif
         go to 10
      endif
c
c     Process any stacked sub-lists
c
      if (stktop.ne.0) then
         ll = lstack(stktop)
         lr = rstack(stktop)
         stktop = stktop - 1
         go to 10
      endif
c
      return
      end

      
      
!     Find max/min of atomistic region      
      subroutine findAtomRegionSize(atomCoord,dx, dy, tol, 
     $		innerRegion, atomRegion)
      use mod_global
      implicit none
      
      type Region
      	double precision xmin, xmax, ymin, ymax
      end type Region
      TYPE(Region) atomRegion, innerRegion      
      double precision atomCoord(ndf,*), dx, dy, tol

!	local variables      
      double precision xmin, xmax, ymin, ymax, x, y, large  
      integer iNode           
!	functions      
      logical insideRegion  
      
      large = 1.e30
            
      xmax=-large
      xmin= large
      
      ymax=-large
      ymin= large
      
      do iNode=1,numnp
      
         x=atomCoord(1,iNode)
         y=atomCoord(2,iNode)
	 
         if(insideRegion(atomCoord(1:2,iNode),innerRegion)) then	    
            xmax=max(xmax,x)
            ymax=max(ymax,y)
            xmin=min(xmin,x)
            ymin=min(ymin,y)	    
         endif
	 
      enddo
      
      xmax=xmax-dx+tol
      ymax=ymax+tol
      xmin=xmin+dx-tol
      ymin=ymin-tol
      
      xmax=xmax-dx
      xmin=xmin+dx
      ymax=ymax-dy
      ymin=ymin+dy
            
      atomRegion%xmin = xmin
      atomRegion%xmax = xmax
      atomRegion%ymin = ymin
      atomRegion%ymax = ymax
            
      return
      end
      
! Qu added begin
!
      subroutine IncAtoms(x)
      use mod_global
      use mod_grain
      implicit none

      double precision x
      dimension x(nxdm,1)

      double precision, allocatable :: x_PAD_TMP(:,:)

      integer i,j,numnp0,minPadID,num_pad

      minPadID=numnp
      do i=1,numnp
        if(IsRelaxed(i).eq.-1)then
	 if(i.lt.minPadID) minPadID=i
	endif
      enddo

      num_pad = numnp-minPadID+1
      allocate(x_PAD_TMP(3,num_pad))
      j=0
      do i=minPadID,numnp
	 j=j+1
	 x_PAD_TMP(1:3,j)=x(1:3,i)
      enddo

      ! Copy over non-pad nodes/atoms
      numnp0=numnp
      do i=2,numperiodz
         do j=1,numnp0
            if(IsRelaxed(j).eq.-1)then
               numnp=numnp+1
               IsRelaxed(numnp)=IsRelaxed(j)
               x(1:2,numnp)=x(1:2,j)
               x(3,numnp)=x(3,j)+(i-1)*grains(1)%dcell(3)
            endif
         enddo
       enddo
       numnp0=numnp
       do i=2,numperiodz
         do j=1,numnp0
            if((IsRelaxed(j).ne.0).and.(IsRelaxed(j).ne.-1))then
               numnp=numnp+1
               IsRelaxed(numnp)=IsRelaxed(j)
               x(1:2,numnp)=x(1:2,j)
               x(3,numnp)=x(3,j)+(i-1)*grains(1)%dcell(3)
            endif
         enddo
       enddo

      deallocate(x_PAD_TMP)

      return
      end     
!
! Qu added ends  
      
************************************************************************
      subroutine PerturbMesh(x,nxdm,numnp,perturb)
      end
      
      
! ************************************************************************
!       logical function inside(x,x1,x2,y1,y2)
!       implicit none
!       double precision x(2),x1,x2,y1,y2
!       inside=(x(1).gt.x1).and.(x(1).lt.x2).and.(x(2).gt.y1).and.(x(2).lt
!      $     .y2)
!       end
!             
