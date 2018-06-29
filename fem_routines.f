c
c $Id: fem_routines.f,v 1.1.1.1 2003/03/12 20:09:00 shastry Exp
	subroutine fem_setup
     &    (numnp, numel, x, id, is_relaxed, ix, itx, cc,
     &     inverse_map, inverse_el_map)
c
c123456789012345678901234567890123456789012345678901234567890123456789012
c         1         2         3         4         5         6         7

	implicit none
	include 'fem_parameters.par'

	integer numnp, numel
	double precision x(3,*)
	integer  id(3,*), is_relaxed(*), ix(4,*), itx(3,*)
	double precision cc(6,6)
	integer inverse_map(*), inverse_el_map(*)

	integer fe_locate

	integer i, j, ij, k, ier
	character*80 error_message

	data i_flag /0/

	if(i_flag.eq.0) then
	  i_flag = 1
	else
	  error_message = 'Second call to fem_setup'
	  call error_handler(error_message)
	endif

	call fe_elastic(cc)

	nnodes=0
	nfixed=0
	npad = 0
	do i=1,numnp
	  if((is_relaxed(i).eq.2).or.(is_relaxed(i).eq.0)) then
	    nnodes = nnodes+1
	    if(nnodes.gt.maxgeo) then
	      error_message =
     &          'maxgeo needs to be increased (is_relaxed=(2/0))'
	      call error_handler(error_message)
	    endif
	    imap(nnodes)=i
	    inverse_map(i)=nnodes
	    do j=1,ndof
	      x0(j,nnodes) = x(j,i)
	    enddo
c z-component (not needed)
	    x0(ndof+1,nnodes) = x(ndof+1,i)

	    do j=1,ndof
	      if((is_relaxed(i).eq.2).or.(id(j,i).eq.1)) then
	        nfixed=nfixed+1
	        if(nfixed.gt.maxfixed) then
	          error_message = 'maxfixed too small (id=1 + 3*n_int)'
	          call error_handler(error_message)
	        endif
	        ifixed(nfixed)=ndof*(nnodes-1)+j
	        ifix_hold(1,nfixed)=j
	        ifix_hold(2,nfixed)=nnodes
	      endif
	    enddo
	  else
	    inverse_map(i)=0
c
	    if(is_relaxed(i).eq.-1) then
	      npad = npad + 1
	      if(npad.gt.maxpad) then
	        error_message = 'maxpad too small'
	        call error_handler(error_message)
	      endif
	      padmap(npad) = i
	    endif
c
	  endif
	enddo
c
	nelm=0
	nsegm = 0
c
	do i=1,numel
	  if(ix(4,i).eq.0) then
	    if((inverse_map(ix(1,i)).eq.0).or.
     &         (inverse_map(ix(2,i)).eq.0).or.
     &         (inverse_map(ix(3,i)).eq.0)) then
	      print *, 'Element number: ', i
	      error_message =
     &          'Inconsistency between ix and id in this element'
	      call error_handler(error_message)
	    endif
	    nelm = nelm+1
	    if(nelm.gt.maxlmn) then
	      error_message =
     &          'maxlmn (els with ix=0) needs to be increased'
	      call error_handler(error_message)
	    endif
	    inverse_el_map(i) = nelm
	    do j=1,knode
	      iconn(j,nelm)=inverse_map(ix(j,i))
	      iadj(j, nelm) = itx(j, i)
	    enddo
c
	    do j=1, knode
	      if((itx(j,i).eq.0).or.(ix(4, itx(j,i)).ne.0)) then
	        nsegm = nsegm+1
	        if(nsegm.gt.maxsegm) then
	          error_message = 'maxsegm needs to be increased'
	          call error_handler(error_message)
	        endif
	        isegm(1,nsegm)=iconn(j,nelm)
	        ij=j+1
	        if(ij.gt.knode) ij=ij-knode
	        isegm(2,nsegm)=iconn(ij,nelm)
	        call linecoef(x0(1,isegm(1,nsegm)),x0(1,isegm(2,nsegm))
     $               ,absegm(1,nsegm))
	      endif
	    enddo
	  else
	    inverse_el_map(i) = 0
	  endif
	enddo

	do i = 1, nelm
	  do j = 1, knode
	    if(iadj(j, i).ne.0) then
	      iadj(j, i) = inverse_el_map(iadj(j, i))
	    endif
	  enddo
	enddo

	if (npad.ne.0) then
	  j = fe_locate(x(1,padmap(1)),0)
	  if(j.eq.0) then
	    print *, 'Atom #', padmap(1)
	    print *, (x(j, padmap(1)), j=1,ndof)
	    error_message = 
     &        'fem_setup: Pad atom outside the continuum region'
	    call error_handler(error_message)
	  endif
	  padelmnt(1) = j
	  call fe_tricoord(x0(1,iconn(1,j)),x0(1,iconn(2,j)),
     &      x0(1,iconn(3,j)),x(1,padmap(1)), padtricoord(1,1))
	  do i = 2, npad
	    j = fe_locate(x(1, padmap(i)),padelmnt(i-1))
	    if(j.eq.0) then
	      print *, 'Atom #', padmap(i)
	      print *, (x(j, padmap(i)), j=1,ndof)
	      error_message = 
     &          'fem_setup: Pad atom outside the continuum region'
	      call error_handler(error_message)
	    endif
	    padelmnt(i) = j
	    call fe_tricoord(x0(1,iconn(1,j)),x0(1,iconn(2,j)),
     &        x0(1,iconn(3,j)),x(1,padmap(i)), padtricoord(1,i))
	  enddo
	endif
c
	call fe_makemat()
	call fe_fixdsp()
	call fe_chol(ier)
	!print *,'ier: ', ier
	!print *
	if(ier.ne.0) then
	  error_message = 'Matrix inversion failed!'
	  call error_handler(error_message)
	endif
	return
	end



	subroutine
     &    fd_solve(b, bc, prop, cz, id, is_relaxed, rhs, forces, e0)
c
c123456789012345678901234567890123456789012345678901234567890123456789012
c         1         2         3         4         5         6         7
c
c if the node is on the interface u=b
c if the node is on the boundary u=prop*bc
c
c  Returns rhs, forces, e0
c  This subroutine should never be called!!!
c
	implicit none
	include 'fem_parameters.par'
	double precision bc(3,*), b(3, *)
	integer id(3,*), is_relaxed(*)
	double precision prop, cz

	double precision rhs(*), forces(*), e0
	double precision u_bc(maxfixed), f_segm(ndof, maxsegm)
	double precision presv(maxfixed)
	integer i, j

	call fd_update_u_tilde_bc(u_bc)
	call fd_update_f_tilde_bc(f_segm)

	do 1 i=1,nequ
 1        rhs(i)=0.0d0

c body forces and boundary tractions

	do 2 i=1,nnodes
	do 2 j=1,ndof

c add body forces to all unconstrained d.o.f. including the interface
c fe_substitute() overwrites interfacial nodes

	  if(id(j,imap(i)).eq.0) then
	    rhs((i-1)*ndof+j) = prop*bc(j, imap(i))/cz
	    !rhs((i-1)*ndof+j) = bc(j, imap(i))/cz
	  endif
 2      continue

	do i=1,nfixed
	  if (is_relaxed(imap(ifix_hold(2,i))).eq.2) then
	    presv(i)=b(ifix_hold(1,i), imap(ifix_hold(2,i)))
	  else
	    presv(i)=prop*bc(ifix_hold(1,i), imap(ifix_hold(2,i)))
	  endif
	enddo

c Set up the b.c. for the ^ field

	call fd_corrective_displ_bc(presv, u_bc)
	call fd_corrective_force_bc(rhs, f_segm)

	call fe_substitute(rhs, presv)

	return
	end

	subroutine fe_force_energy(rhs, forces, e0)
	implicit none
	include 'fem_parameters.par'
	double precision rhs(*), forces(*), e0
	integer i, jmin, jmax, j

	do i=1,nequ
	  forces(i) = a_stiff(1, i)*rhs(i)
	  jmax=mbandw
	  if(jmax+i-1.gt.nequ) jmax=nequ-i+1
	  do j=2, jmax
	    forces(i)=forces(i)+a_stiff(j, i)*rhs(i+j-1)
	  enddo

	  jmin=mbandw-1
	  if((i+1-mbandw).le.0) jmin=i-1
	  do j=1, jmin
	    forces(i)=forces(i)+a_stiff(j+1, i-j)*rhs(i-j)
	  enddo
	enddo

	e0=0.0d0
	do i=1,nequ
	  e0=e0+rhs(i)*forces(i)/2.0d0
	enddo
	return
	end


	subroutine error_handler(message)
	character*80 message
	print '(A,A)', 'Error: ', message
	stop
	end
************************************************************************
      subroutine FindEntryPoint(b,r,xd)
      implicit none
      include 'fem_parameters.par'
      double precision b(3),r(2),xd(3),xint(2),dist2,rdist
      integer i
      logical between
      double precision a1,a2,a3,bot,tol
      parameter(tol=1.e-6)
      integer ism
c
c     eqn of line on which disl moves: a1*x+a2*y+a3=0
c
      if(b(1).eq.0.d0.and.r(1).eq.0.d0) then
	 a1=1.d0
	 a2=0.d0
	 a3=0.d0
      else if(b(2).eq.0.d0.and.r(2).eq.0.d0) then
	 a1=0.d0
	 a2=1.d0
	 a3=0.d0
      else
	 a1=b(2)/(b(1)*r(2)-b(2)*r(1))
	 a2=-b(1)/(b(1)*r(2)-b(2)*r(1))
	 a3=1.d0
      endif
c
c     find all intersections between boundary segments and this line,
c     keep the closest one as the entry point
c     
      rdist=1.e30
      do ism=1,nsegm
	 bot=(a1*absegm(2,ism)-absegm(1,ism)*a2)
	 if(abs(bot).gt.tol) then
	    xint(1)=(absegm(3,ism)*a2-absegm(2,ism)*a3)/bot
	    xint(2)=(-absegm(3,ism)*a1+absegm(1,ism)*a3)/bot
c
c found an intersection:
c
	    if(between(x0(1,isegm(1,ism)),x0(1,isegm(2,ism)),xint)) then
	       dist2=(r(1)-xint(1))**2+(r(2)-xint(2))**2
	       if(dist2.lt.rdist) then
		  rdist=dist2
		  xd(1:2)=xint(1:2)
	       endif
	    endif
	 endif
      enddo
      end      
************************************************************************
      subroutine sliprange(b,r,range,index)
      implicit none
      include 'fem_parameters.par'
      double precision b(3),r(2),range(2),xint(2),RANGETOL
      parameter(RANGETOL=0.5)
      integer index,i
      logical between
      double precision a1,a2,a3,bot,tol
      parameter(tol=1.e-6)
      integer ism
c
c     eqn of line on which disl moves: a1*x+a2*y+a3=0
c
      if(b(1).eq.0.d0.and.r(1).eq.0.d0) then
	 a1=1.d0
	 a2=0.d0
	 a3=0.d0
      else if(b(2).eq.0.d0.and.r(2).eq.0.d0) then
	 a1=0.d0
	 a2=1.d0
	 a3=0.d0
      else
	 a1=b(2)/(b(1)*r(2)-b(2)*r(1))
	 a2=-b(1)/(b(1)*r(2)-b(2)*r(1))
	 a3=1.d0
      endif
c
c     depending on orientation of b, use x or y for comparisons
c      
      if(abs(b(1)).gt.abs(b(2))) then
	 index=1
      else
	 index=2
      endif
c
c     find all intersections between boundary segments and this line,
c     limit dislocation motion accordingly
c     
      range(1)=-1.e30
      range(2)=1.e30
      do ism=1,nsegm
	 bot=(a1*absegm(2,ism)-absegm(1,ism)*a2)
	 if(abs(bot).gt.tol) then
	    xint(1)=(absegm(3,ism)*a2-absegm(2,ism)*a3)/bot
	    xint(2)=(-absegm(3,ism)*a1+absegm(1,ism)*a3)/bot
c
c found an intersection:
c
	    if(between(x0(1,isegm(1,ism)),x0(1,isegm(2,ism)),xint)) then
	       if(r(index).lt.xint(index)) then
	          range(2)=min(xint(index),range(2))
	       else
	          range(1)=max(xint(index),range(1))               
	       endif
	    endif
	 endif
      enddo

c     set the closeness to the interface allowed.

      range(2)=range(2)-RANGETOL
      range(1)=range(1)+RANGETOL
      if(index.eq.1) then
	 write(*,*) 'X-range for the dislocation',range(1),range(2)
      else
	 write(*,*) 'Y-range for the dislocation',range(1),range(2)
      endif
      end

