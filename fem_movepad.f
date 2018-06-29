c
c $Id: fem_movepad.f,v 1.1.1.1 2003/03/12 20:09:00 shastry Exp $

	subroutine
     &    fem_move_pad(x,b,ix,bc,prop,cz,id,is_relaxed,e_out,f_out,
     &    FullField, MoveDisl,strainE0,numel,avevirst,rank)

c123456789012345678901234567890123456789012345678901234567890123456789012
c         1         2         3         4         5         6         7
c
c if the node is on the interface u=b
c if the node is on the boundary u=prop*bc
c changes b 
c e_out for informational purposes only

	implicit none
	include 'fem_parameters.par'
	include 'disl_parameters.par'
	double precision x(3,*), b(3,*), bc(3,*), f_out(3, *)
	double precision avevirst(3,3,*)
	integer numel,rank
	integer id(3,*), ix(4,*),is_relaxed(*)
	double precision prop, cz
	double precision e_out,strainE0
	logical FullField, MoveDisl
c
	double precision rhs(maxeqs), forces(maxeqs), e0, eplus, eminus
	integer i
	integer fe_locate, elem_old, elem_plus, elem_minus
	character*80 error_message
	logical CentralDiff,ForwardDiff,BackDiff
!
	if(i_flag.ne.1) then
	  error_message = 'fem_solve: call fem_setup first!'
	  call error_handler(error_message)
	endif
	if(i_disl.ne.1) then
	  error_message = 'fem_solve: call disl_setup first!'
	  call error_handler(error_message)
	endif

!       solve FEM and DD problem
	call
     &    fd_solve(b, bc, prop, cz, id, is_relaxed, rhs, forces, e0)
!
!  compute the P.-K. force on dislocations
	if(MoveDisl) then
	  call fd_peach_koeller(rhs)
	endif
!
	do i=1,nfixed
	  f_out(ifix_hold(1,i), imap(ifix_hold(2,i))) = 
     &      f_out(ifix_hold(1,i), imap(ifix_hold(2,i)))
     &      -forces(ifixed(i))*cz
	enddo
	e_out = e0*cz
!
!       move dislocations based on P.-K. force
        if(MoveDisl) call move_dis(10.0)
!
!       move pad atoms according to tilda and hat fields
	call fd_movepad(x, rhs, b)
!
!       update b so b=u_hat+u_tilda
	if(FullField) then
	  call fd_full_field(rhs, b)
	endif
!
!       compute total strain energy of fem region
        call fe_stress_check(rhs,strainE0)

	return
	end



	subroutine fd_movepad(x, rhs, b)
	implicit none
	include 'fem_parameters.par'
	double precision x(3, *), rhs(*), b(3,*)
c
	double precision u_out(3)
	integer i, j, k, lmn, n
c
	do i = 1, npad
	  n = padmap(i)
	  lmn = padelmnt(i)
	  do j = 1, ndof
	    b(j,n) = 0.d0
	    do k = 1, knode
	      b(j,n) = b(j,n) + 
     &          rhs( (iconn(k,lmn)-1)*ndof+j ) * padtricoord(k,i)
	    enddo
	  enddo
	  call disl_displ(x(1, n), u_out)
	  do j = 1, ndof
	    b(j, n) = b(j, n) + u_out(j)
	  enddo
	enddo
	return
	end
	


	subroutine fe_tricoord(r1, r2, r3, p, coord)
c
c  Yet another wrapper for intri from feaplib
c
	implicit none
	double precision r1(3), r2(3), r3(3), p(3), coord(3)
c
	double precision x1(2), x2(2), x3(2), pt(2)
	logical intri, ontri
	integer i
	character*80 error_message
c
	do i = 1, 2
	  x1(i) = r1(i)
	  x2(i) = r2(i)
	  x3(i) = r3(i)
	  pt(i) = p(i)
	enddo
	if(.not.intri(x1,x2,x3,pt,coord,ontri)) then
	  error_message = 'fe_tricoord: the atom is not in the element'
	  call error_handler(error_message)
	endif
	return
	end

c
c $Log: fem_movepad.f,v $
c Revision 1.1.1.1  2003/03/12 20:09:00  shastry
c vijay-   Initial import.
c
c


!!!!!!!!dw added subroutine!!!!!!!!!!!!!!!!!!
        subroutine move_dis(alpha)
        use MPI
	! use mod_global


        implicit none
        include 'disl_parameters.par'

        !integer nprocs
        integer rank,ierr
        double precision alpha, mobility, max_vel
        double precision sf_f(10000), aa, bb, friction_f
        double precision min_pos, temperature, time_step_con
        integer fe_locate, i, elem_old
        double precision b
        double precision ddis
        character*80 error_message
        integer ndisl_old, kris_counter, ndisl_kris, kris_time
        character *13 plot_name
        character *5 plot_num

!!!! MZ friction force parameters

	double precision JAMES_friction


!!!!    hacked parameters
        data ndisl_old /0/
        data kris_time /0/


        !call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

        min_pos=-20.0
        temperature=100.0
        time_step_con=5.0e-4
        mobility=time_step_con/(6.242e-2*5.0e-8*temperature)
        ! 3rd parameter is damping coef from Olmstead paper
        max_vel=time_step_con*2000.0
    !    sf_f=0.0*.077*6.242e-2 ! sf energy in J/m2
    !    friction_f=110*6.241506e-6*burg_length(i) 
         ! (flow stress/3) = 283MPa/3 = 94.33MPa = 5.887821e-4eV/ang^3  
         !or 75.76 for GPzones ! or 110 for peak
        if(ndisl.ne.ndisl_old)then
          ndisl_old=ndisl
          !open(88,FILE='sf_data',STATUS='old')
          !read(88,*) ndisl_kris
	  ndisl_kris=1000
          do kris_counter=1,ndisl_kris
            !read(88,*) sf_f(kris_counter)
	    sf_f(kris_counter)=0.d0
          enddo
          !close(88)
        endif
        
!!!!    end of hacked parameters

        ndisl_kris=0
        do i = 1, ndisl
          if(elem_disl(i).ne.0) then
            ndisl_kris=ndisl_kris+1



!	if(r_disl(1, i).gt.70)then
!            if(r_disl(1,i).ge.0.and.
!     & r_disl(2,i).gt.1.63299316185*r_disl(1,i))then
!               friction_f=500*6.241506e-6*burg_length(i) ! (flow stress/3) = 283MPa/3 = 94.33MPa = 5.887821e-4eV/ang^3  !or 75.76 for GPzones ! or 110 for peak
!            else
!               friction_f=75.76*6.241506e-6*burg_length(i)
!            endif
!-            print*,'friction_f',friction_f



            !! SET DISLOCATION FRICTION FORCE !!


!!MZ: added route to set different glide resistances for each simulation

		if(rank.eq.0) then		
	 		open(unit=200, file='md.inp', status='old')
	 		read(200,*) ! mdTemp
         		read(200,*) ! timestep
	 		read(200,*) ! FEMStepMin
	 		read(200,*) ! Nsteps 
	 		read(200,*) ! o2_conc
	 		read(200,*) ! h2o_conc
	 		read(200,*) ! stepsPerFE
	 		read(200,*) ! dumpcg 
	 		read(200,*) ! EAMdebug, Monolayer, monoThick, Adhes
	 		read(200,*) ! Ftol ! Should be in [eV/Ang] units
	 		read(200,*) ! Sbarrier
	 		read(200,*) ! scaleChg
	 		read(200,*) ! K_max, fstart
	 		read(200,*) JAMES_friction  ! in MPa
	 		close(200)
		endif

      		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      		CALL MPI_BCAST(JAMES_friction,1,MPI_DOUBLE,
     &  		0,MPI_COMM_WORLD,ierr)

		! Print James_friction to file 'output'	
		if(rank.eq.0) then
	    		print *,'JAMES - friction',JAMES_friction
		endif	

            friction_f = 1.0 * JAMES_friction * 
     & 			6.241506e-6 * burg_length(i) 

!!! MZ: ended

            !friction_f = 6.241506e-6 * burg_length(i) * JAMES_friction
            !friction_f = 2.0 * 110.0 * 6.241506e-6 * burg_length(i)
            !friction_f = 5.0 * 110.0 * 6.241506e-6 * burg_length(i)

            pk_f(i) = (pk_force(1,i)*burgers(1,i)
     &                + pk_force(2,i)*burgers(2,i))/burg_length(i)
!-            write(*,*) i,' non sf disl force = ',pk_f(i)
            ddis=sqrt(r_disl(1,i)**2+r_disl(2,i)**2)
            pk_f(i) = pk_f(i) + 0.0 !sf_f(i)          !WARNING kris removed all sf_f!
!-            write(*,*) i,'with sf disl force = ',pk_f(i)
            
!!!!        add friction force to resist movment
            if(abs(pk_f(i))-abs(friction_f).le.0.0)then
              pk_f(i)=0.0
              !print*,'disl',i,'is stuck'
            elseif(pk_f(i).ge.0.0)then
              pk_f(i)=pk_f(i)-friction_f
            elseif(pk_f(i).lt.0.0)then
              pk_f(i)=pk_f(i)+friction_f
            else
              print*,'ERROR with friction force!'
              print*,'pk_force',pk_force(1,i),pk_force(2,i)
              print*,'pk_f',pk_f(i)
              !print*,'friction_f',friction_f
              !print*,'burg',burgers(1,i),burgers(2,i)
              print*,'burg_length',burg_length(i)
              pk_f(i)=0.0
              !stop
            endif
             
!           if(pk_f(i).ge.0.0)then
!              if(pk_f(i).gt.friction_f)then
!                pk_f(i)=pk_f(i)-friction_f
!              else
!                pk_f(i)=0.0
!              endif
!           else !(pk_f(i).lt.0.0)then
!              if(pk_f(i).lt.-1*friction_f)then
!                 pk_f(i)=pk_f(i)+friction_f
!              else
!                 pk_f(i)=0.0
!              endif
!!           else
!!              print*,'ERROR with friction force!'
!!              print*,'pk_force',pk_force(1,i),pk_force(2,i)
!!              print*,'pk_f',pk_f(i)
!!              !print*,'friction_f',friction_f
!!              !print*,'burg',burgers(1,i),burgers(2,i)
!!              print*,'burg_length',burg_length(i)
!!              stop
!           endif
!-            print*,i,'disl force with fric',pk_f(i) 

!!!!        upper limit on velocities
            if (abs(pk_f(i)*mobility).gt.max_vel) then
               pk_f(i)=max_vel/mobility*pk_f(i)/abs(pk_f(i))
!-               print*,'vel > max vel'
            endif
!-            write(*,*) 'old disl pos = ',r_disl(1, i), r_disl(2, i),ddis
!-            write(*,*) i,' total disl force = ',pk_f(i)
!            if(r_disl(2,i).gt.-100.0.or.pk_f(i).gt.0.0) then
!            if(r_disl(2,i).le.min_pos.or.pk_f(i).lt.0.0) then
            if(1.eq.1)then
!               if(i.eq.1.and.pk_f(i).le.0.0)then
!               else
             r_disl(1, i) = r_disl(1,i) + 
     &         mobility*pk_f(i)*burgers(1,i)/burg_length(i)
             r_disl(2, i) = r_disl(2,i) + 
     &         mobility*pk_f(i)*burgers(2,i)/burg_length(i)
!               endif
            endif
!            endif
!-            write(*,*) 'new disl pos = ',r_disl(1, i), r_disl(2, i)
c            write(*,*) ' '
!            print*,'DISL',i,r_disl(1,i),r_disl(2,i)
            elem_old = elem_disl(i)
            elem_disl(i)=fe_locate(r_disl(1,i), elem_disl(i))
            if(elem_disl(i).eq.0) then
              if(elem_old.gt.0) then
                elem_disl(i) = -elem_old
              else
                elem_disl(i) = elem_old
              endif
            endif
          endif
        enddo

!        if(ndisl.ne.ndisl_old)then
!          ndisl_old=ndisl
!          open(92,FILE='restart_disl',STATUS='replace')
!          write(92,*),ndisl,' Number of dislocations'
!          do kris_counter=1,ndisl
!            write(92,*) r_disl(1,kris_counter),r_disl(2,kris_counter)
!     &,burgers(1,kris_counter),burgers(2, kris_counter),
!     &burgers(3,kris_counter)
!          enddo
!          close(92) 
!        endif
        kris_time=kris_time+1
        if(mod(kris_time,2000).eq.0)then
        plot_name='restart_dislo'   !print dislocation positions at time of restart
         if(rank.eq.0)then
	  write(plot_num,'(i4)') kris_time/2000+1000
          open(93,FILE=plot_name//plot_num,STATUS='new')
          write(93,*) ndisl_kris,' number of disl '
	 endif
         do kris_counter=1,ndisl
           if(elem_disl(kris_counter).ne.0)then
             if(rank.eq.0) write(93,'(7E27.17)') 
     & r_disl(1,kris_counter),r_disl(2,kris_counter),
     & burgers(1,kris_counter),burgers(2,kris_counter),
     & burgers(3,kris_counter)
     & ,theta_e(kris_counter),theta_s(kris_counter)
           endif
         enddo
         if(rank.eq.0) close(93)
	endif

        return
        end
