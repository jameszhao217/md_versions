c
c $Id: fem_alan.f,v 1.2 2004/04/21 14:29:39 shastry Exp $

        subroutine fe_makemat()

c123456789012345678901234567890123456789012345678901234567890123456789012
c         1         2         3         4         5         6         7

        implicit none
        integer maxsts, maxstn
        parameter (maxsts=4)
        parameter (maxstn=3)
        include 'fem_parameters.par'
        double precision d(maxsts, maxstn)
        double precision xl(knode*ndof), s(knode*ndof, knode*ndof)

        integer nodedf, mdif, i1, i2, mtest, i, j
        character*80 error_message

        call fe_dmat(d)

c compute bandwidth

        nequ =ndof*nnodes

        nodedf=0
        do i=1,nelm
          mdif=0
          do 3 i1=1,knode
          do 3 i2=1,knode
            mtest=abs(iconn(i1,i)-iconn(i2,i))
            if (mtest.ge.mdif) mdif=mtest
 3        continue
          if (mdif.ge.nodedf) nodedf=mdif
        enddo
        mbandw=ndof*(nodedf+1)
        !write(*,10) nequ,mbandw
 10     format(' no. equations, bandwidth= ', 2i8)
        if (mbandw.gt.maxbnd) then
          error_message = ' *** not enough storage *** '
          call error_handler(error_message)
        end if
c
        do i=1,nequ
          do j=1,mbandw
          a_stiff(j,i)=0.0d0
          ad_stiff(j,i)=0.0d0
          enddo
        enddo
c
        do i=1,nelm
          call fe_find(i, xl)
          call fe_stiff(xl, s, d)
          call fe_place(i, s)
        enddo
        return
        end


        subroutine fe_find(lmn, xl)
        implicit none
        include 'fem_parameters.par'
        integer lmn
        double precision xl(knode*ndof)
        integer j, k
c
        do 3 j=1,knode
        do 3 k=1,ndof
 3      xl(ndof*(j-1)+k)=x0(k, iconn(j,lmn))
        return
        end



        subroutine fe_dmat(d)
        implicit none
        integer maxsts, maxstn
        parameter(maxsts=4)
        parameter(maxstn=3)
        double precision d(maxsts,maxstn)
        character*80 error_message
c
        double precision cc(6,6)
        double precision xe, xnu, xlambda, xmu
        integer i_elas
c
        common /elastic/ xe, xnu, xlambda, xmu, cc, i_elas
c
c  1<=> eps(1,1)
c  2<=> eps(2,2)
c  3<=> 2*eps(1,2)
c
c  1<=> s(1,1)
c  2<=> s(2,2)
c  3<=> s(1,2)
c  4<=> s(3,3)
c
        if(i_elas.ne.1) then
          error_message = 'fe_dmat: call fe_elastic first!'
          call error_handler(error_message)
        endif
c
        d(1,1)=cc(1,1)
        d(1,2)=cc(1,2)
        d(1,3)=cc(1,6)
        d(2,1)=cc(2,1)
        d(2,2)=cc(2,2)
        d(2,3)=cc(2,6)
        d(3,1)=cc(6,1)
        d(3,2)=cc(6,2)
        d(3,3)=cc(6,6)
        d(4,1)=cc(3,1)
        d(4,2)=cc(3,2)
        d(4,3)=cc(3,6)
        return
        end



        subroutine fe_stiff(xl, s, d)
        implicit none
        integer maxsts, maxstn
        parameter (maxsts=4)
        parameter (maxstn=3)
        include 'fem_parameters.par'
        double precision d(maxsts, maxstn)
        double precision xl(knode*ndof), s(knode*ndof, knode*ndof)
        double precision pn(knode), qn(knode)
        double precision bpq(ndof, ndof, knode*ndof),
     &                   b(maxstn, knode*ndof)
        double precision xjac(ndof, ndof), xji(ndof, ndof)
c
        integer i, j, i1, i2, j1, j2
        double precision det, area
c
c Everything is hard-coded, since I don't know where I may need
c the derivatives of shape functions
c
c  1<=> eps(1,1)
c  2<=> eps(2,2)
c  3<=> 2*eps(1,2)
c
c Shape functions are: N_1 = 1 - p - q; N_2 = p; N_3 = q;
c
        pn(1)=-1.0d0
        pn(2)=1.0d0
        pn(3)=0.0d0
        qn(1)=-1.0d0
        qn(2)=0.0d0
        qn(3)=1.0d0
c
        do 110 i=1,knode*ndof
        do 110 i2=1,ndof
        do 110 i1=1,ndof
 110      bpq(i1,i2,i)=0.0d0
c
        do i=1,knode
          bpq(1,1,ndof*(i-1)+1) = pn(i)
          bpq(1,2,ndof*(i-1)+1) = qn(i)
          bpq(2,1,ndof*(i-1)+ndof) = pn(i)
          bpq(2,2,ndof*(i-1)+ndof) = qn(i)
        enddo
c
        do 2 i=1,ndof
        do 2 j=1,ndof
2         xjac(i,j)=0.0d0
c
        do i=1, knode
          xjac(1,1)=xjac(1,1)+pn(i)*xl((i-1)*ndof+1)
          xjac(1,2)=xjac(1,2)+qn(i)*xl((i-1)*ndof+1)
          xjac(2,1)=xjac(2,1)+pn(i)*xl((i-1)*ndof+2)
          xjac(2,2)=xjac(2,2)+qn(i)*xl((i-1)*ndof+2)
        enddo
c
        det = xjac(1,1)*xjac(2,2) - xjac(1,2)*xjac(2,1)
c
        xji(1,1) = xjac(2,2)/det
        xji(2,2) = xjac(1,1)/det
        xji(1,2) = -xjac(1,2)/det
        xji(2,1) = -xjac(2,1)/det
c
        area = abs(det)/2.0d0
c
        do i1=1,knode*ndof
          b(1,i1) = bpq(1,1,i1)*xji(1,1)+bpq(1,2,i1)*xji(2,1)
          b(2,i1) = bpq(2,1,i1)*xji(1,2)+bpq(2,2,i1)*xji(2,2)
          b(3,i1) = bpq(1,1,i1)*xji(1,2)+bpq(1,2,i1)*xji(2,2)
     &      +bpq(2,1,i1)*xji(1,1)+bpq(2,2,i1)*xji(2,1)
        enddo
c
        do 10 i1=1, ndof*knode
        do 10 i2=1, ndof*knode
          s(i1,i2)=0.0d0
          do 15 j1=1,maxstn
          do 15 j2=1,maxstn
  15        s(i1,i2)=s(i1,i2)+b(j1,i1)*d(j1,j2)*b(j2,i2)
  10      s(i1,i2) = s(i1,i2)*area
c
        return
        end



        subroutine fe_place(lmn,s)
        implicit none
        include 'fem_parameters.par'
        integer lmn
        double precision  s(knode*ndof, knode*ndof)
        integer i1, i2, i, k, j1, j2, j, l
c
      do 1 i1=1,ndof
      do 1 i2=1,knode
      i=ndof*(i2-1)+i1
      k=ndof*(iconn(i2,lmn)-1)+i1
      do 1 j1=1,ndof
      do 1 j2=1,knode
      j=ndof*(j2-1)+j1
      l=ndof*(iconn(j2,lmn)-1)+j1
      if (k.ge.l) then
      a_stiff(k-l+1, l)=a_stiff(k-l+1, l)+s(i,j)
      ad_stiff(k-l+1, l)=ad_stiff(k-l+1, l)+s(i,j)
      end if
 1    continue
      return
      end


      subroutine fe_fixdsp()
        implicit none
        include 'fem_parameters.par'
        integer i, ic, j
        double precision hold_stiff
c
      hold_stiff=0.0d0
      do 3 i=1,nequ
      if (abs(ad_stiff(1, i)).gt.hold_stiff) then
        hold_stiff=abs(ad_stiff(1, i))
      endif
 3    continue
      do 1 i=1,nfixed
      ic=ifixed(i)
      do 2 j=1,nequ
      if ((j.lt.ic).and.(ic-j+1.le.mbandw)) then
      ad_stiff(ic-j+1, j)=0.0d0
      else if ((j.gt.ic).and.(j-ic+1.le.mbandw)) then
      ad_stiff(j-ic+1, ic)=0.0d0
      end if
 2    continue
      ad_stiff(1, ic)=hold_stiff
 1    continue
      return
      end


      subroutine fe_substitute(rhs, presv)
        implicit none
        include 'fem_parameters.par'
        double precision rhs(*), presv(*)
        integer i, ic, j, k, kmin, kmax
c
      do 1 i=1,nfixed
      ic=ifixed(i)
      do 2 j=1,nequ
      if ((j.lt.ic).and.(ic-j+1.le.mbandw)) then
      rhs(j)=rhs(j)-presv(i)*a_stiff(ic-j+1, j)
      else if ((j.gt.ic).and.(j-ic+1.le.mbandw)) then
      rhs(j)=rhs(j)-presv(i)*a_stiff(j-ic+1, ic)
      end if
 2    continue
 1    continue
      do 3 i=1,nfixed
      ic=ifixed(i)
 3    rhs(ic)=presv(i)*ad_stiff(1,ic)
c
c substitution
c
      do 6 i=2,nequ
      kmin=i+1-mbandw
      if (kmin.lt.1) kmin=1
      do 6 k=kmin,(i-1)
 6    rhs(i)=rhs(i)-ad_stiff(i-k+1, k)*rhs(k)
c
c divide by diagonals
c
      do 7 i=1,nequ
 7    rhs(i)=rhs(i)/ad_stiff(1, i)
c
c back substitution
c
      do 8 i=(nequ-1),1,-1
      kmax=i-1+mbandw
      if (kmax.gt.nequ) kmax=nequ
      do 8 k=(i+1),kmax
 8    rhs(i)=rhs(i)-ad_stiff(k-i+1, i)*rhs(k)
      return
      end



      subroutine fe_chol(ier)
        implicit none
        include 'fem_parameters.par'
        integer ier, i, k, j, kmin, jmax
        double precision piv, test
c
      ier=0
      test=0.0d0
      do 10 i=1,nequ
      if (ad_stiff(1, i).gt.test) test=ad_stiff(1, i)
 10   continue
      test=1.d-8*test
c
c factorization
c
      if (ad_stiff(1,1).lt.test) then
      ier=1
      return
      end if
      piv=1.0d0/ad_stiff(1,1)
      do 1 k=2,mbandw
      ad_stiff(k, 1)=piv*ad_stiff(k, 1)
 1    continue
      do 2 i=2,nequ
      kmin=i+1-mbandw
      if (kmin.lt.1) kmin=1
      do 3 k=kmin,(i-1)
 3    ad_stiff(1, i)=ad_stiff(1, i)
     &            -ad_stiff(i-k+1, k)*ad_stiff(1, k)*ad_stiff(i-k+1, k)
      if (i.lt.nequ) then
      if (ad_stiff(1, i).lt.test) then
      ier=i
      return
      end if
      piv=1.0d0/ad_stiff(1, i)
      jmax=i-1+mbandw
      if (jmax.gt.nequ) jmax=nequ
      do 4 j=(i+1),jmax
      kmin=j+1-mbandw
      if (kmin.lt.1) kmin=1
      do 5 k=kmin,(i-1)
 5    ad_stiff(j-i+1, i)=ad_stiff(j-i+1, i)
     &  -ad_stiff(i-k+1, k)*ad_stiff(1, k)*ad_stiff(j-k+1, k)
 4    ad_stiff(j-i+1, i)=piv*ad_stiff(j-i+1, i)
      end if
 2    continue
      return
      end



        subroutine fe_strain(lmn, rhs, ee_out)
        implicit none
        integer maxsts, maxstn
        parameter (maxsts=4)
        parameter (maxstn=3)
        include 'fem_parameters.par'
c
        integer lmn
        double precision rhs(*), ee_out(maxstn)
c
        double precision xl(knode*ndof), s(knode*ndof, knode*ndof)
        double precision pn(knode), qn(knode)
        double precision bpq(ndof, ndof, knode*ndof),
     &                   b(maxstn, knode*ndof)
        double precision xjac(ndof, ndof), xji(ndof, ndof)
c
        integer i, j, i1, i2, j1, j2
        double precision det

        character*80 error_message

        if((lmn.lt.1).or.(lmn.gt.nelm)) then
          print *, lmn
          error_message = 'fe_strain: element is out of range!'
          call error_handler(error_message)
        endif

        call fe_find(lmn, xl)
c
c Everything is hard-coded, since I don't know where I may need
c the derivatives of shape functions
c
c  1<=> eps(1,1)
c  2<=> eps(2,2)
c  3<=> 2*eps(1,2)
c
c Shape functions are: N_1 = 1 - p - q; N_2 = p; N_3 = q;
c
        pn(1)=-1.0d0
        pn(2)=1.0d0
        pn(3)=0.0d0
        qn(1)=-1.0d0
        qn(2)=0.0d0
        qn(3)=1.0d0
c
        do 110 i=1,knode*ndof
        do 110 i2=1,ndof
        do 110 i1=1,ndof
 110      bpq(i1,i2,i)=0.0d0
c
        do i=1,knode
          bpq(1,1,ndof*(i-1)+1) = pn(i)
          bpq(1,2,ndof*(i-1)+1) = qn(i)
          bpq(2,1,ndof*(i-1)+ndof) = pn(i)
          bpq(2,2,ndof*(i-1)+ndof) = qn(i)
        enddo
c
        do 2 i=1,ndof
        do 2 j=1,ndof
2         xjac(i,j)=0.0d0
c
        do i=1, knode
          xjac(1,1)=xjac(1,1)+pn(i)*xl((i-1)*ndof+1)
          xjac(1,2)=xjac(1,2)+qn(i)*xl((i-1)*ndof+1)
          xjac(2,1)=xjac(2,1)+pn(i)*xl((i-1)*ndof+2)
          xjac(2,2)=xjac(2,2)+qn(i)*xl((i-1)*ndof+2)
        enddo
c
        det = xjac(1,1)*xjac(2,2) - xjac(1,2)*xjac(2,1)
c
        xji(1,1) = xjac(2,2)/det
        xji(2,2) = xjac(1,1)/det
        xji(1,2) = -xjac(1,2)/det
        xji(2,1) = -xjac(2,1)/det
c
        do i=1,knode*ndof
          b(1,i) = bpq(1,1,i)*xji(1,1)+bpq(1,2,i)*xji(2,1)
          b(2,i) = bpq(2,1,i)*xji(1,2)+bpq(2,2,i)*xji(2,2)
          b(3,i) = bpq(1,1,i)*xji(1,2)+bpq(1,2,i)*xji(2,2)
     &      +bpq(2,1,i)*xji(1,1)+bpq(2,2,i)*xji(2,1)
        enddo
c
        call fe_find_rhs(lmn, xl, rhs)
        do i = 1, maxstn
          ee_out(i) = 0.0d0
          do j=1, knode*ndof
            ee_out(i) = ee_out(i) + b(i,j)*xl(j)
          enddo
        enddo
        return
        end



        subroutine fe_find_rhs(lmn, xl, rhs)
        implicit none
        include 'fem_parameters.par'
        integer lmn
        double precision xl(knode*ndof), rhs(*)
        integer j, k
c
        do 3 j=1,knode
        do 3 k=1,ndof
 3      xl(ndof*(j-1)+k)=rhs((iconn(j,lmn)-1)*ndof+k)
        return
        end



        subroutine fe_stress(lmn, rhs, s_out)
        implicit none
        integer maxsts, maxstn
        parameter (maxsts=4)
        parameter (maxstn=3)
        include 'fem_parameters.par'
c
        integer lmn
        double precision rhs(*), s_out(3)
c
        double precision d(maxsts, maxstn), ee_out(maxstn)
        integer i, j
        character*80 error_message
c
c  1<=> s(1,1)
c  2<=> s(2,2)
c  3<=> s(1,2)
c
        if((lmn.lt.1).or.(lmn.gt.nelm)) then
          print *, lmn
          error_message = 'fe_strain: element is out of range!'
          call error_handler(error_message)
        endif
c
        call fe_dmat(d)
        call fe_strain(lmn, rhs, ee_out)
        do i = 1, 3
          s_out(i) = 0
          do j=1,maxstn
            s_out(i) = s_out(i) + d(i,j)*ee_out(j)
          enddo
        enddo
c
        return
        end



        function fe_locate(r,iel0)
c
c  Returns the element number corresponding to position r or zero.
c  Breadth first search in FORTRAN
c
        implicit none
        include 'fem_parameters.par'
        integer fe_locate, iel0
        double precision r(ndof+1)
c
        integer Q_MAX
        parameter(Q_MAX=512)
        integer queue(Q_MAX)
        integer idummy, iel, ptr_put, ptr_get, i, lmn
        logical todo(maxlmn), fe_in_tri
c
        character*80 error_message
c
        if(i_flag.ne.1) then
          error_message = 'fe_locate: call fem_setup first!'
          call error_handler(error_message)
        endif
c
        iel = iel0
        if(abs(iel).gt.nelm) then
          print *, 'iel: ', iel, ' nelm: ', nelm
          error_message = 'fe_locate: initial element out of boundary'
          call error_handler(error_message)
        endif
        if( iel.eq.0 ) iel = 1
        if( iel.lt.0 ) iel = -iel
c
c  Checking iel in order not to waste time on array initialization
c
        if(fe_in_tri(x0(1,iconn(1,iel)),x0(1,iconn(2,iel)),
     &    x0(1,iconn(3,iel)),r)) then
          fe_locate=iel
          return
        endif
c
c initialize the fringe
c
        do i=1,nelm
          todo(i)=.true.
        enddo
c
c Breadth-first search checking first iel again to simplify the code
c Circular queue of size Q_MAX
c ptr_put pointer to where to put an element
c ptr_get where the last element was taken (i.e., needs to be incremented
c before taking something out of the queue).
c
        ptr_get = Q_MAX
        ptr_put = 2
        queue(1) = iel
        todo(iel) = .false.
c
        do idummy = 1, nelm
          ptr_get = ptr_get + 1
          if(ptr_get.eq.Q_MAX+1) ptr_get = 1
          if(ptr_get.eq.ptr_put) then
c  Empty queue. The search over the connected component is exhasted
c  Probably we are done
            goto 100
          endif
c
          iel = queue(ptr_get)
          if(fe_in_tri(x0(1,iconn(1,iel)),x0(1,iconn(2,iel)),
     &        x0(1,iconn(3,iel)),r)) then
            fe_locate=iel
            return
          endif
          do i = 1, 3
            lmn = iadj(i,iel)
            if((lmn.ne.0).and.todo(lmn)) then
              if(ptr_put.eq.ptr_get) then
                print *, 'Q_MAX: ', Q_MAX
                error_message = 'fe_locate: queue overflow'
                call error_handler(error_message)
              endif
              queue(ptr_put) = lmn
              todo(lmn) = .false.
              ptr_put = ptr_put + 1
              if(ptr_put.eq.Q_MAX+1) ptr_put = 1
            endif
          enddo
        enddo
c
 100    fe_locate=0
        return
        end



        function fe_in_tri(r1, r2, r3, p)
c
c  Wrapper for intri from feaplib
c
        implicit none
        double precision r1(3), r2(3), r3(3), p(3)
        logical fe_in_tri
c
        double precision x1(2), x2(2), x3(2), pt(2), s(3)
        logical intri, ontri
        integer i
c
        do i = 1, 2
          x1(i) = r1(i)
          x2(i) = r2(i)
          x3(i) = r3(i)
          pt(i) = p(i)
        enddo
        fe_in_tri = intri(x1,x2,x3,pt,s,ontri)
        return
        end


c
c $Log: fem_alan.f,v $
c Revision 1.2  2004/04/21 14:29:39  shastry
c vijay-    fem_parameters.par fem_alan.f: increased storage.
c
c Revision 1.1.1.1  2003/03/12 20:09:00  shastry
c vijay-   Initial import.
c
c Revision 1.7  2001/12/13 07:31:24  shilkrot
c Implemented breadth first search to find the element number for
c a dislocation. Changed the interface of fe_locate to include the starting
c element for the search. Old fe_locate is in fem_services.
c Changed the interface of fem_setup. Now two arrays used as temp space are
c passed from outside as the last two parameters.
c
c Revision 1.6  2001/11/13 03:32:44  shilkrot
c Added functions computing strain and stress for a given element.
c Wrote a function fe_locate which must be replaced in the future.
c
c Revision 1.5  2001/08/22 03:18:35  shilkrot
c Fixed the expression for the energy and polished fem_alan a little bit.
c This wersion works with dislocation passing.
c
c Revision 1.4  2001/07/12 06:36:49  shilkrot
c Elastic constants are now passed to fem_setup and further to fe_elastic
c
c Revision 1.3  2001/06/25 18:55:53  shilkrot
c Changed stiffness matrix for the transposed one
c Otherwise, no major changes
c
c Revision 1.2  2001/06/18 01:19:17  shilkrot
c a_stiff becomes ad_stiff and vice versa
c
c Revision 1.1  2001/06/18 00:17:25  shilkrot
c The routines initially written by A. Needleman
c
c
