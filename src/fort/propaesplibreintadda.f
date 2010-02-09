c     Authors: P. C. Chaumet and A. Rahmani
c     Date: 03/10/2009

c     Purpose: This routine compute the integration of the Green Tensor
c     (field susceptibility tensor) of free space to improve the
c     convergence rate of the discrete dipole approximation method.

c     Reference: if you use this routine in your research, please
c     reference, as appropriate: P. C. Chaumet, A. Sentenac and
c     A. Rahmani "Coupled dipole method for scatterers with large
c     permittivity", Phys. Rev. E 70(3), 036606-6 (2004).

c     license: GNU GPL

      subroutine propaespacelibreintadda(Rij,k0a,arretecube,relreq,
     $                                   result)
      implicit none
      integer i,j
c     definition of the position of the dipole, observation, wavenumber
c     ,wavelength, spacing lattice
      double precision k0a,arretecubem
      double precision x,y,z,arretecube,k0,xx0,yy0,zz0
      double precision Rij(3),result(12)

c     Variables needs for the integration
      integer  KEY, N, NF, NDIM, MINCLS, MAXCLS, IFAIL, NEVAL, NW
      parameter (nw=4000000,ndim=3,nf=12)
      double precision A(NDIM), B(NDIM), WRKSTR(NW)
      double precision  ABSEST(NF), FINEST(NF), ABSREQ, RELREQ,err
      
      double precision Id(3,3),Rab,Rvect(3)

      external fonctionigtadda

      common/k0xyz/k0,x,y,z,xx0,yy0,zz0

      x=Rij(1)
      y=Rij(2)
      z=Rij(3)
      k0=k0a
      arretecubem=arretecube*0.1d0

c     We perform the integration of the tensor
c     definition for the integration
      MINCLS = 1000
      MAXCLS = 1000000
      KEY = 0
      ABSREQ = 0.d0
      
      A(1)=-arretecube/2.d0
      A(2)=-arretecube/2.d0
      A(3)=-arretecube/2.d0
      B(1)=+arretecube/2.d0
      B(2)=+arretecube/2.d0
      B(3)=+arretecube/2.d0
      
      xx0=1.d0
      yy0=1.d0
      zz0=1.d0
      if (dabs(z).le.arretecubem) then
         zz0=0.d0
      endif
      if (dabs(x).le.arretecubem) then
         xx0=0.d0
      endif
      if (dabs(y).le.arretecubem) then
         yy0=0.d0
      endif

      call  DCUHRE(NDIM,NF,A,B, MINCLS, MAXCLS, fonctionigtadda,
     $      ABSREQ,RELREQ,KEY,NW,0,finest,ABSEST,NEVAL,IFAIL, WRKSTR) 
      
      do N = 1,NF
         FINEST(N)=FINEST(N)/arretecube/arretecube/arretecube
      enddo

      if (ifail.ne.0) then
         write(*,*) 'IFAIL in IGT routine',IFAIL
      endif

      do i = 1,6
        result(2*i-1)=finest(i)
        result(2*i)=finest(i+6)
      enddo
    
      end
c*************************************************************
      subroutine fonctionigtadda(ndim,zz,nfun,f)
      implicit none
      integer n,ndim,nfun
      double precision zz(ndim),f(nfun)
      
      integer i,j
      double precision x,y,z,x0,y0,z0,k0,Id(3,3),Rab,Rtenseur(3,3)
     $     ,Rvect(3),xx0,yy0,zz0
      double complex propaesplibre(3,3),const1,const2
      common/k0xyz/k0,x,y,z,xx0,yy0,zz0

      x0=zz(1)
      y0=zz(2)
      z0=zz(3)

      Rab=0.d0
      Rvect(1)=(x-x0)
      Rvect(2)=(y-y0)
      Rvect(3)=(z-z0)

      do i=1,3
         do j=1,3
            Id(i,j)=0.d0
            if (i.eq.j) Id(i,i)=1.d0
            Rtenseur(i,j)=Rvect(i)*Rvect(j)
         enddo
         Rab=Rab+Rvect(i)*Rvect(i)
      enddo
      Rab=dsqrt(Rab)

c     normalise pour avoir le vecteur unitaire
      do i=1,3
         do j=1,3
            Rtenseur(i,j)=Rtenseur(i,j)/(Rab*Rab)
         enddo
      enddo
    
      const1=(Rab*(1.d0,0.d0))**(-3.d0)-(0.d0,1.d0)*k0*(Rab**(-2.d0))
      const2=k0*k0/Rab*(1.d0,0.d0)
      do i=1,3
         do j=1,3
            propaesplibre(i,j)=((3.d0*Rtenseur(i,j)-Id(i,j))*const1+
     *           (Id(i,j)-Rtenseur(i,j))*const2)*
     *           cdexp((0.d0,1.d0)*k0*Rab)
         enddo
      enddo

      f(1)=dreal(propaesplibre(1,1))
      f(2)=dreal(propaesplibre(1,2))*xx0*yy0
      f(3)=dreal(propaesplibre(1,3))*xx0*zz0
      f(4)=dreal(propaesplibre(2,2))
      f(5)=dreal(propaesplibre(2,3))*yy0*zz0
      f(6)=dreal(propaesplibre(3,3))

      f(7)=dimag(propaesplibre(1,1))
      f(8)=dimag(propaesplibre(1,2))*xx0*yy0
      f(9)=dimag(propaesplibre(1,3))*xx0*zz0
      f(10)=dimag(propaesplibre(2,2))
      f(11)=dimag(propaesplibre(2,3))*yy0*zz0
      f(12)=dimag(propaesplibre(3,3))

      end
