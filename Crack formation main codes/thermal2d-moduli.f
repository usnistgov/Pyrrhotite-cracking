c                    Copyright Notice
c
c This software was developed at the National Institute of Standards
c and Technology (NIST) by employees of the Federal Government in the
c course of their official duties. Pursuant to title 17 Section 105
c of the United States Code this software is not subject to copyright
c protection and is in the public domain.
c
c This is an experimental system.  NIST assumes no responsibility
c whatsoever for its use by other parties, and makes no guarantees,
c expressed or implied, about its quality, reliability, or any other
c characteristic.
c
c We would appreciate acknowledgment if the software is used.
c
c
c
c Contact: Osamah H. Dehwah  osamah.dehwah@nist.gov
c          NIST
c          Engineering Laboratory
c          Materials and Structural Systems Division
c          Infrastructure Materials Group

c This version does not allow macrostrains to change
c macrostrains are set at beginning and used to compute
c the elastic moduli, bulk and shear at the same time
c all eigenstrains = 0

c This version also puts some randomness into matrix, to see if that
c makes the cracks have more of a random direction rather than following
c the i-j directions.

c PERIODIC BOUNDARY CONDITIONS - no displacements held at constant
c  ************************  thermal2d.f  ************************
c  BACKGROUND

c  Program adjusts dimensions of unit cell,
c  [(1 + macrostrain) times dimension],
c  in response to phases that have a non-zero eigenstrain and
c  arbitrary elastic moduli tensors.
c  All three macrostrains can adjust their values (2-d program), and are
c  stored in the last two positions in the displacement vector u,
c  as listed below. Periodic boundaries are maintained.
c  In the comments below, (USER) means that this is a section of code
c  that the user might have to change for  his particular problem.
c  Therefore the user is encouraged to search for this string.

c  PROBLEM AND VARIABLE DEFINITION

c  The problem being solved is the minimization of the elastic energy
c  1/2 uAu + bu + C + Tu + Y, where b and C are also functions of the
c  macrostrains.
c  The small array zcon computes the thermal strain energy associated
c  with macrostrains (C term), T is the thermal energy term linear in the
c  displacements (built from ss), b is the regular energy term linear in the
c  displacements, u is the displacements including the macrostrains, gb
c  is the energy gradient vector, h,Ah are auxiliary vectors,
c  dk is the single pixel stiffness matrix, pix is the phase
c  identification vector (pix=1 for phase 1, etc.), and ib is the
c  integer matrix for mapping labels from the 1-9 nearest neighbor
c  labelling to the 1-d system labelling.
c  The array prob(i) contains the volume fractions of the i'th phase,
c  strxx, etc. are the three independent area averaged stresses
c  (Voigt notation), sxx, etc. are the three independent area
c  averaged strains (Voigt notation, not counting the thermal strains),
c  cmod(i,3,3) gives the elastic modulus tensor
c  of the i'th phase, and dk(i,4,2,4,2) is the stiffness matrix of the i'th
c  phase.  The parameter nphase gives the number of phases being considered
c  in the problem, and is set by the user.

c  DIMENSIONS

c  The main arrays of the problem, u, gb, h, Ah, b, and T, are dimensioned
c  as (nx*ny)+2, which is the number of nodal displacements plus two for
c  the macrostrains.
c  The program currently assumes that the number of different phases is
c  100, since phasemod, eigen, and ss (the moduli, eigenstrains, and
c  auxiliary eigenstrain variable for each phase)
c  and dk are dimensioned to have at most 100 different kinds.  This
c  is easily changed.  The parameter nphase gives the number of phases
c  expected in the problem, and is user-specified.  All major arrays
c  are passed to subroutines via simple common statements.

c  NOTE ON USE OF PROGRAM:  Program is set up to allow the macrostrains,
c  which control the overall size of the system, to be dynamic
c  variables, which are adjusted in order to minimize the overall
c  energy.  That means that if there are no eigenstrains specified
c  for any of the phases, the overall strain will always relax to
c  zero.  If it is desired to simply apply a fixed strain, with no
c  eigenstrains, then in subroutines Energy and Dembx, one must
c  zero out the elements of gb (in energy and in dembx) that
c  correspond to the macrostrains.  This is easily done.
c  This will fix the gradients of the macrostrains to always to be
c  zero, so that they will not change, so the applied strain (initial
c  values of the macrostrains) will remain fixed.

c  STRONGLY SUGGESTED:  READ MANUAL BEFORE USING PROGRAM!!!

c  (USER)  Change these dimensions and in other subroutines at same time.
c  For example, search and replace all occurrences throughout the
c  program of "(402" by "(1602", to go from a 20 x 20 system to
c  a 40 x 40 system.
c The u and similar arrays are dimensioned (nx times ny) + 2.

        real u(546914,2),gb(546914,2),b(546914,2)
        real h(546914,2),Ah(546914,2),T(546914,2)
	real C,dk(100,4,2,4,2)
        real cmod(100,3,3),ss(100,4,2),eigen(100,3)
	real zcon(2,2,2,2),pk(3,4,2)
	real phasemod(100,2),prob(100)
	integer in(9),jn(9),kn(9)
	integer*4 ib(546914,9)
	integer*2 pix(546914)

        common/list1/strxx,stryy,strxy
        common/list2/h,Ah
	common/list3/ib
        common/list4/pix
	common/list5/dk,b,C,zcon,Y
	common/list6/u
	common/list7/gb
        common/list8/cmod,T,eigen
        common/list10/phasemod,nphase,ss
        common/list11/sxx,syy,sxy
        common/list12/pmax
        common/list13/ncrack

c  (USER)  Unit 9 is the microstructure input file, unit 7
c  is the results output file.
c	open (9,file='micro.in')
	open (7,file='ab4-elas-mod.out')

c (USER) nx and ny are the size of the lattice
c note that the dimension of the arrays is (nx times ny) +2
        nx=844
        ny=648
c       ns=total number of sites
        ns=nx*ny
        write(7,9010) nx,ny,ns
9010  format(' nx= ',i4,'  ny= ',i4,' ns = ',i8)
        call srand(-189)

c  Add two more entries in the displacement vector for the 3 macrostrains,
c  u(ns+1,1) = exx, u(ns+1,2) = eyy, u(ns+2,1) = exy, u(ns+2,2) = never used
        nss=ns+2

c  (USER) nphase is the number of phases being considered in the problem.
c  The values of pix(m) will run from 1 to nphase.
c  add extra zero phase for cracks = 4.
	nphase=3+1+10

        call flush(7)
c  (USER) gtest is the stopping criterion, compared to gg=gb*gb.
c  If gtest=abc*ns, when gg < gtest, the rms value per pixel
c  of gb is less than sqrt(abc).
c  The value of gtest must be adjusted for each microstructure to
c  enable valid results.
        gtest=1.e-12*(nx*ny)
        write(7,*) 'relaxation criterion gtest = ',gtest

c  (USER)
c  The parameter phasemod(i,j) is the bulk (i,1)and shear (i,2) moduli of
c  the i'th phase. These can be
c  input in terms of Young's modulus E (i,1) and Poisson's ratio nu (i,2).
c  The program, in do 1144 loop, then changes them to bulk and shear
c  moduli, using relations for isotropic elastic moduli.
c  For anisotropic moduli tensors, one can directly input the whole tensor
c  cmod in subroutine femat, and skip this part.
c  If you wish to input in terms of bulk (i,1) and shear (i,2) moduli,
c  then simply comment out do 1144 loop.

cc mortar (moduli are in GPa). Poisson's ratio is unitless.
c	phasemod(1,1)=50.0
c	phasemod(1,2)=0.25
cc chert (expansive aggregate)
c	phasemod(2,1)=80.
c	phasemod(2,2)=0.2
cc regular aggregate
c	phasemod(3,1)=80.
c        phasemod(3,2)=0.2
cc cracked material
c	phasemod(4,1)=0.0
c        phasemod(4,2)=0.2
c        do 178 i=5,14
c	phasemod(i,1)=50.0+(rand(0)-0.5)*16.0
c	phasemod(i,2)=0.25
c178     continue

cc convert E and nu to K and G
c        do 1144 i=1,nphase
c        save=phasemod(i,1)
c        phasemod(i,1)=phasemod(i,1)/2./(1.-phasemod(i,2))
c        phasemod(i,2)=save/2./(1.+phasemod(i,2))
c1144    continue

      open(777,file="phasemod_values")
        do i=1,nphase
        read(777,*)phasemod(i,1),phasemod(i,2)
        end do
      close(777)

c  (USER) input eigen (thermal) strains for each phase.
c  (1=xx, 2=yy, 3=xy).
c   A positive value means that the phase wants to expand, a negative
c   value means the phase wants to shrink.
        eigen(1,1)=0.0
        eigen(1,2)=0.0
        eigen(1,3)=0.0
        eigen(2,1)=0.0
        eigen(2,2)=0.0
        eigen(2,3)=0.0
        eigen(3,1)=0.0
        eigen(3,2)=0.0
        eigen(3,3)=0.0
        eigen(4,1)=0.0
        eigen(4,2)=0.0
        eigen(4,3)=0.0
        do 47 i=5,14
        eigen(i,1)=eigen(1,1)
        eigen(i,2)=eigen(1,2)
        eigen(i,3)=eigen(1,3)
47      continue

c  Construct the 9 neighbor table, ib(m,n)

c  First construct the 9 neighbor table in terms of delta i and delta j
c  information (see Table 3 in manual)
      in(1)=0
      in(2)=1
      in(3)=1
      in(4)=1
      in(5)=0
      in(6)=-1
      in(7)=-1
      in(8)=-1
      in(9)=0

      jn(1)=1
      jn(2)=1
      jn(3)=0
      jn(4)=-1
      jn(5)=-1
      jn(6)=-1
      jn(7)=0
      jn(8)=1
      jn(9)=0

c  Now construct neighbor table according to 1-d labels
c  Matrix ib(m,n) gives the 1-d label of the n'th neighbor (n=1,9) of
c  the node labelled m.
      do 1020 j=1,ny
      do 1020 i=1,nx
      m=nx*(j-1)+i
      do 1004 n=1,9
      i1=i+in(n)
      j1=j+jn(n)
      if(i1.lt.1) i1=i1+nx
      if(i1.gt.nx) i1=i1-nx
      if(j1.lt.1) j1=j1+ny
      if(j1.gt.ny) j1=j1-ny
      m1=nx*(j1-1)+i1
      ib(m,n)=m1
1004  continue
1020  continue

c  Compute the average stress and strain, as well as the macrostrains (overall
c  system size and shape) in each microstructure.
c  (USER) npoints is the number of microstructures to use.
        npoints=1

c  Read in an initial microstructure in subroutine ppixel, and set up pix(m)
c  with the appropriate phase assignments.
        call ppixel(nx,ny,ns,nphase)
        call assig(ns,nphase,prob)
	do 8055 i=1,nphase
	write(7,9065) i,prob(i)
8055	continue
c  output elastic moduli (bulk and shear) for each phase

c make image (00.pgm) of uncracked microstructure
        call image(ns,0,nx,ny)
        do 8000 micro=1,npoints
c  Count and output the area fractions of the different phases
        call assig(ns,nphase,prob)
	do 8050 i=1,nphase
	write(7,9065) i,prob(i)
9065	format(' Area fraction of phase ',i3,'  is ',f10.8)
8050	continue
c  output elastic moduli (bulk and shear) for each phase
        write(7,*) '  Phase Moduli'
        do 111 i=1,nphase
        write(7,9020) i,phasemod(i,1),phasemod(i,2)
9020	format(' Phase ',i3,' bulk = ',f12.6,' shear = ',f12.6)
111	continue
c  output thermal strains for each phase
        write(7,*) '  Thermal Strains'
        do 119 i=1,nphase
        write(7,9029) i,eigen(i,1),eigen(i,2),eigen(i,3)
9029    format('Phase ',i3,'  ',3f6.2)
119	continue
        call flush(7)

c  (USER) Set inital macrostrains of computational cell
      u(ns+1,1)=0.01
      u(ns+1,2)=0.01
      u(nss,1)=0.01
      u(nss,2)=0.0
c Apply homogeneous macroscopic strain as the initial condition
c to displacement variables
	do 1050 j=1,ny
        do 1050 i=1,nx
		m=nx*(j-1)+i
		x=float(i-1)
		y=float(j-1)
                u(m,1)=x*u(ns+1,1)+y*u(nss,1)
                u(m,2)=x*u(nss,1)+y*u(ns+1,2)
1050	continue

c  Set up the finite element stiffness matrices,the constant, C,
c  the vector, b, required for the energy. b and C depend on the macrostrains.
c  When they are updated, the values of b and C are updated too via
c  calling subroutine femat.
c  Only compute the thermal strain terms the first time femat is called,
c  (iskip=0) as they are unaffected by later changes (iskip=1) in
c  displacements and macrostrains.
c  Compute initial value of gradient gb and gg=gb*gb.
        iskip=0
	call femat(nx,ny,ns,iskip)
        call energy(nx,ny,ns,utot)
 	gg=0.0
        do 100 m2=1,2
        do 100 m=1,nss
        gg=gg+gb(m,m2)*gb(m,m2)
100     continue
	write(7,9042) utot,gg
9042	format(' energy = ',e15.8,'  gg=  ',e15.8)
        call flush(7)

c  Relaxation loop
c  (USER) kmax is the maximum number of times that dembx will be called,
c  with ldemb conjugate gradient steps performed during each call.
c  The total number of conjugate gradient steps allowed for a given elastic
c  computation is kmax*ldemb.
        kmax=400
        ldemb=1000
        ltot=0
        ncrack=200

        do 5000 kkk=1,kmax

c  Call dembx to implement conjugate gradient routine
        write(7,*) 'Going into dembx, call no. ',kkk
        call dembx(nx,ny,ns,Lstep,gg,gtest,ldemb,kkk)
        ltot=ltot+Lstep
c  Call energy to compute energy after dembx call. If gg < gtest, this
c  will be the final energy.  If gg is still larger than gtest, then this
c  will give an intermediate energy with which to check how the
c  relaxation process is coming along. The call to energy does not
c  change the gradient or the value of gg.
c  Need to first call femat to update the vector b, as the value of the
c  components of b depend on the macrostrains.
        iskip=1
	call femat(nx,ny,ns,iskip)
        call energy(nx,ny,ns,utot)
	write(7,9043) utot,gg,ltot
9043    format(' energy = ',e15.8,' gg=  ',e15.8,' ltot = ',i6)
        call flush(7)

c  If relaxation process is finished, jump out of loop
        if(gg.lt.gtest) goto 444
c  Output stresses, strains, and macrostrains as an additional aid in judging
c  how well the relaxation process is proceeding.
      iswitch=0
      call stress(nx,ny,ns,iswitch)
      write(7,*) ' stresses:  xx,yy,xy'
      write(7,*) strxx,stryy,strxy
      write(7,*) ' strains:  xx,yy,xy'
      write(7,*) sxx,syy,sxy

      write(7,*) ' macrostrains in same order'
2745  format(2i8,3f15.8)
      call flush(7)
      call flush(8)
      write(7,*) 'avg = ',(u(ns+1,1)+u(ns+1,2))/2.

5000    continue

444   iswitch=1
      write(7,*) ' going into stress iswitch=1'
      call flush(7)
      write(7,*) u(ns+1,1),u(ns+1,2),u(ns+2,1)
      call stress(nx,ny,ns,iswitch)
      write(7,*) ' out of stress iswitch=1'
      call flush(7)
      write(7,*) ' stresses:  xx,yy,xy'
      write(7,*) strxx,stryy,strxy
      write(7,*) ' strains:  xx,yy,xy'
      write(7,*) sxx,syy,sxy
      call flush(7)

      write(7,*) ' macrostrains in same order'
      write(7,*) u(ns+1,1),u(ns+1,2),u(ns+2,1)
      write(7,*) 'avg = ',(u(ns+1,1)+u(ns+1,2))/2.
      call flush(7)

        bulk=(strxx+stryy)/(sxx+syy)/2.
        shear=strxy/sxy
        young=4.*bulk*shear/(bulk+shear)
        poisson=(bulk-shear)/(bulk+shear)
        write(7,*) ' bulk modulus = ',bulk
        write(7,*) ' shear modulus = ',shear
        write(7,*) ' Young modulus = ',young
        write(7,*) ' Poisson ratio = ',poisson
8000    continue

        end

c  Subroutine sets up the stiffness matrices, the linear term in the
c  regular displacements, b, and the constant term, C, which come from
c  the periodic boundary conditions, the term linear in the displacments,
c  T, that comes from the thermal strains, and the constant term Y.

      subroutine femat(nx,ny,ns,iskip)
      real u(546914,2),b(546914,2),T(546914,2)
      real dk(100,4,2,4,2),phasemod(100,2),dndx(4),dndy(4)
      real g(3,3),econ,ck(3,3),cmu(3,3),cmod(100,3,3)
      real es(3,4,2),zcon(2,2,2,2),ss(100,4,2)
      real eigen(100,3),delta(4,2)
      integer is(4),iskip
      integer*4 ib(546914,9)
      integer*2 pix(546914)

	common/list3/ib
        common/list4/pix
	common/list5/dk,b,C,zcon,Y
	common/list6/u
        common/list8/cmod,T,eigen
        common/list10/phasemod,nphase,ss

      nss=ns+2
c  Generate dk, zcon, T, and Y on first pass. After that they are
c  constant, since they are independent of the macrostrains.  Only b gets
c  upgraded as the macrostrains change.
c  Line number 1221 is the routine for b.
      if(iskip.eq.1) goto 1221

c  initialize stiffness matrices
      do 40 m=1,nphase
      do 40 l=1,2
      do 40 k=1,2
      do 40 j=1,4
      do 40 i=1,4
      dk(m,i,k,j,l)=0.0
40    continue
c  initialize zcon matrix (gives C term for arbitrary macrostrains)
      do 42 i=1,2
      do 42 j=1,2
      do 42 mi=1,2
      do 42 mj=1,2
      zcon(i,mi,j,mj)=0.0
42    continue
c  (USER) An anisotropic elastic moduli tensor could be input at this point,
c  bypassing the following code, which assumes isotropic elasticity
c  (only two independent numbers making up the elastic moduli
c  tensor, the bulk modulus K and the shear modulus G).

c  Set up elastic moduli matrices for each kind of element
c  ck and cmu are the bulk modulus and shear modulus matrices, which
c  need to multiplied by the actual bulk and shear moduli in each phase.

      ck(1,1)=1.0
      ck(1,2)=1.0
      ck(1,3)=0.0
      ck(2,1)=1.0
      ck(2,2)=1.0
      ck(2,3)=0.0
      ck(3,1)=0.0
      ck(3,2)=0.0
      ck(3,3)=0.0

      cmu(1,1)=1.0
      cmu(1,2)=-1.0
      cmu(1,3)=0.0
      cmu(2,1)=-1.0
      cmu(2,2)=1.0
      cmu(2,3)=0.0
      cmu(3,1)=0.0
      cmu(3,2)=0.0
      cmu(3,3)=1.0

      do 31 k=1,nphase
      do 21 j=1,3
      do 21 i=1,3
      cmod(k,i,j)=phasemod(k,1)*ck(i,j)+phasemod(k,2)*cmu(i,j)
21    continue
31    continue
c  Set up Simpson's integration rule weight vector
      do 30 j=1,3
      do 30 i=1,3
      nm=0
      if(i.eq.2) nm=nm+1
      if(j.eq.2) nm=nm+1
      g(i,j)=4.0**nm
30    continue

c  Loop over the nphase kinds of pixels and
c  Simpson's rule quadrature points in order to compute the stiffness
c  matrices.  Stiffness matrices of bilinear finite elements are quadratic
c  in x and y, so that Simpson's rule quadrature gives exact results.
      do 4000 ijk=1,nphase
      do 3000 j=1,3
      do 3000 i=1,3
      x=float(i-1)/2.0
      y=float(j-1)/2.0
c  dndx means the negative derivative with respect to x, of the shape
c  matrix N (see manual, Sec. 2.2), dndy is similar.
      dndx(1)=-(1.0-y)
      dndx(2)=(1.0-y)
      dndx(3)=y
      dndx(4)=-y
      dndy(1)=-(1.0-x)
      dndy(2)=-x
      dndy(3)=x
      dndy(4)=(1.0-x)
c  now build strain matrix
      do 2799 n1=1,3
      do 2799 n2=1,4
      do 2799 n3=1,2
      es(n1,n2,n3)=0.0
2799  continue
      do 2797 n=1,4
      es(1,n,1)=dndx(n)
      es(2,n,2)=dndy(n)
      es(3,n,1)=dndy(n)
      es(3,n,2)=dndx(n)
2797  continue
c  now do matrix multiply to determine value at (x,y), multiply by
c  proper weight, and sum into dk, the stiffness matrix
      do 900 mm=1,2
      do 900 nn=1,2
      do 900 ii=1,4
      do 900 jj=1,4
c  define sum over strain matrices and elastic moduli matrix for
c  stiffness matrix
      sum=0.0
      do 890 kk=1,3
      do 890 ll=1,3
      sum=sum+es(kk,ii,mm)*cmod(ijk,kk,ll)*es(ll,jj,nn)
890   continue
      dk(ijk,ii,mm,jj,nn)=dk(ijk,ii,mm,jj,nn)+g(i,j)*sum/36.
900   continue
3000  continue
4000  continue

c  Now compute the ss matrices, which give the thermal strain terms
c  for the i'th phase, single pixel.

      dndx(1)=-0.5
      dndx(2)=0.5
      dndx(3)=0.5
      dndx(4)=-0.5
      dndy(1)=-0.5
      dndy(2)=-0.5
      dndy(3)=0.5
      dndy(4)=0.5
c  now build average strain matrix
      do 3799 n1=1,3
      do 3799 n2=1,4
      do 3799 n3=1,2
      es(n1,n2,n3)=0.0
3799  continue
      do 3797 n=1,4
      es(1,n,1)=dndx(n)
      es(2,n,2)=dndy(n)
      es(3,n,1)=dndy(n)
      es(3,n,2)=dndx(n)
3797  continue
      do 3598 mmm=1,nphase
      do 3798 nn=1,2
      do 3798 mm=1,4
      sum=0.0
      do 3698 nm=1,3
      do 3698 n=1,3
      sum=sum+cmod(mmm,n,nm)*es(n,mm,nn)*eigen(mmm,nm)
3698  continue
      ss(mmm,mm,nn)=sum
3798  continue
3598  continue

c  now call subroutine const to generate zcon
c  zcon is a (2,2) x (2,2) matrix
      call const(dk,ns,zcon,nx,ny)

c  Now set up linear term, T, for thermal energy. It does not depend
c  on the macrostrains or displacements, so there is no need to update it
c  as the macrostrains change. T is built up out of the ss matrices.

      nss=ns+2
      do 6066 m2=1,2
      do 6066 m=1,nss
      T(m,m2)=0.0
6066  continue

c  For all cases, the correspondence between 1-4 finite element node
c  labels and the 1-9 neighbor labels is (see Table 4 in manual):
c  1:ib(m,9), 2:ib(m,3),3:ib(m,2),4:ib(m,1)
      is(1)=9
      is(2)=3
      is(3)=2
      is(4)=1
c  Do all points, but no macrostrain terms
c  note:  factor of 2 on linear thermal term is cancelled
c  by factor of 1/2 out in front of total energy term
      do 6601 j=1,ny
      do 6601 i=1,nx
      m=nx*(j-1)+i
      do 6600 mm=1,4
      do 6600 nn=1,2
      T(ib(m,is(mm)),nn)=T(ib(m,is(mm)),nn)-ss(pix(m),mm,nn)
6600  continue
6601  continue

c  now need to pick up and sum in all terms multiplying macrostrains
      do 7788 ipp=1,2
      do 7788 jpp=1,2
      if(ipp.eq.2.and.jpp.eq.2) goto 7788
      exx=0.0
      eyy=0.0
      exy=0.0
      if(ipp.eq.1.and.jpp.eq.1) exx=1.0
      if(ipp.eq.1.and.jpp.eq.2) eyy=1.0
      if(ipp.eq.2.and.jpp.eq.1) exy=1.0

c  x=nx face
      do 6001 i2=1,2
      do 6001 i4=1,4
      delta(i4,i2)=0.0
      if(i4.eq.2.or.i4.eq.3) then
      delta(i4,1)=exx*nx
      delta(i4,2)=exy*nx
      end if
6001  continue

      do 6000 j=1,ny-1
      m=j*nx
      do 6900 nn=1,2
      do 6900 mm=1,4
      T(ns+ipp,jpp)=T(ns+ipp,jpp)-ss(pix(m),mm,nn)*delta(mm,nn)
6900  continue
6000  continue

c  y=ny face
      do 6011 i2=1,2
      do 6011 i4=1,4
      delta(i4,i2)=0.0
      if(i4.eq.3.or.i4.eq.4) then
      delta(i4,1)=exy*ny
      delta(i4,2)=eyy*ny
      end if
6011  continue
      do 6010 i=1,nx-1
      m=nx*(ny-1)+i
      do 6901 nn=1,2
      do 6901 mm=1,4
      T(ns+ipp,jpp)=T(ns+ipp,jpp)-ss(pix(m),mm,nn)*delta(mm,nn)
6901  continue
6010  continue

c  x=nx y=ny corner
      do 6031 i2=1,2
      do 6031 i4=1,4
      delta(i4,i2)=0.0
      if(i4.eq.2) then
      delta(i4,1)=exx*nx
      delta(i4,2)=exy*nx
      end if
      if(i4.eq.4) then
      delta(i4,1)=exy*ny
      delta(i4,2)=eyy*ny
      end if
      if(i4.eq.3) then
      delta(i4,1)=exy*ny+exx*nx
      delta(i4,2)=eyy*ny+exy*nx
      end if
6031  continue
      m=nx*ny
      do 6903 nn=1,2
      do 6903 mm=1,4
      T(ns+ipp,jpp)=T(ns+ipp,jpp)-ss(pix(m),mm,nn)*delta(mm,nn)
6903  continue
6030  continue
7788  continue

c  now compute Y, the 0.5(eigen)Cij(eigen) energy, doesn't ever change
c  with macrostrain or displacements
      Y=0.0
      do 8811 m=1,ns
      do 8811 n=1,3
      do 8811 nn=1,3
      Y=Y+0.5*eigen(pix(m),n)*cmod(pix(m),n,nn)*eigen(pix(m),nn)
8811  continue

c  Following needs to be run after every change in macrostrain
c  when energy is recomputed.

1221  continue
c  Use auxiliary variables (exx, etc.) instead of u() variable, for
c  convenience, and to make the following code easier to read.
      exx=u(ns+1,1)
      eyy=u(ns+1,2)
      exy=u(nss,1)
c  Now set up vector for linear term that comes from periodic boundary
c  conditions.  Notation and conventions same as for T term.
c  This is done using the stiffness matrices, and the periodic terms
c  in the macrostrains.  It is easier to set up b in this way than to
c  analytically write out all the terms involved.

      do 5000 m2=1,2
      do 5000 m=1,ns
      b(m,m2)=0.0
5000  continue
c  For all cases, the correspondence between 1-4 finite element node
c  labels and the 1-9 neighbor labels is (see Table 4 in manual):
c  1:ib(m,9), 2:ib(m,3),3:ib(m,2),4:ib(m,1)
      is(1)=9
      is(2)=3
      is(3)=2
      is(4)=1

      C=0.0
c  x=nx face
      do 2001 i2=1,2
      do 2001 i4=1,4
      delta(i4,i2)=0.0
      if(i4.eq.2.or.i4.eq.3) then
      delta(i4,1)=exx*nx
      delta(i4,2)=exy*nx
      end if
2001  continue

      do 2000 j=1,ny-1
      m=j*nx
      do 1900 nn=1,2
      do 1900 mm=1,4
      sum=0.0
      do 1899 m2=1,2
      do 1899 m4=1,4
      sum=sum+delta(m4,m2)*dk(pix(m),m4,m2,mm,nn)
1899  continue
      b(ib(m,is(mm)),nn)=b(ib(m,is(mm)),nn)+sum
1900  continue
2000  continue
c  y=ny face
      do 2011 i2=1,2
      do 2011 i4=1,4
      delta(i4,i2)=0.0
      if(i4.eq.3.or.i4.eq.4) then
      delta(i4,1)=exy*ny
      delta(i4,2)=eyy*ny
      end if
2011  continue
      do 2010 i=1,nx-1
      m=nx*(ny-1)+i
      do 1901 nn=1,2
      do 1901 mm=1,4
      sum=0.0
      do 2099 m2=1,2
      do 2099 m4=1,4
      sum=sum+delta(m4,m2)*dk(pix(m),m4,m2,mm,nn)
2099  continue
      b(ib(m,is(mm)),nn)=b(ib(m,is(mm)),nn)+sum
1901  continue
2010  continue

c  x=nx y=ny corner
      do 2031 i2=1,2
      do 2031 i4=1,4
      delta(i4,i2)=0.0
      if(i4.eq.2) then
      delta(i4,1)=exx*nx
      delta(i4,2)=exy*nx
      end if
      if(i4.eq.4) then
      delta(i4,1)=exy*ny
      delta(i4,2)=eyy*ny
      end if
      if(i4.eq.3) then
      delta(i4,1)=exy*ny+exx*nx
      delta(i4,2)=eyy*ny+exy*nx
      end if
2031  continue
      m=nx*ny
      do 1903 nn=1,2
      do 1903 mm=1,4
      sum=0.0
      do 2029 m2=1,2
      do 2029 m4=1,4
      sum=sum+delta(m4,m2)*dk(pix(m),m4,m2,mm,nn)
2029  continue
      b(ib(m,is(mm)),nn)=b(ib(m,is(mm)),nn)+sum
1903  continue
2030  continue

      return
      end

c  Subroutine that computes derivatives of the b-vector with respect
c  to the macrostrains.  Since b is linear in the macrostrains, the
c  derivative with respect to any one of them can be computed simply
c  by letting that macrostrain, within the subroutine, be equal to one,
c  and all the other macrostrains to be zero.
c  Very similar to 1221 loop in femat for b.

      subroutine bgrad(nx,ny,ns,exx,eyy,exy)
      real b(546914,2)
      real dk(100,4,2,4,2),delta(4,2),zcon(2,2,2,2)
      integer is(4)
      integer*4 ib(546914,9)
      integer*2 pix(546914)

      common/list3/ib
      common/list4/pix
      common/list5/dk,b,C,zcon,Y

c  exx, eyy, and exy are the artificial macrostrains used
c  to get the gradient terms (appropriate combinations of 1's and 0's).

c  Set up vector for linear term

      do 5000 m2=1,2
      do 5000 m=1,ns
      b(m,m2)=0.0
5000  continue
      is(1)=9
      is(2)=3
      is(3)=2
      is(4)=1

c  x=nx face
      do 2001 i2=1,2
      do 2001 i4=1,4
      delta(i4,i2)=0.0
      if(i4.eq.2.or.i4.eq.3) then
      delta(i4,1)=exx*nx
      delta(i4,2)=exy*nx
      end if
2001  continue

      do 2000 j=1,ny-1
      m=j*nx
      do 1900 nn=1,2
      do 1900 mm=1,4
      sum=0.0
      do 1899 m2=1,2
      do 1899 m4=1,4
      sum=sum+delta(m4,m2)*dk(pix(m),m4,m2,mm,nn)
1899  continue
      b(ib(m,is(mm)),nn)=b(ib(m,is(mm)),nn)+sum
1900  continue
2000  continue
c  y=ny face
      do 2011 i2=1,2
      do 2011 i4=1,4
      delta(i4,i2)=0.0
      if(i4.eq.3.or.i4.eq.4) then
      delta(i4,1)=exy*ny
      delta(i4,2)=eyy*ny
      end if
2011  continue
      do 2010 i=1,nx-1
      m=nx*(ny-1)+i
      do 1901 nn=1,2
      do 1901 mm=1,4
      sum=0.0
      do 2099 m2=1,2
      do 2099 m4=1,4
      sum=sum+delta(m4,m2)*dk(pix(m),m4,m2,mm,nn)
2099  continue
      b(ib(m,is(mm)),nn)=b(ib(m,is(mm)),nn)+sum
1901  continue
2010  continue

c  x=nx y=ny corner
      do 2031 i2=1,2
      do 2031 i4=1,4
      delta(i4,i2)=0.0
      if(i4.eq.2) then
      delta(i4,1)=exx*nx
      delta(i4,2)=exy*nx
      end if
      if(i4.eq.4) then
      delta(i4,1)=exy*ny
      delta(i4,2)=eyy*ny
      end if
      if(i4.eq.3) then
      delta(i4,1)=exy*ny+exx*nx
      delta(i4,2)=eyy*ny+exy*nx
      end if
2031  continue
      m=nx*ny
      do 1903 nn=1,2
      do 1903 mm=1,4
      sum=0.0
      do 2029 m2=1,2
      do 2029 m4=1,4
      sum=sum+delta(m4,m2)*dk(pix(m),m4,m2,mm,nn)
2029  continue
      b(ib(m,is(mm)),nn)=b(ib(m,is(mm)),nn)+sum
1903  continue
2030  continue

      return
      end

c  Subroutine computes the quadratic term in the macrostrains, that comes
c  from the periodic boundary conditions, and sets it up as a
c  (2,2) x (2,2) matrix that couples to the three macrostrains

      subroutine const(dk,ns,zcon,nx,ny)
      real dk(100,4,2,4,2),zcon(2,2,2,2),delta(4,2)
      real pp(4,4),s(4,4)
      integer*2 pix(546914)
      common/list4/pix

c  routine to set up 4 x 4 matrix for energy term involving macro-strains
c  only, pulled out of femat (really 3 x 3, as ns+2,2 term is not used).
c  Idea is to evaluate the quadratic term in the macrostrains, C,
c  repeatedly for choices of strain like
c  exx=1, exy=1, all others = 0, build up 10 choices, then recombine
c  to get matrix elements by themselves

      nss=ns+2

      do 1111 i=1,4
      do 1111 j=1,4
      s(i,j)=0.0
      pp(i,j)=0.0
1111  continue

      do 5000 ii=1,4
      do 5000 jj=ii,4
      econ=0.0
      exx=0.0
      eyy=0.0
      exy=0.0
      if(ii.eq.1.and.jj.eq.1) exx=1.0
      if(ii.eq.2.and.jj.eq.2) eyy=1.0
      if(ii.eq.3.and.jj.eq.3) exy=1.0
      if(ii.eq.1.and.jj.eq.2) then
      exx=1.0
      eyy=1.0
      end if
      if(ii.eq.1.and.jj.eq.3) then
      exx=1.0
      exy=1.0
      end if
      if(ii.eq.2.and.jj.eq.3) then
      eyy=1.0
      exy=1.0
      end if
c  x=nx face
      do 2001 i2=1,2
      do 2001 i4=1,4
      delta(i4,i2)=0.0
      if(i4.eq.2.or.i4.eq.3) then
      delta(i4,1)=exx*nx
      delta(i4,2)=exy*nx
      end if
2001  continue

      do 2000 j=1,ny-1
      m=j*nx
      do 1900 nn=1,2
      do 1900 mm=1,4
      do 1899 m2=1,2
      do 1899 m4=1,4
      econ=econ+0.5*delta(m4,m2)*dk(pix(m),m4,m2,mm,nn)*delta(mm,nn)
1899  continue
1900  continue
2000  continue
c  y=ny face
      do 2011 i2=1,2
      do 2011 i4=1,4
      delta(i4,i2)=0.0
      if(i4.eq.3.or.i4.eq.4) then
      delta(i4,1)=exy*ny
      delta(i4,2)=eyy*ny
      end if
2011  continue
      do 2010 i=1,nx-1
      m=nx*(ny-1)+i
      do 1901 nn=1,2
      do 1901 mm=1,4
      do 2099 m2=1,2
      do 2099 m4=1,4
      econ=econ+0.5*delta(m4,m2)*dk(pix(m),m4,m2,mm,nn)*delta(mm,nn)
2099  continue
1901  continue
2010  continue

c  x=nx y=ny corner
      do 2031 i2=1,2
      do 2031 i4=1,4
      delta(i4,i2)=0.0
      if(i4.eq.2) then
      delta(i4,1)=exx*nx
      delta(i4,2)=exy*nx
      end if
      if(i4.eq.4) then
      delta(i4,1)=exy*ny
      delta(i4,2)=eyy*ny
      end if
      if(i4.eq.3) then
      delta(i4,1)=exy*ny+exx*nx
      delta(i4,2)=eyy*ny+exy*nx
      end if
2031  continue
      m=nx*ny
      do 1903 nn=1,2
      do 1903 mm=1,4
      do 2029 m2=1,2
      do 2029 m4=1,4
      econ=econ+0.5*delta(m4,m2)*dk(pix(m),m4,m2,mm,nn)*delta(mm,nn)
2029  continue
1903  continue
2030  continue
      pp(ii,jj)=econ*2.

5000  continue
      do 6000 i=1,4
      do 6000 j=i,4
      if(i.eq.j) s(i,j)=pp(i,j)
      if(i.ne.j) then
      s(i,j)=pp(i,j)-pp(i,i)-pp(j,j)
      end if
6000  continue
      do 7000 i=1,4
      do 7000 j=1,4
      pp(i,j)=0.5*(s(i,j)+s(j,i))
7000  continue
c  now map pp(i,j) into zcon(2,2,2,2)
      do 7200 i=1,2
      do 7200 j=1,2
      do 7200 mi=1,2
      do 7200 mj=1,2
      if(i.eq.1) ii=i+mi-1
      if(i.eq.2) ii=i+mi
      if(j.eq.1) jj=j+mj-1
      if(j.eq.2) jj=j+mj
      zcon(i,mi,j,mj)=pp(ii,jj)
      if(i.eq.2.and.mi.eq.2) zcon(i,mi,j,mj)=0.0
      if(j.eq.2.and.mj.eq.2) zcon(i,mi,j,mj)=0.0
7200  continue

      return
      end

c  Subroutine computes the total energy, utot, and the gradient, gb,
c  for the regular displacements as well as for the macrostrains

      subroutine energy(nx,ny,ns,utot)

	real u(546914,2),gb(546914,2)
	real b(546914,2),T(546914,2)
	real cmod(100,3,3),pk(3,4,2),eigen(100,3)
	real dk(100,4,2,4,2),zcon(2,2,2,2),C
 	integer*4 ib(546914,9)
 	integer*2 pix(546914)

	common/list3/ib
        common/list4/pix
	common/list5/dk,b,C,zcon,Y
	common/list6/u
	common/list7/gb
        common/list8/cmod,T,eigen

        nss=ns+2

        do 2090 m2=1,2
	do 2090 m=1,nss
	gb(m,m2)=0.0
2090	continue
c  Do global matrix multiply via small stiffness matrices, Ah = (A)(h).
c  The long statement below correctly brings in all the terms from
c  the global matrix A using only the small stiffness matrices.

        do 3000 j=1,2
        do 3000 n=1,2
	do 3000 m=1,ns
      gb(m,j)=gb(m,j)+u(ib(m,1),n)*( dk(pix(ib(m,9)),1,j,4,n)
     &+dk(pix(ib(m,7)),2,j,3,n) )+
     &u(ib(m,2),n)*( dk(pix(ib(m,9)),1,j,3,n) )+
     &u(ib(m,3),n)*( dk(pix(ib(m,9)),1,j,2,n)+dk(pix(ib(m,5)),4,j,3,n))+
     &u(ib(m,4),n)*( dk(pix(ib(m,5)),4,j,2,n) )+
     &u(ib(m,5),n)*( dk(pix(ib(m,6)),3,j,2,n)+dk(pix(ib(m,5)),4,j,1,n))+
     &u(ib(m,6),n)*( dk(pix(ib(m,6)),3,j,1,n) )+
     &u(ib(m,7),n)*( dk(pix(ib(m,6)),3,j,4,n)+dk(pix(ib(m,7)),2,j,1,n))+
     &u(ib(m,8),n)*( dk(pix(ib(m,7)),2,j,4,n) )+
     &u(ib(m,9),n)*( dk(pix(ib(m,9)),1,j,1,n)
     &+dk(pix(ib(m,7)),2,j,2,n)+
     &dk(pix(ib(m,6)),3,j,3,n)+dk(pix(ib(m,5)),4,j,4,n) )
3000	continue

      utot=0.0

      gtot=0.0
      do 3100 m2=1,2
      do 3100 m=1,ns
      utot=utot+0.5*u(m,m2)*gb(m,m2)+b(m,m2)*u(m,m2)
c  this is gradient of energy with respect to normal displacements
      gb(m,m2)=gb(m,m2)+b(m,m2)
3100  continue

c  compute "constant" macrostrain energy term
      C=0.0
      do 7200 i=1,2
      do 7200 j=1,2
      do 7200 mi=1,2
      do 7200 mj=1,2
      C=C+0.5*u(ns+i,mi)*zcon(i,mi,j,mj)*u(ns+j,mj)
7200  continue
      utot=utot+C
c  now add in constant term from thermal energy, Y
      utot=utot+Y
c  now add in linear term in thermal energy
      do 7171 m2=1,2
      do 7171 m=1,nss
      utot=utot+T(m,m2)*u(m,m2)
7171  continue

c  now compute gradient with respect to macrostrains
c  put in piece from first derivative of zcon quadratic term
      do 7300 i=1,2
      do 7300 mi=1,2
      sum=0.0
      do 7250 j=1,2
      do 7250 mj=1,2
      sum=sum+zcon(i,mi,j,mj)*u(ns+j,mj)
7250  continue
      gb(ns+i,mi)=sum
7300  continue
c  add in piece of gradient, for displacements as well as macrostrains,
c  that come from linear term in thermal energy
      do 3150 m2=1,2
      do 3150 m=1,nss
      gb(m,m2)=gb(m,m2)+T(m,m2)
3150  continue
c  now generate part that comes from b . u term
c  do by calling b generation with appropriate macrostrain set to 1 to
c  get that partial derivative, just use bgrad (taken from femat),
c  skip dk and zcon part
      do 8100 ii=1,3
      exx=0.0
      eyy=0.0
      exy=0.0
      if(ii.eq.1) exx=1.0
      if(ii.eq.2) eyy=1.0
      if(ii.eq.3) exy=1.0
      call bgrad(nx,ny,ns,exx,eyy,exy)
      sum=0.0
        do 8200 m2=1,2
	do 8200 m=1,ns
	sum=sum+u(m,m2)*b(m,m2)
8200	continue
	if(ii.eq.1) gb(ns+1,1)=gb(ns+1,1)+sum
	if(ii.eq.2) gb(ns+1,2)=gb(ns+1,2)+sum
	if(ii.eq.3) gb(ns+2,1)=gb(ns+2,1)+sum
8100    continue
c zero out macrostrain part of gradient
        gb(ns+1,1)=0.0
        gb(ns+1,2)=0.0
        gb(ns+2,1)=0.0
        gb(ns+2,2)=0.0
        return
        end

c  Subroutine that carries out the conjugate gradient relaxation process

      subroutine dembx(nx,ny,ns,Lstep,gg,gtest,ldemb,kkk)

      real u(546914,2),gb(546914,2),b(546914,2)
      real h(546914,2),Ah(546914,2)
      real dk(100,4,2,4,2),zcon(2,2,2,2)
      real lambda,gamma
      integer*4 ib(546914,9)
      integer*2 pix(546914)

        common/list2/h,Ah
	common/list3/ib
        common/list4/pix
	common/list5/dk,b,C,zcon,Y
	common/list6/u
	common/list7/gb

      nss=ns+2
c  Initialize the conjugate direction vector on first call to dembx only.
c  For calls to dembx after the first, we want to continue using the value
c  of h determined in the previous call.  Of course, if npoints is greater
c  than 1, then this initialization step will be run each time a new
c  microstructure is used, as kkk will be reset to 1 every time the
c  counter micro is increased.
      if(kkk.eq.1) then
      do 500 m2=1,2
      do 500 m=1,nss
      h(m,m2)=gb(m,m2)
500   continue
      end if
c  Lstep counts the number of conjugate gradient steps taken
c  in each call to dembx
      Lstep=0

      do 800 ijk=1,ldemb
      Lstep=Lstep+1

      do 290 m2=1,2
      do 290 m=1,nss
      Ah(m,m2)=0.0
290   continue
c  Do global matrix multiply via small stiffness matrices, Ah = (A)(h).
c  The long statement below correctly brings in all the terms from
c  the global matrix A using only the small stiffness matrices.
        do 400 j=1,2
        do 400 n=1,2
	do 400 m=1,ns
      Ah(m,j)=Ah(m,j)+h(ib(m,1),n)*( dk(pix(ib(m,9)),1,j,4,n)
     &+dk(pix(ib(m,7)),2,j,3,n) )+
     &h(ib(m,2),n)*( dk(pix(ib(m,9)),1,j,3,n) )+
     &h(ib(m,3),n)*( dk(pix(ib(m,9)),1,j,2,n)+dk(pix(ib(m,5)),4,j,3,n))+
     &h(ib(m,4),n)*( dk(pix(ib(m,5)),4,j,2,n) )+
     &h(ib(m,5),n)*( dk(pix(ib(m,6)),3,j,2,n)+dk(pix(ib(m,5)),4,j,1,n))+
     &h(ib(m,6),n)*( dk(pix(ib(m,6)),3,j,1,n) )+
     &h(ib(m,7),n)*( dk(pix(ib(m,6)),3,j,4,n)+dk(pix(ib(m,7)),2,j,1,n))+
     &h(ib(m,8),n)*( dk(pix(ib(m,7)),2,j,4,n) )+
     &h(ib(m,9),n)*( dk(pix(ib(m,9)),1,j,1,n)
     &+dk(pix(ib(m,7)),2,j,2,n)+
     &dk(pix(ib(m,6)),3,j,3,n)+dk(pix(ib(m,5)),4,j,4,n) )
400    continue

c  The above accurately gives the second derivative matrix with respect
c  to nodal displacements, but fails to give the 2nd derivative terms that
c  include the macrostrains [ du d(strain) and d(strain)d(strain) ].
c  Use repeated calls to bgrad to generate mixed 2nd derivatives terms,
c  plus use zcon in order to correct the matrix multiply and correctly bring
c  in macrostrain terms (see manual, Sec. 2.4).
      do 8100 ii=1,3
      e11=0.0
      e22=0.0
      e12=0.0
      if(ii.eq.1) e11=1.0
      if(ii.eq.2) e22=1.0
      if(ii.eq.3) e12=1.0
      call bgrad(nx,ny,ns,e11,e22,e12)
c  now fill in terms from matrix multiply
c  right hand sides, 1 to ns
      do 3333 m=1,ns
      do 3333 m1=1,2
      if(ii.eq.1) Ah(m,m1)=Ah(m,m1)+b(m,m1)*h(ns+1,1)
      if(ii.eq.2) Ah(m,m1)=Ah(m,m1)+b(m,m1)*h(ns+1,2)
      if(ii.eq.3) Ah(m,m1)=Ah(m,m1)+b(m,m1)*h(nss,1)
3333  continue
c  now do across bottom, 1 to ns
      do 3334 m=1,ns
      if(ii.eq.1) Ah(ns+1,1)=Ah(ns+1,1)+b(m,1)*h(m,1)+
     +b(m,2)*h(m,2)
      if(ii.eq.2) Ah(ns+1,2)=Ah(ns+1,2)+b(m,1)*h(m,1)+
     +b(m,2)*h(m,2)
      if(ii.eq.3) Ah(nss,1)=Ah(nss,1)+b(m,1)*h(m,1)+
     +b(m,2)*h(m,2)
3334  continue
c  now do righthand corner terms, ns+1 to nss
      do 3335 m=1,2
      do 3335 m1=1,2
      if(ii.eq.1) Ah(ns+1,1)=Ah(ns+1,1)+zcon(1,1,m,m1)*h(ns+m,m1)
      if(ii.eq.2) Ah(ns+1,2)=Ah(ns+1,2)+zcon(1,2,m,m1)*h(ns+m,m1)
      if(ii.eq.3) Ah(nss,1)=Ah(nss,1)+zcon(2,1,m,m1)*h(ns+m,m1)
3335  continue

8100    continue

c zero out macrostrain part of Ah
      Ah(ns+1,1)=0.0
      Ah(ns+1,2)=0.0
      Ah(ns+2,1)=0.0
      Ah(ns+2,2)=0.0

      hAh=0.0
      do 530 m2=1,2
      do 530 m=1,nss
      hAh=hAh+h(m,m2)*Ah(m,m2)
530   continue

      lambda=gg/hAh
      do 540 m2=1,2
      do 540 m=1,nss
      u(m,m2)=u(m,m2)-lambda*h(m,m2)
      gb(m,m2)=gb(m,m2)-lambda*Ah(m,m2)
540   continue

      gglast=gg
      gg=0.0
      do 550 m2=1,2
      do 550 m=1,nss
      gg=gg+gb(m,m2)*gb(m,m2)
550   continue
      if(gg.lt.gtest) goto 1000

      gamma=gg/gglast
      do 570 m2=1,2
      do 570 m=1,nss
      h(m,m2)=gb(m,m2)+gamma*h(m,m2)
570   continue

800   continue

1000  continue
      return
      end

c  Subroutine that computes the three average stresses
c  and three average strains.

      subroutine stress(nx,ny,ns,iswitch)

      real u(546914,2),uu(4,2),pmax(546914,2)
      real T(546914,2),eigen(100,3)
      real dndx(4),dndy(4),es(3,4,2),cmod(100,3,3)
      integer*4 ib(546914,9)
      integer*2 pix(546914)

      common/list1/strxx,stryy,strxy
      common/list3/ib
      common/list4/pix
      common/list6/u
      common/list8/cmod,T,eigen
      common/list11/sxx,syy,sxy
      common/list12/pmax

      nss=ns+2
      exx=u(ns+1,1)
      eyy=u(ns+1,2)
      exy=u(nss,1)
      if(iswitch.eq.1) then
c      open(unit=12,file='101bw-stress-strain-allmat.txt')
c      write(12,1134) nx,ny
c1134  format(2i5)
      end if

c set up single pixel strain matrix

      dndx(1)=-0.5
      dndx(2)=0.5
      dndx(3)=0.5
      dndx(4)=-0.5
      dndy(1)=-0.5
      dndy(2)=-0.5
      dndy(3)=0.5
      dndy(4)=0.5

c initialize pmax for iswitch=1
      if(iswitch.eq.1) then
      do 5555 m=1,ns
      pmax(m,1)=0.0
      pmax(m,2)=m
5555  continue
      end if
c  Build average strain matrix, follows code in femat, but for average
c  strain over a pixel, not the strain at a point
      do 2799 n1=1,3
      do 2799 n2=1,4
      do 2799 n3=1,2
      es(n1,n2,n3)=0.0
2799  continue
      do 2797 n=1,4
      es(1,n,1)=dndx(n)
      es(2,n,2)=dndy(n)
      es(3,n,1)=dndy(n)
      es(3,n,2)=dndx(n)
2797  continue
c  now compute average stresses and strains in each pixel
      sxx=0.0
      syy=0.0
      sxy=0.0
      strxx=0.0
      stryy=0.0
      strxy=0.0
      do 470 j=1,ny
      do 470 i=1,nx
      m=(j-1)*nx+i
      if(iswitch.eq.1) then
c      write(14,*) ' i j in stress m ',i,j,m
c      call flush(14)
      end if
c  load in elements of 4-vector using pd. bd. conds.
      do 9898 mm=1,2
      uu(1,mm)=u(m,mm)
      uu(2,mm)=u(ib(m,3),mm)
      uu(3,mm)=u(ib(m,2),mm)
      uu(4,mm)=u(ib(m,1),mm)
9898  continue
c  Correct for periodic boundary conditions, some displacements are wrong
c  for a pixel on a periodic boundary.  Since they come from an opposite
c  face, need to put in applied strain to correct them.
      if(i.eq.nx) then
      uu(2,1)=uu(2,1)+exx*nx
      uu(2,2)=uu(2,2)+exy*nx
      uu(3,1)=uu(3,1)+exx*nx
      uu(3,2)=uu(3,2)+exy*nx
      end if
      if(j.eq.ny) then
      uu(3,1)=uu(3,1)+exy*ny
      uu(3,2)=uu(3,2)+eyy*ny
      uu(4,1)=uu(4,1)+exy*ny
      uu(4,2)=uu(4,2)+eyy*ny
      end if
      str11=0.0
      str22=0.0
      str12=0.0
      s11=0.0
      s22=0.0
      s12=0.0
c********compute average stress and strain tensor in each pixel*************
c  First put thermal strain-induced stresses into stress tensor
      do 465 n=1,3
      str11=str11-cmod(pix(m),1,n)*eigen(pix(m),n)
      str22=str22-cmod(pix(m),2,n)*eigen(pix(m),n)
      str12=str12-cmod(pix(m),3,n)*eigen(pix(m),n)
465   continue
      do 466 n2=1,2
      do 466 n4=1,4
c  compute non-thermal strains in each pixel
      s11=s11+es(1,n4,n2)*uu(n4,n2)
      s22=s22+es(2,n4,n2)*uu(n4,n2)
      s12=s12+es(3,n4,n2)*uu(n4,n2)
      do 466 n=1,3
c  compute stresses in each pixel that include both non-thermal
c  and thermal strains
      str11=str11+cmod(pix(m),1,n)*es(n,n4,n2)*uu(n4,n2)
      str22=str22+cmod(pix(m),2,n)*es(n,n4,n2)*uu(n4,n2)
      str12=str12+cmod(pix(m),3,n)*es(n,n4,n2)*uu(n4,n2)
466   continue
c  Sum local stresses and strains into global stresses and strains
      strxx=strxx+str11
      stryy=stryy+str22
      strxy=strxy+str12
      sxx=sxx+s11
      syy=syy+s22
      sxy=sxy+s12
      if(iswitch.eq.1) then
      sa=0.5*((str11+str22)+sqrt((str11-str22)**2+4.*str12*str12))
      sb=0.5*((str11+str22)-sqrt((str11-str22)**2+4.*str12*str12))
      ea=0.5*((s11+s22)+sqrt((s11-s22)**2+4.*s12*s12))
      eb=0.5*((s11+s22)-sqrt((s11-s22)**2+4.*s12*s12))
      if(ea.gt.eb) pmax(m,1)=ea
      if(eb.gt.ea) pmax(m,1)=eb
      pmax(m,2)=m
c      write(14,*) m,pmax(m,1),pmax(m,2)
c      call flush(14)
c      write(12,1137) i,j,str11,str22,str12,s11,s22,s12
c      write(12,1138) i,j,sa,sb,ea,eb
1137  format(2i4,6f15.9)
1138  format(2i4,4f15.9)
      end if
470   continue

c  Area average global stresses and strains

	strxx=strxx/float(ns)
	stryy=stryy/float(ns)
	strxy=strxy/float(ns)
	sxx=sxx/float(ns)
	syy=syy/float(ns)
	sxy=sxy/float(ns)

      return
      end

c  Subroutine to count volume fractions of various phases

      subroutine assig(ns,nphase,prob)
      integer*2 pix(546914)
      real prob(100)
      common/list4/pix

	do 999 i=1,nphase
	prob(i)=0.0
999	continue

      do 1000 m=1,ns
      do 1000 i=1,nphase
        if(pix(m).eq.i) then
	prob(i)=prob(i)+1
	end if
1000    continue

	do 998 i=1,nphase
	prob(i)=prob(i)/float(ns)
998	continue

          return
          end

c  Subroutine to set up image of microstructure
c  adds border all around of mortar (phase 1)

      subroutine ppixel(nx,ny,ns,nphase)
      integer*2 pix(546914)
      common/list4/pix

      call srand(-147)
c  (USER)  If you want to set up a test image inside the program, instead
c  of reading it in from a file, this should be done inside this subroutine.

      open (9,file='micro.in2')

      do 200 j=1,ny
      do 200 i=1,nx
      m=nx*(j-1)+i
      read(9,*) pix(m)
c turn white border into black matrix
!      if(pix(m).eq.4) pix(m)=1
200   continue
      close(9)

c make phases 5-14 matrix phases but with slightly different moduli

c      do 300 j=1,ny
c     do 300 i=1,nx
c     m=nx*(j-1)+i
c     if(pix(m).eq.1) then
c     mm=5+10*rand(0)
c     if(mm.eq.15) mm=14
c     pix(m)=mm
c     end if
c00   continue
c     close(19)

c   Check for wrong phase labels--less than 1 or greater than nphase
       do 500 m=1,ns
      if(pix(m).lt.1) then
       write(7,*) 'Phase label in pix < 1--error at ',m,pix(m)
      end if
       if(pix(m).gt.nphase) then
        write(7,*) 'Phase label in pix > nphase--error at ',m,pix(m)
       end if
500    continue

      return
      end
      subroutine image(ns,micro,nx,ny)
      integer*2 pix(546914)
      common/list4/pix
      character*6 namef
      character*9 namefff
      character*31 nameff

      namef(1:2)='00'
      namef(3:6)='.pgm'
      if(micro.lt.10) write(namef(2:2),'(I1)') micro
      if(micro.ge.10) write(namef(1:2),'(I2)') micro
      open(unit=14,file=namef)
      write(14,133)
      write(14,134) nx,ny
      write(14,157)
133   format('P2')
134   format(i4,1x,i4)
157   format('4')

      do 741 j=ny,1,-1
      do 741 i=1,nx
      m=nx*(j-1)+i
      mm=pix(m)
c make irregular matrix all look the same in the image
      if(mm.ge.5.and.mm.le.14) mm=1
      write(14,22) mm
22    format(i1)
741   continue
      close(14)

c convert pgm file to gif

      nameff(1:8)='convert '
      nameff(9:14)=namef
      nameff(15:24)=' -depth 8 '
      nameff(25:30)=namef
      nameff(28:30)='gif'
      nameff(31:31)=' '
      call system(nameff)
c rm pgm files
      namefff(1:3)='rm '
      namefff(4:9)=namef
      call system(namefff)

      return
      end
