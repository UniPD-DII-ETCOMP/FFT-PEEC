#include "fintrf.h"      
subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    implicit none
! mwPointer mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
! mwSize stuff for mexing
      mwSize mo,no,siz
! mwPointer stuff for mexing
      mwPointer mxCreateDoubleMatrix
      mwPointer mxGetPr
	  mwPointer L_pr,MM_pr,N_pr,xyz_src_pr,Np_pr,face_pr
	  mwPointer r_f_pr,r_e_pr,u_e_pr,triang_face_pr,idxf_dot_pr
      mwPointer idxf_cross_pr,w_pr,d_pr, N_thread_pr,G_pr
      mwPointer m, n, s
      mwPointer mxGetM, mxGetN, mxGetDimensions
	  mwSize mxGetNumberOfDimensions
      mwPointer mxGetNumberOfElements
      mwPointer pm
!integer (normal not 4/8)
      integer flag,i,j,ii,jj,icount
      integer mxIsNumeric 
      integer*4 ComplexFlag
! fortran subroutine arguments
      real*8,allocatable,dimension(:,:,:,:) ::  triang_face
      real*8,allocatable,dimension(:,:,:) ::  face, r_e, u_e, idxf_cross_re
	  integer*8,allocatable,dimension(:,:,:) :: idxf_cross
	  real*8,allocatable,dimension(:,:) :: xyz_src, r_f, idxf_dot_re
	  integer*8,allocatable,dimension(:,:) :: idxf_dot
	  real*8,allocatable,dimension(:) :: w, G, d
	  real*8 L_re, MM_re, N_re, Np_re, N_thread_re
      integer*8 L, MM, NN, Np, N_thread
	  mwPointer dims(3)
	  mwPointer dims4(4)
	  character*80 msg
	  logical debu
	  debu=.true.
	  !!if(debu) open(unit=66,file='log_Green.txt',status='unknown')
! check for proper number of arguments. 
      if (nrhs .ne. 15) then
        call mexErrMsgIdAndTxt ('MATLAB:Green:nInput', &
                                '15 input argument required.')
      elseif (nlhs .ne. 1) then
        call mexErrMsgIdAndTxt ('MATLAB:Green:nOutput', &
                                '1 output argument required.')
      endif
      !!if(debu) !write(66,*) 'arguments checked'
	  
!    Check to see inputs are numeric.
      do i=1,15
        if (mxIsNumeric(prhs(i)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:Green:NonNumeric', &
                                'Inputs must be numeric.')
        endif
	  enddo
	  !!if(debu) !write(66,*) 'numeric checked' 
	  
	  
!     Check that input #1 is scalar and fetch it
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 1 must be scalar, L.')
      endif	  
      siz = m*n
      L_pr = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(L_pr, L_re, siz)
      L=int(L_re,8)
      !!if(debu) !write(66,*) 'input1 checked, L=',L


!     Check that input #2 is scalar and fetch it
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 1 must be scalar, M.')
      endif	  
      siz = m*n
      MM_pr = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(MM_pr, MM_re, siz)
      MM=int(MM_re,8)
      !!if(debu) !write(66,*) 'input2 checked, M=',MM

!     Check that input #3 is scalar and fetch it
      m = mxGetM(prhs(3))
      n = mxGetN(prhs(3))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 3 must be scalar, N.')
      endif	  
      siz = m*n
      N_pr = mxGetPr(prhs(3))
      call mxCopyPtrToReal8(N_pr, N_re, siz)
      NN=int(N_re,8)
      !if(debu) !write(66,*) 'input3 checked, N=',NN


!     Check that input #5 is scalar and fetch it
      m = mxGetM(prhs(5))
      n = mxGetN(prhs(5))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 5 must be scalar, Np.')
      endif	  
      siz = m*n
      Np_pr = mxGetPr(prhs(5))
      call mxCopyPtrToReal8(Np_pr, Np_re, siz)
      Np=int(Np_re,8)
      !if(debu) !write(66,*) 'input5 checked, Np=',Np


!     Check that input #4and fetch it
      m = mxGetM(prhs(4))
      n = mxGetN(prhs(4))
      if(m .ne. Np .or. n .ne. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 4 must be xyz_src_pr Npx3.')
      endif
      siz = m*n
      xyz_src_pr = mxGetPr(prhs(4))
	  allocate(xyz_src(np,3))
      call mxCopyPtrToReal8(xyz_src_pr, xyz_src, siz)
	  !if(debu) !write(66,*) 'input4 checked xyz_src = ', xyz_src


!     Check that input #6and fetch it
	  !if(debu) !write(66,*) 'mxGetNumberOfDimensions(prhs(6))'
	  !if(debu) !write(66,*) mxGetNumberOfDimensions(prhs(6))	  
	  !if(debu) !write(66,*) 'mxGetDimensions(prhs(6))'
	  !if(debu) !write(66,*) mxGetDimensions(prhs(6))	
	  call mxCopyPtrToPtrArray(mxGetDimensions(prhs(6)), dims ,mxGetNumberOfDimensions(prhs(6)))
	  !if(debu) !write(66,*) 'dims'
	  !if(debu) !write(66,*) dims
	  if(dims(1) .ne. 4 .or. dims(2) .ne. 3 .or. dims(3) .ne. 6) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &  
		 'Input 6 must be face 4x3x6')
      endif
      siz = dims(1)*dims(2)*dims(3)
      face_pr = mxGetPr(prhs(6))
	  allocate(face(4,3,6))
      call mxCopyPtrToReal8(face_pr, face, siz)
	  !if(debu) !write(66,*) 'input6 checked face = ', face


!     Check that input #7and fetch it
      m = mxGetM(prhs(7))
      n = mxGetN(prhs(7))
      if(m .ne. 6 .or. n .ne. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 7 must be r_f 6x3.')
      endif
      siz = m*n
      r_f_pr = mxGetPr(prhs(7))
	  allocate(r_f(6,3))
      call mxCopyPtrToReal8(r_f_pr, r_f, siz)
	  !if(debu) !write(66,*) 'input7 checked r_f = ', r_f

!     Check that input #8and fetch it
	  !if(debu) !write(66,*) 'mxGetNumberOfDimensions(prhs(8))'
	  !if(debu) !write(66,*) mxGetNumberOfDimensions(prhs(8))	  
	  !if(debu) !write(66,*) 'mxGetDimensions(prhs(8))'
	  !if(debu) !write(66,*) mxGetDimensions(prhs(8))	
	  call mxCopyPtrToPtrArray(mxGetDimensions(prhs(8)), dims ,mxGetNumberOfDimensions(prhs(8)))
	  !if(debu) !write(66,*) 'dims'
	  !if(debu) !write(66,*) dims
	  if(dims(1) .ne. 4 .or. dims(2) .ne. 3 .or. dims(3) .ne. 6) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &  
		 'Input 8 must be r_e 4x3x6')
      endif
      siz = dims(1)*dims(2)*dims(3)
      r_e_pr = mxGetPr(prhs(8))
	  allocate(r_e(4,3,6))
      call mxCopyPtrToReal8(r_e_pr, r_e, siz)
	  !if(debu) !write(66,*) 'input8 checked r_e = ', r_e


!     Check that input #9and fetch it
	  !if(debu) !write(66,*) 'mxGetNumberOfDimensions(prhs(9))'
	  !if(debu) !write(66,*) mxGetNumberOfDimensions(prhs(9))	  
	  !if(debu) !write(66,*) 'mxGetDimensions(prhs(9))'
	  !if(debu) !write(66,*) mxGetDimensions(prhs(9))	
	  call mxCopyPtrToPtrArray(mxGetDimensions(prhs(9)), dims ,mxGetNumberOfDimensions(prhs(9)))
	  !if(debu) !write(66,*) 'dims'
	  !if(debu) !write(66,*) dims
	  if(dims(1) .ne. 4 .or. dims(2) .ne. 3 .or. dims(3) .ne. 6) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &  
		 'Input 9 must be u_e 4x3x6')
      endif
      siz = dims(1)*dims(2)*dims(3)
      u_e_pr = mxGetPr(prhs(9))
	  allocate(u_e(4,3,6))
      call mxCopyPtrToReal8(u_e_pr, u_e, siz)
	  !if(debu) !write(66,*) 'input9 checked u_e = ', u_e	  
	  
	  
!     Check that input #10 and fetch it
	  !if(debu) !write(66,*) 'mxGetNumberOfDimensions(prhs(10))'
	  !if(debu) !write(66,*) mxGetNumberOfDimensions(prhs(10))	  
	  !if(debu) !write(66,*) 'mxGetDimensions(prhs(10))'
	  !if(debu) !write(66,*) mxGetDimensions(prhs(10))	
	  call mxCopyPtrToPtrArray(mxGetDimensions(prhs(10)), dims4 ,mxGetNumberOfDimensions(prhs(10)))
	  !if(debu) !write(66,*) 'dims4'
	  !if(debu) !write(66,*) dims4
	  if(dims4(1) .ne. 3 .or. dims4(2) .ne. 3 .or. dims4(3) .ne. 2 .or. dims4(4) .ne. 6) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &  
		 'Input 10 must be u_e 3x3x2x6')
      endif
      siz = dims4(1)*dims4(2)*dims4(3)*dims4(4)
      triang_face_pr = mxGetPr(prhs(10))
	  allocate(triang_face(3,3,2,6))
      call mxCopyPtrToReal8(triang_face_pr, triang_face, siz)
	  !if(debu) !write(66,*) 'input10 checked triang_face = ', triang_face	  	  
	  
	  
!     Check that input #11and fetch it
      m = mxGetM(prhs(11))
      n = mxGetN(prhs(11))
      if(m .ne. 6 .or. n .ne. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 11 must be idxf_dot 6x2.')
      endif
      siz = m*n
      idxf_dot_pr = mxGetPr(prhs(11))
	  allocate(idxf_dot_re(6,2))
	  allocate(idxf_dot(6,2))
      call mxCopyPtrToReal8(idxf_dot_pr, idxf_dot_re, siz)
	  idxf_dot=int(idxf_dot_re,8)
	  !if(debu) !write(66,*) 'input11 checked idxf_dot = ', idxf_dot	  
	  deallocate(idxf_dot_re)
	  
!     Check that input #12and fetch it
	  !if(debu) !write(66,*) 'mxGetNumberOfDimensions(prhs(12))'
	  !if(debu) !write(66,*) mxGetNumberOfDimensions(prhs(12))	  
	  !if(debu) !write(66,*) 'mxGetDimensions(prhs(12))'
	  !if(debu) !write(66,*) mxGetDimensions(prhs(12))	
	  call mxCopyPtrToPtrArray(mxGetDimensions(prhs(12)), dims ,mxGetNumberOfDimensions(prhs(12)))
	  !if(debu) !write(66,*) 'dims'
	  !if(debu) !write(66,*) dims
	  if(dims(1) .ne. 2 .or. dims(2) .ne. 3 .or. dims(3) .ne. 6) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &  
		 'Input 12 must be idxf_cross 2x3x6')
      endif
      siz = dims(1)*dims(2)*dims(3)
      idxf_cross_pr = mxGetPr(prhs(12))
	  allocate(idxf_cross(2,3,6))
	  allocate(idxf_cross_re(2,3,6))
      call mxCopyPtrToReal8(idxf_cross_pr, idxf_cross_re, siz)
	  idxf_cross = int(idxf_cross_re,8)
	  !if(debu) !write(66,*) 'input12 checked idxf_cross = ', idxf_cross	
	  deallocate(idxf_cross_re) 	  
	  
!     Check that input #13and fetch it
      m = mxGetM(prhs(13))
      n = mxGetN(prhs(13))
      if(m .ne. Np .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 13 must be w Npx1.')
      endif
      siz = m*n
      w_pr = mxGetPr(prhs(13))
	  allocate(w(Np))
      call mxCopyPtrToReal8(w_pr, w, siz)
	  !if(debu) !write(66,*) 'input13 checked w = ', w	 	  
	  
	  
!     Check that input #14and fetch it
      m = mxGetM(prhs(14))
      n = mxGetN(prhs(14))
      if(m .ne. 3 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 14 must be w 3x1.')
      endif
      siz = m*n
      d_pr = mxGetPr(prhs(14))
	  allocate(d(3))
      call mxCopyPtrToReal8(d_pr, d, siz)
	  !if(debu) !write(66,*) 'input14 checked d = ', d	  
	  
	  
!     Check that input #15and fetch it
      m = mxGetM(prhs(15))
      n = mxGetN(prhs(15))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 15 must be N_thread 1x1.')
      endif
      siz = m*n
      N_thread_pr = mxGetPr(prhs(15))
      call mxCopyPtrToReal8(N_thread_pr, N_thread_re, siz)
	  N_thread=int(N_thread_re,8)
	  !if(debu) !write(66,*) 'input14 checked N_thread = ', N_thread		  
	  
! call the computational subroutine.
      allocate(G(L*MM*NN))
	  !if(debu) !write(66,*) 'chiamo la subroutine'
      call computeGREEN_f90(L,MM,NN,xyz_src,Np,face,r_f,r_e,u_e,triang_face,idxf_dot,&
                            idxf_cross,w,d,N_thread,G)
      !if(debu) !write(66,*) 'assembled'
	  deallocate(xyz_src,face,r_f,r_e,u_e,triang_face,idxf_dot,idxf_cross,w,d)
! Create a matrix for the return argument
      mo=L*MM*NN
      no=1
	  ComplexFlag = 0
      plhs(1) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
! Load the output into a MATLAB array.
      G_pr = mxGetPr(plhs(1))
      siz=mo*no
      call mxCopyReal8ToPtr(G, G_pr, siz)
      !if(debu) !write(66,*) 'copied G'
	  
	  deallocate(G)
!	  !if(debu) close(66)
      return
end





subroutine computeGREEN_f90(L,M,N,xyz_src,Np,face,r_f,r_e,u_e,triang_face,idxf_dot,&
                            idxf_cross,w,d,N_thread,G)
integer*8 L,M,N,Np
integer*8 mx, my, mz ,ee 
real*8 r_m(3), d(3)
real*8 xyz_src(Np,3)
real*8 xyz_trg(Np,3)
real*8 face(4,3,6), triang_face(3,3,2,6)
real*8 r_f(3),r_e(4,3,6),u_e(4,3,6)
integer*8 idxf_dot(6,2),idxf_cross(2,3,6)
integer*8 N_thread
real*8 w(Np), y
real*8 G(L*M*N), mm(3)

!write(66,*) 'sono qui start'



call omp_set_num_threads(N_thread)
!$OMP PARALLEL SHARED(L,M,N,d,Np,xyz_src,face,r_f,r_e,u_e,triang_face,idxf_dot,idxf_cross,w,G)
!$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(mx,my,mz,mm,ee,r_m,xyz_trg,y)
do mx = 1,L
        do my = 1,M
            do mz = 1,N
			
			    !!write(66,*) 'sono qui 1'
                mm(1) = real(mx,8)
				mm(2) = real(my,8)
				mm(3) = real(mz,8)
				
				!!write(66,*) 'sono qui 2'
				
                !Centre of the observation voxel
				do ee = 1,3
					r_m(ee) = (mm(ee)-1.d0) * d(ee);
				enddo
				
				!!write(66,*) 'sono qui 3'
				
				!!write(66,*) 'r_m', r_m
				
				do ee = 1,Np
					xyz_trg(ee,1:3) = xyz_src(ee,1:3) + r_m !Gauss points on target voxel
				enddo
				
				!!write(66,*) 'sono qui 4'
				
				!!write(66,*) 'xyz_trg', xyz_trg
				
				call Integrate_NumAn3(face,r_f,r_e,u_e,triang_face,idxf_dot,&
				         idxf_cross,xyz_trg,w,Np,y)
						 
				!!write(66,*) 'sono qui 5'		 
						 
				!!write(66,*) 'y', y 		 
                G(mx+(my-1)*L+(mz-1)*L*M) =  y
            enddo
        enddo
enddo
!$OMP END DO NOWAIT 
!$OMP END PARALLEL 

!write(66,*) 'sono qui end'
 
end subroutine computeGREEN_f90




!! Integrale 1/R SUBROUTINE FORTRAN
!************************************************************************
!!Calcolo integrale int_V 1/R d^3r' con R = |r'-r|
!Ref: [1] "Magnetic Flux Density and Vector Potential of Uniform Polyhedral
!Sources" (M. Fabbri)

!INPUT: - reference voxel properties...: containing all source voxel properties (normals,
!vertex, ecc...)
!       - idxf_dot/idxf_cross: indici servono per cross/dot products
!       - P_target: gauss points on target voxel
!       - w: Gauss weigths 
!************************************************************************

subroutine Integrate_NumAn3(face,r_f,r_e,u_e,triang_face,idxf_dot,idxf_cross,P_target,w,dims1,y)
implicit none
real*8 face(4,3,6), triang_face(3,3,2,6)
real*8 r_f(6,3),r_e(4,3,6),u_e(4,3,6)
integer*8 ii,kk,idxf_dot(6,2),idxf_cross(2,3,6)
integer*8 dims1
real*8 P_target(dims1,3),w(dims1)
integer*8 Np, jj, ee
real*8 v1(3), v2(3), v3(3), v4(3)
real*8 v1_1(3),v1_2(3),v2_1(3),v2_2(3),v3_1(3),v3_2(3),v4_1(3),v4_2(3)
real*8 epsil1, epsil2, epsil3, epsil4
real*8 w_e(4)
real*8 omega1, omega2, omega
real*8 sub4(4,3),sub3(3),crn(4,3),dtn,dt,W_f
real*8  GG(6,dims1)
real*8 ftp1_1(3), ftp2_1(3), ftp3_1(3), ftp4_1(3)
real*8 ftp1_2(3), ftp2_2(3), ftp3_2(3), ftp4_2(3)
real*8 D1, D2, pi, I(dims1), y


pi=4.0d0*atan(1.0d0)

Np = dims1
GG(:,:) = 0.d0


do kk = 1,Np
    
    do ii = 1,6
      
        !!Compute we eq. (18)
        v1(1:3) = face(2,1:3,ii)-face(1,1:3,ii)
        v1_1(1:3) = face(2,1:3,ii) - P_target(kk,1:3)
        v1_2(1:3) = face(1,1:3,ii) - P_target(kk,1:3)
        epsil1 = norm2(v1(1:3))/(norm2(v1_1(1:3))+norm2(v1_2(1:3)))
        w_e(1) = log((1.0d0+epsil1)/(1.0d0-epsil1))
        
        v2(1:3) = face(3,1:3,ii)-face(2,1:3,ii)
        v2_1(1:3) = face(3,1:3,ii) - P_target(kk,1:3)
        v2_2(1:3) = face(2,1:3,ii) - P_target(kk,1:3)
        epsil2 = norm2(v2(1:3))/(norm2(v2_1(1:3))+norm2(v2_2(1:3)))
        w_e(2) = log((1.0d0+epsil2)/(1.0d0-epsil2))
        
        v3(1:3) = face(4,1:3,ii)-face(3,1:3,ii)
        v3_1(1:3) = face(4,1:3,ii) - P_target(kk,1:3)
        v3_2(1:3) = face(3,1:3,ii) - P_target(kk,1:3)
        epsil3 = norm2(v3(1:3))/(norm2(v3_1(1:3)) + norm2(v3_2(1:3)))
        w_e(3) = log((1.0d0+epsil3)/(1.0d0-epsil3))
        
        v4(1:3) = face(1,1:3,ii)-face(4,1:3,ii)
        v4_1(1:3) = face(1,1:3,ii) - P_target(kk,1:3)
        v4_2(1:3) = face(4,1:3,ii) - P_target(kk,1:3)
        epsil4 = norm2(v4(1:3))/(norm2(v4_1(1:3)) + norm2(v4_2(1:3)))
        w_e(4) = log((1.0d0+epsil4)/(1.0d0-epsil4))
        

        !!Compute D eq. (21)
        ftp1_1 = triang_face(1,1:3,1,ii)-P_target(kk,1:3)
        ftp2_1 = triang_face(2,1:3,1,ii)-P_target(kk,1:3)
        ftp3_1 = triang_face(3,1:3,1,ii)-P_target(kk,1:3)
        ftp4_1 = triang_face(4,1:3,1,ii)-P_target(kk,1:3)
        
        ftp1_2 = triang_face(1,1:3,2,ii)-P_target(kk,1:3)
        ftp2_2 = triang_face(2,1:3,2,ii)-P_target(kk,1:3)
        ftp3_2 = triang_face(3,1:3,2,ii)-P_target(kk,1:3)
        ftp4_2 = triang_face(4,1:3,2,ii)-P_target(kk,1:3)
        
        D1 =  norm2(ftp1_1(1:3))*norm2(ftp2_1(1:3))*norm2(ftp3_1(1:3)) &
               + norm2(ftp3_1(1:3))*dot_product(ftp1_1(1:3),ftp2_1(1:3)) &
               + norm2(ftp2_1(1:3))*dot_product(ftp1_1(1:3),ftp3_1(1:3)) &
               + norm2(ftp1_1(1:3))*dot_product(ftp2_1(1:3),ftp3_1(1:3))  
        D2 =  norm2(ftp1_2(1:3))*norm2(ftp2_2(1:3))*norm2(ftp3_2(1:3)) &
               + norm2(ftp3_2(1:3))*dot_product(ftp1_2(1:3),ftp2_2(1:3)) &
               + norm2(ftp2_2(1:3))*dot_product(ftp1_2(1:3),ftp3_2(1:3)) &
               + norm2(ftp1_2(1:3))*dot_product(ftp2_2(1:3),ftp3_2(1:3))

        call cross(ftp2_1(1:3),ftp3_1(1:3),v1(1:3))
        omega1 = 2.d0*atan2(dot_product(ftp1_1(1:3),v1(1:3)),D1)
        call cross(ftp2_2(1:3),ftp3_2(1:3),v2(1:3))
        omega2 = 2.d0*atan2(dot_product(ftp1_2(1:3),v2(1:3)),D2)
        omega = omega1 + omega2
        
        !!Compute W_f eq. (17)
        !Compute cross product with normals once for all
		do ee = 1,4
			sub4(ee,1:3) = r_e(ee,1:3,ii)  - P_target(kk,1:3)
		enddo
		
		! crn = sub4(1:4,idxf_cross(1,1:3,ii))*idxf_cross(2,1:3,ii);
        do ee = 1,3
		
			crn(1:4,ee)=sub4(1:4,idxf_cross(1,ee,ii))*real(idxf_cross(2,ee,ii),8)
        enddo
		
        sub3 = r_f(ii,1:3) - P_target(kk,1:3);
        dtn = sub3(idxf_dot(ii,1))*real(idxf_dot(ii,2),8); 
        
        W_f = 0.0d0
        do jj = 1,4
            dt = dot_product(crn(jj,1:3),u_e(jj,1:3,ii))
            W_f = W_f + dt*w_e(jj)
        enddo
        W_f = W_f - dtn*omega

        GG(ii,kk) = GG(ii,kk) + dtn*W_f
        
    enddo    
enddo


I(1:Np)=0.d0
do jj = 1,Np
do ii = 1,6
   I(jj) = I(jj) + 0.5d0*GG(ii,jj)
enddo
enddo

!! Numerical Integration
y = dot_product(w,I)

!! Add Constants
y = y/(4.0d0*pi)

end subroutine Integrate_NumAn3


    
subroutine cross(a, b, c)
implicit none
real*8, DIMENSION(3) :: c
real*8, DIMENSION(3), INTENT(IN) :: a, b
c(1) = a(2) * b(3) - a(3) * b(2)
c(2) = a(3) * b(1) - a(1) * b(3)
c(3) = a(1) * b(2) - a(2) * b(1)
end subroutine cross