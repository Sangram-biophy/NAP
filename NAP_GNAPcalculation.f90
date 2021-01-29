program Non_affine_parameter_framewise
implicit none

real(kind=8),allocatable :: Ref(:,:) , inst(:,:) , dist(:), Ref_boxl(:), Inst_boxl(:)
real(kind=8),allocatable :: dr(:),X(:,:),dist1(:),dist2(:),dist3(:),epsil(:,:)
real(kind=8),allocatable :: Y(:,:),Y_inv(:,:),WORK(:),NAPeps1(:,:),NAPeps2(:,:)
real(kind=8),allocatable :: NAPeps3(:,:),NAPeps(:,:),NAP1(:,:)

real(kind=8) :: Romega,NAP,GNAP,NAP_k

integer :: i, j, k, iatom_protein, natom_protein, nframe, iframe, count1
integer :: l,m,n,o,p,q,r,s,t,z,a,b,c,d,e,f,g,h
integer :: IPIV(3), INFO,mu,alp,neu
integer,allocatable :: list_count(:,:),countf(:)
!character(len=4),allocatable :: protein_atom(:)
character(len=4) :: junk2
character(len=6) :: junk1

open(unit=2,file='try.gro', action='read')
open(unit=3,file='ref.gro', action='read')

! System Information
! Number of Protein Atoms (Lysozyme)
natom_protein = 8108
write(408,*)'Total Protein Atom',natom_protein
! Total Number of Snapshot
nframe= 201
write(408,*)'Total Frame',nframe
! Constants
Romega= 0.45d0

!open(unit=2,file='protein_traj.gro', action='read')
!open(unit=3,file='ref.gro', action='read')


allocate(Ref(natom_protein,3) , inst(natom_protein,3) , dist(natom_protein))
allocate(Ref_boxl(3), Inst_boxl(3),dr(3), list_count(natom_protein,natom_protein), countf(natom_protein))
allocate(dist1(natom_protein), dist2(natom_protein), dist3(natom_protein),X(3,3),epsil(3,3))
allocate(Y(3,3),Y_inv(3,3),WORK(3),NAPeps1(natom_protein,3),NAPeps2(natom_protein,3),NAPeps3(natom_protein,3))
allocate(NAPeps(natom_protein,3),NAP1(natom_protein,3))
!reading reference structure file
do iatom_protein = 1 , natom_protein
         read(3,*) junk1, junk2, i, Ref(iatom_protein,:)
         !write(20,*) junk1, junk2, i, Ref_protein(iatom_protein,:)
end do
read(3,*) Ref_boxl(:)
write(21,*) "reference box length", Ref_boxl(:) !reading box length of reference str

do i=1,natom_protein
        count1=0
     !!!calculation of distance, and counting and listing all atoms inside Romega radius
        do j=1,natom_protein
             dr(:)=Ref(j,:)-Ref(i,:)      
             dr(:)= dr(:)-Ref_boxl(:)*anint(dr(:)/Ref_boxl(:)) !accounting pbc 
             dist(j) = sqrt(dot_product(dr,dr))
            ! write(23,*) i, j, dist(j)
             if (dist(j)/=0.0 .and. dist(j)<Romega) then
                    count1=count1+1
                    list_count(i,count1)=j
             endif
        enddo
        countf(i)=count1
        !write(24,*) i, count1
!        do k=1,count1
 !            write(25,*) i, count1, list_count(i,k)
  !      enddo
enddo
!do i=1,natom_protein
 !       l=countf(i)
  !      do j=1,l
   !             write(35,*) i, j, list_count(i,j)
do iframe=1,nframe
        read(2,*)
        read(2,*)
        do iatom_protein = 1 , natom_protein
                read(2,*) junk1, junk2, z, inst(iatom_protein,:)     !reading simulation file
                !write(22,*) junk1, junk2, i, Inst_protein(iatom_protein,:)
        end do
        read(2,*) Inst_boxl(:)
        NAP_k=0.0d0    !initiating value of sum of NAP in a frame to zero
        do k=1,natom_protein
                
                l=countf(k)
                do s=1,3
                        do t=1,3
                                X(s,t)=0.0d0 !assigning zero values to all elements of X matrix(order:3X3)
                        enddo
                enddo
                !caculating X matrix
                do i=1,3
                        do j=1,3
                                do m=1,l
                                        n=list_count(k,m)
                                        X(i,j)=X(i,j)+del1(inst(n,i),inst(k,i),Inst_boxl(i))*del1(Ref(n,j),Ref(k,j),Ref_boxl(j))
                                       ! write(55,*) k,i,j,X(i,j)
                                 enddo
!                                 write(56,*) k,i,j,X(i,j)
                        enddo
                enddo                        
                !assigning a 3X3 zero matrix to Y 
                do a=1,3
                        do b=1,3
                                Y(a,b)=0.0d0
                        enddo
                enddo
                !calculating Y matrix
                do c=1,3
                        do d=1,3
                                do e=1,l
                                        n=list_count(k,e)
                                        Y(c,d)=Y(c,d)+del1(Ref(n,c),Ref(k,c),Ref_boxl(c))*del1(Ref(n,d),Ref(k,d),Ref_boxl(d))
                                       ! write(55,*) k,i,j,X(i,j)
                                        
                                 enddo
                                 Y_inv(c,d)=0.0d0
                                 epsil(c,d)=0.0d0
                !                 write(57,*) k,c,d,Y(c,d)
                        enddo
                enddo
                Y_inv=Y
                !!!!!!!!!===============================================================
                !!!!!!!!inverting Y matrix assigned to Y_inv using LAPACK module
               ! write(100,*) 'Y=', Y_inv
                call dgetrf(size(Y_inv, 1), size(Y_inv, 2), Y_inv, size(Y_inv, 1), IPIV, INFO)
                call dgetri(size(Y_inv, 1), Y_inv, size(Y_inv, 1), IPIV, WORK, size(WORK), INFO);
                !write(100,*) 'inverse of Y=', Y_inv     
                !!!!!!!==================================================================
                !!!!calculating epsilon 3X3 matrix
                do mu=1,3
                        do alp=1,3
                                do neu=1,3
                                        epsil(mu,alp)=epsil(mu,alp)+X(mu,neu)*Y_inv(neu,alp)
                                enddo
                   !             write(101,*) mu,alp,epsil(mu,alp)
                        enddo
                enddo
                !obtaining two parts of the Falk-Langer NAP equation
                do f=1,l
                        g=list_count(k,f)
                        do h=1,3
                                NAPeps1(g,h)=epsil(h,1)*del1(Ref(g,1),Ref(k,1),Ref_boxl(1));
                                NAPeps2(g,h)=epsil(h,2)*del1(Ref(g,2),Ref(k,2),Ref_boxl(2));
                                NAPeps3(g,h)=epsil(h,3)*del1(Ref(g,3),Ref(k,3),Ref_boxl(3));
                                NAPeps(g,h)=NAPeps1(g,h)+NAPeps2(g,h)+NAPeps3(g,h)
                 !               write(102,*) f,g,h,NAPeps(g,h)
                                NAP1(g,h)=del1(inst(g,h),inst(k,h),Ref_boxl(h))
                  !              write(103,*) f,g,h,NAP1(g,h)
                        enddo
               enddo
               NAP=0.0d0
               !calculating NAP of each of the atom
               do p=1,l
                        q=list_count(k,p)
                        NAP=NAP+((NAP1(q,1)-NAPeps(q,1))**2)+((NAP1(q,2)-NAPeps(q,2))**2)+((NAP1(q,3)-NAPeps(q,3))**2)
                        
               enddo
               write(105,*) iframe,k,NAP
               NAP_k=NAP_k+NAP    !sum of all NAP in a frame
        enddo
        !!!!!GNAP calculation(average of all NAP in a frame)
        GNAP=real(NAP_k/natom_protein)
        write(106,*) iframe,GNAP
enddo
!defining a function del1 for subtracting coordinates along with accounting pbc
contains
real function del1(coor1,coor2,boxlength)
        implicit none
        real(kind=8),intent(in)::coor1,coor2,boxlength
        del1=(coor1-coor2)-(boxlength*anint((coor1-coor2)/boxlength))
end function del1
end program Non_affine_parameter_framewise
