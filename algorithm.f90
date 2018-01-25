    module algorithm

    contains

    subroutine scaling(iMat,oMat,factor)

    !****************************************************************************
    !
    !  Scale input Matrix iMat(n×m) into target resolution oMat(an×am) under three conditions:
    !  1) Expansion. If oMat is certain times larger than iMat in both row and column (a>1),
    !                 the original data is completely copied into vacancies.     
    !  2) Sampling.  If iMat is certain times smaller than oMat in both row and column
    !                 (a<1), the original data is averaged according to the given
    !                 Gaussian Filter factor (0.d0 ~ 1.d0).
    !  3) Identical. If iMat and oMat are of the same size (a=1), the oMat==iMat
    !  #Funtion Call: algorithm::gassianFilter, matAppd, matPro, matTrans      
    !
    !****************************************************************************

    implicit none

    real(kind=double),optional,intent(in):: factor
    real(kind=double),intent(in),dimension(:,:):: iMat
    real(kind=double),intent(out),dimension(:,:):: oMat

    real(kind=double),allocatable,dimension(:):: filter1
    real(kind=double),allocatable,dimension(:,:):: filter2,tmpAvg,tmpTrans,tmpAppd,tmpPro

    integer:: iBnd1,iBnd2,oBnd1,oBnd2,i,j,n,div

    iBnd1=ubound(iMat,1)
    iBnd2=ubound(iMat,2)
    oBnd1=ubound(oMat,1)
    oBnd2=ubound(oMat,2)

    !*********                        Expansion                       *********
    scaleIf: if(oBnd1/iBnd1>1) then
        div=oBnd1/iBnd1
        forall(i=1:iBnd1,j=1:iBnd2) oMat((i-1)*div+1:i*div,(j-1)*div+1:j*div)=iMat(i,j)
    !*********                        Sampling                        *********
    else if (oBnd1/iBnd1<1) then
        if(present(factor)) then
            div=iBnd1/oBnd1
            allocate(filter1(1:2*div+1),filter2(1,1:2*div+1),tmpAvg(1:oBnd1,1:iBnd2))
            call gaussianFilter(2*div+1,factor,filter1D=filter1)

            allocate(tmpPro(1,iBnd2))
            filter2(1,:)=filter1
            latAvg: do i=1,oBnd1
                latIf: if(i==1) then
                    call matPro(filter2(:,1:div+1),iMat((i-1)*div+1:i*div+1,:),tmpPro,reset=0)
                    tmpAvg(i,:)=tmpPro(1,:)*2.d0
                else if(i==oBnd1) then
                    call matPro(filter2(:,1:div+1),iMat((i-1)*div:i*div,:),tmpPro,reset=0)
                    tmpAvg(i,:)=tmpPro(1,:)*2.d0
                else
                    call matPro(filter2,iMat((i-1)*div:(i+1)*div,:),tmpPro)
                    tmpAvg(i,:)=tmpPro(1,:)
                end if latIf
            end do latAvg
            deallocate(tmpPro)

            allocate(tmpAppd(oBnd1,2*div+1),tmpTrans(2*div+1,1),tmpPro(oBnd1,1))
            call matTrans(filter2,tmpTrans,reset=0)
            deallocate(filter2)
            call matAppd(tmpAvg(:,1:div),tmpAvg(:,iBnd2-div:iBnd2),tmpAppd,reset=0)
            longAvg: do j=1,oBnd2
                longIf:    if(j==1 .or. j==oBnd2)then
                    call matPro(tmpAppd,tmpTrans,tmpPro)
                    oMat(:,j)=tmpPro(:,1)
                else
                    call matPro(tmpAvg(:,(j-1)*div:(j+1)*div),tmpTrans,tmpPro)
                    oMat(:,j)=tmpPro(:,1)
                end if longIf
            end do longAvg
            deallocate(tmpTrans,tmpPro,tmpAppd)
        else
            write(*,"('Standard Deviation Parameter Missing!!!')")
        end if
    else
    !*********                        Identical                         *********
        oMat=iMat
    end if scaleIf

    end subroutine scaling

    subroutine gaussianFilter(size,sigma,filter1D,filter2D)
    !****************************************************************************
    ! 
    !  Generation of Guassian Filter for convolution in either 1D or 2D with 
    !   standard deviation sigma
    !
    !****************************************************************************

    implicit none

    integer,intent(in):: size
    real(kind=double),intent(in):: sigma
    real(kind=double),allocatable,dimension(:),optional,intent(out):: filter1D
    real(kind=double),allocatable,dimension(:,:),optional,intent(out):: filter2D

    integer:: i,j,r,even
    real(kind=double):: filtersum

    r=size/2
    even=mod(size,2)
    even=1-even

    if(present(filter1D)) then
        allocate(filter1D(1:size))
        forall(i=1:size) filter1D(i)=exp(-1*(int(i+0.5*even-1-r)**2.d0)/(2.d0*sigma**2.d0))
        filtersum=sum(filter1D)
        filter1D=filter1D/filtersum
    else if(present(filter2D)) then
        allocate(filter2D(1:size,1:size))
        forall(i=1:size,j=1:size) filter2D(i,j)=exp(-1*(int(i+0.5*even-1-r)**2.d0+ &
        	int(i+0.5*even-1-r)**2.d0)/(2.d0*sigma**2.d0))
        filtersum=sum(filter2D)
        filter2D=filter2D/filtersum
    else
        write(*,"('Filter dimension not specified!')")
    end if

    end subroutine gaussianFilter
    
    subroutine matAppd(iMat1,iMat2,oMat,flag,reset)
    !****************************************************************************
    !
    !  Append two 2-dimensional matrixs into one.
    !   if flag==0 append in column
    !      flag/=0 append in row 
    !  reset: recorder
    !
    !****************************************************************************

    implicit none

    integer,optional,intent(in):: flag,reset
    real(kind=double),dimension(:,:),intent(in):: iMat1,iMat2
    real(kind=double),allocatable,dimension(:,:),intent(out):: oMat

    integer:: i,j,i1s1,i1s2,i2s1,i2s2
    integer,save:: count

    i1s1=size(iMat1,1) ;    i1s2=size(iMat1,2)
    i2s1=size(iMat2,1) ;    i2s2=size(iMat2,2)

    resetIf: if(present(reset)) then
        count=0
    else
        count=count+1
    end if resetIf

    flagIf: if(present(flag) .and. flag/=0) then

        columnIf: if(i1s2==i2s2) then
            allocate(oMat(i1s1+i2s1,i1s2))
            oMat(1:i1s1,:)=iMat1
            oMat(i1s1+1:i1s1+i2s1,:)=iMat2
        else
            write(*,"(1x,i8,'   matAppd: Matirx size or squence dismatches in column !')") count
        end if columnIf
    else
        rowIf: if(i1s1==i2s1) then
            allocate(oMat(i1s1,i1s2+i2s2))
            oMat(:,1:i1s2)=iMat1
            oMat(:,i1s2+1:i1s2+i2s2)=iMat2
        else
            write(*,"(1x,i8,'matAppd: Matirx size or squence dismatches in row !')") count
        end if rowIf
    end if flagIf

    end subroutine matAppd

    subroutine matTrans(iMat,oMat,reset)
    !****************************************************************************
    !
    !  Transpose iMat into oMat. 
    !   reset: recorder.
    !
    !****************************************************************************
    implicit none

    integer,optional,intent(in):: reset
    real(kind=double),dimension(:,:),intent(in):: iMat
    real(kind=double),dimension(:,:),intent(out):: oMat

    integer:: i,j,is1,is2,os1,os2
    integer,save:: count

    is1=size(iMat,1) ;    is2=size(iMat,2)
    os1=size(oMat,1) ;    os2=size(oMat,2)

    resetIf: if(present(reset)) then
        count=0
    else
        count=count+1
    end if resetIf

    if(is1==os2 .and. is2==os1) then
        forall(i=1:is1,j=1:is2) oMat(j,i)=iMat(i,j)
    else
        write(*,"(1x,i8,'   matTrans: Matirx size or squence dismatches!')") count
    end if

    end subroutine matTrans

    subroutine matPro(iMat1,iMat2,oMat,reset)
    !****************************************************************************
    !
    !  Perform the matrix product, and the matrix size must match.
    !   reset: recorder.
    !
    !****************************************************************************
    implicit none

    integer,optional,intent(in):: reset
    real(kind=double),dimension(:,:),intent(in):: iMat1,iMat2
    real(kind=double),allocatable,dimension(:,:),intent(out):: oMat

    integer:: i,j,i1s1,i1s2,i2s1,i2s2
    integer,save:: count

    i1s1=size(iMat1,1) ;    i2s1=size(iMat2,1)
    i1s2=size(iMat1,2) ;    i2s2=size(iMat2,2)

    resetIf: if(present(reset)) then
        count=0
    else
        count=count+1
    end if resetIf

    allocate(oMat(i1s1,i2s2))

    if(i1s2==i2s1)then
        forall(i=1:i1s1,j=1:i2s2)
            oMat(i,j)=sum(iMat1(i,:)*iMat2(:,j))
        end forall
    else
        write(*,"(1x,i8,'   matPro: Matirx size or squence dismatches!')") count
    end if

    end subroutine matPro

    subroutine binarization(iMat,oMat)
    !****************************************************************************
    !
    !  Binarize the input Matrix into 1.d0 and 0.d0. If the element>0.d0, the
    !   output=1.d0; if the elememt<=0.d0, the output=0.d0 
    !
    !****************************************************************************

    implicit none

    real(kind=double),dimension(:,:),intent(in):: iMat
    real(kind=double),allocatable,dimension(:,:),intent(out):: oMat

    integer:: i,j,lat,lon
    lat=size(iMat,1)
    lon=size(iMat,2)
    allocate(oMat(lat,lon))
    
    do i=1,lat
        do j=1,lon
            if(iMat(i,j)>0.d0) then
                oMat(i,j)=1.d0
            else
                oMat(i,j)=0.d0
            end if
        end do
    end do
    end subroutine
    
    real(kind=double) function nchoosek(k,n)
    !****************************************************************************
    !
    !  Calculate the binomial coefficients of n choose k. 
    !
    !****************************************************************************
    implicit none

    integer,intent(in)::k,n

    integer::i

    nchoosek=1.d0
    if(k>n) then
        write(*,"('Fuction nchoosek error with k>n!')")
    else if(k==0 .or. k==n) then
        nchoosek=1.d0
    else
        Judge:if(k>n/2) then
            do  i=1,n-k
                nchoosek=nchoosek*(i+k)/i
            end do
        else
            do  i=1,k
                nchoosek=nchoosek*(i-k+n)/i
            end do
        end if Judge
    end if

    end function nchoosek
    real(kind=double) function pochInt(x,n)
    !****************************************************************************
    !
    !  Perform the Pochhammer function.
    !
    !****************************************************************************
    real(kind=double),intent(in)::x
    integer,intent(in)::n

    integer:: i
    real(kind=double):: tmp

    tmp=1.d0
    if(n<1) then
        write(*,"('Pochhammer function parameter n isnt a positive integer!')")
    else
        do i=0,n-1
            tmp=tmp*(x+i)
        end do
    end if
    pochInt=tmp

    end function pochInt

    end module algorithm