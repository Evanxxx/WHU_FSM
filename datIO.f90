    module datIO

    use para
    use algorithm
    use SHTOOLS

    contains

    subroutine crustIn(filename, mat)
    !****************************************************************************
    !
    !  Crust1.0 Model Longtitude starts at inverse prime meridian which should
    !  transform into prime meridian coordinate
    !
    !****************************************************************************

    implicit none

    integer :: status, i, j, l, k, n, m
    character(*),intent(in):: filename
    real(kind=double),dimension(np,crsLa,crsLo), intent(out) :: mat

50  open (unit=50, file=trim(crsDir)//filename, status='old', action='read', iostat=status)
    fileopen: if(status == 0) then
        Latitude: do i=1,crsLa
            Longtitude: do j=1,crsLo
                read(50, *) (mat(k,i,j), k=1, np)
            end do Longtitude
        end do Latitude
        close(50)
        write(*,"(1x, a24, ' input completed!')") trim(filename)
        mat=mat*1000.d0 !Unit convert.
    else
        write(*, 100) trim(filename)
100     format(1x, a24, '-open failed')
    end if fileopen
    end subroutine crustIn

    subroutine realBin(filename, oMat)
    !****************************************************************************
    !  Read especially the ETOPO1 topograpy files.
    !  Etopo1:: cell-registered file
    !          ncols         21600
    !          nrows         10800
    !          xllcenter     -180.00
    !          yllcenter     -90.00
    !          cellsize      0.01666666667
    !          NODATA_value  -99999
    !          byteorder     LSBFIRST
    !          NUMBERTYPE    4_BYTE_FLOAT
    !
    !****************************************************************************
    implicit none

    character(*),intent(in):: filename
    real(kind=double),allocatable,dimension(:,:),intent(out):: oMat

    integer :: status, i, j, count
    character(120)::f
    real(kind=binDTM),allocatable,dimension(:,:):: tmpMat
    real(kind=double),allocatable,dimension(:,:):: scalMat
51  open (unit=51, file=filename, status='old',iostat=status,form='unformatted',access='stream')

    allocate(tmpMat(dtmLa,dtmLo),scalMat(dtmLa,dtmLo))

    fileOpen: if( status == 0) then
        Latitude: do i=1, dtmLa
            read(51) (tmpMat( i, j),j=1,dtmLo)
        end do Latitude
        close(51)
        write(*,"(1x, a24, ' input completed!')") trim(filename)
    else
        write(*, 100) filename
100     format(1x, a24, '-open failed')
    end if fileopen
    scalMat=dble(tmpMat)

    deallocate(tmpMat)
    !Scale imported dtm data into target resolution
    if(uniLa<dtmLa) then
        allocate(oMat(uniLa,uniLo))
        call scaling(scalMat,oMat,0.2d0) !In meter
    else
        allocate(oMat(dtmLa,dtmLo))
        oMat=scalMat
    end if

    end subroutine realBin

    subroutine gemmaIn(matTB,matRho,polynm,prov,ratio)
    !****************************************************************************
    !
    ! Read GEMMA files of upper/lower boundaries, density variation, density 
    !  polynomial funciton, geological provinces and global scale factor
    !  matTB:  upper/lower boundaries listed as Upper:(2i-1,:,:) and Lower: (2i,:,:)
    !  matRho: density variation, especially in crustal layer the averaged density
    !  polynm: density polynomial parameters associated with geological provinces
    !  ratio:  global scale factor 
    !
    !****************************************************************************
    integer,intent(out)::  prov(gmaLa,gmaLo)
    real(kind=double),intent(out):: polynm(12,3),matTB(14,gmaLa,gmaLo),matRho(7,gmaLa,gmaLo),&
        ratio(gmaLa,gmaLo)

    integer:: i,j,n,m,k,status
    real(kind=double):: tmp1,tmp2,thk,&
        continTab(101, 6), oceanTab(5,2),oceanTmp(5)

    do n=1,7
52      open(unit=52,file=trim(gmaDir)//trim(gmaLyr(n))//'_top.asc',status='old',action='read',iostat=status)
        io1if:if(status==0) then
            skip1do: do i=1,5
                read(52,*)
            end do skip1do
            read1do: do j=1,gmaLa
                read(52,*)  matTB(n*2-1,j,:)
            end do read1do
            write(*,"(a8,'_top input completed!')") trim(gmaLyr(n))
        else
            write(*,"(a8, '_top input failed!')") trim(gmaLyr(n))
        end if io1if
        close(52)

        open(unit=52,file=trim(gmaDir)//trim(gmaLyr(n))//'_bottom.asc',status='old',action='read',iostat=status)
        io2if:if(status==0) then
            skip2do: do i=1,5
                read(52,*)
            end do skip2do
            read2do: do j=1,gmaLa
                read(52,*)  matTB(n*2,j,:)
            end do read2do
            write(*,"(a8,'_bottom input completed!')") trim(gmaLyr(n))
        else
            write(*,"(a8, '_bottom input failed!')") trim(gmaLyr(n))
        end if io2if
        close(52)
        crustif: if(n/=6) then
            open(unit=52,file=trim(gmaDir)//trim(gmaLyr(n))//'_density.asc',status='old',action='read',iostat=status)
            io3if:if(status==0) then
                skip3do: do i=1,5
                    read(52,*)
                end do skip3do
                read3do: do j=1,gmaLa
                    read(52,*)  matRho(n,j,:)
                end do read3do
                write(*,"(a8,'_density input completed!')") trim(gmaLyr(n))
            else
                write(*,"(a8, '_density input failed!')") trim(gmaLyr(n))
            end if io3if
            close(52)
        else
            open(unit=52,file=trim(gmaDir)//'crust_Geo_prov.asc',status='old',action='read',iostat=status)
            if(status==0) then
                do i=1,5
                    read(52,*)
                end do
                do j=1,gmaLa
                    read(52,*)  prov(j,:)
                end do
                write(*,"('crust_Geo_prov input completed!')")
            else
                write(*,"('crust_Geo_prov input failed!')")
            end if
            close(52)

            open(unit=52,file=trim(gmaDir)//'crust_cal_par.asc',status='old',action='read',iostat=status)
            if(status==0) then
                do i=1,5
                    read(52,*)
                end do
                do j=1,gmaLa
                    read(52,*)  ratio(j,:)
                end do
                write(*,"('crust_cal_par.asc input completed!')")
            else
                write(*,"('crust_cal_par.asc input failed!')")
            end if
            close(52)

            open(unit=52,file=trim(gmaDir)//'gemma_density_v02.txt',status='old',action='read',iostat=status)
            if(status==0) then
                do i=1,32 !Head Skipping
                    read(52,*)
                end do
                do j=1,101
                    read(52,*)  tmp1,tmp2,continTab(j,:)
                end do
                do i=1,6 !Gap Skipping
                    read(52,*)
                end do
                do j=1,5
                    read(52,*) tmp1,tmp2,oceanTab(j,:)
                end do
                write(*,"('gemma_density_v02.txt input completed!')")
            else
                write(*,"('gemma_density_v02.txt input failed!')")
            end if
            close(52)

            open(unit=52,file=trim(gmaDir)//'crust_polynomials.txt',status='old',action='read',iostat=status)

            if(status==0) then
                do i=1,13 !Head Skipping
                    read(52,*)
                end do
                do j=1,12
                    read(52,*)  polynm(j,:)
                end do
                write(*,"('crust_polynomials.txt input completed!')")
            else
                write(*,"('crust_polynomials.txt input failed!')")
            end if
            close(52)

            do i=1,gmaLa
                do j=1,gmaLo
                    thk=matTB(11,i,j)-matTB(12,i,j)
                    !T=matTB(5,:,:) B=matTB(10,:,:) M=matTB(12,:,:)
                    if(prov(i,j)==1) then
                        if(thk<polynm(1,1)) then
                            matRho(6,i,j)=polynm(1,2)+polynm(1,3)*thk/2.d0
                        else if(thk<polynm(2,1)) then
                            matRho(6,i,j)=polynm(2,2)+(polynm(1,2)-polynm(2,2))*polynm(1,1)/thk+&
                                          (polynm(1,3)-polynm(2,3))*polynm(1,1)**2/thk/2.d0+&
                                          polynm(2,3)*thk/2.d0
                        else 
                            matRho(6,i,j)=polynm(3,2)+(polynm(1,2)-polynm(2,2))*polynm(1,1)/thk+&
                                          (polynm(2,2)-polynm(3,2))*polynm(2,1)/thk+&
                                          (polynm(1,3)-polynm(2,3))*polynm(1,1)**2/thk/2.d0+&
                                          (polynm(2,3)-polynm(3,3))*polynm(2,1)**2/thk/2.d0+&
                                          polynm(3,3)*thk/2.d0
                        end if
                    else if(prov(i,j)==8) then
                        matRho(6,i,j)=2600.d0
                    else
                        matRho(6,i,j)=(sum(continTab(1:floor(thk),prov(i,j)-1))+&
                            (thk-floor(thk))*continTab(floor(thk)+1,prov(i,j)-1))/thk
                    end if
                end do
            end do
            matRho(6,:,:)=matRho(6,:,:)*ratio
        end if crustif
    end do
    polynm(:,1)=polynm(:,1)*1.d3 ! km to m
    polynm(:,3)=polynm(:,3)*1.d-3 ! km to m
    matTB=matTB*1000.d0 !unit convert
    end subroutine

    subroutine distGravity(filename,lmax,distGrav)
    !****************************************************************************
    !
    !  Read in gravity disturbance files calculated at ICGEM
    !
    !****************************************************************************

    character(*),intent(in):: filename
    integer,intent(in):: lmax
    real(kind=double),allocatable,dimension(:,:),intent(out):: distGrav

    integer:: lmax0,samples,status,i,j,k
    real(kind=double):: tmp1,tmp2,tmp3

    allocate(distGrav(361,721))

53  open(unit=53,file=trim(ggmDir)//filename//".gdf",status='old',action='read',iostat=status)
    fileopen: if(status==0) then
        do k=1,35
            read(53,*)
        end do
        do i=1,361
            do j=1,721
                read(53,*) tmp1,tmp2,tmp3,distGrav(i,j)
            end do
        end do
        close(53)
    end if fileopen
    end subroutine distGravity
    
    subroutine xyzOut( filename,iMat2, iMat3)
    !****************************************************************************
    !
    ! Output .xyz files convenient for GMT 
    !  Flag=.TRUE. for 3-Dimension input array
    !  Flag=.FALSE. for 2-Dimension input array
    !
    !****************************************************************************

    implicit none

    integer:: status, i, j, n ,nCol,nRow,nLayer
    real:: div

    character(*),intent(in):: filename
    real(kind=double),dimension(:, :),intent(in),optional::iMat2
    real(kind=double),dimension(:, :, :),intent(in),optional::iMat3

    Dimensional:if(present(iMat3)) then
        nRow=size(iMat3,1)
        nCol=size(iMat3,2)
        nLayer=size(iMat3,3)
        div=360.d0/nCol

        tLayer: do n=1,nLayer

54          open(unit=54, file=trim(oDir)//trim(oriLyr(n))//filename, status='replace', action='write', iostat=status)

            multiFileopen: if(status==0) then
                rWrite: do i=1, nRow
                    write(54,101) (-180+j*div-div/2, 90-i*div+div/2,iMat3(i,j,n), j=1,nCol)
101                 format(1x,2f10.4,g24.16)
                end do rwrite
                close(54)
                write(*,102) trim(oriLyr(n))//trim(filename)
102             format(1x,a24,'-output completed!')
            else
                write(*,103) trim(oriLyr(n))//trim(filename)
103             format(1x,a24,'-replace falied')
            end if multiFileopen
        end do tLayer

    else if( present(iMat2)) then
        nRow=size(iMat2,1)
        nCol=size(iMat2,2)
        div=360.d0/nCol

        open(unit=54, file=trim(oDir)//filename, status='replace', action='write', iostat=status)

        fileopen: if(status==0) then
            do i=1,nRow
                write(54,101) (-180+j*div-div/2, 90-i*div+div/2,iMat2(i,j), j=1,nCol)
            end do
            close(54)
            write(*,102) filename
        else
            write(*,103) filename
        end if fileopen
    else
        write(*,"('Input Array dismatches!')")
    end if Dimensional

    end subroutine xyzOut

    subroutine xyzIn( filename,oMat,flag,rescale)

    !****************************************************************************
    !
    !  Read in .xyz files
    !   If flag=0 .xyz file formats in (latitude,longitude,z) order
    !     flag/=0 .xyz file formats in (longitude,latitude,z) order
    !   and
    !   If rescale=0 rescale output arrays into (uniLa,uniLo)
    !     rescale/=0 do not rescale
    !
    !****************************************************************************

    implicit none

    integer,optional,intent(in):: flag,rescale
    character(*),intent(in):: filename
    real(kind=double),allocatable,dimension(:,:),intent(out):: oMat

    integer:: i,j,n,status,interval
    real(kind=double):: tmp1,tmp2
    real(kind=double),dimension(2):: interArry
    real(kind=double),allocatable,dimension(:,:):: tmpMat

55  open(unit=55,file=filename,iostat=status,action='read')

    fileopen:if(status==0) then
        flagif:if(present(flag) .and. flag==0) then
            do i=1,2
                read(55,*) tmp1,interArry(i)
            end do
        elseif (present(flag) .and. flag/=0) then
            do i=1,2
                read(55,*) interArry(i)
            end do
        end if flagif
        interval=nint(1.d0/(interArry(2)-interArry(1)))
        write(*,"(a24,' interval=  ',i12)") filename,interval

        allocate(oMat(interval*180,interval*360),tmpMat(uniLa,uniLo))

        rewind(unit=55,iostat=status)
        filerewind:if(status==0) then
            do i=1,interval*180
                read(55,*) (tmp1,tmp2,oMat(i,j),j=1,interval*360)
            end do
            write(*,"(a24, ' input completed!')") filename
            rescalif: if(present(rescale) .and. rescale==0) then
                call scaling(oMat,tmpMat,factor=0.2d0)
                write(*,"(a24, ' rescaled!')") filename
                oMat=tmpMat
            else
                write(*,"(a24, ' didnt rescale!')") filename                
            end if rescalif
        else
            write(*,"(a24, ' rewind failed')") filename
        end if filerewind
    else
        write(*,"(a24, ' input failed!')") filename
    end if fileopen
    close(55)

    end subroutine xyzIn

    subroutine coeffiOut( filename, iMat, flag)
    !****************************************************************************
    !
    ! Output sphercial coefficicents or degree variances
    !  Flag=0 Normal coefficients
    !  Flag=1 Dimensionless degree variances: Degree Only
    !  Flag=2 2-Dimension degree variances: Degree and Order
    !
    !****************************************************************************

    implicit none

    integer:: status,i,j,power

    integer,optional,intent(in):: flag
    character(*),intent(in):: filename
    real(kind=double),dimension(:,:,:),intent(in):: iMat

    real(kind=double),allocatable,dimension(:):: dlv
    real(kind=double),allocatable,dimension(:,:):: d2v
    power=size(iMat,2)

    flagif: if(present(flag) .and. flag==1) then
        allocate(dlv(power))
        forall(i=1:power) dlv(i)=sum(iMat(:,i,:)**2.d0)
56      open(unit=56,file=trim(oDir)//'raw/'//filename//'.dlv',status='replace',action='write',iostat=status)
        dlvfileopen:if(status==0) then
            write(56,105) (i-1,dlv(i),i=1,power)
105         format(1x,i6,g24.16)
            write(*,"('Dimensionless degree variances output completed!')")
        else
            write(*,"('Dimensionless degree variances output failed!')")
        end if dlvfileopen
        close(56)
        deallocate(dlv)
    else if(present(flag) .and. flag==2) then
        allocate(d2v(power,power))
        forall(i=1:power,j=1:power) d2v(i,j)=iMat(1,i,j)**2.d0+iMat(2,i,j)**2.d0
        open(unit=56,file=trim(oDir)//'raw/'//filename//'.d2v',status='replace',action='write',iostat=status)
        d2vfileopen:if(status==0) then
            do i=1,power
                write(56,106) (i-1,j-1,d2v(i,j),j=1,i)
            end do
106         format(1x,2i6,g24.16)
            write(*,"('2-Dimension degree variances output completed!')")
        else
            write(*,"('2-Dimension degree variances output failed!')")
        end if d2vfileopen
        close(56)
        deallocate(d2v)
    else
        open(unit=56,file=trim(oDir)//'raw/'//filename//'.effi',status='replace',action='write',iostat=status)
        fileopen:if(status==0) then
            do i=1,power
                write(56,107) (i-1,j-i,iMat(2,i,i-j+1),j=1,i-1)
                write(56,107) (i-1,j-1,iMat(1,i,j),j=1,i)
            end do
107         format(1x,2i6,g24.16)
            write(*,"('Coeffients output completed!')")
        else
            write(*,"('Coeffients output failed!')")
        end if fileopen
        close(56)
    end if flagif
    write(*,"('Filename: ',a12,3x,'Degree: ', i6,3x)") filename, power
    end subroutine coeffiOut

    subroutine admitCorrOut(filename, admit, corr)
    !****************************************************************************
    !
    !  Output admittance and Coefficients Correaltion calculated from SHAdimtCorr
    !
    !****************************************************************************
    character(*),intent(in):: filename
    real(kind=double),dimension(:),intent(in):: admit, corr

    integer:: i,lmax

    lmax=size(admit)
57  open(unit=57,file=trim(oDir)//'raw/'//filename//'.corr',status='replace',action='write')
    write(57,108) (i-1,corr(i),i=1,lmax)
108 format(1x,i6,g12.4)
    close(57)
    write(*,"('Correlation output completed!')")
    end subroutine admitCorrOut
    
    subroutine scatterOut(filename,iMat)
    !****************************************************************************
    !
    ! Output the scatter point coodinates with gravity disturbance and depth files
    !
    !****************************************************************************
     implicit none 
     
     character(*),intent(in):: filename
     real(kind=double),dimension(:,:),intent(in):: iMat
     
     integer:: i,n
     
     n=size(iMat,2)
     
58   open(unit=58,file=trim(oDir)//'raw/'//filename//'.stm',status='replace',action='write')
     if(iMat(2,i)/=0) then
         write(58,109) (iMat(1,i),iMat(2,i)/1.d+3,i=1,n)
109      format(1x,2g12.4)
     end if
     close(58)
     write(*,"('Scatter map output comleted!')")
     end subroutine scatterOut
    
    end module datIO
