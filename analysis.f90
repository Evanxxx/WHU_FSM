    module analysis

    use para
    use algorithm
    use datIO
    use SHTOOLS

    contains

    !Truncated error omited!

    subroutine kernal(nchoosek,rhoRatio,low,up,lmax,power,effi,flag)
    !****************************************************************************
    ! Calculate the kernal (Height Function) in FSM.
    !  
    !  rhoRatio: spatial variation of density devided by the averaged Earth 
    !            density
    !  up,low:   upper/lower boundary function
    !  lmax:     the maxium of calculated order/degree
    !  power:    the truncated degree of height binomial expansion
    !  effi:     the resultant spherical coefficients
    !  flag:     demonstating the very FSM is to perform, normal FSM for lateral
    !            FSM computation or higher order FSM engaged in radial FSM 
    !            computation. If flag==0 normal FSM, flag/=0 higher order FSM
    ! #Function call: algorithm:: binarization, nchoosek, scaling
    !                 para::      radius,
    !                 SHTOOLS::   SHExpandDH
    !                  
    !****************************************************************************
    integer,intent(in):: lmax,power,flag
    real(kind=double),external:: nchoosek
    real(kind=double),dimension(:,:),intent(in)::rhoRatio,low,up
    real(kind=double),allocatable,dimension(:,:,:),intent(out):: effi

    integer:: i,j,nlat,nlon,index,lmaxDH
    real(kind=double)::ref
    real(kind=double),allocatable,dimension(:,:):: fBnd,eqRatio,msk
    real(kind=double),allocatable,dimension(:,:,:):: tmp
    ! Change of reference surface
    ! ref=sum(low)
    ! ref=sum(up)
    ref=(sum(low)+sum(up))/16.d0/(lmax+1.d0)**2.d0+radius
    ! ----------------------------------------------------------------------------
    nlat=2*(lmax+1)
    nlon=4*(lmax+1)
    
    allocate(effi(2,lmax+1,lmax+1),tmp(2,lmax+1,lmax+1),fBnd(nlat,nlon),eqRatio(nlat,nlon),&
        msk(nlat,nlon)) 
    call binarization(abs(up),msk) !Mask
    effi=0.d0
    call scaling(rhoRatio,eqRatio)
    
    if(flag==0) then ! Normal FSM
        index=0
    else if(flag/=0) then ! Higher order FSM
        index=1
    end if
   
    do i=1,power+index
        fBnd=((up-(ref-radius)*msk)/ref)**i-((low-(ref-radius)*msk)/ref)**i
        call SHExpandDH(fBnd*eqRatio,nlat,tmp,lmaxDH,sampling=2,lmax_calc=lmax)
        do j=1,lmax+1
            if(i>j+2) then
                tmp(:,j,:)=0.d0
            else
                tmp(:,j,:)=tmp(:,j,:)*nchoosek(i,j+2+index)*3.d0/(j+2.d0+index)/(2.d0*j-1.d0)*(ref/radius)**(j)
            end if
        end do
        effi=effi+tmp
    end do
    end subroutine kernal
    
    subroutine ltrGravity(suffix,Rho,lmax1,lowBnd,upBnd,lmax2,power,ltrGrav)
    !****************************************************************************
    ! Perform the lateral Fast Spectral Forward Modelling with constant density
    !  value in a layer.
    !  suffix:       suffix string for generated files
    !  Rho:          spatial variation of density variation in given layer
    !  lmax1:        the equivalent maxium order/degree of Rho
    !  lowBnd,upBnd: upper/lower boundary function
    !  lmax2:        the equivalent maxium order/degree of boundary function
    !  power:        the truncated degree of height binomial expansion
    !  ltrGrav:      the resultant spatial variation of gravity
    ! #Function call: algorithm:: binarization, nchoosek, scaling
    !                 datIO::     admitCorrOut,coeffiOut
    !                 para::      constM,flatten,pi,radius,radiusWGS
    !                 SHTOOLS::   SHAdmitCorr, SHExpandDH, MakeGravGridDh
    !                  
    !****************************************************************************

    character(*),intent(in):: suffix
    integer,intent(in)::power,lmax1,lmax2
    real(kind=double),dimension(:,:),intent(in)::lowBnd,upBnd,Rho	
    real(kind=double),allocatable,dimension(:,:),intent(out)::ltrGrav

    integer:: lmaxBnd,lmaxRho,samples,nlat,nlon,i,j,k
    real(kind=double)::earthRho
    real(kind=double),allocatable,dimension(:):: admitTmp,corrTmp,admitBnd,corrBnd
    real(kind=double),allocatable,dimension(:,:):: gridPhi,gridTheta,gridTotl,avgRho,tmpGrav
    real(kind=double),allocatable,dimension(:,:,:):: effiRho,effiBnd,effiVnm

    allocate(effiRho(2,lmax1+1,lmax1+1),effiBnd(2,lmax2+1,lmax2+1),effiVnm(2,lmax2+1,lmax2+1))
    earthRho=constM*3.d0/4.d0/pi/radius**3.d0
    
    nlat=(lmax2+1)*2
    nlon=(lmax2+1)*4
    allocate(avgRho(nlat,nlon))
    avgRho=sum(Rho)/(nlon*nlat)
  
    call kernal(nchoosek,Rho/earthRho,lowBnd,upBnd,lmax2,power,effiVnm,0)
    call kernal(nchoosek,avgRho/earthRho,lowBnd,upBnd,lmax2,power,effiBnd,0)
    call SHExpandDH(Rho,nlat,effiRho,lmaxRho,lmax_calc=lmax1,sampling=2)
    ! Degree variances
    call coeffiOut('Rho_'//trim(suffix),effiRho,flag=1)
    call coeffiOut('Vnm_'//trim(suffix),effiVnm,flag=1)
    call coeffiOut('Bnd_'//trim(suffix),effiBnd,flag=1)
    !----------------------------------------------------------------------------
    ! Admittance and Correlation spectra of two real functions.
    allocate(admitTmp(lmax1+1),corrTmp(lmax1+1))
    call SHAdmitCorr(effiRho,effiVnm,lmax1,admitTmp,corrTmp)
    call admitCorrOut('Rho_'//trim(suffix),admitTmp,corrTmp)
    
    call SHAdmitCorr(effiBnd,effiVnm,lmax2,admitTmp,corrTmp)
    call admitCorrOut('Bnd_'//trim(suffix),admitTmp,corrTmp)
    deallocate(admitTmp,corrTmp)
    !----------------------------------------------------------------------------
    deallocate(effiRho,effiBnd)
    allocate(tmpGrav(lmax2*2+2,lmax2*4+4),gridPhi(lmax2*2+2,lmax2*4+4),gridTheta(lmax2*2+2,lmax2*4+4), &
             gridTotl(lmax2*2+2,lmax2*4+4),ltrGrav(nlat,nlon))
    call MakeGravGridDH(effiVnm,lmax2,constGM,radius,radiusWGS, &
        flatten,tmpGrav,gridTheta,gridPhi,gridTotl,samples,sampling=2)
    call scaling(tmpGrav,ltrGrav,factor=0.2d0)
    ! In Gal
    deallocate(gridPhi,gridTheta,gridTotl)
    ltrGrav=ltrGrav*1.d5

    end subroutine ltrGravity

    subroutine rdGravity(suffix,polynm,chii,ratio,lmax1,lowBnd,upBnd,lmax2,power,rdGrav)
    !****************************************************************************
    ! Perform the radial Fast Spectral Forward Modelling with density function
    !  of each individual geological province in crust layer.
    !  suffix:       suffix string for generated files
    !  polynm:       density function table
    !  chii:         geological province indicator
    !  ratio:        global scale facotr indicated in GEMMA
    !  lmax1:        the equivalent maxium order/degree of Rho
    !  lowBnd,upBnd: upper/lower boundary function
    !  lmax2:        the equivalent maxium order/degree of boundary function
    !  power:        the truncated degree of height binomial expansion
    !  rdGrav:       the resultant spatial variation of gravity
    ! #Function call: algorithm:: binarization, nchoosek, scaling
    !                 datIO::     coeffiOut
    !                 para::      constM,flatten,pi,radius,radiusWGS
    !                 SHTOOLS::   SHExpandDH, MakeGravGridDh
    !                  
    !****************************************************************************

    character(*),intent(in):: suffix
    integer,intent(in)::power,lmax1,lmax2,chii(:,:)
    real(kind=double),dimension(:,:),intent(in)::lowBnd,upBnd,polynm,ratio
    real(kind=double),allocatable,dimension(:,:),intent(out)::rdGrav

    integer:: lmaxBnd,lmaxRho,samples,nlat,nlon,stg,i,j,s
    real(kind=double)::earthRho
    real(kind=double),allocatable,dimension(:):: admitTmp,corrTmp,admitBnd,corrBnd,dpt,alp0,alp1
    real(kind=double),allocatable,dimension(:,:):: tmpMsk,prvMsk,eqRatio,gridPhi,gridTheta,gridTotl,&
        bnd0,bnd1,ref,tmpGrav
    real(kind=double),allocatable,dimension(:,:,:):: effiVnm,effiVnm1,effiVnm2,stgMsk

    allocate(effiVnm(2,lmax2+1,lmax2+1),effiVnm1(2,lmax2+1,lmax2+1),effiVnm2(2,lmax2+1,lmax2+1))

    earthRho=constM*3.d0/4.d0/pi/radius**3.d0
    nlat=(lmax2+1)*2
    nlon=(lmax2+1)*4

    allocate(prvMsk(nlat,nlon),tmpMsk(nlat,nlon),eqRatio(nlat,nlon),&
        bnd0(nlat,nlon),bnd1(nlat,nlon),ref(nlat,nlon))
    effiVnm=0.d0
    provloop: do i=1,8 !Geo Province Loop
        stg=multStg(i,1)
        call binarization(dble(abs(chii-i)),prvMsk)
        prvMsk=1.d0-prvMsk !Mask for This Province
        allocate(dpt(stg),alp0(stg),alp1(stg),stgMsk(stg,nlat,nlon))
        do s=1,stg
            dpt(s)=polynm(multStg(i,2)+s-1,1)
            alp0(s)=polynm(multStg(i,2)+s-1,2)
            alp1(s)=polynm(multStg(i,2)+s-1,3)
            if (s==1) then
                call binarization(dpt(s)-upBnd+lowBnd,tmpMsk)!Mask for This stageLayer
                stgMsk(s,:,:)=tmpMsk
            else
                call binarization(dpt(s)-upBnd+lowBnd,tmpMsk)!Mask for This stageLayer
                stgMsk(s,:,:)=tmpMsk
                call binarization(dpt(s-1)-upBnd+lowBnd,tmpMsk) !Mask for the latter stageLayer
                stgMsk(s,:,:)=stgMsk(s,:,:)-tmpMsk
            end if
        end do
        stgloop: do s=1,stg
            if(s==1) then
                bnd0=(lowBnd*stgMsk(s,:,:)+(upBnd-dpt(s))*(1.d0-stgMsk(s,:,:)))*prvMsk
                bnd1=upBnd*stgMsk(s,:,:)*prvMsk
                eqRatio=(alp0(s)+alp1(s)*upBnd)/earthRho
                call kernal(nchoosek,eqRatio*ratio-2670.d0/earthRho,bnd0,bnd1,lmax2,power,effiVnm1,0)
                ref=(bnd1+bnd0)/2.d0
                eqRatio=(-alp1(s)*ref)/earthRho
                call kernal(nchoosek,eqRatio*ratio,bnd0,bnd1,lmax2,power,effiVnm2,1)
            else if(s==stg) then
                bnd0=lowBnd*stgMsk(s,:,:)*prvMsk
                bnd1=(upBnd-dpt(s-1))*stgMsk(s,:,:)*prvMsk
                eqRatio=(alp0(s)+alp1(s)*upBnd)/earthRho
                call kernal(nchoosek,eqRatio*ratio-2670.d0/earthRho,bnd0,bnd1,lmax2,power,effiVnm1,0)
                ref=(bnd1+bnd0)/2.d0
                eqRatio=(-alp1(s)*ref)/earthRho
                call kernal(nchoosek,eqRatio*ratio,bnd0,bnd1,lmax2,power,effiVnm2,1)
            else                
                tmpMsk=0.d0
                do j=1,s
                    tmpMsk=tmpMsk+stgMsk(j,:,:)
                end do
                bnd0=(lowBnd*stgMsk(s,:,:)+(upBnd-dpt(s))*(1.d0-tmpMsk))*prvMsk
                bnd1=(upBnd-dpt(s-1))*(1.d0-tmpMsk+stgMsk(s,:,:))*prvMsk
                eqRatio=(alp0(s)+alp1(s)*upBnd)/earthRho
                call kernal(nchoosek,eqRatio*ratio-2670.d0/earthRho,bnd0,bnd1,lmax2,power,effiVnm1,0)
                ref=(bnd1+bnd0)/2.d0
                eqRatio=(-alp1(s)*ref)/earthRho
                call kernal(nchoosek,eqRatio*ratio,bnd0,bnd1,lmax2,power,effiVnm2,1)
            end if
            effiVnm=effiVnm+effiVnm1+effiVnm2
        end do stgloop
        deallocate(dpt,alp0,alp1,stgMsk)
    end do provloop
    deallocate(effiVnm1,effiVnm2,tmpMsk,prvMsk,bnd0,bnd1,ref)
    !Degree variances
    call coeffiOut('Vnm_'//trim(suffix),effiVnm,flag=1)
    allocate(gridPhi(lmax2*2+2,lmax2*4+4), gridTheta(lmax2*2+2,lmax2*4+4), &
        & gridTotl(lmax2*2+2,lmax2*4+4),tmpGrav(lmax2*2+2,lmax2*4+4),rdGrav(nlat,nlon))
    call MakeGravGridDH(effiVnm,lmax2,constGM,radius,radiusWGS, &
        flatten,tmpGrav,gridTheta,gridPhi,gridTotl,samples,sampling=2)
    call scaling(tmpGrav,rdGrav,factor=0.2d0)
    deallocate(gridPhi,gridTheta,gridTotl,tmpGrav)
    
    rdGrav=rdGrav*1.d5 ! m/s2 to mGal
    end subroutine rdGravity
    
    subroutine distCorr(suffix,dist,moho,topo)
    !****************************************************************************
    !
    ! Correlation Coefficients between graivity disturbance and Moho Depth/Solid 
    !  Topography
    !  dist: gravity disturbance
    !  moho: moho depth
    !  topo: solid topography
    ! #Function call: algorithm:: scaling
    !                 datIO::     admitCorrOut,coeffiOut,scatterOut
    !                 para::      constM,flatten,pi,radius,radiusWGS
    !                 SHTOOLS::   SHAdmitCorr, SHExpandDH, MakeGravGridDh
    !
    !****************************************************************************
    character(*),intent(in):: suffix
    real(kind=double),dimension(:,:),intent(in):: dist,moho,topo
    
    integer:: nlon,nlat,lmax,i,j
    real(kind=double),allocatable,dimension(:):: admitTmp,corrTmp,zeros,wi
    real(kind=double),allocatable,dimension(:,:):: mask,tmpTrans,sctTmp
    real(kind=double),allocatable,dimension(:,:,:):: effiDist,effiTmp,scatter
    
    nlat=size(dist,1)
    nlon=size(dist,2)
    lmax=nlat/2-1
        
    allocate(admitTmp(lmax+1),corrTmp(lmax+1),effiDist(2,lmax+1,lmax+1),&
             effiTmp(2,lmax+1,lmax+1),tmpTrans(lmax+1,2*lmax+1))
    
    call SHExpandDH(dist,nlat,effiDist,lmax,sampling=2)
    call coeffiOut('distGravity_'//trim(suffix),effiDist,flag=1)
    call SHExpandDH(moho,nlat,effiTmp,lmax,sampling=2)
    call SHAdmitCorr(effiTmp,effiDist,lmax,admitTmp,corrTmp)
    call admitCorrOut('moho_'//trim(suffix),admitTmp,corrTmp)      
    call SHExpandDH(topo,nlat,effiTmp,lmax,sampling=2)
    call SHAdmitCorr(effiTmp,effiDist,lmax,admitTmp,corrTmp)
    call admitCorrOut('topo_'//trim(suffix),admitTmp,corrTmp)
   
    allocate(mask(nlat,nlon),scatter(2,2,nlat*nlon),sctTmp(2,nlat*nlon))
    
    call binarization(topo,mask)
    ! Scatter map coordianted with gravity disturbance 
    !                         and Moho Depth
    forall(i=1:nlat, j=1:nlon)  
    	scatter(1,1,(j-1)*nlat+i)=dist(i,j)*mask(i,j)
    	scatter(2,1,(j-1)*nlat+i)=dist(i,j)*(1.d0-mask(i,j))
    	scatter(1,2,(j-1)*nlat+i)=moho(i,j)*mask(i,j)
    	scatter(2,2,(j-1)*nlat+i)=moho(i,j)*(1.d0-mask(i,j))
    end forall
    sctTmp=scatter(1,:,:)
    call scatterOut('mohoC_'//trim(suffix),sctTmp)
    sctTmp=scatter(2,:,:)
    call scatterOut('mohoO_'//trim(suffix),sctTmp)   
    ! Scatter map coordianted with gravity disturbance 
    !                         and Solid Topography
    forall(i=1:nlat, j=1:nlon)  
    	scatter(1,1,(j-1)*nlat+i)=dist(i,j)*mask(i,j)
    	scatter(2,1,(j-1)*nlat+i)=dist(i,j)*(1.d0-mask(i,j))
    	scatter(1,2,(j-1)*nlat+i)=topo(i,j)*mask(i,j)
    	scatter(2,2,(j-1)*nlat+i)=topo(i,j)*(1.d0-mask(i,j))
    end forall
    sctTmp=scatter(1,:,:)
    call scatterOut('topoC_'//trim(suffix),sctTmp)
    sctTmp=scatter(2,:,:)
    call scatterOut('topoO_'//trim(suffix),sctTmp) 
    
   deallocate(admitTmp,corrTmp,scatter,mask,sctTmp)
    end subroutine
    
    subroutine rscRho(rho,bnds,modRho,flag)
    !****************************************************************************
    !
    !  If flag=0,  return the topographic mass and geoid below crust
    !     flag/=0, return the Sediments-Consolidated Crust
    !  Mean density is not appropriate when calculating gravity which
    !   isn't proportional to radius.
    !****************************************************************************
    real(kind=double),dimension(:,:,:),intent(in):: rho,bnds
    real(kind=double),allocatable,dimension(:,:,:),intent(out):: modRho
    integer,optional,intent(in):: flag

    integer:: n,lat,lon,i,j,k,l
    real(kind=double),allocatable,dimension(:,:):: tmp1,tmp2

    n=size(rho,1)
    lat=size(rho,2)
    lon=size(rho,3)
    allocate(modRho(2,lat,lon),tmp1(lat,lon),tmp2(lat,lon))

    flagif: if (present(flag) .and. flag==0) then
        ido: do i=1,lat
            jdo: do j=1,lon
                tmp1(i,j)=1.d-280
                tmp2(i,j)=0.d0
                l=3
                ukdo: do k=3,n
                    zerojudg: if(bnds(k,i,j)>=0) then
                        if(k>3) then
                            tmp1(i,j)=tmp1(i,j)+(bnds(k-1,i,j)-bnds(k,i,j))*rho(k-1,i,j)
                        end if
                        l=l+1
                    else
                        if(k==n) then
                            cycle
                        else
                            tmp2(i,j)=tmp2(i,j)+(bnds(k,i,j)-bnds(k+1,i,j))*rho(k,i,j)
                        end if
                    end if zerojudg
                end do ukdo
                tmp1(i,j)=tmp1(i,j)+bnds(l-1,i,j)*rho(l-1,i,j)
                if(l>4) then
                    tmp2(i,j)=tmp2(i,j)-bnds(l,i,j)*rho(l-1,i,j)
                end if
            end do jdo
        end do ido
        modRho(1,:,:)=(bnds(3,:,:)+abs(bnds(3,:,:)))/2.d0/tmp1
        where(modRho(1,:,:)>0)
            modRho(1,:,:)=1.d0/modRho(1,:,:)
        end where
        modRho(2,:,:)=tmp2/((bnds(3,:,:)-abs(bnds(3,:,:)))/2.d0-bnds(9,:,:))
    else
        do i=1,2
            tmp1=1.d-280
            do j=1,3
                tmp1=tmp1+(bnds(i*3+j,:,:)-bnds(i*3+j-1,:,:))*rho(i*3+j-1,:,:)
                if(i==1) then
                    modRho(i,:,:)=(bnds(i*3+3,:,:)-bnds(i*3,:,:))/tmp1
                else
                    modRho(i,:,:)=tmp1/(bnds(i*3+3,:,:)-bnds(i*3,:,:))
                end if
            end do
        end do
        do i=1,lat
            do j=1,lon
                if (modRho(1,i,j)/=0.d0) then
                    modRho(1,i,j)=1.d0/modRho(1,i,j)
                end if
            end do
        end do
    end if flagif
    end subroutine

    end module analysis