!======================================================================-
!Create dataset for CNN
!
!
!======================================================================-
program FitzHugh_Nagumo_param
    !======================================================================-
    ! 宣言
    !======================================================================-
    implicit none
    integer,parameter                   ::  size=64, dx2=16384
    double precision,parameter          ::  dt=6.103515625
    double precision                    ::  u(size,size), u(size:size),a,b,c

    function f(u,v)
        implicit none
        double precision    :: u, v
        f=u-u**3-v
    end function f

    function g(u,v)
        implicit none
        double precision    :: u, v, a, b, c
        g=c*(u-a*v-b)
    end function g

    function NL_f(u,v,f) result()
        implicit none
        integer             ::  i,j
        double precision    ::  u, v
        double precision,parameter  ::  size=size(u)
        double precision    ::  out(size,size)

        do i = 1, size
            do j= 1, size
                out(i,j)=f(u(i,j),v(i,j))
            enddo
        enddo 
        
    end function NL_f

    function NL_f(u::Array{Float64,2},v::Array{Float64,2})
        out=similar(u)
        for i in  eachindex(u)  i:: Int64
            out[i]=f(u[i],v[i])
        end
        return out
    end
    
    function NL_g(u::Array{Float64,2},v::Array{Float64,2}, model)
        out=similar(u)
        for i in  eachindex(u)  i:: Int64
            out[i]=g(u[i],v[i],model)
        end
        return out
    end
    
    function diffusion(system,distribution)
        out=similar(distribution)
        @. @inbounds out=-4*distribution
        #diff x
            @. @inbounds out[1:system.size-1,:]+=@view distribution[2:system.size,:]
            @. @inbounds out[system.size,:]+=@view distribution[1,:]
            @. @inbounds out[2:system.size,:]+=@view distribution[1:system.size-1,:]
            @. @inbounds out[1,:]+=@view distribution[system.size,:]
        #diff y
            @. @inbounds out[:,1:system.size-1]+=@view distribution[:,2:system.size]
            @. @inbounds out[:,system.size]+=@view distribution[:,1]
            @. @inbounds out[:,2:system.size]+=@view distribution[:,1:system.size-1]
            @. @inbounds out[:,1]+=@view distribution[:,system.size]
        return  out*system.invdx
    end
    
    function onestep_u(model,system,variable)
        @inbounds variable.u+=system.dt*(model.Du*diffusion(system,variable.u)+NL_f(variable.u,variable.v))
        return variable.u
    end
    
    function onestep_v(model,system,variable)
        @inbounds variable.v+=system.dt*(model.Dv*diffusion(system,variable.v)+NL_g(variable.u,variable.v,model))
        return variable.v
    end
    
    function onestep(model,system,variable)
        variable.u=onestep_u(model,system,variable)
        variable.v=onestep_v(model,system,variable)
        return variable
    end
    
































































    
    if((lsize.gt.lmax).or.(lsize.lt.1)) then
      write(*,*) 'error in lsize:',lsize
      stop
    endif
    if((inist.lt.-1).or.(inist.gt.1)) then
      write(0,*) 'error in inist:',inist
      stop
    endif
    if((iseed.le.0).or.(mod(iseed,2).eq.0)) then
      write(0,*) 'error in iseed:',iseed
      stop
    endif
    if(imcs.lt.0) then
      write(0,*) 'error in imcs:',imcs
      stop
    endif
    if(nmcs.lt.0) then
      write(0,*) 'error in nmcs:',nmcs
      stop
    endif
    if(nblock.lt.0) then
      write(0,*) 'error in nblock:',nblock
      stop
    endif
    !==========================================================-
    !   shutsuryoku
    !==========================================================-
    write(*,'(a10,i12)')    '   lsize= ',lsize
    write(*,'(a10,e24.16)') '    temp= ',temp
    write(*,'(a10,i12)')    '   inist= ',inist
    write(*,'(a10,i12)')    '    imcs= ',imcs
    write(*,'(a10,i12)')    '    nmcs= ',nmcs
    write(*,'(a10,i12)')    '  nblock= ',nblock
    write(*,'(a10,i12)')    '   iseed= ',iseed
    !======================================================================-
    ! shokisettei
    !======================================================================-
    !   ransu shokika
    !==========================================================-
    call rndini(iseed)
    !==========================================================-
    !   tonari no hyo shokika
    !==========================================================-
    no=lsize*lsize
    lsize1=lsize-1
    do i=0,lsize1
      ip(i)=i+1
      im(i)=i-1
    enddo
    ip(lsize1)=0
    im(0)=lsize1
    !==========================================================-
    !   spin haichi shokika
    !==========================================================-
    do ix=0,lsize1
      do iy=0,lsize1
        if(inist.ge.1) then
          ispin(ix,iy)=1
        else if(inist.le.-1) then
          ispin(ix,iy)=-1
        else
          call rndu(rnd)
          if(rnd.ge.0.5d0) then
            ispin(ix,iy)=1
          else
            ispin(ix,iy)=-1
          endif
        endif
      enddo
    enddo
    !==========================================================-
    !   sen'i-kakuritsu shokika
    !==========================================================-
    do n=-nnno,nnno
      de=2.0d0*dble(n)/temp
      if(de.gt.0.0d0) then
        trprob(n)=dexp(-de)
      else
        trprob(n)=1.0d0
      endif
    enddo
    !======================================================================-
    ! shoki loop
    !======================================================================-
    do nm=1,imcs
      !==========================================================-
      !   1 Monte Carlo step
      !==========================================================-
      do n=1,no
        call rndu(rnd)
        ix=lsize*rnd
        call rndu(rnd)
        iy=lsize*rnd
        ie=ispin(ix,iy)*(ispin(ip(ix),   iy )+ispin(im(ix),   iy ) &
             +ispin(   ix ,ip(iy))+ispin(   ix ,im(iy)))
        call rndu(rnd)
        if(rnd.lt.trprob(ie)) then
          ispin(ix,iy)=-ispin(ix,iy)
        endif
      enddo
    enddo
    !======================================================================-
    ! block loop
    !======================================================================-
    do nb=1,nblock
      !==========================================================-
      !   block-wa shokika
      !==========================================================-
      amag =0.0d0
      amag2=0.0d0
      aene =0.0d0
      aene2=0.0d0
      !==========================================================-
      !   Monte Carlo step loop
      !==========================================================-
      do nm=1,nmcs
        !==============================================-
        !     1 Monte Carlo step
        !==============================================-
        do n=1,no
          call rndu(rnd)
          ix=lsize*rnd
          call rndu(rnd)
          iy=lsize*rnd
          ie=ispin(ix,iy)*(ispin(ip(ix),   iy )+ispin(im(ix),   iy ) &
               +ispin(   ix ,ip(iy))+ispin(   ix ,im(iy)))
          call rndu(rnd)
          if(rnd.lt.trprob(ie)) then
            ispin(ix,iy)=-ispin(ix,iy)
          endif
        enddo
        !==============================================-
        !     sokutei
        !==============================================-
        iene=0
        imag=0
        do ix=0,lsize1
          do iy=0,lsize1
            imag=imag+ispin(ix,iy)
            iene=iene+ispin(ix,iy)*(ispin(ip(ix),iy)+ispin(ix,ip(iy)))
          enddo
        enddo
        !==============================================-
        !     block-wa kasan
        !==============================================-
        aene =aene +dble(iene)
        aene2=aene2+dble(iene)*dble(iene)
        amag =amag +dble(imag)
        amag2=amag2+dble(imag)*dble(imag)
      enddo
      !==========================================================-
      !   block heikin
      aene =aene /dble(nmcs)
      aene2=aene2/dble(nmcs)
      amag =amag /dble(nmcs)
      amag2=amag2/dble(nmcs)
      write(*,'(a10,i12)')    '   block= ',nb
      write(*,'(a10,e24.16)') '    aene= ',aene
      write(*,'(a10,e24.16)') '   aene2= ',aene2
      write(*,'(a10,e24.16)') '    amag= ',amag
      write(*,'(a10,e24.16)') '   amag2= ',amag2
      dene =aene2-aene*aene
      ave  =-aene/dble(no)
      avc  =(dene/dble(no))/(temp*temp)
      dmag =amag2-amag*amag
      avm  =amag/dble(no)
      avx  =dmag/dble(no)
      write(*,'(a10,e24.16)') '       e= ',ave
      write(*,'(a10,e24.16)') '       c= ',avc
      write(*,'(a10,e24.16)') '       m= ',avm
      write(*,'(a10,e24.16)') '       x= ',avx
    enddo
    !
    stop
  end program ising
  