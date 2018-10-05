!******************************************************************************************************************************************************************************************************
!***********************************************************************************VARIABLES**********************************************************************************************************
!******************************************************************************************************************************************************************************************************
!CONTAINS ALL THE GLOBAL VARIABLES
!MOST OF THE NAMES ARE SELF-EXPLAINATORY, IF NOT, SEE THE SUBROUTINES THEY ARE USED IN TO UNDERSTAND THEIR MEANINGS

module variables
    
    implicit none
    save

    integer                                  :: final_time,iter,interval,vis_flag,finish,periodic_x,periodic_y,periodic_z,duc_count
    integer                                  :: firstx,lastx,firsty,lasty,firstz,lastz,startx,starty,startz,endx,endy,endz
    integer                                  :: startx_vis,starty_vis,startz_vis,endx_vis,endy_vis,endz_vis,startx_print,starty_print,startz_print,endx_print,endy_print,endz_print
    integer                                  :: myrank,mperr
    integer,parameter                        :: sizex=6,sizey=100,sizez=100,zprocs=2,yprocs=2
    !each processor must have at least 6 points in the directions in which tthe boundary conditions are not periodic(this might be obvious or redundant but still) 
    !Not sure what is the thinking behind above comment. I see it now, 3 points per reocessor are needed irrespective of BCs
    integer,parameter                        :: sz=sizez/zprocs+6,sy=sizey/yprocs+6,sx=sizex+6,nprocs=zprocs*yprocs
    integer,parameter                        :: n_sp=5,n_di=3,chem_flag=1,vib_flag=1 !number of species,diatomic species and the two flags
    !n_di is only used in the case of vibrational terms. It can be incorporated in Cp and Cv but they are specified explicitly
    !keep n_sp=1 and n_di=1 when ideal gas case is used, or any other model which doesn't involve dealing with seperate species
    !make sure n_di<=n_sp
    integer,parameter                        :: neq=n_sp+4+vib_flag
    double precision,parameter               :: pi=2.0d0*asin(1.0d0),cl_para=1000d0!cl_para is the clustering parameter in y direction
    double precision                         :: alpha_x_global(1,sy,sz,neq),alpha_y(sx,0:(zprocs-1),sz,neq),alpha_y_global(sx,0:(zprocs-1),sz,neq),alpha_z(sx,sy,0:(yprocs-1),neq),alpha_z_global(sx,sy,0:(yprocs-1),neq)
    double precision,dimension(sx,sy,sz,neq) :: q,q0,ddx_fluxf,ddy_fluxg,ddz_fluxh,visf,visg,vish,fluxf,fluxg,fluxh
    double precision,dimension(sx,sy,sz)     :: x,y,z,y_ph,dy_dy_ph,T,Tv,ro,p,u,v,w,ducros,vis,cond,cond_v,vortcty,heatx,heaty,heatz!,m_check
    double precision,dimension(n_di)         :: theta_v
    double precision,dimension(n_sp)         :: Mw,Cp,Cv,Cv_tr,h0,At,Bt,Ct
    double precision,dimension(17)           :: cf,eta,theta,cm,a1,a2,a3,a4,a5
    double precision                         :: gamma_ref,gamma_star,R_u,R_star,R_ref,Re,Mach,prandtl,dx,dy,dz,dxin,dyin,dzin,dt,cfl,cptim,walltim1,walltim2,lenc,lenx,leny,lenz,lenz_in,lo_star,co_star,co_ref,ro_ref,ro_star,vis_ref,vis_star,T_ref,T_star,Lewis
    double precision                         :: vel_star,p_star,T_wall,Tv_star,E_star,Ev_star,del_h0,intial_turb_KE,t_turb,k_eq_star_eq
    
end module variables

!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************










!******************************************************************************************************************************************************************************************************
!********************************************************************************MAIN PROGRAM**********************************************************************************************************
!******************************************************************************************************************************************************************************************************
  
program heisenberg

    use variables
    implicit none
    include 'mpif.h'

    call MPI_INIT(mperr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, mperr)
    walltim1 = MPI_Wtime()

    call constants()
    call intialise()
    call data_in()
    call compute()
    
    call MPI_FINALIZE(mperr)

end program heisenberg

!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************










!******************************************************************************************************************************************************************************************************
!***********************************************************************SUBROUTINE FOR TIME INTEGRATION OF NS******************************************************************************************
!******************************************************************************************************************************************************************************************************
!q CONTAINS CONSERVATIVE VARIABLES
!THE 3D NS EQUATIONS IN THE FORM dq/dt + dF/dx + dG/dy + dH/dz = dFv/dx + dGv/dy + dHv/dz  SOLVED
!IN heisenberg, F=fluxf,ddx_fluxf=dF/dx,Fv=visf AND SO ON 
!q0 TEMPORARILY STORES q OF PRESENT ITERATION IN TIME INTEGRATION
!MODIFIED, 4TH ORDER EXPLICIT RUNGE-KUTTA TIME INEGRATION IS USED HERE. REFER 1ST VOLUME OF HOFFMANN FOR DETAILS

subroutine compute()

    use variables
    implicit none
    include 'mpif.h'
    integer                                    :: i,j,k,m,l,n,flag=0,flag_global=0,test_time
    double precision                           :: d6dx_vect,d6dy_vect,d6dz_vect,d6dx,d6dy,d6dz
    double precision, dimension(sx,sy,sz,n_sp) :: ws
    double precision, dimension(sx,sy,sz)      :: wv
    double precision, dimension(sx,sy,sz,neq)  :: convergence=0d0

    iter=0
    l=0
  
    q0(:,:,:,:) = q(:,:,:,:)
    ws(:,:,:,:) = 0d0
    wv(:,:,:) = 0d0

    call bc(l)
    call message_pass(l)
    call update(l)
    call bc_viscous()
    call message_pass_viscous(l)
    if (chem_flag==1) call source(ws)
    if (vib_flag==1) call vib_source(ws,wv)
    ! call stats()
    call output(iter,ws)

    test_time = final_time - 1
    time:do iter=1,final_time
        flag = 0
        do l=4,1,-1
            ! m_check(:,:,:) = 0d0
            do m=1,neq
                do k=startz,endz
                    do j=starty,endy
                        do i=startx,endx
                            if(m<=n_sp) then
                                q(i,j,k,m) = q0(i,j,k,m) - 1.0d0/l*dt*(ddx_fluxf(i,j,k,m) + ddy_fluxg(i,j,k,m) + ddz_fluxh(i,j,k,m)) + &
                                                           1.0d0/l*dt*(d6dx_vect(visf,i,j,k,m) + d6dy_vect(visg,i,j,k,m) + d6dz_vect(vish,i,j,k,m)) + &
                                                           1.0d0/l*dt*ws(i,j,k,m)*chem_flag
                                ! m_check(i,j,k) = m_check(i,j,k) + d6dx_vect(visf,i,j,k,m) + d6dy_vect(visg,i,j,k,m) + d6dz_vect(vish,i,j,k,m) + ws(i,j,k,m)
                            else if(m<=n_sp+4) then
                                q(i,j,k,m) = q0(i,j,k,m) - 1.0d0/l*dt*(ddx_fluxf(i,j,k,m) + ddy_fluxg(i,j,k,m) + ddz_fluxh(i,j,k,m)) + &
                                                           1.0d0/l*dt*(d6dx_vect(visf,i,j,k,m) + d6dy_vect(visg,i,j,k,m) + d6dz_vect(vish,i,j,k,m))
                            else
                                q(i,j,k,m) = q0(i,j,k,m) - 1.0d0/l*dt*(ddx_fluxf(i,j,k,m) + ddy_fluxg(i,j,k,m) + ddz_fluxh(i,j,k,m)) + &
                                                           1.0d0/l*dt*(d6dx_vect(visf,i,j,k,m) + d6dy_vect(visg,i,j,k,m) + d6dz_vect(vish,i,j,k,m)) + &
                                                           1.0d0/l*dt*wv(i,j,k)
                            end if
                        end do
                    end do
                end do
            end do

            call bc(l)
            call message_pass(l)
            call update(l)
            if (chem_flag==1) call source(ws)
            if (vib_flag==1) call vib_source(ws,wv)
            call bc_viscous()
            call message_pass_viscous(l)
        end do

        ! if (mod(iter,100)==0) then
        !     ! if(norm2((q(:,:,:,:)-q0(:,:,:,:))/q0(:,:,:,:))<1e-6) flag = 1
        !     if(maxval(abs((q(:,:,:,:)-q0(:,:,:,:))/q0(:,:,:,:)))<1e-6) flag = 1
        !     CALL MPI_ALLREDUCE(flag,flag_global,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mperr)
        !     if (flag_global==4) then
        !         write(*,*)'CONVERGED'
        !         ! call output(iter)                        
        !     end if
        ! end if

        if(myrank==0) write(*,*) iter
        
        q0(:,:,:,:) = q(:,:,:,:)
     
        ! call stats()
        ! if(iter<4000) then
            if (mod(iter,interval)==0) call output(iter,ws)
        ! else if(iter<8000) then
        !     if (mod(iter,100)==0) call output(iter,ws)
        ! else
        !     if (mod(iter,1000)==0) call output(iter,ws)
        ! end if
        ! call output(iter,ws)

        ! if(test_time>final_time) stop   !if stop or even pause is used following if condition, one of the processors always interrupts output command of others 
        ! if(flag_global==4) test_time = final_time + 1

    end do time

end subroutine compute

!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************









!******************************************************************************************************************************************************************************************************
!*****************************************************************SUPPLY THE INPUT PARAMETERS INCLUDING INITIAL CONDITIONS*****************************************************************************
!******************************************************************************************************************************************************************************************************
!star quantities are all quantities at the initial stage
!ref quantities are reference quantities based on 300K temperature or something else that depends on the quantity

subroutine data_in()   

    use variables
    implicit none
    character(len=40)                 :: filename
    integer                           :: i,j,k,i1,j1,k1,m,range1,range2,numrc1,numrc2,numrc3,numrc4,statu,ierr,restart
    double precision                  :: p1,ratio,vortex,radius,delta,dum,ignore,ignore1,add,R_s,Z_s,Cp_star,cs1_star_eq,cs2_star_eq
    double precision, dimension(n_sp) :: cs_star,ro_s_star
    double precision, dimension(n_di) :: ev_s_star,ev_s

    
    !********** reference quantities **********
    R_ref = 287.15d0
    gamma_ref = 1.4d0
    ro_ref = 1.20d0
    T_ref = 300d0
    vis_ref = 0.000019830d0
    co_ref = sqrt(gamma_ref*R_ref*T_ref)
    !******************************************

    

    
    !********** flow parameters **********
    Mach = 1.1
    prandtl = 0.7d0
    Re = 515d0 !can be changed to match to match a paricular initial Re_lambda in isoturb
    Lewis = 1d0
    !*************************************




    !********** initial quantities **********
    T_star = 5000d0
    ro_star = ro_ref
    cs_star(1) = 0.768d0
    cs_star(2) = 0.232d0
    cs_star(3) = 0d0
    cs_star(4) = 0d0
    cs_star(5) = 0d0
    !****************************************




    !******************** quantities dependent on initial and/or reference quantities ********************
    R_star = R_u*add(cs_star,Mw,n_sp,-1)
    Cp_star = add(cs_star,Cp,n_sp,1)
   
    p_star = ro_star*R_star*T_star

    ro_s_star(1) = cs_star(1)*ro_star
    ro_s_star(2) = cs_star(2)*ro_star
    ro_s_star(3) = cs_star(3)*ro_star
    ro_s_star(4) = cs_star(4)*ro_star
    ro_s_star(5) = cs_star(5)*ro_star

    gamma_star = 1d0 + ro_star*R_star/add(ro_s_star,Cv,n_sp,1)!depending on the case you are going to run, gamma can be a function of space.
    
    co_star = sqrt(gamma_star*R_star*T_star)    
    vis_star = vis_ref*(T_star/T_ref)**0.670d0
    lo_star = Re*vis_star/ro_star/co_star
    vel_star = Mach*co_star
    Tv_star = T_star
    T_wall = T_star
    
    ! Z_s = 10000d0/T_star
    ! k_eq_star_eq = cm(1)*exp(a1(1)/Z_s + a2(1) + a3(1)*log(Z_s) + a4(1)*Z_s + a5(1)*Z_s**2)
    ! cs1_star_eq = 1d0/(1d0 + k_eq_star_eq)
    ! cs2_star_eq = k_eq_star_eq/(1d0 + k_eq_star_eq)

    ev_s_star(:) = vib_flag*((R_u/Mw(1:n_di))*theta_v(:)/(exp(theta_v(:)/Tv_star)-1))
    !*****************************************************************************************************
    


    
    !********** time quantities **********
    dt = 0.020d0*lo_star/co_star
    final_time = 2000
    interval = 50
    !*************************************




    !********** domain size and grid **********
    lenc = 1d0
    leny = 2d0*pi*lo_star
    lenx = 2d0*pi*lo_star
    lenz = 2d0*pi*lo_star
    
    dxin = lenx/(sizex-1)
    dyin = leny/(sizey-1) 
    dzin = lenz/(sizez-1)

    dx = dxin
    dy = lenc/(sizey-1)
    dz = dzin

    call gridin()
    !******************************************




    !**************************************** periodicity and limits for do loops ****************************************
    vis_flag = 1
    finish = 0
    periodic_x = 1 !1 FOR PERIODIC BCs AND 0 FOR NON-PERIODIC
    periodic_y = 1
    periodic_z = 1
    restart = 0    

    startx = 4
    endx = sx - 3
    starty = 4
    endy = sy - 3
    startz = 4
    endz = sz - 3

    firstx = 1       !THESE WILL CHANGE ONLY IF BOUNDARY CONDITIONS ARE NOT PERIODIC
    lastx = sx
    firsty = 1
    lasty = sy
    firstz = 1
    lastz = sz
    
    if(periodic_x==0) then
        firstx = startx - 1
        lastx = endx + 1
    end if
    if(periodic_y==0) then
        if(myrank/zprocs==0) firsty = starty - 1
        if(myrank/zprocs==yprocs-1) lasty = endy + 1        
    end if
    if(periodic_z==0) then
        if(mod(myrank,zprocs)==0) firstz = startz - 1
        if(mod(myrank,zprocs)==zprocs-1) lastz = endz + 1
    end if

    startx_vis = startx       !THESE WILL CHANGE ONLY IF BOUNDARY CONDITIONS ARE NOT PERIODIC
    endx_vis = endx
    starty_vis = starty
    endy_vis = endy
    startz_vis = startz
    endz_vis = endz
    
    if(periodic_x==0) then
        startx_vis = firstx 
        endx_vis = lastx
    end if
    if(periodic_y==0) then
        if(myrank/zprocs==0) starty_vis = firsty
        if(myrank/zprocs==yprocs-1) endy_vis = lasty        
    end if
    if(periodic_z==0) then
        if(mod(myrank,zprocs)==0) startz_vis = firstz
        if(mod(myrank,zprocs)==zprocs-1) endz_vis = lastz
    end if
  
    startx_print = startx 
    endx_print = endx
    starty_print = starty
    endy_print = endy
    startz_print = startz
    endz_print = endz
    !*********************************************************************************************************************




    !****************************************** initialisation of flow variables *****************************************
    vortex=0.20d0
    radius=0.000020d0

    Tv(:,:,:) = Tv_star
    do k=firstz,lastz
        do j=firsty,lasty
            do i=firstx,lastx
                ratio = dsqrt((z(i,j,k) - lenz/2d0)**2 + (y_ph(i,j,k) - leny/2d0)**2)/radius
                v(i,j,k) = vel_star*(vortex*(z(i,j,k) - lenz/2d0)/radius*exp(-ratio**2/2))
                w(i,j,k) = vel_star*(1d0 - vortex*(y_ph(i,j,k) - leny/2d0)/radius*exp(-ratio**2/2))
                u(i,j,k) = 0d0
                T(i,j,k) = T_star - 0.50d0*(vortex*vel_star*exp(-ratio**2/2))**2/Cp_star
                ro(i,j,k) = ro_star*(T(i,j,k)/T_star)**(1/(gamma_star - 1d0))
                p(i,j,k) = ro(i,j,k)*R_star*T(i,j,k)
                Tv(i,j,k) = T(i,j,k)
                ev_s(:) = vib_flag*((R_u/Mw(1:n_di))*theta_v(:)/(exp(theta_v(:)/Tv(i,j,k))-1))
                if(n_sp/=1) q(i,j,k,1:n_sp)  = cs_star(1:n_sp)*ro(i,j,k)
                if(n_sp==1) q(i,j,k,n_sp)  = ro(i,j,k)
                q(i,j,k,n_sp+1)  = ro(i,j,k)*u(i,j,k)
                q(i,j,k,n_sp+2)  = ro(i,j,k)*v(i,j,k)
                q(i,j,k,n_sp+3)  = ro(i,j,k)*w(i,j,k)
                q(i,j,k,n_sp+4)  = add(q(i,j,k,1:n_sp),Cv,n_sp,1)*T(i,j,k) + 0.50d0*ro(i,j,k)*(u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2) + add(q(i,j,k,1:n_di),ev_s,n_di,1)*vib_flag + add(q(i,j,k,1:n_sp),h0,n_sp,1)*chem_flag
                if (vib_flag==1) q(i,j,k,n_sp+5) = add(q(i,j,k,1:n_di),ev_s,n_di,1)
            end do
        end do
    end do
    !****************************************** initialisation of flow variables *****************************************

   
end subroutine data_in

!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************










!******************************************************************************************************************************************************************************************************
!************************************************************************************OUTPUT TO DAT FILE************************************************************************************************
!******************************************************************************************************************************************************************************************************

subroutine output(nth,ws)

    use variables
    implicit none
    include 'mpif.h'
    integer            :: i,j,k,uni,ierr,digit1,digit2,digit3,digit4,digit5,digit6,digit7,digit8,digit9,digit10,digit11,digit12,numrc1,numrc2,numrc3,numrc4,criteria
    double precision, dimension(sx,sy,sz,n_sp), intent(in) :: ws
    integer,intent(in) :: nth 
    character(len=40)  :: filename
    character(len=8)   :: date
    character(len=10)  :: time
    double precision   :: check,vortctyz,vortctyy,vortctyx,d6dx,d6dy,d6dz
    double precision   :: drvtx,drvty,drvtz,total_vel,local_a,local_mach

    digit1 = nth/1000000
    digit2 = mod(nth,1000000)
    digit3 = digit2/100000
    digit4 = mod(digit2,100000)
    digit5 = digit4/10000
    digit6 = mod(digit4,10000)
    digit7 = digit6/1000
    digit8 = mod(digit6,1000)
    digit9 = digit8/100
    digit10 = mod(digit8,100)
    digit11 = digit10/10
    digit12 = mod(digit10,10)

    numrc1 = myrank/100
    numrc2 = mod(myrank,100)
    numrc3 = numrc2/10
    numrc4 = mod(numrc2,10)

    filename='neq'//char(nprocs+48)//'_'//char(digit1+48)//char(digit3+48)//char(digit5+48)//char(digit7+48)//char(digit9+48)//char(digit11+48)//char(digit12+48)//'p'//char(numrc1+48)//char(numrc3+48)//char(numrc4+48)//'.dat'
    open(unit=40,file=filename,STATUS='REPLACE',ACTION='WRITE',IOSTAT=ierr)
    if(ierr/=0) write(*,*)filename,ierr
    write(40,'(a180)') 'variables="z","y","x","ro1","ro2","ro3","ro4","ro5","ro","p","T","Tv","u","v","w"'
    write(40,*) 'zone I=',endz_print-startz_print+1, 'J=',endy_print-starty_print+1, 'K=',endx_print-startx_print+1,' F=POINT'
    do i=startx_print,endx_print
        do j=starty_print,endy_print
            do k=startz_print,endz_print
                write(40,'(15E20.10)') z(i,j,k)/lo_star,y_ph(i,j,k)/lo_star,x(i,j,k)/lo_star,q(i,j,k,1:n_sp)/ro_star,ro(i,j,k)/ro_star,p(i,j,k)/p_star,T(i,j,k)/T_star,Tv(i,j,k)/Tv_star,u(i,j,k)/co_star,v(i,j,k)/co_star,w(i,j,k)/co_star
            end do
        end do
    end do
    close(40)

    filename='derivs'//char(nprocs+48)//'_'//char(digit1+48)//char(digit3+48)//char(digit5+48)//char(digit7+48)//char(digit9+48)//char(digit11+48)//char(digit12+48)//'p'//char(numrc1+48)//char(numrc3+48)//char(numrc4+48)//'.dat'
    open(unit=40,file=filename,STATUS='REPLACE',ACTION='WRITE',IOSTAT=ierr)
    if(ierr/=0) write(*,*)filename,ierr
    write(40,'(a180)') 'variables="z","y","x","dtdx","dtdy","dtdz","dtvdx","dtvdy","dtvdz"'
    write(40,*) 'zone I=',endz_print-startz_print+1, 'J=',endy_print-starty_print+1, 'K=',endx_print-startx_print+1,' F=POINT'
    do i=startx_print,endx_print
        do j=starty_print,endy_print
            do k=startz_print,endz_print
                write(40,'(9E20.10)') z(i,j,k)/lo_star,y_ph(i,j,k)/lo_star,x(i,j,k)/lo_star,d6dx(T,i,j,k),d6dy(T,i,j,k),d6dz(T,i,j,k),d6dx(Tv,i,j,k),d6dy(Tv,i,j,k),d6dz(Tv,i,j,k)
            end do
        end do
    end do
    close(40)

    if(iter==0 .and. myrank==0) then
        open(unit=30,file='1_data',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=ierr)
        write(30,*)Mw(1),Mw(2),'Mw(1),Mw(2)'
        write(30,*)dt,'dt'
        write(30,*)Mach,'Mach'
        write(30,*)R_star,'R_star'
        write(30,*)co_star,'co_star'
        write(30,*)vel_star,'vel_star'
        write(30,*)p_star,'p_star'
        write(30,*)T_star,'T_star'
        write(30,*)ro_star,'ro_star'
        write(30,*)n_sp,'n_sp'
        write(30,*)chem_flag,'chem_flag'
        ! write(30,*)'chem_flag=1 but source is commented to see if shock remains stationary'
        close(30)
    end if


    ! filename='Time.dat'
    ! open(unit=30,file=filename,STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=ierr)
    ! call CPU_TIME(cptim)
    ! write(30,*) 'CPU Time taken =',cptim
    ! walltim2=MPI_Wtime()
    ! write(30,*) 'Time taken =',walltim2-walltim1
    ! close(30)
    
    ! call CPU_TIME(cptim)
    ! write(*,*) 'CPU Time taken =',cptim
    walltim2=MPI_Wtime()
    if(myrank==0) write(*,*) 'Time taken =',walltim2-walltim1
    
    ! call DATE_AND_TIME(DATE=date,TIME=time)
    ! if(myrank==0) write(*,*)'output for',nth,'      on  ',date(7:8),'-',date(5:6),'-',date(1:4),'  at  ',time(1:2),':',time(3:4),':',time(5:6),':',time(8:10)

    ! if(myrank==0) write(*,*)'output for',nth

end subroutine output

!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************










!******************************************************************************************************************************************************************************************************
!************************************************************************************OUTPUT TO DAT FILE************************************************************************************************
!******************************************************************************************************************************************************************************************************

! subroutine stats()   

!     use variables
!     implicit none
!     include 'mpif.h'
!     character(len=40)           :: filename
!     integer                     :: i,j,k,i1,j1,k1,m,range1,range2,numrc1,numrc2,numrc3,numrc4,statu,ierr,restart,WENO
!     real*8                      :: ro1,p1,ratio,vortx,radius,delta,dum,global_average,d6dx,d6dy,d6dz
!     real*8                      :: u_bar,v_bar,w_bar,ro_bar,ro1_bar,ro2_bar,p_bar,T_bar,vis_bar,u_p_sq_bar,v_p_sq_bar,w_p_sq_bar,ro_p_sq_bar,p_p_sq_bar,T_p_sq_bar,cs1_bar,cs2_bar
!     real*8                      :: lambda,lambda_x,lambda_y,lambda_z,skew,skew_x,skew_y,skew_z,turb_KE,Re_lambda,M_t
!     real*8, dimension(sx,sy,sz) :: u_p,v_p,w_p,ro_p,p_p,T_p,dudx,dvdy,dwdz,cs1,cs2

!     cs1(:,:,:) = q(:,:,:,1)/ro(:,:,:)
!     cs2(:,:,:) = q(:,:,:,2)/ro(:,:,:)
!     u_bar = global_average(u,1)
!     v_bar = global_average(v,1)
!     w_bar = global_average(w,1)
!     ro_bar = global_average(ro,1)
!     ro1_bar = global_average(q(:,:,:,1),1)
!     ro2_bar = global_average(q(:,:,:,2),1)
!     cs1_bar = global_average(cs1,1)
!     cs2_bar = global_average(cs2,1)
!     p_bar = global_average(p,1)
!     T_bar = global_average(T,1)
!     vis_bar = global_average(vis,1)

!     do k=startz,endz
!         do j=starty,endy
!             do i=startx,endx
!                 u_p(i,j,k) = u(i,j,k) - u_bar
!                 v_p(i,j,k) = v(i,j,k) - v_bar
!                 w_p(i,j,k) = w(i,j,k) - w_bar
!                 ro_p(i,j,k) = ro(i,j,k) - ro_bar
!                 p_p(i,j,k) = p(i,j,k) - p_bar
!                 T_p(i,j,k) = T(i,j,k) - T_bar
!                 dudx(i,j,k) = d6dx(u,i,j,k)
!                 dvdy(i,j,k) = d6dy(v,i,j,k)
!                 dwdz(i,j,k) = d6dz(w,i,j,k) 
!             end do
!         end do
!     end do

!     u_p_sq_bar = global_average(u_p,2)
!     v_p_sq_bar = global_average(v_p,2)
!     w_p_sq_bar = global_average(w_p,2)
!     ro_p_sq_bar = global_average(ro_p,2)
!     p_p_sq_bar = global_average(p_p,2)
!     T_p_sq_bar = global_average(T_p,2) 

!     if(iter==0) intial_turb_KE = u_p_sq_bar + v_p_sq_bar + w_p_sq_bar
!     turb_KE = u_p_sq_bar + v_p_sq_bar + w_p_sq_bar
!     M_t = dsqrt(turb_KE)/co_star

!     lambda_x = dsqrt(u_p_sq_bar/global_average(dudx,2))
!     lambda_y = dsqrt(v_p_sq_bar/global_average(dvdy,2))
!     lambda_z = dsqrt(w_p_sq_bar/global_average(dwdz,2))
!     lambda = (lambda_x + lambda_y + lambda_z)/3d0

!     skew_x = global_average(dudx,3)/(global_average(dudx,2))**1.5d0
!     skew_y = global_average(dvdy,3)/(global_average(dvdy,2))**1.5d0
!     skew_z = global_average(dwdz,3)/(global_average(dwdz,2))**1.5d0
!     skew = (skew_x + skew_y + skew_z)/3d0

!     Re_lambda = lambda*sqrt(turb_KE/3d0)*ro_bar/vis_bar     !Re*(lambda/lo_star)*(ro_bar/ro_star)*M_t/(vis_bar/vis_star)/sqrt(3d0)    !(both these expressions are equivalent)
!     ! Re_lambda = lambda*sqrt(u_p_sq_bar)*ro_bar/vis_bar

!     if (iter==0) then
!         t_turb = lambda/sqrt(turb_KE/3d0)
!         write(*,*)'turbulent time scale (non-dim)',t_turb/(lo_star/co_star)
!         write(*,*)'Re_lambda',Re_lambda
!     end if

!     CALL MPI_ALLREDUCE(duc_count,WENO,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mperr)

!     if(myrank==0) then
!         open(unit=40,file='stats.dat',STATUS='UNKNOWN',ACTION='WRITE',POSITION='APPEND',IOSTAT=ierr)
!         write(40,'(1x,28E30.18)')iter*dt/t_turb,lambda_x/lo_star,lambda_y/lo_star,lambda_z/lo_star,lambda/lo_star,skew_x,skew_y,skew_z,skew,Re_lambda,turb_KE/intial_turb_KE,M_t,sqrt(u_p_sq_bar)/co_star &
!                                  ,sqrt(v_p_sq_bar)/co_star,sqrt(w_p_sq_bar)/co_star,sqrt(ro_p_sq_bar)/ro_star,sqrt(p_p_sq_bar)/p_star,sqrt(T_p_sq_bar)/T_star,u_bar/co_star,v_bar/co_star,w_bar/co_star, &
!                                  ro_bar/ro_star,ro1_bar/ro_star,ro2_bar/ro_star,cs1_bar,cs2_bar,p_bar/p_star,T_bar/T_star
!         close(40)
!     end if

! end subroutine stats

!******************************************************************************************************************************************************************************************************





!***********************************************************************************CALCULATES EIGENVECTORS********************************************************************************************

real*8 function global_average(var,power)

use variables
implicit none
include 'mpif.h'
real*8, dimension(sx,sy,sz), intent(in) :: var
integer, intent(in)                     :: power
integer                                 :: i,j,k
real*8                                  :: sum

sum = 0d0

do k=startz,endz
    do j=starty,endy
        do i=startx,endx
            sum = sum + var(i,j,k)**power 
        end do
    end do
end do

CALL MPI_ALLREDUCE(sum,global_average,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mperr)
global_average = global_average/(sizex*sizey*sizez*1d0)

return

end function global_average

!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************










!******************************************************************************************************************************************************************************************************
!*********************************************************************************THE BOUNDARY CONDITIONS**********************************************************************************************
!******************************************************************************************************************************************************************************************************
!SEE message_pass FOR q_side and q_top

subroutine bc(rkl)
  
    use variables
    implicit none
    include 'mpif.h'
    integer                                      :: i,j,k,m,tag1,tag2,tag3,tag4
    integer                                      :: status(mpi_status_size),q_top,q_side
    integer,intent(in)                           :: rkl
    double precision                             :: fextrapolate,bextrapolate,vel_square,gamma,add,R_s,vib_temp
    double precision, dimension(sx,sy,sz,1:n_sp) :: cs
    double precision, dimension(1:n_sp)          :: ro_s,one
    double precision, dimension(1:n_di)          :: ev_s

    call MPI_TYPE_VECTOR(sz*neq,3*sx,sx*sy,MPI_REAL8,q_top,mperr)
    call MPI_TYPE_COMMIT(q_top,mperr)
    call MPI_TYPE_VECTOR(neq,sx*sy*3,sx*sy*sz,MPI_REAL8,q_side,mperr)
    call MPI_TYPE_COMMIT(q_side,mperr)

    tag1=10
    tag2=20
    tag3=30
    tag4=40

    
    !PERIODIC BOUNDARY CONDITION IN Y DIRECTION
    if (periodic_y==1) then
        if (myrank/zprocs == 0) then
            call MPI_SEND(q(1,4,1,1),1,q_top,myrank+zprocs*(yprocs-1),tag3,MPI_COMM_WORLD,mperr)
        end if
        if (myrank/zprocs == yprocs-1) then
            call MPI_RECV(q(1,sy-2,1,1),1,q_top,myrank-zprocs*(yprocs-1),tag3,MPI_COMM_WORLD,status,mperr)
            call MPI_SEND(q(1,sy-5,1,1),1,q_top,myrank-zprocs*(yprocs-1),tag4,MPI_COMM_WORLD,mperr)
        end if
        if (myrank/zprocs == 0) then
            call MPI_RECV(q(1,1,1,1),1,q_top,myrank+zprocs*(yprocs-1),tag4,MPI_COMM_WORLD,status,mperr)
        end if
    end if
    
    !PERIODIC BOUNDARY CONDITION IN Z DIRECTION
    if (periodic_z==1) then
        if (mod(myrank,zprocs) == 0) then
            call MPI_SEND(q(1,1,4,1),1,q_side,myrank+zprocs-1,tag1,MPI_COMM_WORLD,mperr)
        end if
        if (mod(myrank,zprocs) == zprocs-1) then
            call MPI_RECV(q(1,1,sz-2,1),1,q_side,myrank-zprocs+1,tag1,MPI_COMM_WORLD,status,mperr)
            call MPI_SEND(q(1,1,sz-5,1),1,q_side,myrank-zprocs+1,tag2,MPI_COMM_WORLD,mperr)
        end if
        if (mod(myrank,zprocs) == 0) then
            call MPI_RECV(q(1,1,1,1),1,q_side,myrank+zprocs-1,tag2,MPI_COMM_WORLD,status,mperr)
        end if
    end if
    
    !PERIODIC BOUNDARY CONDITION IN X DIRECTION
    if (periodic_x==1) then
        q(sx-2,:,:,:) = q(4,:,:,:)
        q(sx-1,:,:,:) = q(5,:,:,:)
        q(sx,:,:,:) = q(6,:,:,:)
        q(1,:,:,:) = q(sx-5,:,:,:)
        q(2,:,:,:) = q(sx-4,:,:,:)
        q(3,:,:,:) = q(sx-3,:,:,:)
    end if

end subroutine bc

!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************










!******************************************************************************************************************************************************************************************************
!*****************************************************************************************THE MESH*****************************************************************************************************
!******************************************************************************************************************************************************************************************************
!x,y,z are the computational computational cordinates, x_ph,y_ph,z_ph is the physical coordinates
!d()_d() are the ttransformation metrics

subroutine gridin()   

    use variables
    implicit none
    integer      :: i,j,k,ierr,numrc1,numrc2,numrc3,numrc4
    double precision       :: rmin,rmax,radius,thet,delta
    character*40 :: filename

    do k=1,sz
        do j=1,sy 
            do i=1,sx
                x(i,j,k) = (i-4)*dx
                z(i,j,k) = mod(myrank,zprocs)*dz*(sz-6)+(k-4)*dz
                y(i,j,k) = (myrank/zprocs)*dy*(sy-6)+(j-4)*dy
                y_ph(i,j,k) = leny*((cl_para + 1d0) - (cl_para - 1d0)*((cl_para + 1d0)/(cl_para - 1d0))**(1d0 - y(i,j,k)))/(((cl_para + 1d0)/(cl_para - 1d0))**(1d0 - y(i,j,k)) + 1d0)
                dy_dy_ph(i,j,k) = 2d0*cl_para/(leny*(cl_para**2d0 - (1d0 - y_ph(i,j,k)/leny)**2)*log((cl_para + 1d0)/(cl_para - 1d0)))
            end do
        end do
    end do

end subroutine gridin

!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************










!******************************************************************************************************************************************************************************************************
!**********************************************************************COMPUTES ALL THE CONVECTIVE AND VISCOUS FLUXES**********************************************************************************
!******************************************************************************************************************************************************************************************************

subroutine update(rkl)

    use variables
    implicit none
    include 'mpif.h'
    integer,intent(in)                                :: rkl
    integer                                           :: i,j,k,m,n,l
    double precision                                  :: energy,i_engy,div_duc,omg_duc
    double precision                                  :: towxx,towyy,towzz,towxy,towxz,towyz,towyx,towzx,towzy
    double precision, allocatable, dimension(:,:,:,:) :: flux_iph,flux_jph,flux_kph
    double precision,dimension(sx,sy,sz,n_sp)         :: cs 
    double precision,dimension(sx,sy,sz)              :: local_mach,tot_enthalpy,ducros2,a
    double precision, dimension(neq)                  :: eigen
    double precision, dimension(n_sp)                 :: mu,kap,X_t,phi_t,dum1,dum2,f_diff,g_diff,h_diff,hs,ev_s,one
    double precision, dimension(n_di)                 :: kap_v,dum_v
    double precision                                  :: d6dx,d6dy,d6dz,d6dx_vect,d6dy_vect,d6dz_vect,e1,e2,vel_square,add,vib_temp
    double precision                                  :: D,R_s,Mw_mix !molecular weight of the mixture    

    one(:) = 1d0

    do k=firstz,lastz
        do j=firsty,lasty
            do i=firstx,lastx
                ro(i,j,k) = add(q(i,j,k,1:n_sp),one,n_sp,1)
                cs(i,j,k,1:n_sp) = q(i,j,k,1:n_sp)/ro(i,j,k)    
                u(i,j,k) = q(i,j,k,n_sp+1)/ro(i,j,k)
                v(i,j,k) = q(i,j,k,n_sp+2)/ro(i,j,k)
                w(i,j,k) = q(i,j,k,n_sp+3)/ro(i,j,k)
                if (ISNAN(u(i,j,k)) == .TRUE.) then
                    ! if(myrank==0) write(*,*) iter,rkl,i,j,k!,ro(i,j,k),p(i,j,k),T(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k)
                    write(*,*) iter,rkl,i,j,k!,ro(i,j,k),p(i,j,k),T(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k)
                    STOP
                end if
                vel_square = u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)
                if (vib_flag==1) T(i,j,k) = (q(i,j,k,n_sp+4) - q(i,j,k,n_sp+5) - 0.50d0*ro(i,j,k)*vel_square - chem_flag*add(q(i,j,k,1:n_sp),h0,n_sp,1))/add(q(i,j,k,1:n_sp),Cv,n_sp,1)
                if (vib_flag==0) T(i,j,k) = (q(i,j,k,n_sp+4) - 0.50d0*ro(i,j,k)*vel_square - chem_flag*add(q(i,j,k,1:n_sp),h0,n_sp,1))/add(q(i,j,k,1:n_sp),Cv,n_sp,1)
                R_s = R_u*add(q(i,j,k,1:n_sp),Mw,n_sp,-1)/ro(i,j,k)
                p(i,j,k) = ro(i,j,k)*R_s*T(i,j,k)
                a(i,j,k) = sqrt(R_s*T(i,j,k)*(1 + R_s*ro(i,j,k)/add(q(i,j,k,1:n_sp),Cv,n_sp,1)))
                local_mach(i,j,k) = dsqrt(vel_square)/a(i,j,k)
                tot_enthalpy(i,j,k) = p(i,j,k) + q(i,j,k,n_sp+4)
                if (vib_flag==1) Tv(i,j,k) = vib_temp(i,j,k)
                if (vib_flag==0) Tv(i,j,k) = T(i,j,k)

                ! if (vis_flag == 1 .and. n_sp==1) then
                    ! vis(i,j,k) = vis_ref*(T(i,j,k)/T_ref)**(3d0/2d0)*((T_ref+110d0)/(T(i,j,k)+110d0))
                    ! cond(i,j,k) = vis(i,j,k)*R_star*gamma_star/(gamma_star - 1.0d0)/prandtl  
                    ! vis(i,j,k) = vis_ref*(T(i,j,k)/T_ref)**0.670d0
                    ! cond(i,j,k) = vis(i,j,k)*R_star*gamma_star/(gamma_star - 1.0d0)/prandtl        
                    ! cond_v(i,j,k) = 0d0
                ! elseif (vis_flag == 1) then
                    mu(:) = 0.1d0*exp((At(:)*log(T(i,j,k)) + Bt(:))*log(T(i,j,k)) + Ct(:))
                    kap(:) = mu(:)*(Cv(:) + 3d0*Cv_tr(:)/2d0)
                    kap_v(:) = vib_flag*(mu(1:n_di)*(R_u/Mw(1:n_di))*(theta_v(:)/Tv(i,j,k))**2*(exp(theta_v(:)/Tv(i,j,k))/(exp(theta_v(:)/Tv(i,j,k))-1)**2))
                    Mw_mix = add(q(i,j,k,1:n_sp),Mw,n_sp,-1)/ro(i,j,k)
                    X_t(:) = Mw_mix*q(i,j,k,1:n_sp)/(Mw(:)*ro(i,j,k))
                    do m=1,n_sp
                        dum1(:) = X_t(:)*(1 + (mu(m)/mu(:))**(1/2d0)*(Mw(:)/Mw(m))**(1/4d0))**2
                        dum2(:) = sqrt(8*(1 + Mw(m)/Mw(:)))
                        phi_t(m) = add(dum1,dum2,n_sp,-1)
                    end do
                    dum1(:) = X_t(:)*mu(:)
                    dum2(:) = X_t(:)*kap(:)
                    dum_v(:) = X_t(1:n_di)*kap_v(:)
                    vis(i,j,k) = add(dum1,phi_t,n_sp,-1)
                    cond(i,j,k) = add(dum2,phi_t,n_sp,-1)
                    cond_v(i,j,k) = add(dum_v,phi_t(1:n_di),n_di,-1)
                ! elseif (vis_flag==0) then
                !     vis(i,j,k) = 0d0
                !      cond(i,j,k) = 0d0
                !     cond_v(i,j,k) = 0d0
                ! end if
            end do
        end do
    end do

    

    ! duc_count = 0
    ! ducros(:,:,:) = 0.0d0
    ! ducros2(:,:,:) = 0d0
    ! e1 = 0.01d0
    ! e2 = 0.5d0
    ! do k=startz,endz
    !     do j=starty,endy
    !         do i=startx,endx               
    !             if ((local_mach(i,j,k)<1d0 .and. (((ro(i,j,k)-ro(i-1,j,k))/ro_star)>e2 .or. ((ro(i,j,k)-ro(i,j-1,k))/ro_star)>e2 .or. ((ro(i,j,k)-ro(i,j,k-1))/ro_star)>e2)) .or. &
    !                 (local_mach(i,j,k)>=1d0 .and. (((ro(i,j,k)-ro(i-1,j,k))/ro_star)>e1 .or. ((ro(i,j,k)-ro(i,j-1,k))/ro_star)>e1 .or. ((ro(i,j,k)-ro(i,j,k-1))/ro_star)>e1))  ) then              
    !                 duc_count = duc_count+1
    !                 do n=k-3,k+3
    !                     ducros(i,j,n)=0.6
    !                 end do
    !                 do n=j-3,j+3
    !                     ducros(i,n,k)=0.6
    !                 end do
    !                 do n=i-3,i+3
    !                     ducros(n,j,k)=0.6
    !                 end do
    !             end if
    !         end do
    !     end do
    ! end do
    ducros(:,:,:) = 0d0!.6

    do k=startz,endz
        do j=starty,endy
            do i=startx,endx
                if (ducros(i,j,k)>0.5d0) then
                    ducros2(i,j,k) = 0.6
                    ducros2(i-1,j,k) = 0.6
                    ducros2(i,j-1,k) = 0.6
                    ducros2(i,j,k-1) = 0.6
                end if
            end do
        end do
    end do

    
    
    do k=firstz,lastz
        do j=firsty,lasty
            do i=firstx,lastx
                fluxf(i,j,k,:)  = q(i,j,k,:)*u(i,j,k)
                fluxf(i,j,k,n_sp+4)  = fluxf(i,j,k,n_sp+4) + p(i,j,k)*u(i,j,k)
                
                fluxg(i,j,k,:)  = q(i,j,k,:)*v(i,j,k)
                fluxg(i,j,k,n_sp+4)  = fluxg(i,j,k,n_sp+4) + p(i,j,k)*v(i,j,k)

                fluxh(i,j,k,:)  = q(i,j,k,:)*w(i,j,k)
                fluxh(i,j,k,n_sp+4)  = fluxh(i,j,k,n_sp+4) + p(i,j,k)*w(i,j,k)
            end do
        end do
    end do

    do k=startz_vis,endz_vis
        do j=starty_vis,endy_vis
            do i=startx_vis,endx_vis
                if(n_sp/=1) then
                    D = Lewis*cond(i,j,k)/add(q(i,j,k,1:n_sp),Cp,n_sp,1)  !the formula says conductivity but doesnt specify whther translational or rotational or what. so tr+rot cond is taken here.
                    do m=1,n_sp
                        f_diff(m) = -ro(i,j,k)*D*d6dx(cs(:,:,:,m),i,j,k)                
                        g_diff(m) = -ro(i,j,k)*D*d6dy(cs(:,:,:,m),i,j,k)
                        h_diff(m) = -ro(i,j,k)*D*d6dz(cs(:,:,:,m),i,j,k)
                    end do
                else
                    f_diff(1:n_sp) = 0d0
                    g_diff(1:n_sp) = 0d0
                    h_diff(1:n_sp) = 0d0
                end if

                ev_s(1:n_di) = vib_flag*((R_u/Mw(1:n_di))*theta_v(:)/(exp(theta_v(:)/Tv(i,j,k))-1))
                ev_s(n_di+1:n_sp) = 0d0

                hs(:) = Cv(:)*T(i,j,k) + R_u*T(i,j,k)/Mw(:) + vib_flag*ev_s(:) + chem_flag*h0(:)

                towxx = vis(i,j,k)*(4*d6dx(u,i,j,k)-2*(d6dy(v,i,j,k)+d6dz(w,i,j,k)))/3
                towyy = vis(i,j,k)*(4*d6dy(v,i,j,k)-2*(d6dx(u,i,j,k)+d6dz(w,i,j,k)))/3
                towzz = vis(i,j,k)*(4*d6dz(w,i,j,k)-2*(d6dy(v,i,j,k)+d6dx(u,i,j,k)))/3
                towxy = vis(i,j,k)*(d6dy(u,i,j,k)+d6dx(v,i,j,k))
                towyz = vis(i,j,k)*(d6dy(w,i,j,k)+d6dz(v,i,j,k))
                towzx = vis(i,j,k)*(d6dz(u,i,j,k)+d6dx(w,i,j,k))

                heatx(i,j,k) = -cond(i,j,k)*d6dx(T,i,j,k) - cond_v(i,j,k)*d6dx(Tv,i,j,k)*vib_flag
                heaty(i,j,k) = -cond(i,j,k)*d6dy(T,i,j,k) - cond_v(i,j,k)*d6dy(Tv,i,j,k)*vib_flag
                heatz(i,j,k) = -cond(i,j,k)*d6dz(T,i,j,k) - cond_v(i,j,k)*d6dz(Tv,i,j,k)*vib_flag

                visf(i,j,k,1:n_sp)  = -f_diff(1:n_sp)
                visf(i,j,k,n_sp+1)  = towxx
                visf(i,j,k,n_sp+2)  = towxy
                visf(i,j,k,n_sp+3)  = towzx
                visf(i,j,k,n_sp+4)  = u(i,j,k)*towxx + v(i,j,k)*towxy + w(i,j,k)*towzx - heatx(i,j,k) - add(f_diff,hs,n_sp,1)
                if (vib_flag==1) visf(i,j,k,n_sp+5) = -(-cond_v(i,j,k)*d6dx(Tv,i,j,k)) - add(f_diff,ev_s,n_di,1)
                
                visg(i,j,k,1:n_sp)  = -g_diff(1:n_sp)
                visg(i,j,k,n_sp+1)  = towxy
                visg(i,j,k,n_sp+2)  = towyy
                visg(i,j,k,n_sp+3)  = towyz
                visg(i,j,k,n_sp+4)  = u(i,j,k)*towxy + v(i,j,k)*towyy + w(i,j,k)*towyz - heaty(i,j,k) - add(g_diff,hs,n_sp,1)
                if (vib_flag==1) visg(i,j,k,n_sp+5) = -(-cond_v(i,j,k)*d6dx(Tv,i,j,k)) - add(g_diff,ev_s,n_di,1)
                
                vish(i,j,k,1:n_sp)  = -h_diff(1:n_sp) 
                vish(i,j,k,n_sp+1)  = towzx
                vish(i,j,k,n_sp+2)  = towyz
                vish(i,j,k,n_sp+3)  = towzz
                vish(i,j,k,n_sp+4)  = u(i,j,k)*towzx + v(i,j,k)*towyz + w(i,j,k)*towzz - heatz(i,j,k) - add(h_diff,hs,n_sp,1)
                if (vib_flag==1) vish(i,j,k,n_sp+5) = -(-cond_v(i,j,k)*d6dx(Tv,i,j,k)) - add(h_diff,ev_s,n_di,1)
            end do
        end do
    end do

    
    alpha_x_global(1,:,:,:) = 0d0
    do k=firstz,lastz
        do j=firsty,lasty
            do i=firstx,lastx
                call eigen_val(u(i,j,k),a(i,j,k),neq,eigen)
                do l=1,neq
                    if(abs(eigen(l))>alpha_x_global(1,j,k,l)) alpha_x_global(1,j,k,l) = abs(eigen(l))
                end do
            end do
        end do
    end do

    alpha_y(:,:,:,:) = 0d0
    alpha_y_global(:,:,:,:) = 0d0
    do i=firstx,lastx
        do k=firstz,lastz
            do j=firsty,lasty
                call eigen_val(v(i,j,k),a(i,j,k),neq,eigen)
                do l=1,neq
                    if(abs(eigen(l))>alpha_y(i,mod(myrank,zprocs),k,l)) alpha_y(i,mod(myrank,zprocs),k,l) = abs(eigen(l))
                end do
            end do
        end do
    end do
    do i=firstx,lastx
        do k=firstz,lastz
            do l=1,neq
                do m=0,zprocs-1
                    CALL MPI_ALLREDUCE(alpha_y(i,m,k,l),alpha_y_global(i,m,k,l),1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mperr)
                end do
            end do
        end do
    end do

    alpha_z(:,:,:,:) = 0d0
    alpha_z_global(:,:,:,:) = 0d0
    do j=firsty,lasty
        do i=firstx,lastx
            do k=firstz,lastz
                call eigen_val(w(i,j,k),a(i,j,k),neq,eigen)
                do l=1,neq
                    if(abs(eigen(l))>alpha_z(i,j,myrank/zprocs,l)) alpha_z(i,j,myrank/zprocs,l) = abs(eigen(l))
                end do
            end do
        end do
    end do
    do j=firsty,lasty
        do i=firstx,lastx
            do l=1,neq
                do m=0,yprocs-1
                    CALL MPI_ALLREDUCE(alpha_z(i,j,m,l),alpha_z_global(i,j,m,l),1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mperr)
                end do
            end do
        end do
    end do

    
    allocate(flux_iph(sx,sy,sz,neq),flux_jph(sx,sy,sz,neq),flux_kph(sx,sy,sz,neq))
    do k=startz-1,endz
        do j=starty-1,endy
            do i=startx-1,endx
                if (ducros2(i,j,k)>0.5d0) then
                    call wenox(i,j,k,flux_iph(i,j,k,:)) 
                    call wenoy(i,j,k,flux_jph(i,j,k,:)) 
                    call wenoz(i,j,k,flux_kph(i,j,k,:))
                end if
            end do
        end do
    end do

    do k=startz,endz
        do j=starty,endy
            do i=startx,endx
                if (ducros(i,j,k)>0.5d0) then
                    ddx_fluxf(i,j,k,:) = (flux_iph(i,j,k,:)-flux_iph(i-1,j,k,:))/dx
                    ddy_fluxg(i,j,k,:) = dy_dy_ph(i,j,k)*(flux_jph(i,j,k,:)-flux_jph(i,j-1,k,:))/dy
                    ddz_fluxh(i,j,k,:) = (flux_kph(i,j,k,:)-flux_kph(i,j,k-1,:))/dz
                else                
                    do m=1,n_sp
                        ddx_fluxf(i,j,k,m)  = 0.50d0*d6dx_vect(fluxf,i,j,k,m) + 0.50d0*u(i,j,k)*d6dx_vect(q,i,j,k,m) + 0.50d0*q(i,j,k,m)*d6dx(u,i,j,k)
                    end do
                    ddx_fluxf(i,j,k,n_sp+1)  = 0.50d0*d6dx_vect(fluxf,i,j,k,n_sp+1) + 0.50d0*u(i,j,k)*d6dx_vect(q,i,j,k,n_sp+1) + 0.50d0*q(i,j,k,n_sp+1)*d6dx(u,i,j,k) + d6dx(p,i,j,k)
                    ddx_fluxf(i,j,k,n_sp+2)  = 0.50d0*d6dx_vect(fluxf,i,j,k,n_sp+2) + 0.50d0*u(i,j,k)*d6dx_vect(q,i,j,k,n_sp+2) + 0.50d0*q(i,j,k,n_sp+2)*d6dx(u,i,j,k)
                    ddx_fluxf(i,j,k,n_sp+3)  = 0.50d0*d6dx_vect(fluxf,i,j,k,n_sp+3) + 0.50d0*u(i,j,k)*d6dx_vect(q,i,j,k,n_sp+3) + 0.50d0*q(i,j,k,n_sp+3)*d6dx(u,i,j,k)
                    ddx_fluxf(i,j,k,n_sp+4)  = 0.50d0*d6dx_vect(fluxf,i,j,k,n_sp+4) + 0.50d0*u(i,j,k)*d6dx(tot_enthalpy,i,j,k) + 0.50d0*tot_enthalpy(i,j,k)*d6dx(u,i,j,k)
                    if(vib_flag==1) ddx_fluxf(i,j,k,n_sp+5)  = 0.50d0*d6dx_vect(fluxf,i,j,k,n_sp+5) + 0.50d0*u(i,j,k)*d6dx_vect(q,i,j,k,n_sp+5) + 0.50d0*q(i,j,k,n_sp+5)*d6dx(u,i,j,k)

                    do m=1,n_sp
                        ddy_fluxg(i,j,k,m)  = 0.50d0*d6dy_vect(fluxg,i,j,k,m) + 0.50d0*v(i,j,k)*d6dy_vect(q,i,j,k,m) + 0.50d0*q(i,j,k,m)*d6dy(v,i,j,k)
                    end do
                    ddy_fluxg(i,j,k,n_sp+1)  = 0.50d0*d6dy_vect(fluxg,i,j,k,n_sp+1) + 0.50d0*v(i,j,k)*d6dy_vect(q,i,j,k,n_sp+1) + 0.50d0*q(i,j,k,n_sp+1)*d6dy(v,i,j,k)
                    ddy_fluxg(i,j,k,n_sp+2)  = 0.50d0*d6dy_vect(fluxg,i,j,k,n_sp+2) + 0.50d0*v(i,j,k)*d6dy_vect(q,i,j,k,n_sp+2) + 0.50d0*q(i,j,k,n_sp+2)*d6dy(v,i,j,k) + d6dy(p,i,j,k)
                    ddy_fluxg(i,j,k,n_sp+3)  = 0.50d0*d6dy_vect(fluxg,i,j,k,n_sp+3) + 0.50d0*v(i,j,k)*d6dy_vect(q,i,j,k,n_sp+3) + 0.50d0*q(i,j,k,n_sp+3)*d6dy(v,i,j,k)
                    ddy_fluxg(i,j,k,n_sp+4)  = 0.50d0*d6dy_vect(fluxg,i,j,k,n_sp+4) + 0.50d0*v(i,j,k)*d6dy(tot_enthalpy,i,j,k) + 0.50d0*tot_enthalpy(i,j,k)*d6dy(v,i,j,k)
                    if(vib_flag==1) ddy_fluxg(i,j,k,n_sp+5)  = 0.50d0*d6dy_vect(fluxg,i,j,k,n_sp+5) + 0.50d0*v(i,j,k)*d6dy_vect(q,i,j,k,n_sp+5) + 0.50d0*q(i,j,k,n_sp+5)*d6dy(v,i,j,k)

                    do m=1,n_sp
                        ddz_fluxh(i,j,k,m)  = 0.50d0*d6dz_vect(fluxh,i,j,k,m) + 0.50d0*w(i,j,k)*d6dz_vect(q,i,j,k,m) + 0.50d0*q(i,j,k,m)*d6dz(w,i,j,k)
                    end do
                    ddz_fluxh(i,j,k,n_sp+1)  = 0.50d0*d6dz_vect(fluxh,i,j,k,n_sp+1) + 0.50d0*w(i,j,k)*d6dz_vect(q,i,j,k,n_sp+1) + 0.50d0*q(i,j,k,n_sp+1)*d6dz(w,i,j,k)
                    ddz_fluxh(i,j,k,n_sp+2)  = 0.50d0*d6dz_vect(fluxh,i,j,k,n_sp+2) + 0.50d0*w(i,j,k)*d6dz_vect(q,i,j,k,n_sp+2) + 0.50d0*q(i,j,k,n_sp+2)*d6dz(w,i,j,k)
                    ddz_fluxh(i,j,k,n_sp+3)  = 0.50d0*d6dz_vect(fluxh,i,j,k,n_sp+3) + 0.50d0*w(i,j,k)*d6dz_vect(q,i,j,k,n_sp+3) + 0.50d0*q(i,j,k,n_sp+3)*d6dz(w,i,j,k) + d6dz(p,i,j,k)
                    ddz_fluxh(i,j,k,n_sp+4)  = 0.50d0*d6dz_vect(fluxh,i,j,k,n_sp+4) + 0.50d0*w(i,j,k)*d6dz(tot_enthalpy,i,j,k) + 0.50d0*tot_enthalpy(i,j,k)*d6dz(w,i,j,k)
                    if(vib_flag==1) ddz_fluxh(i,j,k,n_sp+5)  = 0.50d0*d6dz_vect(fluxh,i,j,k,n_sp+5) + 0.50d0*w(i,j,k)*d6dz_vect(q,i,j,k,n_sp+5) + 0.50d0*q(i,j,k,n_sp+5)*d6dz(w,i,j,k)
                end if
            end do
        end do
    end do
    deallocate(flux_iph,flux_jph,flux_kph)

end subroutine update



!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************










!******************************************************************************************************************************************************************************************************
!******************************************************SUBROUTINES AND FUNCTIONS NEEDED FOR CALCULATING CONVECTIVE FLUX DERIVATIVES USING WENO*********************************************************
!******************************************************************************************************************************************************************************************************
!REFER ICASE REPORT ON ENO AND WENO BY CHI-WANG SHU
!CONSTRUCTS FLUX AT i+1/2 USING FIFTH ORDER WENO. LOCAL LAX-FRIEDRICH FLUX SPLITTING IS USED.
!THE VARIABLE NAMES ARE SIMILAR TO WHAT USED IN THE REFERENCE
!_l=LEFT, _r=RIGHT, _pr=PRIMITIVE, et=TOTAL ENERGY, ht=TOTAL ENTHALPY


!******************************************************************************FLUX CONSTRUCTION IN X DIRECTION****************************************************************************************

subroutine wenox(i,j,k,f_ih)
    
    use variables
    implicit none
    integer, intent(in)                          :: i,j,k
    double precision, dimension(neq),intent(out) :: f_ih
    double precision, dimension(neq)             :: plus0,plus1,plus2,minus0,minus1,minus2,beta0,beta1,beta2,beta0_dash,beta1_dash,beta2_dash,temp                               
    double precision, dimension(neq)             :: alpha0,alpha1,alpha2,alpha0_dash,alpha1_dash,alpha2_dash
    double precision, dimension(neq)             :: wt0,wt1,wt2,wt0_dash,wt1_dash,wt2_dash
    double precision, dimension(neq,neq)         :: right_eigen,left_eigen
    double precision, dimension(neq)             :: flux_minus,flux_plus
    double precision, dimension(n_sp)            :: ro_r, ro_l, ro_avg
    double precision                             :: f_charac(i-2:i+3,neq),f_charac_plus(i-2:i+3,neq),f_charac_minus(i-2:i+3,neq)
    double precision                             :: q_charac_x(i-2:i+3,neq)
    double precision                             :: epsilon=1e-6, d0=3d0/10d0, d1=3d0/5d0, d2=1d0/10d0, d0_dash=1d0/10d0, d1_dash=3d0/5d0, d2_dash=3d0/10d0
    double precision                             :: u_r, u_l, v_r, v_l, w_r, w_l, ht_r, ht_l, a_r, a_l, u_avg,Tv_r,Tv_l, v_avg, w_avg, ht_avg, a_avg, p_pr, et_pr,Tv_avg,gamma
    integer                                      :: l, ll, entropy

    
    call primitive(i,j,k,ro_l,u_l,v_l,w_l,ht_l,a_l,Tv_l)
    call primitive(i+1,j,k,ro_r,u_r,v_r,w_r,ht_r,a_r,Tv_r)
    call averages(ro_r,ro_l,u_r,u_l,v_r,v_l,w_r,w_l,ht_r,ht_l,Tv_r,Tv_l,ro_avg,u_avg,v_avg,w_avg,ht_avg,a_avg,Tv_avg,gamma)
    call eigenvector(1,0,n_sp,neq,ro_avg(1:n_sp),u_avg,v_avg,w_avg,a_avg,ht_avg,h0(1:n_sp),(gamma-1d0),right_eigen)
    call eigenvector(1,1,n_sp,neq,ro_avg(1:n_sp),u_avg,v_avg,w_avg,a_avg,ht_avg,h0(1:n_sp),(gamma-1d0),left_eigen)
    

    if (i==firstx) then
        do l=i,i+3,1
            temp(:) = fluxf(l,j,k,:)
            temp(n_sp+1) = fluxf(l,j,k,n_sp+1)+p(l,j,k)
            f_charac(l,:) = MATMUL(left_eigen,temp)
        end do
        do l=i-2,i-1,1
            temp(:) = (fluxf(firstx,j,k,:))*(10d0*abs(firstx-l))**12d0
            temp(n_sp+1) = (fluxf(firstx,j,k,n_sp+1)+p(firstx,j,k))*(10d0*abs(firstx-l))**12d0
            f_charac(l,:) = MATMUL(left_eigen,temp)
        end do
    else if(i==firstx+1) then
        do l=i-1,i+3,1
            temp(:) = fluxf(l,j,k,:)
            temp(n_sp+1) = fluxf(l,j,k,n_sp+1)+p(l,j,k)
            f_charac(l,:) = MATMUL(left_eigen,temp)
        end do
        l=i-2
            temp(:) = (fluxf(firstx,j,k,:))*(10d0*(abs(firstx-l)-1))**12d0
            temp(n_sp+1) = (fluxf(firstx,j,k,n_sp+1)+p(firstx,j,k))*(10d0*(abs(firstx-l)-1))**12d0
            f_charac(l,:) = MATMUL(left_eigen,temp)
    else if(i==lastx) then
        do l=i-2,i,1
            temp(:) = fluxf(l,j,k,:)
            temp(n_sp+1) = fluxf(l,j,k,n_sp+1)+p(l,j,k)
            f_charac(l,:) = MATMUL(left_eigen,temp)
        end do
        do l=i+1,i+3,1
            temp(:) = (fluxf(lastx,j,k,:))*(10d0*abs(lastx-l))**12d0
            temp(n_sp+1) = (fluxf(lastx,j,k,n_sp+1)+p(lastx,j,k))*(10d0*abs(lastx-l))**12d0
            f_charac(l,:) = MATMUL(left_eigen,temp)
        end do
    else if(i==lastx-1) then
        do l=i-2,i+1,1
            temp(:) = fluxf(l,j,k,:)
            temp(n_sp+1) = fluxf(l,j,k,n_sp+1)+p(l,j,k)
            f_charac(l,:) = MATMUL(left_eigen,temp)
        end do
        do l=i+2,i+3,1
            temp(:) = (fluxf(lastx,j,k,:))*(10d0*(abs(lastx-l)-1))**12d0
            temp(n_sp+1) = (fluxf(lastx,j,k,n_sp+1)+p(lastx,j,k))*(10d0*(abs(lastx-l)-1))**12d0
            f_charac(l,:) = MATMUL(left_eigen,temp)
        end do
    else if(i==lastx-2) then
        do l=i-2,i+2,1
            temp(:) = fluxf(l,j,k,:)
            temp(n_sp+1) = fluxf(l,j,k,n_sp+1)+p(l,j,k)
            f_charac(l,:) = MATMUL(left_eigen,temp)
        end do
        l=i+3
            temp(:) = (fluxf(lastx,j,k,:))*(10d0*(abs(lastx-l)-2))**12d0
            temp(n_sp+1) = (fluxf(lastx,j,k,n_sp+1)+p(lastx,j,k))*(10d0*(abs(lastx-l)-2))**12d0
            f_charac(l,:) = MATMUL(left_eigen,temp)
    else
        do l=i-2,i+3,1
            temp(:) = fluxf(l,j,k,:)
            temp(n_sp+1) = fluxf(l,j,k,n_sp+1)+p(l,j,k)
            f_charac(l,:) = MATMUL(left_eigen,temp)
        end do
    end if


    do l=i-2,i+3,1
        if(l<firstx) then
            q_charac_x(l,:)  = MATMUL(left_eigen,q(firstx,j,k,:))
        else if(l>lastx) then
            q_charac_x(l,:)  = MATMUL(left_eigen,q(lastx,j,k,:))
        else
            q_charac_x(l,:)  = MATMUL(left_eigen,q(l,j,k,:))
        end if
    end do


    do l=i-2,i+3,1
        f_charac_plus(l,:) = 0.5d0*(f_charac(l,:) + alpha_x_global(1,j,k,:)*q_charac_x(l,:))
        f_charac_minus(l,:) = 0.5d0*(f_charac(l,:) - alpha_x_global(1,j,k,:)*q_charac_x(l,:))
    end do


    minus0(:) = f_charac_plus(i,:)/3d0 + 5d0*f_charac_plus(i+1,:)/6d0 - f_charac_plus(i+2,:)/6d0
    minus1(:) = -f_charac_plus(i-1,:)/6d0 + 5d0*f_charac_plus(i,:)/6d0 + f_charac_plus(i+1,:)/3d0
    minus2(:) = f_charac_plus(i-2,:)/3d0 - 7d0*f_charac_plus(i-1,:)/6d0 + 11d0*f_charac_plus(i,:)/6d0

    plus0(:) = 11d0*f_charac_minus(i+1,:)/6d0 - 7d0*f_charac_minus(i+2,:)/6d0 + f_charac_minus(i+3,:)/3d0
    plus1(:) = f_charac_minus(i,:)/3d0 + 5d0*f_charac_minus(i+1,:)/6d0 - f_charac_minus(i+2,:)/6d0
    plus2(:) = -f_charac_minus(i-1,:)/6d0 + 5d0*f_charac_minus(i,:)/6d0 + f_charac_minus(i+1,:)/3d0

    beta0(:) = 13d0*(f_charac_plus(i,:) - 2d0*f_charac_plus(i+1,:) + f_charac_plus(i+2,:))**2/12d0 + (3d0*f_charac_plus(i,:) - 4d0*f_charac_plus(i+1,:) + f_charac_plus(i+2,:))**2/4d0
    beta1(:) = 13d0*(f_charac_plus(i-1,:) - 2d0*f_charac_plus(i,:) + f_charac_plus(i+1,:))**2/12d0 + (f_charac_plus(i-1,:) - f_charac_plus(i+1,:))**2/4d0
    beta2(:) = 13d0*(f_charac_plus(i-2,:) - 2d0*f_charac_plus(i-1,:) + f_charac_plus(i,:))**2/12d0 + (f_charac_plus(i-2,:) - 4d0*f_charac_plus(i-1,:) + 3d0*f_charac_plus(i,:))**2/4d0

    beta0_dash(:) = 13d0*(f_charac_minus(i+1,:) - 2d0*f_charac_minus(i+2,:) + f_charac_minus(i+3,:))**2/12d0 + (3d0*f_charac_minus(i+1,:) - 4d0*f_charac_minus(i+2,:) + f_charac_minus(i+3,:))**2/4d0
    beta1_dash(:) = 13d0*(f_charac_minus(i,:) - 2d0*f_charac_minus(i+1,:) + f_charac_minus(i+2,:))**2/12d0 + (f_charac_minus(i,:) - f_charac_minus(i+2,:))**2/4d0
    beta2_dash(:) = 13d0*(f_charac_minus(i-1,:) - 2d0*f_charac_minus(i,:) + f_charac_minus(i+1,:))**2/12d0 + (f_charac_minus(i-1,:) - 4d0*f_charac_minus(i,:) + 3d0*f_charac_minus(i+1,:))**2/4d0

    alpha0(:) = d0/(epsilon + beta0(:))**2
    alpha1(:) = d1/(epsilon + beta1(:))**2
    alpha2(:) = d2/(epsilon + beta2(:))**2
    alpha0_dash(:) = d0_dash/(epsilon + beta0_dash(:))**2
    alpha1_dash(:) = d1_dash/(epsilon + beta1_dash(:))**2
    alpha2_dash(:) = d2_dash/(epsilon + beta2_dash(:))**2
    wt0(:) = alpha0(:)/(alpha0(:) + alpha1(:) + alpha2(:))
    wt1(:) = alpha1(:)/(alpha0(:) + alpha1(:) + alpha2(:))
    wt2(:) = alpha2(:)/(alpha0(:) + alpha1(:) + alpha2(:))
    wt0_dash(:) = alpha0_dash(:)/(alpha0_dash(:) + alpha1_dash(:) + alpha2_dash(:))
    wt1_dash(:) = alpha1_dash(:)/(alpha0_dash(:) + alpha1_dash(:) + alpha2_dash(:))
    wt2_dash(:) = alpha2_dash(:)/(alpha0_dash(:) + alpha1_dash(:) + alpha2_dash(:))

    temp(:) = wt0(:)*minus0(:) + wt1(:)*minus1(:) + wt2(:)*minus2(:)
    flux_plus = MATMUL(right_eigen,temp)
    temp(:) = wt0_dash(:)*plus0(:) + wt1_dash(:)*plus1(:) + wt2_dash(:)*plus2(:)
    flux_minus = MATMUL(right_eigen,temp)

    f_ih(:) = flux_plus(:) + flux_minus(:)

end subroutine wenox

!******************************************************************************************************************************************************************************************************





!******************************************************************************FLUX CONSTRUCTION IN Y DIRECTION****************************************************************************************

subroutine wenoy(i,j,k,g_ih)
    
    use variables
    implicit none
    integer, intent(in)                          :: i,j,k
    double precision, dimension(neq),intent(out) :: g_ih
    double precision, dimension(neq)             :: plus0,plus1,plus2,minus0,minus1,minus2,beta0,beta1,beta2,beta0_dash,beta1_dash,beta2_dash,temp                               
    double precision, dimension(neq)             :: alpha0, alpha1, alpha2, alpha0_dash, alpha1_dash, alpha2_dash
    double precision, dimension(neq)             :: wt0, wt1, wt2, wt0_dash, wt1_dash, wt2_dash
    double precision, dimension(neq,neq)         :: right_eigen, left_eigen
    double precision, dimension(neq)             :: flux_minus, flux_plus
    double precision, dimension(n_sp)            :: ro_r, ro_l, ro_avg
    double precision                             :: g_charac(j-2:j+3,neq), g_charac_plus(j-2:j+3,neq), g_charac_minus(j-2:j+3,neq)           !fluxes in characteristic form
    double precision                             :: q_charac_y(j-2:j+3,neq)
    double precision                             :: epsilon=1e-6, d0=3d0/10d0, d1=3d0/5d0, d2=1d0/10d0, d0_dash=1d0/10d0, d1_dash=3d0/5d0, d2_dash=3d0/10d0
    double precision                             :: u_r, u_l, v_r, v_l, w_r, w_l, ht_r, ht_l, a_r, a_l,Tv_r,Tv_l, u_avg, v_avg, w_avg, ht_avg, a_avg, p_pr, et_pr,Tv_avg,gamma
    integer                                      :: l, ll, entropy

    
    call primitive(i,j,k,ro_l,u_l,v_l,w_l,ht_l,a_l,Tv_l)
    call primitive(i,j+1,k,ro_r,u_r,v_r,w_r,ht_r,a_r,Tv_r)
    call averages(ro_r,ro_l,u_r,u_l,v_r,v_l,w_r,w_l,ht_r,ht_l,Tv_r,Tv_l,ro_avg,u_avg,v_avg,w_avg,ht_avg,a_avg,Tv_avg,gamma)
    call eigenvector(2,0,n_sp,neq,ro_avg(1:n_sp),u_avg,v_avg,w_avg,a_avg,ht_avg,h0(1:n_sp),(gamma-1d0),right_eigen)
    call eigenvector(2,1,n_sp,neq,ro_avg(1:n_sp),u_avg,v_avg,w_avg,a_avg,ht_avg,h0(1:n_sp),(gamma-1d0),left_eigen)


    if(j==firsty) then
        do l=j,j+3,1
            temp(:) = fluxg(i,l,k,:)
            temp(n_sp+2) = fluxg(i,l,k,n_sp+2)+p(i,l,k)
            g_charac(l,:) = MATMUL(left_eigen,temp)
        end do
        do l=j-2,j-1,1
            temp(:) = (fluxg(i,firsty,k,:))*(10d0*abs(firsty-l))**12d0
            temp(n_sp+2) = (fluxg(i,firsty,k,n_sp+2)+p(i,firsty,k))*(10d0*abs(firsty-l))**12d0
            g_charac(l,:) = MATMUL(left_eigen,temp)
        end do
    else if(j==firsty+1) then
        do l=j-1,j+3,1
            temp(:) = fluxg(i,l,k,:)
            temp(n_sp+2) = fluxg(i,l,k,n_sp+2)+p(i,l,k)
            g_charac(l,:) = MATMUL(left_eigen,temp)
        end do
        l=j-2
            temp(:) = (fluxg(i,firsty,k,:))*(10d0*(abs(firsty-l)-1))**12d0
            temp(n_sp+2) = (fluxg(i,firsty,k,n_sp+2)+p(i,firsty,k))*(10d0*(abs(firsty-l)-1))**12d0
            g_charac(l,:) = MATMUL(left_eigen,temp)
    else if(j==lasty) then
        do l=j-2,j,1
            temp(:) = fluxg(i,l,k,:)
            temp(n_sp+2) = fluxg(i,l,k,n_sp+2)+p(i,l,k)
            g_charac(l,:) = MATMUL(left_eigen,temp)
        end do
        do l=j+1,j+3,1
            temp(:) = (fluxg(i,lasty,k,:))*(10d0*abs(lasty-l))**12d0
            temp(n_sp+2) = (fluxg(i,lasty,k,n_sp+2)+p(i,lasty,k))*(10d0*abs(lasty-l))**12d0
            g_charac(l,:) = MATMUL(left_eigen,temp)
        end do
    else if(j==lasty-1) then
        do l=j-2,j+1,1
            temp(:) = fluxg(i,l,k,:)
            temp(n_sp+2) = fluxg(i,l,k,n_sp+2)+p(i,l,k)
            g_charac(l,:) = MATMUL(left_eigen,temp)
        end do
        do l=j+2,j+3,1
            temp(:) = (fluxg(i,lasty,k,:))*(10d0*(abs(lasty-l)-1))**12d0
            temp(n_sp+2) = (fluxg(i,lasty,k,n_sp+2)+p(i,lasty,k))*(10d0*(abs(lasty-l)-1))**12d0
            g_charac(l,:) = MATMUL(left_eigen,temp)
        end do
    else if(j==lasty-2) then
        do l=j-2,j+2,1
            temp(:) = fluxg(i,l,k,:)
            temp(n_sp+2) = fluxg(i,l,k,n_sp+2)+p(i,l,k)
            g_charac(l,:) = MATMUL(left_eigen,temp)
        end do
        l=j+3
            temp(:) = (fluxg(i,lasty,k,:))*(10d0*(abs(lasty-l)-2))**12d0
            temp(n_sp+2) = (fluxg(i,lasty,k,n_sp+2)+p(i,lasty,k))*(10d0*(abs(lasty-l)-2))**12d0
            g_charac(l,:) = MATMUL(left_eigen,temp)
    else
        do l=j-2,j+3,1
            temp(:) = fluxg(i,l,k,:)
            temp(n_sp+2) = fluxg(i,l,k,n_sp+2)+p(i,l,k)
            g_charac(l,:) = MATMUL(left_eigen,temp)
        end do
    end if


    do l=j-2,j+3,1
        if(l<firsty) then
            q_charac_y(l,:) = MATMUL(left_eigen,q(i,firsty,k,:))
        else if(l>lasty) then
            q_charac_y(l,:) = MATMUL(left_eigen,q(i,lasty,k,:))
        else
            q_charac_y(l,:) = MATMUL(left_eigen,q(i,l,k,:))
        end if        
    end do


    do l=j-2,j+3,1
        g_charac_plus(l,:) = 0.5d0*(g_charac(l,:) + alpha_y_global(i,mod(myrank,zprocs),k,:)*q_charac_y(l,:))
        g_charac_minus(l,:) = 0.5d0*(g_charac(l,:) - alpha_y_global(i,mod(myrank,zprocs),k,:)*q_charac_y(l,:))
    end do


    minus0(:) = g_charac_plus(j,:)/3d0 + 5d0*g_charac_plus(j+1,:)/6d0 - g_charac_plus(j+2,:)/6d0
    minus1(:) = -g_charac_plus(j-1,:)/6d0 + 5d0*g_charac_plus(j,:)/6d0 + g_charac_plus(j+1,:)/3d0
    minus2(:) = g_charac_plus(j-2,:)/3d0 - 7d0*g_charac_plus(j-1,:)/6d0 + 11d0*g_charac_plus(j,:)/6d0

    plus0(:) = 11d0*g_charac_minus(j+1,:)/6d0 - 7d0*g_charac_minus(j+2,:)/6d0 + g_charac_minus(j+3,:)/3d0
    plus1(:) = g_charac_minus(j,:)/3d0 + 5d0*g_charac_minus(j+1,:)/6d0 - g_charac_minus(j+2,:)/6d0
    plus2(:) = -g_charac_minus(j-1,:)/6d0 + 5d0*g_charac_minus(j,:)/6d0 + g_charac_minus(j+1,:)/3d0

    beta0(:) = 13d0*(g_charac_plus(j,:) - 2d0*g_charac_plus(j+1,:) + g_charac_plus(j+2,:))**2/12d0 + (3d0*g_charac_plus(j,:) - 4d0*g_charac_plus(j+1,:) + g_charac_plus(j+2,:))**2/4d0
    beta1(:) = 13d0*(g_charac_plus(j-1,:) - 2d0*g_charac_plus(j,:) + g_charac_plus(j+1,:))**2/12d0 + (g_charac_plus(j-1,:) - g_charac_plus(j+1,:))**2/4d0
    beta2(:) = 13d0*(g_charac_plus(j-2,:) - 2d0*g_charac_plus(j-1,:) + g_charac_plus(j,:))**2/12d0 + (g_charac_plus(j-2,:) - 4d0*g_charac_plus(j-1,:) + 3d0*g_charac_plus(j,:))**2/4d0

    beta0_dash(:) = 13d0*(g_charac_minus(j+1,:) - 2d0*g_charac_minus(j+2,:) + g_charac_minus(j+3,:))**2/12d0 + (3d0*g_charac_minus(j+1,:) - 4d0*g_charac_minus(j+2,:) + g_charac_minus(j+3,:))**2/4d0
    beta1_dash(:) = 13d0*(g_charac_minus(j,:) - 2d0*g_charac_minus(j+1,:) + g_charac_minus(j+2,:))**2/12d0 + (g_charac_minus(j,:) - g_charac_minus(j+2,:))**2/4d0
    beta2_dash(:) = 13d0*(g_charac_minus(j-1,:) - 2d0*g_charac_minus(j,:) + g_charac_minus(j+1,:))**2/12d0 + (g_charac_minus(j-1,:) - 4d0*g_charac_minus(j,:) + 3d0*g_charac_minus(j+1,:))**2/4d0

    alpha0(:) = d0/(epsilon + beta0(:))**2
    alpha1(:) = d1/(epsilon + beta1(:))**2
    alpha2(:) = d2/(epsilon + beta2(:))**2
    alpha0_dash(:) = d0_dash/(epsilon + beta0_dash(:))**2
    alpha1_dash(:) = d1_dash/(epsilon + beta1_dash(:))**2
    alpha2_dash(:) = d2_dash/(epsilon + beta2_dash(:))**2
    wt0(:) = alpha0(:)/(alpha0(:) + alpha1(:) + alpha2(:))
    wt1(:) = alpha1(:)/(alpha0(:) + alpha1(:) + alpha2(:))
    wt2(:) = alpha2(:)/(alpha0(:) + alpha1(:) + alpha2(:))
    wt0_dash(:) = alpha0_dash(:)/(alpha0_dash(:) + alpha1_dash(:) + alpha2_dash(:))
    wt1_dash(:) = alpha1_dash(:)/(alpha0_dash(:) + alpha1_dash(:) + alpha2_dash(:))
    wt2_dash(:) = alpha2_dash(:)/(alpha0_dash(:) + alpha1_dash(:) + alpha2_dash(:))

    temp(:) = wt0(:)*minus0(:) + wt1(:)*minus1(:) + wt2(:)*minus2(:)
    flux_plus = MATMUL(right_eigen,temp)
    temp(:) = wt0_dash(:)*plus0(:) + wt1_dash(:)*plus1(:) + wt2_dash(:)*plus2(:)
    flux_minus = MATMUL(right_eigen,temp)

    g_ih(:) = flux_plus(:) + flux_minus(:)

end subroutine wenoy

!******************************************************************************************************************************************************************************************************






!******************************************************************************FLUX CONSTRUCTION IN Z DIRECTION****************************************************************************************

subroutine wenoz(i,j,k,h_ih)
    
    use variables
    implicit none
    integer, intent(in)                          :: i,j,k
    double precision, dimension(neq),intent(out) :: h_ih
    double precision, dimension(neq)             :: plus0,plus1,plus2,minus0,minus1,minus2,beta0,beta1,beta2,beta0_dash,beta1_dash,beta2_dash,temp                           
    double precision, dimension(neq)             :: alpha0, alpha1, alpha2, alpha0_dash, alpha1_dash, alpha2_dash
    double precision, dimension(neq)             :: wt0, wt1, wt2, wt0_dash, wt1_dash, wt2_dash
    double precision, dimension(neq,neq)         :: right_eigen, left_eigen
    double precision, dimension(neq)             :: flux_minus, flux_plus
    double precision, dimension(n_sp)            :: ro_r, ro_l, ro_avg
    double precision                             :: h_charac(k-2:k+3,neq), h_charac_plus(k-2:k+3,neq), h_charac_minus(k-2:k+3,neq)             !fluxes in characteristic form 
    double precision                             :: q_charac_z(k-2:k+3,neq)
    double precision                             :: epsilon=1e-6, d0=3d0/10d0, d1=3d0/5d0, d2=1d0/10d0, d0_dash=1d0/10d0, d1_dash=3d0/5d0, d2_dash=3d0/10d0
    double precision                             :: u_r, u_l, v_r, v_l, w_r, w_l, ht_r, ht_l, a_r, a_l,Tv_r,Tv_l, u_avg, v_avg, w_avg,ht_avg, a_avg, p_pr, et_pr,Tv_avg,gamma
    integer                                      :: l, ll, entropy

    
    call primitive(i,j,k,ro_l,u_l,v_l,w_l,ht_l,a_l,Tv_l)
    call primitive(i,j,k+1,ro_r,u_r,v_r,w_r,ht_r,a_r,Tv_r)
    call averages(ro_r,ro_l,u_r,u_l,v_r,v_l,w_r,w_l,ht_r,ht_l,Tv_r,Tv_l,ro_avg,u_avg,v_avg,w_avg,ht_avg,a_avg,Tv_avg,gamma)
    call eigenvector(3,0,n_sp,neq,ro_avg(1:n_sp),u_avg,v_avg,w_avg,a_avg,ht_avg,h0(1:n_sp),(gamma-1d0),right_eigen)
    call eigenvector(3,1,n_sp,neq,ro_avg(1:n_sp),u_avg,v_avg,w_avg,a_avg,ht_avg,h0(1:n_sp),(gamma-1d0),left_eigen)

    
    if(k==firstz) then
        do l=k,k+3,1
            temp(:) = fluxh(i,j,l,:)
            temp(n_sp+3) = fluxh(i,j,l,n_sp+3)+p(i,j,l)
            h_charac(l,:) = MATMUl(left_eigen,temp)
        end do
        do l=k-2,k-1,1
            temp(:) = (fluxh(i,j,firstz,:))*(10d0*abs(firstz-l))**12d0
            temp(n_sp+3) = (fluxh(i,j,firstz,n_sp+3)+p(i,j,firstz))*(10d0*abs(firstz-l))**12d0
            h_charac(l,:) = MATMUl(left_eigen,temp)
        end do
    else if(k==firstz+1) then
        do l=k-1,k+3,1
            temp(:) = fluxh(i,j,l,:)
            temp(n_sp+3) = fluxh(i,j,l,n_sp+3)+p(i,j,l)
            h_charac(l,:) = MATMUl(left_eigen,temp)
        end do
        l=k-2
            temp(:) = (fluxh(i,j,firstz,:))*(10d0*(abs(firstz-l)-1))**12d0
            temp(n_sp+3) = (fluxh(i,j,firstz,n_sp+3)+p(i,j,firstz))*(10d0*(abs(firstz-l)-1))**12d0
            h_charac(l,:) = MATMUl(left_eigen,temp)
    else if(k==lastz) then
        do l=k-2,k,1
            temp(:) = fluxh(i,j,l,:)
            temp(n_sp+3) = fluxh(i,j,l,n_sp+3)+p(i,j,l)
            h_charac(l,:) = MATMUl(left_eigen,temp)
        end do
        do l=k+1,k+3,1
            temp(:) = (fluxh(i,j,lastz,:))*(10d0*abs(lastz-l))**12d0
            temp(n_sp+3) = (fluxh(i,j,lastz,n_sp+3)+p(i,j,lastz))*(10d0*abs(lastz-l))**12d0
            h_charac(l,:) = MATMUl(left_eigen,temp)
        end do
    else if(k==lastz-1) then
        do l=k-2,k+1,1
            temp(:) = fluxh(i,j,l,:)
            temp(n_sp+3) = fluxh(i,j,l,n_sp+3)+p(i,j,l)
            h_charac(l,:) = MATMUl(left_eigen,temp)
        end do
        do l=k+2,k+3,1
            temp(:) = (fluxh(i,j,lastz,:))*(10d0*(abs(lastz-l)-1))**12d0
            temp(n_sp+3) = (fluxh(i,j,lastz,n_sp+3)+p(i,j,lastz))*(10d0*(abs(lastz-l)-1))**12d0
            h_charac(l,:) = MATMUl(left_eigen,temp)
        end do
    else if(k==lastz-2) then
        do l=k-2,k+2,1
            temp(:) = fluxh(i,j,l,:)
            temp(n_sp+3) = fluxh(i,j,l,n_sp+3)+p(i,j,l)
            h_charac(l,:) = MATMUl(left_eigen,temp)
        end do
        l=k+3
            temp(:) = (fluxh(i,j,lastz,:))*(10d0*(abs(lastz-l)-2))**12d0
            temp(n_sp+3) = (fluxh(i,j,lastz,n_sp+3)+p(i,j,lastz))*(10d0*(abs(lastz-l)-2))**12d0
            h_charac(l,:) = MATMUl(left_eigen,temp)
    else
        do l=k-2,k+3,1
            temp(:) = fluxh(i,j,l,:)
            temp(n_sp+3) = fluxh(i,j,l,n_sp+3)+p(i,j,l)
            h_charac(l,:) = MATMUl(left_eigen,temp)
        end do
    end if


    do l=k-2,k+3,1
        if(l<firstz) then
            q_charac_z(l,:) = MATMUl(left_eigen,q(i,j,firstz,:))
        else if(l>lastz) then
            q_charac_z(l,:) = MATMUl(left_eigen,q(i,j,lastz,:))
        else
            q_charac_z(l,:) = MATMUl(left_eigen,q(i,j,l,:))
        end if
    end do


    do l=k-2,k+3,1
        h_charac_plus(l,:) = 0.5d0*(h_charac(l,:) + alpha_z_global(i,j,myrank/zprocs,:)*q_charac_z(l,:))
        h_charac_minus(l,:) = 0.5d0*(h_charac(l,:) - alpha_z_global(i,j,myrank/zprocs,:)*q_charac_z(l,:))
    end do


    minus0(:) = h_charac_plus(k,:)/3d0 + 5d0*h_charac_plus(k+1,:)/6d0 - h_charac_plus(k+2,:)/6d0
    minus1(:) = -h_charac_plus(k-1,:)/6d0 + 5d0*h_charac_plus(k,:)/6d0 + h_charac_plus(k+1,:)/3d0
    minus2(:) = h_charac_plus(k-2,:)/3d0 - 7d0*h_charac_plus(k-1,:)/6d0 + 11d0*h_charac_plus(k,:)/6d0

    plus0(:) = 11d0*h_charac_minus(k+1,:)/6d0 - 7d0*h_charac_minus(k+2,:)/6d0 + h_charac_minus(k+3,:)/3d0
    plus1(:) = h_charac_minus(k,:)/3d0 + 5d0*h_charac_minus(k+1,:)/6d0 - h_charac_minus(k+2,:)/6d0
    plus2(:) = -h_charac_minus(k-1,:)/6d0 + 5d0*h_charac_minus(k,:)/6d0 + h_charac_minus(k+1,:)/3d0

    beta0(:) = 13d0*(h_charac_plus(k,:) - 2d0*h_charac_plus(k+1,:) + h_charac_plus(k+2,:))**2/12d0 + (3d0*h_charac_plus(k,:) - 4d0*h_charac_plus(k+1,:) + h_charac_plus(k+2,:))**2/4d0
    beta1(:) = 13d0*(h_charac_plus(k-1,:) - 2d0*h_charac_plus(k,:) + h_charac_plus(k+1,:))**2/12d0 + (h_charac_plus(k-1,:) - h_charac_plus(k+1,:))**2/4d0
    beta2(:) = 13d0*(h_charac_plus(k-2,:) - 2d0*h_charac_plus(k-1,:) + h_charac_plus(k,:))**2/12d0 + (h_charac_plus(k-2,:) - 4d0*h_charac_plus(k-1,:) + 3d0*h_charac_plus(k,:))**2/4d0

    beta0_dash(:) = 13d0*(h_charac_minus(k+1,:) - 2d0*h_charac_minus(k+2,:) + h_charac_minus(k+3,:))**2/12d0 + (3d0*h_charac_minus(k+1,:) - 4d0*h_charac_minus(k+2,:) + h_charac_minus(k+3,:))**2/4d0
    beta1_dash(:) = 13d0*(h_charac_minus(k,:) - 2d0*h_charac_minus(k+1,:) + h_charac_minus(k+2,:))**2/12d0 + (h_charac_minus(k,:) - h_charac_minus(k+2,:))**2/4d0
    beta2_dash(:) = 13d0*(h_charac_minus(k-1,:) - 2d0*h_charac_minus(k,:) + h_charac_minus(k+1,:))**2/12d0 + (h_charac_minus(k-1,:) - 4d0*h_charac_minus(k,:) + 3d0*h_charac_minus(k+1,:))**2/4d0

    alpha0(:) = d0/(epsilon + beta0(:))**2
    alpha1(:) = d1/(epsilon + beta1(:))**2
    alpha2(:) = d2/(epsilon + beta2(:))**2
    alpha0_dash(:) = d0_dash/(epsilon + beta0_dash(:))**2
    alpha1_dash(:) = d1_dash/(epsilon + beta1_dash(:))**2
    alpha2_dash(:) = d2_dash/(epsilon + beta2_dash(:))**2
    wt0(:) = alpha0(:)/(alpha0(:) + alpha1(:) + alpha2(:))
    wt1(:) = alpha1(:)/(alpha0(:) + alpha1(:) + alpha2(:))
    wt2(:) = alpha2(:)/(alpha0(:) + alpha1(:) + alpha2(:))
    wt0_dash(:) = alpha0_dash(:)/(alpha0_dash(:) + alpha1_dash(:) + alpha2_dash(:))
    wt1_dash(:) = alpha1_dash(:)/(alpha0_dash(:) + alpha1_dash(:) + alpha2_dash(:))
    wt2_dash(:) = alpha2_dash(:)/(alpha0_dash(:) + alpha1_dash(:) + alpha2_dash(:))

    temp(:) = wt0(:)*minus0(:) + wt1(:)*minus1(:) + wt2(:)*minus2(:)
    flux_plus = MATMUL(right_eigen,temp)
    temp(:) = wt0_dash(:)*plus0(:) + wt1_dash(:)*plus1(:) + wt2_dash(:)*plus2(:)
    flux_minus = MATMUL(right_eigen,temp)
    
    h_ih(:) = flux_plus(:) + flux_minus(:)

end subroutine wenoz

!******************************************************************************************************************************************************************************************************




!***********************************************************************************CALCULATES EIGENVECTORS********************************************************************************************

subroutine eigen_val(u,a,num,eigen)
    
    implicit none
    integer, intent(in)                           :: num
    double precision, intent(in)                  :: u,a
    double precision, dimension(num), intent(out) :: eigen
    
    !ideal gas
    if(num==5) then
        eigen(1) = u - a
        eigen(2) = u
        eigen(3) = u + a
        eigen(4) = u
        eigen(5) = u
    else if(num==6) then
        eigen(1) = u
        eigen(2) = u
        eigen(3) = u
        eigen(4) = u
        eigen(5) = u + a
        eigen(6) = u - a
    end if

end subroutine eigen_val

!******************************************************************************************************************************************************************************************************





!******************************************************************CALCULATES PRIMITIVE VARIABLES GIVEN CONSERVATIVE VARIABLES*************************************************************************

subroutine primitive(i,j,k,ro_pr,u_pr,v_pr,w_pr,ht_pr,a_pr,Tv_pr)
    
    use variables
    implicit none
    integer, intent(in)                            :: i,j,k
    double precision, intent(out)                  :: u_pr,v_pr,w_pr,ht_pr,a_pr,Tv_pr
    double precision, dimension(n_sp), intent(out) :: ro_pr
    double precision, dimension(n_sp)              :: one=1d0
    double precision                               :: ro_tot,Trtr_pr,e_pr,h_pr,p_pr,et_pr,add,R_s

    ro_pr(:) = q(i,j,k,1:n_sp)
    u_pr = u(i,j,k)
    v_pr = v(i,j,k)
    w_pr = w(i,j,k)
    ht_pr = q(i,j,k,n_sp+4)/ro(i,j,k) + p(i,j,k)/ro(i,j,k)
    R_s = R_u*add(ro_pr(:),Mw(:),n_sp,-1)/ro(i,j,k)
    a_pr = sqrt(R_s*T(i,j,k)*(1 + R_s*ro(i,j,k)/add(ro_pr(:),Cv(:),n_sp,1)))

    Tv_pr = Tv(i,j,k)    

end subroutine primitive

!******************************************************************************************************************************************************************************************************





!********************************************************************************CALCULATES ROE AVERAGES***********************************************************************************************

subroutine averages(ro_r,ro_l,u_r,u_l,v_r,v_l,w_r,w_l,ht_r,ht_l,Tv_r,Tv_l,ro_avg,u_avg,v_avg,w_avg,ht_avg,a_avg,Tv_avg,gamma)
    
    use variables
    implicit none
    double precision, dimension(n_sp), intent(in)  :: ro_r,ro_l
    double precision, intent(in)                   :: u_r,u_l,v_r,v_l,w_r,w_l,ht_r,ht_l,Tv_r,Tv_l
    double precision, dimension(n_sp), intent(out) :: ro_avg
    double precision, intent(out)                  :: u_avg,v_avg,w_avg,ht_avg,a_avg,Tv_avg,gamma
    double precision, dimension(n_sp)              :: one=1d0
    double precision, dimension(n_di)              :: ev_s   
    double precision                               :: k,average,ro_tot_avg,add,R_s

    ro_avg(:) = (ro_r(:)+ro_l(:))/2d0
    ro_tot_avg = add(ro_avg,one,n_sp,1)
    u_avg = (u_r+u_l)/2d0
    v_avg = (v_r+v_l)/2d0
    w_avg = (w_r+w_l)/2d0
    ht_avg = (ht_r+ht_l)/2d0
    Tv_avg = (Tv_r+Tv_l)/2d0
    
    k = (u_avg*u_avg + v_avg*v_avg + w_avg*w_avg)/2d0
    ev_s(1:n_di) = vib_flag*((R_u/Mw(1:n_di))*theta_v(:)/(exp(theta_v(:)/Tv_avg)-1))
    
    R_s = R_u*add(ro_avg(:),Mw(:),n_sp,-1)/ro_tot_avg
    gamma = (1 + R_s*ro_tot_avg/add(ro_avg(:),Cv(:),n_sp,1))    
    a_avg = sqrt((gamma-1d0)*(ht_avg - k - chem_flag*add(ro_avg,h0,n_sp,1)/ro_tot_avg - vib_flag*add(ro_avg(1:n_di),ev_s(1:n_di),n_di,1)/ro_tot_avg))
    
end subroutine averages

!******************************************************************************************************************************************************************************************************





!roe averages I think will not be valid for the reactive system, so will use simple mean here
! !********************************************************************************CALCULATES ROE AVERAGES***********************************************************************************************

! subroutine roe_averages(ro_r,ro_l,u_r,u_l,v_r,v_l,w_r,w_l,ht_r,ht_l,u_avg,v_avg,w_avg,ro_avg,ht_avg,a_avg)
    
!     use variables
!     implicit none
!     double precision, intent(in)  :: ro_r, ro_l, u_r, u_l, v_r, v_l, w_r, w_l, ht_r, ht_l
!     double precision, intent(out) :: u_avg, v_avg, w_avg, ro_avg, ht_avg, a_avg
!     double precision              :: vt, average

!     ro_avg = sqrt(ro_r*ro_l)
!     u_avg = average(ro_r,ro_l,u_r,u_l)
!     v_avg = average(ro_r,ro_l,v_r,v_l)
!     w_avg = average(ro_r,ro_l,w_r,w_l)
!     ht_avg = average(ro_r,ro_l,ht_r,ht_l)
!     vt = sqrt((u_avg*u_avg + v_avg*v_avg + w_avg*w_avg))
!     a_avg = sqrt((gamma-1)*(ht_avg-0.5d0*vt*vt))

! end subroutine roe_averages

! !******************************************************************************************************************************************************************************************************





! !******************************************************************************function USED IN ROE AVERAGES*******************************************************************************************

! double precision function average(ro_r,ro_l,fr,fl)
    
!     implicit none
!     double precision, intent(in) :: ro_r, ro_l, fr, fl
    
!     average = (sqrt(ro_l)*fl + sqrt(ro_r)*fr)/(sqrt(ro_r) + sqrt(ro_l))
!     return

! end function average

! !******************************************************************************************************************************************************************************************************


!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************










!******************************************************************************************************************************************************************************************************
!**************************************************************** *COMPUTES CHEMICAL SOURCE TERMS FOR THE DENSITY EQUATIONS****************************************************************************
!******************************************************************************************************************************************************************************************************

subroutine source(ws)
  
    use variables
    implicit none
    double precision, dimension(sx,sy,sz,n_sp), intent(out) :: ws
    double precision, dimension(17)                         :: k_eq,k_f,k_b
    double precision, dimension(n_sp)                       :: ro_s
    double precision                                        :: Z_s,R1,R2,R3,R4,R5
    double precision                                        :: sum, T_eff
    integer                                                 :: i,j,k,m

    do k=startz,endz
        do j=starty,endy
            do i=startx,endx
                T_eff = sqrt(T(i,j,k)*Tv(i,j,k))

                do m=1,17
                    if(vib_flag==1 .and. m<=15) then
                        k_f(m) = cf(m)*T_eff**eta(m)*exp(-theta(m)/T_eff)
                    else
                        k_f(m) = cf(m)*T(i,j,k)**eta(m)*exp(-theta(m)/T(i,j,k))
                    end if
                    Z_s = 10000d0/T(i,j,k)
                    k_eq(m) = cm(m)*exp(a1(m)/Z_s + a2(m) + a3(m)*log(Z_s) + a4(m)*Z_s + a5(m)*Z_s**2)
                    k_b(m) = (cf(m)*T(i,j,k)**eta(m)*exp(-theta(m)/T(i,j,k)))/k_eq(m)
                end do

                R1 = 0d0
                R2 = 0d0
                R3 = 0d0
                do m = 1,5
                    R1 = R1 + k_b(m)*(ro_s(4)/Mw(4))*(ro_s(4)/Mw(4))*(ro_s(m)/Mw(m)) - k_f(m)*(ro_s(1)/Mw(1))*(ro_s(m)/Mw(m))
                    R2 = R2 + k_b(m+5)*(ro_s(5)/Mw(5))*(ro_s(5)/Mw(5))*(ro_s(m)/Mw(m)) - k_f(m+5)*(ro_s(2)/Mw(2))*(ro_s(m)/Mw(m))
                    R3 = R3 + k_b(m+10)*(ro_s(4)/Mw(4))*(ro_s(5)/Mw(5))*(ro_s(m)/Mw(m)) - k_f(m+10)*(ro_s(3)/Mw(3))*(ro_s(m)/Mw(m))
                end do
                R4 = k_b(16)*(ro_s(3)/Mw(3))*(ro_s(4)/Mw(4)) - k_f(16)*(ro_s(1)/Mw(1))*(ro_s(5)/Mw(5))
                R5 = k_b(17)*(ro_s(2)/Mw(2))*(ro_s(4)/Mw(4)) - k_f(17)*(ro_s(3)/Mw(3))*(ro_s(5)/Mw(5))

                ws(i,j,k,1) = Mw(1)*(R1 + R4)
                ws(i,j,k,2) = Mw(2)*(R2 - R5)
                ws(i,j,k,3) = Mw(3)*(R3 - R4 + R5)
                ws(i,j,k,4) = Mw(4)*(-2d0*R1 - R3 - R4 - R5)
                ws(i,j,k,5) = Mw(5)*(-2d0*R2 - R3 + R4 + R5)
            end do
        end do
    end do

end subroutine source

!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************










!******************************************************************************************************************************************************************************************************
!*****************************************************************************COMPUTES VIBRATIONAL SOURCE TERM*****************************************************************************************
!******************************************************************************************************************************************************************************************************

subroutine vib_source(ws,wv)
  
    use variables
    implicit none
    double precision, dimension(sx,sy,sz,n_sp), intent(in) :: ws
    double precision, dimension(sx,sy,sz), intent(out)     :: wv
    double precision, dimension(n_sp)                      :: ro_s,X_t,tow_sr,A_sr,mu_sr,one
    double precision, dimension(n_di)                      :: ev_s,ev_s_star,tow_s,Q_tv
    double precision                                       :: sum,Mw_mix,add
    integer                                                :: i,j,k,r,s

    one(:) = 1d0

    do k=startz,endz
        do j=starty,endy
            do i=startx,endx
                ro_s(:) = q(i,j,k,1:n_sp)
                
                Mw_mix = add(ro_s,Mw,n_sp,-1)/ro(i,j,k)

                do s=1,3
                    ev_s(s)       = (R_u/Mw(s))*(theta_v(s)/(exp(theta_v(s)/Tv(i,j,k))-1))
                    ev_s_star(s)  = (R_u/Mw(s))*(theta_v(s)/(exp(theta_v(s)/T(i,j,k))-1))
                
                    do r=1,5
                        X_t(r) = ro_s(r)*Mw_mix/(Mw(r)*ro(i,j,k))
                        mu_sr(r) = Mw(s)*Mw(r)/(Mw(s) + Mw(r))
                        A_sr(r) = 1.16e-3*sqrt(mu_sr(r))*theta_v(s)**(4d0/3)
                        tow_sr(r) = exp(A_sr(r)*(T(i,j,k)**(-1d0/3) - 0.015*mu_sr(r)**(1d0/4)) - 18.42)/(p(i,j,k)/101325d0)
                    end do

                    tow_s(s) = add(X_t,one,n_sp,1)/add(X_t,tow_sr,n_sp,-1)
                
                    Q_tv(s) = ro_s(s)*(ev_s_star(s) - ev_s(s))/tow_s(s)
                end do

                wv(i,j,k) = Q_tv(1) + Q_tv(2) + Q_tv(3) + ws(i,j,k,1)*ev_s(1) + ws(i,j,k,2)*ev_s(2) + ws(i,j,k,3)*ev_s(3)
            end do
        end do
    end do

end subroutine vib_source

!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************










!******************************************************************************************************************************************************************************************************
!**************************************************************************BOUNDARY CONDITIONS FOR VISCOUS FLUXES**************************************************************************************
!******************************************************************************************************************************************************************************************************
!This subroutine is to avoid computing viscous stresses and heat fluxes using one-sided derivatives when the same could be obtained using central derivatives when boundary conditions are periodic
!SEE message_pass FOR q_side and q_top
!WE don't have to send all three visf,visg,vish in all three direction as done in older codes visf will be used only in x BC, likewise for visg ang vish

subroutine bc_viscous(rkl)
  
    use variables
    implicit none
    include 'mpif.h'
    integer            :: i,j,k,m,tag1,tag2,tag3,tag4
    integer            :: status(mpi_status_size),q_top,q_side
    integer,intent(in) :: rkl
    
    call MPI_TYPE_VECTOR(sz*neq,3*sx,sx*sy,MPI_REAL8,q_top,mperr)
    call MPI_TYPE_COMMIT(q_top,mperr)
    call MPI_TYPE_VECTOR(neq,sx*sy*3,sx*sy*sz,MPI_REAL8,q_side,mperr)
    call MPI_TYPE_COMMIT(q_side,mperr)

    tag1=100
    tag2=200
    tag3=300
    tag4=400

    !PERIODIC BOUNDARY CONDITION IN Z DIRECTION
    if(periodic_z==1) then
        if (mod(myrank,zprocs) == 0) then
            call MPI_SEND(vish(1,1,4,1),1,q_side,myrank+zprocs-1,tag1,MPI_COMM_WORLD,mperr)
        end if
        if (mod(myrank,zprocs) == zprocs-1) then
            call MPI_RECV(vish(1,1,sz-2,1),1,q_side,myrank-zprocs+1,tag1,MPI_COMM_WORLD,status,mperr)
            call MPI_SEND(vish(1,1,sz-5,1),1,q_side,myrank-zprocs+1,tag2,MPI_COMM_WORLD,mperr)
        end if
        if (mod(myrank,zprocs) == 0) then
            call MPI_RECV(vish(1,1,1,1),1,q_side,myrank+zprocs-1,tag2,MPI_COMM_WORLD,status,mperr)
        end if
    end if

    !PERIODIC BOUNDARY CONDITION IN Y DIRECTION
    if(periodic_y==1) then
        if (myrank/zprocs == 0) then
            call MPI_SEND(visg(1,4,1,1),1,q_top,myrank+zprocs*(yprocs-1),tag3,MPI_COMM_WORLD,mperr)
        end if
        if (myrank/zprocs == yprocs-1) then
            call MPI_RECV(visg(1,sy-2,1,1),1,q_top,myrank-zprocs*(yprocs-1),tag3,MPI_COMM_WORLD,status,mperr)
            call MPI_SEND(visg(1,sy-5,1,1),1,q_top,myrank-zprocs*(yprocs-1),tag4,MPI_COMM_WORLD,mperr)
        end if
        if (myrank/zprocs == 0) then
            call MPI_RECV(visg(1,1,1,1),1,q_top,myrank+zprocs*(yprocs-1),tag4,MPI_COMM_WORLD,status,mperr)
        end if
    end if
    
    !PERIODIC BOUNDARY CONDITION IN X DIRECTION
    if(periodic_x==1) then
        visf(sx-2,:,:,:) = visf(4,:,:,:)
        visf(sx-1,:,:,:) = visf(5,:,:,:)
        visf(sx,:,:,:) = visf(6,:,:,:) 
        visf(1,:,:,:) = visf(sx-5,:,:,:)
        visf(2,:,:,:) = visf(sx-4,:,:,:)
        visf(3,:,:,:) = visf(sx-3,:,:,:)
    end if

end subroutine bc_viscous

!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************









!******************************************************************************************************************************************************************************************************
!****************************************************************FOR INTER-PROCESSOR COMMUNICATION OF CONSERVATIVE VARIABLES***************************************************************************
!******************************************************************************************************************************************************************************************************
!q_top AND q_side ARE NEW DATA TYPES USED FOR SENDING DATA HERE AND FOR PERIODIC BC'SENDING
!q OR visf ETC. WILL HAVE neq*SZ*SY*SX ELEMENTS
!q_side IS A NEW DATA TYPE WHICH WILL HAVE neq*3*SY*SX ELEMENTS WHICH CORRESPOND TO ALL THE ninE COMPONENTS OF q OR visf ETC. FOR ALL THE POINTS IN 3 CONSECUTIVE X-Y PLANES
!q_top IS A NEW DATA TYPE WHICH WILL HAVE neq*SZ*3*SX ELEMENTS WHICH CORRESPOND TO ALL THE ninE COMPONENTS OF q OR visf ETC. FOR ALL THE POINTS IN 3 CONSECUTIVE X-Z PLANES
!GOOGLE "MPI_TYPE_VECTOR" AND "HOW IS DATA STORED IN MEMORY IN FORTRAN" FOR MORE

subroutine message_pass(rkl)
    use variables
    implicit none
    include 'mpif.h'
    integer            :: i,j,k,m,dest,tag1,tag2,tag3,tag4,uni
    integer            :: status(mpi_status_size),q_top,q_side
    integer,intent(in) :: rkl

    call MPI_TYPE_VECTOR(sz*neq,3*sx,sx*sy,MPI_REAL8,q_top,mperr)
    call MPI_TYPE_COMMIT(q_top,mperr)
    call MPI_TYPE_VECTOR(neq,sx*sy*3,sx*sy*sz,MPI_REAL8,q_side,mperr)
    call MPI_TYPE_COMMIT(q_side,mperr)

    tag1=10
    tag2=20
    tag3=30
    tag4=40
  
    !MESSAGE PASSING IN Y DIRECTION
    !SEND TO BOTTOM,RECV FROM BOTTOM
    if (myrank/zprocs > 0) then
        call MPI_SEND(q(1,4,1,1),1,q_top,myrank-zprocs,tag1,MPI_COMM_WORLD,mperr)
        call MPI_RECV(q(1,1,1,1),1,q_top,myrank-zprocs,tag2,MPI_COMM_WORLD,status,mperr)
    end if
    !SEND TO TOP,RECV FROM TOP
    if(myrank/zprocs < yprocs-1) then
        call MPI_RECV(q(1,sy-2,1,1),1,q_top,myrank+zprocs,tag1,MPI_COMM_WORLD,status,mperr)
        call MPI_SEND(q(1,sy-5,1,1),1,q_top,myrank+zprocs,tag2,MPI_COMM_WORLD,mperr)     
    end if


    !MESSAGE PASSING IN Z DIRECTION
    !SEND TO LEFT
    if (mod(myrank,zprocs) /= 0) then
        call MPI_SEND(q(1,1,4,1),1,q_side,myrank-1,tag3,MPI_COMM_WORLD,mperr)
        call MPI_RECV(q(1,1,1,1),1,q_side,myrank-1,tag4,MPI_COMM_WORLD,status,mperr)
    end if
    !SEND TO RIGHT
    if (mod(myrank,zprocs) /= zprocs-1) then
        call MPI_RECV(q(1,1,sz-2,1),1,q_side,myrank+1,tag3,MPI_COMM_WORLD,status,mperr)
        call MPI_SEND(q(1,1,sz-5,1),1,q_side,myrank+1,tag4,MPI_COMM_WORLD,mperr)       
    end if

end subroutine message_pass

!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************









!******************************************************************************************************************************************************************************************************
!*******************************************************************FOR INTER-PROCESSOR COMMUNICATION OF VISCOUS FLUXES********************************************************************************
!******************************************************************************************************************************************************************************************************
!This subroutine is to avoid computing viscous stresses and heat fluxes using one-sided derivatives at the inter-processor boundaries
!WE don't have to send all three visf,visg,vish in all three direction as done in older codes. only visf will be used in ddx, likewise for visg ang vish

subroutine message_pass_viscous(rkl)

    use variables
    implicit none
    include 'mpif.h'
    integer            :: i,j,k,m,dest,tag1,tag2,tag3,tag4,uni
    integer            :: status(mpi_status_size),q_top,q_side
    integer,intent(in) :: rkl

    call MPI_TYPE_VECTOR(sz*neq,3*sx,sx*sy,MPI_REAL8,q_top,mperr)
    call MPI_TYPE_COMMIT(q_top,mperr)
    call MPI_TYPE_VECTOR(neq,sx*sy*3,sx*sy*sz,MPI_REAL8,q_side,mperr)
    call MPI_TYPE_COMMIT(q_side,mperr)

    tag1=10
    tag2=20
    tag3=30
    tag4=40
  
    !MESSAGE PASSING IN Y DIRECTION
    if (myrank/zprocs > 0) then
        call MPI_SEND(visg(1,4,1,1),1,q_top,myrank-zprocs,tag1,MPI_COMM_WORLD,mperr)
        call MPI_RECV(visg(1,1,1,1),1,q_top,myrank-zprocs,tag2,MPI_COMM_WORLD,status,mperr)
    end if
    if(myrank/zprocs < yprocs-1) then
        call MPI_RECV(visg(1,sy-2,1,1),1,q_top,myrank+zprocs,tag1,MPI_COMM_WORLD,status,mperr)
        call MPI_SEND(visg(1,sy-5,1,1),1,q_top,myrank+zprocs,tag2,MPI_COMM_WORLD,mperr)     
    end if

    !MESSAGE PASSING IN Z DIRECTION
    if (mod(myrank,zprocs) /= 0) then
        call MPI_SEND(vish(1,1,4,1),1,q_side,myrank-1,tag3,MPI_COMM_WORLD,mperr)
        call MPI_RECV(vish(1,1,1,1),1,q_side,myrank-1,tag4,MPI_COMM_WORLD,status,mperr)
    end if
    if (mod(myrank,zprocs) /= zprocs-1) then
        call MPI_RECV(vish(1,1,sz-2,1),1,q_side,myrank+1,tag3,MPI_COMM_WORLD,status,mperr)
        call MPI_SEND(vish(1,1,sz-5,1),1,q_side,myrank+1,tag4,MPI_COMM_WORLD,mperr)       
    end if  

end subroutine message_pass_viscous

!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************





















!******************************************************************************************************************************************************************************************************
!*************************************************************************INITIALISE EVERYTHING TO ZERO TO BE SAFE*************************************************************************************
!******************************************************************************************************************************************************************************************************

subroutine intialise()

    use variables
    implicit none
    
    ddx_fluxf = 0d0
    ddy_fluxg = 0d0
    ddz_fluxh = 0d0
    visf = 0d0
    visg = 0d0
    vish = 0d0
    fluxf = 0d0
    fluxg = 0d0
    fluxh = 0d0
    x = 0d0
    y = 0d0
    z = 0d0
    y_ph = 0d0
    dy_dy_ph = 0d0

end subroutine intialise

!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************










!******************************************************************************************************************************************************************************************************
!*****************************************************************************COMPUTES VIBRATIONAL SOURCE TERM*****************************************************************************************
!******************************************************************************************************************************************************************************************************

double precision function vib_temp(i,j,k)
  
    use variables
    implicit none
    integer, intent(in)     :: i,j,k
    integer                 :: m
    double precision        :: func,func_der,eps,Ev,x_i,x_ip1
    double precision, dimension(n_sp) :: ro_s

    ro_s(:) = q(i,j,k,1:n_sp)

    x_i = Tv(i,j,k)
    Ev = q(i,j,k,n_sp+5)
    m=0
    eps = 10d0
    
    do while (abs(eps) > 1e-6)
    ! do m=1,4 !This is what shankar did in his code
        func = ro_s(1)*(R_u/Mw(1))*(theta_v(1)/(exp(theta_v(1)/x_i)-1)) + ro_s(2)*(R_u/Mw(2))*(theta_v(2)/(exp(theta_v(2)/x_i)-1)) + ro_s(3)*(R_u/Mw(3))*(theta_v(3)/(exp(theta_v(3)/x_i)-1)) - Ev
        func_der = ro_s(1)*(R_u/Mw(1))*(theta_v(1)/x_i)**2*(exp(theta_v(1)/x_i)/(exp(theta_v(1)/x_i)-1)**2) + ro_s(2)*(R_u/Mw(2))*(theta_v(2)/x_i)**2*(exp(theta_v(2)/x_i)/(exp(theta_v(2)/x_i)-1)**2) + ro_s(3)*(R_u/Mw(3))*(theta_v(3)/x_i)**2*(exp(theta_v(3)/x_i)/(exp(theta_v(3)/x_i)-1)**2)
        
        x_ip1 = x_i - func/func_der
        eps = (x_ip1 - x_i)/x_i
        x_i = x_ip1

        m = m + 1
        if(m>1000000) then
            write(*,*)'convergence problem in "vib_temp"'
            write(*,*)iter,i,j,k,m,func,eps
            stop
        end if
    end do                
    
    vib_temp = x_i              
    
end function vib_temp

!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************










!******************************************************************************************************************************************************************************************************
!************************************************************FOR VARIOUS SUMS INVOLVING PRODUCTS OF DENSITIES AND SOME PROPERTY************************************************************************
!******************************************************************************************************************************************************************************************************

double precision function add(var1,var2,sp,flag)

    use variables
    implicit none
    integer,intent(in)              :: sp,flag !1 for sum of products, -1 for sum of divisions
    double precision,dimension(sp),intent(in) :: var1,var2
    integer                         :: i

    if(abs(flag)/=1) write(*,*)'problem with flag in function add'

    add = 0d0
    do i=1,sp
        add = add + var1(i)*var2(i)**flag
    end do    

end function add

!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************










!******************************************************************************************************************************************************************************************************
!***********************************************************************************COMPUTES DERIVATIVES***********************************************************************************************
!******************************************************************************************************************************************************************************************************



!******************************************************************************FOR DERIVATIVES IN i DIRECTION******************************************************************************************

double precision function d6dx(var,i,j,k)

    use variables
    implicit none
    integer,intent(in)                    :: i,j,k
    double precision,dimension(sx,sy,sz),intent(in) :: var

    if(i==firstx) then
        d6dx = (-10d0*var(i+6,j,k)+72d0*var(i+5,j,k)-225d0*var(i+4,j,k)+400d0*var(i+3,j,k)-450d0*var(i+2,j,k)+360d0*var(i+1,j,k)-147d0*var(i,j,k))/(60d0*dx)
    else if(i==firstx+1) then
        d6dx = (2d0*var(i+5,j,k)-15d0*var(i+4,j,k)+50d0*var(i+3,j,k)-100d0*var(i+2,j,k)+150d0*var(i+1,j,k)-77d0*var(i,j,k)-10d0*var(i-1,j,k))/(60d0*dx)
    else if(i==firstx+2) then
        d6dx = (-var(i+4,j,k)+8d0*var(i+3,j,k)-30d0*var(i+2,j,k)+80d0*var(i+1,j,k)-35d0*var(i,j,k)-24d0*var(i-1,j,k)+2d0*var(i-2,j,k))/(60d0*dx)
    else if(i==lastx) then
        d6dx = -(-10d0*var(i-6,j,k)+72d0*var(i-5,j,k)-225d0*var(i-4,j,k)+400d0*var(i-3,j,k)-450d0*var(i-2,j,k)+360d0*var(i-1,j,k)-147d0*var(i,j,k))/(60d0*dx)
    else if(i==lastx-1) then
        d6dx = -(2d0*var(i-5,j,k)-15d0*var(i-4,j,k)+50d0*var(i-3,j,k)-100d0*var(i-2,j,k)+150d0*var(i-1,j,k)-77d0*var(i,j,k)-10d0*var(i+1,j,k))/(60d0*dx)
    else if(i==lastx-2) then
        d6dx = -(-var(i-4,j,k)+8d0*var(i-3,j,k)-30d0*var(i-2,j,k)+80d0*var(i-1,j,k)-35d0*var(i,j,k)-24d0*var(i+1,j,k)+2d0*var(i+2,j,k))/(60d0*dx)
    else
        d6dx = (var(i+3,j,k)-9d0*var(i+2,j,k)+45d0*var(i+1,j,k)-45d0*var(i-1,j,k)+9d0*var(i-2,j,k)-var(i-3,j,k))/(60d0*dx)
    end if

end function d6dx

!******************************************************************************************************************************************************************************************************





!******************************************************************************FOR DERIVATIVES IN j DIRECTION******************************************************************************************

double precision function d6dy(var,i,j,k)

    use variables
    implicit none
    integer,intent(in)                    :: i,j,k
    double precision,dimension(sx,sy,sz),intent(in) :: var

    if(j==firsty) then
        d6dy = dy_dy_ph(i,j,k)*(-10d0*var(i,j+6,k)+72d0*var(i,j+5,k)-225d0*var(i,j+4,k)+400d0*var(i,j+3,k)-450d0*var(i,j+2,k)+360d0*var(i,j+1,k)-147d0*var(i,j,k))/(60d0*dy)
    else if(j==firsty+1) then
        d6dy = dy_dy_ph(i,j,k)*(2d0*var(i,j+5,k)-15d0*var(i,j+4,k)+50d0*var(i,j+3,k)-100d0*var(i,j+2,k)+150d0*var(i,j+1,k)-77d0*var(i,j,k)-10d0*var(i,j-1,k))/(60d0*dy)
    else if(j==firsty+2) then
        d6dy = dy_dy_ph(i,j,k)*(-var(i,j+4,k)+8d0*var(i,j+3,k)-30d0*var(i,j+2,k)+80d0*var(i,j+1,k)-35d0*var(i,j,k)-24d0*var(i,j-1,k)+2d0*var(i,j-2,k))/(60d0*dy)
    else if(j==lasty) then
        d6dy = -dy_dy_ph(i,j,k)*(-10d0*var(i,j-6,k)+72d0*var(i,j-5,k)-225d0*var(i,j-4,k)+400d0*var(i,j-3,k)-450d0*var(i,j-2,k)+360d0*var(i,j-1,k)-147d0*var(i,j,k))/(60d0*dy)
    else if(j==lasty-1) then
        d6dy = -dy_dy_ph(i,j,k)*(2d0*var(i,j-5,k)-15d0*var(i,j-4,k)+50d0*var(i,j-3,k)-100d0*var(i,j-2,k)+150d0*var(i,j-1,k)-77d0*var(i,j,k)-10d0*var(i,j+1,k))/(60d0*dy)
    else if(j==lasty-2) then
        d6dy = -dy_dy_ph(i,j,k)*(-var(i,j-4,k)+8d0*var(i,j-3,k)-30d0*var(i,j-2,k)+80d0*var(i,j-1,k)-35d0*var(i,j,k)-24d0*var(i,j+1,k)+2d0*var(i,j+2,k))/(60d0*dy)
    else
        d6dy = dy_dy_ph(i,j,k)*(var(i,j+3,k)-9d0*var(i,j+2,k)+45d0*var(i,j+1,k)-45d0*var(i,j-1,k)+9d0*var(i,j-2,k)-var(i,j-3,k))/(60d0*dy)
    end if

end function d6dy

!******************************************************************************************************************************************************************************************************





!******************************************************************************FOR DERIVATIVES IN k DIRECTION******************************************************************************************

double precision function d6dz(var,i,j,k)

    use variables
    implicit none
    integer,intent(in) :: i,j,k
    double precision,dimension(sx,sy,sz),intent(in) :: var 

    if(k==firstz) then
        d6dz = (-10d0*var(i,j,k+6)+72d0*var(i,j,k+5)-225d0*var(i,j,k+4)+400d0*var(i,j,k+3)-450d0*var(i,j,k+2)+360d0*var(i,j,k+1)-147d0*var(i,j,k))/(60d0*dz)
    else if(k==firstz+1) then
        d6dz = (2d0*var(i,j,k+5)-15d0*var(i,j,k+4)+50d0*var(i,j,k+3)-100d0*var(i,j,k+2)+150d0*var(i,j,k+1)-77d0*var(i,j,k)-10d0*var(i,j,k-1))/(60d0*dz)
    else if(k==firstz+2) then
        d6dz = (-var(i,j,k+4)+8d0*var(i,j,k+3)-30d0*var(i,j,k+2)+80d0*var(i,j,k+1)-35d0*var(i,j,k)-24d0*var(i,j,k-1)+2d0*var(i,j,k-2))/(60d0*dz)
    else if(k==lastz) then
        d6dz = -(-10d0*var(i,j,k-6)+72d0*var(i,j,k-5)-225d0*var(i,j,k-4)+400d0*var(i,j,k-3)-450d0*var(i,j,k-2)+360d0*var(i,j,k-1)-147d0*var(i,j,k))/(60d0*dz)
    else if(k==lastz-1) then
        d6dz = -(2d0*var(i,j,k-5)-15d0*var(i,j,k-4)+50d0*var(i,j,k-3)-100d0*var(i,j,k-2)+150d0*var(i,j,k-1)-77d0*var(i,j,k)-10d0*var(i,j,k+1))/(60d0*dz)
    else if(k==lastz-2) then
        d6dz = -(-var(i,j,k-4)+8d0*var(i,j,k-3)-30d0*var(i,j,k-2)+80d0*var(i,j,k-1)-35d0*var(i,j,k)-24d0*var(i,j,k+1)+2d0*var(i,j,k+2))/(60d0*dz)
    else
        d6dz = (var(i,j,k+3)-9d0*var(i,j,k+2)+45d0*var(i,j,k+1)-45d0*var(i,j,k-1)+9d0*var(i,j,k-2)-var(i,j,k-3))/(60d0*dz)
    end if

end function d6dz

!******************************************************************************************************************************************************************************************************


!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************










!******************************************************************************************************************************************************************************************************
!***********************************************************************************COMPUTES DERIVATIVES***********************************************************************************************
!******************************************************************************************************************************************************************************************************



!******************************************************************************FOR DERIVATIVES IN i DIRECTION******************************************************************************************

double precision function d6dx_vect(var,i,j,k,m)!differs from d6dx. just that var here is a 5D vector

    use variables
    implicit none
    integer,intent(in)                        :: i,j,k,m
    double precision,dimension(sx,sy,sz,neq),intent(in) :: var

    if(i==firstx) then
        d6dx_vect = (-10d0*var(i+6,j,k,m)+72d0*var(i+5,j,k,m)-225d0*var(i+4,j,k,m)+400d0*var(i+3,j,k,m)-450d0*var(i+2,j,k,m)+360d0*var(i+1,j,k,m)-147d0*var(i,j,k,m))/(60d0*dx)
    else if(i==firstx+1) then
        d6dx_vect = (2d0*var(i+5,j,k,m)-15d0*var(i+4,j,k,m)+50d0*var(i+3,j,k,m)-100d0*var(i+2,j,k,m)+150d0*var(i+1,j,k,m)-77d0*var(i,j,k,m)-10d0*var(i-1,j,k,m))/(60d0*dx)
    else if(i==firstx+2) then
        d6dx_vect = (-var(i+4,j,k,m)+8d0*var(i+3,j,k,m)-30d0*var(i+2,j,k,m)+80d0*var(i+1,j,k,m)-35d0*var(i,j,k,m)-24d0*var(i-1,j,k,m)+2d0*var(i-2,j,k,m))/(60d0*dx)
    else if(i==lastx) then
        d6dx_vect = -(-10d0*var(i-6,j,k,m)+72d0*var(i-5,j,k,m)-225d0*var(i-4,j,k,m)+400d0*var(i-3,j,k,m)-450d0*var(i-2,j,k,m)+360d0*var(i-1,j,k,m)-147d0*var(i,j,k,m))/(60d0*dx)
    else if(i==lastx-1) then
        d6dx_vect = -(2d0*var(i-5,j,k,m)-15d0*var(i-4,j,k,m)+50d0*var(i-3,j,k,m)-100d0*var(i-2,j,k,m)+150d0*var(i-1,j,k,m)-77d0*var(i,j,k,m)-10d0*var(i+1,j,k,m))/(60d0*dx)
    else if(i==lastx-2) then
        d6dx_vect = -(-var(i-4,j,k,m)+8d0*var(i-3,j,k,m)-30d0*var(i-2,j,k,m)+80d0*var(i-1,j,k,m)-35d0*var(i,j,k,m)-24d0*var(i+1,j,k,m)+2d0*var(i+2,j,k,m))/(60d0*dx)
    else
        d6dx_vect = (var(i+3,j,k,m)-9d0*var(i+2,j,k,m)+45d0*var(i+1,j,k,m)-45d0*var(i-1,j,k,m)+9d0*var(i-2,j,k,m)-var(i-3,j,k,m))/(60d0*dx)
    end if

end function d6dx_vect

!******************************************************************************************************************************************************************************************************





!******************************************************************************FOR DERIVATIVES IN j DIRECTION******************************************************************************************

double precision function d6dy_vect(var,i,j,k,m)!differs from d6dy. just that var here is a 5D vector

    use variables
    implicit none
    integer,intent(in)                        :: i,j,k,m
    double precision,dimension(sx,sy,sz,neq),intent(in) :: var

    if(j==firsty) then
        d6dy_vect = dy_dy_ph(i,j,k)*(-10d0*var(i,j+6,k,m)+72d0*var(i,j+5,k,m)-225d0*var(i,j+4,k,m)+400d0*var(i,j+3,k,m)-450d0*var(i,j+2,k,m)+360d0*var(i,j+1,k,m)-147d0*var(i,j,k,m))/(60d0*dy)
    else if(j==firsty+1) then
        d6dy_vect = dy_dy_ph(i,j,k)*(2d0*var(i,j+5,k,m)-15d0*var(i,j+4,k,m)+50d0*var(i,j+3,k,m)-100d0*var(i,j+2,k,m)+150d0*var(i,j+1,k,m)-77d0*var(i,j,k,m)-10d0*var(i,j-1,k,m))/(60d0*dy)
    else if(j==firsty+2) then
        d6dy_vect = dy_dy_ph(i,j,k)*(-var(i,j+4,k,m)+8d0*var(i,j+3,k,m)-30d0*var(i,j+2,k,m)+80d0*var(i,j+1,k,m)-35d0*var(i,j,k,m)-24d0*var(i,j-1,k,m)+2d0*var(i,j-2,k,m))/(60d0*dy)
    else if(j==lasty) then
        d6dy_vect = -dy_dy_ph(i,j,k)*(-10d0*var(i,j-6,k,m)+72d0*var(i,j-5,k,m)-225d0*var(i,j-4,k,m)+400d0*var(i,j-3,k,m)-450d0*var(i,j-2,k,m)+360d0*var(i,j-1,k,m)-147d0*var(i,j,k,m))/(60d0*dy)
    else if(j==lasty-1) then
        d6dy_vect = -dy_dy_ph(i,j,k)*(2d0*var(i,j-5,k,m)-15d0*var(i,j-4,k,m)+50d0*var(i,j-3,k,m)-100d0*var(i,j-2,k,m)+150d0*var(i,j-1,k,m)-77d0*var(i,j,k,m)-10d0*var(i,j+1,k,m))/(60d0*dy)
    else if(j==lasty-2) then
        d6dy_vect = -dy_dy_ph(i,j,k)*(-var(i,j-4,k,m)+8d0*var(i,j-3,k,m)-30d0*var(i,j-2,k,m)+80d0*var(i,j-1,k,m)-35d0*var(i,j,k,m)-24d0*var(i,j+1,k,m)+2d0*var(i,j+2,k,m))/(60d0*dy)
    else
        d6dy_vect = dy_dy_ph(i,j,k)*(var(i,j+3,k,m)-9d0*var(i,j+2,k,m)+45d0*var(i,j+1,k,m)-45d0*var(i,j-1,k,m)+9d0*var(i,j-2,k,m)-var(i,j-3,k,m))/(60d0*dy)
    end if

end function d6dy_vect

!******************************************************************************************************************************************************************************************************





!******************************************************************************FOR DERIVATIVES IN k DIRECTION******************************************************************************************

double precision function d6dz_vect(var,i,j,k,m)!differs from d6dz. just that var here is a 5D vector

    use variables
    implicit none
    integer,intent(in)                          :: i,j,k,m
    double precision,dimension(sx,sy,sz,neq),intent(in) :: var 

    if(k==firstz) then
        d6dz_vect = (-10d0*var(i,j,k+6,m)+72d0*var(i,j,k+5,m)-225d0*var(i,j,k+4,m)+400d0*var(i,j,k+3,m)-450d0*var(i,j,k+2,m)+360d0*var(i,j,k+1,m)-147d0*var(i,j,k,m))/(60d0*dz)
    else if(k==firstz+1) then
        d6dz_vect = (2d0*var(i,j,k+5,m)-15d0*var(i,j,k+4,m)+50d0*var(i,j,k+3,m)-100d0*var(i,j,k+2,m)+150d0*var(i,j,k+1,m)-77d0*var(i,j,k,m)-10d0*var(i,j,k-1,m))/(60d0*dz)
    else if(k==firstz+2) then
        d6dz_vect = (-var(i,j,k+4,m)+8d0*var(i,j,k+3,m)-30d0*var(i,j,k+2,m)+80d0*var(i,j,k+1,m)-35d0*var(i,j,k,m)-24d0*var(i,j,k-1,m)+2d0*var(i,j,k-2,m))/(60d0*dz)
    else if(k==lastz) then
        d6dz_vect = -(-10d0*var(i,j,k-6,m)+72d0*var(i,j,k-5,m)-225d0*var(i,j,k-4,m)+400d0*var(i,j,k-3,m)-450d0*var(i,j,k-2,m)+360d0*var(i,j,k-1,m)-147d0*var(i,j,k,m))/(60d0*dz)
    else if(k==lastz-1) then
        d6dz_vect = -(2d0*var(i,j,k-5,m)-15d0*var(i,j,k-4,m)+50d0*var(i,j,k-3,m)-100d0*var(i,j,k-2,m)+150d0*var(i,j,k-1,m)-77d0*var(i,j,k,m)-10d0*var(i,j,k+1,m))/(60d0*dz)
    else if(k==lastz-2) then
        d6dz_vect = -(-var(i,j,k-4,m)+8d0*var(i,j,k-3,m)-30d0*var(i,j,k-2,m)+80d0*var(i,j,k-1,m)-35d0*var(i,j,k,m)-24d0*var(i,j,k+1,m)+2d0*var(i,j,k+2,m))/(60d0*dz)
    else
        d6dz_vect = (var(i,j,k+3,m)-9d0*var(i,j,k+2,m)+45d0*var(i,j,k+1,m)-45d0*var(i,j,k-1,m)+9d0*var(i,j,k-2,m)-var(i,j,k-3,m))/(60d0*dz)
    end if

end function d6dz_vect

!******************************************************************************************************************************************************************************************************


!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************










!******************************************************************************************************************************************************************************************************
!*********************************************************************************POLYNOMIAL EXTRAPOLATION*********************************************************************************************
!******************************************************************************************************************************************************************************************************

double precision function fextrapolate(var,var_ind,i,s)
   
    implicit none
    integer, intent(in) :: s
    double precision,dimension(s),intent(in) :: var,var_ind
    integer,intent(in) :: i

    ! fextrapolate = 5d0*var(i+1) - 10d0*var(i+2) + 10d0*var(i+3) - 5d0*var(i+4) + var(i+5)
    ! fextrapolate = 4d0*var(i+1) - 6d0*var(i+2) + 4d0*var(i+3) - var(i+4)
    ! fextrapolate = 3d0*var(i+1) - 3d0*var(i+2) + var(i+3)
    fextrapolate = ((var_ind(i) - var_ind(i+2))/(var_ind(i+1) - var_ind(i+2)))*var(i+1) + ((var_ind(i) - var_ind(i+1))/(var_ind(i+2) - var_ind(i+1)))*var(i+2)
    ! fextrapolate = var(i+1)   

end function fextrapolate

!******************************************************************************************************************************************************************************************************





!******************************************************************************************************************************************************************************************************

double precision function bextrapolate(var,var_ind,i,s)

   implicit none
   integer, intent(in) :: s
   double precision,dimension(s),intent(in) :: var,var_ind
   integer,intent(in) :: i

   ! bextrapolate = 5d0*var(i-1) - 10d0*var(i-2) + 10d0*var(i-3) - 5d0*var(i-4) + var(i-5)
   ! bextrapolate = 4d0*var(i-1) - 6d0*var(i-2) + 4d0*var(i-3) - var(i-4)
   ! bextrapolate = 3d0*var(i-1) - 3d0*var(i-2) + var(i-3)
   bextrapolate = ((var_ind(i) - var_ind(i-2))/(var_ind(i-1) - var_ind(i-2)))*var(i-1) + ((var_ind(i) - var_ind(i-1))/(var_ind(i-2) - var_ind(i-1)))*var(i-2)
   ! bextrapolate = var(i-1)

end function bextrapolate

!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************
!******************************************************************************************************************************************************************************************************

