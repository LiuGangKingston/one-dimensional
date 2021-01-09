! This is the FORTRAN 90 code for the one-dimensional simplified model system in the paper:
! Liu, G. A new equation for period vectors of crystals under external stress and temperature
! in statistical physics: mechanical equilibrium condition and equation of state.
! Eur. Phys. J. Plus 136, 48 (2021). https://doi.org/10.1140/epjp/s13360-020-01010-6
!
! In this code, almost all routine and variable names are self-explained, except temporary ones.
!
! Only one input file called in.dat is expected. If it runs successfully, many data files *.dat
! will be generated. For details about all the data files, please check the readme.md file.
!
! Please send all questions, comments, and other inquiries to Gang Liu: gl.cell@outlook.com .
! January 08th, 2021.




module data
    implicit none

    real*8,  parameter :: boltzmann_constant = 1.380649d-23
    real*8,  parameter :: absolute_zero = 1.0d-100
    real*8,  parameter :: neg_absolute_zero = - absolute_zero
    real*8,  parameter :: epsilon = 5.8000d-20
    real*8,  parameter :: lambda = 2.5000d0
    real*8,  parameter :: sigma = lambda
    real*8,  parameter :: coefficient = 4.0d0 * epsilon
    real*8,  parameter :: force_accuracy = 0.01d-10
    real*8,  parameter :: equtn_accuracy = 0.01d-10
    real*8,  parameter :: identified_useless_real = -0.7777777777d20
    integer, parameter :: identified_useless_integer = -7777777
    integer, parameter :: pow12 = 12, pow6 = 6
    integer, parameter :: pow6more  = pow6  * (pow6  + 1)
    integer, parameter :: pow12more = pow12 * (pow12 + 1)
    integer, parameter :: max_lines_in_a_figure = 26
    integer, parameter :: general_check_file_unit = 80
    integer, parameter :: one_time_input_file_unit = 11
    integer, parameter :: one_time_output_file_unit = 12
    integer, parameter :: must_cells = 10


    integer :: total_external_temperature_lines = identified_useless_integer
    real*8  :: external_temperature = identified_useless_real

    integer :: total_external_force_lines = identified_useless_integer
    real*8  :: external_force = identified_useless_real
    real*8  :: external_force_begin = identified_useless_real

    real*8  :: kinetic_internal_stress = identified_useless_real
    real*8  :: line_external_force_features(1:max_lines_in_a_figure) = identified_useless_real
    real*8  :: line_external_temperature_features(1:max_lines_in_a_figure) = identified_useless_real

    real*8  :: rigrous_zero_single_potential = identified_useless_real
    real*8  :: rigrous_zero_single_force = identified_useless_real
    real*8  :: rigrous_p_min_force_l_on_r =  identified_useless_real

    real*8  :: the_period_min_force_l_on_r = identified_useless_real
    real*8  :: the_min_force_l_on_r = identified_useless_real
    real*8  :: the_period = identified_useless_real
    real*8  :: force_of_left_on_right = identified_useless_real
    real*8  :: force_of_l_on_r_derivative = identified_useless_real
    real*8  :: the_cell_potential   = identified_useless_real
    integer :: cell_number = identified_useless_integer

end module data




module basics

    use data
    implicit none


contains

    function single_potential(r)
        implicit none
        real*8 :: r
        real*8 :: single_potential
        real*8 :: an_inter_v
        an_inter_v = sigma/r
        single_potential = coefficient * (an_inter_v**pow12 - an_inter_v**pow6)
        return
    end function single_potential


    subroutine single_potential_check()
        implicit none
        integer :: i
        real*8  :: r
        open(one_time_output_file_unit, file='single.potentials.dat')
        do i = 1, 1000
            r = 2.2 + 0.005 * i
            write(one_time_output_file_unit, *) r, single_potential(r), 0.0d0
        end do
        close(one_time_output_file_unit)
        return
    end subroutine single_potential_check


    subroutine calculate_cell_potential()
        implicit none
        real*8  :: acc, a_temp
        cell_number = 0
        the_cell_potential = 0.0d0
        potential_accumulated: do
            cell_number = cell_number + 1
            a_temp = the_cell_potential
            the_cell_potential   = the_cell_potential + single_potential(cell_number * the_period)

            if (cell_number .gt. must_cells) then
                if((abs(a_temp) .lt. absolute_zero) .and. (abs(the_cell_potential) .lt. absolute_zero)) then
                  print*,                          'Zero cell potential!', a_temp, the_cell_potential
                  write(general_check_file_unit,*) 'Zero cell potential!', a_temp, the_cell_potential
                  exit potential_accumulated
                end if
                acc = abs( (the_cell_potential - a_temp) / a_temp)
                if(acc .lt. force_accuracy) exit potential_accumulated
            end if
        end do potential_accumulated
        return
    end subroutine calculate_cell_potential


    function single_force(r)
        implicit none
        real*8 :: r
        real*8 :: single_force
        real*8 :: an_inter_f
        an_inter_f = sigma/r
        single_force = coefficient * (pow12 * an_inter_f**pow12 - pow6* an_inter_f**pow6) / r
        return
    end function single_force


    function force_unit_to_newton(f)
        implicit none
        real*8 :: force_unit_to_newton, f
        force_unit_to_newton = f * 1.0d10
        return
    end function force_unit_to_newton


    subroutine single_force_check()
        implicit none
        integer :: i
        real*8  :: r, sf
        open(one_time_output_file_unit, file='single.forces.dat')
        do i = 1, 1000
            r = 2.4 + 0.005 * i
            sf = single_force(r)
            write(one_time_output_file_unit, *) r, force_unit_to_newton(sf), 0.0d0
        end do
        close(one_time_output_file_unit)
        return
    end subroutine single_force_check


    subroutine calculate_force_of_left_on_right()
        implicit none
        real*8  :: acc, a_temp
        cell_number = 0
        force_of_left_on_right = 0.0d0
        force_accumulated: do
            cell_number = cell_number + 1
            a_temp = force_of_left_on_right
            force_of_left_on_right = force_of_left_on_right + cell_number * single_force(cell_number * the_period)

            if (cell_number .gt. must_cells) then
                if((abs(a_temp) .lt. absolute_zero) .and. (abs(force_of_left_on_right) .lt. absolute_zero))        then
                  print*,                          'Zero force from left on the right!', a_temp, force_of_left_on_right
                  write(general_check_file_unit,*) 'Zero force from left on the right!', a_temp, force_of_left_on_right
                  exit force_accumulated
                end if
                acc = abs( (force_of_left_on_right - a_temp) / a_temp)
                !write(general_check_file_unit,*) acc, force_of_left_on_right, a_temp
                if(acc .lt. force_accuracy) exit force_accumulated
            end if
        end do force_accumulated
        return
    end subroutine calculate_force_of_left_on_right


    subroutine c_force_of_l_on_r_derivative()
        implicit none
        real*8  :: acc, an_inter_f, a_temp
        cell_number = 0
        force_of_l_on_r_derivative = 0.0d0
        derivative_accumulated: do
            cell_number = cell_number + 1
            an_inter_f = sigma / (cell_number * the_period)
            a_temp = force_of_l_on_r_derivative
            force_of_l_on_r_derivative = force_of_l_on_r_derivative - &
                                       & coefficient * (pow12more * an_inter_f**pow12 - pow6more * an_inter_f**pow6) / &
                                       & (the_period * the_period)
            if (cell_number .gt. must_cells) then
                if((abs(a_temp) .lt. absolute_zero) .and. (abs(force_of_l_on_r_derivative) .lt. absolute_zero)) then
                  print*, 'Zero derivative force from left on the right!', a_temp, force_of_l_on_r_derivative
                  write(general_check_file_unit,*) 'Zero derivative force from left on the right!', &
                                                  & a_temp, force_of_l_on_r_derivative
                  exit derivative_accumulated
                end if
                acc = abs( (force_of_l_on_r_derivative - a_temp) / a_temp)
                !write(general_check_file_unit,*) acc, force_of_l_on_r_derivative, a_temp
                if(acc .lt. force_accuracy) exit derivative_accumulated
            end if
        end do derivative_accumulated
        return
    end subroutine c_force_of_l_on_r_derivative


    function expansion_rate_of_temperature()
        implicit none
        real*8 :: expansion_rate_of_temperature
        real*8 :: an_inter_f
        call c_force_of_l_on_r_derivative()
        an_inter_f = boltzmann_constant * external_temperature / (the_period * the_period)
        expansion_rate_of_temperature = boltzmann_constant / &
                                      & ( (an_inter_f - force_of_l_on_r_derivative) * the_period * the_period )
        return
    end function expansion_rate_of_temperature


    subroutine check_force_of_left_on_right()
        implicit none
        integer :: i
        open(one_time_output_file_unit, file='left.on.right.forces.dat')
        do i = 1, 1000
            the_period = 2.4 + 0.005 * i
            call calculate_force_of_left_on_right()
            write(one_time_output_file_unit, *) the_period, force_unit_to_newton(force_of_left_on_right),   &
                                                          & force_of_left_on_right, 0.0d0, cell_number
        end do
        close(one_time_output_file_unit)
        return
    end subroutine check_force_of_left_on_right


    subroutine save_find_the_min_l_on_r(loop)
        implicit none
        integer :: loop
        the_period_min_force_l_on_r = the_period
        call calculate_force_of_left_on_right()
        the_min_force_l_on_r = force_of_left_on_right
        print*,                           "The the_period_min_force_l_on_r: ", the_period_min_force_l_on_r
        write(general_check_file_unit, *) "The the_period_min_force_l_on_r: ", the_period_min_force_l_on_r
        write(general_check_file_unit, *) "The the_min_force_l_on_r: ", the_min_force_l_on_r, " with loops of ", loop
        return
    end subroutine save_find_the_min_l_on_r


    subroutine find_the_min_l_on_r()
        implicit none
        integer :: loop
        logical :: r1positive
        real*8  :: r1, r2, rz, t1, t2

        t1 = 0.01d0
        the_period = t1
        call c_force_of_l_on_r_derivative()
        r1 =   force_of_l_on_r_derivative
        if(r1 .gt. absolute_zero) then
            r1positive = .true.
        else if(r1 .lt. neg_absolute_zero) then
            r1positive = .false.
        else
            call save_find_the_min_l_on_r(0)
            return
        end if

        t2 = 200.d0
        the_period = t2
        call c_force_of_l_on_r_derivative()
        r2 =   force_of_l_on_r_derivative
        if(abs(r2) .lt. absolute_zero) then
            call save_find_the_min_l_on_r(0)
            return
        end if
        if(r1positive) then
            if(r2 .gt. absolute_zero) then
                print*, "Sorry in find_the_min_l_on_r subroutine all positive: ", r1, r2
                stop
            end if
        else
            if(r2 .lt. neg_absolute_zero) then
                print*, "Sorry in find_the_min_l_on_r subroutine all negtive: ", r1, r2
                stop
            end if
        end if

        loop = 0
        finding_loop: do
            loop = loop + 1
            the_period = (t1 + t2)/2
            call c_force_of_l_on_r_derivative()
            rz =   force_of_l_on_r_derivative
            if(abs(rz)    .lt.  absolute_zero)  exit finding_loop
            if((t2-t1)/t2 .lt. force_accuracy)  exit finding_loop
            if(r1positive) then
                 if(rz .gt. absolute_zero) then
                       t1 = the_period
                       r1 = rz
                 else
                       t2 = the_period
                       r2 = rz
                 end if
            else
                 if(rz .gt. absolute_zero) then
                       t2 = the_period
                       r2 = rz
                 else
                       t1 = the_period
                       r1 = rz
                 end if
            end if
        end do finding_loop

        call save_find_the_min_l_on_r(loop)

        return
    end subroutine find_the_min_l_on_r


    subroutine rigrous_mininum_points()
        implicit none
        integer :: i
        real*8  :: r1, r2
        rigrous_zero_single_potential = sigma
        rigrous_zero_single_force = sigma * (1.0d0 * pow12 / pow6) ** (1.0d0 / (pow12 - pow6))
        write(general_check_file_unit, *)
        write(general_check_file_unit, *)
        write(general_check_file_unit, *) "Rigorous zero potential, force, and minimum force left on right: "
        write(general_check_file_unit, *) "Rig zero single potential: ", rigrous_zero_single_potential
        write(general_check_file_unit, *) "Rig zero single force: ", rigrous_zero_single_force
        r1 = 0.0d0
        r2 = 0.0d0
        do i = 1, 8
            r1 = r1 + 1.0d0/i**pow12
            r2 = r2 + 1.0d0/i**pow6
            rigrous_p_min_force_l_on_r = sigma * (r1 * pow12more/ (r2*pow6more)) ** (1.0d0 / (pow12 - pow6))
            write(general_check_file_unit, *) "Rig period for min force from L on R: ", rigrous_p_min_force_l_on_r, "i=", i
        end do
        write(general_check_file_unit, *) "Cal period for min force from L on R: ", the_period_min_force_l_on_r
        write(general_check_file_unit, *)
        return
    end subroutine rigrous_mininum_points


    function the_equation()
        implicit none
        real*8  :: the_equation
        kinetic_internal_stress = boltzmann_constant * external_temperature / the_period
        call calculate_force_of_left_on_right()
        the_equation = force_of_left_on_right + kinetic_internal_stress + external_force
        return
    end function the_equation


    function the_kinetic_energy()
        implicit none
        real*8  :: the_kinetic_energy
        the_kinetic_energy = 3.0d0 * boltzmann_constant * external_temperature / 2.0d0
        return
    end function the_kinetic_energy


    function get_the_external_force()
        implicit none
        real*8  :: get_the_external_force
        kinetic_internal_stress = boltzmann_constant * external_temperature / the_period
        call calculate_force_of_left_on_right()
        get_the_external_force = -1.0d0 * (force_of_left_on_right + kinetic_internal_stress)
        return
    end function get_the_external_force


    function get_the_external_temperature()
        implicit none
        real*8  :: get_the_external_temperature
        call calculate_force_of_left_on_right()
        get_the_external_temperature = -1.0d0 * (force_of_left_on_right + external_force) * &
                                                the_period / boltzmann_constant
        return
    end function get_the_external_temperature


    subroutine find_the_period()
        implicit none
        integer :: loop
        real*8  :: t1, r1, t2, r2, rz
        t1 = 1.0d0
        loop = 0
        loop1: do
            loop = loop + 1
            if(loop .gt. 100) then
               print*, "Sorry in loop1 of find_the_period subroutine: ", loop, the_period, r1
               stop
            end if
            t1 = t1*0.1d0
            the_period = t1
            r1 = the_equation()
            if(r1 .gt. absolute_zero) exit loop1
        end do loop1

        t2 = the_period_min_force_l_on_r
        the_period = t2
        r2 = the_equation()
        if(r2 .gt. absolute_zero) then
            print*, "Sorry for t2 in find_the_period subroutine: ", t1, t2, r1, r2, &
                    & force_of_left_on_right, kinetic_internal_stress, external_force
            write(general_check_file_unit, *) "Sorry for t2 in find_the_period subroutine: ", t1, t2, r1, r2, &
                    & force_of_left_on_right, kinetic_internal_stress, external_force
            stop
        else if(r2 .lt. neg_absolute_zero) then
            loop = 0
            finding_loop: do
                loop = loop + 1
                the_period = (t1 + t2)/2
                rz = the_equation()
                !write(general_check_file_unit, *) loop, the_period, rz, t1, r1
                if(rz .gt. absolute_zero) then
                   t1 = the_period
                   r1 = rz
                else
                   t2 = the_period
                   r2 = rz
                end if
                if(abs(rz)    .lt.  absolute_zero)  exit finding_loop
                if((t2-t1)/t2 .lt. equtn_accuracy)  exit finding_loop
            end do finding_loop
        end if

        call calculate_force_of_left_on_right()
        return
    end subroutine find_the_period


    subroutine initial()
        implicit none
        integer :: i
        open(general_check_file_unit, file='general.checking.dat')
        write(general_check_file_unit, *) "For checking purpose ..."
        print*, "Beginning ..."
        open(one_time_input_file_unit, file='in.dat')

        read(one_time_input_file_unit,*) total_external_force_lines
        if(total_external_force_lines .gt. 26) then
           print*, "Too many lines of period vs force required in the in.dat: ", total_external_force_lines
           stop
        else if(total_external_force_lines .gt. 0) then
           read(one_time_input_file_unit,*) (line_external_temperature_features(i), i=1,total_external_force_lines)
           read(one_time_input_file_unit,*) external_force_begin
        end if

        read(one_time_input_file_unit,*) total_external_temperature_lines
        if(total_external_temperature_lines .gt. 26) then
           print*, "Too many lines of period vs temperature required in the in.dat: ", total_external_temperature_lines
           stop
        else if(total_external_temperature_lines .gt. 0) then
           read(one_time_input_file_unit,*) (line_external_force_features(i), i=1,total_external_temperature_lines)
        end if

        close(one_time_input_file_unit)

        return
    end subroutine initial


    subroutine finalize()
        implicit none
        print*, "Done."
        write(general_check_file_unit, *) "Finalizing done."
        close(general_check_file_unit)
        return
    end subroutine finalize


end module basics




subroutine calculate_period_vs_force()
    use basics
    implicit none
    integer :: line_i, i, total_points=20000
    real*8  :: step
    do line_i = 1, total_external_force_lines
          write(general_check_file_unit, *) "Working in calculate_period_vs_force: ", line_i
          open(one_time_output_file_unit, file="period.vs.force."//char(64+line_i)//".dat")
          external_temperature = line_external_temperature_features(line_i)

          step = (the_min_force_l_on_r * (-1.1d0) - external_force_begin) / (total_points - 1)
          external_force_loop: do i = 1, total_points
             external_force = external_force_begin + step * (i - 1)
             the_period = the_period_min_force_l_on_r
             if(the_equation() .gt. absolute_zero) cycle external_force_loop
             call find_the_period()
             write(one_time_output_file_unit,*) external_temperature, force_unit_to_newton(external_force), &
                                              & the_period
          end do external_force_loop

          write(one_time_output_file_unit,*)
          close(one_time_output_file_unit)
    end do
    write(general_check_file_unit, *) "Done in calculate_period_vs_force. "
    write(general_check_file_unit, *)
    return
end subroutine calculate_period_vs_force




subroutine hookes_elastics()
    use basics
    implicit none
    integer :: line_i
    real*8  :: lower_force = 0.0d0, upper_force = 1.0d-20, a1, a2, h
    open(one_time_output_file_unit, file="hookes.elastics.dat")
    if(total_external_force_lines .gt. 0) then
       do line_i = 1, total_external_force_lines
          external_temperature = line_external_temperature_features(line_i)
          print*, "Hooke's elastics: ", line_i, external_temperature
          external_force = lower_force
          call find_the_period()
          a1 = the_period
          external_force = upper_force
          the_period = the_period_min_force_l_on_r
          if(the_equation() .gt. absolute_zero) cycle
          call find_the_period()
          a2 = the_period
          h = (force_unit_to_newton(upper_force) - force_unit_to_newton(lower_force)) / ((a2-a1)/a1)
          write(one_time_output_file_unit,*) external_temperature, force_unit_to_newton(lower_force), &
                                              & force_unit_to_newton(upper_force), a1, a2, a2-a1, h
          print*, "Hooke's elastics result: " , a1, a2, a2-a1, h
       end do
    end if
    write(one_time_output_file_unit,*)
    close(one_time_output_file_unit)
    return
end subroutine hookes_elastics




subroutine calculate_work_and_heat()
    use basics
    implicit none
    integer :: line_i, i, total_points=10000
    real*8  :: lower_force = 0.0d0, step, a1, work, heat, kinetic, energy, energy_0
    do line_i = 1, total_external_force_lines
          write(general_check_file_unit, *) "Working in calculate_work_and_heat: ", line_i
          open(one_time_output_file_unit, file="work.and.heat."//char(64+line_i)//".dat")
          external_temperature = line_external_temperature_features(line_i)

          external_force = lower_force
          call find_the_period()
          a1 = the_period
          call calculate_cell_potential()
          energy_0 = the_cell_potential + the_kinetic_energy()
          step = (the_period_min_force_l_on_r - a1) / (total_points - 1)
          work = 0.0d0

          external_force_loop: do i = 1, total_points
             the_period = a1 + step * (i-1)
             external_force = get_the_external_force()
             !if(external_force .lt. absolute_zero) exit external_force_loop
             work = work + external_force * step
             call calculate_cell_potential()
             kinetic = the_kinetic_energy()
             energy = the_cell_potential + kinetic
             heat = energy - energy_0 - work
             write(one_time_output_file_unit,*) external_temperature, force_unit_to_newton(external_force), &
                                              & the_period, work, the_cell_potential, kinetic, energy, heat
          end do external_force_loop
          write(one_time_output_file_unit,*)
          close(one_time_output_file_unit)
          print*, "Work and heat: " , line_i, work, the_cell_potential, kinetic, energy, heat
    end do
    write(general_check_file_unit, *) "Done in calculate_work_and_heat."
    write(general_check_file_unit, *)
    return
end subroutine calculate_work_and_heat




subroutine calculate_period_vs_temperature()
    use basics
    implicit none
    integer :: line_i, i, total_points = 1000
    real*8  :: ini_t=0.0d0, step, end_t, a1
    print*, "total_external_temperature_lines: ", total_external_temperature_lines
    if(total_external_temperature_lines .gt. 0) then
       do line_i = 1, total_external_temperature_lines
          open(one_time_output_file_unit, file="period.vs.temperature."//char(64+line_i)//".dat")
          external_force = line_external_force_features(line_i)
          write(general_check_file_unit, *) "Working in calculate_period_vs_temperature: ", line_i

          external_temperature = ini_t
          call find_the_period()
          a1 = the_period
          the_period = the_period_min_force_l_on_r
          end_t = get_the_external_temperature() * 0.98d0
          end_t = get_the_external_temperature() * 0.999999999999d0
          if(end_t .lt. absolute_zero) cycle
          step = (end_t - ini_t) / (total_points - 1)
          do i = 1, total_points
             external_temperature = ini_t + step * (i - 1)
             call find_the_period()
             write(one_time_output_file_unit,*) external_temperature, force_unit_to_newton(external_force), &
                                              & the_period, (the_period - a1)/(step*a1), &
                                              & expansion_rate_of_temperature()
             a1 = the_period
          end do

          close(one_time_output_file_unit)
       end do
    end if
    write(general_check_file_unit, *) "Done in calculate_period_vs_temperature: "
    write(general_check_file_unit, *)
    return
end subroutine calculate_period_vs_temperature




subroutine work_and_heat_for_temperature()
    use basics
    implicit none
    integer :: line_i, i, total_points=10000
    real*8  :: ini_t=0.0d0, step, end_t, a1, work, heat, kinetic, energy, energy_0, p_energy_0, p_heat
    print*, "total_external_temperature_lines: ", total_external_temperature_lines
    do line_i = 1, total_external_temperature_lines
          open(one_time_output_file_unit, file="work.for.temperature."//char(64+line_i)//".dat")
          external_force = line_external_force_features(line_i)
          write(general_check_file_unit, *) "Working in work_and_heat_for_temperature: ", line_i

          external_temperature = ini_t
          call find_the_period()
          a1 = the_period
          call calculate_cell_potential()
          p_energy_0 = the_cell_potential
          energy_0   = the_cell_potential + the_kinetic_energy()
          step = (the_period_min_force_l_on_r  - a1) / (total_points - 1)
          work = 0.0d0

          do i = 1, total_points
             the_period = a1 + step * (i-1)
             external_temperature = get_the_external_temperature()
             work = work + external_force * step
             call calculate_cell_potential()
             kinetic = the_kinetic_energy()
             energy = the_cell_potential + kinetic
             heat = energy - energy_0 - work
             p_heat = the_cell_potential - p_energy_0 - work
             write(one_time_output_file_unit,*) external_temperature, force_unit_to_newton(external_force), &
                                              & the_period, work, the_cell_potential, kinetic, energy, heat, p_heat
          end do

          write(one_time_output_file_unit,*)
          close(one_time_output_file_unit)
          print*, "Work and heat for temperature: " , line_i, work, the_cell_potential, kinetic, energy, heat
    end do
    write(general_check_file_unit, *) "Done in work_and_heat_for_temperature."
    write(general_check_file_unit, *)
    return
end subroutine work_and_heat_for_temperature




subroutine checks_and_special_points()
    use basics
    implicit none
    call single_potential_check()
    call single_force_check()
    call check_force_of_left_on_right()
    call find_the_min_l_on_r()
    call rigrous_mininum_points()
    write(general_check_file_unit, *) "Done in checks_and_special_points."
    write(general_check_file_unit, *)
    return
end subroutine checks_and_special_points




program one_dim
    use basics
    implicit none
    call initial()
    call checks_and_special_points()
    call calculate_period_vs_force()
    call hookes_elastics()
    call calculate_work_and_heat()
    call calculate_period_vs_temperature()
    call work_and_heat_for_temperature()
    call finalize()
    stop
end program one_dim



