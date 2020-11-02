program center_probability

        Real :: d, cos_theta
        Real :: x_max, y_max, x_max_down, x_max_up, y_max_up, y_max_down, x_max_tail_up
        Real :: P_center, P_tail
        Real :: nO_tot
        Real :: nO_center, nO_tail
        Integer :: io, io2
        Character(100) :: max_x, max_y

        call get_command_argument(1,max_x)
        call get_command_argument(2,max_y)
        !Convert max_x max_y to real
        read(max_x,*) x_max
        read(max_y,*) y_max

        open(1,file='HBDONOR_network_tot_dL1')
        open(2,file='Probability_center')
        
        nO_tot = 0
        nO_tail = 0
        nO_center = 0
        x_max_down = x_max - 0.1
        x_max_up = x_max + 0.1
        y_max_down = y_max -0.2
        y_max_up = y_max + 0.2
        x_max_tail_up = x_max + 0.2

        do
            read(1,*,iostat = io2) d, cos_theta
            if (io2 /= 0) exit
            nO_tot = nO_tot + 1
                if (d >= x_max_down .and. d <= x_max_up .and. cos_theta >= y_max_down .and. cos_theta <= y_max_up) then
                    nO_center = nO_center + 1
                else if (d > x_max_tail_up .and. cos_theta >= -0.4 .and. cos_theta <= 0.4) then
                        nO_tail = nO_tail + 1
                end if
        end do
        print "(F8.0)", nO_tot
        print "(2F8.0)", nO_center, nO_tail
        P_center = nO_center/nO_tot
        P_tail = nO_tail/nO_tot
        print "(2F6.3)", P_center, P_tail
        close (1)
        write(2,FMT="(A,F8.0)") "Total number of Oxygen =" , nO_tot
        write(2,FMT="(A,F8.0)") "Number of Oxygen in the center region of 3D-plot = " , nO_center
        write(2,FMT="(A,F8.0)") "Number of Oxygen in tail of 3D-plot = " , nO_tail
        write(2,FMT="(A,F5.3)") "Probability of Oxygen in the center region of 3D-plot = " , P_center
        write(2,FMT="(A,F5.3)") "Probability of Oxygen in the tail region of 3D-plot = " , P_tail
        close (2)
end program
