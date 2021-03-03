       subroutine rfftw_f77_create_plan()
       write(6,*) 'error: rfftw_f77_create_plan: FFTW is not included '
     &           ,'in this version.'
       stop 'NOFFTW'
       end

       subroutine rfftw_f77_one()
       write(6,*) 'error: rfftw_f77_one: FFTW is not included '
     &           ,'in this version.'
       stop 'NOFFTW'
       end

       subroutine rfftw_f77_destroy_plan()
       write(6,*) 'error: rfftw_f77_destroy_plan: FFTW is not included '
     &           ,'in this version.'
       stop 'NOFFTW'
       end
