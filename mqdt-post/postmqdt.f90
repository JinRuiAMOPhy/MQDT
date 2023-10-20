program main
   use type_def
   use constants
   use file_cmdl_io
   use interpolate
   use darray
   use search_root
   use envset
   implicit none
   type(mqdt) :: T
   type(Smat) :: S
   type(fname) :: FN
   type(control) :: CTR
   type(grid) :: G
   type(physobserv) :: phys
   type(spectrodata) :: spec
   call environmentset()
   call initiate(FN, T, S, CTR, G)
   !if(T%nop< T%nchan) then
     call findlev(FN, T, CTR, S, phys, spec)
   !end if
   call finalize(T, S, phys, spec)
end program
