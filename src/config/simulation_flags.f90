module simulation_flags

  implicit none

  logical :: useSerialIO,     &
       useViscosity,          &
       useSplitConvection,    &
       useSplitViscosity,     &
       useUpwinding,          &
       useIBM,                &
       useShockCapturing,     &
       useTargetState,        &
       useConstantCFL,        &
       useGravity,            &
       useAdjoint,            &
       useContinuousAdjoint,  &
       useBaseline,           &
       useGradient,           &
       useParticles,          &
       usePTKE,               &
       useLES,                &
       useMassCorrection

  ! Old flags
  logical :: predictionOnly        = .true.,  &
       isBaselineAvailable         = .false., &
       isGradientAvailable         = .false., &
       twoWayCoupling              = .false., &
       boundScalars                = .false.

end module simulation_flags


! =============================== !
! Initialize the simulation flags !
! =============================== !
subroutine simulation_flags_setup

  ! Internal modules
  use simulation_flags

  ! External modules
  use parser

  implicit none

  call parser_read('use serial io', useSerialIO, .false.)
  call parser_read('include viscous terms', useViscosity, .false.)
  call parser_read('use adjoint solver', useAdjoint, .false.)
  call parser_read('disable adjoint solver', predictionOnly, .true.)
  call parser_read('use upwinding', useUpwinding, .false.)
  call parser_read('convective splitting', useSplitConvection, .not.useUpwinding)
  call parser_read('viscous splitting', useSplitViscosity, .false.)
  call parser_read('bound scalars', boundScalars, .false.)
  call parser_read('use immersed boundary', useIBM, .false.)
  call parser_read('use shock capturing', useShockCapturing, .false.)
  call parser_read('use target state', useTargetState, .true.)
  call parser_read('use constant CFL mode', useConstantCfl, .true.)
  call parser_read('baseline prediction available', isBaselineAvailable, .false.)
  call parser_read('adjoint gradient available', isGradientAvailable, .false.)
  call parser_read('use continuous adjoint', useContinuousAdjoint, .false.)
  call parser_read('include gravity', useGravity, .false.)
  call parser_read('include particles', useParticles, .false.)
  call parser_read('use ptke', usePTKE, .false.)
  call parser_read('use viscosity model', useLES, .false.)
  call parser_read('use mass correction', useMassCorrection, .false.)
  if (useParticles) then
     call parser_read('two way coupling', twoWayCoupling, .false.)
  else
     twoWayCoupling = .false.
  end if

  return
end subroutine simulation_flags_setup
