! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Implements the device kernel for the `actuator_line_source_term_t` type.

module actuator_line_source_term_device
  use num_types, only : rp
  use field_list, only : field_list_t
  use field, only : field_t
  use device, only : device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use device_math, only : device_rzero, device_add2s2
  use registry, only : neko_registry
  implicit none
  private

  public :: actuator_line_source_term_compute_device

contains

  !> Computes the actuator line source term on the device.
  !! @param kernel Gaussian convolution kernel.
  !! @param kernel_j Field for holding kernel values and performing operations on the GPU.
  !! @param F_2D Idealized point source on each actuator line segment.
  !! @param fields The right-hand side, which should be the velocity components.
  subroutine actuator_line_source_term_compute_device(kernel, kernel_j, F_2D, fields)
    real(kind=rp), intent(in) :: kernel(:,:)
    type(field_t), intent(inout), pointer :: kernel_j
    real(kind=rp), intent(in) :: F_2D(size(kernel, 2), 3)
    type(field_list_t), intent(inout) :: fields

    integer :: j, n, n_alm
    type (field_t), pointer :: u, v, w, u_rhs, v_rhs, w_rhs
    real(kind=rp) :: kernel_temp(size(kernel, 1))

    n = fields%item_size(1)
    n_alm = size(kernel, 2)

    ! Get pointers to velocity fields
    u => neko_registry%get_field_by_name('u')
    v => neko_registry%get_field_by_name('v')
    w => neko_registry%get_field_by_name('w')

    ! Get pointers to momentum equations right hand sides
    u_rhs => fields%get_by_index(1)
    v_rhs => fields%get_by_index(2)
    w_rhs => fields%get_by_index(3)

    ! Clear old source terms
    call device_rzero(u_rhs%x_d, n)
    call device_rzero(v_rhs%x_d, n)
    call device_rzero(w_rhs%x_d, n)

    ! ! Update source terms
    do j = 1, n_alm
      kernel_temp = kernel(:,j)
      call device_memcpy(kernel_temp, kernel_j%x_d, n, HOST_TO_DEVICE, sync=.true.)
      call device_add2s2(u_rhs%x_d, kernel_j%x_d, -F_2D(j,1), n)
      call device_add2s2(v_rhs%x_d, kernel_j%x_d, -F_2D(j,2), n)
      call device_add2s2(w_rhs%x_d, kernel_j%x_d, -F_2D(j,3), n)
    end do

  end subroutine actuator_line_source_term_compute_device

end module actuator_line_source_term_device
