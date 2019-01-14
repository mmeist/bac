module circular_mesh
implicit none
private :: delaunay_condition, connect_prisms, wrap_idx, wrap_idx_inplace
public :: calc_mesh, calc_points, calc_n_tetras, calc_n_verts

contains

subroutine calc_points(verts_per_ring, n_slices, R0, Rmax, Z0,  points)
    double precision, parameter :: pi = 3.14159265358979d0

    double precision, intent(in) :: R0, Z0, Rmax
    integer, intent(in) :: n_slices
    integer, dimension(:), intent(in) :: verts_per_ring ! without venter vert; e.g. (/6, 8, 10/)
    double precision, dimension(:, :), intent(out) :: points !(r, phi, z)

    integer :: vert, ring, slice, n_rings, n_verts, verts_per_slice, vert_idx
    double precision :: r, theta, phi

    n_rings = size(verts_per_ring)
    n_verts = calc_n_verts(verts_per_ring, n_slices)
    verts_per_slice = n_verts / n_slices

    points(:, 1) = 0.d0 ! center vertex

    vert_idx = 2
    do ring = 1, n_rings
        r = (Rmax * ring) / n_rings
        do vert = 1, verts_per_ring(ring)
            theta = (2.d0 * pi * (vert - 1)) / verts_per_ring(ring)
            points(:, vert_idx) = r*(/cos(theta), 0.d0, sin(theta)/)
            vert_idx = vert_idx + 1
        end do
    end do

    points(1, :verts_per_slice) = points(1, :verts_per_slice) + R0
    points(3, :verts_per_slice) = points(3, :verts_per_slice) + Z0

    ! copy and rotate slice
    do slice = 2, n_slices
        vert_idx = (slice- 1) * verts_per_slice + 1
        phi = (2.d0 * pi * (slice - 1)) / n_slices
        points(:, vert_idx:vert_idx + verts_per_slice - 1) = points(:, 1:verts_per_slice)
        points(2, vert_idx:vert_idx + verts_per_slice - 1) = phi
    end do

end subroutine calc_points

subroutine calc_mesh(verts_per_ring, n_slices, points, verts, neighbours, neighbour_faces)
    integer, intent(in) :: n_slices
    integer, dimension(:), intent(in) :: verts_per_ring! without venter vert; e.g. (/6, 8, 10/)
    double precision, dimension(:, :), intent(in) :: points !points in first (r, phi, z) 
    integer, dimension(:, :), intent(out) :: verts, neighbours, neighbour_faces

    integer, dimension(4 , 6) :: tetra_conf, mask_theta, mask_phi, mask_r
    integer, dimension(4, 3) :: slice_offset, ring_offset, segment_offset
    integer, dimension(size(verts_per_ring)) :: prisms_per_ring
    integer, allocatable, dimension(:) :: bottom_prisms

    integer :: n_tetras, n_rings, n_verts, verts_per_slice, tetras_per_slice, &
               mask_idx, prev_ring_verts, prism_orientation, &
               slice, ring, segment, base_idx, prism_idx, tetra_idx, lower_off, upper_off, i

    double precision, dimension(2) :: u, v, p, q

    n_rings = size(verts_per_ring)
    
    if (any(verts_per_ring < 3) .or. n_rings < 1 .or. n_slices < 3) then
        print *, "Invalid parameters for function calc_mesh"
        stop
    end if
    
    verts_per_slice = sum(verts_per_ring) + 1 ! verts per slice
    n_verts = verts_per_slice * n_slices
    
    n_tetras = calc_n_tetras(verts_per_ring, n_slices)
    if (.not. size(verts, dim=2) == n_tetras .and. size(neighbours, dim=2) == n_tetras &
            .and. size(neighbour_faces, dim=2) == n_tetras) then
        print *, "Invalid size of output arrays for function calc_mesh"
        stop
    end if

    tetras_per_slice = n_tetras / n_slices
    
    prisms_per_ring = verts_per_ring
    prisms_per_ring(2:) = prisms_per_ring(2:) + verts_per_ring(:n_rings - 1)

    allocate(bottom_prisms(verts_per_slice - 1))

    print *, "verts_per_slice", verts_per_slice
    print *, "n_verts", verts_per_slice * n_slices
    print *, "tetras per slice", tetras_per_slice
    print *, "prisms per slice", tetras_per_slice / 3
    print *, "n tetras", n_tetras
    print *, "tetras_per_ring", prisms_per_ring * 3
    print *, "prisms_per_ring", prisms_per_ring * 3
    print *, "---"

    ! 
    !     7------6
    !    /|     /|
    !   3------2 |
    !   | 5----|-4           r     
    !   |/     |/            ^ ,phi
    !   1------0 *           |/    
    !               theta <--*

    ! top facing prism \/
    ! top edge from 2 - 7
    tetra_conf(:, :3) = reshape((/ &
        0, 2, 3, 7,  & !1 
        0, 2, 6, 7,  & !1
        0, 4, 6, 7  & !1
        /), (/4, 3/))
    
    ! calculate the bottom facing prism tetraeder configurations from the top facing one one
    ! bottom facing prism /\
    tetra_conf(:, 4:) = tetra_conf(:, :3)
    where (iand(tetra_conf(:, 4:) , 2) > 0)
        tetra_conf(:, 4:) = ieor(tetra_conf(:, 4:) , 1)
    end where
    where (iand(tetra_conf(:, 4:) , 1) > 0)
        tetra_conf(:, 4:) = ieor(tetra_conf(:, 4:) , 2)
    end where

    do i = 1, 6
        print *, tetra_conf(:, i)
        if (mod(i, 3) == 0) print *,
    end do

    ! extract masks from tetra_conf
    mask_theta = iand(tetra_conf, 1)
    mask_r = ishft(iand(tetra_conf, 2), -1)
    mask_phi = ishft(iand(tetra_conf, 4), -2)

    ! prefill neighbour_faces
    neighbour_faces(1, :) = 4
    neighbour_faces(2, :) = - 1
    neighbour_faces(3, :) = - 1
    neighbour_faces(4, :) = 1

    ! first slice
    prism_idx = 1
    tetra_idx = 1
    base_idx = 1
    do ring = 1, n_rings
        upper_off = 0
        lower_off = 0
        
        if (ring == 1) then
            prev_ring_verts = 1
        else
            prev_ring_verts = verts_per_ring(ring - 1)
        end if

        do segment = 1, prisms_per_ring(ring)

            ! --- calculate mask stuff ---
            ! u, v, p ,q = 0, 2, 3, 1

            if (verts_per_ring(ring) == prev_ring_verts) then
                prism_orientation = mod((segment - 1), 2)
            else
                u = points((/1, 3/), base_idx + modulo(lower_off, prev_ring_verts))
                v = points((/1, 3/), base_idx + prev_ring_verts + modulo(upper_off, verts_per_ring(ring)))
                p = points((/1, 3/), base_idx + prev_ring_verts + modulo(upper_off + 1, verts_per_ring(ring)))
                q = points((/1, 3/), base_idx + modulo(lower_off + 1, prev_ring_verts))

                if (ring == 1 .or. delaunay_condition(u, v, p, q, .true.)) then
                    ! triangle 0, 2, 3 satisfies delaunay condition => prism is up facing
                    prism_orientation  = 0
                else
                    ! prism is down facing
                    prism_orientation = 1
                end if
            end if
            
            mask_idx = 1 + prism_orientation * 3

            ! --- with masks calculate verts ---
            slice_offset = mask_phi(:, mask_idx:mask_idx + 2) * verts_per_slice
            ring_offset = mask_r(:, mask_idx:mask_idx + 2)
            segment_offset = mask_theta(:, mask_idx:mask_idx + 2)

            where (ring_offset /= 0) segment_offset = modulo(segment_offset + upper_off, verts_per_ring(ring))
            where (ring_offset == 0) segment_offset = modulo(segment_offset + lower_off, prev_ring_verts)
            ring_offset = ring_offset * prev_ring_verts

            verts(:, tetra_idx:tetra_idx + 2) = base_idx + slice_offset + ring_offset + segment_offset

            !print *, prism_idx, mask_idx, "|", verts(:, tetra_idx:tetra_idx + 2)

            ! --- connect ---
            ! connect neighbours tetras in this prism
            call connect_prisms(prism_idx, prism_idx, verts, neighbours, neighbour_faces)

            ! connect with prisms in neighbouring slices
            neighbours(1, tetra_idx) = tetra_idx + tetras_per_slice
            neighbours(4, tetra_idx + 2) = tetra_idx - tetras_per_slice

            ! connect with previous prism on ring
            if (.not. segment == 1) then
                call connect_prisms(prism_idx, prism_idx - 1, verts, neighbours, neighbour_faces)
            end if

            ! --- increment ---
            if (prism_orientation == 0) then ! top facing prism \/
                bottom_prisms(base_idx + (prev_ring_verts - 1) + upper_off) = prism_idx
                upper_off = upper_off + 1
            else ! bottom facing prism /\
                ! connect with prism on previous ring
                call connect_prisms(prism_idx, bottom_prisms(base_idx - 1 + lower_off), &
                                        verts, neighbours, neighbour_faces)
                lower_off = lower_off + 1 
            end if

            prism_idx = prism_idx + 1
            tetra_idx = (prism_idx - 1) * 3 + 1

        end do

        ! connect first and last prism in ring
        call connect_prisms(prism_idx - 1, prism_idx - prisms_per_ring(ring), verts, neighbours, neighbour_faces)

        base_idx = base_idx + prev_ring_verts

    end do

    deallocate(bottom_prisms)

    ! for all the other slices we can calculate the verts by incrimenting the
    ! vert indices of the first(idx = 0) slice by slice * verts_per_slice
    ! and the neighbours by incrementing the neighbours of the first slice by slice * tetras_per_slice
    ! neighbour_faces can just be copied
    do slice = 1, n_slices - 1 ! slice is the 0-based index of the slice
        tetra_idx = 1 + slice * tetras_per_slice

        verts(:, tetra_idx:tetra_idx + tetras_per_slice - 1) = verts(:, 1:tetras_per_slice) + slice * verts_per_slice

        ! wrap indices around in the last slice
        if (slice == n_slices - 1) then
            call wrap_idx_inplace(verts(:, tetra_idx:tetra_idx + tetras_per_slice - 1), n_verts)
        end if
    end do

    do slice = 1, n_slices - 1 ! slice is the 0-based index of the slice
        tetra_idx = 1 + slice * tetras_per_slice

        neighbours(:, tetra_idx:tetra_idx + tetras_per_slice - 1) = &
                neighbours(:, 1:tetras_per_slice) + slice * tetras_per_slice

        ! wrap indices around in the last slice
        if (slice == n_slices - 1) then
            call wrap_idx_inplace(neighbours(:, tetra_idx:tetra_idx + tetras_per_slice - 1), n_tetras)
        end if
    end do

    ! we dont wrap the neighbours index of the first slice earlier otherwise we would need to wrap the 
    ! index of every incrementally calculated slice
    call wrap_idx_inplace(neighbours(:, 1:tetras_per_slice), n_tetras)

    do slice = 1, n_slices - 1 ! slice is the 0-based index of the slice
        neighbour_faces(:, tetra_idx:tetra_idx + tetras_per_slice - 1) = neighbour_faces(:, 1:tetras_per_slice)
    end do
    
    ! remove outer faces after modulo operation
    where (neighbour_faces == -1)
        neighbours = -1
    end where

    print *, "finished"
    print *, "tetra_idx after finishing:", tetra_idx

end subroutine calc_mesh

pure function calc_n_verts(verts_per_ring, n_slices) result(n_verts)
    ! calculate the number of vertices in a torus given verts_per_ring and n_slices
    integer, dimension(:), intent(in) :: verts_per_ring
    integer, intent(in) :: n_slices
    integer :: n_verts

    n_verts = (sum(verts_per_ring) + 1) * n_slices
end function calc_n_verts

pure function calc_n_tetras(verts_per_ring, n_slices) result(n_tetras)
    ! calculate the number of tetraeders in a torus given verts_per_ring, n_rings and n_slices
    integer, dimension(:), intent(in) :: verts_per_ring
    integer, intent(in) :: n_slices
    integer :: n_tetras

    n_tetras = (sum(verts_per_ring) + sum(verts_per_ring(:size(verts_per_ring) - 1))) * 3 * n_slices
end function calc_n_tetras

pure function delaunay_condition(u, v, p, q, fixed_order) result(valid)
    ! check if the triangle [u, v, p] satisfies the delaunay condition with respect to point q
    ! returns: 
    !   valid: .true. the triangle a-b-c satisfies the delaunay condition 
    !          .false. otherwise

    double precision, dimension(2), intent(in) :: u, v, p, q
    logical, intent(in), optional :: fixed_order
    logical :: valid

    double precision :: a, b, c, d, e, f, g, h, i
    double precision :: delta, gamma

    a = u(1) - q(1)
    b = u(2) - q(2)
    c = (u(1) - q(1)) ** 2.d0 + (u(2) - q(2)) ** 2.d0
    d = v(1) - q(1)
    e = v(2) - q(2)
    f = (v(1) - q(1)) ** 2.d0 + (v(2) - q(2)) ** 2.d0
    g = p(1) - q(1)
    h = p(2) - q(2)
    i = (p(1) - q(1)) ** 2.d0 + (p(2) - q(2)) ** 2.d0

    delta = a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h

    if (.not.present(fixed_order) .or. fixed_order .eqv. .false.) then
        ! needed for result to be invariant to order of points
        gamma = p(1)*u(2) - p(1)*v(2) - p(2)*u(1) + p(2)*v(1) + u(1)*v(2) - u(2)*v(1)
    else 
        ! if the order of the points is counterclockwise we dont need to calculate gamma
        gamma = 1.d0
    end if

    valid = (delta * gamma < 0.0d0)

end function delaunay_condition

subroutine connect_prisms(prism_1_idx, prism_2_idx, verts, neighbours, neighbour_faces)
    integer, intent(in) :: prism_1_idx, prism_2_idx
    integer, dimension(:, :), intent(in) :: verts
    integer, dimension(:, :), intent(out) :: neighbours, neighbour_faces

    integer :: tetra_1_base, tetra_2_base, tetra_1_idx, tetra_2_idx, tetra_1_off, tetra_2_off, &
               face_idx_1, face_idx_2, i
    
    logical, dimension(4) :: same_vert_1, same_vert_2
    integer :: deleteme

    tetra_1_base = (prism_1_idx - 1) * 3 + 1
    tetra_2_base = (prism_2_idx - 1) * 3 + 1

    deleteme = 0
    do tetra_1_off = 0, 2
        tetra_1_idx = tetra_1_base + tetra_1_off
        do tetra_2_off = 0, 2
            tetra_2_idx = tetra_2_base + tetra_2_off

            if (tetra_1_idx == tetra_2_idx) cycle ! so we can use this method for inner connections of a prism

            do i = 1, 4
                same_vert_1(i) = any(verts(:, tetra_2_idx) == verts(i, tetra_1_idx))
                same_vert_2(i) = any(verts(i, tetra_2_idx) == verts(:, tetra_1_idx))
            end do

            if (count(same_vert_1) == 3) then
                face_idx_1 = minloc(transfer(same_vert_1 , 1, size=4), dim=1)
                face_idx_2 = minloc(transfer(same_vert_2 , 1, size=4), dim=1)

                neighbours(face_idx_1, tetra_1_idx) = tetra_2_idx
                neighbours(face_idx_2, tetra_2_idx) = tetra_1_idx

                neighbour_faces(face_idx_1, tetra_1_idx) = face_idx_2
                neighbour_faces(face_idx_2, tetra_2_idx) = face_idx_1

                deleteme = deleteme + 1
            end if
        end do
    end do

    if (deleteme < 2) then
        print *, "--- <2 match!!!", tetra_1_base, tetra_2_base, "prism", prism_1_idx, prism_2_idx, tetra_1_idx
        print *, "                ", verts(:, tetra_1_base:tetra_1_base + 2)
        print *, "                ", verts(:, tetra_2_base:tetra_2_base + 2)
        print *, "---"
        
    end if

end subroutine connect_prisms

elemental function wrap_idx(index, period) result(wrapped_index)
    ! wrap around the 1-based periodic index which may be larger or smaller than the period
    integer, intent(in) :: index
    integer, intent(in) :: period
    integer :: wrapped_index
    wrapped_index = modulo(index - 1, period) + 1
end function wrap_idx

elemental subroutine wrap_idx_inplace(index, period)
    integer, intent(inout) :: index
    integer, intent(in) :: period
    index = wrap_idx(index, period)
end subroutine wrap_idx_inplace

end module circular_mesh

program test
    use circular_mesh, only: calc_mesh, calc_points, calc_n_tetras, calc_n_verts
    implicit none
      
    integer, allocatable, dimension(:, :) :: v, n, nf
    double precision, allocatable, dimension(:, :) :: points
    integer :: n_tetras, n_points, i
    
    integer :: n_slices = 3
    !integer, dimension(3) :: vps = (/6, 9, 12/)
    integer, dimension(3) :: vps
    vps = (/(i, i=6, 6 + size(vps) * 2 - 1, 2)/)

    n_tetras = calc_n_tetras(vps, n_slices)
    n_points = calc_n_verts(vps, n_slices)

    allocate(v(4, n_tetras))
    allocate(n(4, n_tetras))
    allocate(nf(4, n_tetras))
    allocate(points(3, n_points))

    call calc_points(vps, n_slices, 171.d0, 96.d0, 0.d0, points)
    call calc_mesh(vps, n_slices, points(:, :n_points / n_slices), v, n, nf)
    !do i = 1, size(points, 2)
    !    print *, i, points(:, i)
    !end do
    !do i = 1, size(v, 2)
    !    print *, v(:, i)
    !    if (mod(i, 3) == 0) print *,
    !end do
end program test