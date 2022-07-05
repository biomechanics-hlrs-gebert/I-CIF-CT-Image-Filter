REAL(rk)      , DIMENSION(:,:,:), INTENT(IN) :: kernel
INTEGER(ik)   , DIMENSION(6),     INTENT(IN) :: srb ! subarray_reduced_bndaries

INTEGER(ik) :: accumulator, ii, jj, kk, ll, mm, nn, border
INTEGER(ik), DIMENSION(3) :: kernel_size

kernel_size = SIZE(kernel, DIM=1)

border = (kernel_size(1)-1)/2

DO kk = srb(3), srb(6)
DO jj = srb(2), srb(5)
DO ii = srb(1), srb(4)
    accumulator = 0
    DO ll = -border, border
    DO mm = -border, border
    DO nn = -border, border
        accumulator = accumulator + (kernel( ll+border+1_ik, mm+border+1_ik, nn+border+1_ik) * &
                subarray(ii + ll, jj + mm, kk + nn))
    END DO
    END DO
    END DO
    result_subarray(ii - border, jj - border, kk - border) = accumulator
END DO
END DO
END DO
