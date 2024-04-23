SUBROUTINE pentadiagonal_solver(x, a, b, c, d, e, N)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    REAL(KIND=8), INTENT(IN) :: a(N), b(N), c(N), d(N), e(N)
    REAL(KIND=8), INTENT(OUT) :: x(N)
    REAL(KIND=8) :: alpha(N), beta(N)
    INTEGER :: i

    ! Compute intermediate coefficients
    alpha(1) = b(1)
    beta(1) = c(1) / alpha(1)
    DO i = 2, N
        alpha(i) = b(i) - a(i) * beta(i-1) * e(i-1)
        beta(i) = c(i) / alpha(i)
    ENDDO

    ! Forward sweep
    x(1) = d(1) / alpha(1)
    DO i = 2, N
        x(i) = (d(i) - a(i) * x(i-1)) / alpha(i)
    ENDDO

    ! Backward substitution
    DO i = N-1, 1, -1
        x(i) = x(i) - beta(i) * x(i+1)
    ENDDO
END SUBROUTINE pentadiagonal_solver
