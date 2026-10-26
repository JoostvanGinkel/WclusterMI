subroutine w1calc(data, scluster, n, p, nimp, ncluster, w1)

implicit none

integer :: n, p, nimp, ncluster
integer :: scluster(n)
double precision :: data(n, p, nimp)
double precision :: w1

integer :: c, m, i, j, k
integer :: ck
double precision :: wck, d2

w1 = 0.0d0

do c = 1, ncluster

   ck = 0
   do i = 1, n
      if (scluster(i) .eq. c) ck = ck + 1
   end do

   if (ck .gt. 1) then

      wck = 0.0d0

      do m = 1, nimp

         do i = 1, n-1
            if (scluster(i) .eq. c) then
               do j = i+1, n
                  if (scluster(j) .eq. c) then
                     d2 = 0.0d0
                     do k = 1, p
                        d2 = d2 + (data(i,k,m) - data(j,k,m))**2
                     end do
                     wck = wck + d2
                  end if
               end do
            end if
         end do

      end do

      wck = wck / (dble(nimp) * dble(ck))
      w1 = w1 + wck

   end if

end do

end subroutine w1calc

subroutine w_bar(data, centroids, n, p, m, c, cluster)

implicit none

integer, intent(in) :: n, p, m, c
double precision, intent(in) :: data(n,p,m)
double precision, intent(in) :: centroids(c,p,m)
integer, intent(out) :: cluster(n)

integer :: i, j, k, cl
double precision :: ssq, dist, min_dist
double precision :: ssquares(c)

do i = 1, n

    do cl = 1, c
        ssq = 0.0d0

        do k = 1, m
            do j = 1, p
                ssq = ssq + (data(i,j,k) - centroids(cl,j,k))**2
            end do
        end do

        ssquares(cl) = ssq / dble(m)
    end do

    min_dist = ssquares(1)
    cluster(i) = 1

    do cl = 2, c
        if (ssquares(cl) .lt. min_dist) then
            min_dist = ssquares(cl)
            cluster(i) = cl
        end if
    end do

end do

end subroutine w_bar

subroutine w_star(data, centroids, n, p, m, c, cluster)

implicit none

integer, intent(in) :: n, p, m, c
double precision, intent(in) :: data(n,p,m)
double precision, intent(in) :: centroids(c,p)
integer, intent(out) :: cluster(n)

integer :: i, j, k, cl
double precision :: ssq, dist, min_dist
double precision :: ssquares(c)

do i = 1, n

    do cl = 1, c
        ssq = 0.0d0

        do k = 1, m
            do j = 1, p
                ssq = ssq + (data(i,j,k) - centroids(cl,j))**2
            end do
        end do

        ssquares(cl) = ssq / dble(m)
    end do

    min_dist = ssquares(1)
    cluster(i) = 1

    do cl = 2, c
        if (ssquares(cl) .lt. min_dist) then
            min_dist = ssquares(cl)
            cluster(i) = cl
        end if
    end do

end do

end subroutine w_star