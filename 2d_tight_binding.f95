program main
    use iso_fortran_env
    implicit none

    INTEGER(KIND=8) :: matrix_nn(100_8,100_8), p, x1,y1,x2,y2,x3,y3,x4,y4
    INTEGER(KIND=8) :: matrix_nnn(100_8,100_8)
    !INTEGER(kind=8)  :: state(10), val, k, hh_list(4), 

    

    matrix_nn = nn_matrix(10_8,10_8,1_8)
    write(1,'(100(",",G0,:))') matrix_nn
    print "(100i2)", matrix_nn
    print *, "---------"
    matrix_nnn = nnn_matrix(10_8,10_8,1_8)
    write(2,'(100(",",G0,:))') matrix_nnn
    print "(100i2)", matrix_nnn
    !call nn_hop_list(5_8,2_8,0_8,hh_list)
    !print *, hh_list
    !call Inv_Map(9_8,2_8,p,q)
    !call pbc_nn(p,5_8,q,2_8)
    !PRINT *, p, q
    !call nnn(3_8,2_8,4_8,3_8,x1,y1,x2,y2,x3,y3,x4,y4)
    










    contains


    !THE SUBROUTINES...



    !Swapping 2 elements (of given index i, j) of an array.
    subroutine swap(arr, i, j)
        implicit none
        INTEGER(kind=int64), INTENT(IN) :: i, j
        INTEGER(kind=int64), INTENT(INOUT) :: arr(:)
        INTEGER(kind=int64) :: temp
        temp = arr(i)
        arr(i) = arr(j)
        arr(j) = temp
        
    end subroutine swap


    !MAP 2d-->1d for i={0,...,m-1}, j={0,...,n-1} gives val={0,...,m*n-1}
    function Map(i, j, n)
        implicit none
        INTEGER(kind=int64) :: Map
        INTEGER(kind=int64) :: i, j, n
        Map = i*n + j
        RETURN
    end function

    !Inverse MAP 1d-->2d for val={0,...,m*n-1} gives i={0,...,m-1}, j={0,...,n-1}
    subroutine Inv_Map(val, n, i, j)
        INTEGER(kind=int64), INTENT(IN) :: val, n
        INTEGER(kind=int64), INTENT(OUT) :: i, j
        j = MOD(val, n)
        i = (val - j)/n
    end subroutine Inv_Map
    
    
    !reverse of a 1d array
    !USE WITH CAUTION. HAS SIDEFFECTS (MUTATES INPUT ARRAY)
    subroutine reverse(arr, len)
        INTEGER(kind=int64), INTENT(IN) :: len
        INTEGER(kind=int64) :: arr(len), reverse_arr(len), i
        do i = 1, len
            reverse_arr(i) = arr(len-i+1)
        end do
        arr = reverse_arr
    end subroutine

    !FACTORIAL FUNCTION
    INTEGER(kind=int64) function fact(n1)
        INTEGER(kind=int64) :: n1, fact1, i1

        if (n1 == 0) then
            fact = 1
            RETURN
        end if
    
        fact1=1
        do i1=1, n1
            fact1=fact1*i1
        end do
    
        fact=fact1
        return
    end function
    !end of factorial function
    
    !COMBINATION FUNCTION
    INTEGER(kind=int64) function comb(n,r)
        INTEGER(kind=int64) :: n , r
        if (n .eq. 0) then
            comb = 0
            return
        end if
        comb = (fact(n))/(fact(r)*fact(n-r))
        return 
    end function
    !end of combbination function

    !Function that converts POB representation of state to a unique INTEGER(kind=int64) "pob2integer" = 1,...,comb(mn,mn/2)
    !NOT USED IN THIS PROGRAM.
    INTEGER(kind=int64) function pob2number(arr, len)
        INTEGER(kind=int64) :: len, temp, i, j
        INTEGER(kind=int64) :: arr(0:len-1)
        INTEGER(kind=int64) :: temp_arr(0:len-1)
        temp_arr = arr
        call reverse(temp_arr, len)
        pob2number = 0
        do i = 0, len-1
            temp = 0
            do j = 0, i
                temp = temp + temp_arr(j)
            end do
            pob2number = pob2number + temp_arr(i)*comb(i, temp)
        end do
        RETURN
    end function


    !Function that converts unique INTEGER(kind=int64) value v to pob representation as a 1d array of size=n, with r 1s...
    function number2pob(n, r, v) result(arr)
        implicit none

        INTEGER(kind=int64), intent (in) :: n
        INTEGER(kind=int64), intent (in) :: r
        INTEGER(kind=int64), intent (in) :: v

        INTEGER(kind=int64), dimension(n) :: arr

        INTEGER(kind=int64) :: j
        INTEGER(kind=int64) :: temp, temp1
        INTEGER(kind=int64) :: p

        INTEGER(kind=int64) :: k

        j = n
        temp = v

        do k = r, 1, -1
            do while (.true.)
                j = j - 1
                p = comb(j,k)
                if (temp >= p) then
                    temp = temp - p
                    arr(j+1) = 1 !.true.
                    
                else
                    arr(j+1) = 0 !.false.
                end if
                if (arr(j+1) == 1) then
                    exit
                end if
            end do
        end do
        
        if (j>0) then
            do k = j, 1, -1
                arr(k) = 0
            end do
        end if
        temp1 = SIZE(arr)
        CALL reverse(arr, temp1)
        RETURN
    end function number2pob










    !given m,n and i,j returns coordinates enforcing pbc hops...
    subroutine pbc_nn(i,m,j,n)
        implicit none
        INTEGER(kind=int64), INTENT(IN) :: m, n
        INTEGER(kind=int64), INTENT(INOUT) :: i, j
        i = MOD((MOD(i,m) + m),m)
        j = MOD((MOD(j,n) + n),n)
    end subroutine
        
    !given m,n and i,j gives its NN coordinates with pbc
    subroutine nn(i,j,m,n,x1,y1,x2,y2,x3,y3,x4,y4)
        implicit none
        INTEGER(kind=int64), INTENT(IN) :: i, j, m, n
        INTEGER(kind=int64), INTENT(INOUT) :: x1,y1,x2,y2,x3,y3,x4,y4
        x1 = i
        y1 = j+1
        call pbc_nn(x1,m,y1,n)
        !print "(i0,a1,i0)", x1,",",y1
        x2 = i+1
        y2 = j
        call pbc_nn(x2,m,y2,n)
        !print "(i0,a1,i0)", x2,",",y2
        x3 = i
        y3 = j-1
        call pbc_nn(x3,m,y3,n)
        !print "(i0,a1,i0)", x3,",",y3
        x4 = i-1
        y4 = j
        call pbc_nn(x4,m,y4,n)
        !print "(i0,a1,i0)", x4,",",y4
    end subroutine

    !pob_val = {0,...,m*n-1}
    subroutine nn_hop_list(m,n,val,h_list)
        implicit none
        INTEGER(kind=int64), INTENT(IN) :: m,n,val
        INTEGER(kind=int64), INTENT(OUT) :: h_list(4)
        INTEGER(kind=int64) :: state(m*n), temp_state(m*n), k,a,b,c,d,loc
        INTEGER(kind=int64) :: x1,y1,x2,y2,x3,y3,x4,y4, i, j
        state = 0
        state(val+1) = 1
        !PRINT "(10i0)", state
        !h_list = -1
        !count = 1
        k = val
        call Inv_Map(k,n,i,j)
        !PRINT "(3i0)", i, j, k
        call nn(i,j,m,n,x1,y1,x2,y2,x3,y3,x4,y4)
        !print "(8i2)",x1,y1,x2,y2,x3,y3,x4,y4
        a = Map(x1,y1,n)+1
        b = Map(x2,y2,n)+1
        c = Map(x3,y3,n)+1
        d = Map(x4,y4,n)+1

       
        h_list(1) = a

        h_list(2) = b

        h_list(3) = c

        h_list(4) = d
        
    end subroutine


    !generates hamiltonian matrix of NN hopping for given m,n of the lattice
    function nn_matrix(m,n, t) result(arr)
        implicit none
        INTEGER(kind=int64), INTENT(IN) :: m,n, t
        INTEGER(kind=int64) :: dim, i, j, temp_list(4), k
        REAL(kind=real64), ALLOCATABLE :: arr(:,:)
        dim = m*n
        ALLOCATE (arr(dim,dim))
        arr = 0 
        do i  = 1, dim
            CALL nn_hop_list(m,n,i-1,temp_list)
            !print "(4i0)", temp_list
            do j = 1, 4
                arr(i,temp_list(j)) = t
            end do
            
        end do

        arr = 1_8*arr

        do i = 1, dim
            arr(i,i) = 0
        end do
        return

    end function














    !given m,n and i,j gives its NNN coordinates with pbc
    subroutine nnn(i,j,m,n,x1,y1,x2,y2,x3,y3,x4,y4)
        implicit none
        INTEGER(kind=int64), INTENT(IN) :: i, j, m, n
        INTEGER(kind=int64), INTENT(INOUT) :: x1,y1,x2,y2,x3,y3,x4,y4
        INTEGER(kind=int64) :: temp
        x1 = i+1
        y1 = j+1
        !print "(i0,a1,i0)", x1,",",y1
        if (x1>(m-1) .or. y1>(n-1)) then
            do while (x1.ne.0 .and. y1.ne.0)
                x1 = x1 -1
                y1 = y1 -1
                
            end do
        end if


        !print "(i0,a1,i0)", x1,",",y1

        x2 = i + 1
        y2 = j - 1
        if (x2 > m-1 .or. y2 < 0) then
            do while (x2.ne.0 .and. y2.ne.(n-1))
                x2 = x2 -1
                y2 = y2 +1
            end do
        end if
        
        !print "(i0,a1,i0)", x2,",",y2

        x3 = i-1
        y3 = j-1
        if (x3 < 0 .or. y3 < 0) then
            do while (x3.ne.(m-1) .and. y3.ne.(n-1))
                x3 = x3 +1
                y3 = y3 +1
            end do
        end if
        !print "(i0,a1,i0)", x3,",",y3

        x4 = i-1
        y4 = j+1
        if (x4 < 0 .or. y4 > n-1) then
            do while (x4.ne.(m-1) .and. y4.ne.0)
                x4 = x4 +1
                y4 = y4 -1
            end do
        end if
           

        !print "(i0,a1,i0)", x4,",",y4
    end subroutine

    !pob_val = {0,...,m*n-1}
    subroutine nnn_hop_list(m,n,val,h_list)
        implicit none
        INTEGER(kind=int64), INTENT(IN) :: m,n,val
        INTEGER(kind=int64), INTENT(OUT) :: h_list(4)
        INTEGER(kind=int64) :: state(m*n), temp_state(m*n), k,a,b,c,d,loc
        INTEGER(kind=int64) :: x1,y1,x2,y2,x3,y3,x4,y4, i, j
        state = 0
        state(val+1) = 1
        !PRINT "(10i0)", state
        !h_list = -1
        !count = 1
        k = val
        call Inv_Map(k,n,i,j)
        !PRINT "(3i0)", i, j, k
        call nnn(i,j,m,n,x1,y1,x2,y2,x3,y3,x4,y4)
        !print "(8i2)",x1,y1,x2,y2,x3,y3,x4,y4
        a = Map(x1,y1,n)+1
        b = Map(x2,y2,n)+1
        c = Map(x3,y3,n)+1
        d = Map(x4,y4,n)+1

       
        h_list(1) = a

        h_list(2) = b

        h_list(3) = c

        h_list(4) = d
        
    end subroutine

    !generates hamiltonian matrix of NNN hopping for given m,n of the lattice
    function nnn_matrix(m,n, t) result(arr)
        implicit none
        INTEGER(kind=int64), INTENT(IN) :: m,n, t
        INTEGER(kind=int64) :: dim, i, j, temp_list(4), k
        REAL(kind=real64), ALLOCATABLE :: arr(:,:)
        dim = m*n
        ALLOCATE (arr(dim,dim))
        arr = 0 
        do i  = 1, dim
            CALL nnn_hop_list(m,n,i-1,temp_list)
            !print "(4i0)", temp_list
            do j = 1, 4
                arr(i,temp_list(j)) = t
            end do
            
        end do
        arr = 1_8*arr
        do i = 1, dim
            arr(i,i) = 0
        end do
        return


    end function





end program