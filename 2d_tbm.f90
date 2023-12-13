



program main
    use iso_fortran_env
    implicit none
    INTEGER(kind=8) :: arr0(5) = [2,4,6,8,10] , i, j, n = 3, m = 10
    INTEGER(kind=int64) :: arr1(12)  = [1,1,1,1,1,1,0,0,0,0,0,0]
    INTEGER(kind=int64)  :: state(10), val, k, hh_list(48), matrix(100,100)
    INTEGER(kind=int64) :: x0,y0,x1,y1,x2,y2,x3,y3,x4,y4

    !print "(6I2)", arr0
    !CALL swap(arr0, 1_8, 2_8)
    !print "(6I2)", arr0


    !print "(4I0)", Map(1, 0, 3)


    !i = 0,...,m-1 and j = 0,...,n-1
    !CALL Inv_Map(Map(0_8, 1_8, n), n, i, j)
    !print "(1a,g0,1a,1a,g0)","i=", i, ",", "j=", j
    !call Number2POB(4, 2, 0, arr1)
    !do k = 0 , comb(10_8,5_8)-1
        !state = number2pob(10,5,k)
        !val = pob2number(state, SIZE(state))
        !print "(10i0,1a,g0)", state,"~", val
    !end do
    !val = pob2number(arr1,12)
    !!print *, val
    !call nnn_hop_list(5_8,2_8,251_8, hh_list)
    !call nnn(4_8,1_8,5_8,2_8,x1,y1,x2,y2,x3,y3,x4,y4)
    !PRINT *,Map(3,0,2)
    !call nn_hop_list(4,3,0,hh_list)
    !state = number2pob(10,5,0)
    !print *, pob2number(state,10)
    !print "(10i0)",state
    !call swap(state,1,10)
    !print *, pob2number(state,10)
    !print "(10i0)",state
    matrix = nn_matrix(3_8,3_8,1_8)
    write(1,'(100(",",G0,:))') matrix
    matrix = nnn_matrix(3_8,3_8,1_8)
    write(2,'(100(",",G0,:))') matrix
    !print "(924i2)", matrix






    
   

  
    















    contains


    !the subroutines
    subroutine swap(arr, i, j)
        implicit none
        INTEGER(kind=int64), INTENT(IN) :: i, j
        INTEGER(kind=int64), INTENT(INOUT) :: arr(:)
        INTEGER(kind=int64) :: temp
        temp = arr(i)
        arr(i) = arr(j)
        arr(j) = temp
        
    end subroutine swap


    !i={0,...,m-1}, j={0,...,n-1}
    function Map(i, j, n)
        implicit none
        INTEGER(kind=int64) :: Map
        INTEGER(kind=int64) :: i, j, n
        Map = i*n + j
        RETURN
    end function
    !val={0,...,m*n-1}, i={0,...,m-1}, j={0,...,n-1}
    subroutine Inv_Map(val, n, i, j)
        INTEGER(kind=int64), INTENT(IN) :: val, n
        INTEGER(kind=int64), INTENT(OUT) :: i, j
        j = MOD(val, n)
        i = (val - j)/n
    end subroutine Inv_Map
    
    
    !reverse of a 1d array
    subroutine reverse(arr, len)
        INTEGER(kind=int64), INTENT(IN) :: len
        INTEGER(kind=int64) :: arr(len), reverse_arr(len)
        do i = 1, len
            reverse_arr(i) = arr(len-i+1)
        end do
        arr = reverse_arr
    end subroutine

    !this function cleans an array marked with -1
    function clean_array(arr,len) 
        INTEGER(kind=int64) :: len
        INTEGER(kind=int64) :: arr(len), i, count, temp(len)
        INTEGER(kind=int64), ALLOCATABLE :: clean_array(:)
        temp = arr
        do i = 1, len
            if (temp(i)==-1) then 
                count = i-1
                exit
            end if
        end do
        ALLOCATE (clean_array(count))
        clean_array = temp(1:count)
        return
    end function


    
    
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
    
    
    !append an element to the end of an 1d array
    function append(arr,len, elem)
        INTEGER(kind=int64) :: len, elem
        INTEGER(kind=int64) :: arr(len), append(len)
        append(1:len) = arr
        append(len+1) = elem
        RETURN
    end function
    
    !Function that converts POB representation of state to a unique INTEGER(kind=int64) "pob2integer" = 1,...,comb(mn,mn/2)
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


    !Function that converts unique INTEGER(kind=int64) value v t0 pob representation as a 1d array of size=n, with r 1s...
    function number2pob(n, r, v) 
        implicit none

        INTEGER(kind=int64), intent (in) :: n
        INTEGER(kind=int64), intent (in) :: r
        INTEGER(kind=int64), intent (in) :: v

        INTEGER(kind=int64), dimension(n) :: number2pob

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
                    number2pob(j+1) = 1 !.true.
                    
                else
                    number2pob(j+1) = 0 !.false.
                end if
                if (number2pob(j+1) == 1) then
                    exit
                end if
            end do
        end do
        
        if (j>0) then
            do k = j, 1, -1
                number2pob(k) = 0
            end do
        end if
        temp1 = SIZE(number2pob)
        CALL reverse(number2pob, temp1)
        RETURN
    end function number2pob



    !given m,n and i,j returns coordinates enforcing pbc hops...
    subroutine pbc_nn(i,m,j,n)
        implicit none
        INTEGER(kind=int64), INTENT(IN) :: m, n
        INTEGER(kind=int64), INTENT(INOUT) :: i, j
        if (i<0) then
            i = m-1
        else if (i>m-1) then
            i  = 0
        else
            CONTINUE
        end if

        if (j<0) then
            j= n-1
        else if (j>n-1) then
            j  = 0
        else
            CONTINUE
        end if
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

    !given m,n and i,j gives its NNN coordinates with pbc
    !some serious issues in enforcing pbc
    !hence its is non-pbc currently...
    subroutine nnn(i,j,m,n,x1,y1,x2,y2,x3,y3,x4,y4)
        implicit none
        INTEGER(kind=int64), INTENT(IN) :: i, j, m, n
        INTEGER(kind=int64), INTENT(INOUT) :: x1,y1,x2,y2,x3,y3,x4,y4
        INTEGER(kind=int64) :: temp
        x1 = i+1
        y1 = j+1
        if (x1>m-1 .and. y1>n-1) then
            x1=i
            y1=j
        end if
        if (y1>n-1) then
            temp = x1
            x1 = i
            y1 = j
        else if (x1>m-1) then
            temp = y1
            y1 = j
            x1 = i
        end if


        !print "(i0,a1,i0)", x1,",",y1

        x2 = i + 1
        y2 = j - 1
        if (x2 > m-1 .and. y2 < 0) then
            x2 = i
            y2 = j
        end if
        if (y2 < 0) then
            temp = x2
            x2 = i
            y2 = j
        else if (x2 > m-1) then
            temp = y2
            y2 = j
            x2 = i
        end if
        
        !print "(i0,a1,i0)", x2,",",y2

        x3 = i-1
        y3 = j-1
        if (x3 < 0 .and. y3 < 0) then
            x3 = i
            y3 = j
        end if
        if (y3 < 0) then
            temp = x3
            x3 = i
            y3 = j
        else if (x3 < 0) then
            temp = y3
            y3 = j
            x3 = i
        end if



        !print "(i0,a1,i0)", x3,",",y3

        x4 = i-1
        y4 = j+1
        if (x4 < 0 .and. y4 > n-1) then
            x4 = i
            y4 = j
        end if
        if (y4 > n-1) then
            temp = x4
            x4 = i
            y4 = j
        else if (x4 < 0) then
            temp = y4
            y4 = j
            x4 = i
        end if

        !print "(i0,a1,i0)", x4,",",y4
    end subroutine

    !pob_val = {0,...,comb(m*n,m*n/2)-1}
    subroutine nn_hop_list(m,n,pob_val,h_list)
        implicit none
        INTEGER(kind=int64), INTENT(IN) :: m,n,pob_val
        INTEGER(kind=int64), INTENT(OUT) :: h_list(4*m*n)
        INTEGER(kind=int64) :: state(m*n), temp_state(m*n), k,a,b,c,d,count
        INTEGER(kind=int64) :: x1,y1,x2,y2,x3,y3,x4,y4, i, j
        state = number2pob(m*n,m*n/2,pob_val)
        h_list = -1
        count = 1
        do k = 1,m*n
            if (state(k)==1) then
                call Inv_Map(k-1,n,i,j)
                !print "(12i0)",state
                !call reverse(state, m*n)
                call nn(i,j,m,n,x1,y1,x2,y2,x3,y3,x4,y4)
                a = Map(x1,y1,n)+1
                b = Map(x2,y2,n)+1
                c = Map(x3,y3,n)+1
                d = Map(x4,y4,n)+1
                !print "(2i2)", i, j
                !print "(8i2)", x1,y1,x2,y2,x3,y3,x4,y4

                if (state(a)==0) then
                    temp_state = state
                    call swap(temp_state,a,k)
                    h_list(count) = pob2number(temp_state,m*n)
                    count = count + 1
                    !print *, "a"
                    !print "(g0,1a,g0)",a,",",k
                    !print "(10i0)",temp_state
                    !print "(g0)", pob2number(temp_state,m*n)
                end if
                if (state(b)==0) then
                    temp_state = state
                    call swap(temp_state,b,k)
                    h_list(count) = pob2number(temp_state,m*n)
                    count = count + 1
                    !print *, "b"
                    !print "(g0,1a,g0)",b,",",k
                    !print "(10i0)",temp_state
                    !print "(g0)", pob2number(temp_state,m*n)
                end if
                if (state(c)==0) then
                    temp_state = state
                    call swap(temp_state,c,k)
                    h_list(count) = pob2number(temp_state,m*n)
                    count = count + 1
                    !print *, "c"
                    !print "(g0,1a,g0)",c,",",k
                    !print "(10i0)",temp_state
                    !print "(g0)", pob2number(temp_state,m*n)
                end if
                if (state(d)==0) then
                    temp_state = state
                    call swap(temp_state,d,k)
                    h_list(count) = pob2number(temp_state,m*n)
                    count = count + 1
                    !print *, "d"
                    !print "(g0,1a,g0)",d,",",k
                    !print "(10i0)",temp_state
                    !print "(g0)", pob2number(temp_state,m*n)
                end if
            end if
        end do
        !print *, clean_array(h_list,4*m*n)

        
    end subroutine


    
    !pob_val = {0,...,comb(m*n,m*n/2)-1}
    subroutine nnn_hop_list(m,n,pob_val,h_list)
        implicit none
        INTEGER(kind=int64), INTENT(IN) :: m,n,pob_val
        INTEGER(kind=int64), INTENT(OUT) :: h_list(4*m*n)
        INTEGER(kind=int64) :: state(m*n), temp_state(m*n), k,a,b,c,d,count
        INTEGER(kind=int64) :: x1,y1,x2,y2,x3,y3,x4,y4, i, j
        state = number2pob(m*n,m*n/2,pob_val)
        h_list = -1
        count = 1
        do k = 1,m*n
            if (state(k)==1) then
                call Inv_Map(k-1,n,i,j)
                !print "(2i0)",state
                !print *, "-----"
                !call reverse(state, m*n)
                call nnn(i,j,m,n,x1,y1,x2,y2,x3,y3,x4,y4)
                a = Map(x1,y1,n)+1
                b = Map(x2,y2,n)+1
                c = Map(x3,y3,n)+1
                d = Map(x4,y4,n)+1
                !print "(2i2)", i, j
                !print "(8i2)", x1,y1,x2,y2,x3,y3,x4,y4

                if (state(a)==0) then
                    temp_state = state
                    call swap(temp_state,a,k)
                    h_list(count) = pob2number(temp_state,m*n)
                    count = count + 1
                    !print *, "a"
                    !print "(g0,1a,g0)",a,",",k
                    !print "(10i0)",temp_state
                    !print "(g0)", pob2number(temp_state,m*n)
                end if
                if (state(b)==0) then
                    temp_state = state
                    call swap(temp_state,b,k)
                    h_list(count) = pob2number(temp_state,m*n)
                    count = count + 1
                    !print *, "b"
                    !print "(g0,1a,g0)",b,",",k
                    !print "(10i0)",temp_state
                    !print "(g0)", pob2number(temp_state,m*n)
                end if
                if (state(c)==0) then
                    temp_state = state
                    call swap(temp_state,c,k)
                    h_list(count) = pob2number(temp_state,m*n)
                    count = count + 1
                    !print *, "c"
                    !print "(g0,1a,g0)",c,",",k
                    !print "(10i0)",temp_state
                    !print "(g0)", pob2number(temp_state,m*n)
                end if
                if (state(d)==0) then
                    temp_state = state
                    call swap(temp_state,d,k)
                    h_list(count) = pob2number(temp_state,m*n)
                    count = count + 1
                    !print *, "d"
                    !print "(g0,1a,g0)",d,",",k
                    !print "(10i0)",temp_state
                    !print "(g0)", pob2number(temp_state,m*n)
                end if
            end if
        end do
        !print *, clean_array(h_list,4*m*n)

        
    end subroutine




    function nn_matrix(m,n, t) result(arr)
        implicit none
        INTEGER(kind=int64), INTENT(IN) :: m,n, t
        INTEGER(kind=int64) :: dim, i, j, temp_list(4*m*n), k
        INTEGER(kind=int64), ALLOCATABLE :: arr(:,:), clean_list(:)
        dim = comb(m*n,m*n/2)
        ALLOCATE (arr(dim,dim))
        arr = 0 
        do i  = 1, dim
            CALL nn_hop_list(m,n,i-1,temp_list)
            k = size(clean_array(temp_list, 4*m*n))
            ALLOCATE (clean_list(k))
            clean_list = clean_array(temp_list, 4*m*n)
            !print *, clean_list
            do j = 1, 4*m*n
                if (temp_list(j) ==-1) then
                    exit
                else
                    arr(i,temp_list(j)+1) = t
                end if
            end do
            DEALLOCATE (clean_list)
        end do
        return

    end function



    



        


    function nnn_matrix(m,n, t) result(arr)
        implicit none
        INTEGER(kind=int64), INTENT(IN) :: m,n, t
        INTEGER(kind=int64) :: dim, i, j, temp_list(4*m*n), k
        INTEGER(kind=int64), ALLOCATABLE :: arr(:,:), clean_list(:)
        dim = comb(m*n,m*n/2)
        ALLOCATE (arr(dim,dim))
        arr = 0 
        do i  = 1, dim
            CALL nnn_hop_list(m,n,i-1,temp_list)
            k = size(clean_array(temp_list, 4*m*n))
            ALLOCATE (clean_list(k))
            clean_list = clean_array(temp_list, 4*m*n)
            !print *, clean_list
            do j = 1, 4*m*n
                if (temp_list(j) ==-1) then
                    exit
                else
                    arr(i,temp_list(j)+1) = t
                end if
            end do
            DEALLOCATE (clean_list)
        end do
        return
        end function

        



    
      
        
      
      
      
      
      



end program main