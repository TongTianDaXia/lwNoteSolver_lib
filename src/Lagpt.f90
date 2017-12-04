        module  Laguerre
        
    contains
    
            pure subroutine lagpts(n)
            integer,intent(inout)::    n
            n = 2*n
            end subroutine lagpts
            
        end module 
    
        