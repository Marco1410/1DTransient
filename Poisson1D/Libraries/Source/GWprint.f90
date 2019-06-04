Module Print_All
    Use tools, Only: dp, sp
    Implicit none

    Interface PrintAll
        Module Procedure Print_IntVec, Print_IntMat
        Module Procedure Print_Real4Vec, Print_Real4Mat
        Module Procedure Print_Real8Vec, Print_Real8Mat
    End Interface PrintAll


Contains

    Subroutine Print_IntVec(vec,dim)
        Implicit none
        Integer,Intent(IN), Dimension(:):: vec
        Integer,Intent(In), Optional::dim
        Integer i
        Write(*,"(/,2x,A)") "Vector information"
        Write(*,"(2x,A)") "Integer type"
        Write(*,"(2x,A,x,I3)") "Total components:",size(vec)
        If(Present(dim))Then
            Write(*,"(2x,A,x,I3)") "Components used:",dim
            Do i=1, dim
                Write(*,"(*(I5))") vec(i)
             End Do
             Print*,
        Else
            Write(*,"(2x,A,x,I3)") "Components used:",size(vec)
            Do i=lbound(vec,1), ubound(vec,1)
                Write(*,"(*(I5))") vec(i)
             End Do
             Print*,
        End If
    End Subroutine Print_IntVec

    Subroutine Print_IntMat(mat, row, col)
        Implicit none
        Integer,Intent(IN), Dimension(:,:):: mat
        Integer,Intent(In), Optional:: row, col
        Integer i, j
        Write(*,"(/,2x,A)") "Matrix information"
        Write(*,"(2x,A)") "Integer type"
        Write(*,"(2x,A,x,I3,I3)") "Total Shape (row, col):",shape(mat)
        If(Present(row) .and. Present(col)) Then
            Write(*,"(2x,A,x,I3,I3)") "Shape used (row, col):",row, col
            Do i=1,row
                Write(*,"(*(I5))")( mat(i,j), j=1, col)
             End Do
             Print*, 
        Else if(Present(row) .or. Present(col)) Then
            Write(*,"(2x,A)") "Matrix assumed square, row=col"
            Do i=1,row
                Write(*,"(*(I5))")( mat(i,j), j=1, row)
             End Do
             Print*, 
        Else
            Write(*,"(2x,A,x,I3,I3)") "Shape used (row, col):",shape(mat)
            Do i=lbound(mat,1),ubound(mat,1)
                Write(*,"(*(I5))")( mat(i,j), j=lbound(mat,2),ubound(mat,2))
             End Do
             Print*, 
        End If
    End Subroutine Print_IntMat

    Subroutine Print_Real4Vec(vec, dim)
        Implicit none
        Real(sp),Intent(IN), Dimension(:):: vec
        Integer, Intent(In), Optional:: dim
        Integer i
        Write(*,"(/,2x,A)") "Vector information"
        Write(*,"(2x,A)") "Real4 type"
        Write(*,"(2x,A,x,I3)") "Total components:",size(vec)
        If(Present(dim))Then
            Write(*,"(2x,A,x,I3)") "Components used:",dim
            Do i=1, dim
                Write(*,"(*(e13.5))") vec(i)
             End Do
             Print*, 
        Else
            Write(*,"(2x,A,x,I3)") "Components used:",size(vec)
            Do i=lbound(vec,1), ubound(vec,1)
                Write(*,"(*(e13.5))") vec(i)
             End Do
             Print*, 
        End If
    End Subroutine Print_Real4Vec

    Subroutine Print_Real4Mat(mat, row, col)
        Implicit none
        Real(sp),Intent(IN), Dimension(:,:):: mat
        Integer, Intent(In), Optional:: row, col
        Integer i, j
        Write(*,"(/,2x,A)") "Matrix information"
        Write(*,"(2x,A)") "Real4 type"
        Write(*,"(2x,A,x,I3,I3)") "Total Shape (row, col):",shape(mat)
        If(Present(row) .and. Present(col)) Then
            Do i=1,row
                Write(*,"(*(e13.5))")( mat(i,j), j=1,col)
             End Do
             Print*,
        Else if(Present(row) .or. Present(col)) Then
            Write(*,"(2x,A)") "Please set the correct shape, row & col"
            Do i=1,row
                Write(*,"(*(e13.5))")( mat(i,j), j=1,row)
             End Do
             Print*,
        Else
            Do i=lbound(mat,1),ubound(mat,1)
                Write(*,"(*(e13.5))")( mat(i,j), j=lbound(mat,2),ubound(mat,2))
             End Do
             Print*,
        End If
    End Subroutine Print_Real4Mat

    Subroutine Print_Real8Vec(vec, dim)
        Implicit none
        Real(dp),Intent(IN), Dimension(:):: vec
        Integer, Intent(In), Optional:: dim
        Integer i
        Write(*,"(/,2x,A)") "Vector information"
        Write(*,"(2x,A)") "Real8 type"
        Write(*,"(2x,A,x,I3)") "Total components:",size(vec)
        If(Present(dim))Then
            Write(*,"(2x,A,x,I3)") "Components used:",dim
            Do i=1, dim
                Write(*,"(*(e13.5))") vec(i)
             End Do
             Print*,
        Else
            Write(*,"(2x,A,x,I3)") "Components used:",size(vec)
            Do i=lbound(vec,1), ubound(vec,1)
                Write(*,"(*(e13.5))") vec(i)
             End Do
             Print*,
        End If
    End Subroutine Print_Real8Vec

    Subroutine Print_Real8Mat(mat, row, col)
        Implicit none
        Real(dp),Intent(IN), Dimension(:,:):: mat
        Integer, Intent(In), Optional:: row, col
        Integer i, j
        Write(*,"(/,2x,A)") "Matrix information"
        Write(*,"(2x,A)") "Real8 type"
        Write(*,"(2x,A,x,I3,I3)") "Total Shape (row, col):",shape(mat)
        If(Present(row) .and. Present(col)) Then
            Do i=1,row
                Write(*,"(*(e13.5))")( mat(i,j), j=1,col)
             End Do
             Print*,
        Else if(Present(row) .or. Present(col)) Then
            Write(*,"(2x,A)") "Matrix assumed square, row=col"
            Do i=1,row
                Write(*,"(*(e13.5))")( mat(i,j), j=1,row)
             End Do
             Print*,
        Else
            Do i=lbound(mat,1),ubound(mat,1)
                Write(*,"(*(e13.5))")( mat(i,j), j=lbound(mat,2),ubound(mat,2))
             End Do
             Print*,
        End If
    End Subroutine Print_Real8Mat

End Module Print_All




