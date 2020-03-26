!all of shape is assumed convex
module GeometryLib
use constants
use arrayOpsLib
implicit none

    private
    
    public:: segmentLength
    public:: geoArea, trilateralArea, quadrilateralArea, polygonArea
    public:: tetrahedronVolume
    !--
    public:: geometryShape
    public:: segment, trilateral, quadrilateral, tetrahedron, hexahedron, prism
    
    !-----------------------------------------------------
    !enum for common geommetry
    integer(ip),parameter:: segment = -999
    integer(ip),parameter:: trilateral = -1000
    integer(ip),parameter:: quadrilateral = -1001
    integer(ip),parameter:: tetrahedron = -1002
    integer(ip),parameter:: hexahedron = -1003
    integer(ip),parameter:: prism = -1004
    
    
    interface geoArea
        procedure:: trilateralArea
        procedure:: quadrilateralArea
        procedure:: polygonArea
    end interface geoArea
    
    interface geoVol
        procedure:: tetrahedronVolume
        procedure:: hexVol
    end interface geoVol
    
    
contains


    !-----------------------------------------------
    !geometry in 1d and above
    !-----------------------------------------------
    pure real(rp) function segmentLength(s,e)
    real(rp),dimension(:),intent(in)::  s,e
        segmentLength = norm2(s-e)
    end function segmentLength

    
    !-----------------------------------------------
    !geometry in 2d and above
    !-----------------------------------------------
    pure real(rp) function trilateralArea(p1,p2,p3)
    real(rp),dimension(:),intent(in)::  p1,p2,p3
    real(rp)::                          s1,s2,s3,a
        s1 = segmentLength(p1,p2)
        s2 = segmentLength(p2,p3)
        s3 = segmentLength(p1,p3)
        a = (s1+s2+s3)*0.5_rp
        trilateralArea = sqrt(a*(a-s1)*(a-s2)*(a-s3))
    end function trilateralArea
    
    pure real(rp) function quadrilateralArea(p1,p2,p3,p4)
    real(rp),dimension(:),intent(in)::  p1,p2,p3,p4
        quadrilateralArea = trilateralArea(p1,p2,p3) + trilateralArea(p1,p3,p4)
    end function quadrilateralArea
    
    !the order of coordinate is important, p(3,n) is the coordinate
    pure real(rp) function polygonArea(p)
    real(rp),dimension(:,:),intent(in)::        p 
    integer(ip)::                               n,i
        n = size(p,2); polygonArea = 0._rp
        do i=2,n-1
            polygonArea = polygonArea + trilateralArea(p(:,1), p(:,i), p(:,i+1))
        enddo
        polygonArea = abs(polygonArea)
    end function polygonArea
    
    
    !---------------------------------------------------
    !geometry in 3d
    !---------------------------------------------------
    !free of order of (p1,p2,p3,p4)
    !v = 1/6 det{p1-p4,p2-p4,p3-p4} = |(p1-p4) .ip. ((p2-p4) .cp. p3-p4)|/6
    pure real(rp) function tetrahedronVolume(p1,p2,p3,p4)
    real(rp),dimension(:),intent(in)::  p1,p2,p3,p4
        tetrahedronVolume = (p1-p4) .ip. ((p2-p4) .cpv. (p3-p4))
        tetrahedronVolume = abs(tetrahedronVolume)/6._rp
    end function tetrahedronVolume
    
    
    !--
    !p1-p2-p3-p4
    !p5-p6-p7-p8
    pure real(rp) function hexVol(p1,p2,p3,p4,p5,p6,p7,p8)
    real(rp),dimension(:),intent(in)::  p1,p2,p3,p4,p5,p6,p7,p8
        call disableprogram()
        hexVol = 0._rp
    end function hexVol
    
    
    !------------------------------------------------------
    pure function geometryShape(enumShape)
    integer(ip),intent(in)::    enumshape
    character(cl)::             geometryShape
        select case(enumshape)
        case(segment)
            geometryShape = 'segment'
        case(trilateral)
            geometryShape = 'trilateral'
        case(quadrilateral)
            geometryShape = 'quadrilateral'
        case(tetrahedron)
            geometryShape = 'tetrahedron'
        case(hexahedron)
            geometryShape = 'hexahedron'
        case(prism)
            geometryShape = 'prism'
        case default
            call disableprogram
        end select
    end function geometryShape
    
end module GeometryLib